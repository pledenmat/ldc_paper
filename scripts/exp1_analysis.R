#' @title Analysis of the data from experiment 1
#' @author Pierre Le Denmat
#' @description Analysis of the data from experiment 1 in the context of the research article: 
#' "A low-dimensional approximation of optimal confidence".
#' @details The script is divided into several sections:
#' - Preprocessing: The data is loaded and preprocessed.
#' - Model fitting: Load the model fits performed on the cluster (or fit the model if not already done).
#' - Simulation from the model: Simulate data from the best fitting parameters
#' - Parameter analysis: Look at condition effect on the estimated parameters
#' - Simulation analysis: Reproduce the behavioral analysis on simulated data
#' - Plots: Generate the plots for the paper


rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(Rcpp) 
source('../functions/quantilefit_function_fullRTconf.R')
library(DEoptim)
library(MALDIquant)
library(reshape)
library(OneR)
library(scales)
library(lmerTest)
library(fields) 
library(prob)
library(car)
library(multcomp)
library(effectsize)
setwd(curdir)

stat_test <- F

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Preprocessing -----------------------------------------------------------

N <- 50
for(i in 1:N){
  if(i == 1){
    Data1 <- read.csv(paste0('../data/exp1/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('../data/exp1/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data1 <- rbind(Data1,temp)
  }
}

head(Data1)

Data1 <- subset(Data1, running == "main")

Data1['response'] <- 0
Data1$response[Data1$resp == "['n']"] <- 1

Data1$rt <- Data1$rt/1000
Data1$RTconf <- Data1$RTconf/1000

Data1 <- subset(Data1,rt>.1&rt<4) #There was no trial below 200ms

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data1$sub)); subs <- unique(Data1$sub); exclusion <- c()
tasks <- unique(Data1$task)
par(mfrow=c(2,2))
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data1,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (max(tempDat$rt)>3000) {
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT")
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
    }
  }
}

Data1 <- subset(Data1,!(sub %in% exclusion)) #Sub 10 and 49 removed

Data1$response[Data1$response==0] <- -1

Data1 <- Data1[,c("sub","task","selfconf","difflevel","rt","response","cor","cj","RTconf","block")]
names(Data1) <- c("sub","task","selfconf","coh","rt","resp","cor","cj","RTconf","block")
write.csv(df,"../data/dataexp1_prior_belief.csv",row.names = FALSE) # Export data to fit the model on the cluster

condLab <- unique(Data1$selfconf); Ncond <- length(condLab) #Change column name according to condition label in dataset
subs <- unique(Data1$sub);N<-length(subs)

difficulty <- sort(unique(Data1$coh));
# Fits Load ---------------------------------------------------------------
bound <- matrix(NA,N,Ncond);v <- matrix(NA,N,Ncond);ter <- matrix(NA,N,Ncond);conf_rt <- matrix(NA,N,Ncond);alpha <- matrix(NA,N,Ncond);beta <- matrix(NA,N,Ncond); resid <- matrix(NA,N,Ncond)
v2 <- matrix(NA,N,Ncond);v3 <- matrix(NA,N,Ncond); #Adjust the number of drift parameters to the model loaded
vratio <- matrix(NA,N,Ncond)

for(i in 1:N){
  for(c in 1:Ncond){
    print(paste('Running participant',i,'from',N,"condition",c))
    load(paste0('../fits/exp1/results_sub_',subs[i],'_',condLab[c],'.Rdata'))
    bound[i,c] <- results$optim$bestmem[1]
    ter[i,c] <- results$optim$bestmem[2]
    vratio[i,c] <- results$optim$bestmem[7]
    alpha[i,c] <-   results$optim$bestmem[8]
    beta[i,c] <-   results$optim$bestmem[9]
    v[i,c] <- results$optim$bestmem[10]
    v2[i,c] <-   results$optim$bestmem[11]
    v3[i,c] <-   results$optim$bestmem[12]
    resid[i,c] <- results$optim$bestval
  }
}

param_1 <- data.frame(drift = c(v,v2,v3),bound=rep(bound,3),ter=rep(ter,3),
                      alpha=rep(alpha,3),beta=rep(beta,3),sub=rep(subs,3*Ncond),
                      condition=rep(condLab,each=N,length.out=N*Ncond*3),
                      difflevel=rep(difficulty,each=N*Ncond),resid=rep(resid,3),
                      vratio=rep(vratio,3))

# Simulation from model ---------------------------------------------------

#Generate model simulations
rm(Simuls1);nsim <- 1

for(i in 1:N){
  print(paste('simulating',i,'from',N))
  for(c in 1:Ncond){
    temp_dat <- subset(Data1,sub==subs[i]&selfconf==condLab[c])
    ntrials <- with(temp_dat,aggregate(sub,by=list(difflevel=coh),table))
    temp <- data.frame(chi_square_optim(c(bound[i,c],ter[i,c],0,nsim,.1,.001,vratio[i,c],alpha[i,c],beta[i,c],v[i,c],v2[i,c],v3[i,c]),
                                        temp_dat,0))
    #Match to the number of trials per cell in the data
    temp <- rbind(subset(temp,drift==v[i,c])[sample(ntrials[ntrials$difflevel=="average","x"]),],
                  subset(temp,drift==v2[i,c])[sample(ntrials[ntrials$difflevel=="easy","x"]),],
                  subset(temp,drift==v3[i,c])[sample(ntrials[ntrials$difflevel=="hard","x"]),])
    if(!exists('Simuls1')){ Simuls1 <- cbind(temp,rep(condLab[c],nsim),rep(subs[i],nsim))
    }else{ Simuls1 <- rbind(Simuls1,cbind(temp,rep(condLab[c],nsim),rep(subs[i],nsim)))
    }
  }
}
names(Simuls1) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','cj_1','condition','sub')
Simuls1$exp <- 1

cj1 <- with(Data1,aggregate(cj,by=list(sub=sub,condition=selfconf,coh=coh),mean)); cj1 <- cast(cj1, sub~condition+coh)
cj1 <- c(cj1$lowSC_average,cj1$highSC_average,cj1$mediumSC_average,cj1$lowSC_easy,cj1$highSC_easy,cj1$mediumSC_easy,cj1$lowSC_hard,cj1$highSC_hard,cj1$mediumSC_hard)
param_1$cj <- cj1

coherences <- sort(unique(Data1$coh))
Simuls1$coh <- 0
for (i in 1:N) {
  for(d in 1:length(coherences)) Simuls1$coh[Simuls1$sub==subs[i] & Simuls1$drift %in% c(unique(subset(Simuls1,sub==subs[i])$drift)[d],unique(subset(Simuls1,sub==subs[i])$drift)[d+3],unique(subset(Simuls1,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}

# Parameter analysis ------------------------------------------------------
if (stat_test) {
  ##Drift rate
  m <- lmer(drift ~ condition*difflevel + (1|sub),data=param_1, REML = F)
  mcondition <- lmer(drift ~ condition*difflevel + (1 + condition|sub),data=param_1, REML = F)
  anova(m,mcondition)
  mdifflevel <- lmer(drift ~ condition*difflevel + (1 + difflevel|sub),data=param_1, REML = F) #Singular
  anova(mcondition)
  
  ##Bound
  m <- lmer(bound ~ condition + (1 |sub),data=param_1)
  anova(m) #There's a main effect here
  post.hoc <- glht(m, linfct = mcp(condition = 'Tukey'))
  summary(post.hoc)
  
##Non-decision time
  m <- lmer(ter ~ condition + (1 |sub),data=param_1)
  anova(m)
  
  ##alpha
  m <- lmer(alpha ~ condition + (1 |sub),data=param_1)
  anova(m)
  post.hoc <- glht(m, linfct = mcp(condition = 'Tukey'))
  summary(post.hoc)
  
  ##beta
  m <- lmer(beta ~ condition + (1 |sub),data=param_1)
  anova(m)
  post.hoc <- glht(m, linfct = mcp(condition = 'Tukey'))
  summary(post.hoc)
}
# Behavior analysis of prediction -------------------------------------------------
if (stat_test) {
  m.int <- lmer(rt~condition*coh+(1|sub),data=Simuls1,control = lmerControl(optimizer = "bobyqa"),REML = F)
  m.cond <- lmer(rt~condition*coh+(1+condition|sub),data=Simuls1,control = lmerControl(optimizer = "bobyqa"),REML = F)
  anova(m.int,m.cond)
  m.diff <- lmer(rt~condition*coh+(1+coh|sub),data=Simuls1,control = lmerControl(optimizer = "bobyqa"),REML = F)
  anova(m.diff,m.cond)
  m.both <- lmer(rt~condition*coh+(1+coh+condition|sub),data=Simuls1,control = lmerControl(optimizer = "bobyqa"),REML = F)
  anova(m.diff,m.both)
  plot(resid(m.both),Simuls1$rt) #Linearity
  leveneTest(residuals(m.both) ~ Simuls1$condition) #Homogeneity of variance
  qqmath(m.both) #Normality
  anova(m.both)
  
  m.cond <- lmer(cj~condition*coh+(1+condition|sub),data=Simuls1,REML = F)
  m.both <- lmer(cj~condition*coh+(1+condition+coh|sub),data=Simuls1,control = lmerControl(optimizer = "bobyqa"),REML = F)
  anova(m.cond,m.both)
  m.full <- lmer(cj~condition*coh+(condition*coh|sub),data=Simuls1,
                 control = lmerControl(optimizer = "bobyqa"),REML = F)
  anova(m.both,m.full)
  plot(resid(m.full),Simuls1$cj) #Linearity
  leveneTest(residuals(m.full) ~ Simuls1$condition) #Homogeneity of variance
  qqmath(m.full) #Normality
  anova(m.full)
  post.hoc <- glht(m.full, linfct = mcp(condition = 'Tukey'))
  summary(post.hoc)
  
  m.full <- lmer(cj~selfconf*coh*cor+(selfconf+coh+cor|sub),data=Data1,
                 control = lmerControl(optimizer = "bobyqa"),REML = F)
  anova(m.full)
  m.cor <- lmer(data=subset(Data1,cor==1),cj~selfconf*coh+(selfconf|sub),REML = F)
  anova(m.cor)
  m.error <- lmer(data=subset(Data1,cor==0),cj~selfconf*coh+(selfconf|sub),REML = F)
  anova(m.error)
  m.int <- glmer(cor~condition*coh+(1|sub),data=Simuls1,family="binomial",
                 control=glmerControl(optimizer="bobyqa"))
  m.cond <- glmer(cor~condition*coh+(1+condition|sub),data=Simuls1,family="binomial",
                  control=glmerControl(optimizer="bobyqa"))
  m.diff <- glmer(cor~condition*coh+(1+coh|sub),data=Simuls1,family="binomial",
                  control=glmerControl(optimizer="bobyqa"))
  anova(m.cond,m.diff)
  plot(resid(m.cond),Simuls1$cor) #Linearity
  leveneTest(residuals(m.cond) ~ Simuls1$cor) #Homogeneity of variance
  qqmath(m.cond) #Normality
  Anova(m.cond)
}
# Plot Layout -------------------------------------------------------------
setwd("../plots")
cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
cex.layout <- .66 
ratio <- .66/cex.layout # To get comparable sizes with figure 4
cex.lab <- cex_size(10,cex.layout) 
cex.lab.mtext <- cex.lab*cex.layout # Set with mtext so ignore mfrow
cex.ax <- cex_size(8,cex.layout)
cex.leg <- cex_size(8,cex.layout)
cex.datdot <- cex_size(8,cex.layout);
ylab_orientation <- 1 # 1 for horizontal, 0 for vertical
lwd.dat <- 1.5*ratio
mar.rt <- c(1,4,0,0)+.1
mar.acc <- c(5,4,0,0)+.1
mar.cj <- c(5,4,0,2)+.1
col_minus <- rgb(205,51,51,maxColorValue = 255)
col_minus_shade <- rgb(205,51,51,51,maxColorValue = 255)
col_baseline <- rgb(0,139,139,maxColorValue = 255)
col_baseline_shade <- rgb(0,139,139,51,maxColorValue = 255)
col_plus <- rgb(205,149,12,maxColorValue = 255)
col_plus_shade <- rgb(205,149,12,51,maxColorValue = 255)
tiff(file = "exp1_results.tif",width = 13,height=8,units="cm",compression="lzw",res = 1200)
font <- "Arial"
windowsFonts(A = windowsFont(font))
par(family="A")
layout(matrix(c(1,2,3,1,4,4),ncol=2),widths = c(1,2),heights = c(.1,.4,.4))
par(mar=c(0,0,0,0))
plot.new()
legend("top",legend=c("Positive","Average","Negative"),lty=1,bty = "n",horiz=T,
       title = "Fake Feedback condition",pch=rep(16,3),inset=.1, cex = cex.leg,
       col=c("darkgoldenrod3","cyan4","brown3"),y.intersp = .9)

# RT ----------------------------------------------------------------------

par(mar=mar.rt)

#Aggregate RT for data
rtlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(rt,by=list(sub,coh),mean));names(rtlow) <- c('sub','coh','rt')
rtlow <- cast(rtlow,sub~coh,)
rtmed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(rt,by=list(sub,coh),mean));names(rtmed) <- c('sub','coh','rt')
rtmed <- cast(rtmed,sub~coh)
rthigh <- with(subset(Data1,selfconf=="highSC"),aggregate(rt,by=list(sub,coh),mean));names(rthigh) <- c('sub','coh','rt')
rthigh <- cast(rthigh,sub~coh)

x <- rtlow[,c(2:4)];xmed <- rtmed[,c(2:4)];xhigh <- rthigh[,c(2:4)]
n <- length(x)

#Aggregate RT for model
simDat_low <- subset(Simuls1,condition=="lowSC");simDat_low <- simDat_low[c('sub','coh','rt')]
simDat_med <- subset(Simuls1,condition=="mediumSC");simDat_med <- simDat_med[c('sub','coh','rt')]
simDat_high <- subset(Simuls1,condition=="highSC");simDat_high <- simDat_high[c('sub','coh','rt')]

snrrtlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrrtmedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrrthighSim <- aggregate(.~sub+coh,simDat_high,mean)

x_sim = cast(snrrtlowSim,sub~coh,value='rt')
xmed_sim = cast(snrrtmedSim,sub~coh,value='rt')
xhigh_sim = cast(snrrthighSim,sub~coh,value='rt')

# Reorder difficulty from hard to easy
x <- x[,c("hard","average","easy")]
xmed <- xmed[,c("hard","average","easy")]
xhigh <- xhigh[,c("hard","average","easy")]
x_sim <- x_sim[,c("hard","average","easy")]
xmed_sim <- xmed_sim[,c("hard","average","easy")]
xhigh_sim <- xhigh_sim[,c("hard","average","easy")]

stripchart(x, ylim=c(.6,1.1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=cex.ax,family="A",las=ylab_orientation)
mtext("Response time (s)",2,at=.85,line=3,cex=cex.lab.mtext,family="A");
axis(1,at=0:(n-1),labels=c("","",""), cex.axis=cex.ax);
# Empirical data
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="brown3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwd.dat,length=0,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="cyan4",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwd.dat,length=0,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="darkgoldenrod3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwd.dat,length=0,col="darkgoldenrod3")
#Model predictions
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N)),(colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N)),(colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))
# ACC ---------------------------------------------------------------------

par(mar=mar.acc)

#Aggregate accuracy for data
corlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(cor,by=list(sub,coh),mean));names(corlow) <- c('sub','coh','cor')
corlow <- cast(corlow,sub~coh,)
cormed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(cor,by=list(sub,coh),mean));names(cormed) <- c('sub','coh','cor')
cormed <- cast(cormed,sub~coh)
corhigh <- with(subset(Data1,selfconf=="highSC"),aggregate(cor,by=list(sub,coh),mean));names(corhigh) <- c('sub','coh','cor')
corhigh <- cast(corhigh,sub~coh)

x <- corlow[,c(2:4)];xmed <- cormed[,c(2:4)];xhigh <- corhigh[,c(2:4)]
n <- length(x)


#Aggregate accuracy for model
simDat_low <- subset(Simuls1,condition=="lowSC");simDat_low <- simDat_low[c('sub','coh','cor')]
simDat_med <- subset(Simuls1,condition=="mediumSC");simDat_med <- simDat_med[c('sub','coh','cor')]
simDat_high <- subset(Simuls1,condition=="highSC");simDat_high <- simDat_high[c('sub','coh','cor')]

snrcorlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrcormedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrcorhighSim <- aggregate(.~sub+coh,simDat_high,mean)

x_sim = cast(snrcorlowSim,sub~coh,value='cor')
xmed_sim = cast(snrcormedSim,sub~coh,value='cor')
xhigh_sim = cast(snrcorhighSim,sub~coh,value='cor')

# Reorder difficulty from hard to easy
x <- x[,c("hard","average","easy")]
xmed <- xmed[,c("hard","average","easy")]
xhigh <- xhigh[,c("hard","average","easy")]
x_sim <- x_sim[,c("hard","average","easy")]
xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")]

# Start plotting
stripchart(x, ylim=c(.6,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=cex.ax,las=ylab_orientation)
mtext("Accuracy",2,at=.8,line=3,cex=cex.lab.mtext,family="A");
axis(1,at=0:(n-1),labels=c("Hard","Average","Easy"), cex.axis=cex.ax,family="A");
mtext("Trial difficulty",1,3,at=1,cex=cex.lab.mtext,family="A")
# Model predictions
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N)),(colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N)),(colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))
# Empirical data
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="brown3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwd.dat,length=0,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="cyan4",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwd.dat,length=0,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="darkgoldenrod3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwd.dat,length=0,col="darkgoldenrod3")
# CJ ----------------------------------------------------------------------

par(mar=mar.cj)

## Experiment 1
cjlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(cj,by=list(sub,coh,cor),mean));
names(cjlow) <- c('sub','coh','cor','cj')
cjlow_cor <- subset(cjlow,cor==1); cjlow_err <- subset(cjlow,cor==0)
cjlow_cor <- cast(cjlow_cor,sub~coh); cjlow_err <- cast(cjlow_err,sub~coh)
cjmed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(cj,by=list(sub,coh,cor),mean));
names(cjmed) <- c('sub','coh','cor','cj')
cjmed_cor <- subset(cjmed,cor==1); cjmed_err <- subset(cjmed,cor==0)
cjmed_cor <- cast(cjmed_cor,sub~coh); cjmed_err <- cast(cjmed_err,sub~coh)
cjhigh <- with(subset(Data1,selfconf=="highSC"),aggregate(cj,by=list(sub,coh,cor),mean));
names(cjhigh) <- c('sub','coh','cor','cj')
cjhigh_cor <- subset(cjhigh,cor==1); cjhigh_err <- subset(cjhigh,cor==0)
cjhigh_cor <- cast(cjhigh_cor,sub~coh); cjhigh_err <- cast(cjhigh_err,sub~coh)

#Simulations
simDat_low <- with(subset(Simuls1,condition=="lowSC"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_low) <-c('sub','coh','cor','cj')
simDat_low_cor <- subset(simDat_low,cor==1);simDat_low_err <- subset(simDat_low,cor==0)
simDat_low_cor <- cast(simDat_low_cor,sub~coh);simDat_low_err <- cast(simDat_low_err,sub~coh)
simDat_med <- with(subset(Simuls1,condition=="mediumSC"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_med) <- c('sub','coh','cor','cj')
simDat_med_cor <- subset(simDat_med,cor==1);simDat_med_err <- subset(simDat_med,cor==0)
simDat_med_cor <- cast(simDat_med_cor,sub~coh);simDat_med_err <- cast(simDat_med_err,sub~coh)
simDat_high <- with(subset(Simuls1,condition=="highSC"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_high) <- c('sub','coh','cor','cj')
simDat_high_cor <- subset(simDat_high,cor==1);simDat_high_err <- subset(simDat_high,cor==0)
simDat_high_cor <- cast(simDat_high_cor,sub~coh);simDat_high_err <- cast(simDat_high_err,sub~coh)


#aggregate cj for model
xsim_cor <- simDat_low_cor[,c(2:4)];xsimmed_cor <- simDat_med_cor[,c(2:4)];xsimhigh_cor <- simDat_high_cor[,c(2:4)]
xsim_err <- simDat_low_err[,c(2:4)];xsimmed_err <- simDat_med_err[,c(2:4)];xsimhigh_err <- simDat_high_err[,c(2:4)]
xminus_cor <- cjlow_cor[,c(2:4)];xbaseline_cor <- cjmed_cor[,c(2:4)];xplus_cor <- cjhigh_cor[,c(2:4)]
xminus_err <- cjlow_err[,c(2:4)];xbaseline_err <- cjmed_err[,c(2:4)];xplus_err <- cjhigh_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","average","easy")];
xbaseline_cor <- xbaseline_cor[,c("hard","average","easy")];
xplus_cor <- xplus_cor[,c("hard","average","easy")]
xminus_err <- xminus_err[,c("hard","average","easy")];
xbaseline_err <- xbaseline_err[,c("hard","average","easy")];
xplus_err <- xplus_err[,c("hard","average","easy")]
xsim_cor <- xsim_cor[,c("hard","average","easy")];
xsimmed_cor <- xsimmed_cor[,c("hard","average","easy")];
xsimhigh_cor <- xsimhigh_cor[,c("hard","average","easy")]
xsim_err <- xsim_err[,c("hard","average","easy")];
xsimmed_err <- xsimmed_err[,c("hard","average","easy")];
xsimhigh_err <- xsimhigh_err[,c("hard","average","easy")]

stripchart(xminus_cor,ylim=c(2.5,5.7), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "",
           main = NULL,family="A")
axis(1,at=0:(n-1),labels=c("Hard","Average","Easy"), cex.axis=cex.ax,family="A");
axis(2, seq(2.5,5.5,.5), cex.axis=cex.ax,las=ylab_orientation,family="A")
title(ylab = "Confidence", xlab= "Trial difficulty",line = 2.5,cex.lab=cex.lab)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsim_cor,na.rm=T) + (colSds(as.matrix(xsim_cor))/sqrt(N)),
          (colMeans(xsim_cor,na.rm=T) - colSds(as.matrix(xsim_cor))/sqrt(N))[3:1]),
        border=F,col=col_minus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimmed_cor,na.rm=T) + (colSds(as.matrix(xsimmed_cor))/sqrt(N)),
          (colMeans(xsimmed_cor,na.rm=T) - colSds(as.matrix(xsimmed_cor))/sqrt(N))[3:1]),
        border=F,col=col_baseline_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimhigh_cor,na.rm=T) + (colSds(as.matrix(xsimhigh_cor))/sqrt(N)),
          (colMeans(xsimhigh_cor,na.rm=T) - colSds(as.matrix(xsimhigh_cor))/sqrt(N))[3:1]),
        border=F,col=col_plus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsim_err,na.rm=T) + (colSds(as.matrix(xsim_err))/sqrt(N)),
          (colMeans(xsim_err,na.rm=T) - colSds(as.matrix(xsim_err))/sqrt(N))[3:1]),
        border=F,col=col_minus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimmed_err,na.rm=T) + (colSds(as.matrix(xsimmed_err))/sqrt(N)),
          (colMeans(xsimmed_err,na.rm=T) - colSds(as.matrix(xsimmed_err))/sqrt(N))[3:1]),
        border=F,col=col_baseline_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimhigh_err,na.rm=T) + (colSds(as.matrix(xsimhigh_err))/sqrt(N)),
          (colMeans(xsimhigh_err,na.rm=T) - colSds(as.matrix(xsimhigh_err))/sqrt(N))[3:1]),
        border=F,col=col_plus_shade)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N),lwd=lwd.dat,col=col_minus,length=0)
means <- sapply(xbaseline_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_cor),na.rm=T)/sqrt(N),lwd=lwd.dat,col=col_baseline,length=0)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N),lwd=lwd.dat,col=col_plus,length=0)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N),lwd=lwd.dat,col=col_minus,length=0)
means <- sapply(xbaseline_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_err),na.rm=T)/sqrt(N),lwd=lwd.dat,col=col_baseline,length=0)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N),lwd=lwd.dat,col=col_plus,length=0)
legend("bottomleft",legend=c("Correct trials","Error trials"),
       title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
       cex=cex.leg,lwd=lwd.dat,seg.len=1.5)
legend("topleft",legend=c("Empirical data","Model prediction"),
       title = "",bty = "n",inset=0,lty=c(1,NA),pch=c(16,15),pt.cex=c(1,2),
       cex=cex.leg,lwd=lwd.dat,seg.len=1.5,col="lightgrey")
dev.off()

# Plot estimated parameters -------------------------------------------------------------------
jpeg(filename = "estimated_par_exp1.tif",
     width = 19,
     height = 14,
     units = 'cm',
     res = 1200)
windowsFonts(A = windowsFont("Arial"))
layout(matrix(c(4,1,4,1,4,2,5,2,5,3,5,3),ncol=6))
diff_order <- c("hard","average","easy")
cond_order <- c("lowSC","mediumSC","highSC")
cex.layout <- .66
cexax <- cex_size(8,cex.layout)
cexlab <- cex_size(10,cex.layout)
cexmain <- cex_size(8,cex.layout)
cexleg <- cex_size(8,cex.layout)
linelab <- 2.5
lwdmean <- 3
col_indiv_points <- rgb(.7,.7,.7,.5)
jit_size <- .3
col_hard <- rgb(247,104,161,maxColorValue = 255)
col_med <- rgb(197,27,138,maxColorValue = 255)
col_easy <- rgb(122,1,119,maxColorValue = 255)
leg_squeeze_within <- .5
leg_squeeze_between <- 1
leg_text_width <- .5*c(1,1.5,1)
par(family="A")

# DDM parameters ----------------------------------------------------------


##Non-decision time
par(mar=c(5,3.5,1,2)+0.1)
plot_ter <- with(param_1,aggregate(ter,by=list(sub=sub,condition=condition),mean))
plot_ter <- cast(plot_ter,sub~condition)
plot_ter <- plot_ter[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_ter),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_ter),max(plot_ter)),xaxt='n',main=NULL,family="A")
title(ylab="Non-decision time",xlab="Feedback condition", line = linelab, 
      cex.lab = cexlab,family="A")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=cexax,family="A")
for(i in 1:N) lines(jitter(1:Ncond,jit_size),plot_ter[i,1:Ncond],type='b',lty=2,col=col_indiv_points,pch=19)
points(colMeans(plot_ter),type='b',lwd=lwdmean)
error.bar(1:Ncond,colMeans(plot_ter),colSds(plot_ter,na.rm=T)/sqrt(N),lwd=lwdmean,length=0)

##Bound
plot_bound <- with(param_1,aggregate(bound,by=list(sub=sub,condition=condition),mean))
plot_bound <- cast(plot_bound,sub~condition)
plot_bound <- plot_bound[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_bound),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_bound),max(plot_bound)),xaxt='n',family="A");
title(ylab="Bound",xlab="Feedback condition", line = linelab, cex.lab = cexlab,family="A")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=cexax,family="A")
for(i in 1:N) lines(jitter(1:Ncond,jit_size),plot_bound[i,1:Ncond],type='b',family="A",lty=2,col=col_indiv_points,pch=19)
points(colMeans(plot_bound),type='b',lwd=lwdmean)
error.bar(1:Ncond,colMeans(plot_bound),colSds(plot_bound,na.rm=T)/sqrt(N),lwd=lwdmean,length=0)

##Drift interaction
plot_drift_minus <- with(subset(param_1,difflevel=="hard"),
                         aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_minus <- cast(plot_drift_minus,sub~condition)
plot_drift_minus <- plot_drift_minus[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_drift_minus),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,xlim=c(.8,Ncond+.2),
     ylim=c(min(plot_drift_minus),.31),xaxt='n',family="A");
title(ylab="Drift rate",xlab="Feedback condition",family="A", line = linelab, cex.lab = cexlab)
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=cexax,family="A")
points(colMeans(plot_drift_minus),type='b',lwd=lwdmean,col=col_hard)
error.bar(1:Ncond,colMeans(plot_drift_minus),
          colSds(plot_drift_minus,na.rm=T)/sqrt(N),lwd=lwdmean,length=0,col=col_hard)

plot_drift_control <- with(subset(param_1,difflevel=="average"),
                           aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_control <- cast(plot_drift_control,sub~condition)
plot_drift_control <- plot_drift_control[,cond_order] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_control),type='b',lwd=lwdmean,col=col_med)
error.bar(1:Ncond,colMeans(plot_drift_control),
          colSds(plot_drift_control,na.rm=T)/sqrt(N),lwd=lwdmean,length=0,col=col_med)

plot_drift_plus <- with(subset(param_1,difflevel=="easy"),
                        aggregate(drift,by=list(sub=sub,condition=condition),mean))
plot_drift_plus <- cast(plot_drift_plus,sub~condition)
plot_drift_plus <- plot_drift_plus[,cond_order] #Reorder columns to have easy -> hard
points(colMeans(plot_drift_plus),type='b',lwd=lwdmean,col=col_easy)
error.bar(1:Ncond,colMeans(plot_drift_plus),
          colSds(plot_drift_plus,na.rm=T)/sqrt(N),lwd=lwdmean,length=0,col=col_easy)
legend("top",border=F,legend=c("Hard","Medium","Easy"),lwd=1,horiz=T,text.width = leg_text_width,
       x.intersp = leg_squeeze_within,col=c(col_hard,col_med,col_easy),bty="n",xjust=0,
       cex=cexleg,title = "Trial Difficulty",seg.len = 1)
# Plot alpha/beta ---------------------------------------------------------
##Alpha
plot_alpha <- with(param_1,aggregate(alpha,by=list(sub=sub,condition=condition),mean))
plot_alpha <- cast(plot_alpha,sub~condition)
plot_alpha <- plot_alpha[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_alpha),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_alpha),max(plot_alpha)),xaxt='n',family="A");
title(ylab="Alpha",xlab="Feedback condition", line = linelab, cex.lab = cexlab,family="A")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=cexax,family="A")
for(i in 1:N) lines(jitter(1:Ncond,jit_size),plot_alpha[i,1:Ncond],type='b',lty=2,col=col_indiv_points,pch=19)
points(colMeans(plot_alpha),type='b',lwd=lwdmean)
error.bar(1:Ncond,colMeans(plot_alpha),colSds(plot_alpha,na.rm=T)/sqrt(N),lwd=lwdmean,length=0)

##Beta
plot_beta <- with(param_1,aggregate(beta,by=list(sub=sub,condition=condition),mean))
plot_beta <- cast(plot_beta,sub~condition)
plot_beta <- plot_beta[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_beta),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_beta),max(plot_beta)),xaxt='n',family="A");
title(ylab="Beta",xlab="Feedback condition", line = linelab, cex.lab = cexlab,family="A")
axis(1,1:Ncond,c("Negative","Average","Positive"),cex.axis=cexax,family="A")
for(i in 1:N) lines(jitter(1:Ncond,jit_size),plot_beta[i,1:Ncond],type='b',lty=2,col=col_indiv_points,pch=19)
points(colMeans(plot_beta),type='b',lwd=lwdmean)
error.bar(1:Ncond,colMeans(plot_beta),colSds(plot_beta,na.rm=T)/sqrt(N),lwd=lwdmean,length=0)

dev.off()
# Confidence prediction separately for corrects and errors ----------------
#' Function to get the right font size in points
cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}

# Feedback condition colors
# Margins
mar.fb <- c(4,3.5,0,2)+0.1
mar.cj <- c(4,3.5,0,2)+0.1
mar.contrast <- c(5,6,0,2)+0.1
mar.parameters <- c(3.5,3.5,0,0)
mar.parameters.overall <- c(3.5,3.5,0,2)
### Adjust sizes and positions
cex.layout <- .66
# Lines and dots
lwd.dat <- 1.5
cex.datdot <- 1.125
cex.legend.square <- 3
# Text
cex.leg <- cex_size(8,cex.layout)
cex.phase <- cex_size(11,cex.layout)*cex.layout
cex.main <- cex_size(12,cex.layout)
cex.axis <- cex_size(8,cex.layout) 
cex.lab <- cex_size(10,cex.layout)*cex.layout
cex.lab.rel <- cex_size(10,cex.layout) # In a big layout, cex for mtext and cex.lab differ
cex.trialtype <- cex_size(8,cex.layout)
adj_ydens <- -.05; 
y_coord <- c(.5,1,1.5)
ylab_orientation <- 1 # 1 for horizontal, 0 for vertical
# Parameter values for feedback generation
fb_par_col <- 'black'
fb_par_pch <- 4

N_temp <- length(unique(Data1$sub))


## Look at the proportion of trials behind each data point
trial_count <- with(Data1,aggregate(resp,list(selfconf,cor,coh,sub),length))
names(trial_count) <- c("selfconf","cor","coh","sub","x")
trial_count <- cast(trial_count,sub+selfconf~cor+coh,value="x")
# Replace NAs with 0s
for (i in 1:ncol(trial_count)) {
  trial_count[is.na(trial_count[,i]),i] <- 0
}
colMeans(trial_count[,3:ncol(trial_count)]) # Average number of trial per cell

trial_count$sum <- rowSums(trial_count[,c(-1,-2)]) #Number of trials per subject and condition
# Now divide each cell by the sum
for (i in 3:(ncol(trial_count)-1)) {
  trial_count[,i] <- trial_count[,i]/trial_count$sum
}
colMeans(trial_count[,3:(ncol(trial_count)-1)]) # Average proportion of trial per cell

# Confidence ~ RT ---------------------------------------------------------
par(mfrow=c(1,1), mar = c(5,4,4,2)+.1)
cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
cex.layout <- .83 
ratio <- .66/cex.layout # To get comparable sizes with figure 4
cex.lab <- cex_size(10,cex.layout)*cex.layout # Set with mtext so ignore mfrow
cex.ax <- cex_size(8,cex.layout)
cex.leg <- cex_size(8,cex.layout)
cex.datdot <- cex_size(8,cex.layout);
lwd.dat <- 1.5*ratio
mar.plot <- c(4,4,2,1)+.1

Data1$rt_bin <- cut(Data1$rt,breaks=c(seq(0,1.5,.5),max(Data1$rt)))
test <- with(Data1,aggregate(cj,by=list(sub,rt_bin),mean))
names(test) <- c("sub","rt_bin","cj")
test <- cast(test,sub~rt_bin)
test <- test[,-1]


Simuls1$rt_bin <- cut(Simuls1$rt,breaks=c(seq(0,1.5,.5),max(Simuls1$rt)))
test_sim <- with(Simuls1,aggregate(cj,by=list(sub,rt_bin),mean))
names(test_sim) <- c("sub","rt_bin","cj")
test_sim <- cast(test_sim,sub~rt_bin)
test_sim <- test_sim[,-1]


n <- length(test)

na_count_dat <- sapply(test, function(y) sum(length(which(!is.na(y)))))
na_count_sim <- sapply(test_sim, function(y) sum(length(which(!is.na(y)))))

cex.datdot <- 1

tiff("cj_rt_exp1.tiff",width=8,height=8,units="cm",res=600,compression="lzw")
par(mar=mar.plot)

stripchart(test,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(3,6), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = "Experiment 1")
mtext("Confidence",2,at=4.5,line=2.5,cex=cex.lab);
mtext("RT bin (s)",1,at=n/2,line=2.5,cex=cex.lab);
axis(1,at=0:(n-1),labels=c(names(test)[-4],expression(paste("(1.5,",infinity,"]"))), cex.axis=cex.ax);
axis(2, seq(3,6,.5), cex.axis=cex.ax, las = 2)
means <- sapply(test, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=1)
error.bar(0:(n-1),means,colSds(as.matrix(test),na.rm=T)/sqrt(na_count_dat),lwd=lwd.dat)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim,na.rm=T) + (colSds(as.matrix(test_sim),na.rm=T)/sqrt(na_count_sim)),
          (colMeans(test_sim,na.rm=T) - colSds(as.matrix(test_sim),na.rm=T)/sqrt(na_count_sim))[n:1]),
        border=F,col = rgb(0,0,0,.5))
# Add trial accuracy as a legend
legend("topright",c("Empirical data","Model predictions"),lty=c(1,NA),
       lwd=1.5*ratio,col=c("black",rgb(0,0,0,.5)),border=NA,pch=c(16,15),pt.cex=c(1,2),
       cex=cex.leg,bty="n")
dev.off()

## Correlation between confidence and RT
# Behavior
correlation_data <- c()
for (subj in unique(Data1$sub)) {
  correlation_data <- c(correlation_data,cor(Data1[Data1$sub==subj,]$cj,Data1[Data1$sub==subj,]$rt,method = "spearman"))
}
t.test(correlation_data)

# Model predictions
correlation_sim <- c()
for (subj in unique(Simuls1$sub)) {
  correlation_sim <- c(correlation_sim,cor(Simuls1[Simuls1$sub==subj,]$cj,Simuls1[Simuls1$sub==subj,]$rt,method = "spearman"))
}
t.test(correlation_sim)

