rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage)
go_to("functions")
library(Rcpp) 
source('quantilefit_function_fullRTconf.R')
library(DEoptim)
library(MALDIquant)
library(reshape)
library(OneR)
library(scales)
library(lattice) # qqmath 
library(lmerTest)
library(fields) 
library(prob)
library(car)
library(multcomp)
setwd(curdir)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Data Load EXP 1 ---------------------------------------------------------------
go_to("data")
go_to("exp1")
# setwd(paste0(curdir,"/exp1"))
Data1 <- read.csv('dataexp1_prior_belief.csv')
Data1 <- subset(Data1, rt > .1 & rt < 4) # Same procedure as prior belief paper

condLab <- unique(Data1$selfconf); Ncond <- length(condLab) #Change column name according to condition label in dataset
subs <- unique(Data1$sub);N<-length(subs)

difficulty <- sort(unique(Data1$coh));
# Fits Load ---------------------------------------------------------------
go_to("fits")
go_to("exp1")

bound <- matrix(NA,N,Ncond);v <- matrix(NA,N,Ncond);ter <- matrix(NA,N,Ncond);conf_rt <- matrix(NA,N,Ncond);alpha <- matrix(NA,N,Ncond);beta <- matrix(NA,N,Ncond); resid <- matrix(NA,N,Ncond)
v2 <- matrix(NA,N,Ncond);v3 <- matrix(NA,N,Ncond); #Adjust the number of drift parameters to the model loaded
vratio <- matrix(NA,N,Ncond)

for(i in 1:N){
  for(c in 1:Ncond){
    print(paste('Running participant',i,'from',N,"condition",c))
    load(paste0('results_sub_',subs[i],'_',condLab[c],'.Rdata'))
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
##V-ratio
m <- lmer(vratio ~ condition + (1 |sub),data=param_1)
anova(m)
# Behavior analysis of prediction -------------------------------------------------
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
# Plot Layout -------------------------------------------------------------
go_to("plots")

ratio <- .66/.83 # To get comparable sizes with figure 4
cex.ax <- 1*ratio
cex.lab <- 1*ratio
cex.ax <- 1*ratio
cex.leg <- 1*ratio
cex.datdot <- 1*ratio;lwd.dat <- 1.5*ratio
mar.rt <- c(1,4,2,2)+.1
mar.acc <- c(5,4,0,2)+.1
mar.cj <- c(5,4,2,2)+.1
tiff(file = "exp1_results.tif",width = 16,height=8,units="cm",compression="lzw",res = 1200)
layout(matrix(c(1,2,3,3),ncol=2),widths = c(1,2))
par(mar=c(0,0,0,0))
par(mar=c(5,4,4,2)+.1)

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
           main="",cex.axis=cex.ax)
mtext("Response time (s)",2,at=.85,line=3,cex=cex.lab);
axis(1,at=0:(n-1),labels=c("","",""), cex.axis=cex.ax);
# Empirical data
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="brown3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="cyan4",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="darkgoldenrod3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="darkgoldenrod3")
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
           main="",cex.axis=cex.ax)
mtext("Accuracy",2,at=.8,line=3,cex=cex.lab);
axis(1,at=0:(n-1),labels=c("Hard","Average","Easy"), cex.axis=cex.ax);
mtext("Trial difficulty",1,3,at=1,cex=cex.lab)
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
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="cyan4",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="darkgoldenrod3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="darkgoldenrod3")
# CJ ----------------------------------------------------------------------

par(mar=mar.cj)

## Experiment 1
#Aggregate conf for data
CJlow <- with(subset(Data1,selfconf=="lowSC"),aggregate(cj,by=list(sub,coh),mean));names(CJlow) <- c('sub','coh','cj')
CJlow <- cast(CJlow,sub~coh,)
CJmed <- with(subset(Data1,selfconf=="mediumSC"),aggregate(cj,by=list(sub,coh),mean));names(CJmed) <- c('sub','coh','cj')
CJmed <- cast(CJmed,sub~coh)
CJhigh <- with(subset(Data1,selfconf=="highSC"),aggregate(cj,by=list(sub,coh),mean));names(CJhigh) <- c('sub','coh','cj')
CJhigh <- cast(CJhigh,sub~coh)

# Remove subject column
x <- CJlow[,c(2:4)];xmed <- CJmed[,c(2:4)];xhigh <- CJhigh[,c(2:4)]
n <- length(x)


#Aggregate conf for model
simDat_low <- subset(Simuls1,condition=="lowSC");simDat_low <- simDat_low[c('sub','coh','cj')]
simDat_med <- subset(Simuls1,condition=="mediumSC");simDat_med <- simDat_med[c('sub','coh','cj')]
simDat_high <- subset(Simuls1,condition=="highSC");simDat_high <- simDat_high[c('sub','coh','cj')]

snrCJlowSim <- aggregate(.~sub+coh,simDat_low,mean)
snrCJmedSim <- aggregate(.~sub+coh,simDat_med,mean)
snrCJhighSim <- aggregate(.~sub+coh,simDat_high,mean)

x_sim = cast(snrCJlowSim,sub~coh,value='cj');
xmed_sim = cast(snrCJmedSim,sub~coh,value='cj');
xhigh_sim = cast(snrCJhighSim,sub~coh,value='cj')

# Reorder columns from hard to easy difficulty
x <- x[,c("hard","average","easy")];
xmed <- xmed[,c("hard","average","easy")];
xhigh <- xhigh[,c("hard","average","easy")]
x_sim <- x_sim[,c("hard","average","easy")];
xmed_sim <- xmed_sim[,c("hard","average","easy")];
xhigh_sim <- xhigh_sim[,c("hard","average","easy")];

# Start plotting
stripchart(x, ylim=c(3.5,5.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=cex.ax)
mtext("Confidence",2,at=4.5,line=2.25,cex=cex.lab);
axis(1,at=0:(n-1),labels=c("Hard","Average","Easy"), cex.axis=cex.ax);
mtext("Trial difficulty",1,3,at=1,cex=cex.lab)
legend(.05,5.5,legend=c("Positive","Average","Negative"),lty=1,bty = "n",
       title = "Fake Feedback condition",pch=rep(16,3),inset=.1, cex = cex.leg,
       col=c("darkgoldenrod3","cyan4","brown3"))
legend(1,4,legend=c("Empirical data", "Model prediction"), pch=c(16,15),bty="n",
       col="grey",lty = c(1,NA), border = NA,cex=cex.leg,pt.cex = c(1,2) )

# Plot empirical data
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="brown3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="brown3")
means <- sapply(xmed, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="cyan4",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xmed),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="cyan4")
means <- sapply(xhigh, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col="darkgoldenrod3",lwd=lwd.dat)
error.bar(0:(n-1),means,colSds(as.matrix(xhigh),na.rm=T)/sqrt(N),lwd=lwd.dat,length=.05,col="darkgoldenrod3")

# Add model predictions
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[3:1]),
        border=F,col=rgb(205,51,51,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xmed_sim,na.rm=T) + (colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N)),(colMeans(xmed_sim,na.rm=T) - colSds(as.matrix(xmed_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(0,139,139,51,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xhigh_sim,na.rm=T) + (colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N)),(colMeans(xhigh_sim,na.rm=T) - colSds(as.matrix(xhigh_sim),na.rm=T)/sqrt(N))[3:1]),
        border=F,col=rgb(205,149,12,51,maxColorValue = 255))

dev.off()
