rm(list=ls())
library(reshape);library(effects);library(lmerTest);library(scales);library(prob)
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir) ## change our current working directory

cexkl <- 1.5;cexgr <- 2;lwdgr <- 3;
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Experiment 1 ------------------------------------------------------------
for(i in 1:50){
  if(i == 1){
    Data <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1A/selfconfidence1A_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data <- rbind(Data,temp)
  }
}

head(Data)

Data <- subset(Data, running == "main")

Data['response'] <- 0
Data$response[Data$resp == "['n']"] <- 1

Data <- subset(Data,rt>.2) #There was no trial below 200ms

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data$sub)); subs <- unique(Data$sub); exclusion <- c()
tasks <- unique(Data$task)
par(mfrow=c(2,2))
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data,sub==subs[i]&task==t)
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

Data <- subset(Data,!(sub %in% exclusion)) #Sub 10 and 49 removed

Data$response[Data$response==0] <- -1

df <- Data[,c("sub","task","selfconf","difflevel","rt","response","cor","cj","RTconf","block")]
names(df) <- c("sub","task","selfconf","coh","rt","resp","cor","cj","RTconf","block")
df$rt <- df$rt/1000;df$RTconf <- df$RTconf/1000
# write.csv(df,"dataexp1_helene.csv",row.names = FALSE)

## Per block
for (b in unique(Data$block)) {
  df <- Data[Data$block==b,c("sub","task","selfconf","difflevel","rt","response","cor","cj","RTconf")]
  names(df) <- c("sub","task","selfconf","coh","rt","resp","cor","cj","RTconf")
  df$rt <- df$rt/1000;df$RTconf <- df$RTconf/1000
  write.csv(df,paste0("dataexp1_helene_block_",b,".csv"),row.names = FALSE)
}



# # Behaviour plots: Task ---------------------------------------------------------
# 
# ##Accuracy on task
# dot_colour <- with(subset(df,task=="dotcolour"),aggregate(cor,by=list(coh,sub),mean))
# names(dot_colour) <- c('coh','sub','cor')
# dot_colour <- cast(dot_colour,sub~coh)
# 
# dot_number <- with(subset(df,task=="dotnumber"),aggregate(cor,by=list(coh,sub),mean))
# names(dot_number) <- c('coh','sub','cor')
# dot_number <- cast(dot_number,sub~coh)
# 
# xo <- with(subset(df,task=="xotask"),aggregate(cor,by=list(coh,sub),mean))
# names(xo) <- c('coh','sub','cor')
# xo <- cast(xo,sub~coh)
# 
# 
# par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
# x <- dot_colour[,2:4]; x <- x[,c(2,1,3)]
# xlow <- dot_number[,2:4]; xlow <- xlow[,c(2,1,3)]
# xhigh <- xo[,2:4]; xhigh <- xhigh[,c(2,1,3)]
# 
# means <- sapply(x, mean);n<- length(x)
# stripchart(x,xlim = c(.9,n),ylim = c(.5,1) ,vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# # ,main="Coherence and Feedback effect on Confidence"
# mtext("Accuracy",2,at=.75,line=3,cex=1.75)
# axis(1,at=1:n,labels=names(x),cex.axis=1.75)
# for(i in seq(.5,1,length.out = 6)) abline(h=i,col='lightgrey')
# mtext("Trial difficulty",1,2.5,cex=1.75)
# 
# 
# for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
# lines(1:n,means,type='b',pch=2,col="orange")
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2,col="blue")
# for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0,col="red")
# for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
# legend("left",legend=c("dot colour","dot number","xo task"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5)
# 
# ##RT on task
# dot_colour <- with(subset(df,task=="dotcolour"),aggregate(rt,by=list(coh,sub),mean))
# names(dot_colour) <- c('coh','sub','rt')
# dot_colour <- cast(dot_colour,sub~coh)
# 
# dot_number <- with(subset(df,task=="dotnumber"),aggregate(rt,by=list(coh,sub),mean))
# names(dot_number) <- c('coh','sub','rt')
# dot_number <- cast(dot_number,sub~coh)
# 
# xo <- with(subset(df,task=="xotask"),aggregate(rt,by=list(coh,sub),mean))
# names(xo) <- c('coh','sub','rt')
# xo <- cast(xo,sub~coh)
# 
# 
# par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
# x <- dot_colour[,2:4]; x <- x[,c(2,1,3)]
# xlow <- dot_number[,2:4]; xlow <- xlow[,c(2,1,3)]
# xhigh <- xo[,2:4]; xhigh <- xhigh[,c(2,1,3)]
# 
# means <- sapply(x, mean);n<- length(x)
# stripchart(x,xlim = c(.9,n) ,ylim=c(.5,1.5),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# # ,main="Coherence and Feedback effect on Confidence"
# mtext("Reaction time",2,at=1,line=3,cex=1.75)
# axis(1,at=1:n,labels=names(x),cex.axis=1.75)
# for(i in seq(.6,1.4,length.out = 5)) abline(h=i,col='lightgrey')
# mtext("Trial difficulty",1,2.5,cex=1.75)
# 
# 
# for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
# lines(1:n,means,type='b',pch=2,col="orange")
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2,col="blue")
# for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0,col="red")
# for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
# legend("topleft",legend=c("dot colour","dot number","xo task"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5)
# 
# ##Confidence on task
# dot_colour <- with(subset(df,task=="dotcolour"),aggregate(cj,by=list(coh,sub),mean))
# names(dot_colour) <- c('coh','sub','cj')
# dot_colour <- cast(dot_colour,sub~coh)
# 
# dot_number <- with(subset(df,task=="dotnumber"),aggregate(cj,by=list(coh,sub),mean))
# names(dot_number) <- c('coh','sub','cj')
# dot_number <- cast(dot_number,sub~coh)
# 
# xo <- with(subset(df,task=="xotask"),aggregate(cj,by=list(coh,sub),mean))
# names(xo) <- c('coh','sub','cj')
# xo <- cast(xo,sub~coh)
# 
# 
# par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
# x <- dot_colour[,2:4]; x <- x[,c(2,1,3)]
# xlow <- dot_number[,2:4]; xlow <- xlow[,c(2,1,3)]
# xhigh <- xo[,2:4]; xhigh <- xhigh[,c(2,1,3)]
# 
# means <- sapply(x, mean);n<- length(x)
# stripchart(x,xlim = c(.9,n) ,ylim=c(3,6),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# # ,main="Coherence and Feedback effect on Confidence"
# mtext("Confidence",2,at=4.5,line=3,cex=1.75)
# axis(1,at=1:n,labels=names(x),cex.axis=1.75)
# for(i in seq(3,6,length.out = 7)) abline(h=i,col='lightgrey')
# mtext("Trial difficulty",1,2.5,cex=1.75)
# 
# 
# for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
# lines(1:n,means,type='b',pch=2,col="orange")
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2,col="blue")
# for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0,col="red")
# for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
# legend("topright",legend=c("dot colour","dot number","xo task"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5)
# 
# 
# 
# # Behaviour plots: Condition ----------------------------------------------
# ##Accuracy on task
# traineasy <- with(subset(df,selfconf=="highSC"),aggregate(cor,by=list(coh,sub),mean))
# names(traineasy) <- c('coh','sub','cor')
# traineasy <- cast(traineasy,sub~coh)
# 
# trainaverage <- with(subset(df,selfconf=="mediumSC"),aggregate(cor,by=list(coh,sub),mean))
# names(trainaverage) <- c('coh','sub','cor')
# trainaverage <- cast(trainaverage,sub~coh)
# 
# trainhard <- with(subset(df,selfconf=="lowSC"),aggregate(cor,by=list(coh,sub),mean))
# names(trainhard) <- c('coh','sub','cor')
# trainhard <- cast(trainhard,sub~coh)
# 
# 
# par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
# x <- traineasy[,2:4]; x <- x[,c(2,1,3)]
# xlow <- trainaverage[,2:4]; xlow <- xlow[,c(2,1,3)]
# xhigh <- trainhard[,2:4]; xhigh <- xhigh[,c(2,1,3)]
# 
# means <- sapply(x, mean);n<- length(x)
# stripchart(x,xlim = c(.9,n),ylim = c(.5,1) ,vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# # ,main="Coherence and Feedback effect on Confidence"
# mtext("Accuracy",2,at=.75,line=3,cex=1.75)
# axis(1,at=1:n,labels=names(x),cex.axis=1.75)
# for(i in seq(.5,1,length.out = 6)) abline(h=i,col='lightgrey')
# mtext("Trial difficulty",1,2.5,cex=1.75)
# 
# 
# for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
# lines(1:n,means,type='b',pch=2,col="orange")
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2,col="blue")
# for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0,col="red")
# for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
# legend("left",legend=c("Positive","Average","Negative"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5,title = "Fake feedback condition")
# 
# ##RT on task
# dfcor <- subset(df,cor==1)
# traineasy <- with(subset(dfcor,selfconf=="highSC"),aggregate(rt,by=list(coh,sub),mean))
# names(traineasy) <- c('coh','sub','rt')
# traineasy <- cast(traineasy,sub~coh)
# 
# trainaverage <- with(subset(dfcor,selfconf=="mediumSC"),aggregate(rt,by=list(coh,sub),mean))
# names(trainaverage) <- c('coh','sub','rt')
# trainaverage <- cast(trainaverage,sub~coh)
# 
# trainhard <- with(subset(dfcor,selfconf=="lowSC"),aggregate(rt,by=list(coh,sub),mean))
# names(trainhard) <- c('coh','sub','rt')
# trainhard <- cast(trainhard,sub~coh)
# 
# 
# par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
# x <- traineasy[,2:4]; x <- x[,c(2,1,3)]
# xlow <- trainaverage[,2:4]; xlow <- xlow[,c(2,1,3)]
# xhigh <- trainhard[,2:4]; xhigh <- xhigh[,c(2,1,3)]
# 
# means <- sapply(x, mean);n<- length(x)
# stripchart(x,xlim = c(.9,n) ,ylim=c(.5,1.5),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# # ,main="Coherence and Feedback effect on Confidence"
# mtext("Reaction time",2,at=1,line=3,cex=1.75)
# axis(1,at=1:n,labels=names(x),cex.axis=1.75)
# for(i in seq(.6,1.4,length.out = 5)) abline(h=i,col='lightgrey')
# mtext("Trial difficulty",1,2.5,cex=1.75)
# 
# 
# for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
# lines(1:n,means,type='b',pch=2,col="orange")
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2,col="blue")
# for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0,col="red")
# for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
# legend("topright",legend=c("Positive","Average","Negative"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5,title = "Fake feedback condition")
# 
# 
# ##Confidence on task
# traineasy <- with(subset(df,selfconf=="highSC"),aggregate(cj,by=list(coh,sub),mean))
# names(traineasy) <- c('coh','sub','cj')
# traineasy <- cast(traineasy,sub~coh)
# 
# trainaverage <- with(subset(df,selfconf=="mediumSC"),aggregate(cj,by=list(coh,sub),mean))
# names(trainaverage) <- c('coh','sub','cj')
# trainaverage <- cast(trainaverage,sub~coh)
# 
# trainhard <- with(subset(df,selfconf=="lowSC"),aggregate(cj,by=list(coh,sub),mean))
# names(trainhard) <- c('coh','sub','cj')
# trainhard <- cast(trainhard,sub~coh)
# 
# 
# par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
# x <- traineasy[,2:4]; x <- x[,c(2,1,3)]
# xlow <- trainaverage[,2:4]; xlow <- xlow[,c(2,1,3)]
# xhigh <- trainhard[,2:4]; xhigh <- xhigh[,c(2,1,3)]
# 
# means <- sapply(x, mean);n<- length(x)
# stripchart(x,xlim = c(.9,n) ,ylim=c(3,6),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# # ,main="Coherence and Feedback effect on Confidence"
# mtext("Confidence",2,at=4.5,line=3,cex=1.75)
# axis(1,at=1:n,labels=names(x),cex.axis=1.75)
# for(i in seq(3,6,length.out = 7)) abline(h=i,col='lightgrey')
# mtext("Trial difficulty",1,2.5,cex=1.75)
# 
# 
# for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
# lines(1:n,means,type='b',pch=2,col="orange")
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2,col="blue")
# for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0,col="red")
# for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
# error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
# legend("topright",legend=c("Positive","Average","Negative"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5,title = "Fake feedback condition")
# # Mixed Models ------------------------------------------------------------
# 
# library(car); library(ARTool); library(emmeans); library(lme4); library(multcomp);library(lattice)
# 
# df$coh <- as.factor(df$coh)
# df$selfconf <- as.factor(df$selfconf)
# df$task <- as.factor(df$task)
# # fit <- lmer(cj ~ selfconf * task * coh + (1|sub), data = df,REML = F)
# # fit_a <- lmer(cj ~ selfconf * task * coh + (1+selfconf|sub), data = df,REML = F);anova(fit,fit_a)
# # fit_b <- lmer(cj ~ selfconf * task * coh + (1+task|sub), data = df,REML = F);anova(fit,fit_b)
# # fit_c <- lmer(cj ~ selfconf * task * coh + (1+coh|sub), data = df,REML = F);anova(fit,fit_c)
# # m <- lmer(cj ~ selfconf * task * coh + (selfconf|sub) + (coh|sub) + (task|sub), data = df,REML = F) #Failed to converge
# 
# m <- lmer(cj ~ selfconf * task * coh + (1 + selfconf|sub), data = df,REML = T)
# anova(m)
# # emmeans(m, list(pairwise ~ coh), adjust = "tukey")
# # emmeans(m, list(pairwise ~ selfconf), adjust = "tukey")
# 
# plot(resid(m),df$cj)
# plot(m) #Linearity + Homogeneity check
# qqmath(m) #Normality
# 
# ##RTs
# dfcor <- subset(df,cor==1)
# # fit <- lmer(rt ~ selfconf * task * coh + (1|sub), data = dfcor,REML = F)
# # fit_a <- lmer(rt ~ selfconf * task * coh + (1+selfconf|sub), data = dfcor,REML = F);anova(fit,fit_a)
# # fit_b <- lmer(rt ~ selfconf * task * coh + (1+task|sub), data = dfcor,REML = F);anova(fit,fit_b) #Failed to converge
# # fit_c <- lmer(rt ~ selfconf * task * coh + (1+coh|sub), data = dfcor,REML = F);anova(fit,fit_c) #Failed to converge
# 
# 
# m <- lmer(rt ~ selfconf * task * coh + (1 + selfconf|sub), data = dfcor,REML = T) #Failed to converge
# anova(m)
# # emmeans(m, list(pairwise ~ coh), adjust = "tukey")
# # emmeans(m, list(pairwise ~ selfconf), adjust = "tukey")
# 
# plot(m) #Linearity + Homogeneity check
# qqmath(m) #Normality
# 
# m <- glmer(cor ~ selfconf * task * coh + (1 + selfconf|sub), data = df,family = "binomial",nAGQ = 0) 
# Anova(m)
# Experiment 2 ------------------------------------------------------------

for(i in 1:50){ #123=tinne
  if(i == 1){ 
    Data <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
  }else{
    temp <- read.csv(paste0('RealData_1B/selfconfidence1B_sub',i,'.csv'),fileEncoding="UTF-8-BOM")
    Data <- rbind(Data,temp)
  }
}

head(Data)

Data <- subset(Data, running == "main")

Data['response'] <- 0
Data$response[Data$resp == "['n']"] <- 1

Data <- subset(Data,rt>200) #There was no trial below 200ms

## Diagnostic plot per participant and task + chance performance testing
N <- length(unique(Data$sub)); subs <- unique(Data$sub); exclusion <- c()
tasks <- unique(Data$task)
par(mfrow=c(2,2))
# at_chance <- matrix(NA,nrow=N,ncol=2)
for(i in 1:N){
  for (t in tasks) {
    tempDat <- subset(Data,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    # at_chance[i,] <- c(subs[i], t.test(acc_block$x,mu=.5)$p.value)
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt/1000,frame=F,col=c("black"),main=paste('subject',i,"task :",t),ylab="RT",ylim=c(0,5))
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
  }
}

Data <- subset(Data,!(sub %in% exclusion)) #Sub 12, 32 and 50 removed

Data$response[Data$response==0] <- -1

df <- Data[,c("sub","task","traindiffcond","trialdifflevel","rt","response","cor","cj","RTconf","block")]
names(df) <- c("sub","task","traindiffcond","coh","rt","resp","cor","cj","RTconf","block")
df$rt <- df$rt/1000;df$RTconf <- df$RTconf/1000
# write.csv(df,"dataexp2_helene.csv",row.names = FALSE)

# Per block
for (b in unique(Data$block)) {
  df <- Data[Data$block==b,c("sub","task","traindiffcond","trialdifflevel","rt","response","cor","cj","RTconf")]
  names(df) <- c("sub","task","traindiffcond","coh","rt","resp","cor","cj","RTconf")
  df$rt <- df$rt/1000;df$RTconf <- df$RTconf/1000
  write.csv(df,paste0("dataexp2_helene_block_",b,".csv"),row.names = FALSE)
}

# Mixed Models ------------------------------------------------------------
library(car); library(ARTool); library(emmeans); library(lme4); library(multcomp);library(lattice)

df$coh <- as.factor(df$coh)
df$traindiffcond <- as.factor(df$traindiffcond)
df$task <- as.factor(df$task)
fit <- lmer(cj ~ traindiffcond * task * coh + (1|sub), data = df,REML = F)
fit_a <- lmer(cj ~ traindiffcond * task * coh + (1+traindiffcond|sub), data = df,REML = F);anova(fit,fit_a)
fit_b <- lmer(cj ~ traindiffcond * task * coh + (1+task|sub), data = df,REML = F);anova(fit,fit_b)
fit_c <- lmer(cj ~ traindiffcond * task * coh + (1+coh|sub), data = df,REML = F);anova(fit,fit_c)
m <- lmer(cj ~ traindiffcond * task * coh + (1+traindiffcond+coh+task|sub), data = df,REML = T) #Failed to converge

step_res <- step(m, reduce.fixed=F)

m <- lmer(cj ~ traindiffcond * task * coh + (1 + traindiffcond|sub), data = df,REML = T)
anova(m)
emmeans(m, list(pairwise ~ coh), adjust = "tukey")
emmeans(m, list(pairwise ~ traindiffcond), adjust = "tukey")

plot(resid(m),df$cj)
plot(m) #Linearity + Homogeneity check
qqmath(m) #Normality

##RTs
dfcor <- subset(df,cor==1)
fit <- lmer(rt ~ traindiffcond * task * coh + (1|sub), data = dfcor,REML = F)
fit_a <- lmer(rt ~ traindiffcond * task * coh + (1+traindiffcond|sub), data = dfcor,REML = F);anova(fit,fit_a)
fit_b <- lmer(rt ~ traindiffcond * task * coh + (1+task|sub), data = dfcor,REML = F);anova(fit,fit_b) #Failed to converge
fit_c <- lmer(rt ~ traindiffcond * task * coh + (1+coh|sub), data = dfcor,REML = F);anova(fit,fit_c) #Failed to converge


m <- lmer(rt ~ traindiffcond * task * coh + (1 + traindiffcond|sub), data = dfcor,REML = T) #Failed to converge
anova(m)
emmeans(m, list(pairwise ~ coh), adjust = "tukey")
emmeans(m, list(pairwise ~ traindiffcond), adjust = "tukey")

plot(m) #Linearity + Homogeneity check
qqmath(m) #Normality

##Accuracy
fit <- glmer(cor ~ traindiffcond * task * coh + (1|sub), data = df, family = "binomial",nAGQ = 0)
Anova(fit)


# Behaviour plots: Task ---------------------------------------------------------

##Accuracy on task
dot_colour <- with(subset(df,task=="dotcolour"),aggregate(cor,by=list(coh,sub),mean))
names(dot_colour) <- c('coh','sub','cor')
dot_colour <- cast(dot_colour,sub~coh)

dot_number <- with(subset(df,task=="dotnumber"),aggregate(cor,by=list(coh,sub),mean))
names(dot_number) <- c('coh','sub','cor')
dot_number <- cast(dot_number,sub~coh)

xo <- with(subset(df,task=="xotask"),aggregate(cor,by=list(coh,sub),mean))
names(xo) <- c('coh','sub','cor')
xo <- cast(xo,sub~coh)


par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- dot_colour[,2:4]; x <- x[,c(2,1,3)]
xlow <- dot_number[,2:4]; xlow <- xlow[,c(2,1,3)]
xhigh <- xo[,2:4]; xhigh <- xhigh[,c(2,1,3)]

means <- sapply(x, mean);n<- length(x)
stripchart(x,xlim = c(.9,n),ylim = c(.5,1) ,vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# ,main="Coherence and Feedback effect on Confidence"
mtext("Accuracy",2,at=.75,line=3,cex=1.75)
axis(1,at=1:n,labels=names(x),cex.axis=1.75)
for(i in seq(.5,1,length.out = 6)) abline(h=i,col='lightgrey')
mtext("Trial difficulty",1,2.5,cex=1.75)


for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
lines(1:n,means,type='b',pch=2,col="orange")
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xlow, mean,na.rm=T);
lines(1:n,means,type='b',pch=2,col="blue")
for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
means <- sapply(xhigh, mean,na.rm=T);
lines(1:n,means,type='b',pch=0,col="red")
for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
legend("left",legend=c("dot colour","dot number","xo task"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5)

##RT on task
dot_colour <- with(subset(df,task=="dotcolour"),aggregate(rt,by=list(coh,sub),mean))
names(dot_colour) <- c('coh','sub','rt')
dot_colour <- cast(dot_colour,sub~coh)

dot_number <- with(subset(df,task=="dotnumber"),aggregate(rt,by=list(coh,sub),mean))
names(dot_number) <- c('coh','sub','rt')
dot_number <- cast(dot_number,sub~coh)

xo <- with(subset(df,task=="xotask"),aggregate(rt,by=list(coh,sub),mean))
names(xo) <- c('coh','sub','rt')
xo <- cast(xo,sub~coh)


par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- dot_colour[,2:4]; x <- x[,c(2,1,3)]
xlow <- dot_number[,2:4]; xlow <- xlow[,c(2,1,3)]
xhigh <- xo[,2:4]; xhigh <- xhigh[,c(2,1,3)]

means <- sapply(x, mean);n<- length(x)
stripchart(x,xlim = c(.9,n) ,ylim=c(.5,1.5),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# ,main="Coherence and Feedback effect on Confidence"
mtext("Reaction time",2,at=1,line=3,cex=1.75)
axis(1,at=1:n,labels=names(x),cex.axis=1.75)
for(i in seq(.6,1.4,length.out = 5)) abline(h=i,col='lightgrey')
mtext("Trial difficulty",1,2.5,cex=1.75)


for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
lines(1:n,means,type='b',pch=2,col="orange")
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xlow, mean,na.rm=T);
lines(1:n,means,type='b',pch=2,col="blue")
for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
means <- sapply(xhigh, mean,na.rm=T);
lines(1:n,means,type='b',pch=0,col="red")
for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
legend("topleft",legend=c("dot colour","dot number","xo task"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5)

##Confidence on task
dot_colour <- with(subset(df,task=="dotcolour"),aggregate(cj,by=list(coh,sub),mean))
names(dot_colour) <- c('coh','sub','cj')
dot_colour <- cast(dot_colour,sub~coh)

dot_number <- with(subset(df,task=="dotnumber"),aggregate(cj,by=list(coh,sub),mean))
names(dot_number) <- c('coh','sub','cj')
dot_number <- cast(dot_number,sub~coh)

xo <- with(subset(df,task=="xotask"),aggregate(cj,by=list(coh,sub),mean))
names(xo) <- c('coh','sub','cj')
xo <- cast(xo,sub~coh)


par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- dot_colour[,2:4]; x <- x[,c(2,1,3)]
xlow <- dot_number[,2:4]; xlow <- xlow[,c(2,1,3)]
xhigh <- xo[,2:4]; xhigh <- xhigh[,c(2,1,3)]

means <- sapply(x, mean);n<- length(x)
stripchart(x,xlim = c(.9,n) ,ylim=c(3,6),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# ,main="Coherence and Feedback effect on Confidence"
mtext("Confidence",2,at=4.5,line=3,cex=1.75)
axis(1,at=1:n,labels=names(x),cex.axis=1.75)
for(i in seq(3,6,length.out = 7)) abline(h=i,col='lightgrey')
mtext("Trial difficulty",1,2.5,cex=1.75)


for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
lines(1:n,means,type='b',pch=2,col="orange")
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xlow, mean,na.rm=T);
lines(1:n,means,type='b',pch=2,col="blue")
for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
means <- sapply(xhigh, mean,na.rm=T);
lines(1:n,means,type='b',pch=0,col="red")
for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
legend("topright",legend=c("dot colour","dot number","xo task"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5)



# Behaviour plots: Condition ----------------------------------------------
##Accuracy on task
traineasy <- with(subset(df,traindiffcond=="easy"),aggregate(cor,by=list(coh,sub),mean))
names(traineasy) <- c('coh','sub','cor')
traineasy <- cast(traineasy,sub~coh)

trainaverage <- with(subset(df,traindiffcond=="average"),aggregate(cor,by=list(coh,sub),mean))
names(trainaverage) <- c('coh','sub','cor')
trainaverage <- cast(trainaverage,sub~coh)

trainhard <- with(subset(df,traindiffcond=="hard"),aggregate(cor,by=list(coh,sub),mean))
names(trainhard) <- c('coh','sub','cor')
trainhard <- cast(trainhard,sub~coh)


par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- traineasy[,2:4]; x <- x[,c(2,1,3)]
xlow <- trainaverage[,2:4]; xlow <- xlow[,c(2,1,3)]
xhigh <- trainhard[,2:4]; xhigh <- xhigh[,c(2,1,3)]

means <- sapply(x, mean);n<- length(x)
stripchart(x,xlim = c(.9,n),ylim = c(.5,1) ,vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# ,main="Coherence and Feedback effect on Confidence"
mtext("Accuracy",2,at=.75,line=3,cex=1.75)
axis(1,at=1:n,labels=names(x),cex.axis=1.75)
for(i in seq(.5,1,length.out = 6)) abline(h=i,col='lightgrey')
mtext("Trial difficulty",1,2.5,cex=1.75)


for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
lines(1:n,means,type='b',pch=2,col="orange")
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xlow, mean,na.rm=T);
lines(1:n,means,type='b',pch=2,col="blue")
for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
means <- sapply(xhigh, mean,na.rm=T);
lines(1:n,means,type='b',pch=0,col="red")
for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
legend("left",legend=c("Easy","Average","Hard"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5,title = "Training difficulty")

##RT on task
dfcor <- subset(df,cor==1)
traineasy <- with(subset(dfcor,traindiffcond=="easy"),aggregate(rt,by=list(coh,sub),mean))
names(traineasy) <- c('coh','sub','rt')
traineasy <- cast(traineasy,sub~coh)

trainaverage <- with(subset(dfcor,traindiffcond=="average"),aggregate(rt,by=list(coh,sub),mean))
names(trainaverage) <- c('coh','sub','rt')
trainaverage <- cast(trainaverage,sub~coh)

trainhard <- with(subset(dfcor,traindiffcond=="hard"),aggregate(rt,by=list(coh,sub),mean))
names(trainhard) <- c('coh','sub','rt')
trainhard <- cast(trainhard,sub~coh)


par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- traineasy[,2:4]; x <- x[,c(2,1,3)]
xlow <- trainaverage[,2:4]; xlow <- xlow[,c(2,1,3)]
xhigh <- trainhard[,2:4]; xhigh <- xhigh[,c(2,1,3)]

means <- sapply(x, mean);n<- length(x)
stripchart(x,xlim = c(.9,n) ,ylim=c(.5,1.5),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# ,main="Coherence and Feedback effect on Confidence"
mtext("Reaction time",2,at=1,line=3,cex=1.75)
axis(1,at=1:n,labels=names(x),cex.axis=1.75)
for(i in seq(.6,1.4,length.out = 5)) abline(h=i,col='lightgrey')
mtext("Trial difficulty",1,2.5,cex=1.75)

# stripchart(x,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=alpha('orange',.5),col='white',cex=1.5)
# stripchart(xlow,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg=alpha('blue',.5),col='white')
# stripchart(xhigh,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg=alpha('red',.5),col='white')
for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
lines(1:n,means,type='b',pch=2,col="orange")
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xlow, mean,na.rm=T);
lines(1:n,means,type='b',pch=2,col="blue")
for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
means <- sapply(xhigh, mean,na.rm=T);
lines(1:n,means,type='b',pch=0,col="red")
for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
legend("topleft",legend=c("Easy","Average","Hard"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5,title = "Training difficulty")


##Confidence on task
traineasy <- with(subset(df,traindiffcond=="easy"),aggregate(cj,by=list(coh,sub),mean))
names(traineasy) <- c('coh','sub','cj')
traineasy <- cast(traineasy,sub~coh)

trainaverage <- with(subset(df,traindiffcond=="average"),aggregate(cj,by=list(coh,sub),mean))
names(trainaverage) <- c('coh','sub','cj')
trainaverage <- cast(trainaverage,sub~coh)

trainhard <- with(subset(df,traindiffcond=="hard"),aggregate(cj,by=list(coh,sub),mean))
names(trainhard) <- c('coh','sub','cj')
trainhard <- cast(trainhard,sub~coh)


par(mar=c(5.1, 5.1, 4.1, 4.1),mfrow=c(1,1))
x <- traineasy[,2:4]; x <- x[,c(2,1,3)]
xlow <- trainaverage[,2:4]; xlow <- xlow[,c(2,1,3)]
xhigh <- trainhard[,2:4]; xhigh <- xhigh[,c(2,1,3)]

means <- sapply(x, mean);n<- length(x)
stripchart(x,xlim = c(.9,n) ,ylim=c(3,6),vertical = TRUE, col="white",frame=F,xaxt='n',cex.axis=1.75)
# ,main="Coherence and Feedback effect on Confidence"
mtext("Confidence",2,at=4.5,line=3,cex=1.75)
axis(1,at=1:n,labels=names(x),cex.axis=1.75)
for(i in seq(3,6,length.out = 7)) abline(h=i,col='lightgrey')
mtext("Trial difficulty",1,2.5,cex=1.75)

# stripchart(x,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=alpha('orange',.5),col='white',cex=1.5)
# stripchart(xlow,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg=alpha('blue',.5),col='white')
# stripchart(xhigh,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg=alpha('red',.5),col='white')
for(i in 1:n) points(i,means[i],pch=21,bg="orange",col="white",cex=2)
lines(1:n,means,type='b',pch=2,col="orange")
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="orange")
means <- sapply(xlow, mean,na.rm=T);
lines(1:n,means,type='b',pch=2,col="blue")
for(i in 1:n) points(i,means[i],pch=24,bg="blue",col="white",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="blue")
means <- sapply(xhigh, mean,na.rm=T);
lines(1:n,means,type='b',pch=0,col="red")
for(i in 1:n) points(i,means[i],pch=15,bg="red",col="red",cex=1.5)
error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col="red")
legend("topright",legend=c("Easy","Average","Hard"),pch=c(16,17,15),lty=1,bty = "n", col = c("orange","blue","red"),cex=1.5,title = "Training difficulty")


# Trash -------------------------------------------------------------------

# stripchart(x,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=21,vertical = TRUE,bg=alpha('orange',.5),col='white',cex=1.5)
# stripchart(xlow,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg=alpha('blue',.5),col='white')
# stripchart(xhigh,at=1:n,jitter=.01, method = 'jitter',add=TRUE,pch=24,vertical = TRUE,bg=alpha('red',.5),col='white')
# for(i in 1:n) points(i,means[i],pch=21,bg="black",col="white",cex=2)
# lines(1:n,means,type='b',pch=2)
# # error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
# means <- sapply(xlow, mean,na.rm=T);
# lines(1:n,means,type='b',pch=2)
# # error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
# for(i in 1:n) points(i,means[i],pch=24,bg="black",col="white",cex=1.5)
# means <- sapply(xhigh, mean,na.rm=T);
# lines(1:n,means,type='b',pch=0)
# # error.bar(1:n,means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)
# for(i in 1:n) points(i,means[i],pch=15,bg="black",col="black",cex=1.5)
