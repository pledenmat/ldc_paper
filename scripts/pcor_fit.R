curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)

library(myPackage)
library(reshape)
library(timeSeries)
go_to("functions")
# source("quantilefit_function_DDM.R")
source("rw_createHM_driftsign.r")

font <- "Arial"
windowsFonts(A = windowsFont(font))

#' Function to get the right font size in points
cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}
cex.layout <- 1
cexlab <- cex_size(10,cex.layout)*cex.layout
cexax <- cex_size(8,cex.layout)
cex.leg <- cex_size(8,cex.layout)

# We consider uniformly distributed drift rates for the reference heatmap
Ndrift <- 1000
v_max <- .2
v_min <- .001
mu <- seq(.001,.2,(v_max-v_min)/Ndrift)
mu <- c(-mu,mu)

# Simulation parameters
dt <- .001; nsim <- 100; ev_bound <- .3; ev_window <- dt*10; upperRT <- 5
timesteps <- upperRT/dt; ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)

# Generate the heatmap
go_to("results")
filename <- "hm_ref_unif.Rdata"
if (file.exists(filename)) {
  load(filename)
}else{
  output <- RW_createHM_driftsign(mu, dt=dt, nsim=nsim, ev_bound=ev_bound, ev_window=ev_window, upperRT=upperRT)
  save(output,file=filename)
}

hm_low <- output$lower; hm_up <- output$upper
hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)

colMap <- viridis::viridis(length(hm_up))
go_to("plots")
# Plot reference HM
par(mar=c(3,4,1.5,2))
jpeg(filename = "pcor_hm_axis.jpeg",
     width = 8,
     height = 4,
     units = "cm",
     res = 1200)
par(mar=c(2,2.7,1,1.55))
par(family="A")
fields::image.plot(1:dim(hm_up)[1],1:dim(hm_up)[2],hm_up,zlim=c(0,1),
                   col=colMap,ylab="",xlab='',legend.shrink=.5,
                   main=paste(""),legend.cex=cex.leg,
                   axes=F,cex.main=cexmain,legend.mar=2.5,
                   legend.args = list(text=expression(paste(italic("P"),"(correct)")),cex=cexlab,side = 3,line=.25),
                   axis.args=list(at=seq(0,1,.5),labels=seq(0,1,.5),cex.axis=cexax))
mtext("Evidence",2,at=dim(hm_up)[2]/2,line=2,cex=cexlab);
mtext("Time (s)",1,at=dim(hm_up)[1]/2,line=1,cex=cexlab)
axis(1, at = c(1,dim(hm_up)[1]), labels = c(0,5),cex.axis=cexax)
axis(2, at = c(1,dim(hm_up)[2]/2,dim(hm_up)[2]), labels = c(-ev_bound,0,ev_bound),cex.axis=cexax)
dev.off()
par(mar=c(5,4,4,2)+.1)

# Generate some data from typical DDM parameters
nsub <- 100
bound <- rnorm(nsub,.07,.02)
ter <- rnorm(nsub,.4,.1)
v <- seq(.001,.2,(v_max-v_min)/Ndrift*100)
z <- 0
ntrials <- 1000
sigma <- .1
confRT <- .5
v_ratio <- 1
rm(simdat)
for (subj in 1:nsub) {
  print(paste("Simulating sub",subj,"out of",nsub))
  if(!exists("simdat")){
    simdat <- chi_square_optim_DDM(c(bound[subj],ter[subj],z,ntrials,sigma,dt,confRT,v_ratio,v),NULL,0)
    simdat$sub <- subj
  }else{
    temp <- chi_square_optim_DDM(c(bound[subj],ter[subj],z,ntrials,sigma,dt,confRT,v_ratio,v),NULL,0)
    temp$sub <- subj
    simdat <- rbind(simdat,temp)
  }
}

# Match evidence and time to the heatmap indexing
simdat$closest_evdnc2 <- match.closest(simdat$evidence2,ev_mapping)
simdat$simdatrt2 <- simdat$rt2;
simdat$simdatrt2[simdat$simdatrt2>5] <- 5 #heatmap doesn't go higher
simdat$simdatrt2 <- simdat$simdatrt2/dt #scale with the heatmap

# confidence = p(correct)
simdat[simdat$resp==1,]$cj <- hmvec_up[(simdat[simdat$resp==1,]$closest_evdnc2-1)*timesteps+round(simdat[simdat$resp==1,]$simdatrt2)]
simdat[simdat$resp==-1,]$cj <- hmvec_low[(simdat[simdat$resp==-1,]$closest_evdnc2-1)*timesteps+round(simdat[simdat$resp==-1,]$simdatrt2)]

LDC <- function(a,b,simdat){
  return(( (1+simdat$resp)/2)*(1 /(1 + exp( (1/sqrt(simdat$rt2+.00001))* (-a*simdat$evidence2 - b))))  + ((1-simdat$resp)/2)*(1 /(1 + exp((1/sqrt(simdat$rt2+.00001))* (a*simdat$evidence2 - b))) ))
}

# Fit LDC to data
fit <- nls(cj~( (1+resp)/2)*(1 /(1 + exp( (1/sqrt(rt2+.00001))* (-a*evidence2 - b))))  + 
      ((1-resp)/2)*(1 /(1 + exp((1/sqrt(rt2+.00001))* (a*evidence2 - b))) ),
    data = simdat,
    start = list(a=15,b=0),
    trace = T)
alpha <- coef(fit)["a"]
beta <- coef(fit)["b"]

simdat$cj_model <- LDC(alpha,beta,simdat)

cor.test(simdat$cj,simdat$cj_model,method = "spearman")

# Draw heatmap from best-fitting LDC parameters and then plot it
hm_ldc <- matrix(NA,nrow=timesteps,ncol=length(ev_mapping))
hm_ldc_low <- matrix(NA,nrow=timesteps,ncol=length(ev_mapping))
hm_ldc_up <- matrix(NA,nrow=timesteps,ncol=length(ev_mapping))
i=1;j=1
for(e in ev_mapping){
  for(t in seq(0,5,length.out=timesteps)){
    hm_ldc[i,j] <- 1/(1 + exp(1/sqrt(t+.00001)*(-alpha*e*custom_sign(e) - beta) ) )
    hm_ldc_up[i,j] <- 1/(1+exp((1/sqrt(t+.00001))*(-alpha * e - beta) ))
    hm_ldc_low[i,j] <- 1/(1+exp((1/sqrt(t+.00001))*(alpha * e - beta) ))
    i<-i+1
  }
  i<-1;j<-j+1
}


jpeg(filename = "pcor_ab_axis.jpeg",
     width = 8,
     height = 4,
     units = "cm",
     res = 1200)
par(mar=c(2,2.7,1,1.55))
par(family="A")
fields::image.plot(1:dim(hm_ldc_up)[1],1:dim(hm_ldc_up)[2],hm_ldc_up,zlim=c(0,1),
                   col=colMap,ylab="",xlab='',legend.shrink=.5,
                   main=expression(paste(alpha," = 15.6; ",beta," = ",0)),
                   axes=F,cex.main=1,legend.cex=cex.leg,legend.mar = 2.5,
                   legend.args = list(text="Confidence",side=3,line=.25,cex=cexlab*.9),
                   axis.args=list(at=seq(0,1,.5),labels=seq(0,1,.5),cex.axis=cexax))
mtext("Evidence",2,at=dim(hm_ldc_up)[2]/2,line=2,cex=cexlab);
mtext("Time (s)",1,at=dim(hm_ldc_up)[1]/2,line=1,cex=cexlab)
axis(1, at = c(1,dim(hm_ldc_up)[1]), labels = c(0,5),cex.axis=cexax)
axis(2, at = c(1,dim(hm_ldc_up)[2]/2,dim(hm_ldc_up)[2]), labels = c(-ev_bound,0,ev_bound),cex.axis=cexax)
dev.off()

# Plot functions and parameters -------------------------------------------

cex.layout <- .83
lwdgr <- 2
cexax <- cex_size(8,cex.layout)
cexlab <- cex_size(10,cex.layout)*cex.layout
cexleg <- cex_size(8,cex.layout)
lwdline <- 2
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

# Plot signature of confidence --------------------------------------------

Nbin <- 10
N <- nsub
Ndiff <- length(v)

simdat$cj_bin <- as.numeric(cut(simdat$cj,breaks = seq(0,1,1/Nbin),labels=1:Nbin))/Nbin
Simuls1 <- simdat
Simuls1$cj_bin <- as.numeric(cut(Simuls1$cj_model,breaks = seq(0,1,1/Nbin),labels=1:Nbin))/Nbin
go_to("plots")
jpeg(filename = "signatures_4.jpeg",
     width = 16,
     height = 16,
     units = "cm",
     res = 1200)
par(family="A")
par(mfrow=c(2,2), mar = c(5,4,0,2))
## Signature 1 ====
#Aggregate conf for data
CJcor <- with(simdat,aggregate(cor,by=list(sub,cj_bin),mean));names(CJcor) <- c('sub','cj_bin','cor')
CJcor <- cast(CJcor,sub~cj_bin,value='cor')

#Immediate condition
simDat_cor <- Simuls1[c('sub','cj_bin','cor')]

#aggregate cj_bin for model
snrcj_binCorSim <- aggregate(.~sub+cj_bin,simDat_cor,mean)
x <- CJcor[,c(2:(Nbin+1))]; n <- length(x)
x_sim = cast(snrcj_binCorSim,sub~cj_bin,value='cor')
x_sim <- x_sim[,c(2:(Nbin+1))]

stripchart(x, ylim=c(0,1), xlim=c(-.05,n+-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=cexax)
mtext("Accuracy",2,at=.5,line=2.5,cex=cexlab);
axis(1,at=0:(n-1),labels=names(x), cex.axis=cexax);
mtext("Confidence",1,2.5,at=(Nbin-1)/2,cex=cexlab)
means <- sapply(x, mean, na.rm=T);n<- length(x)
lines(0:(n-1),colMeans(x_sim,na.rm=T),type='l',lty=1,lwd=lwdgr,col=rgb(0,0,0,.5))
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[n:1]),
        border=F,col=rgb(0,0,0,.2))
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05)

## Signature 2 ====
#Aggregate conf for data
CJcor <- with(subset(simdat,cor==1),aggregate(cj,by=list(sub,drift),mean));names(CJcor) <- c('sub','drift','cj')
CJcor <- cast(CJcor,sub~drift,)
CJerr <- with(subset(simdat,cor==0),aggregate(cj,by=list(sub,drift),mean));names(CJerr) <- c('sub','drift','cj')
CJerr <- cast(CJerr,sub~drift)


#Immediate condition
simDat_cor <- subset(Simuls1,cor==1);simDat_cor <- simDat_cor[c('sub','drift','cj')]
simDat_err <- subset(Simuls1,cor==0);simDat_err <- simDat_err[c('sub','drift','cj')]

#aggregate cj for model
snrcjCorSim <- aggregate(.~sub+drift,simDat_cor,mean)
snrcjErrSim <- aggregate(.~sub+drift,simDat_err,mean)
x <- CJcor[,c(2:(Ndiff+1))];xErr <- CJerr[,c(2:(Ndiff+1))]; n <- length(x)
x_sim = cast(snrcjCorSim,sub~drift,value='cj');xErr_sim = cast(snrcjErrSim,sub~drift,value='cj')

stripchart(x, ylim=c(.4,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=cexax)
mtext("Confidence",2,at=.7,line=2.5,cex=cexlab);
axis(1,at=0:(n-1),labels=round(as.numeric(names(x)),digits=2), cex.axis=cexax);
mtext("Drift rate",1,3,at=5,cex=cexlab)
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),colMeans(x_sim,na.rm=T),type='l',lty="dashed",lwd=lwdgr,col=rgb(26,150,65,128,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim))/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim))/sqrt(N))[Ndiff:1]),
        border=F,col=rgb(26,150,65,51,maxColorValue = 255))
lines(0:(n-1),colMeans(xErr_sim,na.rm=T),type='l',lty="dotted",lwd=lwdgr,col=rgb(215,25,28,128,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xErr_sim,na.rm=T) + (colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N)),(colMeans(xErr_sim,na.rm=T) - colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N))[Ndiff:1]),
        border=F,col=rgb(215,25,28,51,maxColorValue = 255))
legend(.05,1,legend=c("Correct","Error"),lty=c(1,1),bty = "n", cex = cexleg,
       col=c(rgb(26,150,65,maxColorValue = 255),rgb(215,25,28,maxColorValue = 255)))
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col = rgb(26,150,65,maxColorValue = 255))
means <- sapply(xErr, mean,na.rm=T)
error.bar(0:(n-1),means,colSds(as.matrix(xErr),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05, col = rgb(215,25,28,maxColorValue = 255))

## Signature 3 ====
#Aggregate conf for data
CJcor <- with(subset(simdat,cj>=median(simdat$cj,na.rm = T)),aggregate(cor,by=list(sub,drift),mean));names(CJcor) <- c('sub','drift','cor')
CJcor <- cast(CJcor,sub~drift,)
CJerr <- with(subset(simdat,cj<median(simdat$cj,na.rm = T)),aggregate(cor,by=list(sub,drift),mean));names(CJerr) <- c('sub','drift','cor')
CJerr <- cast(CJerr,sub~drift)

#Immediate condition
simDat_cor <- subset(Simuls1,cj>=median(simdat$cj,na.rm = T));simDat_cor <- simDat_cor[c('sub','drift','cor')]
simDat_err <- subset(Simuls1,cj<median(simdat$cj,na.rm = T));simDat_err <- simDat_err[c('sub','drift','cor')]

#aggregate cj for model
snrcjCorSim <- aggregate(.~sub+drift,simDat_cor,mean)
snrcjErrSim <- aggregate(.~sub+drift,simDat_err,mean)
x <- CJcor[,c(2:(Ndiff+1))];xErr <- CJerr[,c(2:(Ndiff+1))]; n <- length(x)
x_sim = cast(snrcjCorSim,sub~drift,value='cor');xErr_sim = cast(snrcjErrSim,sub~drift,value='cor')

stripchart(x, ylim=c(.5,1), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           main="",cex.axis=cexax)
mtext("Accuracy",2,at=.75,line=2.5,cex=cexlab);
axis(1,at=0:(n-1),labels=round(as.numeric(names(x)),digits=2), cex.axis=cexax);
mtext("Drift rate",1,3,at=5,cex=cexlab)
means <- sapply(x, mean);n<- length(x)
lines(0:(n-1),colMeans(x_sim,na.rm=T),type='l',lty=1,lwd=lwdgr,col=rgb(102,194,165,128,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(x_sim,na.rm=T) + (colSds(as.matrix(x_sim),na.rm=T)/sqrt(N)),(colMeans(x_sim,na.rm=T) - colSds(as.matrix(x_sim),na.rm=T)/sqrt(N))[Ndiff:1]),
        border=F,col = rgb(102,194,165,51,maxColorValue = 255))
lines(0:(n-1),colMeans(xErr_sim,na.rm=T),type='l',lty=1,lwd=lwdgr,col=rgb(252,141,98,128,maxColorValue = 255))
polygon(c(0:(n-1),(n-1):0),c(colMeans(xErr_sim,na.rm=T) + (colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N)),(colMeans(xErr_sim,na.rm=T) - colSds(as.matrix(xErr_sim),na.rm=T)/sqrt(N))[Ndiff:1]),
        border=F,col = rgb(252,141,98,51,maxColorValue = 255))
legend("bottomright",legend=c("High confidence","Low confidence"),lty = c(1,1),bty = "n", cex = cexleg,
       col = c(rgb(102,194,165,maxColorValue = 255),rgb(252,141,98,maxColorValue = 255)))
error.bar(0:(n-1),means,colSds(as.matrix(x),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05,col = rgb(102,194,165,maxColorValue = 255))
means <- sapply(xErr, mean,na.rm=T)
error.bar(0:(n-1),means,colSds(as.matrix(xErr),na.rm=T)/sqrt(N),lwd=lwdgr,length=.05, col = rgb(252,141,98,maxColorValue = 255))


# Time ~ CJ ---------------------------------------------------------------
simdat$rt_bin <- cut(simdat$rt,breaks=seq(0,5,.5))
test <- with(simdat,aggregate(cj,by=list(sub,rt_bin),mean))
names(test) <- c("sub","rt_bin","cj")
test <- cast(test,sub~rt_bin)
test <- test[,2:length(test)]

test_sim <- with(simdat,aggregate(cj_model,by=list(sub,rt_bin),mean))
names(test_sim) <- c("sub","rt_bin","cj")
test_sim <- cast(test_sim,sub~rt_bin)
test_sim <- test_sim[,2:length(test_sim)]

n <- length(test)

na_count <-sapply(test, function(y) sum(length(which(!is.na(y)))))

cex.datdot <- 1
stripchart(test,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(.6,.8), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = NULL)
mtext("Confidence",2,at=.7,line=2.5,cex=cexlab);
mtext("RT bin",1,at=n/2,line=2.5,cex=cexlab);
axis(1,at=0:(n-1),labels=names(test), cex.axis=cexax);
axis(2, seq(.6,.8,.05), cex.axis=cexax)
means <- sapply(test, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwdline)
error.bar(0:(n-1),means,colSds(as.matrix(test),na.rm=T)/sqrt(na_count),lwd=lwdline)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim,na.rm=T) + (colSds(as.matrix(test_sim),na.rm=T)/sqrt(na_count)),
          (colMeans(test_sim,na.rm=T) - colSds(as.matrix(test_sim),na.rm=T)/sqrt(na_count))[n:1]),
        border=F,col = rgb(0,0,0,.5))


dev.off()


par(mfrow=c(1,1), mar = c(5,4,4,2)+.1)
