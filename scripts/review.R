#' ---
#' Additional code addressing reviewers comments 
#' ---

library(viridis)
library(fields)
# Set working directory to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Reproducing Reviewer #2's simulation ------------------------------------

ldc <- function(side,e,t,a,b) {
  return(1 / (1 + exp(- 1/sqrt(t) * (side*a*e + b))) )
}

ldc_wrong <- function(side,e,t,a,b) {
  return(1 / (1 + exp(- 1/sqrt(t) * (side*a*e - b))) )
}

a <- 1
b <- 0
bpos <- 1
bneg <- -1
e <- seq(0,10,.01)
t <- 1
left <- -1
right <- 1

jpeg("ldc_beta.jpg",width=14,height=14,units = 'cm',res= 300)
par(mfrow=c(2,2),mar=c(4,4,2,1))

plot(e,ldc(left,e,t,a,b),type="l",col="red",bty='n',xlab="",ylab="Confidence",
     main="x = -1, beta = positive",ylim=c(0,1))
lines(e,ldc(left,e,t,a,bpos),col="blue")
abline(h=0.5,lty=2)

plot(e,ldc(right,e,t,a,b),type="l",col="red",bty='n',xlab="",ylab="Confidence",
     main="x = 1, beta = positive",ylim=c(.5,1))
lines(e,ldc(right,e,t,a,bpos),col="blue")
abline(h=0.5,lty=2)

plot(e,ldc(left,e,t,a,b),type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="x = -1, beta = negative",ylim=c(0,.5))
lines(e,ldc(left,e,t,a,bneg),col="blue")
abline(h=0.5,lty=2)

plot(e,ldc(right,e,t,a,b),type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="x = 1, beta = negative",ylim=c(0,1))
lines(e,ldc(right,e,t,a,bneg),col="blue")
abline(h=0.5,lty=2)

dev.off()

jpeg("ldc_wrong_beta.jpg",width=14,height=14,units = 'cm',res= 300)
par(mfrow=c(2,2),mar=c(4,4,2,1))

plot(e,ldc_wrong(left,e,t,a,b),type="l",col="red",bty='n',xlab="",ylab="Confidence",
     main="x = -1, beta = positive",ylim=c(0,.5))
lines(e,ldc_wrong(left,e,t,a,bpos),col="blue")
abline(h=0.5,lty=2)

plot(e,ldc_wrong(right,e,t,a,b),type="l",col="red",bty='n',xlab="",ylab="Confidence",
     main="x = 1, beta = positive",ylim=c(0,1))
lines(e,ldc_wrong(right,e,t,a,bpos),col="blue")
abline(h=0.5,lty=2)

plot(e,ldc_wrong(left,e,t,a,b),type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="x = -1, beta = negative",ylim=c(0,1))
lines(e,ldc_wrong(left,e,t,a,bneg),col="blue")
abline(h=0.5,lty=2)

plot(e,ldc_wrong(right,e,t,a,b),type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="x = 1, beta = negative",ylim=c(.5,1))
lines(e,ldc_wrong(right,e,t,a,bneg),col="blue")
abline(h=0.5,lty=2)

dev.off()

# Simulations with time-independent bias ----------------------------------

ldc_bias <- function(side,evidence,time,a,b,c) {
  return(1 / (1 + exp(- 1/sqrt(time) * (side*a*evidence + b) - c)) )
}

# Parameters for simulation
a <- 15
b <- 0
bpos <- 2
bneg <- -2
cpos <- 2
cneg <- -2
c <- 0

dt <- 0.001
maxRT <- 5
t <- seq(0,maxRT,dt)
bound <- .5
ev_resolution <- 0.005
ev <- seq(-bound,bound,ev_resolution)

# Generate heatmap from simulation
sim_pos_c <- matrix(NA,nrow=length(t),ncol=length(ev))
sim_neg_c <- matrix(NA,nrow=length(t),ncol=length(ev))
sim_c <- matrix(NA,nrow=length(t),ncol=length(ev))
sim_pos_b <- matrix(NA,nrow=length(t),ncol=length(ev))
sim_neg_b <- matrix(NA,length(t),length(ev))
for (i in 1:length(t)) {
  for (j in 1:length(ev)) {
    sim_pos_c[i,j] <- ldc_bias(1,ev[j],t[i],a,b,cpos)
    sim_neg_c[i,j] <- ldc_bias(1,ev[j],t[i],a,b,cneg)
    sim_c[i,j] <- ldc_bias(1,ev[j],t[i],a,b,c)
    sim_pos_b[i,j] <- ldc_bias(1,ev[j],t[i],a,bpos,c)
    sim_neg_b[i,j] <- ldc_bias(1,ev[j],t[i],a,bneg,c)
  }
}

colmap <- viridis(10000)
# Plot heatmap using the viridis color map
jpeg("ldc_pos_bias.jpg",width=14,height=14,units = 'cm',res= 300)
par(mar=c(4,4,2,1))
image.plot(t, ev, sim_pos_c, col=colmap, xlab="Time", ylab="Evidence", 
      main="Positive Time-independent Bias")
dev.off()

jpeg("ldc_neg_bias.jpg",width=14,height=14,units = 'cm',res= 300)
par(mar=c(4,4,2,1))
image.plot(t, ev, sim_neg_c, col=colmap, xlab="Time", ylab="Evidence", 
      main="Negative Time-independent Bias")
dev.off()

jpeg("ldc_no_bias.jpg",width=14,height=14,units = 'cm',res= 300)
par(mar=c(4,4,2,1))
image.plot(t, ev, sim_c, col=colmap, xlab="Time", ylab="Evidence", 
      main="No Time-independent Bias")
dev.off()

# Plot confidence at fixed points in time
jpeg("ldc_bias_c.jpg",width=42,height=14,units = 'cm',res= 300)
par(mfrow=c(1,3),mar=c(4,4,2,1))

plot(ev,sim_pos_c[500,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Time-independent Bias at t = 0.5s",ylim=c(0,1))
lines(ev,sim_neg_c[500,],col="blue")
lines(ev,sim_c[500,],type="l",col="black")
# Add legend
legend("bottomright",legend=c("Positive Bias","Negative Bias","No Bias"),col=c("red","blue","black"),lty=1, bty ='n')

plot(ev,sim_pos_c[1000,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Time-independent Bias at t = 1s",ylim=c(0,1))
lines(ev,sim_neg_c[1000,],col="blue")
lines(ev,sim_c[1000,],type="l",col="black")
# Add legend
legend("bottomright",legend=c("Positive Bias","Negative Bias","No Bias"),col=c("red","blue","black"),lty=1, bty ='n')

plot(ev,sim_pos_c[2000,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Time-independent Bias at t = 2s",ylim=c(0,1))
lines(ev,sim_neg_c[2000,],col="blue")
lines(ev,sim_c[2000,],type="l",col="black")
# Add legend
legend("bottomright",legend=c("Positive Bias","Negative Bias","No Bias"),col=c("red","blue","black"),lty=1, bty ='n')
dev.off()

dev.off()

# Make the same plots with the bias in the beta parameter (all in one figure)
jpeg("ldc_bias_beta.jpg",width=42,height=14,units = 'cm',res= 300)
par(mfrow=c(1,3),mar=c(4,4,2,1))
plot(ev,sim_pos_b[500,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Bias in beta at t = 0.5s",ylim=c(0,1))
lines(ev,sim_neg_b[500,],col="blue")
lines(ev,sim_c[500,],type="l",col="black")
# Add legend
legend("bottomright",legend=c("Positive Bias","Negative Bias","No Bias"),col=c("red","blue","black"),lty=1, bty ='n')

plot(ev,sim_pos_b[1000,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Bias in beta at t = 1s",ylim=c(0,1))
lines(ev,sim_neg_b[1000,],col="blue")
lines(ev,sim_c[1000,],type="l",col="black")
# Add legend
legend("bottomright",legend=c("Positive Bias","Negative Bias","No Bias"),col=c("red","blue","black"),lty=1, bty ='n')

plot(ev,sim_pos_b[2000,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Bias in beta at t = 2s",ylim=c(0,1))
lines(ev,sim_neg_b[2000,],col="blue")
lines(ev,sim_c[2000,],type="l",col="black")
# Add legend
legend("bottomright",legend=c("Positive Bias","Negative Bias","No Bias"),col=c("red","blue","black"),lty=1, bty ='n')
dev.off()

cex.lab <- 1.5
cex.axis <- 1.5
cex.main <- 1.5
cex.legend <- 1.5
# Plot biases in c and beta in the same figure
jpeg("ldc_bias_c_beta.jpg",width=42,height=14,units = 'cm',res= 300)
par(mfrow=c(1,3),mar=c(4,4,2,1))

plot(ev,sim_pos_c[500,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Biases at t = 0.5s",ylim=c(0,1),lty=2, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
lines(ev,sim_pos_b[500,],col="red")
lines(ev,sim_neg_c[500,],col="blue",lty=2)
lines(ev,sim_neg_b[500,],col="blue")
lines(ev,sim_c[500,],type="l",col="black")
# Add legend
legend("bottomright",col=c("red","red","blue","blue","black"),lty=c(2,1,2,1,1), bty ='n',
       legend=c("Positive Bias in c","Positive Bias in beta","Negative Bias in c","Negative Bias in beta","No Bias"),
       cex=cex.legend)

plot(ev,sim_pos_c[1000,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Biases at t = 1s",ylim=c(0,1),lty=2, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
lines(ev,sim_pos_b[1000,],col="red")
lines(ev,sim_neg_c[1000,],col="blue",lty=2)
lines(ev,sim_neg_b[1000,],col="blue")
lines(ev,sim_c[1000,],type="l",col="black")
# Add legend
legend("bottomright",col=c("red","red","blue","blue","black"),lty=c(2,1,2,1,1), bty ='n',
       legend=c("Positive Bias in c","Positive Bias in beta","Negative Bias in c","Negative Bias in beta","No Bias"),
       cex=cex.legend)

plot(ev,sim_pos_c[2000,],type="l",col="red",bty='n',xlab="Evidence",ylab="Confidence",
     main="Biases at t = 2s",ylim=c(0,1),lty=2, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main)
lines(ev,sim_pos_b[2000,],col="red")
lines(ev,sim_neg_c[2000,],col="blue",lty=2)
lines(ev,sim_neg_b[2000,],col="blue")
lines(ev,sim_c[2000,],type="l",col="black")
# Add legend
legend("bottomright",col=c("red","red","blue","blue","black"),lty=c(2,1,2,1,1), bty ='n',
       legend=c("Positive Bias in c","Positive Bias in beta","Negative Bias in c","Negative Bias in beta","No Bias"),
       cex=cex.legend)

dev.off()
