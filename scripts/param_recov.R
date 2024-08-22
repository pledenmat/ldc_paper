#' @title Parameter recovery analysis
#' @author Pierre Le Denmat
#' @description Parameter recovery of the LDC model in the context of the research article: 
#' "A low-dimensional approximation of optimal confidence".
#' @details The script is divided into several sections:
#' - Setup generative parameters: These are randomly sampled from a similar range as the model fits.
#' - Generate data: Generate simulated data from the generative parameters
#' - Model fitting: Load the model fits performed on the cluster (or fit the model if not already done).
#' - Recovery analysis: Look at the correlation between the generative and estimated parameters.
#' + look at the correlation between alpha and beta for possible tradeoff

curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)

library(myPackage)
library(reshape)
library(timeSeries)
setwd("../functions")
source("quantilefit_function_AB.R")
library(gamlss.dist) # Ex-gaussian for generating confRT distribution
setwd("../scripts")


cexmain <- 2;cexax <- 1.5;cexlab <- 2
plots <- F
set.seed(22052023)


# Generative parameters ---------------------------------------------------

nsub <- 200
nsim <- 1
bound <- rnorm(nsub,.07,.02)
ter <- rnorm(nsub,.4,.1)
v1 <- rnorm(nsub,.04,.01)
v2 <- rnorm(nsub,.1,.02)
v3 <- rnorm(nsub,.15,.03)
alpha <- rexGAUS(nsub,15,10,10)
beta <- runif(nsub,-3,2)
z <- 0
ntrials <- 216*3
sigma <- .1
dt <- .001
v_ratio <- 1
t2time <- data.frame(RTconf=rexGAUS(ntrials,mu=.2,sigma=.1,nu=.2)) #' Post-decision accumulation time
t2time[t2time$RTconf<0,"RTconf"] <- .1

# Generate some data from typical DDM parameters ----

setwd("../results")
if (file.exists("data_recov_exp2_ldc.csv")) {
  simdat <- read.csv("data_recov_exp2_ldc.csv")
}else{

  rm(simdat)
  for (subj in 1:nsub) {
    print(paste("Simulating sub",subj,"out of",nsub))
    temp_par <- c(bound[subj],ter[subj],z,nsim,sigma,dt,v_ratio,alpha[subj],beta[subj],v1[subj],v2[subj],v3[subj])
    if(!exists("simdat")){
      simdat <- chi_square_optim_AB(temp_par,t2time,0,condition=1)
      simdat$sub <- subj
    }else{
      temp <- chi_square_optim_AB(temp_par,t2time,0,condition=1)
      temp$sub <- subj
      simdat <- rbind(simdat,temp)
    }
  }
  simdat$difflevel <- rep(c("hard","medium","easy"),each=216)
  simdat$RTconf <- simdat$rt2 - simdat$rt
  # Sanity check
  with(simdat,aggregate(cor,by=list(difflevel),mean))
  
  write.csv(simdat,file = "data_recov_exp2_ldc.csv", row.names = F)
}
simdat2 <- read.csv("hide/data_recov_exp2_ldc.csv")


# Retrieve model fits -----------------------------------------------------

setwd("../fits/recov_exp2")

bound_fit <- rep(NA,nsub)
v_fit <- rep(NA,nsub)
ter_fit <- rep(NA,nsub)
## Adjust the number of drift parameters to the model loaded
v2_fit <- rep(NA,nsub);v3_fit <- rep(NA,nsub)
alpha_fit <- rep(NA,nsub)
beta_fit <- rep(NA,nsub)
resid_fit <- rep(NA,nsub)

for(i in 1:nsub){
  print(paste('Running participant',i,'from',nsub))
  file_name <- paste0('ldc/results_sub_',
                      i,'_AB.Rdata')
  if (file.exists(file_name)) {
    load(file_name)
    if (plots) {
      cost_iter <- results_ab$member$bestvalit[1:results_ab$optim$iter]
      plot(cost_iter, ylab = "Cost function", ylim=c(0, .1), frame = F, 
           type = 'l', main = plot_title)
    }
    bound_fit[i] <- results_ab$optim$bestmem[1]
    ter_fit[i] <- results_ab$optim$bestmem[2]
    alpha_fit[i] <- results_ab$optim$bestmem[8]
    beta_fit[i] <- results_ab$optim$bestmem[9]
    v_fit[i] <- results_ab$optim$bestmem[10]
    v2_fit[i] <-   results_ab$optim$bestmem[11]
    v3_fit[i] <-   results_ab$optim$bestmem[12]
    resid_fit[i] <- results_ab$optim$bestval
  }
}


# Check parameter recovery ------------------------------------------------
# Drift rates (names don't match because difficulty was ordered by alphabetic names
# i.e. easy - hard - medium)
cor.test(v1,v2_fit)
plot(v1,v2_fit,main=paste("r =",round(cor(v1,v2_fit),3)),ylab="Estimated v1",xlab = "Generative v1")
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

cor.test(v2,v3_fit)
plot(v2,v3_fit,main=paste("r =",round(cor(v2,v3_fit),3)),ylab="Estimated v2",xlab = "Generative v2")
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

cor.test(v3,v_fit)
plot(v3,v_fit,main=paste("r =",round(cor(v3,v_fit),3)),ylab="Estimated v3",xlab = "Generative v3")
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

# Bound
cor.test(bound,bound_fit)
plot(bound,bound_fit,main=paste("r =",round(cor(bound,bound_fit),3)),ylab="Estimated Bound",xlab = "Generative Bound")
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

# Non-decision time
cor.test(ter,ter_fit)
plot(ter,ter_fit,main=paste("r =",round(cor(ter,ter_fit),3)),ylab="Estimated Ter",xlab = "Generative Ter")
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

# Alpha
cor.test(alpha[which(alpha>0)],alpha_fit[which(alpha>0)])
plot(alpha[which(alpha>0)],alpha_fit[which(alpha>0)],ylab="Estimated Alpha",xlab = "Generative Alpha",
     main=paste("r =",round(cor(alpha[which(alpha>0)],alpha_fit[which(alpha>0)]),3)))
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

# Beta
cor.test(beta[which(alpha>0)],beta_fit[which(alpha>0)])
plot(beta[which(alpha>0)],beta_fit[which(alpha>0)],ylab="Estimated Beta",xlab = "Generative Beta",
     main=paste("r =",round(cor(beta[which(alpha>0)],beta_fit[which(alpha>0)]),3)))
abline(0,1,col='red')
legend("topleft",bty='n',lty=1,col='red',legend="y = x")

cor.test(alpha,beta)
cor.test(alpha_fit,beta_fit)

plot(alpha[which(alpha>0)],beta[which(alpha>0)],ylab="Generative Beta",xlab = "Generative Alpha",
     main=paste("r =",round(cor(beta[which(alpha>0)],alpha[which(alpha>0)]),3)))

plot(alpha_fit[which(alpha>0)],beta_fit[which(alpha>0)],ylab="Estimated Beta",xlab = "Estimated Alpha",
     main=paste("r =",round(cor(beta_fit[which(alpha>0)],alpha_fit[which(alpha>0)]),3)))
