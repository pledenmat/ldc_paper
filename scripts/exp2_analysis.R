rm(list=ls())
curdir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curdir)
library(myPackage) # Run devtools::install_github("pledenmat/myPackage") to install this custom package
go_to("functions")
library(Rcpp)
library(numDeriv) #Hessian
library(reshape)
library(lmerTest); 
library(emmeans); 
library(lattice) # qqmath 
library(car) # Stat tests
library(timeSeries) #ColSdS
library(ggplot2)
library(MALDIquant)
library(DEoptim)
source("quantilefit_function_AB.R")
source("quantilefit_function_DDM.R")
source("chi_square_optim_hessian.R")
source("build_hm.R")

plots <- F
stat_tests <- F
pcor <- T

# Function ----------------------------------------------------------------
max_count <- function(data){
  return(max(table(data)))
}

error.bar <- function(x, y, upper, lower=upper, length=0,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
# Experiment 2A preprocessing -----------------------------------------------
##First subjects
N <- 50; 
session1_alpha <- 5; session1_beta <- 4 #Data file was different for first participants
go_to("data")
go_to("exp2a")
##Raw data load
first = T
for(i in 1:session1_alpha){
  if (file.exists(paste0('fbmod_alpha_sub',i,'.csv'))) {
    if (first) {
      Data_alpha <- read.csv(paste0('fbmod_alpha_sub',i,'.csv'))
      first = F
    }else{
      temp <- read.csv(paste0('fbmod_alpha_sub',i,'.csv'))
      Data_alpha <- rbind(Data_alpha,temp)
      
    }
  }
}

##Remove participants that did not complete the experiment
check = 3888 #If participant did all trials then trial_check=check 
trial_check = with(Data_alpha,aggregate(block,by=list(sub=sub),sum))
ok = subset(trial_check,x==check)$sub

Data_alpha <- subset(Data_alpha,sub %in% ok)
subs <- unique(Data_alpha$sub); Nsub <- length(subs)
Ntrials <- dim(Data_alpha)[1]/Nsub; Ntask <- 3


##Add the task-independent difficulty labels
hard <- c(44,49,49); medium <- c(47,57,57); easy <- c(53,62,64)
Data_alpha$difflevel <- 0
Data_alpha[Data_alpha$diff %in% hard,'difflevel'] <- 'hard'
Data_alpha[Data_alpha$diff %in% medium,'difflevel'] <- 'medium'
Data_alpha[Data_alpha$diff %in% easy,'difflevel'] <- 'easy'

##Retrieve the condition
c <- 'control'; m <- 'minus'; p <- 'plus'
condition <- c(c,m,p,c,m,p,m,p,c,m,p,c,m,p,c)
condition <- rep(condition,each=Ntrials/Ntask)
Data_alpha$condition <- condition

Data_alpha$manip <- 'alpha'

for(i in (session1_alpha+1):N){
  if (i %in% c(6,7)) {
    temp <- read.csv(paste0('fbmod_alpha_sub',i,'.csv'))
    Data_alpha <- rbind(Data_alpha,temp)
  }
  else if (file.exists(paste0('fbmod_alpha_sub',i,'.csv'))) {
    temp <- read.csv(paste0('fbmod_alpha_sub',i,'.csv'))
    temp <- temp[ , -which(names(temp) %in% c("questionID","questionResp"))] #Remove questionnaire for now
    temp <- temp[temp$task != "",]
    Data_alpha <- rbind(Data_alpha,temp)
  }
}

subs <- unique(Data_alpha$sub); Nsub <- length(subs)
Ntrials <- dim(Data_alpha)[1]/Nsub; Ntask <- 3
Ntrials_task_main <- dim(subset(Data_alpha,running=="main"))[1]/Nsub/Ntask

Data_alpha$response <- 0
Data_alpha[Data_alpha$resp=="['n']",'response'] <- 1

## Diagnostic plot per participant and task + chance performance testing
exclusion <- c()
tasks <- unique(Data_alpha$task)
par(mfrow=c(2,2))
for(i in 1:Nsub){
  for (t in tasks) {
    tempDat <- subset(Data_alpha,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (plots) {
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt,frame=F,col=c("black"),main=paste('subject',subs[i],"task :",t),ylab="RT")
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
    }
  }
}

#' Filter out participants who reported only one confidence level 
#' more than 90% of the time
conf_count <- with(subset(Data_alpha,running=="main"),aggregate(cj,by=list(sub=sub,task=task),max_count))
conf_count$x <- conf_count$x/Ntrials_task_main
exclusion <- c(exclusion, unique(conf_count[conf_count$x>.9,"sub"]))

Data_alpha <- subset(Data_alpha,!(sub %in% exclusion))

age_alpha <- with(Data_alpha,aggregate(age,by=list(sub),mean))
summary(age_alpha$x)
gender_alpha <- table(Data_alpha$gender)

Data_alpha <- subset(Data_alpha,rt>.2 & rt<5) #Trim RTs 

Data_alpha$response[Data_alpha$response==0] <- -1

Data_alpha$condition <- as.factor(Data_alpha$condition)
Data_alpha$difflevel <- as.factor(Data_alpha$difflevel)
Data_alpha$sub <- as.factor(Data_alpha$sub)


# Experiment 2B preprocessing ------------------------------------------------
go_to("exp2b")
##Raw data load
first = T
for(i in 1:session1_beta){
  if (file.exists(paste0('fbmod_beta_sub',i,'.csv'))) {
    if (first) {
      Data_beta <- read.csv(paste0('fbmod_beta_sub',i,'.csv'))
      first = F
    }else{
      temp <- read.csv(paste0('fbmod_beta_sub',i,'.csv'))
      Data_beta <- rbind(Data_beta,temp)
    }
  }
}

##Remove participants that did not complete the experiment
check = 3888 #If participant did all trials then trial_check=check 
trial_check = with(Data_beta,aggregate(block,by=list(sub=sub),sum))
ok = subset(trial_check,x==check)$sub

Data_beta <- subset(Data_beta,sub %in% ok)
subs <- unique(Data_beta$sub); Nsub <- length(subs)
Ntrials <- dim(Data_beta)[1]/Nsub; Ntask <- 3
Ntrials_task_main <- dim(subset(Data_beta,running=="main"))[1]/Nsub/Ntask

##Add the task-independent difficulty labels
hard <- c(44,49,49); medium <- c(47,57,57); easy <- c(53,62,64)
Data_beta$difflevel <- 0
Data_beta[Data_beta$diff %in% hard,'difflevel'] <- 'hard'
Data_beta[Data_beta$diff %in% medium,'difflevel'] <- 'medium'
Data_beta[Data_beta$diff %in% easy,'difflevel'] <- 'easy'

##Retrieve the condition
c <- 'control'; m <- 'minus'; p <- 'plus'
condition <- c(c,m,p,m,p,c,m,p,c)
condition <- rep(condition,each=Ntrials/Ntask)
Data_beta$condition <- condition

Data_beta$manip <- 'beta'

for(i in (session1_beta+1):N){
  if (file.exists(paste0('fbmod_beta_sub',i,'.csv'))) {
    temp <- read.csv(paste0('fbmod_beta_sub',i,'.csv'))
    temp <- temp[ , -which(names(temp) %in% c("questionID","questionResp"))] #Remove questionnaire for now
    temp <- temp[temp$task != "",]
    Data_beta <- rbind(Data_beta,temp)
  }
}


subs <- unique(Data_beta$sub); Nsub <- length(subs)
Ntrials <- dim(Data_beta)[1]/Nsub; Ntask <- 3

Data_beta$response <- 0
Data_beta[Data_beta$resp=="['n']",'response'] <- 1

with(Data_beta,aggregate(diff,by=list(difflevel,task),mean)) #Quick check


## Diagnostic plot per participant and task + chance performance testing
exclusion <- c()
tasks <- unique(Data_beta$task)
par(mfrow=c(2,2))
for(i in 1:Nsub){
  for (t in tasks) {
    tempDat <- subset(Data_beta,sub==subs[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(response,by=list(block=block),mean))
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subs[i],"p =", round(test$p.value,3),"compared to chance"))
    if (plots) {
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt,frame=F,col=c("black"),main=paste('subject',subs[i],"task :",t),ylab="RT")
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subs[i])
    }
  }
}

#' Filter out participants who reported only one confidence level 
#' more than 90% of the time
conf_count <- with(subset(Data_beta,running=="main"),aggregate(cj,by=list(sub=sub,task=task),max_count))
conf_count$x <- conf_count$x/Ntrials_task_main
exclusion <- c(exclusion, unique(conf_count[conf_count$x>.9,"sub"]))

age_beta <- with(Data_beta,aggregate(age,by=list(sub),mean))
age_beta <- age_beta[age_beta$Group.1 != 14,] #Age is 0
summary(age_beta$x)
gender_beta <- table(Data_beta$gender)

Data_beta <- subset(Data_beta,!(sub %in% exclusion))

Data_beta <- subset(Data_beta,rt>.2 & rt<5) #Trim RTs 

Data_beta$response[Data_beta$response==0] <- -1

Data_beta$condition <- as.factor(Data_beta$condition)
Data_beta$difflevel <- as.factor(Data_beta$difflevel)
Data_beta$sub <- as.factor(Data_beta$sub)


# Merge datasets and prepare data ----------------------------------------------------------
Data_alpha$sub <- as.numeric(as.character(Data_alpha$sub))
Data_alpha$sub <- Data_alpha$sub + 50
Data_alpha$sub <- as.factor(Data_alpha$sub)
Data <- rbind(Data_alpha,Data_beta)
Data <- subset(Data, RTconf < 5)

Data$coh <- as.character(Data$difflevel)
train <- subset(Data,running=="training")
train_alpha <- subset(train,manip=="alpha")
train_beta <- subset(train,manip=="beta")
Data <- subset(Data,running=="main")
Data_alpha <- subset(Data,manip=="alpha")
Data_beta <- subset(Data,manip=="beta")

subs <- sort(unique(Data$sub))
cond <- unique(Data$condition)
difficulty <- sort(unique(Data$coh))
manip <- with(Data,aggregate(manip,by=list(sub=sub),unique)) 
cond_ordered <- c("minus","control","plus")

N <- length(subs)
Ncond <- length(cond)
Nmanip <- length(manip)
Ndiff <- length(difficulty)

# Behaviour analysis --------------------------------------------------------------
if (stat_tests) {
  ### RT
  m.int <- lmer(rt~condition*difflevel + (1|sub),data=Data_alpha) #Random intercept model
  m <- lmer(rt~condition*difflevel + (condition|sub),data=Data_alpha,REML = F)
  anova(m.int,m) #Test the added value of the random slope
  m.diff <- lmer(rt~condition*difflevel + (1+difflevel|sub),data=Data_alpha,REML = F,control=lmerControl(optimize="bobyqa"))
  plot(resid(m),Data_alpha$rt) #Linearity
  leveneTest(residuals(m) ~ Data_alpha$condition) #Homogeneity of variance
  qqmath(m) #Normality
  anova(m) #Results
  
  m.int <- lmer(rt~condition*difflevel + (1|sub),data=Data_beta) #Random intercept model
  m <- lmer(rt~condition*difflevel + (condition|sub),data=Data_beta)
  anova(m.int,m) #Test the added value of the random slope
  m.diff <- lmer(rt~condition*difflevel + (difflevel|sub),data=Data_beta,REML = F)
  plot(resid(m),log(Data_beta$rt)) #Linearity
  leveneTest(residuals(m) ~ Data_beta$condition) #Homogeneity of variance
  qqmath(m) #Normality
  anova(m) #Results
  
  ### ACCURACY
  m.int <- glmer(cor~condition*difflevel + (1|sub),data=Data_alpha,family=binomial)
  m <- glmer(cor~condition*difflevel + (condition|sub),data=Data_alpha,family=binomial); 
  anova(m.int,m)
  m.diff <- glmer(cor~condition*difflevel + (difflevel|sub),data=Data_alpha,family=binomial); #singular
  leveneTest(residuals(m) ~ Data_alpha$condition) #Homogeneity of variance
  Anova(m)
  
  m <- glmer(cor~condition*difflevel + (condition|sub),data=Data_beta,family=binomial); Anova(m)
  m.diff <- glmer(cor~condition*difflevel + (difflevel|sub),data=Data_beta,family=binomial); 
  leveneTest(residuals(m) ~ Data_alpha$condition) #Homogeneity of variance
  
  ### CONFIDENCE
  m.int <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_alpha); 
  mcond <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_alpha); #No main effect when taking whole dataset
  mboth <- lmer(cj ~ condition*cor*difflevel + (cor + condition|sub),data = Data_alpha, REML = F)
  anova(mcond,mboth)
  minteraction <- lmer(cj ~ condition*cor*difflevel + (cor * condition|sub),
                       data = Data_alpha, REML = F,control=lmerControl(optimize="bobyqa"))
  anova(mboth,minteraction)
  minteractiondiff <- lmer(cj ~ condition*cor*difflevel + (cor * condition + difflevel|sub),
                           data = Data_alpha, REML = F,control=lmerControl(optimize="bobyqa"))
  plot(resid(minteraction),Data_alpha$cj) #Linearity
  leveneTest(residuals(minteraction) ~ Data_alpha$condition) #Homogeneity of variance
  qqmath(minteraction) #Normality
  anova(minteraction) #Results
  mcor <- lmer(cj ~ condition*difflevel + (condition|sub),data = subset(Data_alpha,cor==1)); 
  merr <- lmer(cj ~ condition*difflevel + (condition|sub),data = subset(Data_alpha,cor==0)); 
  anova(mcor) #Results
  anova(merr) #Results
  
  m.int <- lmer(cj ~ condition*cor*difflevel + (1|sub),data = Data_beta); 
  mcond <- lmer(cj ~ condition*cor*difflevel + (condition|sub),data = Data_beta,REML = F); 
  anova(m.int,mcond)
  mcor <- lmer(cj ~ condition*cor*difflevel + (cor|sub),data = Data_beta,REML = F);
  anova(mcond,mcor) # m1 wins
  mboth <- lmer(cj ~ condition*cor*difflevel + (cor+condition|sub),data = Data_beta,REML = F); 
  anova(mcond,mboth)
  minteraction <- lmer(cj ~ condition*cor*difflevel + (cor*condition|sub),
                       data = Data_beta,REML = F,control=lmerControl(optimize="bobyqa")); 
  anova(mboth,minteraction)
  minteractiondiff <- lmer(cj ~ condition*cor*difflevel + (cor*condition+difflevel|sub),
                       data = Data_beta,REML = F,control=lmerControl(optimize="bobyqa")); 
  plot(resid(minteractiondiff),Data_beta$cj) #Linearity
  leveneTest(residuals(minteractiondiff) ~ Data_beta$condition) #Homogeneity of variance
  qqmath(minteractiondiff) #Normality
  anova(minteractiondiff)
  mcor <- lmer(cj ~ condition*difflevel + (condition|sub),data = subset(Data_beta,cor==1)); 
  merr <- lmer(cj ~ condition*difflevel + (condition|sub),data = subset(Data_beta,cor==0)); 
  anova(mcor) #Results
  anova(merr) #Results
  
  # 3-way interaction
  mint <- lmer(cj ~ condition*cor*manip + (1|sub:manip), data = Data, REML = F)
  m <- lmer(cj ~ condition*cor*manip + (condition|sub:manip), data = Data, REML = F)
  anova(mint,m)
  mcor <- lmer(cj ~ condition*cor*manip + (cor|sub:manip), data = Data, REML = F)
  anova(m,mcor)
  mboth <- lmer(cj ~ condition*cor*manip*difflevel + (cor+condition|sub:manip), data = Data, REML = F)
  anova(m,mboth)
  minteraction <- lmer(cj ~ condition*cor*manip*difflevel + (cor*condition|sub:manip), data = Data, REML = F)
  anova(mboth,minteraction)
  mslopes <- lmer(cj ~ condition*cor*manip*difflevel + (cor+condition+difflevel|sub:manip),
                  data = Data, REML = F, control=lmerControl(optimize="bobyqa")) # Singular fit
  plot(resid(minteraction),Data$cj)
  leveneTest(residuals(minteraction) ~ Data$condition) #Homogeneity of variance
  qqmath(minteraction) #Normality
  anova(minteraction)
  
  ### FEEDBACK exp2A
  train_alpha$cor <- as.factor(train_alpha$cor)
  m <- lmer(fb ~ condition*cor*difflevel + (1+condition|sub),data = train_alpha)
  plot(resid(m),train_alpha$fb)
  leveneTest(residuals(m) ~ train_alpha$condition) #Homogeneity of variance
  qqmath(m) #Normality
  anova(m)
  
  m.err <- lmer(fb ~ condition*difflevel + (condition|sub),data = subset(train_alpha,cor==0))
  anova(m.err)
  m.cor <- lmer(fb ~ condition*difflevel + (1|sub),data = subset(train_alpha,cor==1))
  anova(m.cor)
  
  ### FEEDBACK exp2B
  train_beta$cor <- as.factor(train_beta$cor)
  m <- lmer(fb ~ condition*difflevel*cor + (condition|sub),data = train_beta)
  plot(resid(m),train_beta$fb)
  leveneTest(residuals(m) ~ train_beta$condition) #Homogeneity of variance
  qqmath(m) #Normality
  anova(m)
  
  m.err <- lmer(fb ~ condition*difflevel + (condition|sub),data = subset(train_beta,cor==0))
  anova(m.err)
  m.cor <- lmer(fb ~ condition*difflevel + (1|sub),data = subset(train_beta,cor==1))
  anova(m.cor)
  
  train_alpha$cor <- as.numeric(train_alpha$cor)
  train_beta$cor <- as.numeric(train_beta$cor)
}

# Load model fits --------------------------------------------------------
go_to("results")
if (file.exists("fitted_parameters_2step_fixed.csv")) {
  param <- read.csv("fitted_parameters_2step_fixed.csv")
}else{
  nsim_fit <- 20
  
  bound_afix <- matrix(NA,N,Ncond)
  v_afix <- matrix(NA,N,Ncond)
  ter_afix <- matrix(NA,N,Ncond)
  ## Adjust the number of drift parameters to the model loaded
  v2_afix <- matrix(NA,N,Ncond);v3_afix <- matrix(NA,N,Ncond)
  alpha_afix <- matrix(NA,N,Ncond)
  beta_afix <- matrix(NA,N,Ncond)
  resid_afix <- matrix(NA,N,Ncond)
  
  bound_bfix <- matrix(NA,N,Ncond)
  v_bfix <- matrix(NA,N,Ncond)
  ter_bfix <- matrix(NA,N,Ncond)
  ## Adjust the number of drift parameters to the model loaded
  v2_bfix <- matrix(NA,N,Ncond);v3_bfix <- matrix(NA,N,Ncond)
  alpha_bfix <- matrix(NA,N,Ncond)
  beta_bfix <- matrix(NA,N,Ncond)
  resid_bfix <- matrix(NA,N,Ncond)
  
  bound_bothfix <- matrix(NA,N,Ncond)
  v_bothfix <- matrix(NA,N,Ncond)
  ter_bothfix <- matrix(NA,N,Ncond)
  ## Adjust the number of drift parameters to the model loaded
  v2_bothfix <- matrix(NA,N,Ncond);v3_bothfix <- matrix(NA,N,Ncond)
  alpha_bothfix <- matrix(NA,N,Ncond)
  beta_bothfix <- matrix(NA,N,Ncond)
  resid_bothfix <- matrix(NA,N,Ncond)
  
  bound_nofix <- matrix(NA,N,Ncond)
  v_nofix <- matrix(NA,N,Ncond)
  ter_nofix <- matrix(NA,N,Ncond)
  ## Adjust the number of drift parameters to the model loaded
  v2_nofix <- matrix(NA,N,Ncond);v3_nofix <- matrix(NA,N,Ncond)
  alpha_nofix <- matrix(NA,N,Ncond)
  beta_nofix <- matrix(NA,N,Ncond)
  resid_nofix <- matrix(NA,N,Ncond)
  
  bound_ddm_per_cond <- matrix(NA,N,Ncond)
  v_ddm_per_cond <- matrix(NA,N,Ncond)
  ter_ddm_per_cond <- matrix(NA,N,Ncond)
  ## Adjust the number of drift parameters to the model loaded
  v2_ddm_per_cond <- matrix(NA,N,Ncond);v3_ddm_per_cond <- matrix(NA,N,Ncond)
  alpha_ddm_per_cond <- matrix(NA,N,Ncond)
  beta_ddm_per_cond <- matrix(NA,N,Ncond)
  resid_ddm_per_cond <- matrix(NA,N,Ncond)
  
  go_to("fits")
  go_to("exp2")
  for(i in 1:N){
    for(c in 1:Ncond){
      print(paste('Running participant',i,'from',N,"condition",c))
      file_name <- paste0('alpha/results_sub_',
                          subs[i],'_',manip[i,'x'],'_AB.Rdata')
      if (file.exists(file_name)) {
        load(file_name)
        if (plots) {
          cost_iter <- results_ab$member$bestvalit[1:results_ab$optim$iter]
          plot(cost_iter, ylab = "Cost function", ylim=c(0, .1), frame = F, 
               type = 'l', main = plot_title)
        }
        bound_afix[i,c] <- results_ab$optim$bestmem[1]
        ter_afix[i,c] <- results_ab$optim$bestmem[2]
        alpha_afix[i,c] <- results_ab$optim$bestmem[8]
        beta_afix[i,c] <- results_ab$optim$bestmem[8+c]
        v_afix[i,c] <- results_ab$optim$bestmem[12]
        v2_afix[i,c] <-   results_ab$optim$bestmem[13]
        v3_afix[i,c] <-   results_ab$optim$bestmem[14]
        resid_afix[i,c] <- results_ab$optim$bestval
      }
      
      file_name <- paste0('beta/results_sub_',
                          subs[i],'_',manip[i,'x'],'_AB.Rdata')
      if (file.exists(file_name)) {
        load(file_name)
        if (plots) {
          cost_iter <- results_ab$member$bestvalit[1:results_ab$optim$iter]
          plot(cost_iter, ylab = "Cost function", ylim=c(0, .1), frame = F, 
               type = 'l', main = plot_title)
        }
        bound_bfix[i,c] <- results_ab$optim$bestmem[1]
        ter_bfix[i,c] <- results_ab$optim$bestmem[2]
        alpha_bfix[i,c] <- results_ab$optim$bestmem[7+c]
        beta_bfix[i,c] <- results_ab$optim$bestmem[11]
        v_bfix[i,c] <- results_ab$optim$bestmem[12]
        v2_bfix[i,c] <-   results_ab$optim$bestmem[13]
        v3_bfix[i,c] <-   results_ab$optim$bestmem[14]
        resid_bfix[i,c] <- results_ab$optim$bestval
      }
      
      
      file_name <- paste0('both/results_sub_',
                          subs[i],'_',manip[i,'x'],'_AB.Rdata')
      if (file.exists(file_name)) {
        load(file_name)
        if (plots) {
          cost_iter <- results_ab$member$bestvalit[1:results_ab$optim$iter]
          plot(cost_iter, ylab = "Cost function", ylim=c(0, .1), frame = F, 
               type = 'l', main = plot_title)
        }
        bound_bothfix[i,c] <- results_ab$optim$bestmem[1]
        ter_bothfix[i,c] <- results_ab$optim$bestmem[2]
        alpha_bothfix[i,c] <- results_ab$optim$bestmem[8]
        beta_bothfix[i,c] <- results_ab$optim$bestmem[9]
        v_bothfix[i,c] <- results_ab$optim$bestmem[10]
        v2_bothfix[i,c] <-   results_ab$optim$bestmem[11]
        v3_bothfix[i,c] <-   results_ab$optim$bestmem[12]
        resid_bothfix[i,c] <- results_ab$optim$bestval
      }
      
      
      file_name <- paste0('none/results_sub_',
                          subs[i],'_',manip[i,'x'],'_AB.Rdata')
      if (file.exists(file_name)) {
        load(file_name)
        if (plots) {
          cost_iter <- results_ab$member$bestvalit[1:results_ab$optim$iter]
          plot(cost_iter, ylab = "Cost function", ylim=c(0, .1), frame = F, 
               type = 'l', main = plot_title)
        }
        bound_nofix[i,c] <- results_ab$optim$bestmem[1]
        ter_nofix[i,c] <- results_ab$optim$bestmem[2]
        alpha_nofix[i,c] <- results_ab$optim$bestmem[7+c]
        beta_nofix[i,c] <- results_ab$optim$bestmem[10+c]
        v_nofix[i,c] <- results_ab$optim$bestmem[14]
        v2_nofix[i,c] <-   results_ab$optim$bestmem[15]
        v3_nofix[i,c] <-   results_ab$optim$bestmem[16]
        resid_nofix[i,c] <- results_ab$optim$bestval
      }
      
      file_name <- paste0('ddm_per_condition/results_sub_',subs[i],'_',
                          cond_ordered[c],'_',manip[i,'x'],'_DDM.Rdata')
      if (file.exists(file_name)) {
        load(file_name)
        if (plots) {
          cost_iter <- results_ab$member$bestvalit[1:results_ab$optim$iter]
          plot(cost_iter, ylab = "Cost function", ylim=c(0, .1), frame = F, 
               type = 'l', main = plot_title)
        } 
      } else {
        tempDat <- subset(Data,sub==subs[i]&condition==cond_ordered[c])
        optimal_params <- DEoptim(chi_square_optim_DDM, # function to optimize
                                  # a,ter,z,ntrials,sigma,dt,vratio,drift(s)
                                  lower = c(0, 0, 0, nsim_fit*nrow(tempDat), .1, .001,0,0,0),
                                  upper = c(.2, 2, 0, nsim_fit*nrow(tempDat), .1, .001,.5,.5,.5), 
                                  observations = tempDat,control=c(itermax=1000,steptol=100,reltol=.001,NP=50), returnFit = 1)
        results_ddm <- summary(optimal_params)
        #save individual results
        save(results_ddm, file=file_name)
      }
      bound_ddm_per_cond[i,c] <- results_ddm$optim$bestmem[1]
      ter_ddm_per_cond[i,c] <- results_ddm$optim$bestmem[2]
      v_ddm_per_cond[i,c] <- results_ddm$optim$bestmem[7]
      v2_ddm_per_cond[i,c] <-   results_ddm$optim$bestmem[8]
      v3_ddm_per_cond[i,c] <-   results_ddm$optim$bestmem[9]
      resid_ddm_per_cond[i,c] <- results_ddm$optim$bestval
      
    }
  }
  param_alphafixed <- data.frame(drift = c(v_afix,v2_afix,v3_afix),
                                 bound=rep(bound_afix,Ndiff),ter=rep(ter_afix,Ndiff),
                                 sub=rep(subs,Ndiff*Ncond),
                                 condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                                 alpha = rep(alpha_afix,Ndiff),beta = rep(beta_afix,Ndiff),
                                 manip=rep(manip$x,Ndiff*Ncond),
                                 resid = rep(resid_afix,Ndiff),difflevel=rep(difficulty,each=N*Ncond))
  param_alphafixed$fixed_parameter <- "alpha"
  param_alphafixed$Npar <- 4
  
  param_betafixed <- data.frame(drift = c(v_bfix,v2_bfix,v3_bfix),
                                bound=rep(bound_bfix,Ndiff),ter=rep(ter_bfix,Ndiff),
                                sub=rep(subs,Ndiff*Ncond),
                                condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                                alpha = rep(alpha_bfix,Ndiff),beta = rep(beta_bfix,Ndiff),
                                manip=rep(manip$x,Ndiff*Ncond),
                                resid = rep(resid_bfix,Ndiff),difflevel=rep(difficulty,each=N*Ncond))
  param_betafixed$fixed_parameter <- "beta"
  param_betafixed$Npar <- 4
  
  param_bothfixed <- data.frame(drift = c(v_bothfix,v2_bothfix,v3_bothfix),
                                bound=rep(bound_bothfix,Ndiff),ter=rep(ter_bothfix,Ndiff),
                                sub=rep(subs,Ndiff*Ncond),
                                condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                                alpha = rep(alpha_bothfix,Ndiff),beta = rep(beta_bothfix,Ndiff),
                                manip=rep(manip$x,Ndiff*Ncond),
                                resid = rep(resid_bothfix,Ndiff),difflevel=rep(difficulty,each=N*Ncond))
  param_bothfixed$fixed_parameter <- "both"
  param_bothfixed$Npar <- 2
  
  param_nofixed <- data.frame(drift = c(v_nofix,v2_nofix,v3_nofix),
                              bound=rep(bound_nofix,Ndiff),ter=rep(ter_nofix,Ndiff),
                              sub=rep(subs,Ndiff*Ncond),
                              condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                              alpha = rep(alpha_nofix,Ndiff),beta = rep(beta_nofix,Ndiff),
                              manip=rep(manip$x,Ndiff*Ncond),
                              resid = rep(resid_nofix,Ndiff),difflevel=rep(difficulty,each=N*Ncond))
  param_nofixed$fixed_parameter <- "none"
  param_nofixed$Npar <- 6
  
  param_pcor_cond <- data.frame(drift = c(v_ddm_per_cond,v2_ddm_per_cond,v3_ddm_per_cond),
                           bound=rep(bound_ddm_per_cond,Ndiff),
                           ter=rep(ter_ddm_per_cond,Ndiff),
                           sub=rep(subs,Ndiff*Ncond),
                           condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                           alpha = NA,beta = NA,
                           manip=rep(manip$x,Ndiff*Ncond),
                           resid = rep(resid_ddm_per_cond,Ndiff),difflevel=rep(difficulty,each=N*Ncond))
  param_pcor_cond$fixed_parameter <- "pcor_cond"
  param_pcor_cond$Npar <- 6
  
  param_pcor <- data.frame(drift = c(v_nofix,v2_nofix,v3_nofix),
                                bound=rep(bound_nofix,Ndiff),ter=rep(ter_nofix,Ndiff),
                                sub=rep(subs,Ndiff*Ncond),
                                condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                                alpha = NA,beta = NA,
                                manip=rep(manip$x,Ndiff*Ncond),
                                resid = NA,difflevel=rep(difficulty,each=N*Ncond))
  param_pcor$fixed_parameter <- "pcor"
  param_pcor$Npar <- 0
  
  param <- rbind(param_alphafixed,param_betafixed,param_bothfixed,param_nofixed,
                 param_pcor,param_pcor_cond)
  
  go_to("results")
  write.csv(param,file="fitted_parameters_hm.csv")
  
}
# Simulation from model ---------------------------------------------------
#Generate model simulations
rm(Simuls_afix);rm(Simuls_bfix);rm(Simuls_bothfix);rm(Simuls_nonefix)
rm(Simuls_pcor); rm(Simuls_pcor_cond)
nsim <- 1
# Simulation parameters
dt <- .001; ev_bound <- .5; ev_window <- dt; upperRT <- 5; sigma <- .1
timesteps <- upperRT/dt; ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
dir.create("heatmaps",showWarnings = F)
go_to("heatmaps")
for(i in 1:N){

  print(paste('simulating',i,'from',N))
  temp_dat <- subset(Data,sub==subs[i])
  par_afix <- subset(param,fixed_parameter=="alpha"&sub==subs[i])
  par_afix <- c(mean(par_afix$bound),mean(par_afix$ter),0,nsim,.1,.001,1,mean(par_afix$alpha),
                #Beta for each condition, reordered to have minus,control,plus
                with(par_afix,aggregate(beta,list(condition),mean))[c(2,1,3),"x"],
                with(par_afix,aggregate(drift,list(difflevel),mean))$x) #Drift rate
  par_bfix <- subset(param,fixed_parameter=="beta"&sub==subs[i])
  par_bfix <- c(mean(par_bfix$bound),mean(par_bfix$ter),0,nsim,.1,.001,1,
                #Alpha for each condition, reordered to have minus,control,plus
                with(par_bfix,aggregate(alpha,list(condition),mean))[c(2,1,3),"x"],
                mean(par_bfix$beta),
                with(par_bfix,aggregate(drift,list(difflevel),mean))$x) #Drift rate
  par_bothfix <- subset(param,fixed_parameter=="both"&sub==subs[i])
  par_bothfix <- c(mean(par_bothfix$bound),mean(par_bothfix$ter),0,nsim,.1,.001,1,
                   mean(par_bothfix$alpha),mean(par_bothfix$beta),
                   with(par_bothfix,aggregate(drift,list(difflevel),mean))$x) #Drift rate
  par_nonefix <- subset(param,fixed_parameter=="none"&sub==subs[i])
  par_nonefix <- c(mean(par_nonefix$bound),mean(par_nonefix$ter),0,nsim,.1,.001,1,
                   #Alpha for each condition, reordered to have minus,control,plus
                   with(par_nonefix,aggregate(alpha,list(condition),mean))[c(2,1,3),"x"],
                   #Beta for each condition, reordered to have minus,control,plus
                   with(par_nonefix,aggregate(beta,list(condition),mean))[c(2,1,3),"x"],
                   with(par_nonefix,aggregate(drift,list(difflevel),mean))$x) #Drift rate
  if (pcor) {
    if (i==58) {
      next
    }
    par_pcor <- subset(param,sub==subs[i])
    par_pcor <- c(mean(par_pcor$bound),mean(par_pcor$ter),0,nsim,.1,.001,1,
                  with(par_pcor,aggregate(drift,list(difflevel),mean))$x) #Drift rate
    for (c in 1:Ncond) {
      temp_dat_cond <- subset(temp_dat,condition==cond_ordered[c])
      par_pcor_cond <- subset(param,fixed_parameter=="pcor_cond"&sub==subs[i]&condition==cond_ordered[c])
      par_pcor_cond <- c(mean(par_pcor_cond$bound),mean(par_pcor_cond$ter),0,nsim,.1,.001,1,
                         with(par_pcor_cond,aggregate(drift,list(difflevel),mean))$x) #Drift rate
      
      temp_pcor_cond <- chi_square_optim_DDM_fullconfRT(par_pcor_cond,temp_dat_cond,returnFit = 0)
      file_name <- paste0("hm_",subs[i],"_",cond_ordered[c],".Rdata")
      if (file.exists(file_name)) {
        load(file_name)
      } else {
        hm_up <- build_hm(par_pcor_cond[(length(par_pcor_cond)-Ndiff+1):length(par_pcor_cond)],
                          sigma = sigma, ev_bound = ev_bound, dt = dt, 
                          ev_window = ev_window, upperRT = upperRT)
        save(hm_up,file=file_name)
      }
      hm_low <- 1-hm_up
      hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
      
      temp_pcor_cond$closest_evdnc2 <- match.closest(temp_pcor_cond$evidence2,ev_mapping)
      temp_pcor_cond$rt2[temp_pcor_cond$rt2>5] <- 5 #heatmap doesn't go higher
      temp_pcor_cond$rt2 <- temp_pcor_cond$rt2*timesteps/upperRT #scale with the heatmap
      temp_pcor_cond[temp_pcor_cond$resp==1,]$cj <- hmvec_up[(temp_pcor_cond[temp_pcor_cond$resp==1,]$closest_evdnc2-1)*timesteps+round(temp_pcor_cond[temp_pcor_cond$resp==1,]$rt2)]
      temp_pcor_cond[temp_pcor_cond$resp==-1,]$cj <- hmvec_low[(temp_pcor_cond[temp_pcor_cond$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp_pcor_cond[temp_pcor_cond$resp==-1,]$rt2)]
      
      if(!exists('Simuls_pcor_cond')){ Simuls_pcor_cond <- cbind(temp_pcor_cond,subs[i],manip[i,'x'],cond_ordered[c])
      }else{ Simuls_pcor_cond <- rbind(Simuls_pcor_cond,cbind(temp_pcor_cond,subs[i],manip[i,'x'],cond_ordered[c]))
      }
    }
    temp_pcor <- chi_square_optim_DDM_fullconfRT(par_pcor,temp_dat,returnFit = 0)
    file_name <- paste0("hm_",subs[i],".Rdata")
    if (file.exists(file_name)) {
      load(file_name)
    } else {
      hm_up <- build_hm(par_pcor[(length(par_pcor)-Ndiff+1):length(par_pcor)],
                        sigma = sigma, ev_bound = ev_bound, dt = dt, 
                        ev_window = ev_window, upperRT = upperRT)
      save(hm_up,file=file_name)
    }
    hm_low <- 1-hm_up
    hmvec_low <- as.vector(hm_low); hmvec_up <- as.vector(hm_up)
    
    temp_pcor$closest_evdnc2 <- match.closest(temp_pcor$evidence2,ev_mapping)
    temp_pcor$rt2[temp_pcor$rt2>5] <- 5 #heatmap doesn't go higher
    temp_pcor$rt2 <- temp_pcor$rt2*timesteps/upperRT #scale with the heatmap
    temp_pcor[temp_pcor$resp==1,]$cj <- hmvec_up[(temp_pcor[temp_pcor$resp==1,]$closest_evdnc2-1)*timesteps+round(temp_pcor[temp_pcor$resp==1,]$rt2)]
    temp_pcor[temp_pcor$resp==-1,]$cj <- hmvec_low[(temp_pcor[temp_pcor$resp==-1,]$closest_evdnc2-1)*timesteps+round(temp_pcor[temp_pcor$resp==-1,]$rt2)]
    if(!exists('Simuls_pcor')){ Simuls_pcor <- cbind(temp_pcor,subs[i],manip[i,'x'])
    }else{ Simuls_pcor <- rbind(Simuls_pcor,cbind(temp_pcor,subs[i],manip[i,'x']))
    }
    
  }
  
  
  temp_afix <- chi_square_optim_AB_afix(par_afix,temp_dat,0)
  temp_bfix <- chi_square_optim_AB_bfix(par_bfix,temp_dat,0)
  temp_bothfix <- chi_square_optim_AB(par_bothfix,temp_dat,0)
  temp_nonefix <- chi_square_optim_AB_nofix(par_nonefix,temp_dat,0)
  
  
  if(!exists('Simuls_nonefix')){ Simuls_nonefix <- cbind(temp_nonefix,subs[i],manip[i,'x'])
  }else{ Simuls_nonefix <- rbind(Simuls_nonefix,cbind(temp_nonefix,subs[i],manip[i,'x']))
  }
  if(!exists('Simuls_bothfix')){ Simuls_bothfix <- cbind(temp_bothfix,subs[i],manip[i,'x'])
  }else{ Simuls_bothfix <- rbind(Simuls_bothfix,cbind(temp_bothfix,subs[i],manip[i,'x']))
  }
  if(!exists('Simuls_bfix')){ Simuls_bfix <- cbind(temp_bfix,subs[i],manip[i,'x'])
  }else{ Simuls_bfix <- rbind(Simuls_bfix,cbind(temp_bfix,subs[i],manip[i,'x']))
  }
  if(!exists('Simuls_afix')){ Simuls_afix <- cbind(temp_afix,subs[i],manip[i,'x'])
  }else{ Simuls_afix <- rbind(Simuls_afix,cbind(temp_afix,subs[i],manip[i,'x']))
  }
  
}
Simuls_afix <- data.frame(Simuls_afix);
Simuls_bfix <- data.frame(Simuls_bfix);
Simuls_bothfix <- data.frame(Simuls_bothfix);
Simuls_nonefix <- data.frame(Simuls_nonefix);
names(Simuls_afix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
names(Simuls_bfix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
names(Simuls_bothfix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
names(Simuls_nonefix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')

coherences <- sort(unique(Data$coh))
Simuls_bfix$coh <- 0
for (i in 1:N) {
  for(d in 1:length(coherences)) Simuls_bfix$coh[Simuls_bfix$sub==subs[i] & Simuls_bfix$drift %in% c(unique(subset(Simuls_bfix,sub==subs[i])$drift)[d],unique(subset(Simuls_bfix,sub==subs[i])$drift)[d+3],unique(subset(Simuls_bfix,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}
Simuls_afix$coh <- 0
for (i in 1:N) {
  for(d in 1:length(coherences)) Simuls_afix$coh[Simuls_afix$sub==subs[i] & Simuls_afix$drift %in% c(unique(subset(Simuls_afix,sub==subs[i])$drift)[d],unique(subset(Simuls_afix,sub==subs[i])$drift)[d+3],unique(subset(Simuls_afix,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}
Simuls_nonefix$coh <- 0
for (i in 1:N) {
  for(d in 1:length(coherences)) Simuls_nonefix$coh[Simuls_nonefix$sub==subs[i] & Simuls_nonefix$drift %in% c(unique(subset(Simuls_nonefix,sub==subs[i])$drift)[d],unique(subset(Simuls_nonefix,sub==subs[i])$drift)[d+3],unique(subset(Simuls_nonefix,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}
Simuls_bothfix$coh <- 0
for (i in 1:N) {
  for(d in 1:length(coherences)) Simuls_bothfix$coh[Simuls_bothfix$sub==subs[i] & Simuls_bothfix$drift %in% c(unique(subset(Simuls_bothfix,sub==subs[i])$drift)[d],unique(subset(Simuls_bothfix,sub==subs[i])$drift)[d+3],unique(subset(Simuls_bothfix,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
}
if (pcor) {
  names(Simuls_pcor) <- c("rt","resp","cor","evidence2","rt2","cj","drift","closest_evdnc2","sub","manip")
  names(Simuls_pcor_cond) <- c("rt","resp","cor","evidence2","rt2","cj","drift","closest_evdnc2","sub","manip","condition")
  
  Simuls_pcor$coh <- 0
  for (i in 1:N) {
    for(d in 1:length(coherences)) Simuls_pcor$coh[Simuls_pcor$sub==subs[i] & Simuls_pcor$drift %in% c(unique(subset(Simuls_pcor,sub==subs[i])$drift)[d],unique(subset(Simuls_pcor,sub==subs[i])$drift)[d+3],unique(subset(Simuls_pcor,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
  }
  Simuls_pcor$condition <- cond_ordered
  
  Simuls_pcor <- Simuls_pcor[Simuls_pcor$rt < 5,]
  
  Simuls_pcor$cj_raw <- Simuls_pcor$cj
  Simuls_pcor$cj <- as.numeric(cut(Simuls_pcor$cj,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))
  Simuls_pcor_cond$cj_raw <- Simuls_pcor_cond$cj
  Simuls_pcor_cond$cj <- as.numeric(cut(Simuls_pcor_cond$cj,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))
  Simuls_pcor_alpha <- subset(Simuls_pcor,manip=="alpha")
  Simuls_pcor_beta <- subset(Simuls_pcor,manip=="beta")
  
}
Simuls_afix <- Simuls_afix[Simuls_afix$rt < 5,]
Simuls_bfix <- Simuls_bfix[Simuls_bfix$rt < 5,]
Simuls_nonefix <- Simuls_nonefix[Simuls_nonefix$rt < 5,]
Simuls_bothfix <- Simuls_bothfix[Simuls_bothfix$rt < 5,]

Simuls_bfix_alpha <- subset(Simuls_bfix,manip=="alpha")
Simuls_bfix_beta <- subset(Simuls_bfix,manip=="beta")
Simuls_afix_alpha <- subset(Simuls_afix,manip=="alpha")
Simuls_afix_beta <- subset(Simuls_afix,manip=="beta")
Simuls_nonefix_alpha <- subset(Simuls_nonefix,manip=="alpha")
Simuls_nonefix_beta <- subset(Simuls_nonefix,manip=="beta")
Simuls_bothfix_alpha <- subset(Simuls_bothfix,manip=="alpha")
Simuls_bothfix_beta <- subset(Simuls_bothfix,manip=="beta")
# Compute MSE pcor models -------------------------------------------------
if (pcor) {
  condition <- cond_ordered
  cost_pcor <- data.frame(sub=rep(subs,Ndiff*Ncond),
                          condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                          manip=rep(manip$x,Ndiff*Ncond),
                          resid = NA,difflevel=rep(difficulty,each=N*Ncond))
  cost_pcor_cond <- data.frame(sub=rep(subs,Ndiff*Ncond),
                               condition=rep(cond_ordered,each=N,length.out=N*Ncond*Ndiff),
                               manip=rep(manip$x,Ndiff*Ncond),
                               resid = NA,difflevel=rep(difficulty,each=N*Ncond))
  for (i in 1:N) {
    predictions <- subset(Simuls_pcor,sub==subs[i])
    predictions_cond <- subset(Simuls_pcor_cond,sub==subs[i])
    observations <- subset(Data,sub==subs[i])
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    pred_props_cond <- NULL; pred_props_cond_cj <- NULL
    for (cond in 1:length(condition)) {
      
      # Separate correct and error trials
      c_predicted <- predictions[predictions$cor == 1,]
      e_predicted <- predictions[predictions$cor == 0,]
      
      c_observed <- observations[observations$cor == 1,]
      e_observed <- observations[observations$cor == 0,]
      
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond],]$cj
      e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond],]$cj  
      c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
      e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
      
      c_props_cj <- as.data.frame(table(c_observed[c_observed$condition == condition[cond],]$cj)/dim(observations)[1])
      e_props_cj <- as.data.frame(table(e_observed[e_observed$condition == condition[cond],]$cj)/dim(observations)[1])
      c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
      e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
      obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
      
      c_predicted_cj <- factor(c_predicted_cj,levels=1:6)
      e_predicted_cj <- factor(e_predicted_cj,levels=1:6)
      
      c_pred_proportion_cj <- table(c_predicted_cj)/ dim(predictions)[1]
      
      e_pred_proportion_cj <- table(e_predicted_cj)/ dim(predictions)[1]
      
      pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)

      ### DO THE SAME WITH THE MODEL PER CONDITION
      c_predicted_cond <- predictions_cond[predictions_cond$cor == 1,]
      e_predicted_cond <- predictions_cond[predictions_cond$cor == 0,]
      
      c_predicted_cond_cj <- c_predicted_cond[c_predicted_cond$condition == condition[cond],]$cj
      e_predicted_cond_cj <- e_predicted_cond[e_predicted_cond$condition == condition[cond],]$cj  
      
      
      c_predicted_cond_cj <- factor(c_predicted_cond_cj,levels=1:6)
      e_predicted_cond_cj <- factor(e_predicted_cond_cj,levels=1:6)
      
      c_pred_proportion_cond_cj <- table(c_predicted_cond_cj)/ dim(predictions_cond)[1]
      
      e_pred_proportion_cond_cj <- table(e_predicted_cond_cj)/ dim(predictions_cond)[1]
      
      pred_props_cond_cj <- c(pred_props_cond_cj,c_pred_proportion_cond_cj,e_pred_proportion_cond_cj)
    }
    
    # Combine the quantiles for rts and cj
    obs_props <- c(obs_props,obs_props_cj)
    pred_props <- c(pred_props,pred_props_cj)
    pred_props_cond <- c(pred_props_cond,pred_props_cond_cj)
    
    # Calculate cost
    cost_pcor_cond[cost_pcor_cond$sub==subs[i],"resid"] = sum( (obs_props - pred_props_cond) ^ 2)
    cost_pcor[cost_pcor$sub==subs[i],"resid"] = sum( (obs_props - pred_props) ^ 2)
    param[param$sub==subs[i]&param$fixed_parameter=="pcor_cond","resid"] = sum( (obs_props - pred_props_cond) ^ 2)
    param[param$sub==subs[i]&param$fixed_parameter=="pcor","resid"] = sum( (obs_props - pred_props) ^ 2)
  }
  
}
# Model Comparison --------------------------------------------------------
bic_custom <- function(Residuals,k,n){
  return(log(n)*k+n*log(Residuals/n))
}

param$Ndata_points <- 12*Ncond # 6 CJ quantiles for correct/error
param$bic <- bic_custom(param$resid,param$Npar,param$Ndata_points)

param_nona <- param[complete.cases(param$bic),]
dim(param_nona)==dim(param)
mean_bic <- with(param_nona,aggregate(bic,by=list(fixed_parameter=fixed_parameter,manip=manip),mean))
mean_resid <- with(param_nona,aggregate(resid,by=list(fixed_parameter=fixed_parameter,manip=manip),mean))

mean_bic$delta <- -99
mean_bic[mean_bic$manip=="alpha","delta"] <- 
  mean_bic[mean_bic$manip=="alpha",]$x - 
  min(mean_bic[mean_bic$manip=="alpha",]$x)
mean_bic[mean_bic$manip=="beta","delta"] <- 
  mean_bic[mean_bic$manip=="beta",]$x -
  min(mean_bic[mean_bic$manip=="beta",]$x)
# Compute confidence contrast ------------------------------------------------

#Generate model simulations
nsim <- 1
nrep <- 500
setwd("..")
if (!file.exists("conf_contrast.Rdata")) {
  temp_mat <- matrix(nrow=nrep,ncol=8)
  rep <- 1
  while (rep <= nrep) {
    # Simulation from model ---------------------------------------------------
    rm(Simuls_afix);rm(Simuls_bfix)
    if (file.exists(paste0("simuls/Simuls_afix_",rep,".csv"))) {
      print(paste("Loading prediction n?",rep))
      Simuls_afix <- read.csv(paste0("simuls/Simuls_afix_",rep,".csv"))
      Simuls_bfix <- read.csv(paste0("simuls/Simuls_bfix_",rep,".csv"))
    }else{
      for(i in 1:N){
        print(paste('simulating',i,'from',N))
        temp_dat <- subset(Data,sub==subs[i])
        par_afix <- subset(param,fixed_parameter=="alpha"&sub==subs[i])
        par_afix <- c(mean(par_afix$bound),mean(par_afix$ter),0,nsim,.1,.001,1,mean(par_afix$alpha),
                      #Beta for each condition, reordered to have minus,control,plus
                      with(par_afix,aggregate(beta,list(condition),mean))[c(2,1,3),"x"], 
                      with(par_afix,aggregate(drift,list(difflevel),mean))$x) #Drift rate
        par_bfix <- subset(param,fixed_parameter=="beta"&sub==subs[i])
        par_bfix <- c(mean(par_bfix$bound),mean(par_bfix$ter),0,nsim,.1,.001,1,
                      #Alpha for each condition, reordered to have minus,control,plus
                      with(par_bfix,aggregate(alpha,list(condition),mean))[c(2,1,3),"x"],
                      mean(par_bfix$beta),
                      with(par_bfix,aggregate(drift,list(difflevel),mean))$x) #Drift rate
        
        
        temp_afix <- chi_square_optim_AB_afix(par_afix,temp_dat,0)
        temp_bfix <- chi_square_optim_AB_bfix(par_bfix,temp_dat,0)
        if(!exists('Simuls_bfix')){ Simuls_bfix <- cbind(temp_bfix,subs[i],manip[i,'x'])
        }else{ Simuls_bfix <- rbind(Simuls_bfix,cbind(temp_bfix,subs[i],manip[i,'x']))
        }
        if(!exists('Simuls_afix')){ Simuls_afix <- cbind(temp_afix,subs[i],manip[i,'x'])
        }else{ Simuls_afix <- rbind(Simuls_afix,cbind(temp_afix,subs[i],manip[i,'x']))
        }
      }
      Simuls_afix <- data.frame(Simuls_afix);
      Simuls_bfix <- data.frame(Simuls_bfix);
      names(Simuls_afix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
      names(Simuls_bfix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
      names(Simuls_bothfix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
      names(Simuls_nonefix) <- c('rt','resp','cor','evidence2','rt2', 'cj','drift','condition','cj_1','sub','manip')
      Simuls_afix$nrep <- rep
      Simuls_bfix$nrep <- rep
      
      coherences <- sort(unique(Data$coh))
      Simuls_bfix$coh <- 0
      for (i in 1:N) {
        for(d in 1:length(coherences)) Simuls_bfix$coh[Simuls_bfix$sub==subs[i] & Simuls_bfix$drift %in% c(unique(subset(Simuls_bfix,sub==subs[i])$drift)[d],unique(subset(Simuls_bfix,sub==subs[i])$drift)[d+3],unique(subset(Simuls_bfix,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
      }
      Simuls_afix$coh <- 0
      for (i in 1:N) {
        for(d in 1:length(coherences)) Simuls_afix$coh[Simuls_afix$sub==subs[i] & Simuls_afix$drift %in% c(unique(subset(Simuls_afix,sub==subs[i])$drift)[d],unique(subset(Simuls_afix,sub==subs[i])$drift)[d+3],unique(subset(Simuls_afix,sub==subs[i])$drift)[d+6])] <- coherences[d] #recode drift to coherence
      }
      Simuls_afix <- Simuls_afix[Simuls_afix$rt < 5,]
      Simuls_bfix <- Simuls_bfix[Simuls_bfix$rt < 5,]
      write.csv(Simuls_afix,file = paste0("simuls/Simuls_afix_",rep,".csv"))
      write.csv(Simuls_bfix,file = paste0("simuls/Simuls_bfix_",rep,".csv"))
      
    }
    
    Simuls_bfix_alpha <- subset(Simuls_bfix,manip=="alpha")
    Simuls_bfix_beta <- subset(Simuls_bfix,manip=="beta")
    Simuls_afix_alpha <- subset(Simuls_afix,manip=="alpha")
    Simuls_afix_beta <- subset(Simuls_afix,manip=="beta")
    
    # Model comparison - predictions ------------------------------------------------
    pred_means_a_afix <- with(Simuls_afix_alpha,
                              aggregate(cj,by=list(
                                cor=cor,
                                condition=condition,
                                difflevel=coh,
                                sub=sub),mean))
    pred_means_a_afix <- cast(pred_means_a_afix,sub + difflevel+ cor ~ condition)
    pred_means_a_afix$minus_vs_rest <- rowMeans(pred_means_a_afix[,c("plus","control")]) - pred_means_a_afix$minus
    diff_a_afix <- with(pred_means_a_afix,aggregate(minus_vs_rest,list(cor=cor,sub=sub),mean))
    temp_mat[rep,1:2] <- with(diff_a_afix,aggregate(x,list(cor=cor),mean))$x
    
    pred_means_a_bfix <- with(Simuls_bfix_alpha,aggregate(cj,by=list(cor=cor,condition=condition,difflevel=coh,sub=sub),mean))
    pred_means_a_bfix <- cast(pred_means_a_bfix,sub + difflevel+ cor ~ condition)
    pred_means_a_bfix$minus_vs_rest <- rowMeans(pred_means_a_bfix[,c("plus","control")]) - pred_means_a_bfix$minus
    diff_a_bfix <- with(pred_means_a_bfix,aggregate(minus_vs_rest,list(cor=cor,sub=sub),mean))
    temp_mat[rep,3:4] <- with(diff_a_bfix,aggregate(x,list(cor=cor),mean))$x
    
    pred_means_b_afix <- with(Simuls_afix_beta,aggregate(cj,by=list(cor=cor,condition=condition,difflevel=coh,sub=sub),mean))
    pred_means_b_afix <- cast(pred_means_b_afix,sub + difflevel+ cor ~ condition)
    pred_means_b_afix$minus_vs_rest <- rowMeans(pred_means_b_afix[,c("plus","control")]) - pred_means_b_afix$minus
    diff_b_afix <- with(pred_means_b_afix,aggregate(minus_vs_rest,list(cor=cor,sub=sub),mean))
    temp_mat[rep,5:6] <- with(diff_b_afix,aggregate(x,list(cor=cor),mean))$x
    
    pred_means_b_bfix <- with(Simuls_bfix_beta,aggregate(cj,by=list(cor=cor,condition=condition,difflevel=coh,sub=sub),mean))
    pred_means_b_bfix <- cast(pred_means_b_bfix,sub + difflevel+ cor ~ condition)
    pred_means_b_bfix$minus_vs_rest <- rowMeans(pred_means_b_bfix[,c("plus","control")]) - pred_means_b_bfix$minus
    diff_b_bfix <- with(pred_means_b_bfix,aggregate(minus_vs_rest,list(cor=cor,sub=sub),mean))
    temp_mat[rep,7:8] <- with(diff_b_bfix,aggregate(x,list(cor=cor),mean))$x
    rep <- rep + 1
  }
  conf_contrast <- data.frame(cor=rep(c(0,1),each=dim(temp_mat)[1],length.out=length(temp_mat)),
                            data = rep(c("alpha","beta"),each=dim(temp_mat)[1]*4,length.out=length(temp_mat)),
                            model = rep(c("beta_free","alpha_free"),each=dim(temp_mat)[1]*2,length.out=length(temp_mat)),
                            cj=as.vector(temp_mat), rep = 1:500)
  
  save(conf_contrast,file="conf_contrast.Rdata")
}else{
  load("conf_contrast.Rdata")
}

# Empirical confidence contrast
obs_means_a <- with(Data_alpha,aggregate(cj,by=list(cor=cor,condition=condition,difflevel=difflevel,sub=sub),mean))
obs_means_a <- cast(obs_means_a,sub + difflevel+ cor ~ condition)
obs_means_a <- obs_means_a[complete.cases(obs_means_a),]

obs_means_a$minus_vs_rest <- rowMeans(obs_means_a[,c("plus","control")]) - obs_means_a$minus
diff_a_obs_allcond <- with(obs_means_a,aggregate(minus_vs_rest,list(cor=cor,sub=sub),mean))

obs_means_b <- with(Data_beta,aggregate(cj,by=list(cor=cor,condition=condition,difflevel=difflevel,sub=sub),mean))
obs_means_b <- cast(obs_means_b,sub + difflevel+ cor ~ condition)
obs_means_b <- obs_means_b[complete.cases(obs_means_b),]

obs_means_b$minus_vs_rest <- rowMeans(obs_means_b[,c("plus","control")]) - obs_means_b$minus
diff_b_obs_allcond <- with(obs_means_b,aggregate(minus_vs_rest,list(cor=cor,sub=sub),mean))

# Figure 5 ----------------------------------------------------------------
go_to("plots")

#' Function to get the right font size in points
cex_size <- function(size,cex.layout) {
  return(size/(par()$ps*cex.layout))
}

# Feedback condition colors
col_minus <- rgb(27,158,119,maxColorValue = 255)
col_minus_shade <- rgb(27,158,119,100,maxColorValue = 255)
col_baseline <- rgb(217,95,2,maxColorValue = 255)
col_baseline_shade <- rgb(217,95,2,100,maxColorValue = 255)
col_plus <- rgb(117,112,179,maxColorValue = 255)
col_plus_shade <- rgb(117,112,179,100,maxColorValue = 255)
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
cex.legend <- cex_size(8,cex.layout)
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

tiff(filename = "figure5.tif",units="cm",width=16,height=22,res=600,compression="lzw")
# layout(matrix(c(1,3,5,4,7,9,
#                 2,3,6,4,8,10),ncol=2),height=c(.2,.1,1,.1,1,1))
# layout(matrix(c(1,3,5,4,7,9,11,
#                 1,3,5,4,7,9,11,
#                 1,3,5,4,7,9,11,
#                 1,3,5,4,7,9,12,
#                 2,3,6,4,8,10,13,
#                 2,3,6,4,8,10,13,
#                 2,3,6,4,8,10,13,
#                 2,3,6,4,8,10,14),ncol=8),height=c(.2,.1,1,.1,1,1,1))
layout(matrix(c(1,3,5,4,7,9,11,
                1,3,5,4,7,9,12,
                2,3,6,4,8,10,13,
                2,3,6,4,8,10,14),ncol=4),height=c(.2,.09,1,.09,1,1,1), widths = c(.7,.3,.7,.3))
par(mar=c(0,0,0,0), cex.main = cex.main)
font <- "Arial"
windowsFonts(A = windowsFont(font))
par(family="A")
plot.new()
legend("bottom",legend=c("Minus","Control","Plus"),
       title = NULL,pch=rep(16,3),bty = "n",inset=0,pt.cex=cex.datdot,
       cex = cex.legend,col=c(col_minus,col_baseline,col_plus), horiz = T)
legend("bottom",legend=c("Minus","Control","Plus"),
       title = NULL,pch=rep(15,3),bty = "n",inset=0,
       cex = cex.legend,horiz = T,pt.cex=cex.legend.square,
       col=c(col_minus_shade,col_baseline_shade,col_plus_shade))
title("Experiment 2A",  line = -1)
plot.new()
legend("bottom",legend=c("Minus","Control","Plus"),
       title = NULL,pch=rep(16,3),bty = "n",inset=0,pt.cex=cex.datdot,
       cex = cex.legend,col=c(col_minus,col_baseline,col_plus), horiz = T)
legend("bottom",legend=c("Minus","Control","Plus"),
       title = NULL,pch=rep(15,3),bty = "n",inset=0,
       cex = cex.legend,horiz = T,pt.cex=cex.legend.square,
       col=c(col_minus_shade,col_baseline_shade,col_plus_shade))
title("Experiment 2B", line = -1)
plot.new()
mtext("Training phase", line = -1, cex = cex.phase)
plot.new()
mtext("Testing phase", line = -1, cex = cex.phase)
par(cex.main = 1)
# Feedback presented Exp2A ------------------------------------------------
N_temp <- length(unique(train_alpha$sub))

fbminus <- with(subset(train_alpha,condition=="minus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbminus) <- c('sub','difflevel','cor','fb')
fbminus_cor <- subset(fbminus,cor==1); fbminus_err <- subset(fbminus,cor==0)
fbminus_cor <- cast(fbminus_cor,sub~difflevel); fbminus_err <- cast(fbminus_err,sub~difflevel)

fbbaseline <- with(subset(train_alpha,condition=="control"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbbaseline) <- c('sub','difflevel','cor','fb')
fbbaseline_cor <- subset(fbbaseline,cor==1); fbbaseline_err <- subset(fbbaseline,cor==0)
fbbaseline_cor <- cast(fbbaseline_cor,sub~difflevel); fbbaseline_err <- cast(fbbaseline_err,sub~difflevel)

fbplus <- with(subset(train_alpha,condition=="plus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbplus) <- c('sub','difflevel','cor','fb')
fbplus_cor <- subset(fbplus,cor==1); fbplus_err <- subset(fbplus,cor==0)
fbplus_cor <- cast(fbplus_cor,sub~difflevel); fbplus_err <- cast(fbplus_err,sub~difflevel)


# Drop subject column
xminus_cor <- fbminus_cor[,c(2:4)];xbaseline_cor <- fbbaseline_cor[,c(2:4)];xplus_cor <- fbplus_cor[,c(2:4)]
xminus_err <- fbminus_err[,c(2:4)];xbaseline_err <- fbbaseline_err[,c(2:4)];xplus_err <- fbplus_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","medium","easy")];
xbaseline_cor <- xbaseline_cor[,c("hard","medium","easy")];
xplus_cor <- xplus_cor[,c("hard","medium","easy")]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xbaseline_err <- xbaseline_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]

par(mar=mar.fb)

stripchart(xminus_cor,ylim=c(0,100), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "",
           main = NULL,cex.main=cex.main)
title(ylab = "Feedback", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
axis(1,at=0:(n-1),labels=c("Hard","Medium","Easy"), cex.axis=cex.axis);
axis(2, seq(0,100,20), cex.axis=cex.axis,las=ylab_orientation)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

legend("bottomleft",legend=c("Correct trials","Error trials"),
       title = NULL,lty=c("dashed","dotted"),bty = "n",inset=0,
       cex = cex.legend, lwd=lwd.dat,seg.len=1.5)

# Feedback presented Exp2B ------------------------------------------------

N_temp <- length(unique(train_beta$sub))

fbminus <- with(subset(train_beta,condition=="minus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbminus) <- c('sub','difflevel','cor','fb')
fbminus_cor <- subset(fbminus,cor==1); fbminus_err <- subset(fbminus,cor==0)
fbminus_cor <- cast(fbminus_cor,sub~difflevel); fbminus_err <- cast(fbminus_err,sub~difflevel)
fbbaseline <- with(subset(train_beta,condition=="control"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbbaseline) <- c('sub','difflevel','cor','fb')
fbbaseline_cor <- subset(fbbaseline,cor==1); fbbaseline_err <- subset(fbbaseline,cor==0)
fbbaseline_cor <- cast(fbbaseline_cor,sub~difflevel); fbbaseline_err <- cast(fbbaseline_err,sub~difflevel)
fbplus <- with(subset(train_beta,condition=="plus"),aggregate(fb,by=list(sub,difflevel,cor),mean));
names(fbplus) <- c('sub','difflevel','cor','fb')
fbplus_cor <- subset(fbplus,cor==1); fbplus_err <- subset(fbplus,cor==0)
fbplus_cor <- cast(fbplus_cor,sub~difflevel); fbplus_err <- cast(fbplus_err,sub~difflevel)


#aggregate fb for model
xminus_cor <- fbminus_cor[,c(2:4)];xbaseline_cor <- fbbaseline_cor[,c(2:4)];xplus_cor <- fbplus_cor[,c(2:4)]
xminus_err <- fbminus_err[,c(2:4)];xbaseline_err <- fbbaseline_err[,c(2:4)];xplus_err <- fbplus_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","medium","easy")];
xbaseline_cor <- xbaseline_cor[,c("hard","medium","easy")];
xplus_cor <- xplus_cor[,c("hard","medium","easy")]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xbaseline_err <- xbaseline_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]

par(mar=mar.fb)

stripchart(xminus_cor,ylim=c(0,100), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "",
           main = NULL,cex.main=cex.main)
title(ylab = "Feedback", xlab = "Trial difficulty", line = 2.5,cex.lab=cex.lab.rel)
axis(1,at=0:(n-1),labels=c("Hard","Medium","Easy"), cex.axis=cex.axis);
axis(2, seq(0,100,20), cex.axis=cex.axis,las=ylab_orientation)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

### Exp2A confidence with alpha-free model prediction ====
N_temp <- length(unique(Data_alpha$sub))

cjlow <- with(subset(Data_alpha,condition=="minus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjlow) <- c('sub','difflevel','cor','cj')
cjlow_cor <- subset(cjlow,cor==1); cjlow_err <- subset(cjlow,cor==0)
cjlow_cor <- cast(cjlow_cor,sub~difflevel); cjlow_err <- cast(cjlow_err,sub~difflevel)
cjmed <- with(subset(Data_alpha,condition=="control"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjmed) <- c('sub','difflevel','cor','cj')
cjmed_cor <- subset(cjmed,cor==1); cjmed_err <- subset(cjmed,cor==0)
cjmed_cor <- cast(cjmed_cor,sub~difflevel); cjmed_err <- cast(cjmed_err,sub~difflevel)
cjhigh <- with(subset(Data_alpha,condition=="plus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjhigh) <- c('sub','difflevel','cor','cj')
cjhigh_cor <- subset(cjhigh,cor==1); cjhigh_err <- subset(cjhigh,cor==0)
cjhigh_cor <- cast(cjhigh_cor,sub~difflevel); cjhigh_err <- cast(cjhigh_err,sub~difflevel)

#Simulations
simDat_low <- with(subset(Simuls_bfix_alpha,condition=="minus"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_low) <-c('sub','coh','cor','cj')
simDat_low_cor <- subset(simDat_low,cor==1);simDat_low_err <- subset(simDat_low,cor==0)
simDat_low_cor <- cast(simDat_low_cor,sub~coh);simDat_low_err <- cast(simDat_low_err,sub~coh)
simDat_med <- with(subset(Simuls_bfix_alpha,condition=="control"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_med) <- c('sub','coh','cor','cj')
simDat_med_cor <- subset(simDat_med,cor==1);simDat_med_err <- subset(simDat_med,cor==0)
simDat_med_cor <- cast(simDat_med_cor,sub~coh);simDat_med_err <- cast(simDat_med_err,sub~coh)
simDat_high <- with(subset(Simuls_bfix_alpha,condition=="plus"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_high) <- c('sub','coh','cor','cj')
simDat_high_cor <- subset(simDat_high,cor==1);simDat_high_err <- subset(simDat_high,cor==0)
simDat_high_cor <- cast(simDat_high_cor,sub~coh);simDat_high_err <- cast(simDat_high_err,sub~coh)


#aggregate cj for model
xsim_cor <- simDat_low_cor[,c(2:4)];xsimmed_cor <- simDat_med_cor[,c(2:4)];xsimhigh_cor <- simDat_high_cor[,c(2:4)]
xsim_err <- simDat_low_err[,c(2:4)];xsimmed_err <- simDat_med_err[,c(2:4)];xsimhigh_err <- simDat_high_err[,c(2:4)]
xminus_cor <- cjlow_cor[,c(2:4)];xbaseline_cor <- cjmed_cor[,c(2:4)];xplus_cor <- cjhigh_cor[,c(2:4)]
xminus_err <- cjlow_err[,c(2:4)];xbaseline_err <- cjmed_err[,c(2:4)];xplus_err <- cjhigh_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","medium","easy")];
xbaseline_cor <- xbaseline_cor[,c("hard","medium","easy")];
xplus_cor <- xplus_cor[,c("hard","medium","easy")]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xbaseline_err <- xbaseline_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]
xsim_cor <- xsim_cor[,c("hard","medium","easy")];
xsimmed_cor <- xsimmed_cor[,c("hard","medium","easy")];
xsimhigh_cor <- xsimhigh_cor[,c("hard","medium","easy")]
xsim_err <- xsim_err[,c("hard","medium","easy")];
xsimmed_err <- xsimmed_err[,c("hard","medium","easy")];
xsimhigh_err <- xsimhigh_err[,c("hard","medium","easy")]

par(mar=mar.cj)

stripchart(xminus_cor,ylim=c(3,5.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "",
           main = NULL)
axis(1,at=0:(n-1),labels=c("Hard","Medium","Easy"), cex.axis=cex.axis);
axis(2, seq(3,5.5,.5), cex.axis=cex.axis,las=ylab_orientation)
title(ylab = "Confidence", xlab= "Trial difficulty",line = 2.5,cex.lab=cex.lab.rel)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsim_cor,na.rm=T) + (colSds(as.matrix(xsim_cor))/sqrt(N_temp)),
          (colMeans(xsim_cor,na.rm=T) - colSds(as.matrix(xsim_cor))/sqrt(N_temp))[3:1]),
        border=F,col=col_minus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimmed_cor,na.rm=T) + (colSds(as.matrix(xsimmed_cor))/sqrt(N_temp)),
          (colMeans(xsimmed_cor,na.rm=T) - colSds(as.matrix(xsimmed_cor))/sqrt(N_temp))[3:1]),
        border=F,col=col_baseline_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimhigh_cor,na.rm=T) + (colSds(as.matrix(xsimhigh_cor))/sqrt(N_temp)),
          (colMeans(xsimhigh_cor,na.rm=T) - colSds(as.matrix(xsimhigh_cor))/sqrt(N_temp))[3:1]),
        border=F,col=col_plus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsim_err,na.rm=T) + (colSds(as.matrix(xsim_err))/sqrt(N_temp)),
          (colMeans(xsim_err,na.rm=T) - colSds(as.matrix(xsim_err))/sqrt(N_temp))[3:1]),
        border=F,col=col_minus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimmed_err,na.rm=T) + (colSds(as.matrix(xsimmed_err))/sqrt(N_temp)),
          (colMeans(xsimmed_err,na.rm=T) - colSds(as.matrix(xsimmed_err))/sqrt(N_temp))[3:1]),
        border=F,col=col_baseline_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimhigh_err,na.rm=T) + (colSds(as.matrix(xsimhigh_err))/sqrt(N_temp)),
          (colMeans(xsimhigh_err,na.rm=T) - colSds(as.matrix(xsimhigh_err))/sqrt(N_temp))[3:1]),
        border=F,col=col_plus_shade)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)
### Exp2B confidence with beta-free model prediction ====
N_temp <- length(unique(Data_beta$sub))

cjlow <- with(subset(Data_beta,condition=="minus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjlow) <- c('sub','difflevel','cor','cj')
cjlow_cor <- subset(cjlow,cor==1); cjlow_err <- subset(cjlow,cor==0)
cjlow_cor <- cast(cjlow_cor,sub~difflevel); cjlow_err <- cast(cjlow_err,sub~difflevel)
cjmed <- with(subset(Data_beta,condition=="control"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjmed) <- c('sub','difflevel','cor','cj')
cjmed_cor <- subset(cjmed,cor==1); cjmed_err <- subset(cjmed,cor==0)
cjmed_cor <- cast(cjmed_cor,sub~difflevel); cjmed_err <- cast(cjmed_err,sub~difflevel)
cjhigh <- with(subset(Data_beta,condition=="plus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjhigh) <- c('sub','difflevel','cor','cj')
cjhigh_cor <- subset(cjhigh,cor==1); cjhigh_err <- subset(cjhigh,cor==0)
cjhigh_cor <- cast(cjhigh_cor,sub~difflevel); cjhigh_err <- cast(cjhigh_err,sub~difflevel)

#Simulations
simDat_low <- with(subset(Simuls_afix_beta,condition=="minus"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_low) <-c('sub','coh','cor','cj')
simDat_low_cor <- subset(simDat_low,cor==1);simDat_low_err <- subset(simDat_low,cor==0)
simDat_low_cor <- cast(simDat_low_cor,sub~coh);simDat_low_err <- cast(simDat_low_err,sub~coh)
simDat_med <- with(subset(Simuls_afix_beta,condition=="control"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_med) <- c('sub','coh','cor','cj')
simDat_med_cor <- subset(simDat_med,cor==1);simDat_med_err <- subset(simDat_med,cor==0)
simDat_med_cor <- cast(simDat_med_cor,sub~coh);simDat_med_err <- cast(simDat_med_err,sub~coh)
simDat_high <- with(subset(Simuls_afix_beta,condition=="plus"),aggregate(cj,by=list(sub,coh,cor),mean));
names(simDat_high) <- c('sub','coh','cor','cj')
simDat_high_cor <- subset(simDat_high,cor==1);simDat_high_err <- subset(simDat_high,cor==0)
simDat_high_cor <- cast(simDat_high_cor,sub~coh);simDat_high_err <- cast(simDat_high_err,sub~coh)


#aggregate cj for model
xsim_cor <- simDat_low_cor[,c(2:4)];xsimmed_cor <- simDat_med_cor[,c(2:4)];xsimhigh_cor <- simDat_high_cor[,c(2:4)]
xsim_err <- simDat_low_err[,c(2:4)];xsimmed_err <- simDat_med_err[,c(2:4)];xsimhigh_err <- simDat_high_err[,c(2:4)]
xminus_cor <- cjlow_cor[,c(2:4)];xbaseline_cor <- cjmed_cor[,c(2:4)];xplus_cor <- cjhigh_cor[,c(2:4)]
xminus_err <- cjlow_err[,c(2:4)];xbaseline_err <- cjmed_err[,c(2:4)];xplus_err <- cjhigh_err[,c(2:4)]
n <- length(xminus_err)

xminus_cor <- xminus_cor[,c("hard","medium","easy")];
xbaseline_cor <- xbaseline_cor[,c("hard","medium","easy")];
xplus_cor <- xplus_cor[,c("hard","medium","easy")]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xbaseline_err <- xbaseline_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]
xsim_cor <- xsim_cor[,c("hard","medium","easy")];
xsimmed_cor <- xsimmed_cor[,c("hard","medium","easy")];
xsimhigh_cor <- xsimhigh_cor[,c("hard","medium","easy")]
xsim_err <- xsim_err[,c("hard","medium","easy")];
xsimmed_err <- xsimmed_err[,c("hard","medium","easy")];
xsimhigh_err <- xsimhigh_err[,c("hard","medium","easy")]

par(mar=mar.cj)

stripchart(xminus_cor,ylim=c(3,5.5), xlim=c(-.05,n-1), vertical = TRUE, col="white",frame=F,xaxt='n',
           yaxt = 'n',xlab="",ylab = "",
           main = NULL)
axis(1,at=0:(n-1),labels=c("Hard","Medium","Easy"), cex.axis=cex.axis);
axis(2, seq(3,5.5,.5), cex.axis=cex.axis,las=ylab_orientation)
title(ylab = "Confidence", xlab= "Trial difficulty",line = 2.5,cex.lab=cex.lab.rel)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsim_cor,na.rm=T) + (colSds(as.matrix(xsim_cor))/sqrt(N_temp)),
          (colMeans(xsim_cor,na.rm=T) - colSds(as.matrix(xsim_cor))/sqrt(N_temp))[3:1]),
        border=F,col=col_minus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimmed_cor,na.rm=T) + (colSds(as.matrix(xsimmed_cor))/sqrt(N_temp)),
          (colMeans(xsimmed_cor,na.rm=T) - colSds(as.matrix(xsimmed_cor))/sqrt(N_temp))[3:1]),
        border=F,col=col_baseline_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimhigh_cor,na.rm=T) + (colSds(as.matrix(xsimhigh_cor))/sqrt(N_temp)),
          (colMeans(xsimhigh_cor,na.rm=T) - colSds(as.matrix(xsimhigh_cor))/sqrt(N_temp))[3:1]),
        border=F,col=col_plus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsim_err,na.rm=T) + (colSds(as.matrix(xsim_err))/sqrt(N_temp)),
          (colMeans(xsim_err,na.rm=T) - colSds(as.matrix(xsim_err))/sqrt(N_temp))[3:1]),
        border=F,col=col_minus_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimmed_err,na.rm=T) + (colSds(as.matrix(xsimmed_err))/sqrt(N_temp)),
          (colMeans(xsimmed_err,na.rm=T) - colSds(as.matrix(xsimmed_err))/sqrt(N_temp))[3:1]),
        border=F,col=col_baseline_shade)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(xsimhigh_err,na.rm=T) + (colSds(as.matrix(xsimhigh_err))/sqrt(N_temp)),
          (colMeans(xsimhigh_err,na.rm=T) - colSds(as.matrix(xsimhigh_err))/sqrt(N_temp))[3:1]),
        border=F,col=col_plus_shade)
means <- sapply(xminus_cor, mean)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_cor, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dashed")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_cor),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)

means <- sapply(xminus_err, mean, na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_minus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xminus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_minus)
means <- sapply(xbaseline_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_baseline,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xbaseline_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_baseline)
means <- sapply(xplus_err, mean,na.rm=T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,col=col_plus,lwd=lwd.dat,lty = "dotted")
error.bar(0:(n-1),means,colSds(as.matrix(xplus_err),na.rm=T)/sqrt(N_temp),lwd=lwd.dat,col=col_plus)
# Plot confidence contrast ------------------------------------------------
# Experiment 2A -----------------------------------------------------------
a_free <- subset(conf_contrast,data=="alpha"&model=="alpha_free")
a_free <- cast(a_free,rep~cor,value = "cj")
a_free$interaction <- a_free$`1` - a_free$`0`
b_free <- subset(conf_contrast,data=="alpha"&model=="beta_free")
b_free <- cast(b_free,rep~cor,value = "cj")
b_free$interaction <- b_free$`1` - b_free$`0`
empdat <- cast(diff_a_obs_allcond,sub~cor)
empdat$interaction <- empdat$`1` - empdat$`0`

names(a_free) <- c("rep","error","correct","interaction")
names(b_free) <- c("rep","error","correct","interaction")
names(empdat) <- c("sub","error","correct","interaction")

dat_order <- c("interaction","error","correct")

par(mar=mar.contrast)
plot(a_free$`0`,xlim=c(-.1,.3), ylim=c(0,2), col="white",frame=F,xaxt='n',main="", 
     ylab="",yaxt='n',xlab="");
axis(2,at=y_coord,label=c("Interaction","Error","Correct"),las=1,tick=F,cex.axis=cex.trialtype)
axis(1,at=seq(-.1,.3,.2),label=seq(-.1,.3,.2),cex.axis=cex.axis)
mtext(text = "Trial type",side = 2, line = 5,cex = cex.lab)
mtext(text = "Confidence contrast",side = 1, line = 2.5,cex = cex.lab)
segments(x0=0,y0=0,y1=2,lty=2,col="lightgrey")
for(i in 1:length(y_coord)){
  polygon(density(a_free[,dat_order[i]])$x,(density(a_free[,dat_order[i]])$y/200)+y_coord[i]+adj_ydens,
          col=rgb(240,228,66,maxColorValue = 255),border=F)
  polygon(density(b_free[,dat_order[i]])$x,(density(b_free[,dat_order[i]])$y/200)+y_coord[i]+adj_ydens,
          col=rgb(0,114,178,maxColorValue = 255),border=F)
  points(y=y_coord[i]-.1,mean(empdat[,dat_order[i]]),cex=cex.datdot,pch=19)
  arrows(x0=mean(empdat[,dat_order[i]])+( sd(empdat[,dat_order[i]])/sqrt(length(empdat[,dat_order[i]]))),
         y0=y_coord[i]-.1,x1=mean(empdat[,dat_order[i]])-( sd(empdat[,dat_order[i]])/sqrt(length(empdat[,dat_order[i]]))), 
         angle=90,length=0,lwd=lwd.dat)
}
legend("topleft",border=F,legend=c("Empirical data",expression(paste(alpha,"-free")),expression(paste(beta,"-free"))),
       bty='n',cex=cex.legend,lty = c(1,NA,NA),pch = c(19,15,15),lwd=2,pt.cex=c(1,2,2),
       col=c("black",rgb(240,228,66,maxColorValue = 255),rgb(0,114,178,maxColorValue = 255)))
# Experiment 2B -----------------------------------------------------------

a_free <- subset(conf_contrast,data=="beta"&model=="alpha_free")
a_free <- cast(a_free,rep~cor,value = "cj")
a_free$interaction <- a_free$`1` - a_free$`0`
b_free <- subset(conf_contrast,data=="beta"&model=="beta_free")
b_free <- cast(b_free,rep~cor,value = "cj")
b_free$interaction <- b_free$`1` - b_free$`0`
empdat <- cast(diff_b_obs_allcond,sub~cor)
empdat$interaction <- empdat$`1` - empdat$`0`

names(a_free) <- c("rep","error","correct","interaction")
names(b_free) <- c("rep","error","correct","interaction")
names(empdat) <- c("sub","error","correct","interaction")

dat_order <- c("interaction","error","correct")

par(mar=mar.contrast)
plot(a_free$`0`,xlim=c(-.1,.3), ylim=c(0,2), col="white",frame=F,xaxt="n",main="", 
     cex.lab=cex.lab,cex.axis=cex.axis,ylab="",yaxt='n',xlab="");
axis(2,at=y_coord,label=c("Interaction","Error","Correct"),las=1,tick=F,cex.axis=cex.trialtype)
axis(1,at=seq(-.1,.3,.2),label=seq(-.1,.3,.2),cex.axis=cex.axis)
mtext(text = "Trial type",side = 2, line = 5,cex = cex.lab)
mtext(text = "Confidence contrast",side = 1, line = 2.5,cex = cex.lab)
segments(x0=0,y0=0,y1=2,lty=2,col="lightgrey")
for(i in 1:length(y_coord)){
  polygon(density(a_free[,dat_order[i]])$x,(density(a_free[,dat_order[i]])$y/200)+y_coord[i]+adj_ydens,
          col=rgb(240,228,66,maxColorValue = 255),border=F)
  polygon(density(b_free[,dat_order[i]])$x,(density(b_free[,dat_order[i]])$y/200)+y_coord[i]+adj_ydens,
          col=rgb(0,114,178,maxColorValue = 255),border=F)
  points(y=y_coord[i]-.1,mean(empdat[,dat_order[i]]),cex=cex.datdot,pch=19)
  arrows(x0=mean(empdat[,dat_order[i]])+( sd(empdat[,dat_order[i]])/sqrt(length(empdat[,dat_order[i]]))),
         y0=y_coord[i]-.1,x1=mean(empdat[,dat_order[i]])-( sd(empdat[,dat_order[i]])/sqrt(length(empdat[,dat_order[i]]))), 
         angle=90,length=0,lwd=lwd.dat)
}
# dev.off()

# Plot estimated parameters -------------------------------------------------------------------
# go_to("plots")
# jpeg(filename = "estimated_par_exp2_bestmodels.jpg",
#      width = 16,
#      height = 8,
#      units = 'cm',
#      res = 500)
# # layout(matrix(c(1,3,2,5,1,4,2,6),ncol=2),heights = c(.05,.45,.05,.45))
# layout(matrix(c(1,2,3,4),ncol=4),widths = c(.375,.125,.375,.125))
cond_order <- c("minus","control","plus")
cexax <- cex_size(8,cex.layout)
cexlab <- cex_size(10,cex.layout)
cexmain <- 1.75
linelab <- 2.5
lwdmean <- 2.25
col_indiv_points <- rgb(.7,.7,.7,.5)
jit_size <- .1
col_hard <- rgb(247,104,161,maxColorValue = 255)
col_med <- rgb(197,27,138,maxColorValue = 255)
col_easy <- rgb(122,1,119,maxColorValue = 255)
leg_squeeze_within <- .5
leg_squeeze_between <- 1
leg_text_width <- .5*c(1,1.5,1)
# Plot alpha/beta exp2a ---------------------------------------------------------

n_alpha <- sum(subs %in% Data_alpha$sub)
par(mar=mar.parameters)
##Alpha
plot_alpha <- with(subset(param_betafixed,manip=="alpha"),aggregate(alpha,by=list(sub=sub,condition=condition),mean))
plot_alpha <- cast(plot_alpha,sub~condition)
plot_alpha <- plot_alpha[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_alpha),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_alpha),max(plot_alpha)),xaxt='n',yaxt='n');
title(ylab="Alpha",xlab="Feedback condition", line = linelab, cex.lab=cexlab)
axis(1,1:Ncond,c("Minus","Control","Plus"),cex.axis=cexax)
axis(2,seq(0,60,10),seq(0,60,10),cex.axis=cexax,las=ylab_orientation)
for(i in 1:n_alpha) lines(jitter(1:Ncond,amount=jit_size),plot_alpha[i,1:Ncond],type='b',lty=2,col=col_indiv_points,pch=19)
points(colMeans(plot_alpha),type='b',lwd=lwdmean)
error.bar(1:Ncond,colMeans(plot_alpha),colSds(plot_alpha,na.rm=T)/sqrt(n_alpha),lwd=lwdmean,length=0)
points(c(9,18,36),type='p',pch=fb_par_pch,col=fb_par_col)
legend("topleft",legend=c("Generative values","Empirical fits"),pch=c(fb_par_pch,19),
       cex=cex.legend,bty="n")
par(mar=mar.parameters.overall)
##Beta
plot_beta <- with(subset(param_betafixed,manip=="alpha"),aggregate(beta,by=list(sub=sub),mean))
plot_beta <- plot_beta$x
plot(xlab="",ylab="",plot_beta,frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.9,1.1),ylim=c(min(plot_beta),1),xaxt='n',yaxt='n');
title(ylab="Beta",xlab="", line = linelab, cex.lab=cexlab)
axis(1,1,"Overall",cex.axis=cexax, tick = F)
axis(2,seq(-2.5,1,.5),seq(-2.5,1,.5),cex.axis=cexax,las=ylab_orientation)
points(jitter(rep(1,n_alpha),amount=jit_size/3),plot_beta,col=col_indiv_points,pch=19)
points(1,mean(plot_beta),pch=16,cex=cex.datdot)
points(1,0,type='p',pch=fb_par_pch,col=fb_par_col)
error.bar(1,mean(plot_beta),sd(plot_beta,na.rm=T)/sqrt(n_alpha),lwd=lwdmean,length=0)

# Plot alpha/beta exp2b ---------------------------------------------------------
par(mar=mar.parameters)
n_beta <- sum(subs %in% Data_beta$sub)
##Beta
plot_beta <- with(subset(param_alphafixed,manip=="beta"),aggregate(beta,by=list(sub=sub,condition=condition),mean))
plot_beta <- cast(plot_beta,sub~condition)
plot_beta <- plot_beta[,cond_order] #Reorder columns to have easy -> hard
plot(xlab="",ylab="",colMeans(plot_beta),frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.8,Ncond+.2),ylim=c(min(plot_beta),max(plot_beta)),xaxt='n',yaxt='n');
title(ylab="Beta",xlab="Feedback condition", line = linelab, cex.lab=cexlab)
axis(1,1:Ncond,c("Minus","Control","Plus"),cex.axis=cexax)
axis(2,seq(-2,1,.5),seq(-2,1,.5),cex.axis=cexax,las=ylab_orientation)
for(i in 1:n_beta) lines(jitter(1:Ncond,amount=jit_size),plot_beta[i,1:Ncond],type='b',lty=2,col=col_indiv_points,pch=19)
points(colMeans(plot_beta),type='b',lwd=lwdmean)
points(c(-1,0,1),type='p',pch=fb_par_pch,col=fb_par_col)
error.bar(1:Ncond,colMeans(plot_beta),colSds(plot_beta,na.rm=T)/sqrt(n_beta),lwd=lwdmean,length=0)

par(mar=mar.parameters.overall)

##Alpha
plot_alpha <- with(subset(param_alphafixed,manip=="beta"),aggregate(alpha,by=list(sub=sub),mean))
plot_alpha <- plot_alpha$x
dist <- density(plot_alpha,cut=1)
plot(xlab="",ylab="",plot_alpha,frame=F,type='n',cex.lab=cexlab,cex.axis=cexax,
     xlim=c(.9,1.1),ylim=c(0,35),xaxt='n',yaxt='n');
title(ylab="Alpha",xlab="", line = linelab, cex.lab=cexlab)
axis(1,1,"Overall",cex.axis=cexax, tick = F)
axis(2,seq(0,35,5),seq(0,35,5),cex.axis=cexax,las=ylab_orientation)
points(jitter(rep(1,n_beta),amount=jit_size/3),plot_alpha,col=col_indiv_points,pch=19)
points(1,mean(plot_alpha),pch=16,cex=cex.datdot)
points(1,18,type='p',pch=fb_par_pch,col=fb_par_col)
error.bar(1,mean(plot_alpha),sd(plot_alpha,na.rm=T)/sqrt(n_alpha),lwd=lwdmean,length=0)

dev.off()
# Analyze parameters ------------------------------------------------------
m <- lmer(data=subset(param_betafixed,manip=="alpha"),alpha ~ condition + (1|sub))
anova(m)
emm <- emmeans(m, ~ condition)
pairs(emm)
eta_squared(anova(m),alternative="two")

m <- lmer(data=subset(param_alphafixed,manip=="beta"),beta ~ condition + (1|sub))
anova(m)
eta_squared(anova(m),alternative="two")
emm <- emmeans(m, ~ condition)
pairs(emm)


# Reviewer analysis: confidence*difficulty in error trials ----------------

go_to("plots")
# Alpha-manipulated feedback
N_temp <- length(unique(Data_alpha$sub))

cjlow <- with(subset(Data_alpha,condition=="minus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjlow) <- c('sub','difflevel','cor','cj')
cjlow_cor <- subset(cjlow,cor==1); cjlow_err <- subset(cjlow,cor==0)
cjlow_cor <- cast(cjlow_cor,sub~difflevel); cjlow_err <- cast(cjlow_err,sub~difflevel)
cjmed <- with(subset(Data_alpha,condition=="control"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjmed) <- c('sub','difflevel','cor','cj')
cjmed_cor <- subset(cjmed,cor==1); cjmed_err <- subset(cjmed,cor==0)
cjmed_cor <- cast(cjmed_cor,sub~difflevel); cjmed_err <- cast(cjmed_err,sub~difflevel)
cjhigh <- with(subset(Data_alpha,condition=="plus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjhigh) <- c('sub','difflevel','cor','cj')
cjhigh_cor <- subset(cjhigh,cor==1); cjhigh_err <- subset(cjhigh,cor==0)
cjhigh_cor <- cast(cjhigh_cor,sub~difflevel); cjhigh_err <- cast(cjhigh_err,sub~difflevel)

xminus_err <- cjlow_err[,c(2:4)];xbaseline_err <- cjmed_err[,c(2:4)];xplus_err <- cjhigh_err[,c(2:4)]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xbaseline_err <- xbaseline_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]

delta_minus <- c()
delta_baseline <- c()
delta_plus <- c()
for (i in 1:N_temp) {
  delta_minus <- c(delta_minus,xminus_err[i,3]-xminus_err[i,1])
  delta_baseline <- c(delta_baseline,xbaseline_err[i,3]-xbaseline_err[i,1])
  delta_plus <- c(delta_plus,xplus_err[i,3]-xplus_err[i,1])
}
range_delta <- range(c(delta_minus,delta_plus,delta_baseline),na.rm = T)
jpeg("diff_conf_2B.jpg",units='cm',width=8,height=8,res=300)
par(mar=c(4,4,2,0))
plot(density(delta_minus,na.rm=T),col=col_minus,lwd=2,type='l',main = "Experiment 2B",
     xlab = "Easy - Hard confidence \n in incorrect trials", ylab = 'Density', 
     ylim=c(0,.8),bty='n')
lines(density(delta_minus,na.rm=T),col=col_minus,lwd=2)
lines(density(delta_baseline,na.rm=T),col=col_baseline,lwd=2)
lines(density(delta_plus,na.rm=T),col=col_plus,lwd=2)
legend("topleft",legend=c("Minus","Control","Plus"),bty='n',title = "Feedback condition",
       fill=c(col_minus,col_baseline,col_plus),cex=.8)
dev.off()

# Beta-manipulated feedback
N_temp <- length(unique(Data_beta$sub))

cjlow <- with(subset(Data_beta,condition=="minus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjlow) <- c('sub','difflevel','cor','cj')
cjlow_cor <- subset(cjlow,cor==1); cjlow_err <- subset(cjlow,cor==0)
cjlow_cor <- cast(cjlow_cor,sub~difflevel); cjlow_err <- cast(cjlow_err,sub~difflevel)
cjmed <- with(subset(Data_beta,condition=="control"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjmed) <- c('sub','difflevel','cor','cj')
cjmed_cor <- subset(cjmed,cor==1); cjmed_err <- subset(cjmed,cor==0)
cjmed_cor <- cast(cjmed_cor,sub~difflevel); cjmed_err <- cast(cjmed_err,sub~difflevel)
cjhigh <- with(subset(Data_beta,condition=="plus"),aggregate(cj,by=list(sub,difflevel,cor),mean));
names(cjhigh) <- c('sub','difflevel','cor','cj')
cjhigh_cor <- subset(cjhigh,cor==1); cjhigh_err <- subset(cjhigh,cor==0)
cjhigh_cor <- cast(cjhigh_cor,sub~difflevel); cjhigh_err <- cast(cjhigh_err,sub~difflevel)

xminus_err <- cjlow_err[,c(2:4)];xbaseline_err <- cjmed_err[,c(2:4)];xplus_err <- cjhigh_err[,c(2:4)]
xminus_err <- xminus_err[,c("hard","medium","easy")];
xbaseline_err <- xbaseline_err[,c("hard","medium","easy")];
xplus_err <- xplus_err[,c("hard","medium","easy")]

delta_minus <- c()
delta_baseline <- c()
delta_plus <- c()
for (i in 1:N_temp) {
  delta_minus <- c(delta_minus,xminus_err[i,3]-xminus_err[i,1])
  delta_baseline <- c(delta_baseline,xbaseline_err[i,3]-xbaseline_err[i,1])
  delta_plus <- c(delta_plus,xplus_err[i,3]-xplus_err[i,1])
}
range_delta <- range(c(delta_minus,delta_plus,delta_baseline),na.rm = T)
jpeg("diff_conf_2A.jpg",units='cm',width=8,height=8,res=300)
par(mar=c(4,4,2,0))
plot(density(delta_minus,na.rm=T),col=col_minus,lwd=2,type='l',main = "Experiment 2A",
     xlab = "Easy - Hard confidence \n in incorrect trials", ylab = 'Density', 
     ylim=c(0,.8),bty='n')
lines(density(delta_minus,na.rm=T),col=col_minus,lwd=2)
lines(density(delta_baseline,na.rm=T),col=col_baseline,lwd=2)
lines(density(delta_plus,na.rm=T),col=col_plus,lwd=2)
legend("topleft",legend=c("Minus","Control","Plus"),bty='n',title = "Feedback condition",
       fill=c(col_minus,col_baseline,col_plus),cex=.8)
dev.off()

if (stat_tests) {
  m <- lmer(data = subset(Data_alpha,cor==0), cj ~ difflevel*condition + (difflevel|sub))
  m2 <- lmer(data = subset(Data_alpha,cor==0), cj ~ difflevel*condition + (condition|sub))
  m3 <- lmer(data = subset(Data_alpha,cor==0), cj ~ difflevel*condition + (difflevel+condition|sub))
  anova(m3,m2)
  anova(m3,m)
  anova(m3)
  
  mb <- lmer(data = subset(Data_beta,cor==0), cj ~ difflevel*condition + (condition|sub),
             REML = F,control = lmerControl(optimizer='bobyqa'))
  anova(mb)
}


# Change of mind exploration ----------------------------------------------
# First, we check if indeed the beta feedback manipulation led to varying proportions
# of trials with feedback < 50%
train_beta$fb_com <- as.numeric(train_beta$fb<50)
fb_com <- with(train_beta,aggregate(fb_com,by=list(sub,condition),mean))
names(fb_com) <- c("sub","condition","fb")
fb_com <- cast(fb_com,sub~condition)
colMeans(fb_com)


# Merge the BICs with the main dataframe and add a change of mind column
bfree_bic <- with(subset(param,fixed_parameter=="alpha"),aggregate(bic,list(sub=sub),mean))
names(bfree_bic) <- c("sub","bfree_bic")
afree_bic <- with(subset(param,fixed_parameter=="beta"),aggregate(bic,list(sub=sub),mean))
names(afree_bic) <- c("sub","afree_bic")

# Add BICs to the main dataframe
Data <- merge(Data,afree_bic)
Data <- merge(Data,bfree_bic)

Data$bic_diff <- Data$afree_bic - Data$bfree_bic # Difference in BICs
Data$com <- as.numeric(Data$cj<4) # Change of mind, defined as "guess error" or lower

Data_beta <- subset(Data,manip=="beta")

## Correlation between change of mind and BIC difference
com_bic <- with(Data_beta,aggregate(list(bic_diff,com),by=list(sub),mean))
names(com_bic) <- c("sub","bic_diff","com")

cor.test(com_bic$bic_diff,com_bic$com,method="pearson")
jpeg("com-bic_2b.jpg",width=8,height=8,units="cm",res=300)
par(mar=c(4,4,2,1))
plot(com_bic$bic_diff~com_bic$com,ylab=expression(paste(alpha,"-free - ",beta,"-free BIC")),
     xlab = "Proportion changes of mind", bty='n',pch=1,cex=.8,
     main=paste("r = ",round(cor(com_bic$bic_diff,com_bic$com),3)))
abline(lm(com_bic$bic_diff~com_bic$com),lwd=2)
dev.off()
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1)

## Correlation between the feedback manipulation effect on change of mind and BIC difference
com_cond <- with(Data_beta,aggregate(com,by=list(sub,condition),mean))
names(com_cond) <- c("sub","condition","com")
com_cond <- cast(com_cond, sub~condition)
com_cond$diff <- com_cond$plus - com_cond$minus
mean(com_cond$diff)
sd(com_cond$diff)

jpeg("change_com-bic_2b.jpg",width=8,height=8,units="cm",res=300)
par(mar=c(4,4,2,1))
plot(com_bic$bic_diff~com_cond$diff,ylab=expression(paste(alpha,"-free - ",beta,"-free BIC")),
     xlab = "Plus - Minus Feedback \n proportion changes of mind", bty='n',pch=1,cex=.8,
     main=paste("r = ",round(cor(com_bic$bic_diff,com_cond$diff),3)))
abline(lm(com_bic$bic_diff~com_cond$diff),lwd=2)
dev.off()
cor.test(com_bic$bic_diff,com_cond$diff)

# Check the beta-model-predicted proportions of changes of mind
Simuls_afix_beta$com <- as.numeric(Simuls_afix_beta$cj<4)
beta_com <- with(Simuls_afix_beta,aggregate(com,by=list(sub,condition),mean))
names(beta_com) <- c("sub","condition","com")
beta_com <- cast(beta_com,sub~condition)
beta_com$diff <- beta_com$plus - beta_com$minus
mean(beta_com$diff)
sd(beta_com$diff)

# T-tests
t.test(com_cond$minus,com_cond$plus, paired = T) # Feedback manipulation effect on changes of mind

t.test(beta_com$minus,beta_com$plus, paired = T) # Beta model predicted changes of mind

# Difference in beta between plus and minus conditions
t.test(with(subset(param_alphafixed,manip=="beta"&condition=="minus"),aggregate(beta,list(sub),mean))$x,
       with(subset(param_alphafixed,manip=="beta"&condition=="plus"),aggregate(beta,list(sub),mean))$x,
       paired = T)

# Participants with high proportions of changes of mind --------------

com_sub <- with(Data,aggregate(com,by=list(sub=sub),mean))
hist(com_sub$x) # few participants with high proportions of changes of mind

# Let's investigate further these participants
high_com <- com_sub[com_sub$x>.40,"sub"]
Data_highsub <- subset(Data,sub %in% high_com) # 4 participants

with(Data_highsub,aggregate(cor,by=list(sub),mean)) # One has low accuracy

# To check that these participants still understood properly the confidence scale,
# we look at the correlation between confidence judgments and accuracy
cj_acc <- cast(Data_highsub,value="cor",sub~cj,fun.aggregate = mean)
cor_cj_acc <- function(X) {
  return(cor(X,1:length(X)))
}
cor_highsub <- apply(cj_acc[,2:7],1,cor_cj_acc) # All very high positive correlations except for one

# Just to be sure, we redo the statistical analyses without the one odd participant
Data_beta_com_filtered <- subset(Data_beta, sub != cj_acc[cor_highsub<.8,"sub"])

# RT
m.cond <- lmer(rt~condition*difflevel + (condition|sub),data=Data_beta_com_filtered)
anova(m.cond) #Results

# Accuracy
m.cond <- glmer(cor~condition*difflevel + (condition|sub),data=Data_beta_com_filtered,family=binomial)
Anova(m.cond)

# Confidence
minteraction <- lmer(cj ~ condition*cor*difflevel + (cor*condition|sub),
                     data = Data_beta_com_filtered,REML = F,control=lmerControl(optimize="bobyqa"));
anova(minteraction)

mcor <- lmer(cj ~ condition*difflevel + (condition|sub),data = subset(Data_beta_com_filtered,cor==1)); 
merr <- lmer(cj ~ condition*difflevel + (condition|sub),data = subset(Data_beta_com_filtered,cor==0)); 
anova(mcor) #Results
anova(merr) #Results

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

Data_alpha$rt_bin <- cut(Data_alpha$rt,breaks=c(seq(0,1.5,.5),max(Data_alpha$rt)))
test <- with(Data_alpha,aggregate(cj,by=list(sub,rt_bin),mean))
names(test) <- c("sub","rt_bin","cj")
test <- cast(test,sub~rt_bin)
test <- test[,-1]


Simuls_bfix_alpha$rt_bin <- cut(Simuls_bfix_alpha$rt,breaks=c(seq(0,1.5,.5),max(Simuls_bfix_alpha$rt)))
test_sim <- with(Simuls_bfix_alpha,aggregate(cj,by=list(sub,rt_bin),mean))
names(test_sim) <- c("sub","rt_bin","cj")
test_sim <- cast(test_sim,sub~rt_bin)
test_sim <- test_sim[,-1]


n <- length(test)

na_count_dat <- sapply(test, function(y) sum(length(which(!is.na(y)))))
na_count_sim <- sapply(test_sim, function(y) sum(length(which(!is.na(y)))))

cex.datdot <- 1

go_to("plots")
tiff("cj_rt_dist_exp2a.tiff",width=8,height=8,units="cm",res=600,compression="lzw")
par(mar=mar.plot)

stripchart(test,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(3,6), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = "Experiment 2A")
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


Data_beta$rt_bin <- cut(Data_beta$rt,breaks=c(seq(0,1.5,.5),max(Data_beta$rt)))
test <- with(Data_beta,aggregate(cj,by=list(sub,rt_bin),mean))
names(test) <- c("sub","rt_bin","cj")
test <- cast(test,sub~rt_bin)
test <- test[,-1]


Simuls_afix_beta$rt_bin <- cut(Simuls_afix_beta$rt,breaks=c(seq(0,1.5,.5),max(Simuls_afix_beta$rt)))
test_sim <- with(Simuls_afix_beta,aggregate(cj,by=list(sub,rt_bin),mean))
names(test_sim) <- c("sub","rt_bin","cj")
test_sim <- cast(test_sim,sub~rt_bin)
test_sim <- test_sim[,-1]


n <- length(test)

na_count_dat <- sapply(test, function(y) sum(length(which(!is.na(y)))))
na_count_sim <- sapply(test_sim, function(y) sum(length(which(!is.na(y)))))

cex.datdot <- 1

go_to("plots")
tiff("cj_rt_dist_exp2b.tiff",width=8,height=8,units="cm",res=600,compression="lzw")
par(mar=mar.plot)

stripchart(test,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(3,6), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = "Experiment 2B")
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

## Correlations
# Behavior
correlation_a <- c()
for (subj in unique(Data_alpha$sub)) {
  correlation_a <- c(correlation_a,cor(Data_alpha[Data_alpha$sub==subj,]$cj,Data_alpha[Data_alpha$sub==subj,]$rt,method = "spearman"))
}
t.test(correlation_a)

correlation_b <- c()
for (subj in unique(Data_beta$sub)) {
  correlation_b <- c(correlation_b,cor(Data_beta[Data_beta$sub==subj,]$cj,Data_beta[Data_beta$sub==subj,]$rt,method = "spearman"))
}
t.test(correlation_b)

# Model predictions
correlation_a_sim <- c()
for (subj in unique(Simuls_bfix_alpha$sub)) {
  correlation_a_sim <- c(correlation_a_sim,cor(Simuls_bfix_alpha[Simuls_bfix_alpha$sub==subj,]$cj,Simuls_bfix_alpha[Simuls_bfix_alpha$sub==subj,]$rt,method = "spearman"))
}
t.test(correlation_a_sim)

correlation_b_sim <- c()
for (subj in unique(Simuls_afix_beta$sub)) {
  correlation_b_sim <- c(correlation_b_sim,cor(Simuls_afix_beta[Simuls_afix_beta$sub==subj,]$cj,Simuls_afix_beta[Simuls_afix_beta$sub==subj,]$rt,method = "spearman"))
}
t.test(correlation_b_sim)

# Other Confidence ~ RT analysis ------------------------------------------
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
mar.plot <- c(4,4,0,1)+.1

### Experiment 2A
Data_alpha$rt_bin <- cut(Data_alpha$rt,breaks=c(seq(0,2,.5),max(Data_alpha$rt)))
test <- with(Data_alpha,aggregate(cj,by=list(sub,cor,rt_bin),mean))
names(test) <- c("sub","cor","rt_bin","cj")
test <- cast(test,sub+cor~rt_bin)
test_cor <- test[test$cor==1,c(-1,-2)]
test_err <- test[test$cor==0,c(-1,-2)]


Simuls_bfix_alpha$rt_bin <- cut(Simuls_bfix_alpha$rt,breaks=c(seq(0,2,.5),max(Simuls_bfix_alpha$rt)))
test_sim <- with(Simuls_bfix_alpha,aggregate(cj,by=list(sub,cor,rt_bin),mean))
names(test_sim) <- c("sub","cor","rt_bin","cj")
test_sim <- cast(test_sim,sub+cor~rt_bin)
test_sim_cor <- test_sim[test_sim$cor==1,c(-1,-2)]
test_sim_err <- test_sim[test_sim$cor==0,c(-1,-2)]


n <- length(test_cor)

na_count_dat_cor <- sapply(test_cor, function(y) sum(length(which(!is.na(y)))))
na_count_dat_err <- sapply(test_err, function(y) sum(length(which(!is.na(y)))))
na_count_sim_cor <- sapply(test_sim_cor, function(y) sum(length(which(!is.na(y)))))
na_count_sim_err <- sapply(test_sim_err, function(y) sum(length(which(!is.na(y)))))

cex.datdot <- 1

go_to("plots")
tiff("cj_rt_dist_exp2a.tiff",width=16,height=8,units="cm",res=600,compression="lzw")
layout(matrix(c(1,2,1,3),ncol=2),heights=c(.1,1))
par(mar=c(0,0,0,0))
plot.new()
title(main="Experiment 2A",cex.main=cex.lab,line=-1)
par(mar=mar.plot)

stripchart(test_cor,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(3,6), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = "")
mtext("Confidence",2,at=4.5,line=2.5,cex=cex.lab);
mtext("RT bin (s)",1,at=n/2,line=2.5,cex=cex.lab);
axis(1,at=0:(n-1),labels=c(names(test_cor)[-5],expression(paste("(2,",infinity,"]"))), cex.axis=cex.ax);
axis(2, seq(3,6,.5), cex.axis=cex.ax, las = 2)
# Correct
means <- sapply(test_cor, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=1)
error.bar(0:(n-1),means,colSds(as.matrix(test_cor),na.rm=T)/sqrt(na_count_dat_cor),lwd=lwd.dat)
# Error
means <- sapply(test_err, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=2)
error.bar(0:(n-1),means,colSds(as.matrix(test_err),na.rm=T)/sqrt(na_count_dat_err),lwd=lwd.dat)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim_cor,na.rm=T) + (colSds(as.matrix(test_sim_cor),na.rm=T)/sqrt(na_count_sim_cor)),
          (colMeans(test_sim_cor,na.rm=T) - colSds(as.matrix(test_sim_cor),na.rm=T)/sqrt(na_count_sim_cor))[n:1]),
        border=F,col = rgb(0,0,0,.5))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim_err,na.rm=T) + (colSds(as.matrix(test_sim_err),na.rm=T)/sqrt(na_count_sim_err)),
          (colMeans(test_sim_err,na.rm=T) - colSds(as.matrix(test_sim_err),na.rm=T)/sqrt(na_count_sim_err))[n:1]),
        border=F,col = rgb(0,0,0,.5))
# Add trial accuracy as a legend
legend("topright",c("Empirical correct","Empirical incorrect","Model predictions"),lty=c(1,2,NA),
       lwd=1.5*ratio,col=c("black","black",rgb(0,0,0,.5)),border=NA,pch=c(NA,NA,15),pt.cex=2,
       cex=cex.leg,bty="n")
# dev.off()

# Let's plot the distribution of RT for corrects and errors
go_to("plots")

breaks <- seq(0,5,.25)
Simuls_bfix_alpha <- subset(Simuls_bfix_alpha,rt<4)
# tiff("rt_dist_exp1.tiff",width=8,height=8,units="cm",res=600,compression="lzw")

par(mar=mar.plot)

dist_plot <- hist(Data_alpha$rt[Data_alpha$cor==1],breaks=breaks,col=rgb(0,1,0,.5),border=NA,main="",
                  xlab="",ylab="",cex.lab=cex.lab,cex.axis=cex.ax,xaxt='n')
axis(1,at=seq(0,5,.5),labels=seq(0,5,.5),cex.axis=cex.ax)
mtext("Number of trials",2,at=max(dist_plot$counts)/2,line=2.5,cex=cex.lab);
mtext("RT (s)",1,at=2,line=2.5,cex=cex.lab);
hist(Data_alpha$rt[Data_alpha$cor==0],breaks=breaks,col=rgb(1,0,0,.5),border=NA,add=T)
# Add model predictions as lines
pred_correct <- hist(Simuls_bfix_alpha$rt[Simuls_bfix_alpha$cor==1],breaks=breaks,plot=F)
lines(breaks[-1]-.125,pred_correct$counts,col="darkgreen",lwd=2,lty=1)
pred_error <- hist(Simuls_bfix_alpha$rt[Simuls_bfix_alpha$cor==0],breaks=breaks,plot=F)
lines(breaks[-1]-.125,pred_error$counts,col="darkred",lwd=2,lty=1)
# Legend
legend("topright",c("Empirical correct","Empirical incorrect","Model correct", "Model incorrect"),
       lty=c(NA,NA,1,1),lwd=1.5*ratio,col=c(rgb(0,1,0,.5),rgb(1,0,0,.5),"darkgreen","darkred"),border=NA,
       cex=cex.leg,bty="n",pch=c(15,15,NA,NA),pt.cex=2)
dev.off()

### Experiment 2B

Data_beta$rt_bin <- cut(Data_beta$rt,breaks=c(seq(0,2,.5),max(Data_beta$rt)))
test <- with(Data_beta,aggregate(cj,by=list(sub,cor,rt_bin),mean))
names(test) <- c("sub","cor","rt_bin","cj")
test <- cast(test,sub+cor~rt_bin)
test_cor <- test[test$cor==1,c(-1,-2)]
test_err <- test[test$cor==0,c(-1,-2)]


Simuls_afix_beta$rt_bin <- cut(Simuls_afix_beta$rt,breaks=c(seq(0,2,.5),max(Simuls_afix_beta$rt)))
test_sim <- with(Simuls_afix_beta,aggregate(cj,by=list(sub,cor,rt_bin),mean))
names(test_sim) <- c("sub","cor","rt_bin","cj")
test_sim <- cast(test_sim,sub+cor~rt_bin)
test_sim_cor <- test_sim[test_sim$cor==1,c(-1,-2)]
test_sim_err <- test_sim[test_sim$cor==0,c(-1,-2)]


n <- length(test_cor)

na_count_dat_cor <- sapply(test_cor, function(y) sum(length(which(!is.na(y)))))
na_count_dat_err <- sapply(test_err, function(y) sum(length(which(!is.na(y)))))
na_count_sim_cor <- sapply(test_sim_cor, function(y) sum(length(which(!is.na(y)))))
na_count_sim_err <- sapply(test_sim_err, function(y) sum(length(which(!is.na(y)))))

cex.datdot <- 1

go_to("plots")
tiff("cj_rt_dist_exp2b.tiff",width=16,height=8,units="cm",res=600,compression="lzw")
layout(matrix(c(1,2,1,3),ncol=2),heights=c(.1,1))
par(mar=c(0,0,0,0))
plot.new()
title(main="Experiment 2B",cex.main=cex.lab,line=-1)
par(mar=mar.plot)

stripchart(test_cor,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(3,6), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = "")
mtext("Confidence",2,at=4.5,line=2.5,cex=cex.lab);
mtext("RT bin (s)",1,at=n/2,line=2.5,cex=cex.lab);
axis(1,at=0:(n-1),labels=c(names(test_cor)[-5],expression(paste("(2,",infinity,"]"))), cex.axis=cex.ax);
axis(2, seq(3,6,.5), cex.axis=cex.ax, las = 2)
# Correct
means <- sapply(test_cor, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=1)
error.bar(0:(n-1),means,colSds(as.matrix(test_cor),na.rm=T)/sqrt(na_count_dat_cor),lwd=lwd.dat)
# Error
means <- sapply(test_err, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=2)
error.bar(0:(n-1),means,colSds(as.matrix(test_err),na.rm=T)/sqrt(na_count_dat_err),lwd=lwd.dat)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim_cor,na.rm=T) + (colSds(as.matrix(test_sim_cor),na.rm=T)/sqrt(na_count_sim_cor)),
          (colMeans(test_sim_cor,na.rm=T) - colSds(as.matrix(test_sim_cor),na.rm=T)/sqrt(na_count_sim_cor))[n:1]),
        border=F,col = rgb(0,0,0,.5))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim_err,na.rm=T) + (colSds(as.matrix(test_sim_err),na.rm=T)/sqrt(na_count_sim_err)),
          (colMeans(test_sim_err,na.rm=T) - colSds(as.matrix(test_sim_err),na.rm=T)/sqrt(na_count_sim_err))[n:1]),
        border=F,col = rgb(0,0,0,.5))
# Add trial accuracy as a legend
legend("topright",c("Empirical correct","Empirical incorrect","Model predictions"),lty=c(1,2,NA),
       lwd=1.5*ratio,col=c("black","black",rgb(0,0,0,.5)),border=NA,pch=c(NA,NA,15),pt.cex=2,
       cex=cex.leg,bty="n")
# dev.off()

# Let's plot the distribution of RT for corrects and errors
go_to("plots")

breaks <- seq(0,5,.25)
Simuls_afix_beta <- subset(Simuls_afix_beta,rt<4)
# tiff("rt_dist_exp1.tiff",width=8,height=8,units="cm",res=600,compression="lzw")

par(mar=mar.plot)

dist_plot <- hist(Data_beta$rt[Data_beta$cor==1],breaks=breaks,col=rgb(0,1,0,.5),border=NA,main="",
                  xlab="",ylab="",cex.lab=cex.lab,cex.axis=cex.ax,xaxt='n')
axis(1,at=seq(0,5,.5),labels=seq(0,5,.5),cex.axis=cex.ax)
mtext("Number of trials",2,at=max(dist_plot$counts)/2,line=2.5,cex=cex.lab);
mtext("RT (s)",1,at=2,line=2.5,cex=cex.lab);
hist(Data_beta$rt[Data_beta$cor==0],breaks=breaks,col=rgb(1,0,0,.5),border=NA,add=T)
# Add model predictions as lines
pred_correct <- hist(Simuls_afix_beta$rt[Simuls_afix_beta$cor==1],breaks=breaks,plot=F)
lines(breaks[-1]-.125,pred_correct$counts,col="darkgreen",lwd=2,lty=1)
pred_error <- hist(Simuls_afix_beta$rt[Simuls_afix_beta$cor==0],breaks=breaks,plot=F)
lines(breaks[-1]-.125,pred_error$counts,col="darkred",lwd=2,lty=1)
# Legend
legend("topright",c("Empirical correct","Empirical incorrect","Model correct", "Model incorrect"),
       lty=c(NA,NA,1,1),lwd=1.5*ratio,col=c(rgb(0,1,0,.5),rgb(1,0,0,.5),"darkgreen","darkred"),border=NA,
       cex=cex.leg,bty="n",pch=c(15,15,NA,NA),pt.cex=2)
dev.off()

# Test the relationship between confidence and RT
control <- lmerControl(optimizer="bobyqa")
# Experiment 2A
m <- lmer(data = Data_alpha, cj ~ rt * cor + (rt * cor | sub), REML = F, control = control)
anova(m)

grouped_alpha <- split(Data_alpha[,c('cj','rt','cor')],list(Data_alpha$sub))
correlations_alpha_cor <- sapply(grouped_alpha,function(x) cor(subset(x,cor==1)$cj,subset(x,cor==1)$rt))
correlations_alpha_cor <- correlations_alpha_cor[!is.na(correlations_alpha_cor)]
correlations_alpha_err <- sapply(grouped_alpha,function(x) cor(subset(x,cor==0)$cj,subset(x,cor==0)$rt))
correlations_alpha_err <- correlations_alpha_err[!is.na(correlations_alpha_err)]
grouped_beta <- split(Data_beta[,c('cj','rt','cor')],list(Data_beta$sub))
correlations_beta_cor <- sapply(grouped_beta,function(x) cor(subset(x,cor==1)$cj,subset(x,cor==1)$rt))
correlations_beta_cor <- correlations_beta_cor[!is.na(correlations_beta_cor)]
correlations_beta_err <- sapply(grouped_beta,function(x) cor(subset(x,cor==0)$cj,subset(x,cor==0)$rt))
correlations_beta_err <- correlations_beta_err[!is.na(correlations_beta_err)]
cor_range <- range(c(correlations_alpha_cor,correlations_beta_cor,correlations_alpha_err,correlations_beta_err))

breaks <- seq(-1,1,.1)
# Plot for Experiment 2A
hist(correlations_alpha_cor,breaks=breaks,main="",xlab="Empirical - Model correlation",
     ylab="Nb subjects",xlim=cor_range,col=rgb(0,1,0,.5))
abline(v=mean(correlations_alpha_cor),lwd=2,col="darkgreen")
# Add the distribution for errors
hist(correlations_alpha_err,breaks=breaks,add=T,col=rgb(1,0,0,.5))
abline(v=mean(correlations_alpha_err),lwd=2,col="darkred")

# Plot for Experiment 2B
hist(correlations_beta_cor,breaks=breaks,main="",xlab="Empirical - Model correlation",
     ylab="Nb subjects",xlim=cor_range,col=rgb(0,1,0,.5))
abline(v=mean(correlations_beta_cor),lwd=2,col="darkgreen")
# Add the distribution for errors
hist(correlations_beta_err,breaks=breaks,add=T,col=rgb(1,0,0,.5))
abline(v=mean(correlations_beta_err),lwd=2,col="darkred")

cor.test(subset(Data_alpha,cor==1)$cj,subset(Data_alpha,cor==1)$rt,method="spearman")
cor.test(subset(Data_alpha,cor==0)$cj,subset(Data_alpha,cor==0)$rt,method="spearman")
# Experiment 2B
m <- lmer(data = Data_beta, cj ~ rt * cor + (rt * cor | sub), REML = F, control = control)
anova(m)
cor.test(subset(Data_beta,cor==1)$cj,subset(Data_beta,cor==1)$rt,method="spearman")
cor.test(subset(Data_beta,cor==0)$cj,subset(Data_beta,cor==0)$rt,method="spearman")


Data_beta$rt_bin <- cut(Data_beta$rt,breaks=c(seq(0,2,.5),max(Data_beta$rt)))
test <- with(Data_beta,aggregate(cj,by=list(sub,cor,rt_bin),mean))
names(test) <- c("sub","cor","rt_bin","cj")
test <- cast(test,sub+cor~rt_bin)
test_cor <- test[test$cor==1,c(-1,-2)]
test_err <- test[test$cor==0,c(-1,-2)]


Simuls_afix_beta$rt_bin <- cut(Simuls_afix_beta$rt,breaks=c(seq(0,2,.5),max(Simuls_afix_beta$rt)))
test_sim <- with(Simuls_afix_beta,aggregate(cj,by=list(sub,cor,rt_bin),mean))
names(test_sim) <- c("sub","cor","rt_bin","cj")
test_sim <- cast(test_sim,sub+cor~rt_bin)
test_sim_cor <- test_sim[test_sim$cor==1,c(-1,-2)]
test_sim_err <- test_sim[test_sim$cor==0,c(-1,-2)]


n <- length(test_cor)

na_count_dat_cor <- sapply(test_cor, function(y) sum(length(which(!is.na(y)))))
na_count_dat_err <- sapply(test_err, function(y) sum(length(which(!is.na(y)))))
na_count_sim_cor <- sapply(test_sim_cor, function(y) sum(length(which(!is.na(y)))))
na_count_sim_err <- sapply(test_sim_err, function(y) sum(length(which(!is.na(y)))))

# Adapt dot size to number of trials
max.cex.datdot <- 2
bin_freq <- table(Data_beta$rt_bin)/nrow(Data_beta)
cex.datdot <- bin_freq[as.character(names(test))] * max.cex.datdot/max(bin_freq)
cex.datdot <- 1

go_to("plots")
tiff("cj_rt_dist_exp2b.tiff",width=16,height=8,units="cm",res=600,compression="lzw")
layout(matrix(c(1,2,1,3),ncol=2),heights=c(.1,1))
par(mar=c(0,0,0,0))
plot.new()
title(main="Experiment 2B",cex.main=cex.lab,line=-1)
par(mar=mar.plot)

stripchart(test_cor,vertical = TRUE, col="white",frame=F,xaxt='n',
           ylim=c(3,6), xlim=c(-.05,n-1),
           yaxt = 'n',xlab="",ylab = "",
           main = "")
mtext("Confidence",2,at=4.5,line=2.5,cex=cex.lab);
mtext("RT bin (s)",1,at=n/2,line=2.5,cex=cex.lab);
axis(1,at=0:(n-1),labels=c(names(test_cor)[-5],expression(paste("(2,",infinity,"]"))), cex.axis=cex.ax);
axis(2, seq(3,6,.5), cex.axis=cex.ax, las = 2)
# Correct
means <- sapply(test_cor, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=1)
error.bar(0:(n-1),means,colSds(as.matrix(test_cor),na.rm=T)/sqrt(na_count_dat_cor),lwd=lwd.dat)
# Error
means <- sapply(test_err, mean, na.rm = T)
lines(0:(n-1),means,type='b',pch=16,cex=cex.datdot,lwd=lwd.dat,lty=2)
error.bar(0:(n-1),means,colSds(as.matrix(test_err),na.rm=T)/sqrt(na_count_dat_err),lwd=lwd.dat)
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim_cor,na.rm=T) + (colSds(as.matrix(test_sim_cor),na.rm=T)/sqrt(na_count_sim_cor)),
          (colMeans(test_sim_cor,na.rm=T) - colSds(as.matrix(test_sim_cor),na.rm=T)/sqrt(na_count_sim_cor))[n:1]),
        border=F,col = rgb(0,0,0,.5))
polygon(c(0:(n-1),(n-1):0),
        c(colMeans(test_sim_err,na.rm=T) + (colSds(as.matrix(test_sim_err),na.rm=T)/sqrt(na_count_sim_err)),
          (colMeans(test_sim_err,na.rm=T) - colSds(as.matrix(test_sim_err),na.rm=T)/sqrt(na_count_sim_err))[n:1]),
        border=F,col = rgb(0,0,0,.5))
# Add trial accuracy as a legend
legend("topright",c("Empirical correct","Empirical incorrect","Model predictions"),lty=c(1,2,NA),
       lwd=1.5*ratio,col=c("black","black",rgb(0,0,0,.5)),border=NA,pch=c(NA,NA,15),pt.cex=2,
       cex=cex.leg,bty="n")
# dev.off()

# Let's plot the distribution of RT for corrects and errors
go_to("plots")

breaks <- seq(0,5,.25)
Simuls_afix_beta <- subset(Simuls_afix_beta,rt<4)
# tiff("rt_dist_exp1.tiff",width=8,height=8,units="cm",res=600,compression="lzw")

par(mar=mar.plot)

dist_plot <- hist(Data_beta$rt[Data_beta$cor==1],breaks=breaks,col=rgb(0,1,0,.5),border=NA,main="",
                  xlab="",ylab="",cex.lab=cex.lab,cex.axis=cex.ax,xaxt='n')
axis(1,at=seq(0,5,.5),labels=seq(0,5,.5),cex.axis=cex.ax)
mtext("Number of trials",2,at=max(dist_plot$counts)/2,line=2.5,cex=cex.lab);
mtext("RT (s)",1,at=2,line=2.5,cex=cex.lab);
hist(Data_beta$rt[Data_beta$cor==0],breaks=breaks,col=rgb(1,0,0,.5),border=NA,add=T)
# Add model predictions as lines
pred_correct <- hist(Simuls_afix_beta$rt[Simuls_afix_beta$cor==1],breaks=breaks,plot=F)
lines(breaks[-1]-.125,pred_correct$counts,col="darkgreen",lwd=2,lty=1)
pred_error <- hist(Simuls_afix_beta$rt[Simuls_afix_beta$cor==0],breaks=breaks,plot=F)
lines(breaks[-1]-.125,pred_error$counts,col="darkred",lwd=2,lty=1)
# Legend
legend("topright",c("Empirical correct","Empirical incorrect","Model correct", "Model incorrect"),
       lty=c(NA,NA,1,1),lwd=1.5*ratio,col=c(rgb(0,1,0,.5),rgb(1,0,0,.5),"darkgreen","darkred"),border=NA,
       cex=cex.leg,bty="n",pch=c(15,15,NA,NA),pt.cex=2)
dev.off()

