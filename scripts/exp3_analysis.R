#' @title Analysis of the data from experiment 3
#' @author Pierre Le Denmat
#' @description Analysis of the data from experiment 3 in the context of the research article: 
#' "A low-dimensional approximation of optimal confidence".
#' @details The script is divided into several sections:
#' - Preprocessing: The data is loaded and preprocessed.
#' - Behaviour analysis: The behavioural data is analyzed using (generalized) linear mixed models.

# First, let's clear the workspace
rm(list=ls())

# Set working directory to this script's location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Now, let's load relevant packages
library(readbulk)
library(reshape)
library(lmerTest)
library(lattice)
library(car)
library(emmeans)
library(Rcpp)
library(ggplot2)

plot_all <- FALSE
# Functions ---------------------------------------------------------------

max_count <- function(data){
  return(max(table(data)))
}


# data load and preprocessing ---------------------------------------------
# Load data
data <- readbulk::read_bulk("../data/exp3/")

# Quick look at the data
head(data)
summary(data)

# Remove rows with missing values in the cor column
data <- data[!is.na(data$cor),]

# Remove empty columns
data <- data[,colSums(is.na(data))<nrow(data)]

table(table(data$sub))

# Remove participants who did not come to session 2
data <- data[!(data$sub %in% c(2,8,21,24)),]

# 2 participants entered their Sona ID as subject ID for session 2, let's find who it was
check <- as.numeric(names(table(data$sub)[table(data$sub)<2400]))
with(subset(data,sub %in% check),aggregate(cbind(age,gender,handedness),by=list(sub),unique))
# 57784 -> 30, 58333 -> 26
data[data$sub == 57784,"sub"] <- 30
data[data$sub == 58333,"sub"] <- 26

table(table(data$sub))
# Filter out participants who did not complete the experiment 
complete_session <- 1200 # 2 tasks * (4 blocks * 60 training trials + 6 blocks * 60 test trials)
completed_sessions <- apply(with(data, table(sub, session)) == complete_session, 1, all)
subjects <- names(completed_sessions)[completed_sessions]
subjects <- c(subjects,1) # sub 1 had more trials in session 1
data <- subset(data, sub %in% subjects)
Nsub <- length(subjects)

# Recode the response column to 0s and 1s
data$resp <- ifelse(data$resp == "['c']", 0, 1)

# Diagnostic plot per participant + chance performance testing
exclusion <- c()
tasks <- unique(data$task)
par(mfrow=c(2,2))
for(i in 1:Nsub){
  for (t in tasks) {
    tempDat <- subset(data,sub==subjects[i]&task==t)
    acc_block <- with(tempDat,aggregate(cor,by=list(block=block),mean))
    bias_block <- with(tempDat,aggregate(resp,by=list(block=block),mean))
    if (plot_all) {
      par(mar=c(3,4,2,0))
      plot(acc_block,ylab="Acc (.) and bias (x)",frame=F,ylim=c(0,1));abline(h=.5,lty=2,col="grey")
      points(bias_block,pch=4)
      plot(tempDat$rt,frame=F,col=c("black"),main=paste('subject',subjects[i],t),ylab="RT")
      plot(tempDat$cj,frame=F,col=c("black"),ylim=c(1,6),ylab="conf")
      plot(tempDat$RTconf,frame=F,col=c("black"),ylab="RT_conf")
    }
    test <- binom.test(length(tempDat$cor[tempDat$cor==1]),n=length(tempDat$cor),alternative = "greater")
    print(paste("In t0, sub",subjects[i],"p =", round(test$p.value,3),"compared to chance"))
    if (test$p.value > .05) {
      exclusion <- c(exclusion,subjects[i])
    }
    
  }
}
par(mfrow=c(1,1))

# Remove participants who responded the same confidence level for more than 90% of trials
max_conf_proportion <- with(subset(data,running=='main'),aggregate(cj,by=list(sub=sub),function(x) max_count(x)/length(x)))
exclusion <- c(exclusion, max_conf_proportion[max_conf_proportion$x>.9,1])

age <- with(subset(data,session==1),aggregate(age,list(sub=sub),unique))
summary(age$x)
sd(age$x)
table(data$gender)/2400

data <- subset(data, !(sub %in% exclusion))


# Trim RTs to remove outliers
data <- subset(data, rt < 5 & rt > .1 & RTconf < 5)

# Sort data by subject, session, block, and trial
data <- data[order(data$sub,data$session,data$block,data$trial),]

# Remove training data for analysis
data_training <- subset(data, running == "training")
data <- subset(data, running == "main")

# Change task names
data$task <- as.character(data$task)
data[data$task=='xotask','task'] <- 'LetterTask'
data[data$task=='dotcolour','task'] <- 'ColorTask'
data$task <- factor(data$task,levels = c('LetterTask','ColorTask'))

data$fbtask <- as.character(data$fbtask)
data[data$fbtask=='xotask','fbtask'] <- 'LetterTask'
data[data$fbtask=='dotcolour','fbtask'] <- 'ColorTask'
data$fbtask <- factor(data$fbtask,levels = c('LetterTask','ColorTask'))

# Rename feedback labels for plotting
data[data$feedback=="negativefb","feedback"] <- "Low"
data[data$feedback=="positivefb","feedback"] <- "High"
data$feedback <- factor(data$feedback,levels = c('High','Low'))

data$resp[data$resp==0] <- -1
data$feedback <- as.factor(data$feedback)
data$difflevel <- as.factor(data$difflevel)
# data$sub <- as.factor(data$sub)
data$task <- as.factor(data$task)
data$rt2 <- data$rt + data$RTconf
subjects <- unique(data$sub); Nsub <- length(subjects)
tasks <- unique(data$task); Ntask <- length(tasks)
feedbacks <- unique(data$feedback); Nfeedback <- length(feedbacks)
difflevels <- sort(unique(data$difflevel)); Ndifflevel <- length(difflevels)

data$transfer <- 'Close'
data[data$fbtask == data$task, 'transfer'] <- 'Direct'


# Behavioral analysis -----------------------------------------------------
control <- lmerControl(optimizer="bobyqa")
glmercontrol <- glmerControl(optimizer="bobyqa")

# Set contrast coding
options(contrasts = c("contr.sum","contr.poly"))

# Accuracy
accuracy.same <- glmer(cor ~ feedback*difflevel + (feedback|sub), data=subset(data,transfer=="Direct"), family=binomial, control=glmercontrol)
Anova(accuracy.same)
accuracy.transfer <- glmer(cor ~ feedback*difflevel + (feedback|sub), data=subset(data,transfer=="Close"), family=binomial, control=glmercontrol)
Anova(accuracy.transfer)

## RT
rt.same <- lmer(rt ~ feedback*difflevel + (difflevel+feedback|sub), data=subset(data,transfer=="Direct"), control=control, REML = F)
anova(rt.same)
eta_squared(rt.same, alternative = "two")
rt.transfer <- lmer(rt ~ feedback*difflevel + (difflevel+feedback|sub), data=subset(data,transfer=="Close"), control=control, REML = F)
anova(rt.transfer)
eta_squared(rt.transfer, alternative = "two")

# Confidence
cj.same <- lmer(data = subset(data, transfer=="Direct"), cj ~ feedback*difflevel*cor + (cor+difflevel+feedback|sub), control=control, REML = F)
anova(cj.same)
eta_squared(cj.same, alternative = "two")
post_hoc <- emmeans(cj.same.interaction, pairwise ~ feedback | cor)
pairs(post_hoc)
cj.transfer <- lmer(data = subset(data, transfer=="Close"), cj ~ feedback*difflevel*cor + (cor+difflevel+feedback|sub), control=control, REML = F)
anova(cj.transfer)
eta_squared(cj.transfer, alternative = "two")
post_hoc <- emmeans(cj.transfer, pairwise ~ feedback | cor)
pairs(post_hoc)