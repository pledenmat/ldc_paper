rm(list=ls())
library(optparse)
source("fastmerge.r")
source("build_hm.R")

args = commandArgs(trailingOnly=TRUE)
s = args[1]
# Load dataset --------------------------------------------------------

load("Data_hm.Rdata")
params <- read.csv("fitted_parameters_hm.csv")

subs <- unique(Data$sub)
cond_ordered <- c("minus","control","plus")
Ncond <- length(cond_ordered)
s <- subs[as.numeric(s)]

# Simulation from model ---------------------------------------------------
#Generate model simulations
nsim <- 1

# Simulation parameters
dt <- .001; ev_bound <- .5; ev_window <- dt; upperRT <- 5; sigma <- .1
timesteps <- upperRT/dt; ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)

par_pcor <- subset(param,sub==s)
par_pcor <- c(mean(par_pcor$bound),mean(par_pcor$ter),0,nsim,.1,.001,1,
              with(par_pcor,aggregate(drift,list(difflevel),mean))$x) #Drift rate

for (c in 1:Ncond) {
  par_pcor_cond <- subset(param,fixed_parameter=="pcor"&sub==s&condition==cond_ordered[c])
  par_pcor_cond <- c(mean(par_pcor_cond$bound),mean(par_pcor_cond$ter),0,nsim,.1,.001,1,
                     with(par_pcor_cond,aggregate(drift,list(difflevel),mean))$x) #Drift rate
  
  file_name <- paste0("heatmaps/ldc/hm_",s,"_",cond_ordered[c],".Rdata")
  if (file.exists(file_name)) {
    load(file_name)
  } else {
    hm_up <- build_hm(par_pcor_cond[(length(par_pcor_cond)-Ndiff+1):length(par_pcor_cond)],
                      sigma = sigma, ev_bound = ev_bound, dt = dt, 
                      ev_window = ev_window, upperRT = upperRT)
    save(hm_up,file=file_name)
  }
}

file_name <- paste0("heatmaps/ldc/hm_",s,".Rdata")
if (file.exists(file_name)) {
  load(file_name)
} else {
  hm_up <- build_hm(par_pcor[(length(par_pcor)-Ndiff+1):length(par_pcor)],
                    sigma = sigma, ev_bound = ev_bound, dt = dt, 
                    ev_window = ev_window, upperRT = upperRT)
  save(hm_up,file=file_name)
}