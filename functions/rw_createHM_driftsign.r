# 05/11/2020
#' RW_createHM_driftsign
#' @author Pierre Le Denmat
#' 
#' @description Generates "empirical" heatmaps (HM) by simulating both positive and negative drift rate trajectories. 
#'
#' @param mu : vector of drift rates
#' @param nsim : number of paths generated
#' @param dt : time resolution
#' @param intra_sv : random walk standard deviation
#' @param upperRT : maximum time
#' @param ev_bound : maximum/minimum evidence
#' @param ev_window : evidence resolution
#' @param HMtype : Type of heatmap returned, default all: ('sym' -> symmetrical HM, 'upper' -> upper bound hit heatmap, 'lower' -> lower bound hit)
#'
#' @return Matrix (or list of matrices) describing the Time x Evidence heatmap(s)
library(MALDIquant)
RW_createHM_driftsign <- function(mu, nsim=10000, dt=.0025, intra_sv=.1, upperRT=5, ev_bound=.5, ev_window=.01, HMtype=FALSE){
  
  starting_point <- 0
  ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
  
  timesteps <- upperRT/dt
  step_size <- sqrt(dt)*intra_sv
  time_index <- seq(0,timesteps-1) #To correctly index simulations into the heatmap
  
  passage_count <- rep.int(0,timesteps*length(ev_mapping))
  accuracy_count <- rep.int(0,timesteps*length(ev_mapping))
  passage_count_pos <- rep.int(0,timesteps*length(ev_mapping))
  
  # RW generation
  for (v in 1:length(mu)) {
    print(paste("Running drift rate level",v,"from",length(mu)))
    prob_up <- 0.5*(1+sqrt(dt)/intra_sv*mu[v])
    
    for(i in 1:nsim){
      position <- ((runif(timesteps) < prob_up)*2 - 1) * step_size
      position[1] <- position[1] + starting_point
      evidence <- cumsum(position)
      
      #Match to the heatmap
      mapping <- match.closest(evidence,ev_mapping) 
      value <- mapping + time_index*(length(ev_mapping)) #HM is aggregated per time row    
      
      if(mu[v] > 0){
        accuracy <- evidence > 0
        passage_count_pos[value] <- passage_count_pos[value] + 1 #All trajects where mu > 0
      }else{
        accuracy <- evidence <= 0
      }
      
      passage_count[value] <- passage_count[value] + 1 #All trajects
      accuracy_count[value[accuracy]] <- accuracy_count[value[accuracy]] + 1 #Correct traject
    }
  }
  # HM computation
  hmvector <- accuracy_count/passage_count
  hmvector_pos <- passage_count_pos/passage_count
  
  hm <- matrix(hmvector,nrow=timesteps,ncol=length(ev_mapping),byrow = TRUE)
  hm_pos <- matrix(hmvector_pos,nrow=timesteps,ncol=length(ev_mapping),byrow = TRUE)
  hm_neg <- 1-hm_pos
  
  if(tolower(HMtype)=='sym'){
    return(hm)
  }else if(tolower(HMtype)=='upper'){
    return(hm_pos)
  }else if(tolower(HMtype)=='lower'){
    return(hm_neg)
  }else{
    return(list('sym'= hm, 'upper'= hm_pos, 'lower'= hm_neg))
  }
  
}