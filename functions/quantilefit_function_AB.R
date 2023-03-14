library(Rcpp)
library(myPackage)

model_conf_new_sqrt <- function(a,b,data){
  return(( (1+data$resp)/2)*(1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2 - b))))  + ((1-data$resp)/2)*(1 /(1 + exp((1/sqrt(data$rt2+.00001))* (a*data$evidence2 - b))) ))
}
#' chi_square_optim
#'
#' @param params Parameters of the model to be fitted + fixed parameters
#' @param observations Data to be fitted
#' @param returnFit if 0, then returns predicted data. Else, fit the parameters
#' @param binning True if confidence is binned (from 1 to 6)
#'
#' @return predicted data if returnFit == 0, else returns the prediction error
#'
#' @details This is the function used to simultaneously fit reaction time, accuracy and confidence judgments via quantile optimization.  
#'
#' @examples
#' 
chi_square_optim_AB <- function(params, observations, returnFit,
                                confRT_name = "RTconf", binning = T, 
                                condition = c("minus","control","plus")){
  #First, generate predictions:
  drift <- params[10:length(params)]
  params <- params[1:9]
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio','alpha','beta')
  
  trial = data.frame(matrix(NA,nrow=0,ncol=7));
  names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in drift) {
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=-d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- fastmerge(predictions,predictionsneg)
    predictions['drift'] <- d
    trial <- fastmerge(trial,predictions)
  }
  predictions <- trial
  predictions <- predictions[rep(seq_len(nrow(predictions)), 3), ]
  predictions$condition <- rep(condition,each = nrow(trial))
  
  predictions$cj_1 <- model_conf_new_sqrt(params["alpha"],params["beta"],predictions)

  
  
  if (binning) {
    predictions$cj <- as.numeric(cut(predictions$cj_1,
                                     breaks=seq(0,1,length.out = 7),
                                     include.lowest = TRUE))  
  }else{predictions$cj <- predictions$cj_1}
  
  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations 
  }else{ 
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond],]$cj
      e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond],]$cj  
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      if (binning) {
        c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        
        c_props_cj <- as.data.frame(table(c_observed[c_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        e_props_cj <- as.data.frame(table(e_observed[e_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
        e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
        obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
        
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj == 1),
          sum(c_predicted_cj == 2),
          sum(c_predicted_cj == 3),
          sum(c_predicted_cj == 4),
          sum(c_predicted_cj == 5),
          sum(c_predicted_cj == 6)
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj == 1),
          sum(e_predicted_cj == 2),
          sum(e_predicted_cj == 3),
          sum(e_predicted_cj == 4),
          sum(e_predicted_cj == 5),
          sum(e_predicted_cj == 6)
        ) / dim(predictions)[1]
        
      }else{
        c_quantiles_cj <- quantile(c_observed[c_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        e_quantiles_cj <- quantile(e_observed[e_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        if (any(is.na(e_quantiles_cj))) {
          e_quantiles_cj <- rep(0,5)
        }
        if (any(is.na(c_quantiles_cj))) {
          c_quantiles_cj <- rep(0,5)
        }  
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
          sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
          sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
          sum(c_predicted_cj > c_quantiles_cj[5])
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
          sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
          sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
          sum(e_predicted_cj > e_quantiles_cj[5])
        ) / dim(predictions)[1]
      }
      pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    #Return chiSquare
    return(chiSquare)
  }
}
chi_square_optim_AB_bfix <- function(params, observations, returnFit,
                                     confRT_name = "RTconf", binning = T, 
                                     condition = c("minus","control","plus")){
  #First, generate predictions:
  drift <- params[(length(params)-2):length(params)]
  alphas <- params[(length(params)-6):(length(params)-4)]
  params <- c(params[1:(length(params)-7)],params[length(params)-3])
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio','beta')
  
  trial = data.frame(matrix(NA,nrow=0,ncol=7));names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in drift) {
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=-d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- fastmerge(predictions,predictionsneg)
    predictions['drift'] <- d
    trial <- fastmerge(trial,predictions)
  }
  predictions <- trial
  predictions <- predictions[rep(seq_len(nrow(predictions)), 3), ]
  predictions$condition <- rep(condition,each = nrow(trial))
  predictions$cj_1 <- -99
  
  for (i in 1:length(alphas)) {
    predictions[predictions$condition == condition[i],"cj_1"] <- model_conf_new_sqrt(alphas[i],params["beta"],predictions[predictions$condition == condition[i],])  }
  
  
  if (binning) {predictions$cj <- as.numeric(cut(predictions$cj_1,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
  }else{predictions$cj <- predictions$cj_1}
  
  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations 
  }else{ 
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond],]$cj
      e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond],]$cj  
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      if (binning) {
        c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        
        c_props_cj <- as.data.frame(table(c_observed[c_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        e_props_cj <- as.data.frame(table(e_observed[e_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
        e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
        obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
        
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj == 1),
          sum(c_predicted_cj == 2),
          sum(c_predicted_cj == 3),
          sum(c_predicted_cj == 4),
          sum(c_predicted_cj == 5),
          sum(c_predicted_cj == 6)
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj == 1),
          sum(e_predicted_cj == 2),
          sum(e_predicted_cj == 3),
          sum(e_predicted_cj == 4),
          sum(e_predicted_cj == 5),
          sum(e_predicted_cj == 6)
        ) / dim(predictions)[1]
        
      }else{
        c_quantiles_cj <- quantile(c_observed[c_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        e_quantiles_cj <- quantile(e_observed[e_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        if (any(is.na(e_quantiles_cj))) {
          e_quantiles_cj <- rep(0,5)
        }
        if (any(is.na(c_quantiles_cj))) {
          c_quantiles_cj <- rep(0,5)
        }  
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
          sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
          sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
          sum(c_predicted_cj > c_quantiles_cj[5])
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
          sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
          sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
          sum(e_predicted_cj > e_quantiles_cj[5])
        ) / dim(predictions)[1]
      }
      pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    #Return chiSquare
    return(chiSquare)
  }
}
chi_square_optim_AB_afix <- function(params, observations, returnFit,
                                     confRT_name = "RTconf", binning = T, 
                                     condition = c("minus","control","plus")){
  #First, generate predictions:
  drift <- params[(length(params)-2):length(params)]
  betas <- params[(length(params)-5):(length(params)-3)]
  params <- params[1:(length(params)-6)]
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio','alpha')
  
  trial = data.frame(matrix(NA,nrow=0,ncol=7));
  names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in drift) {
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=-d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- fastmerge(predictions,predictionsneg)
    predictions['drift'] <- d
    trial <- fastmerge(trial,predictions)
  }
  predictions <- trial
  predictions <- predictions[rep(seq_len(nrow(predictions)), 3), ]
  predictions$condition <- rep(condition,each = nrow(trial))
  predictions$cj_1 <- -99
  
  for (i in 1:length(betas)) {
    predictions[predictions$condition == condition[i],"cj_1"] <- model_conf_new_sqrt(params["alpha"],betas[i],predictions[predictions$condition == condition[i],])  }
  
  
  if (binning) {predictions$cj <- as.numeric(cut(predictions$cj_1,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
  }else{predictions$cj <- predictions$cj_1}
  
  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations 
  }else{ 
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond],]$cj
      e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond],]$cj  
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      if (binning) {
        c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        
        c_props_cj <- as.data.frame(table(c_observed[c_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        e_props_cj <- as.data.frame(table(e_observed[e_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
        e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
        obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
        
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj == 1),
          sum(c_predicted_cj == 2),
          sum(c_predicted_cj == 3),
          sum(c_predicted_cj == 4),
          sum(c_predicted_cj == 5),
          sum(c_predicted_cj == 6)
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj == 1),
          sum(e_predicted_cj == 2),
          sum(e_predicted_cj == 3),
          sum(e_predicted_cj == 4),
          sum(e_predicted_cj == 5),
          sum(e_predicted_cj == 6)
        ) / dim(predictions)[1]
        
      }else{
        c_quantiles_cj <- quantile(c_observed[c_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        e_quantiles_cj <- quantile(e_observed[e_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        if (any(is.na(e_quantiles_cj))) {
          e_quantiles_cj <- rep(0,5)
        }
        if (any(is.na(c_quantiles_cj))) {
          c_quantiles_cj <- rep(0,5)
        }  
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
          sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
          sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
          sum(c_predicted_cj > c_quantiles_cj[5])
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
          sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
          sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
          sum(e_predicted_cj > e_quantiles_cj[5])
        ) / dim(predictions)[1]
      }
      pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
      # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
      pred_props_cj[pred_props_cj==0] <- .0000001
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    #Return chiSquare
    return(chiSquare)
  }
}
chi_square_optim_AB_nofix <- function(params, observations, returnFit,
                                      confRT_name = "RTconf", binning = T, 
                                      condition = c("minus","control","plus")){
  #First, generate predictions:
  drift <- params[(length(params)-2):length(params)]
  betas <- params[(length(params)-5):(length(params)-3)]
  alphas <- params[(length(params)-8):(length(params)-6)]
  params <- params[1:(length(params)-9)]
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio')
  
  trial = data.frame(matrix(NA,nrow=0,ncol=7));
  names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in drift) {
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(
      v=-d,a=params['a'],ter=params['ter'],z=params['z'],
      ntrials=params['ntrials']*dim(observations)[1]/2/length(drift)/length(condition),
      s=params['sigma'],dt=params['dt'],
      t2distribution=rep(observations[,confRT_name],times=params['ntrials']),
      postdriftmod=params['vratio']))
    
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- fastmerge(predictions,predictionsneg)
    predictions['drift'] <- d
    trial <- fastmerge(trial,predictions)
  }
  predictions <- trial
  predictions <- predictions[rep(seq_len(nrow(predictions)), 3), ]
  predictions$condition <- rep(condition,each = nrow(trial))
  predictions$cj_1 <- -99
  
  for (i in 1:length(betas)) {
    predictions[predictions$condition == condition[i],"cj_1"] <- model_conf_new_sqrt(alphas[i],betas[i],predictions[predictions$condition == condition[i],])  }
  
  
  if (binning) {predictions$cj <- as.numeric(cut(predictions$cj_1,breaks=seq(0,1,length.out = 7),include.lowest = TRUE))  
  }else{predictions$cj <- predictions$cj_1}
  
  
  
  
  #if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions)
    
    #If we are fitting the model, now compare these predictions to the observations 
  }else{ 
    # again, separate the predections according to the response
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (cond in 1:length(condition)) {
      c_predicted_cj <- c_predicted[c_predicted$condition == condition[cond],]$cj
      e_predicted_cj <- e_predicted[e_predicted$condition == condition[cond],]$cj  
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      if (binning) {
        c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        
        c_props_cj <- as.data.frame(table(c_observed[c_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        e_props_cj <- as.data.frame(table(e_observed[e_observed$condition == condition[cond],]$cj)/dim(observations)[1])
        c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] <- c_obs_proportion_cj[c_obs_proportion_cj$var1 %in% c_props_cj$Var1,"Freq"] + c_props_cj$Freq
        e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] <- e_obs_proportion_cj[e_obs_proportion_cj$var1 %in% e_props_cj$Var1,"Freq"] + e_props_cj$Freq
        obs_props_cj <- c(obs_props_cj,c_obs_proportion_cj$Freq,e_obs_proportion_cj$Freq)
        
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj == 1),
          sum(c_predicted_cj == 2),
          sum(c_predicted_cj == 3),
          sum(c_predicted_cj == 4),
          sum(c_predicted_cj == 5),
          sum(c_predicted_cj == 6)
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj == 1),
          sum(e_predicted_cj == 2),
          sum(e_predicted_cj == 3),
          sum(e_predicted_cj == 4),
          sum(e_predicted_cj == 5),
          sum(e_predicted_cj == 6)
        ) / dim(predictions)[1]
        
      }else{
        c_quantiles_cj <- quantile(c_observed[c_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        e_quantiles_cj <- quantile(e_observed[e_observed$condition == condition[cond],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        if (any(is.na(e_quantiles_cj))) {
          e_quantiles_cj <- rep(0,5)
        }
        if (any(is.na(c_quantiles_cj))) {
          c_quantiles_cj <- rep(0,5)
        }  
        c_pred_proportion_cj <- c(
          sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[2]) - sum(c_predicted_cj <= c_quantiles_cj[1]),
          sum(c_predicted_cj <= c_quantiles_cj[3]) - sum(c_predicted_cj <= c_quantiles_cj[2]),
          sum(c_predicted_cj <= c_quantiles_cj[4]) - sum(c_predicted_cj <= c_quantiles_cj[3]),
          sum(c_predicted_cj <= c_quantiles_cj[5]) - sum(c_predicted_cj <= c_quantiles_cj[4]),
          sum(c_predicted_cj > c_quantiles_cj[5])
        ) / dim(predictions)[1]
        
        e_pred_proportion_cj <- c(
          sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[2]) - sum(e_predicted_cj <= e_quantiles_cj[1]),
          sum(e_predicted_cj <= e_quantiles_cj[3]) - sum(e_predicted_cj <= e_quantiles_cj[2]),
          sum(e_predicted_cj <= e_quantiles_cj[4]) - sum(e_predicted_cj <= e_quantiles_cj[3]),
          sum(e_predicted_cj <= e_quantiles_cj[5]) - sum(e_predicted_cj <= e_quantiles_cj[4]),
          sum(e_predicted_cj > e_quantiles_cj[5])
        ) / dim(predictions)[1]
      }
      pred_props_cj <- c(pred_props_cj,c_pred_proportion_cj,e_pred_proportion_cj)
      # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
      pred_props_cj[pred_props_cj==0] <- .0000001
    }
    
    # Combine the quantiles for rts and cj
    if (binning) {obs_props <- c(obs_props,obs_props_cj)
    }else{obs_props <- c(obs_props,obs_props)}
    
    pred_props <- c(pred_props,pred_props_cj)
    # calculate chi square
    chiSquare = sum( (obs_props - pred_props) ^ 2)
    #Return chiSquare
    return(chiSquare)
  }
}
