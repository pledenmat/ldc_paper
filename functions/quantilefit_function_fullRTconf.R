library(Rcpp)
library(myPackage)

model_conf_new_sqrt <- function(a,b,data){
  return(( (1+data$resp)/2)*(1 /(1 + exp( (1/sqrt(data$rt2+.00001))* (-a*data$evidence2 - b))))  + ((1-data$resp)/2)*(1 /(1 + exp((1/sqrt(data$rt2+.00001))* (a*data$evidence2 - b))) ))
}
model_conf_notime <- function(a,b,data){
  return(( (1+data$resp)/2)*(1 /(1 + exp( (-a*data$evidence2 - b))))  + ((1-data$resp)/2)*(1 /(1 + exp( (a*data$evidence2 - b))) ))
}
#' chi_square_optim
#'
#' @param params Parameters of the model to be fitted + fixed parameters
#' @param observations Data to be fitted
#' @param returnFit if 0, then returns predicted data. Else, fit the parameters
#' @param binning True if confidence is binned (from 1 to 6)
#' @param denominator Divide by predicted values in the goal function if True
#'
#' @return predicted data if returnFit == 0, else returns the prediction error
#'
#' @details This is the function used to simultaneously fit reaction time, accuracy and confidence judgments via quantile optimization.  
#'
#' @examples
#' 
chi_square_optim <- function(params, observations, returnFit,confRT_name = "RTconf", binning = T, denominator = F, model = 1, hm_up = NULL, hm_low = NULL ){
  #First, generate predictions:
  drift <- params[10:length(params)]
  params <- params[1:9]
  names(params) <- c('a','ter','z','ntrials','sigma','dt','vratio','alpha','beta')
  
  trial = data.frame(matrix(NA,nrow=0,ncol=7));names(trial) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  for (d in drift) {
    predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(v=d,a=params['a'],ter=params['ter'],z=params['z'],ntrials=params['ntrials']*dim(observations)[1]/2,s=params['sigma'],dt=params['dt'],t2distribution=rep(observations[,confRT_name],times=params['ntrials']),postdriftmod=params['vratio']))
    predictionsneg <- data.frame(DDM_with_confidence_slow_fullconfRT(v=-d,a=params['a'],ter=params['ter'],z=params['z'],ntrials=params['ntrials']*dim(observations)[1]/2,s=params['sigma'],dt=params['dt'],t2distribution=rep(observations[,confRT_name],times=params['ntrials']),postdriftmod=params['vratio']))
    names(predictions) <- c('rt','resp','cor','evidence2','rt2','cj')
    names(predictionsneg) <- c('rt','resp','cor','evidence2','rt2','cj')
    predictions <- fastmerge(predictions,predictionsneg)
    predictions['drift'] <- d
    trial <- fastmerge(trial,predictions)
  }
  predictions <- trial
  
  if (model==1) {
    predictions$cj_1 <- model_conf_new_sqrt(params["alpha"],params["beta"],predictions)
  }else if (model==2) {
    predictions$cj_1 <- model_conf_notime(params["alpha"],params["beta"],predictions)
  }else if (model==3) {
    predictions$cj_1 <- model_conf_noalpha(params["alpha"],params["beta"],predictions)
  }else if (model==4) {
    predictions$cj_1 <- model_conf_nobeta(params["alpha"],params["beta"],predictions)
  }else if (model==5){
    dt <- params["dt"];ev_bound <- .5; ev_window <- dt*5; upperRT <- 5
    ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window);timesteps <- upperRT/dt
    #match to the heatmap
    predictions$closest_evdnc2 <- match.closest(predictions$evidence2,ev_mapping)
    predictions$closest_cj <- match.closest(predictions$cj,ev_mapping)
    predictions$temprt2 <- predictions$rt2;
    predictions$temprt2[predictions$temprt2>5] <- 5 #heatmap doesn't go higher
    predictions$temprt2 <- predictions$temprt2*timesteps/5 #scale with the heatmap, between 0 and 2000
    
    predictions$cj_1 <- NA
    predictions[predictions$resp==1,]$cj_1 <- hm_up[(predictions[predictions$resp==1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==1,]$temprt2)]
    predictions[predictions$resp==-1,]$cj_1 <- hm_low[(predictions[predictions$resp==-1,]$closest_evdnc2-1)*timesteps+round(predictions[predictions$resp==-1,]$temprt2)]
    predictions <- predictions[complete.cases(predictions$cj_1),]
  }  
  
  

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
    
    coherences <- sort(unique(observations$coh)) # CHANGE ACCORDING TO DATASET
    
    obs_props <- NULL; pred_props <- NULL;obs_props_cj <- NULL; pred_props_cj <- NULL
    for (d in 1:length(drift)) {
      # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
      c_quantiles <- quantile(c_observed[c_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      e_quantiles <- quantile(e_observed[e_observed$coh == coherences[d],]$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
      if (any(is.na(e_quantiles))) {
        e_quantiles <- rep(0,5)
      }
      if (any(is.na(c_quantiles))) {
        c_quantiles <- rep(0,5)
      }      
      # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
      prop_obs_c <- dim(c_observed[c_observed$coh == coherences[d],])[1] / dim(observations)[1]
      prop_obs_e <- dim(e_observed[e_observed$coh == coherences[d],])[1] / dim(observations)[1]
      
      c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
      e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
      obs_props <- c(obs_props,c_obs_proportion,e_obs_proportion)
      
      c_predicted_rt <- sort(c_predicted[c_predicted$drift == drift[d],]$rt)
      e_predicted_rt <- sort(e_predicted[e_predicted$drift == drift[d],]$rt)
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data (scale by N?)
      c_pred_proportion <- c(
        sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[2]) - sum(c_predicted_rt <= c_quantiles[1]),
        sum(c_predicted_rt <= c_quantiles[3]) - sum(c_predicted_rt <= c_quantiles[2]),
        sum(c_predicted_rt <= c_quantiles[4]) - sum(c_predicted_rt <= c_quantiles[3]),
        sum(c_predicted_rt <= c_quantiles[5]) - sum(c_predicted_rt <= c_quantiles[4]),
        sum(c_predicted_rt > c_quantiles[5])
      ) / dim(predictions)[1]
      
      e_pred_proportion <- c(
        sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[2]) - sum(e_predicted_rt <= e_quantiles[1]),
        sum(e_predicted_rt <= e_quantiles[3]) - sum(e_predicted_rt <= e_quantiles[2]),
        sum(e_predicted_rt <= e_quantiles[4]) - sum(e_predicted_rt <= e_quantiles[3]),
        sum(e_predicted_rt <= e_quantiles[5]) - sum(e_predicted_rt <= e_quantiles[4]),
        sum(e_predicted_rt > e_quantiles[5])
      ) / dim(predictions)[1]
      pred_props <- c(pred_props,c_pred_proportion,e_pred_proportion)
      
      # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
      pred_props[pred_props==0] <- .0000001
      
      # Now, do the same for confidence

      
      c_predicted_cj <- c_predicted[c_predicted$drift == drift[d],]$cj
      e_predicted_cj <- e_predicted[e_predicted$drift == drift[d],]$cj  
      # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data
      if (binning) {
        c_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        e_obs_proportion_cj <- data.frame(var1=1:6,Freq=0)
        
        c_props_cj <- as.data.frame(table(c_observed[c_observed$coh == coherences[d],]$cj)/dim(observations)[1])
        e_props_cj <- as.data.frame(table(e_observed[e_observed$coh == coherences[d],]$cj)/dim(observations)[1])
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
        c_quantiles_cj <- quantile(c_observed[c_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
        e_quantiles_cj <- quantile(e_observed[e_observed$coh == coherences[d],]$cj, probs = c(.1,.3,.5,.7,.9), names = FALSE)
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
    if (denominator) {
      chiSquare = sum( ( (obs_props - pred_props) ^ 2) / pred_props )
    }else{chiSquare = sum( (obs_props - pred_props) ^ 2)}
    #Return chiSquare
    return(chiSquare)
  }
}
