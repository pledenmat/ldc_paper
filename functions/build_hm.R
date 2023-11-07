library(fields)
library(mnormt)
library(distr)

#' Fast generation of p(cor) heatmap
#'
#' @param mu Drift rate used to generate the heatmap (only correctly works with 1 drift rate for now)
#' @param sigma 
#' @param dt 
#' @param ev_bound 
#' @param ev_window 
#' @param upperRT 
#'
#' @return
#' @export
#'
#' @examples
build_hm <- function(mu, sigma = .1, dt = .001, ev_bound = .5, ev_window = .01, upperRT = 5){
  ev_mapping <- seq(-ev_bound,ev_bound,by=ev_window)
  time_vec <- seq(0,upperRT,dt)
  if (length(mu)==1) {
    hm <- sapply(time_vec,function(time){
      if (time==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- dnorm(ev_mapping,mean=mu*time,sd=sigma*sqrt(time))
        dist_neg <- dnorm(ev_mapping,mean=-mu*time,sd=sigma*sqrt(time))
        res <- dist_pos / (dist_pos + dist_neg)
      }
      c(res)
    })
  }else{
    hm <- sapply(time_vec,function(time){
      if (time==0) {
        res <- rep(NA,length(ev_mapping))
        res[ceiling(length(res)/2)] <- .5
      }else{
        dist_pos <- lapply(mu, function(x){Norm(x*time,sigma*sqrt(time))})
        dist_neg <- lapply(-mu, function(x){Norm(x*time,sigma*sqrt(time))})
        dist_pos <- UnivarMixingDistribution(Dlist = dist_pos)
        dist_neg <- UnivarMixingDistribution(Dlist = dist_neg)
        dpos <- d(dist_pos)
        dneg <- d(dist_neg)
        res <- dpos(ev_mapping)/(dpos(ev_mapping)+dneg(ev_mapping))
        }
      c(res)
    })
  }
  return(t(hm))
}
