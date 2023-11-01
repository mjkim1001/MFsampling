inverse_linear <- function(y_values, a, b,c) {
  return((y_values - b) / a + c)
}
fit_plr <- function(data, Nbreaks=5){
  o = lm(y~x, data)
  fit.AIC<- sapply(1:Nbreaks, function(x) tryCatch({
    extractAIC(segmented(o, npsi=x))
  }, error = function(e) {
    c(NA,NA)
  }))
  fit.AIC = cbind(extractAIC(o),fit.AIC)
  nbreak = which.min(fit.AIC[2,])-1
  
  if(nbreak){
    o =segmented(o, npsi = nbreak)
  } 
  return(o)
}

sample_plr <- function(o, n=1,xlim){
  
  if(is.null(o$psi)){ # If linear without breaks
    x_breaks = xlim
    slopes = o$coefficients[2]
  }else{
    slopes = slope(o)$x[,1]
    breaks = o$psi[,2]
    intercepts = c(o$coefficients[1], predict(o,data.frame(x=breaks)))
    x_breaks = c(xlim[1],breaks, xlim[2])
  }
  y_breaks <<- predict(o,data.frame(x=x_breaks))
  mhat = function(x){predict(o, data.frame(x=x))}
  if(sum(slopes>=0) %in% c(0,length(slopes))){
    # print("monotone")
    x_intervals = diff(x_breaks)
    normalized_slopes = slopes/sum(slopes * x_intervals)
    pX_hat <- function(x){
      return(normalized_slopes[findInterval(x, x_breaks)])
    }
    cum_prob = cumsum(normalized_slopes * x_intervals)
    count_FinvU = table(findInterval(runif(n), c(0,cum_prob)))
    new_samp=c()
    for(idx in 1:length(cum_prob)){
      if(!is.na(count_FinvU[as.character(idx)])) 
        new_samp= c(new_samp, runif(count_FinvU[as.character(idx)], min = x_breaks[idx], x_breaks[idx+1] ))
    }
  }else{
    mhat_inv <- function(y){
      value = c()
      for( idx in findintervals(y, y_breaks)){
        value = c(value, inverse_linear(y, slopes[idx], intercepts[idx], c(0,breaks)[idx]))
      }
      return(value)
    }
    mprime <- function(x){
      return(slopes[findInterval(x, x_breaks)])
    }
    fYtilde_hat <- build_fYtilde_hat( mprime, mhat_inv, fX)
    pX_hat <- build_pX_hat(mhat, fYtilde_hat, fX, xlim, normalize = F)
    sample_pX <- build_mh_sampler(pX_hat, xlim, n_burn=5e2)
    new_samp=sample_pX(n)
  }
  new_weights = pX_hat(new_samp)/normalize(pX_hat,xlim)
  return(list(new_samp, new_weights))
}
