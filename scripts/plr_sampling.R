# TITLE: noisy_sampling.R
# AUTHOR: Kevin O'Connor
# DATE: 12/9/21

library(ks)
library(rv)
library(pracma)

build_fY_hat <- function(y_vec, weights, h=NA){
  # Builds density estimate for Y from y_vec and weights using kde. 
  
  # Fit KDE
  n <- length(weights)
  weights <- weights*n/sum(weights)
  if (is.na(h)){
    kde_fitted <- kde(x=y_vec, w=weights)
  } else {
    kde_fitted <- kde(x=y_vec, w=weights, h=h)
  }  
  
  # Build function
  fY_hat <- function(y){predict(kde_fitted, x=y)}
  fY_hat <- Vectorize(fY_hat)
  return(fY_hat)
}

build_fYtilde_hat <- function(
    dot_mean,
    mhat_inv, 
    fX, ytlim=NULL
){
  if(!is.null(ytlim)){
    fYtilde_hat <- function(y){
      if(length(ytlim)==1){
        if(y<=ytlim[1]){return(0)}
      }else{
        if(y>ytlim[2] | y<=ytlim[1]){return(0)}
      }
      x_vec <- mhat_inv(y) # can return multiple values
      if(!length(x_vec)){
        print(sprintf("x null at y %f",y))
        return(0)
      }
      derivs <- pmax(1e-10,abs(sapply(x_vec, function(x){dot_mean(x)[1]})))
      fX_vec <- sapply(x_vec, fX)
      return(sum(fX_vec/derivs))
    }
  }else{
    fYtilde_hat <- function(y){
      x_vec <- mhat_inv(y) # can return multiple values
      # print(paste(sprintf("fYtilde called at y %f, x",y),x_vec))
      if(!length(x_vec)){
        print(sprintf("x null at y %f",y))
        return(0)
      }
      derivs <- pmax(1e-10,abs(sapply(x_vec, function(x){dot_mean(x)[1]})))
      fX_vec <- sapply(x_vec, fX)
      return(sum(fX_vec/derivs))
    }
  }
  
  fYtilde_hat = Vectorize(fYtilde_hat)
  return(fYtilde_hat)
}

normalize <- function(gX, xlim) {
  norm_const <- tryCatch({
    integrate(
      gX,
      lower = xlim[1],
      upper = xlim[2],
      subdivisions = 1e3,
      stop.on.error = FALSE
    )$value
  },
  error = function(e) {
    try(norm_const <- area(gX, xlim[1], xlim[2], eps = 1e-3))
    if(is.null(norm_const)){
      try(norm_const = area(gX, xlim[1], xlim[2], eps = 1e-2))
    }
    if(is.null(norm_const)){
      norm_const = integral(gX, xlim[1], xlim[2])
    }
    
    return(norm_const)
  })
  
  return(norm_const)
}

build_pX_hat <- function(
    mhat, 
    fYtilde_hat, 
    fX, 
    xlim=c(0,1),  # Range of X values
    normalize=T
){
  # Builds estimate of optimal proposal density pX from estimate of m and the 
  # estimated density of Ytilde.
  pX_hat_unnorm <- function(x){
    if (x < xlim[1] | x > xlim[2]){
      return(0)
    }
    #  print(sprintf("call fYtilde at x %f, y %f",x,mhat(x)))
    denom <- fYtilde_hat(mhat(x))
    if(is.na(denom)){
      print(sprintf("NA at x %f y %f",x, mhat(x)))
      return(0)}
    if (denom > 0){
      return(fX(x)/denom)
    } else {
      sprintf("denom:%f",denom)
      return(NA)
    }
  }
  pX_hat_unnorm <- Vectorize(pX_hat_unnorm)
  if(normalize){
    pX_hat = function(x){pX_hat_unnorm(x)/normalize(pX_hat_unnorm,xlim)}
    return(Vectorize(pX_hat))
  }else{
    return(pX_hat_unnorm)
  }
}
