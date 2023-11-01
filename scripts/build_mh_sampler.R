
build_mh_sampler <- function(dens, lim, n_burn=100){
  sampler <- function(n){
    propose <- function(x){
      runif(1, min=lim[1], max=lim[2])
    }
    propose_dens <- function(x, y){
      dunif((lim[1]+lim[2])/2, min=lim[1], max=lim[2])
    }
    burn_in <- function(x0){
      x <- x0
      for (i in 1:n_burn){
        x_tmp <- propose(x)
        rat <- dens(x_tmp)*propose_dens(x_tmp, x)/(dens(x)*propose_dens(x, x_tmp))
        a <- min(c(1, rat))
        if(is.na(a)){next}
        if (a >= runif(1)){
          x <- x_tmp
        }
      }
      return(x)
    }
    
    # Sample burn-in
    x0 <- propose(0)
    x <- burn_in(x0)
    
    # Draw "independent" samples
    x_vec <- c()
    for (i in 1:n){
      x_vec <- c(x_vec, burn_in(x))
    }
    return(x_vec)
  }
}
