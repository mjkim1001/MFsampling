set_x_env <- function(type='normal', x_mean=0, x_sd=1){
  if(type=='normal'){
    x_mean <<- x_mean
    x_sd <<- x_sd
    h<<- 0.5
    sample_x <<- function(n){
      x_vec <- rnorm(n,x_mean,x_sd)
      return(x_vec)
    }
    # Density functions
    fX <<- function(x){
      # Truncated normal distribution
      dnorm(x,x_mean,x_sd)
    }
    FX <<- function(x){
      pnorm(x,x_mean,x_sd)
    }
  }else if(type=="exp"){
    h<<- 1
    sample_x <<- function(n){
      x_vec <- rexp(n,1)
      return(x_vec)
    }
    # Density functions
    fX <<- function(x){
      # Truncated normal distribution
      dexp(x,1)
    }
    FX <<- function(x){
      pexp(x,1)
    }
  }
  else if(type=="weibull"){
    h<<- 4
    sample_x <<- function(n){
      x_vec <- rweibull(n,shape=1/2)
      return(x_vec)
    }
    # Density functions
    fX <<- function(x){
      # Truncated normal distribution
      dweibull(x,shape=1/2)
    }
    FX <<- function(x){
      pweibull(x,shape=1/2)
    }
  }
  
  fX <<- Vectorize(fX)
  FX <<- Vectorize(FX)
  
}
set_model<-function(type="p3",a=0,b=0){
  if(type == "p3"){
    m <<- function(x){(x+0.1*x^3)}
    mprime <<- function(x){(1+0.3*x^2)}
    minv <<- function(y){
      uniroot(function(x){m(x)-y}, c(-20,20))[['root']]
    }
    Nbreaks<<-5
  }else if(type=="hetero"){
    m <<- function(x){a*exp(x/b)}
    mprime <<- function(x){a/b*exp(x/b)}
    minv <<- function(y){ b*log(y/a) }
    Nbreaks<<-10
  }else if(type=="p5"){
    m <<- function(x){(x + 0.2*x^3 +0.02*x^5)}
    mprime <<- function(x){(1+0.6*x^2 +0.1*x^4)}
    minv <<- function(y){
      uniroot(function(x){m(x)-y}, c(-20,20))[['root']]
    }
    minv <<- Vectorize(minv)
    Nbreaks<<-5
  }else if(type=="non_monotone"){
    x_breaks=c(-10,-a,a,10)
    slopes=b*c(1.5,-1,1.5)
    m <<- function(x){
      ifelse(x < -a, slopes[1]*(x+a)+b*a, 
             ifelse(x < a, slopes[2]*x, 
                    slopes[3]*(x-a)-b*a))
    }
    y_breaks=m(x_breaks)
    mprime <<- function(x) {
      ifelse(x < -a, slopes[1],
             ifelse(x < a, slopes[2],
                    slopes[3]))
    }
    minv <<- function(y){
      value = c()
      for( idx in findintervals(y, y_breaks)){
        value = c(value, inverse_linear(y, slopes[idx], y_breaks[idx], x_breaks[idx]))
      }
      return(value)
    }
    Nbreaks<<-10
  }else if(type=="norm"){
    m <<- function(x){
      a*exp(-x^2/C^2)
    }
    mprime <<- function(x) {
      -2*a/C^2*x*exp(-x^2/C^2)
    }
    minv <<- function(y){
      if(y==a){return(0)}
      return(c(-1,1)*sqrt(-C^2*log(y/a)))
    }
    Nbreaks<<-10
  }
  m<<-Vectorize(m)
  mprime <<-Vectorize(mprime)
}
set_param <- function(noise_type, noise_sigma=1){
  if (noise_type == "noiseless"){
    sigma <- function(x){0}
    noise_sigma <<- 0
  } else if (noise_type == "homo"){
    sigma <- function(x){1}
    noise_sigma <<- noise_sigma
    h<<- h + noise_sigma/2
  }else if (noise_type == "hetero"){
    print("set_param should come after set_model")
    sigma <- function(x){abs(m(x))}
    h <<- 3*h
    noise_sigma <<- 1/6
  }
  return(sigma)}


set_y_env <- function(noise_type, is_monotone=TRUE, ytlim=NULL, a=0,b=0){
  sample_y <<- function(x){
    m(x) + sigma(x)*rnorm(1, mean=0, sd=noise_sigma)
  }
  sample_y <<- Vectorize(sample_y)
  
  if(!is_monotone){
    if(noise_type == "homo"){
      fYtilde <<- build_fYtilde_hat( mprime, minv, fX, ytlim=ytlim)
      fYtilde <<- Vectorize(fYtilde)
      fY <<- function(y){
        func <- function(z){fYtilde(y-z)*dnorm(z, mean=0, sd=noise_sigma)}
        return(integrate(func, lower=-5*noise_sigma, upper=5*noise_sigma)$value)
      }
      fY <<- Vectorize(fY)
    }else if(noise_type == "noiseless"){
      fY <<- build_fYtilde_hat( mprime, minv, fX)
      fY <<- Vectorize(fY)
    }
    return()
  }
  
  if (noise_type == "noiseless"){
    FY <<- function(y){
      return(FX(minv(y)))
    }
    fY <<- function(y){
      return(fX(minv(y))/mprime(minv(y)))
    }
    FY <<- Vectorize(FY)
    fY <<- Vectorize(fY)
    
  } else if (noise_type == "homo"){
    
    FYtilde <<- function(y){
      return(FX(minv(y)))
    }
    if(!is.null(ytlim)){
      fYtilde <<- function(y){
        if(length(ytlim)==1){
          if(y<=ytlim[1]){return(0)}
        }else{
          if(y>ytlim[2] | y<=ytlim[1]){return(0)}
        }
        return(fX(minv(y))/mprime(minv(y)))
      }
    }else{
      fYtilde <<- function(y){
        return(fX(minv(y))/mprime(minv(y)))
      }
    }
    fY <<- function(y){
      func <- function(z){fYtilde(y-z)*dnorm(z, mean=0, sd=noise_sigma)}
      return(integrate(func, lower=-5*noise_sigma, upper=5*noise_sigma)$value)
    }
    
    fY <<- Vectorize(fY)
    FYtilde <<- Vectorize(FYtilde)
    fYtilde <<- Vectorize(fYtilde)
    
  }else if(noise_type=="hetero"){
    fY <<- function(y){
      func <- function(z){1/(a*(1+z)) * dlnorm(y/(a*(1+z)),meanlog=x_mean/b, sdlog=1/b)*dnorm(z, mean=0, sd=noise_sigma)}
      return(integrate(func, lower=-5*noise_sigma, upper=5*noise_sigma)$value)
    }
    fY <<- Vectorize(fY)
    
    trans <<- function(y){ log(y) }
    sample_z <<- function(x){
      trans(m(x) + sigma(x)*rnorm(1, mean=0, sd=noise_sigma))
    }
    sample_z <<- Vectorize(sample_z) ## Now it is sample_z
    
    mtilde_prime <<- function(x){1/b}
    beta_prime <<- function(y){1/y}
  }
}


sample_tail_x <- function(N, r, xL, is_left=T){
  x_vec <- sample_x(N)
  # 원하는 조건을 만족하는 샘플을 필터링합니다.
  if(is_left) {
    x_vec <- x_vec[x_vec <= xL]
  } else {
    x_vec <- x_vec[x_vec >= xL]
  }
  # 필요한 만큼의 샘플을 얻을 때까지 반복합니다.
  while(length(x_vec) < r) {
    additional_samples <- sample_x(as.integer(N * (r-length(x_vec))/length(x_vec)))
    
    if(is_left) {
      additional_samples <- additional_samples[additional_samples <= xL]
    } else {
      additional_samples <- additional_samples[additional_samples >= xL]
    }
    
    x_vec <- c(x_vec, additional_samples)
  }
  
  return(x_vec[1:r])
}

