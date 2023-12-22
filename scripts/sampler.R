plrSampler <- function(n0, n, xlim,Nbreaks=5){
  x_vec <- runif(n0, min=xlim[1], max=xlim[2])
  y_vec <- sample_y(x_vec)
  new_x <- new_y <- new_w <-c()
  for(sample in 1:n){
    data0 = data.frame(x=c(x_vec,new_x),y=c(y_vec,new_y))
    o = fit_plr(data0, Nbreaks)
    s = sample_plr(o, xlim=xlim)
    
    new_x <- c(new_x, s[[1]])
    new_y <- c(new_y, sample_y(s[[1]]))
    new_w <- c(new_w, fX(s[[1]])/s[[2]])
  }
  return(data.frame(x=new_x,w=new_w,y=new_y))
}
IS_sampler <-function(N,r,n, type="optimal", true_FX=TRUE, data=NULL,x_left=NULL,x_right=NULL){
  w = (n-2*r)/r 
  if(is.null(x_left)){
    if(!is.null(data)){
      data = data %>% arrange(x_vec)
      x_pop = data$x_vec
    }else{
      x_pop <- sort(sample_x(N))
    }
    xL <<- x_pop[r]; xR <<- x_pop[N-r+1]
    x_left = x_pop[1:r]
    x_right = x_pop[(N-r+1):N]
    if(!true_FX){ FXhat <<- function(x) mean(x_pop<=x) }
  }
  if(true_FX){
    PL = FX(xL);  PR = 1-FX(xR) 
  }else{
    PL = FXhat(xL);  PR = 1-FXhat(xR) 
    #If x_left is used, then FXhat should have been defined before.
  }
  
  if(grepl("plr", type)){
    sampled = plrSampler(n0=15, n= n-2*r, xlim= c(xL,xR))
    x_vec = c(x_left,sampled$x, x_right)
    weights = c(rep(w*PL,r),sampled$w,rep(w*PR,r))
    y_vec = c(sample_y(x_left), sampled$y, sample_y(x_right))
  }else{
    if(grepl("optimal", type)){
      norm_const <- normalize(pX_unnorm,c(xL,xR))
      pX_hat <- function(x){pX_unnorm(x)/norm_const}
      sample_pX <- build_mh_sampler(pX_hat, c(xL,xR), n_burn=1e3)
      x_vec = c(x_left,sample_pX(n-2*r), x_right)
      weights = c(rep(w*PL,r),fX(x_vec[(r+1):(n-r)])/pX_hat(x_vec[(r+1):(n-r)]),rep(w*PR,r))
    }else if(type=="uniform"){
      x_vec = c(x_left,runif(n-2*r, min=xL, max=xR), x_right)
      weights = c(rep(w*PL,r),fX(x_vec[(r+1):(n-r)])*(xR-xL),rep(w*PR,r))
      #x_vec <- runif(n, min=-5, max=5)
      #weights <- fX(x_vec)*10 # When true fX is known
    }else if(type=="random"){
      #x_vec <- sample(x_pop,n);
      x_vec <- sample_x(n)
      weights <- rep(1, n)
    }
    y_vec = sample_y(x_vec)
  }
  label= rep(c("L","0","R"),c(r,n-2*r,r))
  return(data.frame(x_vec,weights,y_vec,label))
}
threshold_IS <-function(data, itr=0, N0=N, rr=NULL){
  #When m is monotone
  #If rr is not given, return xL and xR based on the label.
  #If rr is given, then retrun yL and yR based on the (rr)th orderstatistic of Y.
  if(is.null(rr)){
    rL=which(data$x_vec == max((data%>%filter(label=="L"))$x_vec))
    rR=which(data$x_vec == min((data%>%filter(label=="R"))$x_vec))
    return(data.frame(xL=data$x_vec[rL], xR=data$x_vec[rR], yL=data$y_vec[rL], yR=data$y_vec[rR], ymin=min(data$y_vec), ymax=max(data$y_vec) , iteration=itr, N=N0))
  }else{
    data = data %>% arrange(y_vec)
    if(length(rr)>1){
      rL = rr[1];rR=nrow(data)- rr[2] +1
    }else{
      rL = rr
      rR = nrow(data)- rr +1
    }
    return(data.frame(xL=data$x_vec[rL], xR=data$x_vec[rR], yL=data$y_vec[rL], yR=data$y_vec[rR], ymin=min(data$y_vec), ymax=max(data$y_vec) , iteration=itr, N=N0, rr=rL))
  }
}

Fbar_generator<-function(data, top=T,d=1){
  # return `Fbar` function which returns empirical estimate of P(Y>u)
  y_vec = data$y_vec

  Fbar<-function(u){
    weights=data$weights*(nrow(data))/sum(data$weights)
    if(!top){ y_vec = -y_vec; u = -u }
    return( mean((y_vec>u) *(weights)^d) )
  }
  return(Fbar)
}


generate_grid <- function(data,h_support=NULL){
  if(is.null(h_support)){
    return(seq( min(data$y_vec), max(data$y_vec), length.out=500))
  }
  return(seq( min(data$y_vec)-h_support*3.75*noise_sigma, max(data$y_vec)+h_support*3.75*noise_sigma, length.out=500))
}

generate_dens<-function(data, proposal, grid=NULL, h0=h, itr=0, N0=N, rr=50, r=50,hetero=F){
  if(is.null(grid)){
    y_grid = generate_grid(data)}
  else{ y_grid = grid }
  
  if(proposal=="true"){ # When actual sampling density is used.
    return(data.frame("proposal"=proposal, "y"=y_grid, "est_dens"=fY(y_grid), "iteration"=-1, "N"=N0, "h"=NA, is_support=T,rr=rr))
  }
  model <- kde(x=data$y_vec, w=data$weights*(nrow(data))/sum(data$weights), h=h0)
  if(hetero){
    is_support = (trans(y_grid)>=min(data$y_vec) & trans(y_grid) <=max(data$y_vec) )
  }else{
    is_support = (y_grid>=min(data$y_vec) & y_grid <=max(data$y_vec) )
  }
  if(grepl("modified", proposal)){ # When GPD calibration is used
    #if(hetero){
    #  data$y = exp(data$y)
    #}
    
    y_sort = sort(data$y_vec)
    n = nrow(data)
    
    if(length(rr)>1){
      rl=rr[1];rr=rr[2]
      yL =y_sort[rl]; yR=y_sort[n-rr+1]
    }else{
      yL =y_sort[rr]; yR=y_sort[n-rr+1]
    }
    fit.gpd_R = fgpd.weight(data$y_vec, weight=data$weights, u = yR)
    fit.gpd_L = fgpd.weight(-data$y_vec, weight=data$weights, u = -yL)
    
    gpd_param = data.frame(yL = yL, yR=yR, xiL =fit.gpd_L$xi, xiR = fit.gpd_R$xi, sigmaL = fit.gpd_L$sigmau, sigmaR = fit.gpd_R$sigmau, rr = rr )
    
    FbarL = (Fbar_generator(data,top=F))(yL)
    FbarR = (Fbar_generator(data,top=T))(yR)
    xiL = fit.gpd_L$xi; xiR = fit.gpd_R$xi
    sigmaL = fit.gpd_L$sigmau; sigmaR = fit.gpd_R$sigmau
    
    cL=cR=1
    
    if(hetero){
      dens =beta_prime(y_grid) * sapply(trans(y_grid), function(x){
        if(x <= yL){ 
          return(evmix::dgpd(-x,u=-yL, xi = xiL,  sigmau = sigmaL , phiu=FbarL)*cL)
        }else if(x > yR){
          return(evmix::dgpd(x,u=yR, xi = xiR,  sigmau = sigmaR, phiu=FbarR)*cR)
        } else{
          return(predict(model, x=x))
        }
      })
    }else{
      dens =sapply((y_grid), function(x){
        if(x <= yL){ 
          return(evmix::dgpd(-x,u=-yL, xi = xiL,  sigmau = sigmaL , phiu=FbarL)*cL)
        }else if(x > yR){
          return(evmix::dgpd(x,u=yR, xi = xiR,  sigmau = sigmaR, phiu=FbarR)*cR)
        } else{
          return(predict(model, x=x))
        }
      })
    }
    
  }else{ # When importance sampler is used.
    if(hetero){
      dens = sapply(y_grid,function(x) beta_prime(x)*predict(model, x=trans(x)))
    }else{
      dens = sapply(y_grid,function(x) predict(model, x=x))
    }
  }
  
  return(data.frame("proposal"=proposal, "y"=y_grid, "est_dens"=dens, "iteration" = itr, "N"= N0 ,h=h0, is_support=is_support,rr=rr))
}

XY_sampler<-function(N){
  x_vec= sample_x(N)
  y_vec= sample_y(x_vec)
  return(data.frame(x_vec,y_vec, weights=1 ))
}

plot_sampler <-function(data,r=NULL){
  xL = data$x_vec[r]; xR=data$x_vec[nrow(data)-r+1]
  if(is.null(r)){
    return(ggplot(data) + geom_point(aes(x_vec,y_vec)) + ggtitle("Scatter plot of X and Y"))
  }else{
    return(ggplot(data) + geom_point(aes(x_vec,y_vec)) + ggtitle("Scatter plot of sampled X and Y") + 
             geom_vline(data= data.frame(xLR=c(xL,xR)),aes(xintercept = xLR, color="xL/xR")) +
             scale_colour_manual(
               values=c("red"),
               labels=c("xL/xR")) )
  }
}
