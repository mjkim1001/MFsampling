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
IS_sampler <-function(N,r, type="optimal", n=0, true_FX=TRUE, data=NULL,x_left=NULL,x_right=NULL){
  if(n==0){n=3*r; w=1}else{ 
    w = (n-2*r)/r }
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
  }
  
  if(grepl("plr", type)){
    sampled = plrSampler(n0=10, n= n-2*r, xlim= c(xL,xR))
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
  #If rr is given, then retrun yL and yR based on the rrth orderstatistic of Y.
  if(is.null(rr)){
    rL=which(data$x_vec == max((data%>%filter(label=="L"))$x_vec))
    rR=which(data$x_vec == min((data%>%filter(label=="R"))$x_vec))
    return(data.frame(xL=data$x_vec[rL], xR=data$x_vec[rR], yL=data$y_vec[rL], yR=data$y_vec[rR], ymin=min(data$y_vec), ymax=max(data$y_vec) , iteration=itr, N=N0))
  }else{
    data = data %>% arrange(y_vec)
    rL = rr
    rR = nrow(data)- rr +1
    return(data.frame(xL=data$x_vec[rL], xR=data$x_vec[rR], yL=data$y_vec[rL], yR=data$y_vec[rR], ymin=min(data$y_vec), ymax=max(data$y_vec) , iteration=itr, N=N0, rr=rr))
  }
}

Fbar_generator<-function(data, r, top=T){
  x_vec = data$x; y_vec = data$y
  n = nrow(data)
  xL = x_vec[r]; xR = x_vec[n-r+1]
  PL = FX(xL);PR = 1-FX(xR)
  Fbar<-function(u){
    if(top){
      range = ((n-r+1):n)
    }else{
      range = (1:r)
    }
    PP <-ifelse(top,PR,PL)
    if(!top){ y_vec = -y_vec; u = -u }
    return((xR-xL) * mean(((y_vec[(r+1):(n-r)])>u) * fX((x_vec[(r+1):(n-r)]))) * (1-PR-PL) +
             PP * mean((y_vec[range])>u))
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
    is_support = (b(y_grid)>=min(data$y_vec) & b(y_grid) <=max(data$y_vec) )
  }else{
    is_support = (y_grid>=min(data$y_vec) & y_grid <=max(data$y_vec) )
  }
  if(grepl("modified", proposal)){ # When GPD calibration is used
    #if(hetero){
    #  data$y = exp(data$y)
    #}
    
    gpd_param <<- NULL
    y_sort = sort(data$y)
    for( rrr in (rr-5):(rr+5)){
      n=nrow(data)
      yL =y_sort[rrr]; yR=y_sort[n-rrr+1]
      fit.gpd_R = fgpd.weight(data$y, weight=data$w, u = yR)
      fit.gpd_L = fgpd.weight(-data$y, weight=data$w, u = -yL)
      
      gpd_param = rbind(gpd_param, data.frame(yL = yL, yR=yR, xiL =fit.gpd_L$xi, xiR = fit.gpd_R$xi, sigmaL = fit.gpd_L$sigmau, sigmaR = fit.gpd_R$sigmau, rr = rrr ))
    }
    
    yL =y_sort[rr]; yR=y_sort[n-rr+1]
    
    FbarL = (Fbar_generator(data,r,top=F))(yL)
    FbarR = (Fbar_generator(data,r,top=T))(yR)
    sigmaL = median(gpd_param$sigmaL)
    sigmaR = median(gpd_param$sigmaR)
    xiL = median(gpd_param$xiL);xiR = median(gpd_param$xiR)
    cL = predict(model, x=yL-0.01)/evmix::dgpd(-yL+0.01,u=-yL, xi = xiL, sigmau = sigmaL, phiu=FbarL)
    cR = predict(model, x=yR+0.01)/evmix::dgpd(yR+0.01,u=yR, xi = xiR, sigmau = sigmaR, phiu=FbarR)
    cL=cR=1
    # dens = sapply(y_grid, function(x){
    #   if(x <= yL){ 
    #     return(evmix::dgpd(-x,u=-yL, xi = xiL,  sigmau = sigmaL , phiu=FbarL)*cL)
    #   }else if(x > yR){
    #     return(evmix::dgpd(x,u=yR, xi = xiR,  sigmau = sigmaR, phiu=FbarR)*cR)
    #   } else{
    #     if(hetero){
    #       return(beta_prime(x)*predict(model, x=b(x)))
    #     }else{
    #       return(predict(model, x=x))
    #     }
    #   }
    # })
    if(hetero){
      dens =beta_prime(y_grid) * sapply(b(y_grid), function(x){
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
      dens = sapply(y_grid,function(x) beta_prime(x)*predict(model, x=b(x)))
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

plot_data<-function(N, grid, type="optimal",r=50, rr0=40){
  #temporary plotting for myself
  if(type=="all"){
    data = XY_sampler(N)
  }else{
    data= IS_sampler(N,r=50, true_FX = F)
  }
  trhold <- do.call(rbind, lapply(c(10, 20, 40, 50), function(rr){threshold_IS(data, N0=N, rr = rr)}))
  ylim_plot <- summary(data$y_vec)[c(1,6)]
  
  if(type=="optimal"){
    df = generate_dens(data,"optimal", grid, h0=h, N0=N, r=r)
  }else if(type=="modified"){
    df = generate_dens(data,"modified", grid, h0=h, N0=N, rr=rr0, r=r)
  }else if(type=="all"){
    df = generate_dens(data,"random", grid, h0=h, N0=N, r=r)
  }
  
  plot(df$y[df$y>= ylim_plot[1] & df$y <= ylim_plot[2]], log(df$est_dens[df$y>= ylim_plot[1] & df$y <= ylim_plot[2]]),pch='.',ylab="log dens", col=rgb(0,0,0,0.7), ylim = c(-20,0),xlim=c(-25,25), xlab="y")
  points(df$y[df$y < ylim_plot[1] | df$y > ylim_plot[2]], log(df$est_dens[df$y < ylim_plot[1] | df$y > ylim_plot[2]]),pch='.',ylab="log dens", col=rgb(0,1,0,0.7), ylim = c(-20,0),xlim=c(-25,25))
  curve( log(fY(x)), ylim[1],ylim[2], add=T, col="blue")
  
  for(i in 1:9) {
    if(type=="all"){
      data = XY_sampler(N)
    }else{
      data= IS_sampler(N,r=50)
    }
    trhold <- rbind(trhold,do.call(rbind, lapply(c(10, 20, 40, 50), function(rr){threshold_IS(data, N0=N, rr = rr)})))
    ylim_plot <- summary(data$y_vec)[c(1,6)]
    
    if(type=="optimal"){
      df = generate_dens(data,"optimal", grid, h0=h, N0=N, r=r)
    }else if(type=="modified"){
      df = generate_dens(data,"modified", grid, h0=h, N0=N, rr= rr0, r=r)
    }else if(type=="all"){
      df = generate_dens(data,"random", grid, h0=h/2, N0=N, r=r)
    }
    points(df$y[df$y>= ylim_plot[1] & df$y <= ylim_plot[2]], log(df$est_dens[df$y>= ylim_plot[1] & df$y <= ylim_plot[2]]),pch='.',ylab="log dens", col=rgb(0,0,0,0.7), ylim = c(-25,0),xlim=c(-30,30), xlab="y")
    points(df$y[df$y < ylim_plot[1] | df$y > ylim_plot[2]], log(df$est_dens[df$y < ylim_plot[1] | df$y > ylim_plot[2]]),pch='.',ylab="log dens", col=rgb(0,1,0,0.7), ylim = c(-25,0),xlim=c(-30,30))
  }
  
  trhold= as.data.frame(trhold)
  if(type=="modified"){
    r_range<- c(rr0)
  }else{
    r_range<- c(50,20,10)
  }
  # Add vertical lines at certain percentiles of data$y_vec
  for(rrr in r_range){
    for( yy in (trhold %>% filter(rr==rrr))[c('yL','yR')] ){
      abline(v=yy , col=c("red","green4","purple")[which(rrr==r_range)], lty=2)
    }
    text(x=trhold %>% filter(rr==rrr) %>% summarize(yL=min(yL)-2,yR=max(yR)+2), y = -2, srt=90, sprintf("%dth",rrr),col=c("red","green4","purple")[which(rrr==r_range)])
  }
  
}
