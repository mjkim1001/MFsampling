library(foreach)
library(doParallel)
library(beepr)
# This function is for generate_ldens for Figure 1
# We do not fix xL, xR here. effectively sample n points.
# Threshold returns r-th order statistic of sampled Y's
generate_ldens <- function(y_grid,m_idx=1,noise_type,r=50, n=150, N = 10^6, noise_sigma=1){
  ncores= 5
  registerDoParallel(cores=ncores)
  hetero=(noise_type=="hetero")
  combine_list <-function(x, ...){
    mapply(rbind, x, ..., SIMPLIFY = F)
  }
  
  # df contains estimated density value per each grid and method.
  # First, generate true fY values.
  
  df0 = generate_dens(data=NULL,"true", y_grid, N0=N, r=r)
  
  #For 100 iteration, save the following:
  # KDE using all N
  # IS using optimal
  # IS modified via GPD
  # KDE using random n
  # IS using predicted mean
  results <- foreach( i = 1:5, .combine = combine_list,.errorhandling = "remove") %dopar% {
    set.seed(i)
    trhold <- NULL
    df <- NULL
    dataUsed<-NULL
    # Obtain all N Y values.
    dataAll = XY_sampler(N)
    df = rbind(df,generate_dens(dataAll,"total", y_grid, h0=h/2, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(dataAll, N0=N,itr=i, rr= r) %>% mutate(proposal="total"))
    data = IS_sampler(N, r ,n=n,type="optimal", data=dataAll, true_FX = FALSE) 
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="optimal", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"optimal", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="optimal"))
    df = rbind(df,generate_dens(data,"modified", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="modified"))
    
    data = IS_sampler(N, r=r, n=n, type="uniform", data = dataAll, true_FX = FALSE)
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="uniform", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"uniform", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="uniform"))
    data = IS_sampler(N, r=r, type="random", data = dataAll, true_FX = FALSE)
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="random", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"random", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="random"))
    data = IS_sampler(N, r=r,n=n, type="plr", data = dataAll, true_FX = FALSE)
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="plr", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"plr", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero) )
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="plr"))
    df = rbind(df,generate_dens(data,"plr-modified", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="plr-modified"))
    if(hetero){
      trhold = trhold %>% mutate(yL = exp(yL), yR = exp(yR))
    }
    return(list(dataUsed=dataUsed, df=df, trhold=trhold))
  }
  beep()
  results$df = rbind(df0, results$df)
  saveRDS(results,  paste(homedir,sprintf("figures/m%sN%d_%s_noise%d.rds", m_idx,ifelse(log10(N) %% 1,N,log10(N)),noise_type,noise_sigma) ,sep="/") )
  save_png(results, y_grid, m_idx=m_idx, noise_type=noise_type, r=r,N=N, noise_sigma=noise_sigma)
  return(results)
}


save_png <-function(results, y_grid, plotType = c("random","uniform","optimal","modified", "plr","plr-modified","total"),m_idx, noise_type, N0=NULL,r=0,N=0,noise_sigma=1,figure=1){
  if(figure==1){
    p1 = results$df %>% mutate(type="m1", proposal=factor(proposal,levels=plotType)) %>% 
      filter(proposal %in% plotType) %>% 
      ggplot() + 
      geom_line(data= results$df %>% filter(proposal=="true" ) %>% mutate(type="m1")%>% 
                  dplyr::select(-c(proposal,N,rr, iteration)) ,aes(x=y, y=log(est_dens)), col='red') + 
      geom_line(aes(x=y, y=log(est_dens),group=iteration, color=factor(ifelse(is_support,"observed","extended"),levels=c("observed", "extended"))),alpha=0.1) + 
      #geom_point(data=subset(df,!is_support), aes(x=y, y=log(est_dens)),alpha=0.1,shape='.') + 
      facet_wrap(~proposal,ncol=3) +ylim(-30,0)+ xlim(min(y_grid),max(y_grid))+
      scale_color_manual(values=c("observed"="black", "extended"="blue", "yL/yR"='purple'), labels=c("observed", "extended",'yL/yR')) + 
      ylab("Log of Estimated Density") + labs(color="")+xlab("")+
      geom_vline(data= results$trhold %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
                 aes(xintercept = y, color="yL/yR"),linetype="dashed",color='purple') +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      theme(
        legend.position='top', 
        legend.direction="horizontal",
        text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.background = element_rect(fill="transparent"),
        strip.text.x = element_text(size=12)
      )
    
    
    fname = paste(homedir,sprintf("figures/p1_m%sN%d_%s_noise%d.png", m_idx,ifelse(log10(N) %% 1,N,log10(N)),noise_type,noise_sigma) ,sep="/")
    ggsave(filename = fname, plot = p1, width = 10, height = 8)
  }else{
    df_var1 <- subset(results$df,is_support) %>%
      group_by(proposal, y, N) %>% 
      summarise(
        dens_var=var(est_dens,na.rm=T)
      )
    df_var1$dens_var_scaled <- df_var1$dens_var/(fY(df_var1$y)^2)
    
    
    p2=df_var1 %>% 
      filter(proposal %in% c("random", "uniform", "optimal", "modified")) %>%
      mutate(proposal=factor(proposal, levels=c("random", "uniform", "optimal", "modified"))) %>%
      ggplot(aes(x=y, y=log(dens_var_scaled), group=proposal)) +
      geom_line(aes(linetype=proposal)) +
      ylab("Log of Scaled Variance") +
      scale_linetype_manual(
        values=c("dotted", "twodash","longdash", "solid"),
        labels=c("random", "uniform", "optimal", "modified"))+
      #coord_cartesian(ylim=c(-10, 10)) +
      geom_vline(data= results$trhold %>% filter(proposal %in% c( "modified")) %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
                 aes(xintercept = y, color="yL,yR"),linetype="dashed") +
      #facet_grid( ~noise) +
      scale_color_manual(values=c("green4"), labels=c("yL,yR")) + 
      theme(
        # legend.position=c(0.5, 0.9), 
        # legend.direction="horizontal",
        text = element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)
      ) +
      labs(linetype = "Proposal", color= "")+xlab("") + ylim(-6,6)+xlim(min(y_grid),max(y_grid))+
      theme(
        legend.position='top', 
        # legend.direction="horizontal",
        text = element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)
      )
    
    fname = paste(homedir,sprintf("figures/p2_m%sN%d_%s_noise%d.png", m_idx,ifelse(log10(N) %% 1,N,log10(N)),noise_type,noise_sigma) ,sep="/")
    ggsave(filename = fname, plot = p2, width = 10, height = 5)
  }
}

generate_vars <- function(y_grid,m_idx=1,noise_type,r=50, n=150, N = 10^6, noise_sigma=1){
  ncores= 4
  registerDoParallel(cores=ncores)
  hetero=(noise_type=="hetero")
  combine_list <-function(x, ...){
    mapply(rbind, x, ..., SIMPLIFY = F)
  }
  
  data = IS_sampler(N, r ,n=n,type="optimal", true_FX = FALSE) 
  results <- foreach( i = 1:10, .combine = combine_list,.errorhandling = "remove") %dopar% {
    set.seed(i)
    trhold <- NULL
    df <- NULL
    dataUsed<-NULL
    x_left = sample_tail_x(N, r, xL, is_left=TRUE)
    x_right = sample_tail_x(N, r, xR, is_left=FALSE)
    
    data = IS_sampler(N, r ,n=n,type="optimal", x_left=x_left,x_right=x_right, true_FX = FALSE) 
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="optimal", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"optimal", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="optimal"))
    df = rbind(df,generate_dens(data,"modified", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="modified"))
    
    data = IS_sampler(N, r=r, n=n, type="uniform", x_left=x_left,x_right=x_right, true_FX = FALSE)
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="uniform", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"uniform", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="uniform"))
    data = IS_sampler(N, r=r,n=n, type="random", true_FX = FALSE)
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="random", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"random", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="random"))
    
    return(list(dataUsed=dataUsed, df=df, trhold=trhold))
  }
  saveRDS(results,  paste(homedir,sprintf("figures/Varm%sN%d_%s_noise%d.rds", m_idx,ifelse(log10(N) %% 1,N,log10(N)),noise_type,noise_sigma) ,sep="/") )
  save_png(results, y_grid,  m_idx=m_idx, noise_type=noise_type, r=r,N=N, noise_sigma=1, figure=2)
  return(results)
}


generate_unif <- function(y_grid,m_idx=1,noise_type,r=50, n=150, N = 10^6, noise_sigma=1){
  ncores= 4
  registerDoParallel(cores=ncores)
  hetero=(noise_type=="hetero")
  combine_list <-function(x, ...){
    mapply(rbind, x, ..., SIMPLIFY = F)
  }

  
  df0 = generate_dens(data=NULL,"true", y_grid, N0=N, r=r)
  
  results <- foreach( i = 1:100, .combine = combine_list,.errorhandling = "remove") %dopar% {
    set.seed(i)
    trhold <- NULL
    df <- NULL
    dataUsed<-NULL
    # Obtain all N Y values.
    data = IS_sampler(N, r ,n=n,type="optimal", true_FX = FALSE) 
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="optimal", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"optimal", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="optimal"))
    df = rbind(df,generate_dens(data,"modified", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="modified"))
    
    data = IS_sampler(N, r=r, n=n, type="uniform", true_FX = FALSE)
    dataUsed = rbind(dataUsed, data %>% mutate(proposal="uniform", itr=i, N0=N, r=r))
    df = rbind(df,generate_dens(data,"uniform", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
    trhold = rbind(trhold, threshold_IS(data, N0=N,itr=i, rr= r) %>% mutate(proposal="uniform"))
    if(hetero){
      trhold = trhold %>% mutate(yL = exp(yL), yR = exp(yR))
    }
    return(list(dataUsed=dataUsed, df=df, trhold=trhold))
  }
  results$df = rbind(df0, results$df)
  saveRDS(results, sprintf("~/R/PCNN/PCNN/Unifm%sr%dN%d_%s_noise%d.rds", m_idx,r, ifelse(log10(N) %% 1,N,log10(N)), noise_type,noise_sigma) )
  save_png(results, y_grid, m_idx=m_idx, noise_type=noise_type, r=r,N=N, noise_sigma=noise_sigma)
  return(results)
}
