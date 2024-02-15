
h=0.25
df = NULL
for(i in 1:100){
  data = as_tibble(results$dataUsed) %>% filter(proposal=="optimal" & itr==i) %>% dplyr::select(x_vec,weights,y_vec)
  df= rbind(df,generate_dens(data,"optimal", y_grid, h0=h, N0=10^6, itr=i, rr= r, r=r, hetero=F))
}
results$df = rbind(results$df %>% filter(proposal != "optimal"), df)

df = rbind(df,generate_dens(data,"modified", y_grid, h0=h, N0=N, itr=i, rr= r, r=r, hetero=hetero))
table((as_tibble(results$dataUsed) %>% filter(proposal=="optimal"))$N0)
i=1
(as_tibble(results$dataUsed))%>% filter(itr==i & proposal=="uniform")

save_figure <-function(results, plotType = c("random","uniform","optimal","modified", "plr","plr-modified","total"),m_idx, noise_type, N0=NULL,r=0,noise_sigma=0){
  results$trhold = results$trhold %>% 
    mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "random - N0", "modified - N", "PLR - N", "uniform - N","PLR modified - N"))) %>%
    filter( proposal!="uniform - N")
  df0 = results$df %>% filter(proposal=="true") 
  results$df = results$df %>% mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "random - N0", "modified - N", "PLR - N","uniform - N","PLR modified - N","true"))) %>% filter( proposal!="uniform - N")
  
  p1= ggplot() + 
    geom_line(data= df0%>% filter( proposal!="uniform - N") %>%
                dplyr::select(-c(proposal,N,rr, iteration)) ,aes(x=y, y=log10(est_dens)), col='red') + 
    #geom_line(aes(x=y, y=log10(est_dens),group=iteration, color=factor(ifelse(is_support,"observed","extended"),levels=c("observed", "extended"))),alpha=0.1) + 
    geom_point(data=subset(results$df,is_support & proposal!="true" & proposal!="uniform - N"), aes(x=y, y=log10(est_dens), color = "observed"),alpha=0.1,shape='.') +
    geom_point(data=subset(results$df,!is_support& proposal!="true" & proposal!="uniform N "), aes(x=y, y=log10(est_dens), color = "extended"),alpha=0.1,shape='.') +
    facet_wrap(~proposal, ncol=3) +ylim(-10,0)+# xlim(-25,25)+
    ylab("Log10 of Estimated Density") + labs(color="")+xlab("Y")+
    geom_vline(data= results$trhold %>% filter( proposal!="uniform - N") %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
               aes(xintercept = y, color="25th Y"),linetype="dashed",alpha=0.1) +
    scale_color_manual(
      values=c("observed"="black", "extended"="blue", "25th Y"='green4'), 
      labels=c("observed", "extended",'25th Y'),
      breaks=c("observed", "extended", "25th Y")
    )+
    guides(color = guide_legend(override.aes = list(alpha = 1,
                                                    linetype = c("blank", "blank", "dashed"), 
                                                    shape = c(19, 19, NA) 
                                                    ))) +
    theme(
      legend.position=ifelse(m_idx==1,'top','None'), 
      legend.direction="horizontal",
      text=element_text(size=12),
      legend.title=element_text(size=12),
      legend.text=element_text(size=12),
      legend.background = element_rect(fill="transparent"),
      strip.text.x = element_text(size=12)
    )
  
  fname = sprintf("./MFsampling/final2/f1m%dnoise%d.png", m_idx,noise_sigma) 
  ggsave(filename = fname, plot = p1, width = 10, height = ifelse(m_idx==1,5,4.5))
}

results <- readRDS(paste(homedir,"MFsampling/final2/mp1N6000000_homo_noise6.rds",sep = '/'))
save_figure(results, m_idx=1, noise_sigma=6)

results2 <- readRDS(paste(homedir,"MFsampling/final2/mPMN6_homo_noise6.rds",sep = '/'))
save_figure(results2, m_idx=2, noise_sigma=6)


## Figure1 desity for m3
results <- readRDS(paste(homedir,"MFsampling/final2/mhetero-homoN6_homo_noise1.rds",sep = '/'))
results2 <- readRDS(paste(homedir,"MFsampling/final2/mheteroN6_hetero_noise1.rds",sep = '/'))

#hetero
results$df = results$df %>% filter(proposal %in% c("modified", "optimal", "random", "true")) %>% mutate(type="m3")
results2$df = results2$df %>% filter(proposal %in% c("modified", "optimal", "random", "true")) %>% mutate(type="m3 - transformed")

results$trhold = rbind(results$trhold %>% 
                         mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "modified - N"))) %>% filter( proposal!="uniform - N") %>% mutate(type="m3") ,
                       results2$trhold %>% 
                         mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "modified - N"))) %>% filter( proposal!="uniform - N") %>% mutate(type="m3 - transformed"))

df0 = rbind(results$df %>% filter(proposal=="true")%>% mutate(type="m3"), results$df %>% filter(proposal=="true")%>% mutate(type="m3 - transformed"))

results$df = rbind(results$df %>% mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "modified - N","true"))) %>% 
                     filter( proposal!="uniform - N"), results2$df %>% mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "modified - N","true"))) %>% filter( proposal!="uniform - N"))

p= ggplot() + 
  geom_line(data= df0%>% filter( proposal!="uniform - N") %>%
              dplyr::select(-c(proposal,N,rr, iteration)) ,aes(x=y, y=log10(est_dens)), col='red') + 
  #geom_line(aes(x=y, y=log10(est_dens),group=iteration, color=factor(ifelse(is_support,"observed","extended"),levels=c("observed", "extended"))),alpha=0.1) + 
  geom_point(data=subset(results$df,is_support & proposal!="true" & proposal!="uniform - N"), aes(x=y, y=log10(est_dens), color = "observed"),alpha=0.1,shape='.') +
  geom_point(data=subset(results$df,!is_support& proposal!="true" & proposal!="uniform N "), aes(x=y, y=log10(est_dens), color = "extended"),alpha=0.1,shape='.') +
  facet_grid(type~proposal) +ylim(-10,0)+# xlim(-25,25)+
  ylab("Log10 of Estimated Density") + labs(color="")+xlab("Y")+
  geom_vline(data= results$trhold %>% filter( proposal!="uniform - N") %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
             aes(xintercept = y, color="25th Y"),linetype="dashed",alpha=0.1) +
  scale_color_manual(
    values=c("observed"="black", "extended"="blue", "25th Y"='green4'), 
    labels=c("observed", "extended",'25th Y'),
    breaks=c("observed", "extended", "25th Y")
  )+
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  linetype = c("blank", "blank", "dashed"), 
                                                  shape = c(19, 19, NA) ))) +
  theme(
    legend.position='None', 
    legend.direction="horizontal",
    text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.background = element_rect(fill="transparent"),
    strip.text.x = element_text(size=12)
  )



ggsave(filename="./MFsampling/final2/f1m3.png", plot=p, width=10, height = 4.5)





# 6.2.1
df<-NULL
trhold<-NULL
results <- readRDS(paste(homedir,"MFsampling/final2/mp1N6000000_homo_noise6.rds",sep = '/'))
df0 = results$df %>% filter(proposal=="true") 
for(N in c(5,6000000)){
  results <- readRDS(paste(homedir,sprintf("MFsampling/final2/mp1N%d_homo_noise6.rds",N),sep = '/'))
  df=rbind(df, results$df)
  trhold = rbind(trhold, results$trhold)
}
p1= ggplot() + 
  geom_line(data= df0%>% 
              dplyr::select(-c(proposal,N,rr, iteration)) ,aes(x=y, y=log10(est_dens)), col='red') + 
  geom_point(data=subset(df,is_support & proposal%in% c("optimal","modified") )%>% mutate(N=ifelse(log10(N) %% 1,"N0=6*10^6", "N0=10^5"),proposal=factor(proposal, levels=c("optimal","modified"))), aes(x=y, y=log10(est_dens), color = "observed"),alpha=0.1,shape='.') +
  geom_point(data=subset(df,!is_support& proposal%in% c("optimal","modified") )%>% mutate(N=ifelse(log10(N) %% 1,"N0=6*10^6", "N0=10^5"),proposal=factor(proposal, levels=c("optimal","modified"))), aes(x=y, y=log10(est_dens), color = "extended"),alpha=0.1,shape='.') +
  facet_wrap(~N+proposal,ncol=4) +ylim(-10,0)+ xlim(-100,100)+
  ylab("Log10 of Estimated Density") + labs(color="")+xlab("Y")+
  geom_vline(data= trhold%>% filter(proposal%in% c("optimal","modified")) %>% mutate(N=ifelse(log10(N) %% 1,"N0=6*10^6", "N0=10^5"),proposal=factor(proposal, levels=c("optimal","modified"))) %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
             aes(xintercept = y, color="25th Y"),linetype="dashed",alpha=0.1) +
  scale_color_manual(
    values=c("observed"="black", "extended"="blue", "25th Y"='green4'), 
    labels=c("observed", "extended",'25th Y'),
    breaks=c("observed", "extended", "25th Y")
  )+
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  linetype = c("blank", "blank", "dashed"), 
                                                  shape = c(19, 19, NA) ))) +
  theme(
    legend.position='top', 
    legend.direction="horizontal",
    text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.background = element_rect(fill="transparent"),
    strip.text.x = element_text(size=12))

fname = "./MFsampling/final2/f2byN.png"
ggsave(filename = fname, plot = p1, width = 10, height = 2.9)

#6.2.2
calculate_var<-function(results, m_idx, noise_sigma){
  df_var1 <- subset(results$df,is_support) %>%
    group_by(proposal, y, N) %>% 
    summarise(
      dens_var=var(est_dens,na.rm=T)
    )
  df_var1 = df_var1 %>% left_join(results$df %>% filter(proposal =="true") %>% dplyr::select(y, est_dens), by="y")
  df_var1=df_var1%>%mutate(dens_var_scaled = dens_var / est_dens^2)
  return(list(df_var= df_var1%>% mutate(noise=ifelse(noise_sigma, "Homoscedastic", "Noiseless"), type=sprintf("m%d",m_idx)),
              trhold = results$trhold %>% mutate(noise=ifelse(noise_sigma, "Homoscedastic", "Noiseless"), type=sprintf("m%d",m_idx)),
              unif = results$df %>% filter(proposal%in% c("true","uniform","modified"))%>% mutate(noise=ifelse(noise_sigma, "Homoscedastic", "Noiseless"), type=sprintf("m%d",m_idx)) ))
}

df_var = NULL
trhold=NULL
unif = NULL

results <- readRDS(paste(homedir,"MFsampling/final2/mexp-homoN6_homo_noise6.rds",sep = '/'))
vlist <- calculate_var(results,3,6)
vlist <- calculate_var(results,3,0)
vlist <- calculate_var(results,1,6)
vlist <- calculate_var(results,1,0)
# vlist <- calculate_var(results,2,6)
# vlist <- calculate_var(results,2,0)
df_var=rbind(df_var, vlist$df_var)
trhold=rbind(trhold, vlist$trhold)
unif = rbind(unif, vlist$unif)

results <- readRDS(paste(homedir,"MFsampling/final2/mexp-homo-honehalfN6_noiseless_noise0.rds",sep = '/'))
results <- readRDS(paste(homedir,"MFsampling/final2/mp1N6000000_homo_noise6.rds",sep = '/'))
results <- readRDS(paste(homedir,"MFsampling/final2/mp1N6_noiseless_noise0.rds",sep = '/'))
# results <- readRDS(paste(homedir,"MFsampling/final/mPMN6_homo_noise6.rds",sep = '/'))
# results <- readRDS(paste(homedir,"MFsampling/final/mPMN6_noiseless_noise0h1.rds",sep = '/'))

f=3
p1= df_var  %>% filter(type ==sprintf("m%d",f)) %>%
  filter(proposal %in% c("random", "uniform", "optimal", "modified")) %>%
  mutate(proposal=factor(proposal, levels=c("random", "uniform", "optimal", "modified")), noise = factor(noise,levels=c("Noiseless","Homoscedastic"))) %>%
  ggplot(aes(x=y, y=log10(dens_var_scaled), group=proposal)) +
  geom_line(aes(linetype=proposal)) +
  ylab("Log10 of Scaled Variance") +
  scale_linetype_manual(
    values=c("dotted", "twodash","longdash", "solid"),
    labels=c("random", "uniform", "optimal", "modified"))+
  #coord_cartesian(ylim=c(-10, 10)) +
  geom_vline(data= trhold %>% filter(proposal %in% c( "modified") & (type ==sprintf("m%d",f))) %>% 
               pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") %>%
               mutate(noise = factor(noise,levels=c("Noiseless","Homoscedastic"))),
             aes(xintercept = y, color="GPD threshold"),linetype="dashed", alpha=0.1) +
  facet_grid(type~noise) +
  scale_color_manual(values=c("green4"), labels=c("GPD threshold")) + 
  theme(
    # legend.position=c(0.5, 0.9), 
    # legend.direction="horizontal",
    text = element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12)
  ) +xlab(ifelse(f==3,"Y",""))+
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  labs(linetype = "", color= "")+ ylim(-5,5)+ #xlim(-80,80)+ #xlim(ifelse(f==1,-100,-10),ifelse(f==1,100,150))+
  theme(
    legend.position=ifelse(f==1,"top","None"), 
    # legend.direction="horizontal",
    text = element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12)
  )
fname = sprintf("./MFsampling/final2/f2Varp%d.png",f)
ggsave(filename = fname, plot = p1, width = 10, height = ifelse(f==1,2.9,2.5))

trhold2 = trhold %>% filter(proposal%in% c("true","uniform","modified"))

f=1
p=ggplot() + 
  geom_line(data= unif %>% filter(type ==sprintf("m%d",f)) %>% filter(proposal=="true")%>% mutate( noise = factor(noise,levels=c("Noiseless","Homoscedastic"))) %>%
              dplyr::select(-c(proposal,N,rr, iteration)) ,aes(x=y, y=log10(est_dens)), col='red') + 
  geom_point(data=subset(unif,is_support & proposal!="true" )%>% filter(type ==sprintf("m%d",f)) %>%mutate(proposal=factor(proposal, levels=c("modified","uniform")), noise = factor(noise,levels=c("Noiseless","Homoscedastic"))), aes(x=y, y=log10(est_dens), color = "observed"),alpha=0.1,shape='.') +
  geom_point(data=subset(unif,!is_support& proposal!="true" )%>% filter(type ==sprintf("m%d",f)) %>%mutate(proposal=factor(proposal, levels=c("modified","uniform")), noise = factor(noise,levels=c("Noiseless","Homoscedastic"))), aes(x=y, y=log10(est_dens), color = "extended"),alpha=0.1,shape='.') +
  facet_grid(type~noise+proposal) +ylim(-10,0)+ xlim(-100,100)+
  ylab("Log10 of Estimated Density") + labs(color="")+xlab(ifelse(f==3,"Y",""))+
  geom_vline(data= trhold2 %>% filter(type ==sprintf("m%d",f)) %>% mutate(proposal=factor(proposal, levels=c("modified","uniform")), noise = factor(noise,levels=c("Noiseless","Homoscedastic"))) %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
             aes(xintercept = y, color="25th Y"),linetype="dashed",alpha=0.1) +
  scale_color_manual(
    values=c("observed"="black", "extended"="blue", "25th Y"='green4'), 
    labels=c("observed", "extended",'25th Y'),
    breaks=c("observed", "extended", "25th Y")
  )+
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  linetype = c("blank", "blank", "dashed"), 
                                                  shape = c(19, 19, NA) ))) +
  theme(
    legend.position=ifelse(f==1,"top","None"), 
    legend.direction="horizontal",
    text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.background = element_rect(fill="transparent"),
    strip.text.x = element_text(size=12))
fname = sprintf("./MFsampling/final2/f2UnifM%d.png",f)
ggsave(filename = fname, plot = p, width = 10, height = ifelse(f==1,2.9,2.5))

## Figure for poster
results <- readRDS(paste(homedir,"MFsampling/final2/mp1N6000000_homo_noise6.rds",sep = '/'))

trhold = results$trhold %>% 
  filter( proposal %in% c("random", "optimal", "modified")) %>%
  mutate(proposal=factor(proposal, levels=c("random", "optimal", "modified"))) %>% mutate(prob="m1")
  #mutate(proposal = factor(case_when(proposal=="random"~"random - N",proposal=="plr"~"PLR - N", proposal=="total"~"random - N0",proposal=="optimal"~"optimal - N", proposal=="modified"~"modified - N", proposal=="uniform"~"uniform - N", proposal=="plr-modified"~"PLR modified - N",.default = proposal), levels=c("random - N","optimal - N", "random - N0", "modified - N", "PLR - N", "uniform - N","PLR modified - N"))) %>%
  
df0 = results$df %>% filter(proposal=="true") %>% mutate(prob="m1")
df = results$df %>% filter( proposal %in% c("random", "optimal", "modified")) %>%
  mutate(proposal=factor(proposal, levels=c("random", "optimal", "modified"))) %>% mutate(prob="m1")

p1= ggplot() + 
  geom_line(data= df0%>% filter( proposal!="uniform - N") %>%
              dplyr::select(-c(proposal,N,rr, iteration)) ,aes(x=y, y=log10(est_dens)), col='red') + 
  #geom_line(aes(x=y, y=log10(est_dens),group=iteration, color=factor(ifelse(is_support,"observed","extended"),levels=c("observed", "extended"))),alpha=0.1) + 
  geom_point(data=subset(results$df,is_support & proposal!="true" & proposal!="uniform - N"), aes(x=y, y=log10(est_dens), color = "observed"),alpha=0.1,shape='.') +
  geom_point(data=subset(results$df,!is_support& proposal!="true" & proposal!="uniform N "), aes(x=y, y=log10(est_dens), color = "extended"),alpha=0.1,shape='.') +
  facet_wrap(~proposal, ncol=3) +ylim(-10,0)+# xlim(-25,25)+
  ylab("Log10 of Estimated Density") + labs(color="")+xlab("Y")+
  geom_vline(data= results$trhold %>% filter( proposal!="uniform - N") %>% pivot_longer(c("yL","yR"),names_to = "trhold", values_to="y") ,
             aes(xintercept = y, color="25th Y"),linetype="dashed",alpha=0.1) +
  scale_color_manual(
    values=c("observed"="black", "extended"="blue", "25th Y"='green4'), 
    labels=c("observed", "extended",'25th Y'),
    breaks=c("observed", "extended", "25th Y")
  )+
  guides(color = guide_legend(override.aes = list(alpha = 1,
                                                  linetype = c("blank", "blank", "dashed"), 
                                                  shape = c(19, 19, NA) 
  ))) +
  theme(
    legend.position=ifelse(m_idx==1,'top','None'), 
    legend.direction="horizontal",
    text=element_text(size=12),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12),
    legend.background = element_rect(fill="transparent"),
    strip.text.x = element_text(size=12)
  )
