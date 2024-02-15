

File_LM<-"L2_ONRFL_"
Case_LM<-"0020"
Mot_LM<-read.table(paste0(File_LM,Case_LM,".mot"), skip = 2, header=TRUE)

T1 <-1
T2 <-5000
# roll
plot(Mot_LM$Time[T1:T2],Mot_LM$Rot_X[T1:T2],type="l")


File_SC<-"flh_irreg6b-th-"
Case_SC<-"0000020"
Mot_SC<-read.table(paste0(File_SC,Case_SC,".mot"), skip = 2, header=TRUE)
names(Mot_SC)[1] <- "Time"
names(Mot_SC)[2] <- "Xcg"
names(Mot_SC)[3] <- "Ycg"
names(Mot_SC)[4] <- "Zcg"
names(Mot_SC)[5] <- "Rot_X"
names(Mot_SC)[6] <- "Rot_Y"
names(Mot_SC)[7] <- "Rot_Z"

# roll
plot(Mot_SC$Time[T1:T2],Mot_SC$Rot_X[T1:T2],type="l")



File_LSTM<-"lstm_output_for_flh_irreg6b-th-"
Case_LSTM<-"0000020"
Mot_LSTM<-read.table(paste0(File_LSTM,Case_LSTM,".txt"), skip = 1, header=TRUE)
names(Mot_LSTM)[1] <- "Time"
names(Mot_LSTM)[2] <- "Zcg"
names(Mot_LSTM)[3] <- "Rot_X"
names(Mot_LSTM)[4] <- "Rot_Y"



T1 <-500
T2 <-1500

plot(Mot_LM$Time[T1:T2],Mot_LM$Rot_X[T1:T2],type="l")
lines(Mot_SC$Time[T1:T2],Mot_SC$Rot_X[T1:T2],col="red")
lines(Mot_LSTM$Time[T1:T2],Mot_LSTM$Rot_X[T1:T2],col="blue")

plot(Mot_LM$Time[T1:T2],Mot_LM$Rot_Y[T1:T2],type="l")
lines(Mot_SC$Time[T1:T2],Mot_SC$Rot_Y[T1:T2],col="red")
lines(Mot_LSTM$Time[T1:T2],Mot_LSTM$Rot_Y[T1:T2],col="blue")

plot(Mot_LM$Time[T1:T2],Mot_LM$Zcg[T1:T2],type="l")
lines(Mot_SC$Time[T1:T2],Mot_SC$Zcg[T1:T2],col="red")
lines(Mot_LSTM$Time[T1:T2],Mot_LSTM$Zcg[T1:T2],col="blue")




library(ggplot2)
library(reshape2)

Tm1 = 500
Tm2 = 1500
Mot_LM_SC = data.frame(Mot_LM$Time[Tm1:Tm2],Mot_LM$Rot_X[Tm1:Tm2],Mot_SC$Rot_X[Tm1:Tm2],Mot_LSTM$Rot_X[Tm1:Tm2])
names(Mot_LM_SC)[1] <- "time"
names(Mot_LM_SC)[2] <- "LAMP"
names(Mot_LM_SC)[3] <- "SC"
names(Mot_LM_SC)[4] <- "LSTM"

Mot_LM_SC_Melted <- reshape2::melt(Mot_LM_SC, id.var='time')

ggplot(Mot_LM_SC_Melted, aes(x=time, y=value, col=variable)) + 
  geom_line(aes(linetype = variable)) +
  labs(x = "time (s)", y = "Roll (deg)") +
  scale_color_manual(values=c("black", "black","black")) + 
  scale_linetype_manual(values=c("solid", "dotted", "dashed")) +
  theme(legend.position = c(0.975, 0.975),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) 


###

R <- 20
max_all <- matrix(0,R,3)

library(stringr)

for (r in 1:R){
  
  print(r)
  File_LM<-"L2_ONRFL_"
  Case_LM<-str_pad(r, 4, pad = "0")
  Mot_LM<-read.table(paste0(File_LM,Case_LM,".mot"), skip = 2, header=TRUE)
  
  File_SC<-"flh_irreg6b-th-"
  Case_SC<-str_pad(r, 7, pad = "0")
  Mot_SC<-read.table(paste0(File_SC,Case_SC,".mot"), skip = 2, header=TRUE)
  names(Mot_SC)[1] <- "Time"
  names(Mot_SC)[2] <- "Xcg"
  names(Mot_SC)[3] <- "Ycg"
  names(Mot_SC)[4] <- "Zcg"
  names(Mot_SC)[5] <- "Rot_X"
  names(Mot_SC)[6] <- "Rot_Y"
  names(Mot_SC)[7] <- "Rot_Z"
  
  
  File_LSTM<-"lstm_output_for_flh_irreg6b-th-"
  Case_LSTM<-str_pad(r, 7, pad = "0")
  Mot_LSTM<-read.table(paste0(File_LSTM,Case_LSTM,".txt"), skip = 1, header=TRUE)
  names(Mot_LSTM)[1] <- "Time"
  names(Mot_LSTM)[2] <- "Zcg"
  names(Mot_LSTM)[3] <- "Rot_X"
  names(Mot_LSTM)[4] <- "Rot_Y"
  
  max_all[r,] <- c(max(Mot_LM$Rot_X),max(Mot_SC$Rot_X),max(Mot_LSTM$Rot_X))
  #max_all[r,] <- c(max(Mot_LM$Rot_Y),max(Mot_SC$Rot_Y),max(Mot_LSTM$Rot_Y))
  
}



plot(max_all[,1],max_all[,2],xlim=c(34,60),ylim=c(34,60))
points(max_all[,1],max_all[,3],col="red")
abline(a=0,b=1)

cor(max_all[,1],max_all[,2])
cor(max_all[,1],max_all[,3])



mR = data.frame("LAMP"=max_all[,1],"SC"=max_all[,2],"LSTM"=max_all[,3])
mR_Melted <- reshape2::melt(mR, id.var='LAMP')


ggplot(mR_Melted,aes(x=value,y=LAMP,shape=variable))+
  geom_point(size=3)+
  scale_shape_manual(values=c(1,16))+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)+
  labs(x = "SC/LSTM", y = "LAMP") +
  xlim(34, 60) +
  ylim(34, 60) +
  theme(legend.position = c(0.225, 0.975),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) 



###



Rtb <- c(
  109,
  445,
  554,
  568,
  904,
  983,
  1068,
  1215,
  1227,
  1264,
  1302,
  1315,
  1358,
  1422,
  1826,
  1919,
  1938,
  1952,
  1957,
  1997
)
max_all_tb <- matrix(0,length(Rtb),3)

library(stringr)

rr = 1
for (r in Rtb){
  
  print(r)
  File_LM<-"L2_ONRFL_"
  Case_LM<-str_pad(r, 4, pad = "0")
  Mot_LM<-read.table(paste0(File_LM,Case_LM,".mot"), skip = 2, header=TRUE)
  
  File_SC<-"flh_irreg6b-th-"
  Case_SC<-str_pad(r, 7, pad = "0")
  Mot_SC<-read.table(paste0(File_SC,Case_SC,".mot"), skip = 2, header=TRUE)
  names(Mot_SC)[1] <- "Time"
  names(Mot_SC)[2] <- "Xcg"
  names(Mot_SC)[3] <- "Ycg"
  names(Mot_SC)[4] <- "Zcg"
  names(Mot_SC)[5] <- "Rot_X"
  names(Mot_SC)[6] <- "Rot_Y"
  names(Mot_SC)[7] <- "Rot_Z"
  
  
  File_LSTM<-"lstm_output_for_flh_irreg6b-th-"
  Case_LSTM<-str_pad(r, 7, pad = "0")
  Mot_LSTM<-read.table(paste0(File_LSTM,Case_LSTM,".txt"), skip = 1, header=TRUE)
  names(Mot_LSTM)[1] <- "Time"
  names(Mot_LSTM)[2] <- "Zcg"
  names(Mot_LSTM)[3] <- "Rot_X"
  names(Mot_LSTM)[4] <- "Rot_Y"
  
  max_all_tb[rr,] <- c(max(Mot_LM$Rot_X),max(Mot_SC$Rot_X),max(Mot_LSTM$Rot_X))
  #max_all_tb[rr,] <- c(max(Mot_LM$Rot_Y),max(Mot_SC$Rot_Y),max(Mot_LSTM$Rot_Y))
  rr = rr+1
  
}


Rtt <- c(
  0037,
  0306,
  0347,
  0414,
  0500,
  0537,
  0637,
  0767,
  0889,
  0973,
  1081,
  1140,
  1402,
  1417,
  1527,
  1537,
  1666,
  1741,
  1794,
  1879
)
max_all_tt <- matrix(0,length(Rtt),3)

library(stringr)

rr = 1
for (r in Rtt){
  
  print(r)
  File_LM<-"L2_ONRFL_"
  Case_LM<-str_pad(r, 4, pad = "0")
  Mot_LM<-read.table(paste0(File_LM,Case_LM,".mot"), skip = 2, header=TRUE)
  
  File_SC<-"flh_irreg6b-th-"
  Case_SC<-str_pad(r, 7, pad = "0")
  Mot_SC<-read.table(paste0(File_SC,Case_SC,".mot"), skip = 2, header=TRUE)
  names(Mot_SC)[1] <- "Time"
  names(Mot_SC)[2] <- "Xcg"
  names(Mot_SC)[3] <- "Ycg"
  names(Mot_SC)[4] <- "Zcg"
  names(Mot_SC)[5] <- "Rot_X"
  names(Mot_SC)[6] <- "Rot_Y"
  names(Mot_SC)[7] <- "Rot_Z"
  
  
  File_LSTM<-"lstm_output_for_flh_irreg6b-th-"
  Case_LSTM<-str_pad(r, 7, pad = "0")
  Mot_LSTM<-read.table(paste0(File_LSTM,Case_LSTM,".txt"), skip = 1, header=TRUE)
  names(Mot_LSTM)[1] <- "Time"
  names(Mot_LSTM)[2] <- "Zcg"
  names(Mot_LSTM)[3] <- "Rot_X"
  names(Mot_LSTM)[4] <- "Rot_Y"
  
  max_all_tt[rr,] <- c(max(Mot_LM$Rot_X),max(Mot_SC$Rot_X),max(Mot_LSTM$Rot_X))
  #max_all_tt[rr,] <- c(max(Mot_LM$Rot_Y),max(Mot_SC$Rot_Y),max(Mot_LSTM$Rot_Y))
  rr = rr+1
  
}




mR1 = data.frame("LAMP"=max_all[,1],"SC:rd"=max_all[,2],"LSTM:rd"=max_all[,3])
mR1_Melted <- reshape2::melt(mR1, id.var='LAMP')

mR2 = data.frame("LAMP"=max_all_tb[,1],"SC:tb"=max_all_tb[,2],"LSTM:tb"=max_all_tb[,3])
mR2_Melted <- reshape2::melt(mR2, id.var='LAMP')

mR3 = data.frame("LAMP"=max_all_tt[,1],"SC:tt"=max_all_tt[,2],"LSTM:tt"=max_all_tt[,3])
mR3_Melted <- reshape2::melt(mR3, id.var='LAMP')

mR123_Melted <- rbind(mR1_Melted,mR2_Melted,mR3_Melted)


ggplot(mR123_Melted,aes(x=value,y=LAMP,shape=variable))+
  geom_point(size=3)+
  scale_shape_manual(values=c(1,16,2,17,0,15))+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)+
  labs(x = "SC/LSTM", y = "LAMP") +
  xlim(26, 100) +
  ylim(26, 70) +
  theme(legend.position = c(0.225, 0.975),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) 



xa <- c(max_all[,1],max_all_tb[,1],max_all_tt[,1])
ya <- c(max_all[,3],max_all_tb[,3],max_all_tt[,3])
cor(xa,ya)
# not too clear if correlation is meaningful

###


load("all_heave.Rda")
All

library(ggplot2)

ggplot(All,aes(x=SC,y=LAMP,shape=sampling))+
  geom_point(size=3)+
  #ylim(3.5*10^8,1.2*10^9)+
  scale_shape_manual(values=c(1,2,0))+
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)+
  labs(x = "SC heave", y = "LAMP heave") +
  theme(legend.position = c(0.225, 0.975),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) 
  #ggtitle("Heave maxima over random, top and bottom SC records") +
  #theme(plot.title = element_text(hjust = 0.5))




###



library(ks)
library(ggplot2)
library(reshape2)



y_x <- function(x,rho){
  rnorm(length(x),mean=rho*x,sd=sqrt(1-rho^2))
}


rho <- 0.8
#rho <- 0

set.seed(108)
x <- runif(300,-5,5)
y <- y_x(x,rho)
y0 <- rnorm(300,0,1)

w <- dnorm(x)
w <- length(x)*w/sum(w)

fhat <- kde(x=y,h=0.4,w=w)
fhat00 <- kde(x=y0,eval.points=fhat$eval.points,w=rep(1,length(y0)))

PR <- data.frame("Y"=fhat$eval.points,"kde"=log10(fhat$estimate),"true"=log10(dnorm(fhat$eval.points)))
PR_Melted <- reshape2::melt(PR, id.var='Y')

PRall <- data.frame("Y"=fhat$eval.points,"actual"=log10(fhat00$estimate),"uniform"=log10(fhat$estimate),"true"=log10(dnorm(fhat$eval.points)))
PRall_Melted <- reshape2::melt(PRall, id.var='Y')

ggplot(PRall_Melted, aes(x=Y, y=value)) + 
  geom_line(aes(linetype = variable)) +
  labs(x = "Y", y = "log10 PDF") +
  scale_linetype_manual(values=c("dotted", "dashed","solid")) +
  scale_color_manual(values=c("black", "black","black")) + 
  theme(legend.position = c(0.675, 0.4),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) 
  #ggtitle("Random and Uniform Samplings") +
  #theme(plot.title = element_text(hjust = 0.5))


