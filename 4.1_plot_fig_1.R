#####################
##### conceptual plot

##### Vania Lin 
##### coded on 2024-02-05
##### updated on 2024-08-09

#####################
rm(list=ls())

path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)

require(mgcv)
require(tidyverse)
require(MASS)
require(pROC)
require(shape)

#############
##### panel A
rt_truth <- read.csv("est/truth/sim_truth1.csv")
ct0 <- read.csv("est/detected/scenario1.csv")
rt0 <- read.csv("est/rt/rt_obs_scenario1.csv")

daily_count1 <- ct0 %>% group_by(sampled_time) %>% 
        summarise(n=n(),ct_mean=mean(ct.value),
                  ct_skewness=e1071::skewness(ct.value)) %>% ungroup() 

daily_count2 <- daily_count1 %>% full_join(.,rt_truth%>%dplyr::select(step,Rt),
                                           by=c("sampled_time"="step")) %>%
        left_join(.,rt0%>%mutate(days=as.numeric(as.Date(date)-as.Date("2022-01-01"))) %>% 
                          dplyr::select(days,mean,lower_95,upper_95),by=c("sampled_time"="days")) %>%
        mutate(rt_cat=cut(mean,breaks = c(0,0.5,1,1.5,5))) # cut by incidence-based Rt

### get GAM Ct over days 50-100
ct1 <- ct0 %>% filter(sampled_time%in%50:100) %>% mutate(sampled_time=as.numeric(sampled_time))
gam_tmp <- gam(ct.value~s(sampled_time),data=ct1)
daily_count3 <- daily_count2 %>% filter(sampled_time%in%50:100)
pr <- predict(gam_tmp, newdata = daily_count3, type = "response", se = TRUE) # get predicted response values from GAM
daily_count3$ct.gam <- pr$fit
daily_count3$ct.gam.lwr <- pr$fit - qnorm(0.975) * pr$se.fit
daily_count3$ct.gam.upr <- pr$fit + qnorm(0.975) * pr$se.fit

##
vl_traj0 <- read.csv("est/viral_traj/est_vl_scenario1.csv")
set.seed(1)
selected_id <- sample(unique(vl_traj0$i),100,F)
vl_traj1 <- vl_traj0 %>% filter(i%in%selected_id) %>%
        mutate(sampled_delay=onset_time+confirmation_delay-infection_time,
               onset_delay=onset_time-infection_time) 

#pdf("figure1_18jul.pdf",height = 10, width = 16)
### 1) epidemic curve
par(fig=c(0,0.15,0.75,1),mar=c(2,1,1,1)+.1,cex=.8)
plot(NA,xlim=c(0,20),ylim=rev(c(0,4)),axes=F,xlab=NA,ylab=NA)
text(rep(2,3),1:3,c("Specimen\nsampling","Case detection","Epidemic\ncharacteristics"),adj=0)
mtext("A",side=3,font=2,adj=0,cex=2,line=-1)
#
par(fig=c(0.1,0.4,0.75,1),mar=c(2,2,2,1)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,220),ylim=c(0,350),axes=F,xlab=NA,ylab=NA)
for (n in 1:nrow(daily_count2)){
        polygon(daily_count2$sampled_time[n]+rep(c(-.5,.5),each=2),
                c(0,rep(daily_count2$n[n],2),0),col="orange3",border=F)
}       
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=-.7)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Daily number of cases\n(detected and tested)",side=2,line=1.4,cex=.8)
# legend
polygon(rep(c(10,14),each=2),c(330,rep(340,2),330),col="orange3",border=F)
text(14.5,335,"Detected cases",adj=0)
lines(c(10,14),rep(315,2))
text(14.5,315,"Incidence-based Rt",adj=0)
par(new=T)
plot(NA,xlim=c(0,220),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
abline(h=1,lty=2,col="gray70")
daily_rt <- daily_count2 %>% filter(sampled_time %in% 5:200)
polygon(c(daily_rt$sampled_time,rev(daily_rt$sampled_time)),
        c(daily_rt$lower_95,rev(daily_rt$upper_95)),col=alpha("black",.2),border=F) 
lines(daily_rt$sampled_time,daily_rt$mean) # incidence-based Rt
#lines(daily_rt$sampled_time,daily_rt$Rt,col="red3") # simulation truth
axis(4,las=1,hadj=.7,line=-.7)
mtext("Rt",side=4,line=.5,cex=.8)
mtext("Epidemic trajectory",side=3,font=2)

### 2) 
## viral shedding trajectory
par(fig=c(0,0.15,0.6,0.76),mar=c(2,1,1,1)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,20),ylim=rev(c(0,4)),axes=F,xlab=NA,ylab=NA)
text(rep(2,3),1:3,c("Peak viral load","Shedding duration","Time of viral peak"),adj=0)
#
par(fig=c(0.1,0.25,0.6,0.76),mar=c(2,2,2,2)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,35),ylim=rev(c(10,40)),axes=F,xlab=NA,ylab=NA)
for (n in 1:100){
        traj_tmp <- vl_traj1 %>% filter(i==selected_id[n])
        lines(traj_tmp$test.to.infect,traj_tmp$ct.value,col=alpha("gray",.4))
}
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=-.3)
mtext("Days since infection",side=1,line=1,cex=.8)
mtext("Ct value",side=2,line=1.2,cex=.8)
mtext("Individual viral shedding",side=3,font=2)

## population viral load
par(fig=c(0.25,0.4,0.6,0.76),mar=c(2,2,2,2)+.1,cex=.8,new=T)
plot(NA,xlim=c(50,100),ylim=rev(c(24,28)),axes=F,xlab=NA,ylab=NA)
lines(daily_count3$sampled_time,daily_count3$ct.gam)
lines(daily_count3$sampled_time,daily_count3$ct.gam.lwr,lty=2)
lines(daily_count3$sampled_time,daily_count3$ct.gam.upr,lty=2)
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=0)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Ct value",side=2,line=1.5,cex=.8)
mtext("Population viral load",side=3,font=2)

par(fig=c(0.3,0.5,0.75,0.85),mar=c(2,2,2,1)+0.1,new=T)
plot(NA,xlim=c(0,4),ylim=c(0,2),axes=F,xlab=NA,ylab=NA)
Arrows(2.5,1,2.9,1,lwd=2,arr.width=.04)

### 3) correlation plot
par(fig=c(0.47,0.67,0.8,1),mar=c(2,2,2,1)+0.1,new=T)
boxplot(daily_count2$ct_mean~daily_count2$rt_cat,axes=F,ylim=rev(c(20,30)),
        ylab=NA,xlab=NA,boxwex=.15,at=1:4-0.1,whisklty = 1,outpch=16,outcex=.3,staplecol="white",col="gray77")
axis(1,at=1:4,labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5'),line=-.5,padj = -.3,lwd.tick=.5,lwd=.5)
mtext("Rt",side=1,line=1.2,cex=.8)
mtext("Daily mean Ct",side=2,line=2.4,cex=.8)
axis(2,las=1,line=0,hadj=1.1,lwd.tick=.5,lwd=.5)
mtext("Correlation between Rt and\npopulation viral load",side=3,font=2,line=-.2)

par(fig=c(0.47,0.67,0.6,0.8),mar=c(2,2,2,1)+0.1,new=T)
boxplot(daily_count2$ct_skewness~daily_count2$rt_cat,axes=F,ylim=c(-1,1.5),
        ylab=NA,xlab=NA,boxwex=.15,at=1:4-0.1,whisklty = 1,outpch=16,outcex=.3,staplecol="white",col="gray77")
axis(1,at=1:4,labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5'),line=-.5,padj = -.3,lwd.tick=.5,lwd=.5)
mtext("Rt",side=1,line=1.2,cex=.8)
mtext("Daily Ct skewness",side=2,line=2.4,cex=.8)
axis(2,las=1,line=0,hadj=1.1,lwd.tick=.5,lwd=.5)

par(fig=c(0.53,0.73,0.75,0.85),mar=c(2,2,2,1)+0.1,new=T)
plot(NA,xlim=c(0,4),ylim=c(0,2),axes=F,xlab=NA,ylab=NA)
Arrows(2.5,1,2.9,1,lwd=2,arr.width=.04)

### 4) Rt time-series 
par(fig=c(0.7,1,.91,1),mar=c(2,2,2,1)+0.1,new=T)
plot(NA,xlim=c(0,10),ylim=c(0,1),axes=F,xlab=NA,ylab=NA)
lines(c(0,.5),rep(.5,2),lwd=2)
text(.6,.5,"Incidence-based Rt",adj=0)
lines(c(3.5,4),rep(.5,2),lwd=2,col="pink2")
text(4.1,.5,"Ct-based Rt",adj=0)
lines(c(6,6.5),rep(.5,2),lwd=2,col="red3")
text(6.6,.5,"Simulation truth",adj=0)
mtext("Different Rt streams",side=3,font=2)
##
Rt_auc <- function(df){
  df_used <- df %>% mutate(
    rt_true=1*(Rt>=1),
    rt_pred=1*(fit>=1)
  )
  auc_pred <- df_used %>% roc("rt_true","rt_pred")
  auc_out <- round(as.numeric(auc_pred$auc),2)
  #
  return(auc_out)
}
#
dat_train <- daily_count2 %>% filter(sampled_time%in%50:90)
lm_train <- lm(log(mean)~ct_mean+ct_skewness,data=dat_train)
daily_count4 <- daily_count2 %>% filter(sampled_time%in%50:200)
est <- exp(predict(lm_train,daily_count4,interval = "prediction"))
rt_combined <- cbind(daily_count4,est) #%>% filter(sampled_time>=10) #
summary(rt_combined$fit)
auc_used <- c(Rt_auc(rt_combined%>%filter(sampled_time%in%50:90)),Rt_auc(rt_combined%>%filter(sampled_time%in%110:200)))
par(fig=c(0.7,1,0.6,.95),mar=c(2,2,2,1)+0.1,new=T)
plot(NA,xlim=c(50,220),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
polygon(c(rt_combined$sampled_time,rev(rt_combined$sampled_time)),
        c(rt_combined$lower_95,rev(rt_combined$upper_95)),col=alpha("black",.2),border=F) 
lines(rt_combined$sampled_time,rt_combined$mean) # incidence-based Rt
lines(rt_combined$sampled_time,rt_combined$Rt,col="red3") # simulation truth
polygon(c(rt_combined$sampled_time,rev(rt_combined$sampled_time)),
        c(rt_combined$lwr,rev(rt_combined$upr)),col=alpha("pink2",.2),border=F) 
lines(rt_combined$sampled_time,rt_combined$fit,col="pink2") # Ct-based Rt
abline(h=1,lty=2,col="gray70")
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=.5)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Rt",side=2,line=2,cex=.8)
text(150,3.5,paste0("AUC (training): ",auc_used[1]),adj=0)
text(150,3.2,paste0("AUC (testing): ",auc_used[2]),adj=0)
axis(3,at=c(50,90),labels=rep(NA,2),tck=.02)
mtext("Training period",side=3,line=.2,at=70,cex=.8)
axis(3,at=c(110,200),labels=rep(NA,2),tck=.02)
mtext("Testing period",side=3,line=.2,at=155,cex=.8)

###########################################################################################
###########################################################################################
######### panel B (detection examples)
# 1. sampling fraction
ct_scene1 <- read.csv("est/detected/scenario1.csv")
count_scene1 <- ct_scene1 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup() %>%
  full_join(.,ct_scene1%>%sample_frac(.5)%>%group_by(sampled_time)%>%summarise(n_sampled=n())%>%ungroup()) %>%
  replace_na(list(n=0,n_sampled=0))
par(fig=c(0,0.3,0.28,.56),mar=c(2,3,2,2)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,220),ylim=c(0,350),axes=F,xlab=NA,ylab=NA)
for (n in 1:nrow(count_scene1)){
  polygon(count_scene1$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_scene1$n[n],2),0),col="orange3",border=F)
  polygon(count_scene1$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_scene1$n_sampled[n],2),0),col="skyblue2",border=F)
}       
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=-.7)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Daily number of cases",side=2,line=1.2,cex=.8)
# legend
polygon(rep(c(10,14),each=2),c(330,rep(340,2),330),col="orange3",border=F)
text(14.5,335,"Detected cases",adj=0)
polygon(rep(c(10,14),each=2),c(310,rep(320,2),310),col="skyblue2",border=F)
text(14.5,315,"Tested cases",adj=0)
mtext("Varying sampling fractions",side=3,font=2,line=-.5) ### title
par(new=T)
plot(NA,xlim=c(0,220),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
abline(h=1,lty=2,col="gray70")
lines(rt_truth$step,rt_truth$Rt,col="red3") # simulation truth
axis(4,las=1,hadj=.7,line=-.7)
mtext("Rt",side=4,line=.5,cex=.8)
mtext("B",side=3,font=2,at=-27,cex=2)

# 2. biased detection
ct_scene2 <- read.csv("est/detected/scenario10.csv")
count_scene2 <- ct_scene2 %>% group_by(sampled_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild")) %>% ungroup() 
par(fig=c(.3,.6,0.28,.56),mar=c(2,3,2,2)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,220),ylim=c(0,350),axes=F,xlab=NA,ylab=NA)
for (n in 1:nrow(count_scene2)){
  polygon(count_scene2$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_scene2$n[n],2),0),col="orange3",border=F)
  polygon(count_scene2$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_scene2$n_severe[n],2),0),col="red4",border=F)
}       
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=-.7)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Daily number of cases",side=2,line=1.2,cex=.8)
# legend
polygon(rep(c(80,84),each=2),c(330,rep(340,2),330),col="orange3",border=F)
text(84.5,335,"Detected cases with mild infections",adj=0)
polygon(rep(c(80,84),each=2),c(310,rep(320,2),310),col="red4",border=F)
text(84.5,315,"Detected cases with severe infections",adj=0)
mtext("Biased towards severe cases",side=3,font=2,line=-.5) ### title
par(new=T)
plot(NA,xlim=c(0,220),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
abline(h=1,lty=2,col="gray70")
lines(rt_truth$step,rt_truth$Rt,col="red3") # simulation truth
axis(4,las=1,hadj=.7,line=-.7)
mtext("Rt",side=4,line=.5,cex=.8)

# 3. genuine flat testing period
ct_scene3 <- read.csv("est/detected/scenario3.csv")
rt_truth_alt <- read.csv("est/truth/sim_truth2.csv")
count_scene3 <- ct_scene3 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup() 
par(fig=c(0,.3,0,.28),mar=c(2,3,2,2)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,220),ylim=c(0,350),axes=F,xlab=NA,ylab=NA)
for (n in 1:nrow(count_scene3)){
  polygon(count_scene3$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_scene3$n[n],2),0),col="orange3",border=F)
}       
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=-.7)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Daily number of cases",side=2,line=1.2,cex=.8)
# legend
polygon(rep(c(10,14),each=2),c(330,rep(340,2),330),col="orange3",border=F)
text(14.5,335,"Detected cases",adj=0)
lines(c(10,14),rep(315,2),lwd=2,col="red3")
text(14.5,315,"Simulation truth",adj=0)
mtext("Plateaued second wave",side=3,font=2,line=-.5) ### title
par(new=T)
plot(NA,xlim=c(0,220),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
abline(h=1,lty=2,col="gray70")
lines(rt_truth_alt$step,rt_truth_alt$Rt,col="red3") # simulation truth (alternative; flat epidemic)
axis(4,las=1,hadj=.7,line=-.7)
mtext("Rt",side=4,line=.5,cex=.8)

# 4. detection limit
ct_scene4 <- read.csv("est/detected/scenario2.csv")
count_scene4 <- ct_scene4 %>% group_by(sampled_time) %>% summarise(n_detected=n()) %>% ungroup() 
count_combined4 <- full_join(count_scene1%>%dplyr::select(sampled_time,n),count_scene4)
par(fig=c(0.3,.6,0,.28),mar=c(2,3,2,2)+.1,cex=.8,new=T)
plot(NA,xlim=c(0,220),ylim=c(0,350),axes=F,xlab=NA,ylab=NA)
for (n in 1:nrow(count_combined4)){
  polygon(count_combined4$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_combined4$n[n],2),0),col="gray77",border=F)
  polygon(count_combined4$sampled_time[n]+rep(c(-.5,.5),each=2),
          c(0,rep(count_combined4$n_detected[n],2),0),col="orange3",border=F)
}       
axis(1,padj=-1)
axis(2,las=1,hadj=.7,line=-.7)
mtext("Days since outbreak",side=1,line=1,cex=.8)
mtext("Daily number of cases",side=2,line=1.2,cex=.8)
# legend
polygon(rep(c(10,14),each=2),c(330,rep(340,2),330),col="orange3",border=F)
text(14.5,335,"Detected cases",adj=0)
polygon(rep(c(10,14),each=2),c(310,rep(320,2),310),col="gray77",border=F)
text(14.5,315,"Missed cases",adj=0)
mtext("Limited detection",side=3,font=2,line=-.5) ### title
par(new=T)
plot(NA,xlim=c(0,220),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
abline(h=1,lty=2,col="gray70")
lines(rt_truth$step,rt_truth$Rt,col="red3") # simulation truth
axis(4,las=1,hadj=.7,line=-.7)
mtext("Rt",side=4,line=.5,cex=.8)


###########################################################################################
###########################################################################################
######### panel C (trajectory examples)
#####
# vl_full_list <- NULL
# for (k in 1:4){
#   vl_full_list[[k]] <- read.csv(paste0("est/1_relative_time/viral_traj/est_vl_scenario",k,".csv"))
# }
title_used <- c("Peak at day 5\n(irrespective of onset)","Peak before onset","Peak at onset","Peak after onset")
x_left <- rep(c(.64,.82),2)
y_top <- rep(c(.56,.28),each=2)
##
all_keys <- unique(vl_full_list[[1]]$i) 
### fix same individuals (only their shedding trajectory will be varying)
set.seed(615)
selected_keys <- sample(all_keys,200,T)
for (k in 1:4){ ## did not plot out fifth scenario yet
  dat_used <- vl_full_list[[k]]
  ##
  ## trajectory plot
  par(fig=c(x_left[k],x_left[k]+.18,y_top[k]-.28,y_top[k]),mar=c(2,3,2,2)+.1,cex=.8,new=T)
  plot(NA,xlim=c(0,40),ylim=rev(c(10,40)),xlab=NA,ylab=NA,axes=F)
  axis(1,line=0,padj=-.7)
  axis(2,las=1,line=0,hadj=.6)
  mtext("Days since infection",side=1,line=1.1,cex=.8)
  mtext("Ct value",side=2,cex=.8,line=1.5)
  for (n in 1:200){
    traj_tmp <- dat_used %>% filter(i==selected_keys[n])
    lines(traj_tmp$test.to.infect,traj_tmp$ct.value,col=alpha("gray",.5))
  }
  #abline(v=5,lty=2)
  ## relative time plot
  traj_all_tmp <- dat_used %>% filter(i%in%selected_keys) %>% filter(!duplicated(i))
  abline(v=median(traj_all_tmp$incu_period),col="green4",lty=2)
  abline(v=median(traj_all_tmp$time_peak)+0.1,col="blue2",lty=2) # 0.1 is to jagger overlap
  if (k == 1){
    mtext("C",font=2,side=3,at=-7,cex=2)
    ## legends
    lines(c(20,22),rep(12,2),col="green4",lty=3)
    text(22.5,12,"Onset",adj=0)
    lines(c(20,22),rep(15,2),col="blue2",lty=3)
    text(22.5,15,"Viral peak",adj=0)
  }
  mtext(paste0("(",k,") ",title_used[k]),font=2,side=3,line=0)
}

#dev.off()

# fig 1 (10*16)

####