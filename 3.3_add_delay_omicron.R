#####
rm(list=ls())

path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))
source(paste0(path,"func_pathogens.R"))
source("ggplot_default_editing.R")

### only focus on Omicron scenario


rt_truth <- read.csv("est/2_examples/truth/sim_truth3.csv")
start.time <- Sys.time()
## to resample population Ct
vl_traj0 <- read.csv("est/2_examples/viral_traj/est_vl_scenario3.csv")
# get daily case count (by sampling time)
ct0 <- read.csv("est/2_examples/detected/scenario3.csv")
vl_full0 <- read.csv("est/2_examples/viral_detect/est_detect_vl_scenario3.csv")
daily_count <- ct0 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
# also keep the training period same as main
peak_date <- as.Date("2022-01-01") + ct0 %>% group_by(confirmed_time) %>% 
  summarise(n=n()) %>% ### updated: restrict to the first epi wave
  filter(confirmed_time<100) %>%
  filter(n==max(n)) %>% 
  pull(confirmed_time) 
## use the first peak date (if multiple)
train_period <- c(peak_date[1]-21,peak_date[1]+21) 
#
rt0 <- read.csv("est/2_examples/rt/rt_obs_scenario3.csv")


## instead of changing the onset-to-detection delay
## we can change the detection-to-testing delay
## so that it won't affect the incidence-based Rt estimation
set.seed(823)
vl_traj1 <- vl_traj0 %>% group_by(i) %>%
  mutate(detect.to.test=round(runif(1,2,4)),
         delay=floor(incu_period)+confirmation_delay+detect.to.test) %>%
  filter(test.to.infect==delay) %>% ungroup()

### compare old and new detected ct (And delay)
vl_full1 <- vl_full0 %>% left_join(.,vl_traj1 %>% 
                                     ## also attach the time_peak to see 
                                     ## if the delayed detection fall on the declining phase
                                     dplyr::select(i,time_peak,test.to.infect,ct.value),by="i") %>%
  ## relative time in relation to peak
  mutate(test.to.peak.x=test.to.infect.x-time_peak,
         test.to.peak.y=test.to.infect.y-time_peak)
summary(vl_full1$test.to.peak.x)
summary(vl_full1$test.to.peak.y) # fewer were on the negative side

# par(mfrow=c(1,2))
# plot(vl_full1$test.to.peak.x,type="p",col=alpha("gray",.5),ylim=c(-5,25))
# abline(h=0,lty=2)
# plot(vl_full1$test.to.peak.y,type="p",col=alpha("gray",.5),ylim=c(-5,25))
# abline(h=0,lty=2)

##############################
##############################
### update the distribution (for detected Ct)
vl_full_update1 <- vl_full0 %>% 
  left_join(.,vl_traj1 %>% dplyr::select(i,detect.to.test,ct.value),by="i") %>%
  ### update the sampling time (added sampling delay)
  mutate(test.to.infect=test.to.infect+detect.to.test) %>%
  dplyr::select(-ct.value.x) %>%
  rename(ct.value=ct.value.y) 

### check update Ct distribution over time
ct_compare <- ct0 %>% 
  left_join(.,vl_full_update1%>%dplyr::select(i,test.to.infect,ct.value),by="i") %>%
  mutate(sampled_time_y=sampled_time+test.to.infect.y-test.to.infect.x)

ggplot(data=ct_compare)+
  geom_smooth(aes(x=sampled_time,y=ct.value.x),col="red",method="gam")+
  geom_smooth(aes(x=sampled_time_y,y=ct.value.y),col="blue",method="gam")+
  geom_vline(xintercept = 110,linetype = 'dashed',size = .3,color = 'grey') # larger variation for extended testing

####
##
ct_old <- merge_Ct_Rt(rt0,ct0)
ct_new <- merge_Ct_Rt(rt0,ct_compare%>%dplyr::select(-sampled_time)%>%
                        rename(sampled_time=sampled_time_y,ct.value=ct.value.y))

ggplot()+
  geom_line(data=ct_old,aes(x=sampled_time,y=mean),col="red")+
  geom_line(data=ct_new,aes(x=sampled_time,y=mean),col="blue")+
  scale_y_reverse()


### check rt estimate based on updated Ct
##
rt_old <- evaluate_daily_funcs(rt0,ct0,rt_truth,seq(train_period[1],train_period[2],1))
rt_new <- evaluate_daily_funcs(rt0,ct_compare%>%dplyr::select(-sampled_time)%>%
                                 rename(sampled_time=sampled_time_y,ct.value=ct.value.y),
                               rt_truth,seq(train_period[1],train_period[2],1))

ggplot(data=rt_old)+
  geom_line(aes(x=sampled_time,y=Rt))+
  geom_line(aes(x=sampled_time,y=fit),col="red")+
  geom_line(data=rt_new,aes(x=sampled_time,y=fit),col="blue")+
  geom_hline(yintercept = 1,linetype = 'dashed',size = .3,color = 'grey')+
  scale_y_continuous(limits = c(0,5))

Rt_auc(rt_old%>%filter(date>=(as.Date("2022-01-01")+110)))[2]  
Rt_auc(rt_new%>%filter(date>=(as.Date("2022-01-01")+110)))[2]  # check main output


###################
###################
### finalize comparison
omicron_mat <- matrix(NA,2,3)
n_run_time <- 100
auc_tmp1 <- auc_tmp2 <- rep(NA,n_run_time) ## running times
for (i in 1:n_run_time){
  set.seed(i)
  ct_update_old <- sample_ct_only(vl_full0,daily_count)
  ## updated estimates
  rt_tmp1 <- evaluate_daily_funcs(rt0,ct_update_old,rt_truth,
                                 seq(train_period[1],train_period[2],1)) 
  rt_evaluate1 <- rt_tmp1 %>% filter(date>=(as.Date("2022-01-01")+110))
  auc_tmp1[i] <- Rt_auc(rt_evaluate1)[2]
  ###
  daily_count_used <- daily_count 
  ### the 'sampled_time ' has already been updated to include the testing delay
  ct_update <- sample_ct_only(vl_full_update1,daily_count_used,detect_delay = T) 
  ## updated estimates
  rt_tmp2 <- evaluate_daily_funcs(rt0,ct_update,rt_truth,
                                 seq(train_period[1],train_period[2],1)) 
  rt_evaluate2 <- rt_tmp2 %>% filter(date>=(as.Date("2022-01-01")+110))
  auc_tmp2[i] <- Rt_auc(rt_evaluate2)[2]
}

## get point and interval estimates
omicron_mat[1,1:3] <- round(c(median(auc_tmp1),quantile(auc_tmp1,c(.025,.975))),2)
omicron_mat[2,1:3] <- round(c(median(auc_tmp2),quantile(auc_tmp2,c(.025,.975))),2)

omicron_mat

sum(auc_tmp2-auc_tmp1>0)/length(auc_tmp1)*100 # 96%

### check detection and sampling
count_sample <- ct_update %>% group_by(sampled_time) %>% summarize(n=n()) %>% ungroup()
count_case <- ct_update %>% group_by(confirmed_time) %>% summarize(n=n()) %>% ungroup()
## plot sample count and confirm case count over time in a barplot
ggplot(data=count_sample)+
  geom_bar(aes(x=sampled_time,y=n),stat="identity",fill="blue",alpha=.5)+
  geom_bar(data=count_case,aes(x=confirmed_time,y=n),stat="identity",fill="red",alpha=.5)+
  geom_bar(data=daily_count,aes(x=sampled_time,y=n),stat="identity",fill="orange",alpha=.5)+
  scale_y_continuous(sec.axis = sec_axis(~./max(count_sample$n),name="Confirmed case count"))+
  scale_x_continuous(breaks = seq(0,200,20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##
####
set.seed(1)
selected_id <- sample(ct0$i,100,F)
par(fig=c(0,1,0.3,1),mar=c(3,3,2,1)+.1)
plot(NA,xlim=c(0,30),ylim=rev(c(5,40)),axes=F)
axis(1)
mtext("Days since infection",side=1,line=2)
axis(2,las=1,line=-.5)
mtext("Ct value",side=2,line=2)
for (n in 1:200){
  vl_traj_indiv <- vl_traj0 %>% filter(i==selected_id[n])
  lines(vl_traj_indiv$test.to.infect,vl_traj_indiv$ct.value,col=alpha("gray",.4))
}
vl_indiv <- vl_full1 %>% filter(i%in%selected_id)
time_mat <- matrix(NA,3,3) # onset and sampling (before and after)
time_mat[1,] <- quantile(vl_indiv$incu_period,c(.5,.25,.75))
time_mat[2,] <- quantile(vl_indiv$test.to.infect.x,c(.5,.25,.75))
time_mat[3,] <- quantile(vl_indiv$test.to.infect.y,c(.5,.25,.75))
col_used <- c("orange","blue","lightblue3")
text_used <- c("Onset","Original Ct testing","Later Ct testing")
for (t in 1:3){
  ## median and IQR
  lines(time_mat[t,2:3],rep(10-2*t,2),col=col_used[t])
  points(time_mat[t,1],10-2*t,col=col_used[t],pch=16)
  ## legend
  lines(c(19,20),rep(10+2*t,2),col=col_used[t])
  points(19.5,10+2*t,pch=16,col=col_used[t])
  text(20.5,10+2*t,text_used[t],adj=0)
}
vl_flat <- vl_traj0 %>% filter(i%in%selected_id,!duplicated(i))
abline(v=median(vl_flat$time_peak),lty=2)
mtext("A",side=3,font=2,cex=1.2,adj=0)

par(fig=c(0,1,0,0.25),mar=c(3,2,2,1)+.1,new=T)
plot(NA,xlim=c(0.5,1.1),ylim=rev(c(0.5,2.5)),axes=F,xlab=NA,ylab=NA)
axis(1,at=5:10*0.01)
mtext("AUC",side=1,line=.1)
axis(2,at=1:2,labels = paste0(c("Original","Later")," Ct testing"),line=-6,las=1,col="white")
for (k in 1:2){
  lines(omicron_mat[k,2:3],rep(k,2),col=col_used[k+1])
  points(omicron_mat[k,1],k,pch=16,col=col_used[k+1])
  text(1,k,paste0(format(omicron_mat[k,1],nsmall=2),"(",
                  format(omicron_mat[k,2],nsmall=2),",",
                  format(omicron_mat[k,3],nsmall=2),")"),adj=0)
}
mtext("B",side=3,font=2,cex=1.2,adj=0)

# fig_s9 (7*8)

####