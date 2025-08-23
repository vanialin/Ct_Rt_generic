#####
### to investigate the impact of case count

##### only using data from scenario 1
## updated 2024-01-29

rm(list=ls())

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))
source("ggplot_default_editing.R")

rt_truth <- read.csv("est/truth/sim_truth1.csv")
vl_full0 <- read.csv("est/viral_detect/est_detect_vl_scenario4.csv") # delay for training and testing
##
ct0 <- read.csv("est/detected/scenario1.csv")
rt0 <- read.csv("est/rt/rt_obs_scenario1.csv")

###################################
sample_n_num <- 1:20*10
evaluate_mat <- matrix(NA,20,6) 
##
# overall results
daily_count <- ct0 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
# also keep the training period same as main
peak_date <- as.Date("2022-01-01") + daily_count %>% 
  ### updated: restrict to the first epi wave
  filter(sampled_time<110) %>% filter(n==max(n)) %>% pull(sampled_time) 
## use the first peak date (if multiple)
train_period <- c(peak_date[1]-21,peak_date[1]+21) 
#
##### results from sampling all cases (as the reference)
## use the first peak date (if multiple)
rt_tmp0 <- evaluate_daily_funcs(rt0,ct0,rt_truth,seq(train_period[1],train_period[2],1)) 
# make sure days to compare is the same here (across all scenarios)
rt_evaluate0 <- rt_tmp0 %>% filter(date>=(as.Date("2022-01-01")+110))
result_all <- Rt_auc(rt_evaluate0)[2] ## maybe can serve as the dashline (reference for comparison)
###
n_run_time <- 100
auc_mat <- matrix(NA,length(sample_n_num),n_run_time) ## running times
##
##### for varying case count
for (k in 1:20){ # 20 sampling count cut-off
  ## if daily count lower than the cutoff -> sample all cases on that day
  remain_all <- daily_count %>% filter(n<sample_n_num[k]) %>% pull(sampled_time) 
  ct_tmp1 <- ct0 %>% filter(sampled_time%in%remain_all)
  all_days <- unique(daily_count$sampled_time)
  ## if daily count higher than the cutoff -> sample N cases on that day
  resample_days <- all_days[!all_days%in%remain_all]
  #### resampling start here
  for (i in 1:n_run_time){
    ### update here
    set.seed(k*i)
    for (t in 1:length(resample_days)){
      if (t==1){
        ct_daily <- ct0 %>% filter(sampled_time==resample_days[t]) 
        ct_tmp2 <- ct_daily[sample(1:nrow(ct_daily),sample_n_num[k],F),]
      } else {
        ct_daily <- ct0 %>% filter(sampled_time==resample_days[t]) 
        ct_tmp2 <- rbind(ct_tmp2,ct_daily[sample(1:nrow(ct_daily),sample_n_num[k],F),])
      }
    }
    # combine
    ct_update <- rbind(ct_tmp1,ct_tmp2)
    ## updated estimates
    rt_tmp <- evaluate_daily_funcs(rt0,ct_update,rt_truth,
                                   seq(train_period[1],train_period[2],1)) 
    # make sure days to compare is the same here (across all scenarios)
    rt_evaluate <- rt_tmp %>% filter(date>=(as.Date("2022-01-01")+110))
    auc_mat[k,i] <- Rt_auc(rt_evaluate)[2]
  }
  ## get point and interval estimates
  evaluate_mat[k,1:3] <- round(c(median(auc_mat[k,]),quantile(auc_mat[k,],c(.025,.975))),3)
}

for (n in 2:20){ # get relative change
  rel_change <- rep(NA,100)
  for (i in 1:100){
    rel_change[i] <- (auc_mat[n,i]-auc_mat[n-1,i])/auc_mat[n-1,i]
  }
  ##
  evaluate_mat[n,4:6] <- round(c(median(rel_change),quantile(rel_change,c(.025,.975)))*100,2)
}

evaluate_mat

###
par(mfrow=c(2,1),mar=c(3,2,2,3)+.1)
plot(NA,xlim=c(0.5,20.5),ylim=c(.4,1),axes=F,xlab=NA,ylab=NA)
polygon(rep(c(9.5,20.5),each=2),c(.35,rep(1.05,2),.35),col=alpha("gray77",.3),border=F)
axis(1,at=1:20,labels = rep(NA,20))
axis(2,las=1,line=-2)
mtext("Area under the curve",side=2,line=1)
lines(c(.5,20.5),rep(result_all,2),lty=2)
for (n in 1:20){
  lines(rep(n,2),evaluate_mat[n,2:3])
  points(n,evaluate_mat[n,1],pch=16)
}
mtext("A",side=3,font=2,adj=0,cex=1.5)
##
plot(NA,xlim=c(0.5,20.5),ylim=c(-20,30),axes=F,xlab=NA,ylab=NA)
polygon(rep(c(9.5,20.5),each=2),c(-25,rep(35,2),-25),col=alpha("gray77",.3),border=F)
axis(2,las=1,line=-1)
mtext("Relative change",side=2,line=1)
axis(1,at=1:20,labels = sample_n_num)
mtext("Number of sampled cases",side=1,line=2)
lines(c(.5,20.5),rep(0,2),lty=2)
lines(c(.5,20.5),rep(2,2),lty=2,col="gray50")
for (nn in 2:20){
  lines(rep(nn,2),evaluate_mat[nn,5:6])
  points(nn,evaluate_mat[nn,4],pch=16)
} 
mtext("B",side=3,font=2,adj=0,cex=1.5)

# fig_s5 (8*10)

####