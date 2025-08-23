########################

## output for relative difference
## only have ONE incidence-based Rt

## Vania Lin
## updated 2024-06-03
########################

rm(list=ls())

require(ggpubr)
require(vioplot)
require(gridExtra)

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))
source(paste0(path,"func_pathogens.R"))

### get result first
detect_cases <- read.csv("est/1_relative_time/detected/scenario0.csv")
daily_count <- detect_cases %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
# keep the training period same as main
peak_date <- as.Date("2022-01-01") + daily_count %>% 
  ### updated: restrict to the first epi wave
  filter(sampled_time<110) %>% filter(n==max(n)) %>% pull(sampled_time) 
## use the first peak date (if multiple)
train_period <- c(peak_date[1]-21,peak_date[1]+21) 
rt0 <- read.csv("est/1_relative_time/rt/rt_obs_scenario0.csv")
#
#### resampling start here
n_run_time <- 100
auc_mat <- matrix(NA,4,n_run_time) ## running times
rt_truth <- read.csv("est/1_relative_time/truth/sim_truth.csv")
vl_full0 <- read.csv("est/1_relative_time/viral_traj/est_vl_scenario1.csv")# only as a reference
vl_list <- NULL
for (k in 1:4){
  vl_list[[k]] <- read.csv(paste0("est/1_relative_time/viral_detect/est_detect_vl_scenario",k,".csv")) 
}
##
for (ii in 1:n_run_time){
  set.seed(ii)
  ct_update_id <- sample_ct_only(vl_full0,daily_count) %>%
    filter(!duplicated(i)) %>% pull(i) # get detected cases at each scenario
  for (n in 1:4){ ## focus on the first four different scenarios
    ct_update <- vl_list[[n]] %>% filter(i%in%ct_update_id) %>%
      mutate(sampled_time=onset_time+confirmation_delay,
             confirmed_time=sampled_time)
    ## updated estimates
    rt_tmp <- evaluate_daily_funcs(rt0,ct_update,rt_truth,
                                   seq(train_period[1],train_period[2],1)) 
    # make sure days to compare is the same here (across all scenarios)
    rt_evaluate <- rt_tmp%>%filter(date>=(as.Date("2022-01-01")+110))
    auc_mat[n,ii] <- Rt_auc(rt_evaluate)[2]
  }
}

evaluate_mat <- matrix(NA,4,3)
for (k in 1:4){
  evaluate_mat[k,] <- quantile(auc_mat[k,],c(.5,.025,.975))
}
evaluate_mat
#
kruskal.test(list(auc_mat[1,],auc_mat[2,],auc_mat[3,],auc_mat[4,])) # p < 0.001

################################################
########## individual-level variation ########## 
################################################
ct_all <- rbind(vl_list[[1]]%>%mutate(scenario=1),
                vl_list[[2]]%>%mutate(scenario=2),
                vl_list[[3]]%>%mutate(scenario=3),
                vl_list[[4]]%>%mutate(scenario=4)) %>%
  mutate(sampled_time=onset_time+confirmation_delay)

### check the time period for increasing/decreasing epidemic (choose 10 days each as example)
traj_list <- NULL
for (k in 1:4){
  traj_list[[k]] <- read.csv(paste0("est/1_relative_time/viral_traj/est_vl_scenario",k,".csv")) %>%
    filter(i%in%vl_list[[k]]$i) %>% 
    mutate(sampled_delay=onset_time+confirmation_delay-infection_time) %>% 
    # only focus on the trajectory *until* the time of detection
    filter(test.to.infect<=sampled_delay)
}

### 
## plot
par(fig=c(0,1,0.8,1),mar=c(2,2,2,1)+.1,cex=.8)
plot(NA,xlim=c(0,220),ylim=c(0,350),axes=F,xlab=NA,ylab=NA)
axis(1,at=0:22*10,padj=-1)
axis(2,las=1,hadj=1,line=-2)
mtext("Case count",side=2,cex=.8,line=.5)
mtext("Days since outbreak",side=1,cex=.8,line=1.2)
for (n in 1:nrow(daily_count)){
  polygon(rep(daily_count$sampled_time[n]+c(-.5,.5),each=2),
          c(0,rep(daily_count$n[n],2),0),col="darkorange",border="white")
}
par(new=T)
plot(NA,xlim=c(0,220),ylim=c(0,3),axes=F,xlab=NA,ylab=NA)
axis(4,las=1,hadj=.4,line=-2)
mtext("Rt",side=4,cex=.8,line=0)
polygon(rep(c(110,120),each=2),c(0,rep(3.3,2),0),col=alpha("blue3",.3),border=F)
polygon(rep(c(125,135),each=2),c(0,rep(3.3,2),0),col=alpha("red3",.3),border=F)
polygon(rep(c(170,180),each=2),c(0,rep(3.3,2),0),col=alpha("purple",.3),border=F)
# only plot till day 220
epi_rt_used <- rt0 %>% mutate(days=as.numeric(as.Date(date)-as.Date("2022-01-01"))) %>% filter(days<=220)
polygon(c(epi_rt_used$days,rev(epi_rt_used$days)),
        c(epi_rt_used$lower_95,rev(epi_rt_used$upper_95)),
        col=alpha("gray44",.3),border=F) 
lines(epi_rt_used$days,epi_rt_used$median)
abline(h=1,lty=2,col="gray33")
mtext("Period 1",side=3,cex=.8,at=115)
mtext("Period 2",side=3,cex=.8,at=130)
mtext("Period 3",side=3,cex=.8,at=175)
mtext("A",side=3,font=2,line=.8,adj=0)
## check case count
sum(daily_count$n[daily_count$sampled_time%in%110:120])
sum(daily_count$n[daily_count$sampled_time%in%125:130])
sum(daily_count$n[daily_count$sampled_time%in%170:180])

### plot the AUC
a <- rep(NA,4)
for (k in 1:4){
  a[k] <- median(auc_mat[k,])
}
par(fig=c(0,0.4,0,0.8),mar=c(2.5,3,2,1)+.1,new=T)
plot(NA,xlim=c(.4,1),ylim=rev(c(.7,4.3)),axes=F,xlab=NA,ylab=NA)
lines(rep(min(a),2),c(.5,4.5),col="gray77",lty=2)
lines(rep(max(a),2),c(.5,4.5),col="gray77",lty=2)
for (k in 1:4){
  ## change into showing median and lower/upper quantiles instead
  lines(quantile(auc_mat[k,],c(.025,.975)),rep(k,2),lwd=1.3)
  points(median(auc_mat[k,]),k,pch=16,cex=1.2)
}
axis(1,padj = -.5)
axis(2,at=1:4,labels=paste0("Pathogen ",1:4),col="white")
mtext("B",font=2,side=3,at=.36,line=.5)
mtext("AUC",side=1,cex=.8,line=1.5)

###
set.seed(111)
random_id_increasing <- sample(unique(ct_all$i[ct_all$sampled_time%in%110:120&ct_all$scenario==1]),200,F)
random_id_peak <- sample(unique(ct_all$i[ct_all$sampled_time%in%125:135&ct_all$scenario==1]),200,F)
random_id_decreasing <- sample(unique(ct_all$i[ct_all$sampled_time%in%170:180&ct_all$scenario==1]),200,F)
id_list <- list(random_id_increasing,random_id_peak,random_id_decreasing)
#
summary_mat <- matrix(NA,4,6)
tag_period <- c("Increasing epidemic\n(period 1)","Epidemic peak\n(period 2)","Decreasing epidemic\n(period 3)")
for (k in 1:4){ # four pathogen examples
  for (j in 1:3){ # three different periods
    par(fig=c(0.4+0.2*(j-1),0.4+0.2*j,0.8-0.2*k,0.8-0.2*(k-1)),mar=c(2.2,3,2,1)+.1,new=T)
    traj_tmp <- traj_list[[k]] %>% filter(i%in%id_list[[j]])
    plot(NA,xlim=c(-5,30),ylim=rev(c(10,40)),axes=F,xlab=NA,ylab=NA)
    axis(1,padj=-1,at=0:6*5)
    axis(2,hadj=.6,las=1)
    ###
    ### detected post viral peak
    traj_post <- traj_tmp %>% filter(sampled_delay>=time_peak)
    for (n in 1:length(unique(traj_post$i))){
      traj_indiv <- traj_post %>% filter(i==unique(traj_post$i)[n])
      lines(traj_indiv$test.to.infect,traj_indiv$ct.value,col=alpha("gray",.5))
    }
    ### for those tested after viral peak, use pink
    traj_pre <- traj_tmp %>% filter(sampled_delay<time_peak) 
    for (nn in 1:length(unique(traj_pre$i))){
      traj_indiv <- traj_pre %>% filter(i==unique(traj_pre$i)[nn])
      lines(traj_indiv$test.to.infect,traj_indiv$ct.value,col=alpha("pink",.3))
    }
    #abline(v=5,lty=2,col="red")
    values_detected <- vl_list[[k]] %>% filter(i%in%id_list[[j]])
    vioplot(values_detected$ct.value,add=T,at=-2,wex=7,axes=F,xlab=NA,ylab=NA,ylim=rev(c(10,40)))
    summary_mat[k,2*(j-1)+1] <- mean(values_detected$ct.value)
    summary_mat[k,2*j] <- e1071::skewness(values_detected$ct.value)
    if (j == 1){
      median_v <- median(values_detected$ct.value)
      mtext("Ct value",side=2,line=1.6,cex=.8)
      ### panel legend
      mtext(LETTERS[k+2],side=3,font=2,line=.5,adj=0)
    }
    if (k == 4){
      mtext("Days since infection",side=1,line=1.3,cex=.8,at=15)
    } else if (k == 1){
      mtext(tag_period[j],side=3,cex=.8,font=2,line=-.3)
    }
    abline(h=median_v,lty=2)
    #print(skewness(values_detected$ct.value))
  }
} 

round(summary_mat,2)

# fig_3 (9*12)

#save.image(file = "data/relative_time.RData")

####