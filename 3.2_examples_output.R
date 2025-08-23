#####
rm(list=ls())

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))
source(paste0(path,"func_pathogens.R"))
source("ggplot_default_editing.R")

require(ggpubr)

# step 6: get AUC for each example (also for ancestral strain)
evaluate_mat <- matrix(NA,6,3)

### for ancestral strain first
rt_truth_ancestral <- read.csv("est/truth/sim_truth1.csv")
## to resample population Ct
vl_full_ancestral <- read.csv("est/viral_detect/est_detect_vl_scenario1.csv")
# get daily case count (by sampling time)
ct0_ancestral <- read.csv("est/detected/scenario1.csv")
daily_count <- ct0_ancestral %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
# also keep the training period same as main
peak_date <- as.Date("2022-01-01") + ct0_ancestral %>% group_by(confirmed_time) %>% 
  summarise(n=n()) %>% ### updated: restrict to the first epi wave
  filter(confirmed_time<100) %>%
  filter(n==max(n)) %>% 
  pull(confirmed_time) 
## use the first peak date (if multiple)
train_period <- c(peak_date[1]-21,peak_date[1]+21) 
#
rt0_ancestral <- read.csv("est/rt/rt_obs_scenario1.csv")
#### resampling start here
n_run_time <- 100
auc_tmp <- rep(NA,n_run_time) ## running times
for (i in 1:n_run_time){
  set.seed(i)
  ct_update <- sample_ct_only(vl_full_ancestral,daily_count)
  ## updated estimates
  rt_tmp <- evaluate_daily_funcs(rt0_ancestral,ct_update,rt_truth_ancestral,
                                 seq(train_period[1],train_period[2],1)) 
  # make sure days to compare is the same here (across all scenarios)
  rt_evaluate <- rt_tmp%>%filter(date>=(as.Date("2022-01-01")+110))
  auc_tmp[i] <- Rt_auc(rt_evaluate)[2]
}
## get point and interval estimates
evaluate_mat[1,1:3] <- round(c(median(auc_tmp),quantile(auc_tmp,c(.025,.975))),3)

#### results for the five other pathogens
for (n in 1:5){ # 
  rt_truth <- read.csv(paste0("est/2_examples/truth/sim_truth",n,".csv"))
  start.time <- Sys.time()
  ## to resample population Ct
  vl_full0 <- read.csv(paste0("est/2_examples/viral_detect/est_detect_vl_scenario",n,".csv"))
  # get daily case count (by sampling time)
  ct0 <- read.csv(paste0("est/2_examples/detected/scenario",n,".csv")) 
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
  rt0 <- read.csv(paste0("est/2_examples/rt/rt_obs_scenario",n,".csv"))
  #### resampling start here
  n_run_time <- 100
  auc_tmp <- rep(NA,n_run_time) ## running times
  for (i in 1:n_run_time){
    set.seed(i)
    ct_update <- sample_ct_only(vl_full0,daily_count)
    ## updated estimates
    rt_tmp <- evaluate_daily_funcs(rt0,ct_update,rt_truth,
                                   seq(train_period[1],train_period[2],1)) 
    # make sure days to compare is the same here (across all scenarios)
    rt_evaluate <- rt_tmp%>%filter(date>=(as.Date("2022-01-01")+110))
    auc_tmp[i] <- Rt_auc(rt_evaluate)[2]
  }
  ## get point and interval estimates
  evaluate_mat[n+1,1:3] <- round(c(median(auc_tmp),quantile(auc_tmp,c(.025,.975))),3)
  print(paste0("Completed scenario: ",n,", time = ",Sys.time()-start.time,"."))
} 

evaluate_mat

#write.csv(evaluate_mat,"output/evaluate_real_pathogen.csv",row.names = F)


########## plot performance output and example viral trajectories
vl_traj_list <- NULL
vl_traj_list[[1]] <- read.csv("est/viral_traj/est_vl_scenario1.csv")
ct_list <- NULL
ct_list[[1]] <- read.csv("est/detected/scenario1.csv") 
for (k in 1:5){
  vl_traj_list[[k+1]] <- read.csv(paste0("est/2_examples/viral_traj/est_vl_scenario",k,".csv"))
  ct_list[[k+1]] <- read.csv(paste0("est/2_examples/detected/scenario",k,".csv"))
}

### plot
pathogen_used <- c(paste0("SARS-CoV-2\n(",c("Ancestral","Alpha","Delta","Omicron"),")"),"SARS-CoV-1","Influenza A")
## traj 
x_left <- rep(c(0,0.22),3)
y_top <- rep(c(1,0.67,0.34),each=2)
for (k in 1:6){
  vl_traj_tmp <- vl_traj_list[[k]] %>% filter(i%in%ct_list[[k]]$i)
  all_id <- unique(vl_traj_tmp$i)
  set.seed(k)
  selected_id <- sample(all_id,100,F)
  par(fig=c(x_left[k],x_left[k]+0.2,y_top[k]-0.33,y_top[k]),mar=c(3,3,2.7,1)+.1,new=(k!=1))
  if (k%in%1:5){
    plot(NA,xlim=c(0,30),ylim=rev(c(5,40)),axes=F)
  } else {
    plot(NA,xlim=c(0,10),ylim=rev(c(5,40)),axes=F)
  }
  axis(1)
  axis(2,at=1:4*10,las=1)
  for (n in 1:100){
    vl_traj_indiv <- vl_traj_tmp %>% filter(i==selected_id[n])
    lines(vl_traj_indiv$test.to.infect,vl_traj_indiv$ct.value,col=alpha("gray",.4))
  }
  ##
  ct_indiv <- ct_list[[k]] %>% filter(i%in%selected_id) %>%
    mutate(sampled_delay=onset_time+confirmation_delay-infection_time) 
  time_mat <- matrix(NA,2,3) # onset and sampling
  time_mat[1,] <- quantile(ct_indiv$incu_period,c(.5,.25,.75))
  time_mat[2,] <- quantile(ct_indiv$sampled_delay,c(.5,.25,.75))
  col_used <- c("orange","blue")
  for (t in 1:2){
    lines(time_mat[t,2:3],rep(10-2*t,2),col=col_used[t])
    points(time_mat[t,1],10-2*t,col=col_used[t],pch=16)
  }
  ## median time for viral peak
  if (k %in% 1:2){
    abline(v=time_mat[1,1],lty=2)
  } else {
    ct_indiv <- ct_indiv %>% 
      left_join(.,vl_traj_tmp%>%filter(!duplicated(i))%>%dplyr::select(i,time_peak))
    abline(v=median(ct_indiv$time_peak),lty=2) ## median time for viral peak
  }
  mtext(pathogen_used[k],side=3,font=2,adj=0,line=0+.8*(k%in%5:6))
  ##
  if (k==1){
    mtext("A",side=3,font=2,at=-6.5,cex=1.2,line=1)
  } else if (k==2){ # legend
    lines(c(15,17),rep(9,2),lty=2)
    text(17.5,9,"Viral peak",adj=0)
    #
    lines(c(15,17),rep(12,2),col="orange")
    points(16,12,col="orange",pch=16)
    text(17.5,12,"Onset",adj=0)
    lines(c(15,17),rep(16,2),col="blue")
    points(16,16,col="blue",pch=16)
    text(17.5,16,"Detected",adj=0)
  }
  mtext("Ct value",side=2,cex=.9,line=2)
  mtext("Days since infection",side=1,cex=.9,line=2)
}

par(fig=c(.44,1,0,1),mar=c(3.3,2,1.7,1)+.1,new=T)
plot(NA,xlim=c(0.5,6.5),ylim=c(0,1),axes=F,xlab=NA,ylab=NA)
axis(1,at=1:6,labels=c(paste0("SARS-CoV-2\n",c("Ancestral","Alpha","Delta","Omicron")),
                       "  \nSARS-CoV-1","  \nInfluenza A"),padj=.5)
axis(2,las=1,line=-1.5)
mtext("AUC",side=2,line=.9,cex=.9)
for (k in 1:6){
  lines(rep(k,2),evaluate_mat[k,2:3])
  points(k,evaluate_mat[k,1],pch=16)
}
mtext("B",side=3,font=2,adj = 0,cex=1.2)

# fig_4 (8*15)

################################
################################
#### get correlation box plots
###
pathogen_tab <- c(paste0("SARS-CoV-2\n",c("Ancestral","Alpha","Delta","Omicron")),"SARS-CoV-1","Influenza A")
col_used <- c("#073b4c","#118ab2","#06d6a0","#ffd166","#ef476f","#aaaaaa")
rt_list <- rt_truth_list <- NULL
rt_list[[1]] <- rt0_ancestral
rt_truth_list[[1]] <- rt_truth_ancestral

for (k in 1:5){
  rt_list[[k+1]] <- read.csv(paste0("est/2_examples/rt/rt_obs_scenario",k,".csv"))
  rt_truth_list[[k+1]] <- read.csv(paste0("est/2_examples/truth/sim_truth",k,".csv"))
}

rt_est_list <- NULL
for (k in 1:6){
  ct0 <- ct_list[[k]] %>% 
    mutate(sampled_time=onset_time+confirmation_delay,
           confirmed_time=sampled_time)
  rt0 <- rt_list[[k]]
  rt_truth <- rt_truth_list[[k]]
  #
  daily_count <- ct0 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
  # keep the training period same as main
  peak_date <- as.Date("2022-01-01") + daily_count %>% 
    ### updated: restrict to the first epi wave
    filter(sampled_time<110) %>% filter(n==max(n)) %>% pull(sampled_time) 
  ## use the first peak date (if multiple)
  train_period <- c(peak_date[1]-21,peak_date[1]+21) 
  rt_tmp <- evaluate_daily_funcs(rt0,ct0,rt_truth,
                                 seq(train_period[1],train_period[2],1)) %>%
    mutate(rt_cat=cut(Rt,breaks = c(0,0.5,1,1.5,5))) %>% 
    filter(sampled_time>=110) # only include days to compare
  #
  rt_est_list[[k]] <- rt_tmp
}

rt_est_all <- rbind(rt_est_list[[1]]%>%mutate(scenario="Ancestral"),
                    rt_est_list[[2]]%>%mutate(scenario="Alpha"),
                    rt_est_list[[3]]%>%mutate(scenario="Delta"),
                    rt_est_list[[4]]%>%mutate(scenario="Omicron"),
                    rt_est_list[[5]]%>%mutate(scenario="SARS"),
                    rt_est_list[[6]]%>%mutate(scenario="FluA")) %>%
  mutate(scenario=factor(scenario,levels=c("Ancestral","Alpha","Delta","Omicron","SARS","FluA")))
##
pp1 <- ggplot(data = rt_est_all) +
  geom_boxplot(aes(
    x = rt_cat,
    y = mean,fill=scenario),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  scale_x_discrete(name = 'Simulation truth',
                   expand = c(0.05, 0.05),
                   labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  scale_y_reverse(name="Daily Ct mean") +
  geom_hline(yintercept = 25,linetype = 'dashed',size = .5,color = 'gray77')+
  coord_cartesian(ylim=c(35,20))+
  scale_fill_manual(name = 'Pathogen',values = col_used,labels=pathogen_tab)+
  labs(title="A")

pp2 <- ggplot(data = rt_est_all) +
  geom_boxplot(aes(
    x = rt_cat,
    y = skewness,fill=scenario),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  geom_hline(yintercept = 0,linetype = 'dashed',size = .5,color = 'gray77')+
  scale_x_discrete(name = 'Simulation truth',
                   expand = c(0.05, 0.05),
                   labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  scale_y_continuous(name="Daily Ct skewness") +
  scale_fill_manual(name = 'Pathogen',values = col_used,labels=pathogen_tab)+
  labs(title="B")

pp3 <- ggplot(data = rt_est_all) +
  geom_boxplot(aes(
    x = rt_cat,
    y = fit,fill=scenario),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  geom_hline(yintercept = 1,linetype = 'dashed',size = .5,color = 'gray77')+
  scale_x_discrete(name = 'Simulation truth',
                   expand = c(0.05, 0.05),
                   labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  scale_y_continuous(name="Ct-based Rt") +
  coord_cartesian(ylim=c(0,5))+
  scale_fill_manual(name = 'Pathogen',values = col_used,labels=pathogen_tab)+
  labs(title="C")

ggarrange(pp1,pp2,pp3,ncol=3,nrow=1,common.legend = T,legend = "bottom") ## fig_s7 (6*14)

################################################################################
################################################################################
## additional explanation
vl_new_list <- NULL
for (k in 1:6){
  vl_new_list[[k]] <- vl_traj_list[[k]] %>% filter(i%in%ct_list[[k]]$i) %>%
    group_by(i) %>% summarise(peak_time=test.to.infect[which.min(ct.value)]) %>% ungroup() %>%
    left_join(.,ct_list[[k]]) %>%
    mutate(test.to.peak=test.to.infect-peak_time)
}

vl_new_all <-  rbind(vl_new_list[[1]] %>% mutate(pathogen="ancestral"),
                     vl_new_list[[2]] %>% mutate(pathogen="alpha"),
                     vl_new_list[[3]] %>% mutate(pathogen="delta"),
                     vl_new_list[[4]] %>% mutate(pathogen="omicron"),
                     vl_new_list[[5]] %>% mutate(pathogen="sars"),
                     vl_new_list[[6]] %>% mutate(pathogen="flu")) %>%
  mutate(pathogen=factor(pathogen,levels=c("ancestral","alpha","delta","omicron","sars","flu")))


count1 <- vl_new_all %>% group_by(pathogen,test.to.peak) %>% summarise(n=n()) %>% ungroup()

dd <- ggplot(count1)+geom_bar(aes(x=test.to.peak,y=n,fill=pathogen),stat="identity",position="dodge") + 
  scale_fill_brewer(palette="Set1")+facet_wrap(~pathogen,nrow=6,scales="free_y")+
  coord_cartesian(xlim=c(-10,12))


vl_new_all %>% group_by(pathogen) %>% summarise(lower_pt=quantile(test.to.peak,.025),upper_pt=quantile(test.to.peak,.975)) %>% ungroup()

vl_new_all2 <- vl_new_all %>% 
  filter((pathogen=="ancestral"&test.to.peak%in%-1:11)|
           (pathogen=="alpha"&test.to.peak%in%-1:11)|
           (pathogen=="delta"&test.to.peak%in%0:12)|
           (pathogen=="omicron"&test.to.peak%in%-4:8)|
           (pathogen=="sars"&test.to.peak%in%-10:2)|
           (pathogen=="flu"&test.to.peak%in%-2:10))

aa <- ggplot(vl_new_all2)+
  geom_smooth(aes(x=test.to.peak,y=ct.value,color=pathogen))+
  scale_y_reverse()+scale_color_brewer(palette="Set1")

require(cowplot)
plot_grid(plotlist=list(dd,aa), ncol=1, align='v',rel_heights = c(2,1))


####