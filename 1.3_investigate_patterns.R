## CONT'D from 1.2

#############################################
### plot scenario summary (correlation plot)
func_synthesize_scenario <- function(data,scenario="usual",frac){
  truth_tmp <- data[[1]]
  rt_tmp <- data[[2]]
  ct_tmp <- data[[3]] ## detected linelist
  if (scenario=="fractional"){
    ct_tmp1 <- ct_tmp %>% sample_frac(frac)
  } else {
    ct_tmp1 <- ct_tmp
    frac <- NULL
  }
  count_tmp <- ct_tmp1 %>% group_by(confirmed_time) %>% 
    summarise(n=n()) %>% ungroup()
  ######
  peak_date <- as.Date("2022-01-01") + count_tmp %>% filter(confirmed_time<110) %>% 
    filter(n==max(n)) %>% pull(confirmed_time)
  train_period <- c(peak_date[1]-21,peak_date[1]+21)
  rt_out <- evaluate_daily_funcs(rt_tmp,ct_tmp1,truth_tmp,seq(train_period[1],train_period[2],1)) 
  ###
  ###
  return(rt_out)
}

########### three groups of supplementary panels
truth_list <- rt_list <- ct_list <- NULL
truth_seq <- c(rep(1,5),2,rep(1,8))
rt_seq <- c(rep(1,4),2:11)
for (k in 1:14){
  truth_list[[k]] <- read.csv(paste0("est/truth/sim_truth",truth_seq[k],".csv"))
  rt_list[[k]] <- read.csv(paste0("est/rt/rt_obs_scenario",rt_seq[k],".csv"))
  ct_list[[k]] <- read.csv(paste0("est/detected/scenario",rt_seq[k],".csv")) 
}
##
### scenarios 1-4 (sampling frac)
out_list <- NULL
set.seed(1220)
out_list[[1]] <- func_synthesize_scenario(list(truth_list[[1]],rt_list[[1]],ct_list[[1]]),scenario = "usual")
mtext("Scenario 1",side=3,font=2,cex=.9,adj=0)
out_list[[2]] <- func_synthesize_scenario(list(truth_list[[2]],rt_list[[2]],ct_list[[2]]),scenario = "fractional",frac=.8)
mtext("Scenario 2",side=3,font=2,cex=.9,adj=0)
out_list[[3]] <- func_synthesize_scenario(list(truth_list[[3]],rt_list[[3]],ct_list[[3]]),scenario = "fractional",frac=.5)
mtext("Scenario 3",side=3,font=2,cex=.9,adj=0)
out_list[[4]] <- func_synthesize_scenario(list(truth_list[[4]],rt_list[[4]],ct_list[[4]]),scenario = "fractional",frac=.3)

### scenarios 5-6 (flat testing) and scenarios 7-14 (biased testing, with or without delay)
for (i in 1:10){
  out_list[[4+i]] <- func_synthesize_scenario(list(truth_list[[4+i]],rt_list[[4+i]],ct_list[[4+i]]),scenario="usual")
}

################### 
### try to show boxplot for testing periods 
### (restrict to comparison on/after day 110 to reduce uncertainty across scenarios)
### version 1: boxplot (among Rt categories); version 2: residual plot (against Rt)

### sampling fraction
rt_merged1 <- rbind(out_list[[1]]%>%mutate(scenario=1),out_list[[2]]%>%mutate(scenario=2),
                    out_list[[3]]%>%mutate(scenario=3),out_list[[4]]%>%mutate(scenario=4)) %>% 
  ## also make sure all were testing periods
  filter(sampled_time>=110) %>% 
  mutate(rt_cat=cut(Rt,breaks = c(0,0.5,1,1.5,5),labels = 1:4),
         rt_residual=Rt-fit)

### flat epidemic
rt_merged2 <- rbind(out_list[[5]]%>%mutate(scenario=5),out_list[[6]]%>%mutate(scenario=6)) %>% 
  ## also make sure all were testing periods
  filter(sampled_time>=110) %>% 
  mutate(rt_cat=cut(Rt,breaks = c(0,0.5,1,1.5,5),labels = 1:4),
         rt_residual=Rt-fit)

### biased detection
rt_merged3 <- rbind(out_list[[7]]%>%mutate(scenario=7),out_list[[8]]%>%mutate(scenario=8),
                    out_list[[9]]%>%mutate(scenario=9),out_list[[10]]%>%mutate(scenario=10)) %>% 
  ## also make sure all were testing periods
  filter(sampled_time>=110) %>% 
  mutate(rt_cat=cut(Rt,breaks = c(0,0.5,1,1.5,5),labels = 1:4),
         rt_residual=Rt-fit)

### biased detection with delay
rt_merged4 <- rbind(out_list[[11]]%>%mutate(scenario=11),out_list[[12]]%>%mutate(scenario=12),
                    out_list[[13]]%>%mutate(scenario=13),out_list[[14]]%>%mutate(scenario=14)) %>% 
  ## also make sure all were testing periods
  filter(sampled_time>=110) %>% 
  mutate(rt_cat=cut(Rt,breaks = c(0,0.5,1,1.5,5),labels = 1:4),
         rt_residual=Rt-fit)

######
## plot out
pp1 <- ggplot(data = rt_merged1) +
  geom_boxplot(aes(
    x = factor(scenario),
    y = fit,fill=rt_cat),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  geom_hline(yintercept = 1,linetype = 'dashed',size = .3,color = 'grey') +
  scale_x_discrete(name = 'Scenario',
                   expand = c(0.05, 0.05),
                   labels = c("Sample all\ndetected cases\n(1)","Sample 80%\ndetected cases\n(2)",
                              "Sample 50%\ndetected cases\n(3)","Sample 30%\ndetected cases\n(4)")) +
  scale_y_continuous(limits = c(0,3),name="Ct-based Rt") +
  scale_fill_manual(name = 'Simulation truth',
                    values = colorRampPalette(c("lavender","purple3"))(4),
                    labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  labs(title="A") +
  theme(legend.position = c(0.3,0.9))

pp2 <- ggplot(data = rt_merged2) +
  geom_boxplot(aes(
    x = factor(scenario),
    y = fit,fill=rt_cat),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  geom_hline(yintercept = 1,linetype = 'dashed',size = .3,color = 'grey') +
  scale_x_discrete(name = 'Scenario',
                   expand = c(0.05, 0.05),
                   labels = c("Epidemic plateaued\ndue to detection limit\n(5)",
                              "Epidemic with Rt\nfluctuating around 1\n(6)")) +
  scale_y_continuous(limits = c(0,3),name="Ct-based Rt") +
  scale_fill_manual(name = 'Simulation truth',
                    values = colorRampPalette(c("lavender","purple3"))(4),
                    labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  labs(title="B")+
  theme(legend.position = c(0.3,0.9))

pp3 <- ggplot(data = rt_merged3) +
  geom_boxplot(aes(
    x = factor(scenario),
    y = fit,fill=rt_cat),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  geom_hline(yintercept = 1,linetype = 'dashed',size = .3,color = 'grey') +
  scale_x_discrete(name = 'Scenario',
                   expand = c(0.05, 0.05),
                   labels = c("Biased detection\nall time\n(7)",
                              "Biased detection\nduring testing\n(8)",
                              "Biased and increasing\ndetection during testing\n(9)",
                              "Heavily biased\nduring testing\n(10)")) +
  scale_y_continuous(limits = c(0,3),name="Ct-based Rt") +
  scale_fill_manual(name = 'Simulation truth',
                    values = colorRampPalette(c("lavender","purple3"))(4),
                    labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  labs(title="C") +
  theme(legend.position = c(0.3,0.9))

pp4 <- ggplot(data = rt_merged4) +
  geom_boxplot(aes(
    x = factor(scenario),
    y = fit,fill=rt_cat),
    width = 0.3,
    notchwidth = .1,
    outlier.size = .1,
    lwd=.2) + 
  geom_hline(yintercept = 1,linetype = 'dashed',size = .3,color = 'grey') +
  scale_x_discrete(name = 'Scenario',
                   expand = c(0.05, 0.05),
                   labels = c("Biased and delayed\ndetection all time\n(11)",
                              "Biased and delayed\ndetectionduring testing\n(12)",
                              "Biased and increasing\ndetectionduring testing with delay\n(13)",
                              "Heavily biased with delay\nduring testing\n(14)")) +
  scale_y_continuous(limits = c(0,3),name="Ct-based Rt") +
  scale_fill_manual(name = 'Simulation truth',
                    values = colorRampPalette(c("lavender","purple3"))(4),
                    labels = c('<0.5', '0.5-1.0', '1.0-1.5', '>1.5')) +
  labs(title="D") +
  theme(legend.position = c(0.35,0.9))

ggarrange(pp1, pp2, pp3, pp4, 
          ncol = 2, nrow = 2, 
          common.legend = TRUE, legend = "right") # fig_s2 (8*16)

with(rt_merged1%>%filter(scenario==11),table(rt_cat,useNA = 'always'))
with(rt_merged2%>%filter(scenario==6),table(rt_cat,useNA = 'always'))

################### 
### population Ct over time
### provide reasons for the varied performance
set.seed(1)
ct_group1 <- rbind(ct_list[[1]]%>%mutate(scenario=1),
                   ct_list[[2]]%>%sample_frac(.8)%>%mutate(scenario=2),
                   ct_list[[3]]%>%sample_frac(.5)%>%mutate(scenario=3),
                   ct_list[[4]]%>%sample_frac(.3)%>%mutate(scenario=4)) %>%
  mutate(sampled_time=onset_time+confirmation_delay)

p1 <- ggplot(ct_group1,aes(x=sampled_time,y=ct.value,
                           color=factor(scenario))) +
  geom_smooth(method="gam") +
  scale_y_reverse(name="Ct value") +
  coord_cartesian(ylim=c(30,20)) +
  scale_color_manual(name="",values = colorRampPalette(c("royalblue3","lightblue2"))(4),
                     labels = c("Sample all detected cases","Sample 80% detected cases",
                                "Sample 50% detected cases","Sample 30% detected cases")) +
  scale_x_continuous(name="Days")+
  theme(legend.position = c(.8,.8))+labs(title="A")

ct_group2 <- rbind(ct_list[[5]]%>%mutate(scenario=5),
                   ct_list[[6]]%>%mutate(scenario=6)) %>%
  mutate(sampled_time=onset_time+confirmation_delay)

p2 <- ggplot(ct_group2,aes(x=sampled_time,y=ct.value,
                           color=factor(scenario))) +
  geom_smooth(method="gam") +
  scale_y_reverse(name="Ct value") +
  coord_cartesian(ylim=c(30,20)) +
  scale_color_manual(name="",values = c("moccasin","sandybrown"),
                     labels=c("Epidemic plateaued due to detection limit",
                              "Epidemic with Rt fluctuating around 1")) +
  scale_x_continuous(name="Days")+
  theme(legend.position = c(.7,.8))+labs(title="B")

ct_group3 <- rbind(ct_list[[7]]%>%mutate(scenario=7),
                   ct_list[[8]]%>%mutate(scenario=8),
                   ct_list[[9]]%>%mutate(scenario=9),
                   ct_list[[10]]%>%mutate(scenario=10)) %>%
  mutate(sampled_time=onset_time+confirmation_delay)

p3 <- ggplot(ct_group3,aes(x=sampled_time,y=ct.value,
                           color=factor(scenario))) +
  geom_smooth(method="gam") +
  scale_y_reverse(name="Ct value") +
  coord_cartesian(ylim=c(30,20)) +
  scale_color_manual(name="",values = colorRampPalette(c("seashell","red4"))(4),
                     labels=c("Biased detection during both periods",
                              "Biased detection during testing",
                              "Biased and increasing detection during testing",
                              "Heavily biased detection during testing")) +
  scale_x_continuous(name="Days")+
  theme(legend.position = c(.5,.85))+labs(title="C")

ct_group4 <- rbind(ct_list[[11]]%>%mutate(scenario=11),
                   ct_list[[12]]%>%mutate(scenario=12),
                   ct_list[[13]]%>%mutate(scenario=13),
                   ct_list[[14]]%>%mutate(scenario=14)) %>%
  mutate(sampled_time=onset_time+confirmation_delay)

p4 <- ggplot(ct_group4,aes(x=sampled_time,y=ct.value,
                           color=factor(scenario))) +
  geom_smooth(method="gam") +
  scale_y_reverse(name="Ct value") +
  coord_cartesian(ylim=c(30,20)) +
  scale_color_manual(name="",values = colorRampPalette(c("seashell","red4"))(4),
                     labels=c("Biased and delayed detection during both periods",
                              "Biased and delayed detection during testing",
                              "Biased and increasing detection during testing with delay",
                              "Heavily biased detection during testing with delay")) +
  scale_x_continuous(name="Days")+
  theme(legend.position = c(.5,.85))+labs(title="D")

p_all2 <- grid.arrange(grobs=list(p1,p2,p3,p4),
                       heights = c(1, 1), widths = c(1, 1)) # fig_s4 (7*12)


################### 
### performance by epidemic characteristics
truth_list2 <- rt_list2 <- ct_list2 <- NULL
truth_seq <- c(1,1,2)
rt_seq <- 1:3
for (k in 1:3){
  truth_list2[[k]] <- read.csv(paste0("est/truth/sim_truth",truth_seq[k],".csv"))
  rt_list2[[k]] <- read.csv(paste0("est/rt/rt_obs_scenario",rt_seq[k],".csv"))
  ct_list2[[k]] <- read.csv(paste0("est/detected/scenario",rt_seq[k],".csv")) 
}

vl_list2 <- NULL
vl_list2[[1]] <- read.csv("est/viral_detect/est_detect_vl_scenario1.csv") 
vl_list2[[2]] <- read.csv("est/viral_detect/est_detect_vl_scenario2.csv") 

####### define the epidemic periods (based on simulation truth 1)
func_switch <- function(data){
  for (n in 2:length(data)){
    if (data[n-1]*data[n]>0){
      n <- n+1
    } else {
      return(n-1)
    }
  }
}

func_switch(truth_list2[[1]]$Rt-1)
func_switch(truth_list2[[1]]$Rt[64:200]-1)
func_switch(truth_list2[[1]]$Rt[64+46:200]-1)

plot(truth_list2[[1]]$step,truth_list2[[1]]$Rt,type="l")
abline(h=1,lty=2)
#abline(v=c(63,63+45,63+45+43),lty=2)
abline(v=c(100,130),lty=2,col="red")
abline(v=c(130,160),lty=2,col="blue")

## increasing: 100 to 130; declining: 131 to 160
period_start <- c(100,130)

##
evaluate_mat2 <- matrix(NA,7,3)
auc_all <- matrix(NA,7,100)
for (n in 1:3){ # 3 scenarios
  rt_truth <- truth_list2[[n]]
  start.time <- Sys.time()
  ## to resample population Ct
  if (n %in% 1:2){
    vl_full0 <- vl_list2[[1]]
  } else {
    vl_full0 <- vl_list2[[2]]
  }
  # get daily case count (by sampling time)
  ct0 <- ct_list2[[n]]
  daily_count <- ct0 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
  # also keep the training period same as main
  peak_date <- as.Date("2022-01-01") + daily_count %>% 
    ### updated: restrict to the first epi wave
    filter(sampled_time<110) %>% filter(n==max(n)) %>% pull(sampled_time) 
  ## use the first peak date (if multiple)
  train_period <- c(peak_date[1]-21,peak_date[1]+21) 
  #
  rt0 <- rt_list2[[n]]
  #### resampling start here
  n_run_time <- 100
  auc_tmp <- rep(NA,n_run_time) ## running times
  if (n %in% 1:2){ # increasing/decreasing epidemic
    auc_by_time <- matrix(NA,2,n_run_time)
  }
  for (i in 1:n_run_time){
    set.seed(i)
    ct_update <- sample_ct_only(vl_full0,daily_count)
    ## updated estimates
    rt_tmp <- evaluate_daily_funcs(rt0,ct_update,rt_truth,
                                   seq(train_period[1],train_period[2],1)) 
    # make sure days to compare is the same here (across all scenarios)
    rt_evaluate <- rt_tmp%>%filter(date>=(as.Date("2022-01-01")+100)) ## here is 100 instead of 110 as the Rt already increasing on day 100
    # overall
    auc_tmp[i] <- Rt_auc(rt_evaluate)[2]
    if (n%in%1:2){
      for (k in 1:2){
        rt_evaluate_tmp <- rt_tmp%>%filter(date%in%as.character(as.Date("2022-01-01")+period_start[k]:(period_start[k]+30)))
        auc_by_time[k,i] <- Rt_auc(rt_evaluate_tmp)[2]
      }
    }
  }
  #
  if (n %in% 1:2){
    auc_all[3*(n-1)+1,] <- auc_tmp
    auc_all[3*(n-1)+2:3,] <- auc_by_time
  } else {
    auc_all[7,] <- auc_tmp
  }
} 

for (n in 1:7){
  evaluate_mat2[n,] <- round(c(median(auc_all[n,]),quantile(auc_all[n,],c(.025,.975))),3)
}

#####
#pdf("output/fig_s5_29jan.pdf",height = 9,width = 17)
par(fig=c(0,0.5,0,1),mar=c(5,4,4,2)+.1)
plot(NA,xlim=c(.5,7.5),ylim=c(0,1),axes=F,xlab=NA,ylab="Area under the curve")
axis(1,at=1:7,labels = c(rep(c("Overall","Increasing","Decreasing"),2),"Overall"),col.ticks = "white",col = "white")
axis(2,las=1)
lines(c(.5,7.5),rep(evaluate_mat[1,1],2),lty=2,col="gray")
for (k in 1:7){ 
  lines(rep(k,2),evaluate_mat[k,2:3])
  points(k,evaluate_mat[k,1],pch=16)
}
### add more notation for each scenario
axis(3,at=c(.6,3.4),labels=rep(NA,2),tck=.02)
mtext("Scenario 1\n(stable detection)",side=3,line=1,at=2)
axis(3,at=c(3.6,6.4),labels=rep(NA,2),tck=.02)
mtext("Scenario 5\n(limited detection)",side=3,line=1,at=5)
axis(3,at=c(6.6,7.4),labels=rep(NA,2),tck=.02)
mtext("Scenario 6\n(plateaued wave)",side=3,line=1,at=7)

mtext("A",side=3,font=2,cex=1.5,adj=0,line=2)

### epidemic situation (only over testing)
df_used_list <- NULL
scenario_number <- c(1,5,6)
for (k in 1:3){
  df_used_list[[k]] <- 
    func_synthesize_scenario(list(truth_list2[[k]],rt_list2[[k]],ct_list2[[k]])) %>% 
    filter(sampled_time>=100)
}
# plot out
for (k in 1:3){
  df_tmp <- df_used_list[[k]]
  par(fig=c(0.5,1,0.67-0.33*(k-1),1-0.33*(k-1)),new=T,mar=c(3,2,2,3)+.1)
  # case count panel
  plot(NA,xlim=c(100,200),ylim=c(0,355-150*(k%in%2:3)),axes=F,xlab=NA,ylab=NA)
  ## add shaded area for increasing/declining
  if (k %in% 1:2){
    polygon(rep(c(100,130),each=2),c(0,rep(350-150*(k==2),2),0),col=alpha("blue4",.1),border=F)
    polygon(rep(c(130,160),each=2),c(0,rep(350-150*(k==2),2),0),col=alpha("green4",.1),border=F)
    if (k == 1){
      text(115,350,"Increasing epidemic")
      text(145,350,"Decreasing epidemic")
      ## figure legend
      polygon(rep(c(170,173),each=2),c(340,rep(350,2),340),col="orange3",border=F)
      text(174,345,"Reported cases",adj=0)
      lines(c(170,173),rep(315,2),col="red")
      text(174,315,"Simulation truth",adj=0)
      lines(c(170,173),rep(285,2))
      text(174,285,"Incidence-based Rt",adj=0)
      lines(c(170,173),rep(255,2),col="pink2")
      text(174,255,"Ct-based Rt",adj=0)
      ## panel legend
      mtext("B",side=3,font=2,cex=1.5,at=86,line=0)
    }
  }
  for (n in 1:nrow(df_tmp)){
    polygon(rep(c(df_tmp$sampled_time[n]-.5,df_tmp$sampled_time[n]+.5),each=2),
            c(0,rep(df_tmp$count[n],2),0),col="orange3",border=F)
  }
  axis(1,at=5:10*20)
  mtext("Days of reporting",side=1,line=2)
  axis(2,las=1)
  mtext("Daily cases",side=2,line=2.5)
  mtext(paste0("Scenario ",scenario_number[k]),side=3,font=2,adj=0,line=0)
  # Rt panel
  par(new=T)
  plot(NA,xlim=c(100,200),ylim=c(0,4),axes=F,xlab=NA,ylab=NA)
  abline(h=1,lty=2)
  lines(df_tmp$sampled_time,df_tmp$Rt,col="red")
  lines(df_tmp$sampled_time,df_tmp$rt_est)
  lines(df_tmp$sampled_time,df_tmp$fit,col="pink2")
  axis(4,las=1)
  mtext("Rt",side=4,line=1.5,las=1)
}

#dev.off()

####