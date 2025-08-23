########################

## to get model training
## Ct-based Rt

## Vania Lin
## updated 2024-07-09
########################

rm(list=ls())

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))
source("ggplot_default_editing.R")

#### get overall summary pattern for the training period
#### use the same training period to train model (10 days before and 20 days after peak)
#### then use the same testing period for evaluation

####
## first: evaluate the shared characteristics of training period
func_switch <- function(data){
  for (n in 2:length(data)){
    if (data[n-1]*data[n]>0){
      n <- n+1
    } else {
      return(n-1)
    }
  }
}

#
###
training_mat <- matrix(NA,14,7) ## ignore results for rows 2:4
seq_used <- c(rep(1,4),2:11)
frac_used <- c(.8,.5,.3)
for (n in 1:14){ # 14 scenarios
  ct0 <- read.csv(paste0("est/detected/scenario",seq_used[n],".csv"))
  ### sampling small proportion
  if (n%in%2:4){
    set.seed(n)
    ct_tmp <- ct0 %>% sample_frac(frac_used[n-1])
  } else {
    ct_tmp <- ct0
  }
  rt0 <- read.csv(paste0("est/rt/rt_obs_scenario",seq_used[n],".csv"))
  data_tmp <- merge_Ct_Rt(rt0,ct_tmp)
  train_results <- select_training_period(data_tmp%>%filter(sampled_time<110)) # limit to first wave
  train_period <- train_results[[1]]
  # peak date
  training_mat[n,1] <- 
    # peak date defined as the date with maximum case count *by confirmation*
    as.character(as.Date("2022-01-01")+ ct0 %>% group_by(confirmed_time) %>% 
                   ### only limit to first wave 
                   summarise(n=n()) %>% filter(confirmed_time<110) %>%
                   filter(n==max(n)) %>% pull(confirmed_time))
  # switch point
  rt1 <- rt0 %>% mutate(diff=mean-1)
  training_mat[n,2] <- as.character(rt1$date[func_switch(rt1$diff)])
  # start and end
  training_mat[n,3:4] <- as.character(train_period)
  training_mat[n,5] <- as.numeric(as.Date(train_period[1])-as.Date(training_mat[n,1]))
  training_mat[n,6] <- as.numeric(as.Date(train_period[2])-as.Date(training_mat[n,1]))
  # r-squared
  training_mat[n,7] <- train_results[[2]][9]
}
# summary characteristics of training periods
summary(as.numeric(training_mat[c(1,5:14),5])) # -21
summary(as.numeric(training_mat[c(1,5:14),6])) # 20
summary(as.numeric(training_mat[c(1,5:14),7])) # R-squared quite high (>0.5; median=0.69,IQR=0.62-0.76)
### maybe fix the training period as three weeks before and after the peak date for simplification

### summary plot of training periods
seq_plot <- c(1,5:14)
par(fig=c(0,0.75,0,1))
plot(NA,xlim=c(20,110),ylim=rev(c(0.5,length(seq_plot)+.5)),axes=F,xlab=NA,ylab=NA)
axis(1,at=2:11*10)
axis(2,las=1,at=1:length(seq_plot),labels = seq_plot)
mtext("Scenario",side=2,line=2)
mtext("Days since outbreak",side=1,line=2)
for (k in 1:length(seq_plot)){
  lines(as.numeric(as.Date(training_mat[seq_plot[k],3:4])-as.Date("2022-01-01")),rep(k,2))
  points(as.numeric(as.Date(training_mat[seq_plot[k],1])-as.Date("2022-01-01")),k,pch=16,col="red") # case count peak
  points(as.numeric(as.Date(training_mat[seq_plot[k],2])-as.Date("2022-01-01")),k,pch=17,col="blue") # Rt turning point
  ### indicate the time period of 3 weeks before and after
  polygon(as.numeric(as.Date(training_mat[seq_plot[k],1])-as.Date("2022-01-01"))+rep(c(-21,21),each=2),
          c(k-.15,rep(k-.25,2),k-.15),col="pink",border=F)
}
mtext("A",cex=1.5,side=3,font=2,adj=0)
## R-squared
par(fig=c(0.7,1,0,1),new=T)
plot(NA,xlim=c(0,1),ylim=rev(c(0.5,length(seq_plot)+.5)),axes=F,xlab=NA,ylab=NA)
axis(1)
mtext("Adjusted R-squared",side=1,line=2)
for (k in 1:length(seq_plot)){
  polygon(rep(c(0,training_mat[seq_plot[k],7]),each=2),c(k-.5,rep(k+.5,2),k-.5),col="gray",border="white")
}
abline(v=.5,lty=2,col="blue")
mtext("B",cex=1.5,side=3,font=2,adj=0) # fig_s1 (6*10)

##########################
##########################
## next: fix the training period and assess model performance using ROC 
## training period fixed at -21 days and 21 days after date of case count peak
###
evaluate_mat <- matrix(NA,14,9) # row 1-3: AUC; row 4-6: rho; row 7-9: MAPE
seq_detect_vl <- c(rep(1,5),2,rep(3,4),4,rep(5,3)) 
for (n in 1:14){ # 14 scenarios
  if (n %in%c(1:5,7:14)){
    rt_truth <- read.csv("est/truth/sim_truth1.csv")
  } else if (n==6){
    rt_truth <- read.csv("est/truth/sim_truth2.csv")
  } 
  start.time <- Sys.time()
  ## to resample population Ct
  vl_full0 <- read.csv(paste0("est/viral_detect/est_detect_vl_scenario",
                              seq_detect_vl[n],".csv")) # detected viral load for all infected individuals
  # get daily case count (by sampling time)
  ct0 <- read.csv(paste0("est/detected/scenario",seq_used[n],".csv")) 
  #
  daily_count <- ct0 %>% group_by(sampled_time) %>% summarise(n=n()) %>% ungroup() # for choosing peak date
  if (n%in%7:14){
    # for resampling Ct based on severity status 
    daily_count_mild <- ct0 %>% filter(sev_status=="mild") %>%
      group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
    daily_count_severe <- ct0 %>% filter(sev_status!="mild") %>%
      group_by(sampled_time) %>% summarise(n=n()) %>% ungroup()
  } 
  # also keep the training period same as main
  peak_date <- as.Date("2022-01-01") + daily_count %>% 
    ### updated: restrict to the first epi wave
    filter(sampled_time<=110) %>% filter(n==max(n)) %>% pull(sampled_time) 
  ## use the first peak date (if multiple)
  train_period <- c(peak_date[1]-21,peak_date[1]+21) ## updated training period on 2023-12-18
  #
  rt0 <- read.csv(paste0("est/rt/rt_obs_scenario",seq_used[n],".csv"))
  #### resampling start here
  n_run_time <- 100
  auc_tmp <- rho_tmp <- mape_tmp <- rep(NA,n_run_time) ## running times
  for (i in 1:n_run_time){
    set.seed(i)
    ### for scenarios 2-4 (only sample fractional cases)
    if (n %in% 2:4){
      # do not restrict the resampling only among detected cases 
      # (otherwise the uncertainty will not be comparable to other scenarios)
      ct_update <- sample_ct_only(vl_full0, daily_count) %>% # first resample detected cases
        sample_frac(frac_used[n-1]) # then fraction a portion out for Ct values
    } else if (n%in%c(1,5:6)){
      ct_update <- sample_ct_only(vl_full0,daily_count)
    } else {
      # need to sample mild and severe cases separately (otherwise the severe proportion will be around 20%)
      ct_update1 <- sample_ct_only(vl_full0%>%filter(sev_status=="mild"),daily_count_mild)
      ct_update2 <- sample_ct_only(vl_full0%>%filter(sev_status!="mild"),daily_count_severe)
      ct_update <- rbind(ct_update1,ct_update2)
    }
    ## updated estimates
    rt_tmp <- evaluate_daily_funcs(rt0,ct_update,rt_truth,
                                   seq(train_period[1],train_period[2],1)) 
    # make sure days to compare is the same here (across all scenarios)
    rt_evaluate <- rt_tmp%>%filter(date>=(as.Date("2022-01-01")+110))
    #### main evaluation matrix
    auc_tmp[i] <- Rt_auc(rt_evaluate)[2]
    rho_tmp[i] <- cor.test(log(rt_evaluate$fit),log(rt_evaluate$Rt),
                           method = "spearman",use="na.or.complete")$estimate
    #mae_tmp[i] <- mean(abs(log(rt_evaluate$fit)-log(rt_evaluate$Rt)))
    #rmse_tmp[i] <- (log(rt_evaluate$fit)-log(rt_evaluate$Rt))^2 %>% mean %>% sqrt
    mape_tmp[i] <- mean(abs((log(rt_evaluate$Rt)-log(rt_evaluate$fit))/log(rt_evaluate$Rt))) ## the denominator should be the true Rt
  }
  ## get point and interval estimates
  evaluate_mat[n,1:3] <- round(c(median(auc_tmp),quantile(auc_tmp,c(.025,.975))),3)
  evaluate_mat[n,4:6] <- round(c(median(rho_tmp),quantile(rho_tmp,c(.025,.975))),3)
  # evaluate_mat[n,7:9] <- round(c(median(mae_tmp),quantile(mae_tmp,c(.025,.975))),3)
  evaluate_mat[n,7:9] <- round(c(median(mape_tmp),quantile(mape_tmp,c(.025,.975))),3)
  print(paste0("Completed scenario: ",n,", time = ",Sys.time()-start.time,"."))
} 

(evaluate_mat <- round(evaluate_mat,2))
# save the evaluation matrix
# write.csv(evaluate_mat,"output/evaluate_mat_full_8aug.csv",row.names = F)

#############################################
### plot Fig 2 (performance by detection mode)
## 
### panel A: results
par(fig=c(0,1,0.5,1),mar=c(5,10,4,2)+.1)
plot(NA,xlim=c(.5,14.5),ylim=c(0,1),axes=F,xlab=NA,ylab="AUC")
axis(1,at=1:14,labels = 1:14)
mtext("Scenario",side=1,line=2)
axis(2,las=1)
lines(c(.5,14.5),rep(evaluate_mat[1,1],2),lty=2,col="gray")
for (k in 1:14){ 
  lines(rep(k,2),evaluate_mat[k,2:3])
  points(k,evaluate_mat[k,1],pch=16)
}
### add more notation for each scenario
axis(3,at=c(.6,4.4),labels=rep(NA,2),tck=.02,line=.7)
mtext("Varying sampling fractions",side=3,line=1,at=2.5)
axis(3,at=c(4.6,6.4),labels=rep(NA,2),tck=.02,line=.7)
mtext("Plateaued wave",side=3,line=1,at=5.5)
axis(3,at=c(6.6,14.4),labels=rep(NA,2),tck=.02,line=.7)
mtext("Biased detection towards\n severe cases",side=3,line=.7,at=10.7)
axis(3,at=c(6.6,10.4),labels=rep(NA,2),tck=.02,line=-1.5)
mtext("Similar delay\nfor mild and severe cases",side=3,line=-1.5,at=8.5)
axis(3,at=c(10.6,14.4),labels=rep(NA,2),tck=.02,line=-1.5)
mtext("Various delays\nfor mild and severe cases",side=3,line=-1.5,at=12.5)
mtext("A",side=3,font=2,cex=1.5,at=-.35,line=2)

### panel B: scenario summary
col_detect <- colorRampPalette(c("steelblue4","slategray1"))(4)
col_biased <- colorRampPalette(c("mistyrose1","red4"))(3)
label_used <- c("Sampling\nfraction","Epidemic\ncharacteristics","Limited\ndetection",
                "Biased detection \n(similar delay)","Biased detection \n(various delays)")
wave_used <- c("Both\nperiods",rep(c("Training","Testing"),4))
par(fig=c(0,1,0.05,0.55),mar=c(5,10,4,2)+.1,new=T)
plot(NA,xlim=c(.5,14.5),ylim=rev(c(0.5,9.5)),axes=F,xlab=NA,ylab=NA)
axis(1,line=-1,at=1:14,col="white")
mtext("Scenario",line=1,side=1)
axis(2,at=1:9,line=-2,labels = wave_used,col = "white",las=1,cex.axis=.9)
axis(2,at=c(1,2.5,4.5,6.5,8.5),labels = label_used,tck=.03,line=1.6,col = "white",las=1,cex.axis=.9)
axis(2,at=c(1,2,4,6,8,10)-.5,labels = rep(NA,6),tck=.03,line=2)
###
## varying sampling fraction
col_detect_used <- c(col_detect,rep(col_detect[1],10))
for (k in 1:14){
  polygon(rep(c(k-.5,k+.5),each=2),c(.5,rep(1.5,2),.5),col=col_detect_used[k],border=F)
}
# limited detection/flat waves
polygon(rep(c(4.5,5.5),each=2),c(4.5,rep(5.5,2),4.5),col="snow3",border=F)
polygon(rep(c(5.5,6.5),each=2),c(2.5,rep(3.5,2),2.5),col="snow3",border=F)
# biased detection
polygon(rep(c(6.5,7.5),each=2),c(5.5,rep(6.5,2),5.5),col=col_biased[1],border=F)
polygon(rep(c(10.5,11.5),each=2),c(7.5,rep(8.5,2),7.5),col=col_biased[1],border=F)
col_biased_used <- c(col_biased[1],col_biased)
for (k in 1:2){
  for (i in 1:4){
    polygon(rep(c(6.5,7.5)+i-1+4*(k-1),each=2),c(6.5,rep(7.5,2),6.5)+2*(k-1),col=col_biased_used[i],border=F)
    #polygon(rep(c(6.5,7.5),each=2)+2*(k-1),c(6.5,rep(7.5,2),6.5)+i-1+4*(k-1),col=col_biased_used[i],border=F)
  }
}
## box border
for (k in 1:14){
  for (n in 1:9){
    polygon(rep(c(k-.5,k+.5),each=2),c(n-.5,rep(n+.5,2),n-.5))
  }
}
mtext("B",side=3,font=2,at=-.35,cex=1.5,line=1.5) 

####
## legend for color codes
par(fig=c(0,1,0,.15),mar=c(2,2,2,1)+.1,new=T)
plot(NA,xlim=c(0,20),ylim=c(0,1),axes=F,xlab=NA,ylab=NA)
polygon(rep(c(2.5,3.5),each=2),c(0.4,rep(0.6,2),0.4),col="snow3",border=F)
text(3.6,.5,"Plateaued\nwave",adj=0)
for (k in 1:4){
  polygon(rep(c(5.5,6.5),each=2)+k-1,c(0.4,rep(0.6,2),0.4),col=col_detect[k],border=F)
}
text(9.6,.5,"Sampling fraction\n(100% to 30%)",adj=0)
#
for (k in 1:3){
  polygon(rep(c(12.5,13.5),each=2)+k-1,c(.4,rep(.6,2),.4),col=col_biased[k],border=F)
}
text(15.6,.5,"Biased detection towards severe\ncases (moderate to heavy)",adj=0)

# fig 2 (9*12)

#############################################
##### plot Fig S3 (more evaluation matrix)
##
# plot the additional evaluation matrix
par(mfrow=c(1,3),mar=c(3,3,2,1)+.1,cex=1)
# for AUC
plot(NA,xlim=c(0,1.5),ylim=rev(c(0.5,14.5)),axes=F,xlab=NA,ylab=NA)
axis(2,at=1:14,las=1,col ="white")
axis(1,at=0:5*0.2)
mtext("Scenario",side=2,line=1.9)
mtext("AUC",side=1,line=2,at=.5)
for (k in 1:14){
  lines(c(evaluate_mat[k,2],evaluate_mat[k,3]),rep(k,2),lwd=2)
  points(evaluate_mat[k,1],k,pch=19)
  text(1,k,paste0(format(evaluate_mat[k,1],nsmall=2)," (",
                  format(evaluate_mat[k,2],nsmall=2),",",
                  format(evaluate_mat[k,3],nsmall=2),")"),adj=0)
}
mtext("A",side=3,adj=0,font=2,cex=1.2)

#
# for rho
plot(NA,xlim=c(0,1.5),ylim=rev(c(0.5,14.5)),axes=F,xlab=NA,ylab=NA)
axis(2,at=1:14,las=1,col="white")
axis(1,at=0:5*0.2)
mtext("Spearman's rho",side=1,line=2,at=.5)
for (k in 1:14){
  lines(c(evaluate_mat[k,5],evaluate_mat[k,6]),rep(k,2),lwd=2)
  points(evaluate_mat[k,4],k,pch=19)
  text(1,k,paste0(format(evaluate_mat[k,4],nsmall=2)," (",
                  format(evaluate_mat[k,5],nsmall=2),",",
                  format(evaluate_mat[k,6],nsmall=2),")"),adj=0)
}
mtext("B",side=3,adj=0,font=2,cex=1.2)

#
# for MAPE
plot(NA,xlim=rev(c(-5.5,log(16,base = 2))),ylim=rev(c(0.5,14.5)),axes=F,xlab=NA,ylab=NA)
axis(2,at=1:14,las=1,col="white")
axis(1,at=log(c(0.2,0.5,1,2,4,16),base=2),labels = c(0.2,0.5,1,2,4,16))
mtext("MAPE",side=1,line=2,at=1)
for (k in 1:14){
  lines(log(c(evaluate_mat[k,8],evaluate_mat[k,9]),base=2),rep(k-0.1,2),lwd=2)
  points(log(evaluate_mat[k,7],base=2),k-0.1,pch=19)
  text(-2.5,k,paste0(format(evaluate_mat[k,7],nsmall=2)," (",
                  format(evaluate_mat[k,8],nsmall=2),",",
                  format(evaluate_mat[k,9],nsmall=2),")"),adj=0)
}
# add legend (symbols and texts) indicating blue is for MAE and red is for MAPE
mtext("C",side=3,adj=0,font=2,cex=1.2)

### fig_s3_new_matrix (7*15)

####