########################

## to set up scenarios
## for various detection modes

## Vania Lin
## updated 2024-01-30

########################

rm(list=ls())

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))


#############################################################################################
#############################################################################################
# steps 1-2: simulate population-level daily incidence and individual-level infected linelist
simulated_scenes <- NULL
### simulation truth (truth1)
set.seed(1)
simulated_scenes[[1]] <- func_simulate_infected(population_n,times,
                                                par=pars,pars_incu = c(5.2,3.9))
write.csv(simulated_scenes[[1]][[1]],"est/truth/sim_truth1.csv",row.names = F)
#
#### flat second wave (truth2)
pars_alt0 <- read.csv("partab_seir_switch_model_alt.csv")
pars_alt <- pars_alt0$values
names(pars_alt) <- pars_alt0$names

set.seed(2)
simulated_scenes[[2]] <- func_simulate_infected(population_n,times,
                                                par = pars_alt,pars_incu = c(5.2,3.9))
write.csv(simulated_scenes[[2]][[1]],"est/truth/sim_truth2.csv",row.names = F)


# step 3: simulate viral load trajectories for **symptomatic** infected individuals
set.seed(1114)
vl_list <- NULL
for (k in 1:2){ 
  ## only simulate for symptomatic cases
  infect_linelist <- simulated_scenes[[k]][[2]] %>% filter(!is.na(onset_time))
  vl_list[[k]] <- get_indiv_trajectory(infect_linelist)
  # get individual viral trajectory (export)
  write.csv(vl_list[[k]][[1]],paste0("est/viral_traj/est_vl_scenario",k,".csv"),row.names = F)
  # get viral load at detection (export)
  write.csv(vl_list[[k]][[2]],paste0("est/viral_detect/est_detect_vl_scenario",k,".csv"),row.names = F)
}

### for sev
set.seed(1214)
infect_linelist2 <- func_assign_sev(simulated_scenes[[1]][[2]] %>% filter(!is.na(onset_time))) ## linelist with sev_status assigned
##### assuming no delay detection for severe cases, the only difference is the shedding duration
set.seed(1214) 
vl_list[[3]] <- get_indiv_trajectory(infect_linelist2,incl_sev = T)
# get individual viral trajectory (export)
write.csv(vl_list[[3]][[1]],"est/viral_traj/est_vl_scenario3.csv",row.names = F) 
# get viral load at detection (export)
write.csv(vl_list[[3]][[2]],"est/viral_detect/est_detect_vl_scenario3.csv",row.names = F)

### for different delayed detections 
### 1.delayed in both training and testing
set.seed(4)
infect_linelist3 <- infect_linelist2 %>% 
  ### update confirmation delay for severe cases
  mutate(confirmation_delay=ifelse(sev_status=="mild",confirmation_delay,
                                   floor(rgamma(n(), shape = 3,rate = 0.43))))
summary(infect_linelist2$confirmation_delay[infect_linelist2$sev_status=="mild"])
summary(infect_linelist2$confirmation_delay[infect_linelist2$sev_status!="mild"])
summary(infect_linelist3$confirmation_delay[infect_linelist3$sev_status=="mild"])
summary(infect_linelist3$confirmation_delay[infect_linelist3$sev_status!="mild"]) # checked
##### assuming delay detection and prolonged shedding duration for severe cases
set.seed(1214)
vl_list[[4]] <- get_indiv_trajectory(infect_linelist3,incl_sev = T)
# get individual viral trajectory (export)
write.csv(vl_list[[4]][[1]],"est/viral_traj/est_vl_scenario4.csv",row.names = F)
# get viral load at detection (export)
write.csv(vl_list[[4]][[2]],"est/viral_detect/est_detect_vl_scenario4.csv",row.names = F)

### 2.delayed in testing
summary(infect_linelist2$onset_time) # checked
set.seed(5)
infect_linelist4 <- infect_linelist2 %>% 
  ### update confirmation delay for severe cases
  ## individual delay (for severe cases) may not be the same for vl_list[[4]]
  ## since the number for gamma distribution were different in the two scenarios
  mutate(confirmation_delay=ifelse(onset_time<110|sev_status=="mild",confirmation_delay,
                                   floor(rgamma(n(), shape = 3,rate = 0.43))))
summary(infect_linelist3$confirmation_delay[infect_linelist3$sev_status!="mild"&infect_linelist3$onset_time<110])
summary(infect_linelist4$confirmation_delay[infect_linelist4$sev_status!="mild"&infect_linelist4$onset_time<110]) 
summary(infect_linelist3$confirmation_delay[infect_linelist3$sev_status!="mild"&infect_linelist3$onset_time>=110])
summary(infect_linelist4$confirmation_delay[infect_linelist4$sev_status!="mild"&infect_linelist4$onset_time>=110]) # checked
##### assuming delay detection and prolonged shedding duration for severe cases
set.seed(1214)
vl_list[[5]] <- get_indiv_trajectory(infect_linelist4,incl_sev = T)
# get individual viral trajectory (export)
write.csv(vl_list[[5]][[1]],"est/viral_traj/est_vl_scenario5.csv",row.names = F)
# get viral load at detection (export)
write.csv(vl_list[[5]][[2]],"est/viral_detect/est_detect_vl_scenario5.csv",row.names = F)

#### by making the seed similar across the three vl_list
#### can make sure individual trajectories the same across these three scenarios
#### while the detected viral load changes as the detection delay was updated for severe cases

#########################################
#########################################

# step 4: simulate detected cases and their sampled viral load
## fixed probability of detection with certain daily variation
## should use the same line for all scenarios to reduce scenario uncertainty
set.seed(1)
detect_prob <- tibble(t=times_extended,
                      prob=rnorm(length(times_extended),.25,.005),
                      ver="fixed")
head(detect_prob)

####### detected scenarios (N=11) #######
### scenario 1: two waves, regular detection
### scenario 2: two waves, flat detection during testing period
### scenario 3: flat testing wave, regular detection
### scenario 4-6: two waves, biased detection (**severity**) 
###               over both periods (4), or over testing period (5), or increasing during peak (6)
### scenario 7: biased detection -> only severe cases during testing 
### scenarios 8-11: delayed detection for severe cases (flat detection prob during both periods)
###                 over both periods (8), during testing (9), or increasing during peak (10) 
###                 or during testing when only severe would be detected (11)
detected_individuals <- NULL
##
### for scenarios 1-2:
set.seed(1115)
infect_individuals <- vl_list[[1]][[2]] # with their sampled Ct values
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times) 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[1]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario1.csv",row.names=F)  
#
#
detect_count1 <- detect_cases %>% group_by(confirmed_time) %>% summarise(n=n()) %>% ungroup()
with(detect_count1,plot(confirmed_time,n,type="b"))
abline(v=110,lty=2)
abline(h=c(50,100,150),col="blue4")
### for days>110, if case count is over 100, then (randomly) change to only able to detect 100-120 cases
sample_frac_days <- detect_count1 %>% filter(confirmed_time>110,n>100) %>% pull(confirmed_time)
set.seed(11)
sample_n <- sample(100:120,length(sample_frac_days),T)
detect_count <- detect_count1 %>% filter(confirmed_time%in%sample_frac_days) %>% 
  mutate(sample_n=sample_n) %>% rowwise() %>% mutate(actual_n=min(c(sample_n,n))) %>% ungroup()
#
detect_cases_orig <- detect_cases %>% filter(!confirmed_time%in%sample_frac_days)
detect_cases_change1 <- detect_cases %>% filter(confirmed_time%in%sample_frac_days) 
detect_cases_merged <-  detect_cases_orig ## baseline for merging
for (k in 1:length(sample_frac_days)){
  detect_cases_tmp <- detect_cases_change1 %>% filter(confirmed_time==sample_frac_days[k])
  detect_cases_get <- detect_cases_tmp[sample(1:nrow(detect_cases_tmp),detect_count$actual_n[k],F),]
  detect_cases_merged <- bind_rows(detect_cases_merged,detect_cases_get)
}
detect_cases2 <- detect_cases_merged %>% arrange(infection_time)  
##
## check flat testing
count_tmp <- detect_cases2 %>% group_by(confirmed_time) %>% summarise(n=n()) %>% ungroup()
with(count_tmp,plot(confirmed_time,n,type="b"))
#
detected_individuals[[2]] <- detect_cases2  
write.csv(detect_cases2,"est/detected/scenario2.csv",row.names=F) 


### for scenario 3:
set.seed(15)
infect_individuals <- vl_list[[2]][[2]] # with their sampled Ct values
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times) 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[3]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario3.csv",row.names=F)  


### for scenarios 4-7:
set.seed(17)
infect_individuals <- vl_list[[3]][[2]] # with severity status
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "both") 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[4]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario4.csv",row.names=F)  
#
set.seed(71)
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "testing") 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[5]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario5.csv",row.names=F)  

set.seed(6)
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "peak") 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[6]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario6.csv",row.names=F)  

set.seed(7)
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "heavy") 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[7]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario7.csv",row.names=F)  

### for scenario 8-11
## scenario 8
set.seed(8)
infect_individuals <- vl_list[[4]][[2]] # with severity status and *delayed detection* for severe cases during **both** periods
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "both") ## extension of scenario 4
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[8]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario8.csv",row.names=F)  
#
set.seed(9)
infect_individuals <- vl_list[[5]][[2]] # with severity status and *delayed detection* for severe cases during **testing** period
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "testing") ## extension of scenario 5
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[9]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario9.csv",row.names=F)  
#
set.seed(10)
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "peak") # extension of scenario 6
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time) 
detected_individuals[[10]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario10.csv",row.names=F)  

set.seed(11)
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob,
                                 solve_times = times,incl_sev = "heavy") ## extension of scenario 7
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[11]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario11.csv",row.names=F) 


############ add two more scenarios of only detecting severe cases during the testing periods
##### date version: 2023-12-20 (not added to main scenarios yet)
detect_prob2 <- detect_prob
logistic_func <- function(t,start_prob,end_prob, growth_rate, switch_point=100){
  start_prob + (end_prob-start_prob)/(1 + exp(-growth_rate*(t-switch_point)))
}

kk <- logistic_func(101:250,start_prob=.05,end_prob=.5,growth_rate=0.07,
                    # switch point before start of second wave
                    switch_point=100)
for (n in 101:250){
  detect_prob2$prob[detect_prob2$t==n] <- kk[n-100]
}
plot(detect_prob2$t,detect_prob2$prob,type="l")
abline(v=100,lty=2)

## without delay
## scenario added 1
set.seed(12)
infect_individuals <- vl_list[[3]][[2]] 
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob2,
                                 solve_times = times,incl_sev = "heavy") 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[12]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario12.csv",row.names=F)  
#
## with delay
## scenario added 2
set.seed(13)
infect_individuals <- vl_list[[5]][[2]] 
detect_tmp <- simulate_reporting(infect_individuals,
                                 timevarying_prob = detect_prob2,
                                 solve_times = times,incl_sev = "heavy") 
detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
detected_individuals[[13]] <- detect_cases
write.csv(detect_cases,"est/detected/scenario13.csv",row.names=F)  
#

############### check severity scenarios (4-7)
### check severity proportion
df_tmp1 <- detected_individuals[[4]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop1=n_severe/n) %>% ungroup()
df_tmp2 <- detected_individuals[[5]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop2=n_severe/n) %>% ungroup()
df_tmp3 <- detected_individuals[[6]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop3=n_severe/n) %>% ungroup()
df_tmp4 <- detected_individuals[[7]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop4=n_severe/n) %>% ungroup()
df_tmp <- inner_join(df_tmp1%>%dplyr::select(confirmed_time,prop1),
                     df_tmp2%>%dplyr::select(confirmed_time,prop2)) %>%
  inner_join(.,df_tmp3%>%dplyr::select(confirmed_time,prop3)) %>%
  inner_join(.,df_tmp4%>%dplyr::select(confirmed_time,prop4))
with(df_tmp,plot(confirmed_time,prop1,type="l"))
lines(df_tmp$confirmed_time,df_tmp$prop2,col="red") 
lines(df_tmp$confirmed_time,df_tmp$prop3,col="blue")
lines(df_tmp$confirmed_time,df_tmp$prop4,col="green4") # checked

### check scenario 8-11
### check severity proportion
df_tmp1 <- detected_individuals[[8]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop1=n_severe/n) %>% ungroup()
df_tmp2 <- detected_individuals[[9]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop2=n_severe/n) %>% ungroup()
df_tmp3 <- detected_individuals[[10]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop3=n_severe/n) %>% ungroup()
df_tmp4 <- detected_individuals[[11]] %>% group_by(confirmed_time) %>% 
  summarise(n=n(),n_severe=sum(sev_status!="mild"),prop4=n_severe/n) %>% ungroup()
df_tmp <- inner_join(df_tmp1%>%dplyr::select(confirmed_time,prop1),
                     df_tmp2%>%dplyr::select(confirmed_time,prop2)) %>%
  inner_join(.,df_tmp3%>%dplyr::select(confirmed_time,prop3)) %>%
  inner_join(.,df_tmp4%>%dplyr::select(confirmed_time,prop4)) 
with(df_tmp,plot(confirmed_time,prop1,type="l"))
lines(df_tmp$confirmed_time,df_tmp$prop2,col="red") 
lines(df_tmp$confirmed_time,df_tmp$prop3,col="blue") 
lines(df_tmp$confirmed_time,df_tmp$prop4,col="green4") # checked


### check trajectories for mild/severe cases as well
# vl_traj_check <-vl_list[[3]][[1]]
# sev_strings <- c("mild","severe","critical","fatal")
# par(mfrow=c(2,2))
# for (k in 1:4){
#   vl_traj_tmp <- vl_traj_check %>% filter(sev_status==sev_strings[k])
#   plot(NA,xlim=c(0,40),ylim=rev(c(10,40)),axes=F,xlab=NA,ylab=NA)
#   axis(1)
#   axis(2,las=1)
#   abline(v=20,col="red",lty=2)
#   unique_id <- unique(vl_traj_tmp$i)
#   sampled_id <- sample(unique_id,200,F)
#   for (n in 1:200){
#     vl_indiv_tmp <- vl_traj_tmp %>% filter(i==sampled_id[n])
#     lines(vl_indiv_tmp$test.to.infect,vl_indiv_tmp$ct.value,col=alpha("gray",.4))
#   }
#   mtext(sev_strings[k],side=3,font=2,adj=0)
# } # checked

#save.image("data/scenarios.RData")

####