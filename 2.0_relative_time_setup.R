########################

## to set up simulation
## for different relations 
## viral peak vs onset

## Vania Lin
## updated 2024-01-30
########################

rm(list=ls())

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))
source(paste0(path,"func_pathogens.R"))

# steps 1-2: simulate population-level daily incidence and individual-level infected linelist
set.seed(1)
simulated_scene <- func_simulate_infected(population_n,times,pars,pars_incu = c(5.2,3.9))
write.csv(simulated_scene[[1]],"est/1_relative_time/truth/sim_truth.csv",row.names = F)

###########
# step 3: simulate viral load trajectories for **symptomatic** infected individuals (N=4)
infect_linelist <- simulated_scene[[2]] %>% filter(!is.na(onset_time))
##
set.seed(2)
vl_list <- NULL
## viral peak fixed at day 5
vl_list[[1]] <- 
  get_indiv_trajectory_path(infect_linelist,duration_mean=22.6,duration_sd=1.3,
                            vl_mean=22.3,vl_sd=4.2)
write.csv(vl_list[[1]][[1]],"est/1_relative_time/viral_traj/est_vl_scenario1.csv",row.names = F)
# get viral load at detection (export)
write.csv(vl_list[[1]][[2]],"est/1_relative_time/viral_detect/est_detect_vl_scenario1.csv",row.names = F)

### before or after
values <- c(-2,0,2) ## make the time after onset scenario larger difference**
for (k in 1:3){
  set.seed(17+k)
  vl_list[[k+1]] <- 
    get_indiv_trajectory_path(infect_linelist,duration_mean=22.6,duration_sd=1.3,
                              peak_type="varying",peak_time = values[k],vl_mean=22.3,vl_sd=4.2)
  write.csv(vl_list[[k+1]][[1]],paste0("est/1_relative_time/viral_traj/est_vl_scenario",k+1,".csv"),row.names = F)
  # get viral load at detection (export)
  write.csv(vl_list[[k+1]][[2]],paste0("est/1_relative_time/viral_detect/est_detect_vl_scenario",k+1,".csv"),row.names = F)
}

######### check the trajectory for the four scenarios
# par(mfrow=c(4,1))
# for (k in 1:4){
#   vl_traj_tmp <- vl_list[[k]][[1]]
#   plot(NA,xlim=c(0,40),ylim=rev(c(10,40)),axes=F,xlab=NA,ylab=NA)
#   axis(1)
#   axis(2,las=1)
#   unique_id <- unique(vl_traj_tmp$i)
#   sampled_id <- sample(unique_id,200,F)
#   for (n in 1:200){
#     vl_indiv_tmp <- vl_traj_tmp %>% filter(i==sampled_id[n])
#     lines(vl_indiv_tmp$test.to.infect,vl_indiv_tmp$ct.value,col=alpha("gray",.4))
#   }
#   abline(v=5,col="red",lty=2) # onset time
# } # checked

options(mc.cores = 4)

# step 4: simulate detected cases and their sampled viral load
## fixed probability of detection with certain daily variation
## should use the same line for all scenarios to reduce scenario uncertainty
set.seed(1)
detect_prob <- tibble(t=times_extended,
                      prob=rnorm(length(times_extended),.25,.005),
                      ver="fixed")
head(detect_prob)

#### only have one linelist of detected individuals
### only the detected viral loads under each scenario could be different
detected_individuals <- simulate_reporting(infect_linelist,
                                           timevarying_prob = detect_prob,
                                           solve_times = times) 
#
detect_cases <- detected_individuals[[1]] %>% arrange(infection_time)
write.csv(detect_cases,"est/1_relative_time/detected/scenario0.csv",row.names=F)

# step 5: get incidence-based Rt using EpiNow2
set.seed(1)
# delay distributions
incubation_period <- 
  bootstrapped_dist_fit(rlnorm(1000,log(5.2),log(3.9)),dist = "lognormal",max_value = 25)
# GT for simulated data set ***
## Use generation time from SEIR model. 
## Convolution of two exponential distributions 
## mean = to sum of means from both 
## variance = to sum of variances of both
SI_mean <- pars["infectious"] + pars["incubation"]
SI_sd <- sqrt(pars["infectious"]^2 + pars["incubation"]^2)
generation_time <- list(mean=SI_mean,mean_sd=3,sd=SI_sd,sd_sd=3,max=100)
reporting_delay <- 
  bootstrapped_dist_fit(rgamma(1000,shape=1.83,rate=0.43),max_value = 25)
#
get_inc_Rt(detect_cases,n=0,path="1_relative_time/rt")


save.image(file = "data/relative_time.RData")

####