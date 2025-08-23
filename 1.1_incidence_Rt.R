########################

## to get incidence Rt
## five scenarios

## Vania Lin
## updated 2023-12-18
########################

rm(list=ls())

path <- "C:/Users/vanialam/SPH Dropbox/Yun Lin/self/ct_rt_update_program/"
#path <- "~/SPH Dropbox/Yun Lin/self/ct_rt_update_program/" # on mac
setwd(path)
source(paste0(path,"sim_source_general.R"))
source(paste0(path,"func_generic.R"))

#### CONT'D

options(mc.cores = 4)

detected_individuals <- NULL
for (k in 1:13){
  detected_individuals[[k]] <- read.csv(paste0("est/detected/scenario",k,".csv"))
}

# step 5: get incidence-based Rt using EpiNow2
set.seed(777)
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
##
##
for (k in 1:13){ # 13 scenarios for inc-based Rt (updated on 2023-12-18)
  set.seed(k)
  get_inc_Rt(detected_individuals[[k]],n=k)
}

####