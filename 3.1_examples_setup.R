########################

## to set up simulation
## for different pathogens (N=5)

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

pars_example <- read.csv("pathogen_examples.csv") # five real-life examples

# steps 1-2: simulate population-level daily incidence and individual-level infected linelist
simulated_scenes <- NULL
for (k in 1:5){
  pars["incubation"] <- round(pars_example$incubation_mean[k]) # update incu for partab 
  set.seed(k)
  simulated_scenes[[k]] <- func_simulate_infected(population_n,times,
                                                  #log-transformed in the function
                                                  pars,pars_incu = c(pars_example$incubation_mean[k],
                                                                     pars_example$incubation_sd[k]))
  write.csv(simulated_scenes[[k]][[1]],paste0("est/2_examples/truth/sim_truth",k,".csv"),row.names = F)
}


# step 3: simulate viral load trajectories for **symptomatic** infected individuals
set.seed(126)
vl_list <- NULL
for (k in 1:5){ 
  ## only simulate for symptomatic cases
  infect_linelist <- simulated_scenes[[k]][[2]] %>% filter(!is.na(onset_time))
  vl_list[[k]] <- 
    get_indiv_trajectory_path(infect_linelist,
                              # log-transformed in the function
                              pars_example$shedding_mean[k],pars_example$shedding_sd[k],
                              peak_type = "varying", peak_time = pars_example$peak_time[k],
                              vl_mean = pars_example$peak_mean[k],vl_sd = pars_example$peak_sd[k])
  write.csv(vl_list[[k]][[1]],paste0("est/2_examples/viral_traj/est_vl_scenario",k,".csv"),row.names = F)
  # get viral load at detection (export)
  write.csv(vl_list[[k]][[2]],paste0("est/2_examples/viral_detect/est_detect_vl_scenario",k,".csv"),row.names = F)
}


# step 4: simulate detected cases and their sampled viral load
## fixed probability of detection with certain daily variation
## should use the same line for all scenarios to reduce scenario uncertainty
set.seed(1)
detect_prob <- tibble(t=times_extended,
                      prob=rnorm(length(times_extended),.25,.005),
                      ver="fixed")
head(detect_prob)

detected_individuals <- NULL
set.seed(1101)
for (i in 1:5){ # 5 real-life examples
  infect_individuals <- vl_list[[i]][[2]] # with their sampled Ct values
  detect_tmp <- simulate_reporting(infect_individuals,
                                   timevarying_prob = detect_prob,
                                   solve_times = times) 
  #
  detect_cases <- detect_tmp[[1]] %>% arrange(infection_time)
  detected_individuals[[i]] <- detect_cases
  write.csv(detect_cases,paste0("est/2_examples/detected/scenario",i,".csv"),row.names=F)
}

# step 5: get incidence-based Rt using EpiNow2
for (m in 1:5){ 
  set.seed(m)
  # delay distributions
  incubation_period <- 
    bootstrapped_dist_fit(rlnorm(1000,log(pars_example$incubation_mean[m]),
                                 log(pars_example$incubation_sd[m])),
                          dist = "lognormal",max_value = 25)
  # GT for simulated data set ***
  pars["incubation"] <- round(pars_example$incubation_mean[m])
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
  get_inc_Rt(detected_individuals[[m]],n=m,path = "2_examples/rt")
}

####