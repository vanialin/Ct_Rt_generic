####
require(EnvStats)
require(tidyverse)
require(MASS)

rm(list=ls())

###### to align shedding durations from infection
###### also to convert the incubation period into lognormal

approx_sd <- function(x1, x2){
  (x2-x1) / (qnorm(0.975) - qnorm(0.025)) 
}

### numbers in the excel table were all for lognormal distr (for incubation period)
set.seed(38)
incu0 <- rlnormTrunc(10000,log(5.2),log(3.9),min=0,max=25)
summary(incu0)
#
set.seed(1)
approx_sd(15.5,18.6)
dur_moderate0 <- rnormTrunc(10000,17,0.79,min=0,max=25)
dur_moderate1 <- dur_moderate0 + incu0
summary(dur_moderate1)
a <- fitdistr(dur_moderate1,"normal")
b <- fitdistr(dur_moderate1,"lognormal")
AIC(a,b) # lognormal better
(mean_mod_dur <- exp(fitdistr(dur_moderate1,"lognormal")$estimate[1]))
(sd_mod_dur <- exp(fitdistr(dur_moderate1,"lognormal")$estimate[2]))

### 1. alpha
## incubation mean of 5.0, SD 2.3
set.seed(1)
incu_alpha <- rnormTrunc(10000,5.0,2.3,min=0,max=25)
exp(fitdistr(incu_alpha,"lognormal")$estimate[1])
exp(fitdistr(incu_alpha,"lognormal")$estimate[2])
## use the relative difference of pre-alpha and alpha 
## referenced from Sunagawa to get the shedding duration
set.seed(2)
dur_ancestral <- rlnormTrunc(10000,log(22.6),log(1.3),min=0)
dur_alpha <- dur_ancestral - rnorm(10000,3,.25) ## ~3 days shorter (SD=0.25 to add noise)
exp(fitdistr(dur_alpha,"lognormal")$estimate[1])
exp(fitdistr(dur_alpha,"lognormal")$estimate[2])

### 2. delta
## incubation mean of 4.3, SD 2.4
set.seed(3)
incu_delta <- rnormTrunc(10000,4.3,2.4,min=0,max=25)
exp(fitdistr(incu_delta,"lognormal")$estimate[1])
exp(fitdistr(incu_delta,"lognormal")$estimate[2])
##
set.seed(4)
dur_delta <- dur_ancestral + rnormTrunc(10000,6,.25) ## ~ 6 days longer
exp(fitdistr(dur_delta,"lognormal")$estimate[1])
exp(fitdistr(dur_delta,"lognormal")$estimate[2])

### 3. omicron
## incubation mean of 3.2, SD 2.2
set.seed(5)
incu_omicron <- rnormTrunc(10000,3.2,2.2,min=0,max=25)
exp(fitdistr(incu_omicron,"lognormal")$estimate[1])
exp(fitdistr(incu_omicron,"lognormal")$estimate[2])
##
set.seed(6)
## use the incubation period to add on the shedding from onset for ancestral
dur_post_onset <- rnormTrunc(10000,17,0.79,min=0,max=25)
dur_omicron <- incu_omicron + dur_post_onset
exp(fitdistr(dur_omicron,"lognormal")$estimate[1])
exp(fitdistr(dur_omicron,"lognormal")$estimate[2])

### 4. SARS
approx_sd(3.6,4.4)
set.seed(7)
incu_sars <- rnormTrunc(10000,4,0.2,min=0,max=25)
fitdistr(incu_sars,"lognormal")
exp(fitdistr(incu_sars,"lognormal")$estimate[1])
exp(fitdistr(incu_sars,"lognormal")$estimate[2])
set.seed(8)
#incu_sars1 <- rlnormTrunc(10000,1.38,0.05,min=0,max=25)
dur_sars0 <- rnormTrunc(10000,20,1,min=0,max=25)
dur_sars1 <- dur_sars0 + incu_sars
exp(fitdistr(dur_sars1,"lognormal")$estimate[1])
exp(fitdistr(dur_sars1,"lognormal")$estimate[2])

### 5. fluA (Lessler 2008 NEJM)
approx_sd(1,1.8) # caveat: it is not normally distributed
set.seed(9)
incu_flu <- rnormTrunc(10000,1.4,0.2,min=0,max=25)
fitdistr(incu_flu,"lognormal")
exp(fitdistr(incu_flu,"lognormal")$estimate[1])
exp(fitdistr(incu_flu,"lognormal")$estimate[2])
#
approx_sd(4.31,5.29) # sd of shedding duration (wild-type flu; did not differ by subtypes)

######################
##### viral peak #####
######################
set.seed(1)
vl_moderate0 <- rnormTrunc(10000,22.3,4.2,min=10,max=35)
log_vl <- 10^(14.159-vl_moderate0*0.297)
log_vl_higher <- log_vl*10^0.7 # Alpha
vl_higher <- (14.159-log(log_vl_higher,base=10))/0.297
fitdistr(vl_higher,"normal")$estimate[1]
fitdistr(vl_higher,"normal")$estimate[2]

log_vl_higher2 <- log_vl*10^0.6 # Delta
vl_higher2 <- (14.159-log(log_vl_higher2,base=10))/0.297
fitdistr(vl_higher2,"normal")$estimate[1]
fitdistr(vl_higher2,"normal")$estimate[2]

log_vl_sars <- log_vl*1.9
vl_sars <- (14.159-log(log_vl_sars,base=10))/0.297
fitdistr(vl_sars,"normal")$estimate[1]
fitdistr(vl_sars,"normal")$estimate[2]

log_vl_flu <- log_vl*1.6
vl_flu <- (14.159-log(log_vl_flu,base=10))/0.297
fitdistr(vl_flu,"normal")$estimate[1]
fitdistr(vl_flu,"normal")$estimate[2]

####