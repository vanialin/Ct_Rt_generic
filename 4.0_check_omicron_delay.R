#### find out the longer detection delay during wave 5 using Hong Kong data

dat_hk0 <- read.csv("~/SPH Dropbox/Yun Lin/self/data_copy/data_HA/cleaned/data_HA_20230212.csv")

dat_hk1 <- dat_hk0 %>% 
  mutate(date.confirm=as.Date(date.confirm),
         onset.to.confirm=as.numeric(as.Date(date.confirm)-as.Date(date.onset)),
         period=case_when(date.confirm%in%seq(as.Date("2021-12-31"),as.Date("2022-05-22"),1)~"w5",
                          date.confirm%in%seq(as.Date("2022-05-23"),as.Date("2022-09-30"),1)~"w6",
                          date.confirm>=as.Date("2022-10-01")~"w7",
                          date.confirm<as.Date("2020-07-01")~"w1-2",
                          date.confirm>=as.Date("2020-07-01")&date.confirm<as.Date("2021-12-31")~"w3-4"))

dat_hk2 <- dat_hk1 %>% filter(!classification%in%c("Imported"))
aggregate(onset.to.confirm~period, data=dat_hk2, summary)

aaa <- dat_hk2 %>% filter(period=="w1-2",!is.na(date.onset),onset.to.confirm>0) %>% pull(onset.to.confirm)
fitdistr(aaa, "gamma") # likely the original detection delay is fitted using data from waves 1-2

# bbb <- dat_hk2 %>% filter(period%in%c("w6","w7"),!is.na(date.onset),onset.to.confirm>0) %>% pull(onset.to.confirm)
# fitdistr(bbb, "gamma") 

bbb <- dat_hk2 %>% filter(date.confirm<=as.Date("2022-02-06"),!is.na(date.onset),onset.to.confirm>0) %>% pull(onset.to.confirm)
fitdistr(bbb, "gamma") ## maybe use this