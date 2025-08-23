### CONT'D

################################################
########## population-level variation ########## 
################################################
## check Rt and mean Ct over time
detected_id <- unique(detect_cases$i) ## keep same individuals
##
##
rt_est_list <- NULL
for (k in 1:4){
  ct0 <- vl_list[[k]] %>% filter(i%in%detected_id) %>%
    mutate(sampled_time=onset_time+confirmation_delay,
           confirmed_time=sampled_time)
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
    mutate(rt_cat=cut(Rt,breaks = c(0,0.5,1,1.5,5))) # include all days 
  #
  rt_est_list[[k]] <- rt_tmp
}
##
####
p_gam <- ggplot(ct_all,aes(x=sampled_time,y=ct.value,
                        color=factor(scenario))) +
  geom_smooth(method="gam") +
  scale_y_reverse(name="Ct value") +
  coord_cartesian(ylim=c(40,10)) +
  scale_color_manual(name="",values = c("gray55","pink3","orange2","darkgreen"),
                     labels = c("Peak fixed at day 5",
                                "Peak before onset",
                                "Peak at onset",
                                "Peak after onset")) +
  scale_x_continuous(name="Days")+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position = c(.6,.8),legend.spacing.y = unit(.1,'cm'))+
  labs(title="A")

ct_skewness <- ct_all %>% group_by(scenario,sampled_time) %>% 
  summarise(skewness=e1071::skewness(ct.value)) %>% ungroup()
p_skew <- ggplot(ct_skewness,
                 aes(x=sampled_time,y=skewness,color=factor(scenario))) +
  geom_line() +
  scale_y_continuous(name="Daily Ct skewness") +
  scale_color_manual(name="",values = c("gray55","pink3","orange2","darkgreen"),
                     labels = c("Peak fixed at day 5",
                                "Peak before onset",
                                "Peak at onset",
                                "Peak after onset")) +
  scale_x_continuous(name="Days")+
  theme(legend.position = "none")+
  labs(title="")

rt_est_all <- rbind(rt_est_list[[1]]%>%mutate(scenario=1),
                    rt_est_list[[2]]%>%mutate(scenario=2),
                    rt_est_list[[3]]%>%mutate(scenario=3),
                    rt_est_list[[4]]%>%mutate(scenario=4))
p_rt <- ggplot(rt_est_all,aes(x=sampled_time,y=fit,color=factor(scenario))) +
  geom_line() +
  scale_y_continuous(name="Ct-based Rt") +
  scale_color_manual(name="",values = c("gray55","pink3","orange2","darkgreen"),
                     labels = c("Peak fixed at day 5",
                                "Peak before onset",
                                "Peak at onset",
                                "Peak after onset")) +
  geom_hline(yintercept = 1,linetype = 'dashed',size = 1,color = 'gray77')+
  scale_x_continuous(name="Days")+coord_cartesian(ylim=c(0,5))+
  theme(legend.position = "none")+
  labs(title="")

###
category_used <- c("mean","skewness","rt")
panel_title_all <- c("B",rep("",2),"C",rep("",2),"D",rep("",2),"E",rep("",2))
pp_list <- NULL
pp_list[[1]] <- p_gam
pp_list[[2]] <- p_skew
pp_list[[3]] <- p_rt
for (k in 1:4){ # 4 relative relationships
  for (j in 1:3){ # 3 panels
    pp_list[[3*k+j]] <- plot_overview(rt_est_list[[k]],category_used[j],panel_title_all[3*(k-1)+j])
    if (k == 1 & j == 1){
      pp_list[[3*k+j]] <-  pp_list[[3*k+j]] + theme(legend.position = c(.3,.9))
    } else {
      pp_list[[3*k+j]] <-  pp_list[[3*k+j]] + theme(legend.position = "none")
    }
  }
}


pp_all <- grid.arrange(grobs=pp_list,
                       heights = c(1, 1, 1),
                       widths = c(2,1,1,1,1),
                       layout_matrix = rbind(c(1,4,7,10,13),
                                             c(2,5,8,11,14),
                                             c(3,6,9,12,15)))

ggsave("output/fig_s6_5feb.pdf",pp_all,height = 8,width=13)

####