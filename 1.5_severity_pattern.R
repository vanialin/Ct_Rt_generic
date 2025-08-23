compare_mat <- matrix(NA,16,9)

for (k in 1:2){ # severity scenarios
  for (i in 1:2){ # delay scenarios
    number_mat <- matrix(NA,100,9)
    dat_sev1 <- out_list[[4*(i-1)+7+3*(k==2)]] %>% filter(sampled_time>=110)
    for (n in 1:100){
      top_day <- with(dat_sev1,which.max(count))
      days <- dat_sev1$sampled_time[top_day]+c(-21:21)
      model_new <- lm(log(Rt)~mean+skewness,data=dat_sev1%>%filter(sampled_time%in%days))
      dat_sev1$fit_new <- exp(predict(model_new,newdata=dat_sev1))
      
      set.seed(n)
      dat_sev_renew <- dat_sev1[sample(1:nrow(dat_sev1),nrow(dat_sev1),replace = TRUE),]
      #
      number_mat[n,1] <- Rt_auc(dat_sev_renew)[2]
      number_mat[n,2] <- with(dat_sev_renew,cor(log(Rt),log(fit),method = "spearman"))
      number_mat[n,3] <- with(dat_sev_renew,mean(abs(log(fit)-log(Rt))))
      
      ## updated
      number_mat[n,4] <- Rt_auc(dat_sev_renew%>%dplyr::select(-fit)%>%rename(fit=fit_new))[2]
      number_mat[n,5] <- with(dat_sev_renew%>%dplyr::select(-fit)%>%rename(fit=fit_new),cor(log(Rt),log(fit),method = "spearman"))
      number_mat[n,6] <- with(dat_sev_renew%>%dplyr::select(-fit)%>%rename(fit=fit_new),mean(abs(log(fit)-log(Rt))))
      
      number_mat[n,7] <- number_mat[n,4] - number_mat[n,1]
      number_mat[n,8] <- number_mat[n,5] - number_mat[n,2]
      number_mat[n,9] <- number_mat[n,6] - number_mat[n,3]
    }
    ##
    ##
    
    for (d in 1:9){
      compare_mat[9*(i-1)+4*(k-1)+1:3,d] <- quantile(number_mat[,d],c(0.5,0.025,0.975))
    }
  }
}

compare_mat <- round(compare_mat,2)

####
out_mat <- matrix(NA,6,6)
for (n in 1:2){ # delay layer
  for (m in 1:2){ # severity layer
    for (k in 1:3){ # evaluation matrix
      out_mat[3*(n-1)+k,3*(m-1)+1] <- paste0(compare_mat[9*(n-1)+4*(m-1)+1,k]," (",
                                             compare_mat[9*(n-1)+4*(m-1)+2,k],",",
                                             compare_mat[9*(n-1)+4*(m-1)+3,k],")")
      #
      out_mat[3*(n-1)+k,3*(m-1)+2] <- paste0(compare_mat[9*(n-1)+4*(m-1)+1,3+k]," (",
                                             compare_mat[9*(n-1)+4*(m-1)+2,3+k],",",
                                             compare_mat[9*(n-1)+4*(m-1)+3,3+k],")")
      #
      out_mat[3*(n-1)+k,3*(m-1)+3] <- paste0(compare_mat[9*(n-1)+4*(m-1)+1,6+k]," (",
                                             compare_mat[9*(n-1)+4*(m-1)+2,6+k],",",
                                             compare_mat[9*(n-1)+4*(m-1)+3,6+k],")")
    }
  }
  
}


#write.csv(out_mat,"output/table_sev_supp.csv",row.names = F)

