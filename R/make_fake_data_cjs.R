
# data simulation for single release, diff surv for int and seg groups
siml_cjs<- function(big_phi, big_p, mrkd, remv, n_occ, intgr, segr){
  n_occ<- dim(big_phi)[2]+1
  CH<- as.data.frame(matrix(0, ncol= n_occ, nrow= mrkd))
  colnames(CH)<- c('rel','grj','goj','lmj','mcj','jdj','bon','twx')
  CH$brood<- sample(c('CW','AD'), size=mrkd, prob=c(0.5,0.5), replace=TRUE)
  CH$group<- NA
  CH[CH$brood=='CW',]$group<- sample(c('R','T'), size=sum(CH$brood=='CW'), prob=c(0.3,0.7), replace=TRUE)
  CH[CH$brood=='AD',]$group<- sample(c('R','T'), size=sum(CH$brood=='AD'), prob=c(0.3,0.7), replace=TRUE)
  CH$prob<- ifelse(CH[,n_occ+1]=='CW', intgr/sum(CH[,n_occ+1]=='CW'),
                   segr/sum(CH[,n_occ+1]=='AD'))
  CH$age<- NA
  
  for(i in 1:mrkd){
    CH[i,1]<- 1 # first detection
    for(t in 2:n_occ){ # starting loop on the second occ
      ifelse(CH[i,n_occ+1]=='CW', sur<- rbinom(1, 1, big_phi[i,t-1]),
             sur<- rbinom(1, 1, big_phi[i,t-1]-0.07)) # diff surv for int and seg
      if(sur==0) break
      rp<- rbinom(1, 1, big_p[i,t-1]) # detection
      if(rp==1) CH[i,t]<- 1
      if(t==n_occ) next
      rmvd<- rbinom(1, 1, remv) # remove xx% of fish
      if(rmvd==1) CH[i,t]<- sample(c(2,3), size=1, prob=c(0.998,0.002))
      if(CH[i,t]==2|CH[i,t]==3) break
    } # observed fate at t for fish i
    if(sur==0) next
    adu<- rbinom(1, 1, 0.02) # x% adult return (for all)
    #ifelse(sum(CH[i,2:n_occ])==0, adu<- rbinom(1,1,0.05), adu<- rbinom(1,1,0.01))
    ifelse(adu==1, CH$age[i]<- sample(c(0,1,2), size=1, prob=c(0.09,0.38,0.53)),
           CH$age[i]<- NA) # assign age
  } # fish i
  CH$ac0<- ifelse(rowSums(CH[,2:4])==0 & CH$age>1, 1, 0)
  CH$ac0j<- ifelse(rowSums(CH[,2:4])==0 & CH$age>0, 1, 0)
  CH$ac1<- ifelse((CH[,2]==1|CH[,3]==1|CH[,4]==1)&
                    CH[,2]!=2& CH[,3]!=2& CH[,4]!=2&
                    CH[,2]!=3& CH[,3]!=3& CH[,4]!=3& CH$age>1, 1, 0)
  CH$ac1j<- ifelse((CH[,2]==1|CH[,3]==1|CH[,4]==1)&
                     CH[,2]!=2& CH[,3]!=2& CH[,4]!=2&
                     CH[,2]!=3& CH[,3]!=3& CH[,4]!=3& CH$age>0, 1, 0)
  CH$atx<- ifelse((CH[,2]==2|CH[,3]==2|CH[,4]==2)& CH$age>1, 1, 0)
  CH$atxj<- ifelse((CH[,2]==2|CH[,3]==2|CH[,4]==2)& CH$age>0, 1, 0)
  CH$c0type<- ifelse(rowSums(CH[,2:4])==0, 1, 0)
  CH$d2<- ifelse(CH[,2]==2|CH[,2]==3, 1, 0)
  CH$d3<- ifelse(CH[,3]==2|CH[,3]==3, 1, 0)
  CH$d4<- ifelse(CH[,4]==2|CH[,4]==3, 1, 0)
  CH$d50<- ifelse(CH$c0type==1& CH[,5]==2|CH[,5]==3, 1, 0)
  CH$d60<- ifelse(CH$c0type==1& CH[,6]==2|CH[,6]==3, 1, 0)
  CH$d70<- ifelse(CH$c0type==1& CH[,7]==2|CH[,7]==3, 1, 0)
  CH$d51<- ifelse(CH$c0type==0& CH[,5]==2|CH[,5]==3, 1, 0)
  CH$d61<- ifelse(CH$c0type==0& CH[,6]==2|CH[,6]==3, 1, 0)
  CH$d71<- ifelse(CH$c0type==0& CH[,7]==2|CH[,7]==3, 1, 0)
  
  return(CH)
}

