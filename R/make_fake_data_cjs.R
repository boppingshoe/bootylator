
#' data simulation for single release, diff surv for int and seg groups
#'
#' @param big_phi Set reach survivals.
#' @param big_p Set detection probabilities for each reach.
#' @param mrkd Indicates amount of fish released.
#' @param remv Indicates the portion of fish removed (transported or died).
#' @param n_occ Number of events/detection.
#' @param intgr Portion of "intergrated" group. Intergrated + Segregatred = 1.
#' @param surv_diff Difference in survival between the intergrated and segregated groups.
#' @param grp_t Portion of transported group. Group T + Group R = 1
#' @param adu_rtn Adult return rate.
#' @return Simulated detection history and adult counts. Bootstrap ready.
#' @examples
#' n_occ<- 8
#' mrkd<- 5000
#' phi_real<- rep(0.8, n_occ-1) # set survival
#' p_real<- rep(0.45, n_occ-1) # set detection
#' remv<- 0.01 # portion removed at occ 2 and 3
#' intgr<- 0.23
#' surv_diff<- 0.07
#' grp_t<- 0.7
#' adu_rtn<- 0.02
#' big_phi<- matrix(phi_real, ncol= n_occ-1, nrow= mrkd)
#' big_p<- matrix(p_real, ncol= n_occ-1, nrow= mrkd)
#'
#' ch<- siml_cjs(big_phi, big_p, mrkd, remv, n_occ, intgr, surv_diff, grp_t, adu_rtn)
#'
siml_cjs<- function(big_phi, big_p, mrkd, remv, n_occ, intgr, surv_diff, grp_t, adu_rtn){
  # n_occ<- dim(big_phi)[2]+1
  segr<- 1- intgr
  grp_r<- 1- grp_t
  CH<- as.data.frame(matrix(0, ncol= n_occ, nrow= mrkd))
  if(n_occ==8) {colnames(CH)<- c('rel','grj','goj','lmj','mcj','jdj','bon','twx')}
  CH$brood<- sample(c('CW','AD'), size=mrkd, prob=c(0.5,0.5), replace=TRUE)
  CH$group<- NA
  CH[CH$brood=='CW',]$group<- sample(c('R','T'), size=sum(CH$brood=='CW'),
    prob=c(grp_r, grp_t), replace=TRUE)
  CH[CH$brood=='AD',]$group<- sample(c('R','T'), size=sum(CH$brood=='AD'),
    prob=c(grp_r, grp_t), replace=TRUE)
  CH$prob<- ifelse(CH[,n_occ+1]=='CW', intgr/sum(CH[,n_occ+1]=='CW'),
                   segr/sum(CH[,n_occ+1]=='AD'))

  for(i in 1:mrkd){
    CH[i,1]<- 1 # first detection
    for(t in 2:n_occ){ # starting loop on the second occ
      ifelse(CH[i,n_occ+1]=='CW', sur<- rbinom(1, 1, big_phi[i,t-1]),
             sur<- rbinom(1, 1, big_phi[i,t-1]-surv_diff)) # diff surv for int and seg
      if(sur==0) break
      rp<- rbinom(1, 1, big_p[i,t-1]) # detection
      if(rp==1) CH[i,t]<- 1
      if(t==n_occ) next
      rmvd<- rbinom(1, 1, remv) # remove xx% of fish
      if(rmvd==1) CH[i,t]<- sample(c(2,3), size=1, prob=c(0.998,0.002))
      if(CH[i,t]==2|CH[i,t]==3) break
    } # observed fate at t for fish i
    if(sur==0) next
    adu<- rbinom(1, 1, adu_rtn) # x% adult return (for all)
    #ifelse(sum(CH[i,2:n_occ])==0, adu<- rbinom(1,1,0.05), adu<- rbinom(1,1,0.01))
    ifelse(adu==1, CH$age_boa[i]<- sample(c(0,1,2,3), size=1,
      prob=c(0.26, 0.34, 0.39, 0.01)), CH$age_boa[i]<- NA) # assign age
    ifelse(!is.na(CH$age_boa), CH$age_rtn[i]<- sample(c(NA,0,1,2,3), size=1,
      prob=c(0.49, 0.01, 0.22, 0.27, 0.01)), CH$age_rtn[i]<- NA) # assign age
  } # fish i
  # tallying using the corrected data set ----
  # adult counts using age_rtn
  CH$ac0_rtn<- ifelse(CH[,2]==0& CH[,3]==0& CH[,4]==0& CH$age_rtn>1, 1, 0)
  CH$ac0j_rtn<- ifelse(CH[,2]==0& CH[,3]==0& CH[,4]==0& CH$age_rtn>0, 1, 0)

  CH$ac1_rtn<- ifelse((CH[,2]==1|CH[,3]==1|CH[,4]==1)&
      CH[,2]!=2& CH[,3]!=2& CH[,4]!=2&
      CH[,2]!=3& CH[,3]!=3& CH[,4]!=3&
      CH$age_rtn>1, 1, 0)
  CH$ac1j_rtn<- ifelse((CH[,2]==1|CH[,3]==1|CH[,4]==1)&
      CH[,2]!=2& CH[,3]!=2& CH[,4]!=2&
      CH[,2]!=3& CH[,3]!=3& CH[,4]!=3&
      CH$age_rtn>0, 1, 0)

  CH$atx_rtn<- ifelse((CH[,2]==2|CH[,3]==2|CH[,4]==2)& CH$age_rtn>1, 1, 0)
  CH$atxj_rtn<- ifelse((CH[,2]==2|CH[,3]==2|CH[,4]==2)& CH$age_rtn>0, 1, 0)

  CH$at0_rtn<- ifelse((CH[,2]==2|CH[,2]==0)& (CH[,3]==2|CH[,3]==0)
    & (CH[,4]==2|CH[,4]==0)& CH$age_rtn>1, 1, 0)
  CH$at0j_rtn<- ifelse((CH[,2]==2|CH[,2]==0)& (CH[,3]==2|CH[,3]==0)
    & (CH[,4]==2|CH[,4]==0)& CH$age_rtn>0, 1, 0)
  # adult counts age_boa
  CH$ac0_boa<- ifelse(CH[,2]==0& CH[,3]==0& CH[,4]==0& CH$age_boa>1, 1, 0)
  CH$ac0j_boa<- ifelse(CH[,2]==0& CH[,3]==0& CH[,4]==0& CH$age_boa>0, 1, 0)

  CH$ac1_boa<- ifelse((CH[,2]==1|CH[,3]==1|CH[,4]==1)& CH$age_boa>1, 1, 0)
  CH$ac1j_boa<- ifelse((CH[,2]==1|CH[,3]==1|CH[,4]==1)& CH$age_boa>0, 1, 0)

  CH$atx_boa<- ifelse((CH[,2]==2|CH[,3]==2|CH[,4]==2)& CH$age_boa>1, 1, 0)
  CH$atxj_boa<- ifelse((CH[,2]==2|CH[,3]==2|CH[,4]==2)& CH$age_boa>0, 1, 0)

  CH$at0_boa<- ifelse((CH[,2]==2|CH[,2]==0)& (CH[,3]==2|CH[,3]==0)
    & (CH[,4]==2|CH[,4]==0)& CH$age_boa>1, 1, 0)
  CH$at0j_boa<- ifelse((CH[,2]==2|CH[,2]==0)& (CH[,3]==2|CH[,3]==0)
    & (CH[,4]==2|CH[,4]==0)& CH$age_boa>0, 1, 0)

  CH$c0type<- 0
  CH$c0type[CH[,2]==0& CH[,3]==0& CH[,4]==0]<- 1
  CH$d2<- ifelse(CH[,2]==2|CH[,2]==3, 1, 0)
  CH$d3<- ifelse(CH[,3]==2|CH[,3]==3, 1, 0)
  CH$d4<- ifelse(CH[,4]==2|CH[,4]==3, 1, 0)
  CH$d50<- ifelse(CH$c0type==1& CH[,5]==2|CH[,5]==3, 1, 0)
  CH$d60<- ifelse(CH$c0type==1& CH[,6]==2|CH[,6]==3, 1, 0)
  CH$d70<- ifelse(CH$c0type==1& CH[,7]==2|CH[,7]==3, 1, 0)
  CH$d51<- ifelse(CH$c0type==0& CH[,5]==2|CH[,5]==3, 1, 0)
  CH$d61<- ifelse(CH$c0type==0& CH[,6]==2|CH[,6]==3, 1, 0)
  CH$d71<- ifelse(CH$c0type==0& CH[,7]==2|CH[,7]==3, 1, 0)
  # ----

  return(CH)
}

