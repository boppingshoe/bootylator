
#' Import and format data for ready to use by \code{surv_calc()}
#'
#' @param file_name File path where the input csv file is stored.
#' @param wgt 'y' if using weighted sampling probability.
#' @return Capture history and indicators for adult return.
#' @examples
#' detect_data<- format_dat('C:/Users/bobbyhsu/Documents/Temp/SR HCH 2014 MCCA.csv', wgt='n')
#'
format_dat<- function(file_name, wgt){
  # importing data files and select the wanted columns ----
  # the 'burnham' here is actually the capture_di from the original data file
  yomama_in<- read.csv(file=file_name)#, na.strings= c('','NA'))
  yomama <- subset(yomama_in, , c(1,4,16,17,18,23, grep('flag', names(yomama_in))))
  n_col<- ncol(yomama)
  if(n_col==6) {
    names(yomama)<- c("tagId","burnham","twx","boa","return","relDate")
    yomama$group<- 'T'
  } else if(n_col==7) {
    names(yomama)<- c("tagId","burnham","twx","boa","return","relDate","group")
  } else if(n_col==8) {
    names(yomama)<- c("tagId","burnham","twx","boa","return","relDate","group","brood")
  } else stop('Data file is not read properly. Make sure data source is in the correct format.')
  n_occ<- nchar(yomama$burnham[1])+1
  # ----

  # create detection history ----
  fdat<- as.data.frame(matrix(0, nrow=nrow(yomama), ncol=n_occ))
  fdat[,1]<- 1
  for(t in 2:(n_occ-1)){
    fdat[,t]<- as.numeric(substr(yomama$burnham, t,t))
  }
  fdat[,n_occ]<- ifelse(yomama$twx=='',0,1) # adding TWX as the last detection
  if(n_occ==8) colnames(fdat)<- c('rel','grj','goj','lmj','mcj','jdj','bon','twx')
  # ----

  # add columns before original order is altered ----
  fdat$burnham<- yomama$burnham
  fdat$tagId<- yomama$tagId
  fdat$group<- yomama$group

  if(wgt=='y'){
    trim.trailing <- function (x) sub("\\s+$", "", x)
    fdat$brood<- trim.trailing(yomama$brood)
    fdat$prob<- ifelse(fdat$brood=='CW', intgr/sum(fdat$brood=='CW'),
                          segr/sum(fdat$brood=='AD'))
  } else fdat$prob<- 1/nrow(fdat)

  fdat$relDate<- as.Date(substr(yomama$relDate, 1,10))
  yomama$boa[grepl("^ *$",yomama$boa)]<- NA
  yomama$return[grepl("^ *$",yomama$return)]<- NA
  fdat$boa<- as.Date(substr(yomama$boa, 1,10))
  fdat$return<- as.Date(substr(yomama$return, 1,10))
  # age calculated using BOA_OBS (here is named 'boa')
  fdat$age_boa<- as.numeric(format(fdat$boa, '%Y'))-
    as.numeric(format(fdat$relDate, '%Y'))
  # age calculated using GRA_OBS, MCN_OBS, or BOA_OBS2 (here is named 'return')
  fdat$age_rtn<- as.numeric(format(fdat$return, '%Y'))-
    as.numeric(format(fdat$relDate, '%Y'))
  # ----

  # correct records with detection after 2 or 3 ----
  # (order will be altered after correction)
  correct<- function(x){
    pos<- which(x[1:n_occ]==x[length(x)])[1]
    x[(pos+1):n_occ]<- 0
    return(x[-length(x)])
  }
  numr<- function(tmp) as.numeric(as.character(tmp))

  badId3<- fdat[grep('[3]', fdat$burnham), 'tagId']
  if (length(badId3)>0){
    tmp3<- apply(cbind(subset(fdat, tagId%in%badId3, 1:n_occ), 3), 1, correct)
    tmp3<- as.data.frame(cbind(t(tmp3), subset(fdat, tagId%in%badId3,
                                               (n_occ+1):ncol(fdat)) ))
    tmp3[,1:n_occ]<- apply(tmp3[,1:n_occ], 2, numr)
    fdat<- rbind( tmp3, subset(fdat,!(tagId%in%badId3)) )
  }

  badId2<- fdat[grep('[2]', fdat$burnham), 'tagId']
  badId2<- subset(badId2, !(badId2 %in% badId3))
  if (length(badId2)>0){
    tmp2<- apply(cbind(subset(fdat, tagId%in%badId2, 1:n_occ), 2), 1, correct)
    tmp2<- as.data.frame(cbind(t(tmp2), subset(fdat, tagId%in%badId2,
                                               (n_occ+1):ncol(fdat)) ))
    tmp2[,1:n_occ]<- apply(tmp2[,1:n_occ], 2, numr)
    fdat<- rbind( tmp2, subset(fdat,!(tagId%in%badId2)) )
  }
  # ----

  # tallying using the corrected data set ----
  # adult counts using GRA_OBS, MCN_OBS, or BOA_OBS2 (aka 'return')
  fdat$ac0_rtn<- ifelse(fdat[,2]==0& fdat[,3]==0& fdat[,4]==0& fdat$age_rtn>1, 1, 0)
  fdat$ac0j_rtn<- ifelse(fdat[,2]==0& fdat[,3]==0& fdat[,4]==0& fdat$age_rtn>0, 1, 0)

  fdat$ac1_rtn<- ifelse((fdat[,2]==1|fdat[,3]==1|fdat[,4]==1)&
                          fdat[,2]!=2& fdat[,3]!=2& fdat[,4]!=2&
                          fdat[,2]!=3& fdat[,3]!=3& fdat[,4]!=3&
                          fdat$age_rtn>1, 1, 0)
  fdat$ac1j_rtn<- ifelse((fdat[,2]==1|fdat[,3]==1|fdat[,4]==1)&
                           fdat[,2]!=2& fdat[,3]!=2& fdat[,4]!=2&
                           fdat[,2]!=3& fdat[,3]!=3& fdat[,4]!=3&
                           fdat$age_rtn>0, 1, 0)

  fdat$atx_rtn<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)& fdat$age_rtn>1, 1, 0)
  fdat$atxj_rtn<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)& fdat$age_rtn>0, 1, 0)

  fdat$at0_rtn<- ifelse((fdat[,2]==2|fdat[,2]==0)& (fdat[,3]==2|fdat[,3]==0)
                    & (fdat[,4]==2|fdat[,4]==0)& fdat$age_rtn>1, 1, 0)
  fdat$at0j_rtn<- ifelse((fdat[,2]==2|fdat[,2]==0)& (fdat[,3]==2|fdat[,3]==0)
                     & (fdat[,4]==2|fdat[,4]==0)& fdat$age_rtn>0, 1, 0)
  # adult counts using BOA_OBS (aka 'boa')
  fdat$ac0_boa<- ifelse(fdat[,2]==0& fdat[,3]==0& fdat[,4]==0& fdat$age_boa>1, 1, 0)
  fdat$ac0j_boa<- ifelse(fdat[,2]==0& fdat[,3]==0& fdat[,4]==0& fdat$age_boa>0, 1, 0)

  fdat$ac1_boa<- ifelse((fdat[,2]==1|fdat[,3]==1|fdat[,4]==1)& fdat$age_boa>1, 1, 0)
  fdat$ac1j_boa<- ifelse((fdat[,2]==1|fdat[,3]==1|fdat[,4]==1)& fdat$age_boa>0, 1, 0)

  fdat$atx_boa<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)& fdat$age_boa>1, 1, 0)
  fdat$atxj_boa<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)& fdat$age_boa>0, 1, 0)

  fdat$at0_boa<- ifelse((fdat[,2]==2|fdat[,2]==0)& (fdat[,3]==2|fdat[,3]==0)
                        & (fdat[,4]==2|fdat[,4]==0)& fdat$age_boa>1, 1, 0)
  fdat$at0j_boa<- ifelse((fdat[,2]==2|fdat[,2]==0)& (fdat[,3]==2|fdat[,3]==0)
                         & (fdat[,4]==2|fdat[,4]==0)& fdat$age_boa>0, 1, 0)

  fdat$c0type<- 0
  fdat$c0type[fdat[,2]==0& fdat[,3]==0& fdat[,4]==0]<- 1
  fdat$d2<- ifelse(fdat[,2]==2|fdat[,2]==3, 1, 0)
  fdat$d3<- ifelse(fdat[,3]==2|fdat[,3]==3, 1, 0)
  fdat$d4<- ifelse(fdat[,4]==2|fdat[,4]==3, 1, 0)
  fdat$d50<- ifelse(fdat$c0type==1& fdat[,5]==2|fdat[,5]==3, 1, 0)
  fdat$d60<- ifelse(fdat$c0type==1& fdat[,6]==2|fdat[,6]==3, 1, 0)
  fdat$d70<- ifelse(fdat$c0type==1& fdat[,7]==2|fdat[,7]==3, 1, 0)
  fdat$d51<- ifelse(fdat$c0type==0& fdat[,5]==2|fdat[,5]==3, 1, 0)
  fdat$d61<- ifelse(fdat$c0type==0& fdat[,6]==2|fdat[,6]==3, 1, 0)
  fdat$d71<- ifelse(fdat$c0type==0& fdat[,7]==2|fdat[,7]==3, 1, 0)
  # ----

  return(fdat)
}


#' Estimate survivals, detections, and tally adult returns
#'
#' @param ch Input file made by \code{format_dat()} function.
#' @param i Iteration number used by the bootstrap function \code{bootystrapper()}.
#' @param wt Indicates whether to weight the sampling probability.
#' @param wt_i Indicates whether to calculate the original estimates using weighted probability.
#' @return Survivals, detection and returing adult counts
#' @examples
#' for(i in 1:10){
#'   results<- surv_calc(detect_data, i, nocc=8, wt='y', wt_i='y')
#' }
#' results
#'
surv_calc<- function(ch, i, nocc, wt, wt_i){
  # breakdown of int_t, int_r, seg_t, and seg_r ----
  # if comment out, make sure change the output in 'bootystrapper()'
  if(wt=='y') tnr<- unlist(tapply(ch$group, ch$brood, table))
  else tnr<- table(ch$group)
  # ----
  if(i>1) wt_i<- 'n'
  # m-array is constructed using crt group
  if(wt_i=='y'){
    sim_mary<- marray_wtd(ch, nocc)
  } else {sim_mary<- marray(ch, nocc)}
  # elements needed for estimating phi's and p's ----
  m<- c(0, colSums(sim_mary[, 2:nocc]))
  # if(nocc==4) {z<- c(0, sum(sim_mary[1, 3:nocc]), sum(sim_mary[1:2, 4]))}
  # if(nocc==6) {z<- c(0, sum(sim_mary[1, 3:nocc]), sum(sim_mary[1:2, 4:nocc]),
  #                   sum(sim_mary[1:3, 5:nocc]), sum(sim_mary[1:4, 6]))}
  # if(nocc==8) {z<- c(0, sum(sim_mary[1, 3:nocc]), sum(sim_mary[1:2, 4:nocc]),
  #                   sum(sim_mary[1:3, 5:nocc]), sum(sim_mary[1:4, 6:nocc]),
  #                   sum(sim_mary[1:5, 7:nocc]), sum(sim_mary[1:6, nocc]))}

  z<- rep(0, (nocc-1))
  for (i in 1:(nocc-2)){
    z[i+1]<- sum(sim_mary[1:i,(i+2):nocc])
  }
  R<- sim_mary[1:(nocc-1), 1]
  r<- sim_mary[1:(nocc-1), (nocc+1)]
  if (any(r==0)) {
    M<- z*(R+1)/(r+1) + m[1:(nocc-1)] # finite population correction
  } else M<- z*(R)/(r) + m[1:(nocc-1)]
  # ----

  # output params
  # params from crt group ----
  R1<- sim_mary[1,1]
  c0a_rtn  <- sum(ch$ac0_rtn, na.rm=TRUE)
  c0aj_rtn <- sum(ch$ac0j_rtn, na.rm=TRUE)
  c0a_boa  <- sum(ch$ac0_boa, na.rm=TRUE)
  c0aj_boa <- sum(ch$ac0j_boa, na.rm=TRUE)
  d5670<- colSums (cbind(ch$d50, ch$d60, ch$d70))
  m12<- sim_mary[1,2]
  m13<- sim_mary[1,3]
  m14<- sim_mary[1,4]

  phi<- M[2:(nocc-1)]/(M[1:(nocc-2)] - m[1:(nocc-2)] + R[1:(nocc-2)])
  p  <- m[2:(nocc-1)]/M[2:(nocc-1)]
  # ----
  # params from t group ----
  R1t<- sum(ch$group=='T') # t group

  # m12t<- nrow(subset(ch, group=='T'& ch[,2]!=0)) # pre 2006 migration year
  # m13t<- nrow(subset(ch, group=='T'& ch[,2]==0& goj!=0))
  # m14t<- nrow(subset(ch, group=='T'& ch[,2]==0& goj==0& lmj!=0))

  cht<- subset(ch, group=='T')
  x_t<- cbind(nrow(cht[cht[,2]==2,]), nrow(cht[cht[,3]==2,]),
              nrow(cht[cht[,4]==2,]), nrow(cht[cht[,5]==2,])) # t group
  x_0<- cbind(nrow(cht[cht[,2]==0&cht[,3]==2,]),
              nrow(cht[cht[,2:3]==0&cht[,4]==2,]),
              nrow(cht[cht[,2:4]==0&cht[,5]==2,])) # t group
  d234t<- colSums (cbind(cht$d2, cht$d3, cht$d4)) # t group
  d5671t<- colSums (cbind(cht$d51, cht$d61, cht$d71)) # t group

  c0at_rtn <- sum(ch[ch$group=='T', 'ac0_rtn'], na.rm=TRUE)  # t group
  c0ajt_rtn<- sum(ch[ch$group=='T', 'ac0j_rtn'], na.rm=TRUE) #
  c1at_rtn <- sum(ch[ch$group=='T', 'ac1_rtn'], na.rm=TRUE)  #
  c1ajt_rtn<- sum(ch[ch$group=='T', 'ac1j_rtn'], na.rm=TRUE) #
  txat_rtn <- sum(ch[ch$group=='T', 'atx_rtn'], na.rm=TRUE)  #
  txajt_rtn<- sum(ch[ch$group=='T', 'atxj_rtn'], na.rm=TRUE) #
  t0at_rtn <- sum(ch[ch$group=='T', 'at0_rtn'], na.rm=TRUE)  #
  t0ajt_rtn<- sum(ch[ch$group=='T', 'at0j_rtn'], na.rm=TRUE) #
  c0at_boa <- sum(ch[ch$group=='T', 'ac0_boa'], na.rm=TRUE)  #
  c0ajt_boa<- sum(ch[ch$group=='T', 'ac0j_boa'], na.rm=TRUE) #
  c1at_boa <- sum(ch[ch$group=='T', 'ac1_boa'], na.rm=TRUE)  #
  c1ajt_boa<- sum(ch[ch$group=='T', 'ac1j_boa'], na.rm=TRUE) #
  txat_boa <- sum(ch[ch$group=='T', 'atx_boa'], na.rm=TRUE)  #
  txajt_boa<- sum(ch[ch$group=='T', 'atxj_boa'], na.rm=TRUE) #
  t0at_boa <- sum(ch[ch$group=='T', 'at0_boa'], na.rm=TRUE)  #
  t0ajt_boa<- sum(ch[ch$group=='T', 'at0j_boa'], na.rm=TRUE) #

  lgr_atx_rtn<-  sum(ch[ch$group=='T'& ch[,2]==2, 'atx_rtn'], na.rm=TRUE)   # t group
  lgs_atx_rtn<-  sum(ch[ch$group=='T'& ch[,3]==2, 'atx_rtn'], na.rm=TRUE)   #
  lmn_atx_rtn<-  sum(ch[ch$group=='T'& ch[,4]==2, 'atx_rtn'], na.rm=TRUE)   #
  lgr_atxj_rtn<- sum(ch[ch$group=='T'& ch[,2]==2, 'atxj_rtn'], na.rm=TRUE)  #
  lgs_atxj_rtn<- sum(ch[ch$group=='T'& ch[,3]==2, 'atxj_rtn'], na.rm=TRUE)  #
  lmn_atxj_rtn<- sum(ch[ch$group=='T'& ch[,4]==2, 'atxj_rtn'], na.rm=TRUE)  #
  # ----

  if(length(tnr)==1) {
    calc<- cbind(t(phi), t(p), R1, R1t, m12, m13, m14,
                 #m12t, m13t, m14t,
                 x_t, x_0, t(d234t), t(d5671t), t(d5670),
                 c0a_rtn, c0aj_rtn, c0a_boa, c0aj_boa,
                 c0at_rtn, c0ajt_rtn, c1at_rtn, c1ajt_rtn,
                 txat_rtn, txajt_rtn, t0at_rtn, t0ajt_rtn,
                 c0at_boa, c0ajt_boa, c1at_boa, c1ajt_boa,
                 txat_boa, txajt_boa, t0at_boa, t0ajt_boa,
                 lgr_atx_rtn, lgs_atx_rtn, lmn_atx_rtn,
                 lgr_atxj_rtn, lgs_atxj_rtn, lmn_atxj_rtn)
  } else {calc<- cbind(t(phi), t(p), R1, R1t, m12, m13, m14,
               #m12t, m13t, m14t,
               x_t, x_0, t(d234t), t(d5671t), t(d5670),
               c0a_rtn, c0aj_rtn, c0a_boa, c0aj_boa,
               c0at_rtn, c0ajt_rtn, c1at_rtn, c1ajt_rtn,
               txat_rtn, txajt_rtn, t0at_rtn, t0ajt_rtn,
               c0at_boa, c0ajt_boa, c1at_boa, c1ajt_boa,
               txat_boa, txajt_boa, t0at_boa, t0ajt_boa,
               lgr_atx_rtn, lgs_atx_rtn, lmn_atx_rtn,
               lgr_atxj_rtn, lgs_atxj_rtn, lmn_atxj_rtn,
               t(tnr))}
  return(calc)
}


#' Bootstrap using surv_calc and organize output
#'
#' @param d Input file made by \code{format_dat()}.
#' @param fn Function to run the bootstrap on.
#' @param iter Amount of bootstrap iterations.
#' @param n_occ Total detection events including the trawl.
#' @param wgt Indicates whether to weight the sampling probability.
#' @param wgt_int Indicates whether to calculate the original estimates using weighted probability.
#' @return Estimates in a data frame with original estimate as the first row and bootstrap results in the remaining rows.
#' @examples
#' out<- bootystrapper(detect_data, surv_calc, iter=100, n_occ=8, wgt='n', wgt_init='n')
#' head(out)
#'
bootystrapper <- function(d, fn, iter, n_occ, wgt, wgt_init){
  start_time<- Sys.time()

  original <- fn(d, i=1, n_occ, wgt, wgt_init) #run function on original data
  #make an output matrix with NA's
  out <- matrix(data=NA, nrow=(iter+1),ncol=length(original))
  # first row name is original
  rownames(out) <- as.character(c("original",c(1:(iter))))
  # fill first line of output matrix with original run
  out[1,] <- original

  # builds an vector of row numbers in original
  index <- c(1:length(d$prob))
  #starts loop for number of iterations
  pb <- txtProgressBar(min=2, max=(iter+1), char='x', width=50, style = 3)
  for (i in (2:(iter+1))){
    # resample index is "resampled" data rows to use
    if (wgt=='y') sample_index <- sample(index, prob=d$prob, replace=T)
    else sample_index <- sample(index, replace=T)
    # build resampled data from sample.index
    sample_data <- d[sample_index,]
    # run function on resampled data
    out[i,] <- fn(sample_data, i, n_occ, wgt, wgt_init)
    # print booty progress
    # if ((i-1) %in% seq(0,iter, by=50)) cat(i-1, ' ')
    setTxtProgressBar(pb, i)
  } # bootstrap loop
  # build output matrix and return
  close(pb)
  out <- as.data.frame(out)

  if(wgt=='y') {
    colnames(out) <- c('phi1', 'phi2', 'phi3', 'phi4', 'phi5', 'phi6', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'AD_R', 'AD_T', 'CW_R', 'CW_T')
  } else if(n_occ== 4) {
    colnames(out) <- c('phi1', 'phi2', 'p2', 'p3', 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn')
  } else if(n_occ==6) {
    colnames(out) <- c('phi1', 'phi2', 'phi3', 'phi4', 'p2', 'p3', 'p4', 'p5', 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn')
  } else {
    colnames(out) <- c('phi1', 'phi2', 'phi3', 'phi4', 'phi5', 'phi6', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'R group', 'T group')
    } # n_occ= 8 and not weighted

  cat('\n')
  print(Sys.time()- start_time)
  return(out)
}


#' Construct m-ij array
#'
#' @param CH Input file made by \code{format_dat()} function.
#' @param n_occ Total detection events including the trawl.
#' @return m-ij array
#' @examples
#' marray(detect_data, n_occ=8)
#'
marray<- function(CH, n_occ){
  n_ind<- dim(CH)[1]
  m_array<- matrix(data= 0, nrow= n_occ, ncol= n_occ+1)
  for(t in 1:n_occ){
    m_array[t,1]<- sum(CH[CH[,t]==1,t])
  } # R(i)

  for(t in 1:(n_occ-1)){
    m_array[t,(t+1)]<- nrow(CH[CH[,t]==1& CH[,(t+1)]!=0,])
    if (t> n_occ-2) next
    m_array[t,(t+2)]<- nrow(CH[CH[,t]==1& CH[,(t+1)]==0& CH[,(t+2)]!=0,])
    if (t> n_occ-3) next
    for(u in (t+3):n_occ){
      m_array[t,u]<- nrow(CH[CH[,t]==1& rowSums(CH[,(t+1):(u-1)])==0& CH[,u]!=0,])
    }
  } # m(ij)

  for(t in 1:n_occ){
    m_array[t, n_occ+1]<- sum(m_array[t,2:n_occ])
  } # r(i)
  out<- m_array[-n_occ,]
  return(out)
}


#' Construct m-ij array with weighted probability
#'
#' @param CH Input file made by \code{format_dat()} function.
#' @param n_occ Total detection events including the trawl.
#' @return m-ij array
#' @details It needs pre-assigned weights to calculate weighted m-ij. If all weights are the same, the results are the same as the unweighted m-ij.
#' @examples
#' # pre-assign weights for intergrated and segregated populations
#' intgr<- 234012/(234012+813874) # 2014 CW
#' segr <- 813874/(234012+813874) # 2014 AD
#' # make input file, specify weighted sampling
#' wgt_data<- format_dat('C:/Users/bobbyhsu/Documents/Temp/SR HCH 2014 MCCA.csv', wgt='y')
#' marray(wgt_data, n_occ=8)
#'
marray_wtd<- function(CH, n_occ){
  n_ind<- dim(CH)[1]
  wtd<- CH[,'prob']*n_ind
  m_array<- matrix(data= 0, nrow= n_occ, ncol= n_occ+1)
  for(t in 1:n_occ){
    m_array[t,1]<- sum(wtd[CH[,t]==1])
  } # R(i)

  for(t in 1:(n_occ-1)){
    m_array[t,(t+1)]<- sum(wtd[CH[,t]==1& CH[,(t+1)]!=0])
    if (t> n_occ-2) next
    m_array[t,(t+2)]<- sum(wtd[CH[,t]==1& CH[,(t+1)]==0& CH[,(t+2)]!=0])
    if (t> n_occ-3) next
    for(u in (t+3):n_occ){
      m_array[t,u]<- sum(wtd[CH[,t]==1& rowSums(CH[,(t+1):(u-1)])==0& CH[,u]!=0])
    }
  } # m(ij)

  for(t in 1:n_occ){
    m_array[t, n_occ+1]<- sum(m_array[t,2:n_occ])
  } # r(i)
  out<- m_array[-n_occ,]
  return(out)
}

































