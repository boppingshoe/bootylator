
#' Import and format data for ready to use by \code{surv_calc()}
#'
#' @param file_name File path where the input csv file is stored.
#' @param wgt "y" if using weighted sampling probability. User will be prompted to enter the amount of intergrated and segregated fish.
#' @return Capture history and indicators for adult return.
#' @examples
#' detect_data<- format_dat('C:/Users/bobbyhsu/Documents/Temp/SR HCH 2014 MCCA.csv', wgt='n')
#'
format_dat<- function(file_name, wgt){
  # importing data files and select the wanted columns ----
  # the 'burnham' here is actually the capture_di from the original data file
  yomama_in<- read.csv(file=file_name)#, na.strings= c('','NA'))
  if(names(yomama_in)[1]!='tag_id') {
    yomama_in<- read.csv(file=file_name, header = FALSE)
    names(yomama_in)<- c('tag_id', 'burnham_hi', 'length', 'capture_di', 'mort',
      'transport', 'GRJ_OBS',	'GRX_OBS', 'GOJ_OBS',	'LMJ_OBS', 'ICH_OBS',	'MCJ_OBS',
      'JDJ_OBS',	'BON_OBS', 'MCA_OBS',	'TWX_OBS', 'BOA_OBS',	'GRA_OBS', 'srrt',
      'tag_site', 'rel_site', 'coord_id',	'rel_date',	'river_km',	'migr_yr', 'flag',
      'tag_date', 'tag_file', 'wt', 'flags', 'tag_rem')
  }
  # yomama <- subset(yomama_in, , c(1,4,16,17,18,23, grep('flag', names(yomama_in))))
  yomama <- subset(yomama_in, , c(1,4,16,grep('BOA', names(yomama_in))[1],ifelse(grepl('BOA', names(yomama_in)[18]), 17, 18),23, grep('flag', names(yomama_in))))
  n_col<- ncol(yomama)
  if(n_col==6) {
    names(yomama)<- c("tag_id","burnham","twx","boa","return","rel_date")
    yomama$group<- 'T'
  } else if(n_col==7) {
    names(yomama)<- c("tag_id","burnham","twx","boa","return","rel_date","group")
  } else if(n_col==8) {
    names(yomama)<- c("tag_id","burnham","twx","boa","return","rel_date","group","brood")
  } else stop('Data file is not read properly. Make sure data source is in the correct format.')
  n_occ<- nchar(yomama$burnham[1])+1
  # ----

  # create detection history ----
  fdat<- as.data.frame(matrix(0, nrow=nrow(yomama), ncol=n_occ))
  fdat[,1]<- 1
  for(t in 2:(n_occ-1)){
    fdat[,t]<- as.numeric(substr(yomama$burnham, t, t))
  }
  # adding TWX as the last detection
  fdat[,n_occ]<- ifelse(yomama$twx==''|is.na(yomama$twx), 0, 1)
  colnames(fdat)<- paste0('occ',1:n_occ)
  # if(n_occ==8) colnames(fdat)<- c('rel','grj','goj','lmj','mcj','jdj','bon','twx')
  # ----

  # add columns before original order is altered ----
  fdat$burnham<- apply(fdat[,1:n_occ], 1 , function(x) paste0(x, collapse=''))
  fdat$tag_id<- yomama$tag_id
  fdat$group<- yomama$group

  if(wgt=='y'){
    trim.trailing <- function (x) sub("\\s+$", "", x)
    fdat$brood<- trim.trailing(yomama$brood)
    cw<- as.numeric(readline(prompt = 'Integrated: '))
    ad<- as.numeric(readline(prompt = 'Segregated: '))
    intgr<- cw/(cw+ad)
    segr<- ad/(cw+ad)
    fdat$prob<- ifelse(fdat$brood=='CW', intgr/sum(fdat$brood=='CW'),
                          segr/sum(fdat$brood=='AD'))
  } else fdat$prob<- 1/nrow(fdat)

  fdat$rel_date<- as.Date(substr(yomama$rel_date, 1,10))
  yomama$boa[grepl("^ *$",yomama$boa)]<- NA
  yomama$return[grepl("^ *$",yomama$return)]<- NA
  fdat$boa<- as.Date(substr(yomama$boa, 1,10))
  fdat$return<- as.Date(substr(yomama$return, 1,10))
  # grab migration year from file name
  migyr<- as.numeric(regmatches(file_name, regexpr("[0-9]...", file_name)))
  # age calculated using BOA_OBS (here is named 'boa')
  fdat$age_boa<- as.numeric(format(fdat$boa, '%Y'))- migyr
  # age calculated using GRA_OBS, MCN_OBS, or BOA_OBS2 (here is named 'return')
  fdat$age_rtn<- as.numeric(format(fdat$return, '%Y'))- migyr
  # ----

  # correct records with detection after 2 or 3 ----
  # (order will be altered after correction)
  correct<- function(x){
    pos<- which(x[1:n_occ]==2|x[1:n_occ]==3)[1]
    x[(pos+1):n_occ]<- 0
    return(x)
  }

  tempset<- fdat[grep('[23]', fdat$burnham),]
  if(nrow(tempset)> 0){
    posi<- apply(tempset[,1:n_occ], 1, function(x) grep('[23]', x)[1])

    badId23<- tempset[grepl('[123]', substr(tempset$burnham, posi+1, n_occ)), 'tag_id']
    badId23<- badId23[!is.na(badId23)]
    if(length(badId23)> 0){
      tmp23<- apply(subset(fdat, tag_id%in%badId23, 1:n_occ), 1, correct)
      tmp23<- as.data.frame(cbind(t(tmp23), subset(fdat, tag_id%in%badId23,
        (n_occ+1):ncol(fdat)) ))
      fdat<- rbind( tmp23, subset(fdat,!(tag_id%in%badId23)) )

      posi2<- apply(tmp23[,1:n_occ], 1, function(x) grep('[23]', x))
      qnable<- tmp23[as.numeric(substr(tmp23$burnham, posi2+1, n_occ))!=1, ]

      if(nrow(qnable)> 0){
        toshow<- readline(prompt= paste('I found', nrow(qnable),
            'questionable fish. Would you like to see the list (y/shush)? '))
        if(toshow=='y') print(qnable[, c(paste0('occ',1:n_occ),
        'burnham', 'tag_id', 'group', 'rel_date')], max.print=1e+06)
      }
    }
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

  fdat$atx_rtn<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)&
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_rtn>1, 1, 0)
  fdat$atxj_rtn<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)&
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_rtn>0, 1, 0)

  fdat$at0_rtn<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2) &
      fdat[,2]!=1 & fdat[,3]!=1 & fdat[,4]!=1 &
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = '')) &
      fdat$age_rtn>1, 1, 0)
  fdat$at0j_rtn<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2) &
      fdat[,2]!=1 & fdat[,3]!=1 & fdat[,4]!=1 &
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_rtn>0, 1, 0)
  # adult counts using BOA_OBS (aka 'boa')
  fdat$ac0_boa<- ifelse(fdat[,2]==0& fdat[,3]==0& fdat[,4]==0& fdat$age_boa>1, 1, 0)
  fdat$ac0j_boa<- ifelse(fdat[,2]==0& fdat[,3]==0& fdat[,4]==0& fdat$age_boa>0, 1, 0)

  fdat$ac1_boa<- ifelse((fdat[,2]==1|fdat[,3]==1|fdat[,4]==1)&
      fdat[,2]!=2& fdat[,3]!=2& fdat[,4]!=2&
      fdat[,2]!=3& fdat[,3]!=3& fdat[,4]!=3&
      fdat$age_boa>1, 1, 0)
  fdat$ac1j_boa<- ifelse((fdat[,2]==1|fdat[,3]==1|fdat[,4]==1)&
      fdat[,2]!=2& fdat[,3]!=2& fdat[,4]!=2&
      fdat[,2]!=3& fdat[,3]!=3& fdat[,4]!=3&
      fdat$age_boa>0, 1, 0)

  fdat$atx_boa<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)&
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_boa>1, 1, 0)
  fdat$atxj_boa<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2)&
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_boa>0, 1, 0)

  fdat$at0_boa<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2) &
      fdat[,2]!=1 & fdat[,3]!=1 & fdat[,4]!=1 &
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_boa>1, 1, 0)
  fdat$at0j_boa<- ifelse((fdat[,2]==2|fdat[,3]==2|fdat[,4]==2) &
      fdat[,2]!=1 & fdat[,3]!=1 & fdat[,4]!=1 &
      substr(fdat$burnham,2,(n_occ-1))==
      apply(fdat[,2:(n_occ-1)],1, function(x)paste(x,collapse = ''))&
      fdat$age_boa>0, 1, 0)

  fdat$c0type<- 0
  fdat$c0type[fdat[,2]==0& fdat[,3]==0& fdat[,4]==0]<- 1
  fdat$d2<- ifelse(fdat[,2]==2|fdat[,2]==3, 1, 0)
  fdat$d3<- ifelse(fdat[,3]==2|fdat[,3]==3, 1, 0)
  fdat$d4<- ifelse(fdat[,4]==2|fdat[,4]==3, 1, 0)
  fdat$d50<- ifelse(fdat$c0type==1& (fdat[,5]==2|fdat[,5]==3), 1, 0)
  fdat$d60<- ifelse(fdat$c0type==1& (fdat[,6]==2|fdat[,6]==3), 1, 0)
  fdat$d70<- ifelse(fdat$c0type==1& (fdat[,7]==2|fdat[,7]==3), 1, 0)
  fdat$d51<- ifelse(fdat$c0type==0& (fdat[,5]==2|fdat[,5]==3), 1, 0)
  fdat$d61<- ifelse(fdat$c0type==0& (fdat[,6]==2|fdat[,6]==3), 1, 0)
  fdat$d71<- ifelse(fdat$c0type==0& (fdat[,7]==2|fdat[,7]==3), 1, 0)
  # ----

  # add identifiers
  fdat$tag_site<- yomama_in$tag_site
  fdat$rel_site<- yomama_in$rel_site
  fdat$coord_id<- yomama_in$coord_id

  return(fdat)
}


#' Estimate survivals, detections, and tally adult returns
#'
#' @param ch Input file made by \code{format_dat()} function.
#' @param i Iteration number used by the bootstrap function \code{bootystrapper()}.
#' @param nocc Total detection events including the trawl.
#' @param wt Indicates whether to weight the sampling probability.
#' @param wt_i Indicates whether to calculate the original estimates using weighted probability.
#' @param phi_p_only Option to only calculate survivals and detection and not do the adult counts. Default to no ("n") if not specified.
#' @param fpc Option to adapt finite population correction for survival and detection calculations. Default to no ("n") if not specified.
#' @param match_bt4 Option to follow the procedures used by BT4 program. BT4 excludes mini-jacks for both adult and juvenile removal counts. Default to yes ("y") if not specified. If one choosed no ("n"), mini-jacks would be included in the juvenile removal counts.
#' @return Survivals, detection and returing adult counts
#' @examples
#' for(i in 1:10){
#'   results<- surv_calc(detect_data, i, nocc=8, wt='y', wt_i='y')
#' }
#' results
#'
surv_calc<- function(ch, i, nocc, wt, wt_i, phi_p_only, fpc, match_bt4, ...){
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
  z<- rep(0, (nocc-1))
  for (i in 1:(nocc-2)){
    z[i+1]<- sum(sim_mary[1:i,(i+2):nocc])
  }
  R<- sim_mary[1:(nocc-1), 1]
  r<- sim_mary[1:(nocc-1), (nocc+1)]
  if (fpc=='y' & any(r==0)) {
    M<- z*(R+1)/(r+1) + m[1:(nocc-1)] # finite population correction
  } else M<- z*(R)/(r) + m[1:(nocc-1)]

  phi<- M[2:(nocc-1)]/(M[1:(nocc-2)] - m[1:(nocc-2)] + R[1:(nocc-2)])
  p  <- m[2:(nocc-1)]/M[2:(nocc-1)]
  if(phi_p_only=='y') {
    calc<- cbind(t(phi), t(p))
    return(calc)
    stop()
  }
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
  # ----
  # params from t group ----
  R1t<- sum(ch$group=='T') # t group

  # m12t<- nrow(subset(ch, group=='T'& ch[,2]!=0)) # pre 2006 migration year
  # m13t<- nrow(subset(ch, group=='T'& ch[,2]==0& goj!=0))
  # m14t<- nrow(subset(ch, group=='T'& ch[,2]==0& goj==0& lmj!=0))

  cht<- subset(ch, group=='T')
  if(match_bt4=='n'){
    x_t<- cbind(nrow(subset(cht, occ2==2 & as.numeric(substr(burnham,3,nocc-1))==0)),
      nrow(subset(cht, occ3==2 & as.numeric(substr(burnham,4,nocc-1))==0)),
      nrow(subset(cht, occ4==2 & as.numeric(substr(burnham,5,nocc-1))==0)),
      ifelse(nocc>4, nrow(subset(cht, occ5==2 &
          as.numeric(substr(burnham,6,nocc-1))==0)), NA)) # t group
  } else {
    # BT4 doesn't count smolts that return as mini-jacks
    # do this to match BT4 counts
    x_t<- cbind(nrow(subset(cht, occ2==2 &
        as.numeric(substr(burnham,3,nocc-1))==0 &
        (age_rtn!=0|is.na(age_rtn)))),
      nrow(subset(cht, occ3==2 &
          as.numeric(substr(burnham,4,nocc-1))==0 &
          (age_rtn!=0|is.na(age_rtn)))),
      nrow(subset(cht, occ4==2 &
          as.numeric(substr(burnham,5,nocc-1))==0 &
          (age_rtn!=0|is.na(age_rtn)))),
      ifelse(nocc>4, nrow(subset(cht, occ5==2 &
          as.numeric(substr(burnham,6,nocc-1))==0 &
          (age_rtn!=0|is.na(age_rtn)))), NA)) # t group
  }

  x_0<- cbind(nrow(cht[cht[,2]==0 & cht[,3]==2,]),
            nrow(cht[as.numeric(substr(cht$burnham,2,3))==0 & cht[,4]==2,]),
            nrow(cht[as.numeric(substr(cht$burnham,2,4))==0 & cht[,5]==2,])) # t group
  d234t<- colSums (cbind(cht$d2, cht$d3, cht$d4)) # t group
  d5671t<- colSums (cbind(cht$d51, cht$d61, cht$d71)) # t group

  c0at_rtn <- sum(cht[, 'ac0_rtn'], na.rm=TRUE)  # t group
  c0ajt_rtn<- sum(cht[, 'ac0j_rtn'], na.rm=TRUE) #
  c1at_rtn <- sum(cht[, 'ac1_rtn'], na.rm=TRUE)  #
  c1ajt_rtn<- sum(cht[, 'ac1j_rtn'], na.rm=TRUE) #
  txat_rtn <- sum(cht[, 'atx_rtn'], na.rm=TRUE)  #
  txajt_rtn<- sum(cht[, 'atxj_rtn'], na.rm=TRUE) #
  t0at_rtn <- sum(cht[, 'at0_rtn'], na.rm=TRUE)  #
  t0ajt_rtn<- sum(cht[, 'at0j_rtn'], na.rm=TRUE) #
  c0at_boa <- sum(cht[, 'ac0_boa'], na.rm=TRUE)  #
  c0ajt_boa<- sum(cht[, 'ac0j_boa'], na.rm=TRUE) #
  c1at_boa <- sum(cht[, 'ac1_boa'], na.rm=TRUE)  #
  c1ajt_boa<- sum(cht[, 'ac1j_boa'], na.rm=TRUE) #
  txat_boa <- sum(cht[, 'atx_boa'], na.rm=TRUE)  #
  txajt_boa<- sum(cht[, 'atxj_boa'], na.rm=TRUE) #
  t0at_boa <- sum(cht[, 'at0_boa'], na.rm=TRUE)  #
  t0ajt_boa<- sum(cht[, 'at0j_boa'], na.rm=TRUE) #

  lgr_atx_rtn<-  sum(cht[cht[,2]==2, 'atx_rtn'], na.rm=TRUE)   # t group
  lgs_atx_rtn<-  sum(cht[cht[,3]==2, 'atx_rtn'], na.rm=TRUE)   #
  lmn_atx_rtn<-  sum(cht[cht[,4]==2, 'atx_rtn'], na.rm=TRUE)   #
  lgr_atxj_rtn<- sum(cht[cht[,2]==2, 'atxj_rtn'], na.rm=TRUE)  #
  lgs_atxj_rtn<- sum(cht[cht[,3]==2, 'atxj_rtn'], na.rm=TRUE)  #
  lmn_atxj_rtn<- sum(cht[cht[,4]==2, 'atxj_rtn'], na.rm=TRUE)  #
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
#' @param wgt Indicates whether to weight the sampling probability.
#' @param wgt_init Indicates whether to calculate the original estimates using weighted probability.
#' @param phi_p_only Indicate to turn off the phi_p_only option in \code{curv_calc()}. Default is no ("n").
#' @param fpc Indicate to turn off the fpc option in \code{curv_calc()}. Default is yes ("y").
#' @param match_bt4 Indicate to turn off the match_bt4 option in \code{curv_calc()}. The default here is yes ("y").
#' @return Estimates in a data frame with original estimate as the first row and bootstrap results in the remaining rows.
#' @examples
#' out<- bootystrapper(detect_data, surv_calc, iter=100, n_occ=8, wgt='n', wgt_init='n')
#' head(out)
#'
bootystrapper <- function(d, fn, iter, wgt, wgt_init, phi_p_only='n', fpc='y', match_bt4='y', ...){
  start_time<- Sys.time()
  n_occ<- sum(grepl('occ', names(d)))
  original <- fn(d, i=1, n_occ, wgt, wgt_init, phi_p_only, fpc, match_bt4) #run function on original data
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
    out[i,] <- fn(sample_data, i, n_occ, wgt, wgt_init, phi_p_only, fpc, match_bt4)
    # print booty progress
    # if ((i-1) %in% seq(0,iter, by=50)) cat(i-1, ' ')
    setTxtProgressBar(pb, i)
  } # bootstrap loop
  close(pb)
  # build output matrix and return
  out <- as.data.frame(out)
  out$tag_site<- d$tag_site[1]
  out$rel_site<- d$rel_site[1]
  out$coord_id<- d$coord_id[1]

  if(phi_p_only== 'y') {
    colnames(out) <- c(paste0('phi',1:(n_occ-2)), paste0('p',2:(n_occ-1)))
  } else if(wgt== 'y') {
    colnames(out) <- c(paste0('phi',1:(n_occ-2)), paste0('p',2:(n_occ-1)), 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'AD_R', 'AD_T', 'CW_R', 'CW_T', 'tag_site', 'rel_site', 'coord_id')
  } else if(n_occ== 8) { # n_occ= 8 and not weighted
    colnames(out) <- c(paste0('phi',1:(n_occ-2)), paste0('p',2:(n_occ-1)), 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'R group', 'T group', 'tag_site', 'rel_site', 'coord_id')
  } else { # n_occ= 4 or 6
    colnames(out) <- c(paste0('phi',1:(n_occ-2)), paste0('p',2:(n_occ-1)), 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'tag_site', 'rel_site', 'coord_id')
  }

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

































