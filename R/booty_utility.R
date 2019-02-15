
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to bootylator 1.4.4")
}

#' Import and format data for ready to use by \code{surv_calc()}
#'
#' @param file_name File path where the input csv file is stored.
#' @param mig_yr Juvenile migration year.
#' @param wgt "y" if using weighted sampling probability. User will be prompted to enter the amount of intergrated and segregated fish.
#' @return Capture history and indicators for adult return.
#' @examples
#' detect_data<- format_dat('C:/Users/bobbyhsu/Documents/Temp/SR HCH 2014 MCCA.csv', mig_yr= 2014, wgt= 'n')
#'
format_dat<- function(file_name, mig_yr= 'auto', wgt= 'n'){
  # importing data files ----
  yomama_in<- read.csv(file= file_name)#, na.strings= c('','NA'))
  if (sum(grepl('tag_id', names(yomama_in)))!= 1) {
    stop ('Data processing stopped.
      Data column names missing.')
  }  # abort when columns have no names

  # select columns needed ----
  # GRA_OBS for snake, MCA_OBS for others
  yomama<- yomama_in[, c(grep('tag_id', names(yomama_in)),
    grep('capture', names(yomama_in))[1], grep('BOA_OBS', names(yomama_in)),
    grep('MCA_OBS', names(yomama_in)), grep('GRA_OBS', names(yomama_in)),
    grep('flag', names(yomama_in)), grep('rel_date', names(yomama_in)))]
  if (ncol(yomama)== 7) { # should only have 7 columns
  names(yomama)<- c("tag_id","capture","boa","return","group","brood","rel_date")
  } else {
    stop ('Data processing stopped.
      Data file contained both MCA_OBS and GRA_OBS
      (cannot distinguish between Snake and Columbia fish).')
  }
  if (any(names(yomama_in)== 'MCA_OBS')) {
    yomama$group<- 'T'
  }
  if (mig_yr== 'auto') { # grab migration year from file name
    migyr<- as.numeric(regmatches(file_name, regexpr("[0-9]...", file_name)))
  } else { # user set migration year
    migyr<- mig_yr
    if (migyr!= as.numeric(regmatches(file_name, regexpr("[0-9]...", file_name)))) {
      stop ('Data processing stopped.
      Migration year entered did not match the data file name.')
    } # check if migration match file name
  }

  # create detection history ----
  n_occ<- nchar(yomama$capture[1])
  fdat<- as.data.frame(matrix(0, nrow=nrow(yomama), ncol= n_occ))
  fdat[,1]<- 1
  for(t in 2:n_occ){
    fdat[, t]<- as.numeric(substr(yomama$capture, t, t))
  }
  colnames(fdat)<- paste0('occ', 1:n_occ)

  # add columns before original order is altered ----
  # fdat$capture<- apply(fdat[,1:n_occ], 1 , function(x) paste0(x, collapse=''))
  fdat$capture<- yomama$capture
  fdat$tag_id<- yomama$tag_id
  fdat$group<- yomama$group

  if(wgt=='y'){
    fdat$brood<- trimws(yomama$brood)
    cw<- as.numeric(readline(prompt = 'Integrated: '))
    ad<- as.numeric(readline(prompt = 'Segregated: '))
    intgr<- cw/ (cw+ ad)
    segr<- ad/ (cw+ ad)
    fdat$prob<- ifelse(fdat$brood== 'CW', intgr/ sum(fdat$brood== 'CW'),
                          segr/ sum(fdat$brood== 'AD'))
  } else fdat$prob<- 1/ nrow(fdat)

  yomama$boa[grepl("^ *$", yomama$boa)]<- NA # '^' means start of string, '$' end,
  yomama$return[grepl("^ *$", yomama$return)]<- NA # and '*' includes space or empty
  if (grepl('/', substr(yomama$rel_date[1], 1, 4))|
      grepl('-', substr(yomama$rel_date[1], 1, 4))) {
    fdat$rel_date<- as.Date(substr(yomama$rel_date, 1, 10),
      tryFormat= c("%m-%d-%Y", "%m/%d/%Y"))
    fdat$boa<- as.Date(substr(yomama$boa, 1, 10),
      tryFormat= c("%m-%d-%Y", "%m/%d/%Y"))
    fdat$return<- as.Date(substr(yomama$return, 1, 10),
      tryFormat= c("%m-%d-%Y", "%m/%d/%Y"))
  } else {
    fdat$rel_date<- as.Date(substr(yomama$rel_date, 1, 10))
    fdat$boa<- as.Date(substr(yomama$boa, 1, 10))
    fdat$return<- as.Date(substr(yomama$return, 1, 10))
  }
  # age calculated using BOA_OBS (here is named 'boa')
  fdat$age_boa<- as.numeric(format(fdat$boa, '%Y'))- migyr
  # age calculated using GRA_OBS or MCA_OBS (here is named 'return')
  fdat$age_rtn<- as.numeric(format(fdat$return, '%Y'))- migyr

  # correct records with detection after 2 or 3 ----
  # (order will be altered after correction)
  correct<- function(x){
    pos<- which(x[1:n_occ]== 2|x[1:n_occ]== 3)[1]
    x[(pos+1):n_occ]<- 0
    return(x)
  }

  tempset<- fdat[grep('[23]', fdat$capture),]
  if(nrow(tempset)> 0){
    posi<- apply(tempset[,1:n_occ], 1, function(x) grep('[23]', x)[1])

    badId23<- tempset[grepl('[123]', substr(tempset$capture, posi+1, n_occ)), 'tag_id']
    badId23<- badId23[!is.na(badId23)]
    if(length(badId23)> 0){
      tmp23<- apply(subset(fdat, tag_id%in%badId23, 1:n_occ), 1, correct)
      tmp23<- as.data.frame(cbind(t(tmp23), subset(fdat, tag_id%in%badId23,
        (n_occ+ 1):ncol(fdat)) ))
      fdat<- rbind( tmp23, subset(fdat,!(tag_id%in%badId23)) )

      # posi2<- apply(tmp23[, 1:n_occ], 1, function(x) grep('[23]', x))
      # qnable<- tmp23[as.numeric(substr(tmp23$capture, posi2+ 1, n_occ))!= 1, ]
      # if(nrow(qnable)> 0) {
      #   toshow<- readline(prompt= paste('I found', nrow(qnable),
      #       'questionable fish. Would you like to see the list (y/shush)? '))
      #   if(toshow== 'y') print(qnable[, c(paste0('occ', 1:n_occ),
      #   'capture', 'tag_id', 'group', 'rel_date')], max.print= 1e+06)
      # } # Jerry said he's aware of this and requested supressing the message
    }
  }

  # tallying using the corrected data set ----
  # adult counts using GRA_OBS or MCA_OBS (aka 'return')
  fdat$ac0_rtn<- ifelse(fdat[, 2]== 0& fdat[, 3]== 0&
      fdat[, 4]== 0& fdat$age_rtn> 1, 1, 0)
  fdat$ac0j_rtn<- ifelse(fdat[, 2]== 0& fdat[, 3]== 0&
      fdat[, 4]== 0& fdat$age_rtn> 0, 1, 0)

  fdat$ac1_rtn<- ifelse((fdat[, 2]== 1| fdat[, 3]== 1| fdat[, 4]== 1)&
                          fdat[, 2]!= 2& fdat[, 3]!= 2& fdat[, 4]!= 2&
                          fdat[, 2]!= 3& fdat[, 3]!= 3& fdat[, 4]!= 3&
                          fdat$age_rtn> 1, 1, 0)
  fdat$ac1j_rtn<- ifelse((fdat[, 2]==1| fdat[, 3]== 1|fdat[, 4]== 1)&
                           fdat[, 2]!= 2& fdat[, 3]!= 2& fdat[, 4]!= 2&
                           fdat[, 2]!= 3& fdat[, 3]!= 3& fdat[, 4]!= 3&
                           fdat$age_rtn> 0, 1, 0)

  fdat$atx_rtn<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_rtn> 1, 1, 0)
  fdat$atxj_rtn<- ifelse((fdat[ ,2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_rtn> 0, 1, 0)

  fdat$at0_rtn<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      fdat[, 2]!= 1& fdat[, 3]!= 1& fdat[, 4]!= 1&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_rtn> 1, 1, 0)
  fdat$at0j_rtn<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      fdat[, 2]!= 1& fdat[, 3]!= 1& fdat[, 4]!= 1&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_rtn> 0, 1, 0)
  # adult counts using BOA_OBS (aka 'boa')
  fdat$ac0_boa<- ifelse(fdat[, 2]== 0& fdat[, 3]== 0&
      fdat[, 4]== 0& fdat$age_boa> 1, 1, 0)
  fdat$ac0j_boa<- ifelse(fdat[, 2]== 0& fdat[, 3]== 0&
      fdat[, 4]== 0& fdat$age_boa> 0, 1, 0)

  fdat$ac1_boa<- ifelse((fdat[, 2]== 1| fdat[, 3]== 1| fdat[, 4]== 1)&
      fdat[, 2]!= 2& fdat[, 3]!= 2& fdat[, 4]!= 2&
      fdat[, 2]!= 3& fdat[, 3]!= 3& fdat[, 4]!= 3&
      fdat$age_boa> 1, 1, 0)
  fdat$ac1j_boa<- ifelse((fdat[, 2]== 1| fdat[, 3]== 1| fdat[, 4]== 1)&
      fdat[, 2]!= 2& fdat[, 3]!= 2& fdat[, 4]!= 2&
      fdat[, 2]!= 3& fdat[, 3]!= 3& fdat[, 4]!= 3&
      fdat$age_boa> 0, 1, 0)

  fdat$atx_boa<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      substr(fdat$capture, 2, (n_occ-1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_boa> 1, 1, 0)
  fdat$atxj_boa<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_boa> 0, 1, 0)

  fdat$at0_boa<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      fdat[, 2]!= 1& fdat[, 3]!= 1& fdat[, 4]!= 1&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_boa>1, 1, 0)
  fdat$at0j_boa<- ifelse((fdat[, 2]== 2| fdat[, 3]== 2| fdat[, 4]== 2)&
      fdat[, 2]!= 1& fdat[, 3]!= 1& fdat[, 4]!= 1&
      substr(fdat$capture, 2, (n_occ- 1))==
      apply(fdat[, 2:(n_occ- 1)], 1, function(x) paste(x, collapse= ''))&
      fdat$age_boa> 0, 1, 0)

  fdat$c0type<- 0
  fdat$c0type[fdat[, 2]== 0& fdat[, 3]== 0& fdat[, 4]== 0]<- 1
  fdat$d2<- ifelse(fdat[, 2]== 2| fdat[, 2]== 3, 1, 0)
  fdat$d3<- ifelse(fdat[, 3]== 2| fdat[, 3]== 3, 1, 0)
  fdat$d4<- ifelse(fdat[, 4]== 2| fdat[, 4]== 3, 1, 0)
  fdat$d50<- ifelse(fdat$c0type== 1& (fdat[, 5]== 2| fdat[, 5]== 3), 1, 0)
  fdat$d60<- ifelse(fdat$c0type== 1& (fdat[, 6]== 2| fdat[, 6]== 3), 1, 0)
  fdat$d70<- ifelse(fdat$c0type== 1& (fdat[, 7]== 2| fdat[, 7]== 3), 1, 0)
  fdat$d51<- ifelse(fdat$c0type== 0& (fdat[, 5]== 2| fdat[, 5]== 3), 1, 0)
  fdat$d61<- ifelse(fdat$c0type== 0& (fdat[, 6]== 2| fdat[, 6]== 3), 1, 0)
  fdat$d71<- ifelse(fdat$c0type== 0& (fdat[, 7]== 2| fdat[, 7]== 3), 1, 0)

  # add identifiers ----
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
#' @param phi_p_only Option to only calculate survivals and detection and not do the adult counts.
#' @param fpc Option to adapt finite population correction for survival and detection calculations.
#' @return Survivals, detection and returing adult counts
#' @examples
#' surv_calc(detect_data, i= 1, nocc= 8, wt= 'n', wt_i= 'n', phi_p_only= 'y', fpc= 'y')
#'
surv_calc<- function(ch, i, nocc, wt, wt_i, phi_p_only, fpc, ...){
  # breakdown of int_t, int_r, seg_t, and seg_r ----
  # if comment out, make sure change the output in 'bootystrapper()'
  if (wt== 'y') tnr<- unlist(tapply(ch$group, ch$brood, table))
  else tnr<- table(ch$group)
  # ----
  if (i> 1) wt_i<- 'n'
  # m-array is constructed using crt group
  if (wt_i== 'y') sim_mary<- marray_wtd(ch, nocc)
  else sim_mary<- marray(ch, nocc)
  # elements needed for estimating phi's and p's ----
  m<- c(0, colSums(sim_mary[, 2:nocc]))
  z<- rep(0, (nocc- 1))
  for (t in 1:(nocc- 2)) {
    z[t+ 1]<- sum(sim_mary[1:t, (t+2):nocc])
  }
  R<- sim_mary[1:(nocc- 1), 1]
  r<- sim_mary[1:(nocc- 1), (nocc+ 1)]
  if (fpc== 'y' & any(r== 0)) {
    M<- z* (R+ 1)/ (r+ 1)+ m[1:(nocc- 1)] # finite population correction
  } else M<- z* R/ r+ m[1:(nocc- 1)]

  phi<- M[2:(nocc- 1)]/ (M[1:(nocc- 2)]- m[1:(nocc- 2)]+ R[1:(nocc- 2)])
  p  <- m[2:(nocc- 1)]/ M[2:(nocc- 1)]
  if(phi_p_only== 'y') {
    calc<- cbind(t(phi), t(p))
    colnames(calc)<- c(paste0('Phi', 1:(nocc- 2)), paste0('p', 2:(nocc- 1)))
    return(calc)
    stop()
  }
  # ----

  # output params
  # params from crt group ----
  R1<- sim_mary[1, 1]
  c0a_rtn  <- sum(ch$ac0_rtn, na.rm= TRUE)
  c0aj_rtn <- sum(ch$ac0j_rtn, na.rm= TRUE)
  c0a_boa  <- sum(ch$ac0_boa, na.rm= TRUE)
  c0aj_boa <- sum(ch$ac0j_boa, na.rm= TRUE)
  d5670<- colSums (cbind(ch$d50, ch$d60, ch$d70))
  m12<- sim_mary[1, 2]
  m13<- sim_mary[1, 3]
  m14<- sim_mary[1, 4]
  # ----
  # params from t group ----
  R1t<- sum(ch$group== 'T') # t group

  # m12t<- nrow(subset(ch, group== 'T'& ch[, 2]!= 0)) # pre 2006 migration year
  # m13t<- nrow(subset(ch, group== 'T'& ch[, 2]== 0& goj!= 0))
  # m14t<- nrow(subset(ch, group== 'T'& ch[, 2]== 0& goj== 0& lmj!= 0))

  cht<- subset(ch, group== 'T')
  x_t<- cbind(nrow(subset(cht, occ2== 2 &
      as.numeric(substr(capture, 3, nocc- 1))== 0 &
      (age_rtn!= 0| is.na(age_rtn)))),
    nrow(subset(cht, occ3== 2 &
        as.numeric(substr(capture, 4 ,nocc- 1))== 0 &
        (age_rtn!= 0| is.na(age_rtn)))),
    nrow(subset(cht, occ4== 2 &
        as.numeric(substr(capture, 5, nocc- 1))== 0 &
        (age_rtn!= 0| is.na(age_rtn)))),
    ifelse(nocc> 4, nrow(subset(cht, occ5== 2 &
        as.numeric(substr(capture, 6, nocc- 1))== 0 &
        (age_rtn!= 0| is.na(age_rtn)))), NA)) # t group
  # BT4 doesn't count smolts that return as mini-jacks
  # do this to count mini-jacks
  # x_t<- cbind(nrow(subset(cht, occ2== 2 &
  #     as.numeric(substr(capture, 3,nocc- 1))== 0)),
  #   nrow(subset(cht, occ3== 2 & as.numeric(substr(capture, 4,nocc- 1))== 0)),
  #   nrow(subset(cht, occ4== 2 & as.numeric(substr(capture, 5,nocc- 1))== 0)),
  #   ifelse(nocc> 4, nrow(subset(cht, occ5== 2 &
  #       as.numeric(substr(capture, 6,nocc- 1))== 0)), NA)) # t group

  x_0<- cbind(nrow(cht[cht[, 2]== 0 & cht[, 3]== 2,]),
    nrow(cht[as.numeric(substr(cht$capture, 2, 3))== 0 & cht[, 4]== 2,]),
    nrow(cht[as.numeric(substr(cht$capture, 2, 4))== 0 & cht[, 5]== 2,])) # t group
  d234t <- colSums(cbind(cht$d2, cht$d3, cht$d4)) # t group
  d5671t<- colSums(cbind(cht$d51, cht$d61, cht$d71)) # t group

  c0at_rtn <- sum(cht[, 'ac0_rtn'], na.rm= TRUE)  # t group
  c0ajt_rtn<- sum(cht[, 'ac0j_rtn'], na.rm= TRUE) #
  c1at_rtn <- sum(cht[, 'ac1_rtn'], na.rm= TRUE)  #
  c1ajt_rtn<- sum(cht[, 'ac1j_rtn'], na.rm= TRUE) #
  txat_rtn <- sum(cht[, 'atx_rtn'], na.rm= TRUE)  #
  txajt_rtn<- sum(cht[, 'atxj_rtn'], na.rm= TRUE) #
  t0at_rtn <- sum(cht[, 'at0_rtn'], na.rm= TRUE)  #
  t0ajt_rtn<- sum(cht[, 'at0j_rtn'], na.rm= TRUE) #
  c0at_boa <- sum(cht[, 'ac0_boa'], na.rm= TRUE)  #
  c0ajt_boa<- sum(cht[, 'ac0j_boa'], na.rm= TRUE) #
  c1at_boa <- sum(cht[, 'ac1_boa'], na.rm= TRUE)  #
  c1ajt_boa<- sum(cht[, 'ac1j_boa'], na.rm= TRUE) #
  txat_boa <- sum(cht[, 'atx_boa'], na.rm= TRUE)  #
  txajt_boa<- sum(cht[, 'atxj_boa'], na.rm= TRUE) #
  t0at_boa <- sum(cht[, 'at0_boa'], na.rm= TRUE)  #
  t0ajt_boa<- sum(cht[, 'at0j_boa'], na.rm= TRUE) #

  lgr_atx_rtn<-  sum(cht[cht[, 2]== 2, 'atx_rtn'], na.rm= TRUE)  # t group
  lgs_atx_rtn<-  sum(cht[cht[, 3]== 2, 'atx_rtn'], na.rm= TRUE)  #
  lmn_atx_rtn<-  sum(cht[cht[, 4]== 2, 'atx_rtn'], na.rm= TRUE)  #
  lgr_atxj_rtn<- sum(cht[cht[, 2]== 2, 'atxj_rtn'], na.rm= TRUE) #
  lgs_atxj_rtn<- sum(cht[cht[, 3]== 2, 'atxj_rtn'], na.rm= TRUE) #
  lmn_atxj_rtn<- sum(cht[cht[, 4]== 2, 'atxj_rtn'], na.rm= TRUE) #
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


#' Estimate survivals/detection using likelihood method, with logit link, and tally adult returns
#'
#' @param ch Input file made by \code{format_dat()} function.
#' @param i Iteration number used by the bootstrap function \code{bootystrapper()}.
#' @param nocc Total detection events including the trawl.
#' @param wt Indicates whether to weight the sampling probability (bootstrap only).
#' @param phi_p_only Option to only calculate survivals and detection and not do the adult counts.
#' @param logit_link Option to use package "RMark" or "marked" for estimating survival with logit link function. The default is package "RMark."
#' @return Survivals, detection and returing adult counts
#' @examples
#' mark_calc(detect_data, i= 1, nocc= 6, wt= 'n', phi_p_only= 'y')
#'
mark_calc<- function(ch, i, nocc, wt, phi_p_only, logit_link='RMark', ...){
  # breakdown of int_t, int_r, seg_t, and seg_r ----
  # if comment out, make sure change the output in 'bootystrapper()'
  if(wt=='y') tnr<- unlist(tapply(ch$group, ch$brood, table))
  else tnr<- table(ch$group)
  # ----

  m_data<- mark_dat(ch)
  Phi_t<- list(formula= ~time, link= 'logit')
  p_t<- list(formula= ~time, link= 'logit')
  if (logit_link== 'RMark') {
    dat_proc<- RMark::process.data(data= m_data, model= 'CJS')
    dat_ddl<- RMark::make.design.data(dat_proc)
    dat_ddl$Phi$fix= ifelse(dat_ddl$Phi$time== (unique(nchar(m_data$ch))- 1), 1, NA)
    invisible(capture.output(cjs_fit<- RMark::mark(data= dat_proc, ddl= dat_ddl,
      model.parameters= list(Phi= Phi_t, p= p_t),
      output= FALSE,
      delete= TRUE)))
    phi<- summary(cjs_fit)$reals$Phi[[1]]$pim[1, -(nocc- 1)]
    p<- summary(cjs_fit)$reals$p[[1]]$pim[1, -(nocc- 1)]
  } else if (logit_link=='marked') {
    cjs_fit<- marked::crm(m_data, model.parameters= list(Phi= Phi_t, p= p_t))
    phi<- cjs_fit$results$reals$Phi[1:(nocc- 2), 3]
    p<- cjs_fit$results$reals$p[1:(nocc- 2), 3]
  } else {
    stop('Please specifiy either "RMark" or "marked" for logit link.')
  }

  if(phi_p_only=='y') {
    calc<- cbind(t(phi), t(p))
    colnames(calc)<- c(paste0('Phi', 1:(nocc- 2)), paste0('p', 2:(nocc- 1)))
    return(calc)
    stop()
  }
  # ----

  sim_mary<- marray(ch, nocc)
  # output params
  # params from crt group ----
  R1<- sim_mary[1, 1]
  c0a_rtn  <- sum(ch$ac0_rtn, na.rm= TRUE)
  c0aj_rtn <- sum(ch$ac0j_rtn, na.rm= TRUE)
  c0a_boa  <- sum(ch$ac0_boa, na.rm= TRUE)
  c0aj_boa <- sum(ch$ac0j_boa, na.rm= TRUE)
  d5670<- colSums (cbind(ch$d50, ch$d60, ch$d70))
  m12<- sim_mary[1, 2]
  m13<- sim_mary[1, 3]
  m14<- sim_mary[1, 4]
  # ----
  # params from t group ----
  R1t<- sum(ch$group== 'T') # t group

  cht<- subset(ch, group== 'T')
  # BT4 doesn't count smolts that return as mini-jacks
  # there's an option to count mini-jacks in surv_calc (but not here)
  x_t<- cbind(nrow(subset(cht, occ2== 2 &
      as.numeric(substr(capture, 3, nocc- 1))== 0 &
      (age_rtn!=0| is.na(age_rtn)))),
    nrow(subset(cht, occ3== 2 &
        as.numeric(substr(capture, 4,nocc- 1))== 0 &
        (age_rtn!= 0| is.na(age_rtn)))),
    nrow(subset(cht, occ4== 2 &
        as.numeric(substr(capture, 5, nocc- 1))== 0 &
        (age_rtn!= 0| is.na(age_rtn)))),
    ifelse(nocc> 4, nrow(subset(cht, occ5== 2 &
        as.numeric(substr(capture, 6, nocc- 1))== 0 &
        (age_rtn!= 0| is.na(age_rtn)))), NA)) # t group

  x_0<- cbind(nrow(cht[cht[, 2]== 0 & cht[, 3]==2,]),
    nrow(cht[as.numeric(substr(cht$capture, 2, 3))== 0 & cht[, 4]== 2,]),
    nrow(cht[as.numeric(substr(cht$capture, 2, 4))== 0 & cht[, 5]== 2,])) # t group
  d234t<- colSums(cbind(cht$d2, cht$d3, cht$d4)) # t group
  d5671t<- colSums(cbind(cht$d51, cht$d61, cht$d71)) # t group

  c0at_rtn <- sum(cht[, 'ac0_rtn'], na.rm= TRUE)  # t group
  c0ajt_rtn<- sum(cht[, 'ac0j_rtn'], na.rm= TRUE) #
  c1at_rtn <- sum(cht[, 'ac1_rtn'], na.rm= TRUE)  #
  c1ajt_rtn<- sum(cht[, 'ac1j_rtn'], na.rm= TRUE) #
  txat_rtn <- sum(cht[, 'atx_rtn'], na.rm= TRUE)  #
  txajt_rtn<- sum(cht[, 'atxj_rtn'], na.rm= TRUE) #
  t0at_rtn <- sum(cht[, 'at0_rtn'], na.rm= TRUE)  #
  t0ajt_rtn<- sum(cht[, 'at0j_rtn'], na.rm= TRUE) #
  c0at_boa <- sum(cht[, 'ac0_boa'], na.rm= TRUE)  #
  c0ajt_boa<- sum(cht[, 'ac0j_boa'], na.rm= TRUE) #
  c1at_boa <- sum(cht[, 'ac1_boa'], na.rm= TRUE)  #
  c1ajt_boa<- sum(cht[, 'ac1j_boa'], na.rm= TRUE) #
  txat_boa <- sum(cht[, 'atx_boa'], na.rm= TRUE)  #
  txajt_boa<- sum(cht[, 'atxj_boa'], na.rm= TRUE) #
  t0at_boa <- sum(cht[, 'at0_boa'], na.rm= TRUE)  #
  t0ajt_boa<- sum(cht[, 'at0j_boa'], na.rm= TRUE) #

  lgr_atx_rtn<-  sum(cht[cht[, 2]== 2, 'atx_rtn'], na.rm= TRUE)  # t group
  lgs_atx_rtn<-  sum(cht[cht[, 3]== 2, 'atx_rtn'], na.rm= TRUE)  #
  lmn_atx_rtn<-  sum(cht[cht[, 4]== 2, 'atx_rtn'], na.rm= TRUE)  #
  lgr_atxj_rtn<- sum(cht[cht[, 2]== 2, 'atxj_rtn'], na.rm= TRUE) #
  lgs_atxj_rtn<- sum(cht[cht[, 3]== 2, 'atxj_rtn'], na.rm= TRUE) #
  lmn_atxj_rtn<- sum(cht[cht[, 4]== 2, 'atxj_rtn'], na.rm= TRUE) #
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


#' Make an ".inp" file for packages RMark.
#'
#' @param ch Input file containing capture history (from \code{format_dat()} function).
#' @return RMark/marked data file
#' @examples
#' m_data<- mark_dat(ch)
#'
mark_dat<- function(ch) {
  mdat<- as.data.frame(cbind(do.call(paste0,
    as.data.frame(ch[, grep('occ', names(ch))], stringsAsFactors= FALSE)
    ), 1))
  names(mdat)<- c('ch', 'freq')
  mdat$freq<- as.numeric(mdat$freq)
  if (length(mdat[grepl('2', mdat$ch),]$freq)> 0|
      length(mdat[grepl('3', mdat$ch),]$freq)> 0) {
    mdat[grepl('2', mdat$ch)| grepl('3', mdat$ch),]$freq<- (-1)
    mdat$ch<- gsub('2', '1', mdat$ch)
    mdat$ch<- gsub('3', '1', mdat$ch)
  }
  mdat$ch<- as.character(mdat$ch)
  return(mdat) # this is the same str as the '.inp' file
}


#' Bootstrap using surv_calc and organize output
#'
#' @param d Input file made by \code{format_dat()}.
#' @param fn Function to run the bootstrap on.
#' @param iter Amount of bootstrap iterations.
#' @param wgt Indicates whether to weight the sampling probability. Default is no ("n").
#' @param wgt_init Indicates whether to calculate the original estimates using weighted probability. Default is no ("n").
#' @param phi_p_only Indicate to turn off the phi_p_only option in \code{curv_calc()}. Default is no ("n").
#' @param fpc Indicate to turn off the finite population correction option in \code{curv_calc()}. Default is yes ("y").
#' @param logit_link Indicate to use "RMark" or "marked" and estimate using logit link. The default here is none ("n").
#' @param save_name Name to save bootstrap output in CSSOUTPUT in SQL server. No results will be saved if nothing is specified.
#' @return Estimates in a data frame with original estimate as the first row and bootstrap results in the remaining rows.
#' @examples
#' # To conduct standard CSS bootstrap procedures
#' out1<- bootystrapper(detect_data, surv_calc, iter= 1000, fpc= 'n', save_name='SR HCH 2008 CATH')
#' # To conduct weighted bootstrap and produce only survival and detection estimates
#' out2<- bootystrapper(detect_data, surv_calc, iter= 100, wgt= 'y', wgt_init= 'y', phi_p_only= 'y')
#' head(out2)
#'
bootystrapper<- function(d, fn, iter, wgt='n', wgt_init='n', phi_p_only='n', fpc='y', logit_link='n', save_name='none', ...){
  start_time<- Sys.time()
  n_occ<- sum(grepl('occ', names(d)))
  # run function on original data
  if (logit_link=='n') {
    original<- fn(d, i=1, n_occ, wgt, wgt_init, phi_p_only, fpc)
  } else {
    original<- fn(d, i=1, n_occ, wgt, phi_p_only, logit_link)
  }
  #make an output matrix with NA's
  out <- matrix(data= NA, nrow= (iter+ 1), ncol= length(original))
  # first row name is original
  rownames(out)<- as.character(c("original", c(1:(iter))))
  # fill first line of output matrix with original run
  out[1,]<- original

  # builds an vector of row numbers in original
  index<- c(1:length(d$prob))
  # starts loop for number of iterations
  pb<- txtProgressBar(min= 2, max= (iter+1), char= 'x', width= 50, style= 3)
  for (i in 2:(iter+ 1)) {
    # resample index is "resampled" data rows to use
    if (wgt== 'y') sample_index<- sample(index, prob= d$prob, replace= TRUE)
    else sample_index<- sample(index, replace= TRUE)

    # build resampled data from sample.index
    sample_data<- d[sample_index,]

    # run function on resampled data (with quit loop)
    attempt<- 1
    while(is.na(out[i,]) && attempt<= 10) {
      attempt<- attempt+ 1
      try(
        if (logit_link== 'n') {
          out[i,]<- fn(sample_data, i, n_occ, wgt, wgt_init, phi_p_only, fpc)
        } else {
          out[i,]<- fn(sample_data, i, n_occ, wgt, phi_p_only, logit_link)
        }
      )
    }
    # # run function on resampled data
    # if (logit_link=='n') {
    #   out[i,] <- fn(sample_data, i, n_occ, wgt, wgt_init, phi_p_only, fpc)
    # } else {
    #   out[i,] <- fn(sample_data, i, n_occ, wgt, phi_p_only, logit_link)
    # }

    setTxtProgressBar(pb, i) # print booty progress
  } # bootstrap loop
  close(pb)
  # build output matrix and return
  out<- as.data.frame(out)

  if(phi_p_only== 'y') {
    colnames(out)<- c(paste0('phi', 1:(n_occ- 2)), paste0('p', 2:(n_occ- 1)))
  } else if(wgt== 'y') {
    out$tag_site<- d$tag_site[1]
    out$rel_site<- d$rel_site[1]
    out$coord_id<- d$coord_id[1]
    colnames(out) <- c(paste0('phi', 1:(n_occ- 2)), paste0('p', 2:(n_occ- 1)), 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'AD_R', 'AD_T', 'CW_R', 'CW_T', 'tag_site', 'rel_site', 'coord_id')
  } else if(n_occ== 8) { # n_occ= 8 and not weighted
    out$tag_site<- d$tag_site[1]
    out$rel_site<- d$rel_site[1]
    out$coord_id<- d$coord_id[1]
    colnames(out)<- c(paste0('phi', 1:(n_occ- 2)), paste0('p', 2:(n_occ- 1)), 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'R group', 'T group', 'tag_site', 'rel_site', 'coord_id')
  } else { # n_occ= 4 or 6
    out$tag_site<- d$tag_site[1]
    out$rel_site<- d$rel_site[1]
    out$coord_id<- d$coord_id[1]
    colnames(out) <- c(paste0('phi', 1:(n_occ- 2)), paste0('p', 2:(n_occ- 1)), 'R1', 'R1t', 'm12', 'm13', 'm14', 'x12t', 'x1a2t', 'x1aa2t', 'x1aaa2t', 'x102t', 'x1002t', 'x10002t', 'd2t', 'd3t', 'd4t', 'd51t', 'd61t', 'd71t', 'd50', 'd60', 'd70', 'C0adult_rtn', 'C0adultj_rtn', 'C0adult_boa', 'C0adultj_boa', 'C0adult_t_rtn', 'C0adultj_t_rtn', 'C1adult_rtn', 'C1adultj_rtn', 'Txadult_rtn', 'Txadultj_rtn', 'T0adult_rtn', 'T0adultj_rtn', 'C0adult_t_boa', 'C0adultj_t_boa', 'C1adult_boa', 'C1adultj_boa', 'Txadult_boa', 'Txadultj_boa', 'T0adult_boa', 'T0adultj_boa', 'lgradult_rtn', 'lgsadult_rtn', 'lmnadult_rtn', 'lgradultj_rtn', 'lgsadultj_rtn', 'lmnadultj_rtn', 'tag_site', 'rel_site', 'coord_id')
  }

  cat('\n')
  print(Sys.time()- start_time)

  if (save_name!= 'none') {
    # save bootystrap output to sql server CSSOUTPUT
    rand_str <- function(n) {
      do.call(paste0, replicate(7, sample(c(letters,letters,0:9), n, TRUE), FALSE))
    }
    channel <- RODBC::odbcDriverConnect("case=nochange;
      Description=CSSOUTPUT;
      DRIVER=SQL Server;
      SERVER=PITTAG_2016;
      UID=sa;
      PWD=frznool;
      WSID=CUTTHROAT;
      DATABASE=CSSOUTPUT;
      Network=DBMSSOCN")
    RODBC::sqlSave(channel, data.frame(out), tablename=
        paste0('C_T_', save_name, '_bootylator_', rand_str(1),
          format(Sys.time(), '%m%d%Y')))
    RODBC::odbcCloseAll()
  }

  return(out)
}


#' Construct m-ij array
#'
#' @param CH Input file made by \code{format_dat()} function.
#' @param n_occ Total detection events including the trawl.
#' @return m-ij array
#' @examples
#' marray(detect_data, n_occ= 8)
#'
marray<- function(CH, n_occ) {
  n_ind<- dim(CH)[1]
  m_array<- matrix(data= 0, nrow= n_occ, ncol= n_occ+ 1)
  for(t in 1:n_occ){
    m_array[t, 1]<- sum(CH[CH[, t]== 1, t])
  } # R(i)

  for(t in 1:(n_occ- 1)) {
    m_array[t, (t+ 1)]<- nrow(CH[CH[, t]== 1& CH[,(t+ 1)]!= 0,])
    if (t> n_occ- 2) next
    m_array[t,(t+ 2)]<- nrow(CH[CH[, t]== 1& CH[,(t+ 1)]==0& CH[,(t+ 2)]!= 0,])
    if (t> n_occ- 3) next
    for(u in (t+ 3):n_occ) {
      m_array[t, u]<- nrow(CH[CH[, t]== 1&
          rowSums(CH[, (t+ 1):(u- 1)])== 0& CH[, u]!= 0,])
    }
  } # m(ij)

  for(t in 1:n_occ){
    m_array[t, n_occ+ 1]<- sum(m_array[t, 2:n_occ])
  } # r(i)
  out<- m_array[-n_occ,]
  dimnames(out)<- list(1:(n_occ- 1), c('R(i)', 'j=2', 3:n_occ, 'r(i)'))
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
#' marray(wgt_data, n_occ= 8)
#'
marray_wtd<- function(CH, n_occ){
  n_ind<- dim(CH)[1]
  wtd<- CH[,'prob']* n_ind
  m_array<- matrix(data= 0, nrow= n_occ, ncol= n_occ+ 1)
  for(t in 1:n_occ){
    m_array[t, 1]<- sum(wtd[CH[, t]== 1])
  } # R(i)

  for(t in 1:(n_occ- 1)){
    m_array[t, (t+ 1)]<- sum(wtd[CH[, t]== 1& CH[, (t+ 1)]!= 0])
    if (t> n_occ- 2) next
    m_array[t, (t+ 2)]<- sum(wtd[CH[, t]== 1&
        CH[, (t+ 1)]== 0& CH[, (t+ 2)]!= 0])
    if (t> n_occ- 3) next
    for(u in (t+ 3):n_occ){
      m_array[t, u]<- sum(wtd[CH[, t]== 1&
          rowSums(CH[, (t+ 1):(u- 1)])== 0& CH[, u]!= 0])
    }
  } # m(ij)

  for(t in 1:n_occ){
    m_array[t, n_occ+ 1]<- sum(m_array[t, 2:n_occ])
  } # r(i)
  out<- m_array[-n_occ,]
  dimnames(out)<- list(1:(n_occ- 1), c('R(i)', 'j=2', 3:n_occ, 'r(i)'))
  return(out)
}

































