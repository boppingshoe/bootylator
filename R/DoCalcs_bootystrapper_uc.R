

#' Take bootstrap output and return summary data for all annual report/web table for Upper and Mid Columbia River fish
#'
#' @param crt Bootstrap output as input file here.
#' @param target The name of the original input file (eg. SR HCH 2015 MCCA).
#' @param css_group CSS group from the CSSGroups_LookupTable.
#' @param makefile Save bootstrap output in CSSOUTPUT in SQL server, append parameter output in CSSREPORT in SQL server, and make parameter output in csv file in working directory. Default is 'y'.
#' @return Survivals,detection, adult counts, and SARs...
#' @examples
#' ans<- doCalcs_uc(crt, target= 'CR_USK_2015_OKAN_MCJ', css_group= 'OKSR', makefile= 'y')
#' format(ans, scientific= FALSE)


doCalcs_uc <- function(crt, target, css_group, makefile='y', ...){

roundT <- function(x, ...){
  if(is.null(names(list(...)))){
  ans <- floor(x + 0.5)   }else{
  digits <- list(...)[[1]]
  ans <- floor(x*(10^digits) + 0.5)/10^digits
  }
  return(ans)}

#roundT(2.5); roundT(2.25, digits=1)
##vs
#round(2.5); round(2.25, digits=1)

#------------------------------------------------------------------------------
#90%,95% non parametric CI function similar to the old bootstrap program [[modified order, 7-24-2013 JET]]

exactci <- function(x, n, conflev){
  alpha <- (1- conflev)
  if (x == 0) {
    ll <- 0
    ul <- 1 - (alpha/ 2)^ (1/ n) # qbeta(1- alpha/ 2, x+ 1, n- x)
  }
  else if (x == n) {
    ll <- (alpha/ 2)^ (1/ n)
    ul <- 1
  }
  else {
    ll <- 1/ (1+ (n- x+ 1)/ (x* qf(alpha/ 2, 2* x, 2* (n- x+ 1))) )
    ul <- 1/ (1+ (n- x)/
        ((x+ 1) * qf(1- alpha/ 2, 2* (x+ 1), 2* (n- x))) )
  }
  c(ll, ul)
}
# Computes the Clopper/Pearon exact ci
# for a binomial success probability
# for x successes out of n trials with
# confidence coefficient conflev
# from http://users.stat.ufl.edu/~aa/cda/R/one-sample/R1/index.html
# date of visit: 1/14/2019

get.CIs <- function(data, exci='n', exfn, x, n, conflev, ...){
  n_dat<- length(data)- 1
  result <- data.frame(0,0,0,0,0,0,0,0,0,0, NA,NA)
  names(result) <- c("initial", "np_90cill", "np_90ciul", "boots_avg",
    "boots_std", "cv", "p_90cill", "p_90ciul", "np_95cill", "np_95ciul",
    "ex_90cill", "ex_90ciul")
  result$initial    <- data[1]
  result$np_90cill  <- quantile(data, 0.05, na.rm= TRUE)
  result$np_90ciul  <- quantile(data, 0.95, na.rm= TRUE)
  result$boots_avg  <- mean(data[-1], na.rm= TRUE)
  result$boots_std  <- sqrt(var(data[-1], na.rm= TRUE)*(n_dat-1)/n_dat) # pop'n sd
  result$cv         <- result$boots_std/result$boots_avg
  result$p_90cill   <- result$boots_avg - 1.645*result$boots_std
  result$p_90ciul   <- result$boots_avg + 1.645*result$boots_std
  result$np_95cill  <- quantile(data, 0.025, na.rm= TRUE)
  result$np_95ciul  <- quantile(data, 0.975, na.rm= TRUE)
  if (exci== 'y') {
    result$ex_90cill<- exfn(x, n, conflev= 0.9)[1]* 100
    result$ex_90ciul<- exfn(x, n, conflev= 0.9)[2]* 100
  }

  return(result)
}
fill.CIs <- function(data){
  result <- data.frame(data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1], NA,NA)
  names(result) <- c("initial", "np_90cill", "np_90ciul", "boots_avg", "boots_std", "cv", "p_90cill", "p_90ciul", "np_95cill", "np_95ciul", "ex_90cill", "ex_90ciul")
  return(result)
}


#load up the needed pieces------------------------------------------------------
nocc<- sum(grepl('phi', names(crt)))
s1_cjs    <- crt$phi1
s2_cjs    <- crt$phi2
if('phi3' %in% names(crt)) {s3_cjs<- crt$phi3} else {s3_cjs<- (rep(0, 1001))}
if('phi4' %in% names(crt)) {s4_cjs<- crt$phi4} else {s4_cjs<- (rep(0, 1001))}
if('phi5' %in% names(crt)) {s5_cjs<- crt$phi5} else {s5_cjs<- (rep(0, 1001))}
# s4_cjs    <- crt$phi4
p2        <- crt$p2
p3        <- crt$p3
if('p4' %in% names(crt)){p4<- crt$p4} else {p4<- (rep(0, 1001))}
if('p5' %in% names(crt)){p5<- crt$p5} else {p5<- (rep(0, 1001))}
if('p6' %in% names(crt)){p6<- crt$p6} else {p6<- (rep(0, 1001))}
r1        <- crt$R1

totaladult_b<- round((crt$C0adult_t_boa + crt$C1adult_boa + crt$Txadult_boa),0)

totaladult_bj<- round((crt$C0adultj_t_boa + crt$C1adultj_boa + crt$Txadultj_boa),0)

totaladult_m<- round((crt$C0adult_t_rtn + crt$C1adult_rtn + crt$Txadult_rtn),0)

totaladult_mj<- round((crt$C0adultj_t_rtn + crt$C1adultj_rtn + crt$Txadultj_rtn),0)
# ----



#Smolt populations--------------------------------------------------------------
#smolts
# popula_cjs<- r1* s1_cjs

#overall SARs by location and with and without jacks
sar_tws_cr_b <-  totaladult_b/(s1_cjs * r1)
sar_tws_cr_bj <- totaladult_bj/(s1_cjs * r1)
sar_rel_b <-  totaladult_b/r1
sar_rel_bj <-  totaladult_bj/r1
sar_tws_cr_m <-  totaladult_m/(s1_cjs * r1)
sar_tws_cr_mj <- totaladult_mj/(s1_cjs * r1)
sar_rel_m <-  totaladult_m/r1
# ----

#load up answers----------------------------------------------------------------

S1 <- get.CIs(s1_cjs)
if('phi2' %in% names(crt)) {S2<- get.CIs(s2_cjs)} else {S2 <- get.CIs(rep(0, 1001))}
if('phi3' %in% names(crt)) {S3<- get.CIs(crt$phi3)} else {S3 <- get.CIs(rep(0, 1001))}
if('phi4' %in% names(crt)) {S4<- get.CIs(crt$phi4)} else {S4 <- get.CIs(rep(0, 1001))}
if('phi5' %in% names(crt)) {S5<- get.CIs(crt$phi5)} else {S5 <- get.CIs(rep(0, 1001))}

SR<- get.CIs(apply(crt[grep('phi', names(crt))][-1], 1, prod))
R1<- get.CIs(r1)
p2<- get.CIs(p2)
if('p3' %in% names(crt)){p3<- get.CIs(p3)} else {p3 <- get.CIs(rep(0, 1001))}
if('p4' %in% names(crt)){p4<- get.CIs(p4)} else {p4 <- get.CIs(rep(0, 1001))}
if('p5' %in% names(crt)){p5<- get.CIs(p5)} else {p5 <- get.CIs(rep(0, 1001))}
if('p6' %in% names(crt)){p6<- get.CIs(p6)} else {p6 <- get.CIs(rep(0, 1001))}

popula_cjs<- get.CIs(r1* s1_cjs)
adults_b<- get.CIs(totaladult_b)
adults_bj<- get.CIs(totaladult_bj)
adults_m<- get.CIs(totaladult_m)
adults_mj<- get.CIs(totaladult_mj)

releaseSAR_b<- get.CIs(sar_rel_b* 100, exci='y',
  exactci, x= totaladult_b[1], n= r1[1])
releaseSAR_bj<- get.CIs(sar_rel_bj* 100, exci='y',
  exactci, x= totaladult_bj[1], n= r1[1])
overallSAR_b<- get.CIs(sar_tws_cr_b* 100, exci='y',
  exactci, x= totaladult_b[1], n= (r1* s1_cjs)[1])
overallSAR_bj<- get.CIs(sar_tws_cr_bj* 100, exci='y',
  exactci, x= totaladult_bj[1], n= (r1* s1_cjs)[1])

releaseSAR_m<- get.CIs(sar_rel_m* 100, exci='y',
  exactci, x= totaladult_m[1], n= r1[1])
overallSAR_m<- get.CIs(sar_tws_cr_m* 100, exci='y',
  exactci, x= totaladult_m[1], n= (r1* s1_cjs)[1])
overallSAR_mj<- get.CIs(sar_tws_cr_mj* 100, exci='y',
  exactci, x= totaladult_mj[1], n= (r1* s1_cjs)[1])
# McNpop       <- get.CIs(r1* s1_cjs)

#wrangle format and output------------------------------------------------------
#choose parameters:
parm <- c(
          # general stuff for evaluation; has various reach survival versions
          # 'S1', 'S2', 'S3', 'S4', 'SR', 'R1'
          paste0('S', 1:nocc)
          , 'SR', 'R1'
          , paste0('p', 2:(nocc+1))
          , 'popula_cjs'
          #,'p2', 'p3', 'p4','p5', 'popula_cjs'

          ,'adults_b','adults_bj'
          ,'adults_m','adults_mj'

          #overall SARs
          ,'releaseSAR_b','releaseSAR_bj','overallSAR_b', 'overallSAR_bj'
          ,'releaseSAR_m','overallSAR_m', 'overallSAR_mj'
  )
#
ansDF <- get(parm[1])
for(i in parm[-1]){ansDF <- rbind(ansDF, get(i))}
ans <- cbind(parm, ansDF)


# option to save output to working directory and sql server
if (makefile == 'y' | makefile == 'Y') {
  # save bootystrap output to sql server
  rand_str <- function(n) {
    do.call(paste0, replicate(7, sample(c(letters, letters, 0:9), n, TRUE), FALSE))
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
  RODBC::sqlSave(channel, data.frame(crt), tablename=
      paste0('C_T_', target, '_bootylator_', rand_str(1), format(Sys.time(), '%m%d%Y')))
  RODBC::odbcCloseAll()

  # write doCalcs output in csv file to working directory
  nm <- paste("doCalcs_uc ", target, format(Sys.time(), '%Y-%m-%d %H%M%S'),
    " ", ".csv", sep = "")
  # ans$css_group<-css_group
  write.table(ans, paste(nm, sep = "\\"),
    col.names = T, row.names = F, sep = ',', quote = F)

  # also send to clipboard
  # write.table(ans,
  #             'clipboard',
  #             col.names = T, row.names = F, sep = ',',
  #             quote = F)

  # save to sql server CSS report
  channel2 <- RODBC::odbcDriverConnect("case=nochange;
                                Description=CSSREPORT;
                                DRIVER=SQL Server;
                                SERVER=PITTAG_2016;
                                UID=sa;
                                PWD=frznool;
                                WSID=CUTTHROAT;
                                DATABASE=CSSREPORT;
                                Network=DBMSSOCN")
  # get rid of things that sql doesn't like
  ans[ans=="Inf"]<- NA
  ans[ans=="NaN"]<- NA
  # add things to output that will identify data in the sql table
  migr_yr<- as.numeric(regmatches(target, regexpr("[0-9]...", target)))

  ans$input_sqlfile<- target
  ans$migr_year<- migr_yr
  ans$css_group<- css_group # needs manual input
  ans$rel_site<- crt$rel_site[1]
  ans$tag_site<- crt$tag_site[1]
  ans$coord_id<- crt$coord_id[1]
  ans$flag<- as.character('')
  # IS THIS DATA ALREADY IN THE TABLE?
  inthesqltable<- RODBC::sqlQuery(channel2,
    paste("select css_group, migr_year, parm, thecount = count(initial)
      from BOOTSTRAP_RESULTS
      where css_group = ", "'", css_group, "'",
      " and migr_year = ", "'", migr_yr, "'",
      " and parm = 'overallSAR'
      group by css_group, migr_year, parm", sep="")
    )

  # EITHER UPDATE OR SAVE DEPENDING ON ANSWER
  if(length(inthesqltable$thecount) != '0') {
    RODBC::sqlUpdate(channel2, data.frame(ans), tablename="BOOTSTRAP_RESULTS",
      index = c('css_group', 'migr_year', 'parm'),
      verbose = FALSE, test = FALSE, nastring = NULL, fast = TRUE)
  } else {
    RODBC::sqlSave(channel2, data.frame(ans), tablename="BOOTSTRAP_RESULTS",
      safer=TRUE, append=TRUE)
  }

  RODBC::odbcCloseAll()
}

# finally print in R
# print(format(ans, scientific = F))
rownames(ans)<- parm
return(ans)
}
#------------------------------------------------------------------------------
