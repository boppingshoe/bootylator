#Function to take BT4 output and return summary data for all annual report/web tables.
#Last modified 5/19/2014 by JM
#makefile = 'Y'
#reaches =  6
#target  =  'SR HCH 2010 RAPH'
#migr_yr = 2010

# doCalcs_v3_Bootstrap_in_R <- function(makefile, target, location, migr_yr, css_group,rel_site,tag_site,coord_id){
# library(RODBC)
# channel <- odbcDriverConnect("case=nochange;
#                               Description=CSSOUTPUT;
#                               DRIVER=SQL Server;
#                               SERVER=PITTAG_SQL6;
#                               UID=sa;
#                               PWD=frznool;
#                               WSID=CUTTHROAT;
#                               DATABASE=CSSOUTPUT;
#                               Network=DBMSSOCN")
#
# #target<-'HCH 2011 WSPH'
# tables<-sqlQuery(channel, paste("select * from sys.tables where name like " , "'%", target, "%'", " order by name asc", sep=""))
# tables$name
#
# crtnm<-factor(tables$name[1])
# #rnm<-factor(tables$name[2])
# #tnm<-factor(tables$name[3])
#
# crtfile<-sqlFetch(channel, crtnm, colnames = FALSE, rownames = TRUE)
# #tfile<-sqlFetch(channel, tnm, colnames = FALSE, rownames = TRUE)
# #rfile<-sqlFetch(channel, rnm, colnames = FALSE, rownames = TRUE)
#
#
# crt<-crtfile
# #t<-tfile
# #r<-rfile

#
#------------------------------------------------------------------------------
#Code written by Jack Tuomik
#May 2012
#this function rounds things in similar way to Excel and I think Foxpro
#by default, R has it's own rounding method that is different from the other two
doCalcs_uc <- function(crt, ...){

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
get.CIs <- function(data){
  n<- length(data)- 1
  result <- data.frame(0,0,0,0,0,0,0,0,0,0)
  names(result) <- c("initial", "np_90cill", "np_90ciul", "boots_avg",
    "boots_std", "cv", "p_90cill", "p_90ciul",
    "np_95cill", "np_95ciul")
  result$initial    <- data[1]
  result$np_90cill  <- quantile(data, 0.05, na.rm= TRUE)
  result$np_90ciul  <- quantile(data, 0.95, na.rm= TRUE)
  result$boots_avg  <- mean(data[-1], na.rm= TRUE)
  result$boots_std  <- sqrt(var(data[-1], na.rm= TRUE)*(n-1)/n) # pop'n sd
  result$cv         <- result$boots_std/result$boots_avg
  result$p_90cill   <- result$boots_avg - 1.645*result$boots_std
  result$p_90ciul   <- result$boots_avg + 1.645*result$boots_std
  result$np_95cill  <- quantile(data, 0.025, na.rm= TRUE)
  result$np_95ciul  <- quantile(data, 0.975, na.rm= TRUE)

  return(result)
}
fill.CIs <- function(data){

            result <- data.frame(data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1])
            names(result) <- c("initial", "np_90cill", "np_90ciul", "boots_avg", "boots_std",
                                "cv", "p_90cill", "p_90ciul", "np_95cill", "np_95ciul")
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
popula_cjs  <- r1 * s1_cjs

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

SR<- apply(crt[grep('phi', names(crt))][-1], 1, prod)
R1<- get.CIs(r1)
p2<- get.CIs(p2)
if('p3' %in% names(crt)){p3<- get.CIs(p3)} else {p3 <- get.CIs(rep(0, 1001))}
if('p4' %in% names(crt)){p4<- get.CIs(p4)} else {p4 <- get.CIs(rep(0, 1001))}
if('p5' %in% names(crt)){p5<- get.CIs(p5)} else {p5 <- get.CIs(rep(0, 1001))}
if('p6' %in% names(crt)){p6<- get.CIs(p6)} else {p6 <- get.CIs(rep(0, 1001))}

popula_cjs<- get.CIs(r1 * s1_cjs)
adults_b<- get.CIs(totaladult_b)
adults_bj<- get.CIs(totaladult_bj)
adults_m<- get.CIs(totaladult_m)
adults_mj<- get.CIs(totaladult_mj)

releaseSAR_b<- get.CIs(sar_rel_b) * 100
releaseSAR_bj<- get.CIs(sar_rel_bj) * 100
overallSAR_b<- get.CIs(sar_tws_cr_b) * 100
overallSAR_bj<- get.CIs(sar_tws_cr_bj) * 100
releaseSAR_m<- get.CIs(sar_rel_m) * 100
overallSAR_m<- get.CIs(sar_tws_cr_m) * 100
overallSAR_mj<- get.CIs(sar_tws_cr_mj) * 100

#McNpop       <- get.CIs(popula_cjs)

#wrangle format and output------------------------------------------------------
#choose parameters:
parm <- c(
          #general stuff for evaluation; has various reach survival versions
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
#
ansDF <- get(parm[1])
for(i in parm[-1]){ansDF <- rbind(ansDF, get(i))}

ans <- cbind(parm, ansDF)
#makefile = 'Y'
#write file
# if(makefile == 'Y'){
# nm <- paste("doCalcs_v3_bstrp_R ", target, format(Sys.time(), '%Y-%m-%d %H%M%S')
#                   ," ",".csv", sep = "")
# write.table(ans,
#             paste(nm, sep = "\\"),
#             col.names = T, row.names = F, sep = ',',
#             quote = F)}
# #also send to clipboard
# write.table(ans,
#             'clipboard',
#             col.names = T, row.names = F, sep = ',',
#             quote = F)
#
# # save to sql server CSS report
# channel3 <- odbcDriverConnect("case=nochange;
#                               Description=CSSREPORT;
#                               DRIVER=SQL Server;
#                               SERVER=PITTAG_SQL6;
#                               UID=sa;
#                               PWD=frznool;
#                               WSID=CUTTHROAT;
#                               DATABASE=CSSREPORT;
#                               Network=DBMSSOCN")
#add things to output that will identify data in the sql table
#for practice
#migr_yr<-2011
#css_group<-"WSPH"
#rel_site<-"WSPH"
#tag_site<-"WSPH"
#coord_id<-"JAR"

# ans$input_sqlfile<-target
# ans$migr_year<-migr_yr
# ans$css_group<-css_group
# ans$rel_site<-rel_site
# ans$tag_site<-tag_site
# ans$coord_id<-coord_id
# ans$flag<-as.character('')
##IS THIS DATA ALREADY IN THE TABLE?
# inthesqltable<-sqlQuery(channel3, paste("select css_group, migr_year, parm, thecount = count(initial) from BOOTSTRAP_RESULTS_COL where css_group = " , "'", css_group, "'", " and migr_year = " , "'", migr_yr, "'", " and parm = 'overallSAR_b' group by css_group, migr_year, parm", sep=""))

#new version updates if existing data otherwise inserts if not
# if(length(inthesqltable$thecount) != '0'){
# {sqlUpdate(channel3, data.frame(ans),tablename="BOOTSTRAP_RESULTS_COL", index = c('css_group', 'migr_year', 'parm'), verbose = FALSE, test = FALSE, nastring = NULL,fast = TRUE)}}else{
# sqlSave(channel3,data.frame(ans),tablename="BOOTSTRAP_RESULTS_COL",safer=TRUE,append=TRUE)}
#old version just inserts data
#sqlSave(channel3,data.frame(ans),tablename="BOOTSTRAP_RESULTS_COL",safer=TRUE,append=TRUE)
#finally print in R
# print(format(ans, scientific = F))
rownames(ans)<- parm
return(ans)
}
#------------------------------------------------------------------------------
