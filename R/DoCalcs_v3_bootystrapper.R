
#' Take bootstrap output and return summary data for all annual report/web table for Snake River fish
#'
#' @param crt Bootstrap output as input file here.
#' @param reaches Number of reach expansion. Default is 6.
#' @param species "CH" if the species was Chinook. "ELSE" for all others.
#' @param target The name of the original input file (eg. SR HCH 2015 MCCA).
#' @param css_group CSS group from the CSSGroups_LookupTable.
#' @param makefile Save bootstrap output in CSSOUTPUT in SQL server, append parameter output in CSSREPORT in SQL server, and make parameter output in csv file in working directory. Default is 'y'.
#' @return Survivals,detection, adult counts, and SARs...
#' @examples
#' ans<- doCalcs_v3(crt, reaches= 6, species= 'CH', target= 'SR HCH 2015 MCCA', css_group= 'MCCA', makefile= 'y')
#' format(ans, scientific= FALSE)


doCalcs_v3 <- function(crt, reaches=6, species, target, css_group, makefile='y', ...) {

#------------------------------------------------------------------------------
# Modified from doCalcs_v3 (by Jack Tuomikoski)

# R rounds to the nearest even number, roundT rounds up

roundT <- function(x, ...){
  if(is.null(names(list(...)))) {
    ans <- floor(x + 0.5)   }
  else {
    digits <- list(...)[[1]]
    ans <- floor(x*(10^digits) + 0.5)/10^digits}
  return(ans)
}

# roundT(2.5); roundT(2.25, digits=1)
# vs.
# round(2.5); round(2.25, digits=1)

# 90%,95% non parametric CI function similar to the old bootstrap program
# [[modified order, 7-24-2013 JET]]

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
    result$ex_90cill<- exfn(x, n, conflev=0.9)[1]
    result$ex_90ciul<- exfn(x, n, conflev=0.9)[2]
  }

  return(result)
}
fill.CIs <- function(data){
  result <- data.frame(data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1],data[1], NA,NA)
  names(result) <- c("initial", "np_90cill", "np_90ciul", "boots_avg", "boots_std", "cv", "p_90cill", "p_90ciul", "np_95cill", "np_95ciul", "ex_90cill", "ex_90ciul")
  return(result)
}

#load up the needed pieces------------------------------------------------------
s1_cjs    <- crt$phi1

# This is not correct if s2 drifts above one, may need an s2mod field here
s2_cjs    <- crt$phi2
s3_cjs    <- crt$phi3
s4_cjs    <- crt$phi4
s5_cjs    <- crt$phi5
s6_cjs    <- crt$phi6

#Modify s3*s3
#  for 2010 and a few other cases where s2*s3 is over 1
#  Affects :
#    c0_cjs, c1_cjs, c1_cjsNEW, t0_cjs, t1_cjs
#    C0_SAR, C1_SAR, C1_SARNEW, T0_SAR, Tx_SAR, TIR, D
s2s3_cjsMod <- ifelse(s2_cjs * s3_cjs > 1, 1, s2_cjs * s3_cjs)

p2        <- crt$p2
p3        <- crt$p3
p4        <- crt$p4
p5        <- crt$p5
p6        <- crt$p6
p7        <- crt$p7


x12.       <- crt$x12t # t group
x102.      <- crt$x102t #
x1002.     <- crt$x1002t #
x10002.    <- crt$x10002t #
x1a2.      <- crt$x1a2t #
x1aa2.     <- crt$x1aa2t #
x1aaa2.    <- crt$x1aaa2t #

delta5_0  <- crt$d50
delta6_0  <- crt$d60
delta7_0  <- crt$d70
delta5_1.  <- crt$d51t # t gruop
delta6_1.  <- crt$d61t #
delta7_1.  <- crt$d71t #

delta2.   <- crt$d2t # t group
delta3.   <- crt$d3t #
delta4.   <- crt$d4t #

R1        <- crt$R1
R1.       <- crt$R1t # t group

m12       <- crt$m12
m13       <- crt$m13
m14       <- crt$m14

# m12.      <- crt$m12t # t group
# m13.      <- crt$m13t #
# m14.      <- crt$m14t #

# 1/15/2019 added ifelse statements for chinooka and "others" counting methods
if (species== 'CH') { # chinooka don't count jack
  c0adults  <- round(crt$C0adult_rtn, 0)
  c0adults. <- round(crt$C0adult_t_rtn, 0)
  c1adults. <- round(crt$C1adult_rtn, 0)
  t0adults. <- round(crt$T0adult_rtn, 0)
  t1adults. <- round(crt$Txadult_rtn, 0)
} else {
  c0adults  <- round(crt$C0adultj_rtn, 0)
  c0adults. <- round(crt$C0adultj_t_rtn, 0)
  c1adults. <- round(crt$C1adultj_rtn, 0)
  t0adults. <- round(crt$T0adultj_rtn, 0)
  t1adults. <- round(crt$Txadultj_rtn, 0)
}

### May not need to explicitly label these since the output would be the overall SAR by location and with and without jacks assigned in next section
c0adults_gj.  <- round(crt$C0adultj_t_rtn, 0)
c1adults_gj. <- round(crt$C1adultj_rtn, 0)
t0adults_gj. <- round(crt$T0adultj_rtn, 0)
t1adults_gj. <- round(crt$Txadultj_rtn, 0)

if (species == 'CH') { # chinooka don't count jack
  c0adults_b. <- round(crt$C0adult_t_boa, 0)
  c1adults_b. <- round(crt$C1adult_boa, 0)
  t0adults_b. <- round(crt$T0adult_boa, 0)
  t1adults_b. <- round(crt$Txadult_boa, 0)
} else {
  c0adults_b. <- round(crt$C0adultj_t_boa, 0)
  c1adults_b. <- round(crt$C1adultj_boa, 0)
  t0adults_b. <- round(crt$T0adultj_boa, 0)
  t1adults_b. <- round(crt$Txadultj_boa, 0)
}

c0adults_bj. <- round(crt$C0adultj_t_boa, 0)
c1adults_bj. <- round(crt$C1adultj_boa, 0)
t0adults_bj. <- round(crt$T0adultj_boa, 0)
t1adults_bj. <- round(crt$Txadult_boa, 0)

# deal with steelhead and sockeye not having 'jacks'
if (species== 'CH') {
  totaladult.<- round(crt$C0adult_t_rtn+ crt$C1adult_rtn+ crt$Txadult_rtn)
  totaladult_gj.<- round(crt$C0adultj_t_rtn+ crt$C1adultj_rtn+ crt$Txadultj_rtn)
  totaladult_b.<- round(crt$C0adult_t_boa+ crt$C1adult_boa+ crt$Txadult_boa)
  totaladult_bj.<- round(crt$C0adultj_t_boa+ crt$C1adultj_boa+ crt$Txadultj_boa)
} else {
  totaladult.<- totaladult_gj.<-
    round(crt$C0adultj_t_rtn+ crt$C1adultj_rtn+ crt$Txadultj_rtn)
  # totaladult_gj.<- round(crt$C0adultj_t_rtn+ crt$C1adultj_rtn+ crt$Txadultj_rtn)
  totaladult_b.<- totaladult_bj.<-
    round(crt$C0adultj_t_boa+ crt$C1adultj_boa+ crt$Txadultj_boa)
  # totaladult_bj.<- round(crt$C0adultj_t_boa+ crt$C1adultj_boa+ crt$Txadultj_boa)
}

# make sure for Tx SARs to LGS and LMM use the lgsadults2_gra and lmnadults2_gra for the Xa2 and xaa2 seen or unseen above post-2005.
if (species== 'CH') {
  lgradults. <- round(crt$lgradult_rtn, 0)
  lgsadults. <- round(crt$lgsadult_rtn, 0)
  lmnadults. <- round(crt$lmnadult_rtn, 0)
} else {
  lgradults. <- round(crt$lgradultj_rtn, 0)
  lgsadults. <- round(crt$lgsadultj_rtn, 0)
  lmnadults. <- round(crt$lmnadultj_rtn, 0)
}


#reaches <- 6 for an example plug in
#juvenile survival--------------------------------------------------------------
if(reaches == 6){
  vc_mcn  <- s2_cjs * s3_cjs * s4_cjs
  vc_jda  <- s2_cjs * s3_cjs * s4_cjs * s5_cjs
  vc_cjs  <- s2_cjs * s3_cjs * s4_cjs * s5_cjs * s6_cjs
}else if(reaches == 5){
  vc_mcn  <-  s2_cjs * s3_cjs * s4_cjs
  vc_jda  <-  s2_cjs * s3_cjs * s4_cjs
  vc_cjs  <- (s2_cjs * s3_cjs * s4_cjs * s5_cjs)^(285.9 / 216.4)
}else if(reaches == 4){
  vc_mcn  <-  s2_cjs * s3_cjs * s4_cjs
  vc_jda  <- (s2_cjs * s3_cjs * s4_cjs)^(216.4 / 139.8)
  vc_cjs  <- (s2_cjs * s3_cjs * s4_cjs)^(285.9 / 139.8)
}else if(reaches == 3){
  vc_mcn  <- (s2_cjs * s3_cjs)^(139.8 / 65.86)
  vc_jda  <- (s2_cjs * s3_cjs)^(216.4 / 65.86)
  vc_cjs  <- (s2_cjs * s3_cjs)^(285.9 / 65.86)
}

#do several types at once
vc_cjs_0expan   <- (s2_cjs * s3_cjs * s4_cjs * s5_cjs * s6_cjs)
vc_cjs_1expan   <- (s2_cjs * s3_cjs * s4_cjs * s5_cjs)^(285.9/216.4)
vc_cjs_2expan   <- (s2_cjs * s3_cjs * s4_cjs)^(285.9/139.8)
vc_cjs_3expan   <- (s2_cjs * s3_cjs)^(285.9/65.86)
# vc_jda_01expan  <- (s2_cjs * s3_cjs * s4_cjs * s5_cjs)
# vc_jda_2expan   <- (s2_cjs * s3_cjs * s4_cjs ) ^ (216.4/139.8)
# vc_jda_3expan   <- (s2_cjs * s3_cjs ) ^ (216.4/65.86)
# vc_mcn_012expan <- (s2_cjs * s3_cjs * s4_cjs)
# vc_mcn_3expan   <- (s2_cjs * s3_cjs) ^ (139.8/65.86)

#Barge survival for smolts - mechanically this is a harmonic mean---------------
#No rounding is added to match old bootstrap code
vt_cjs_NEW <- 0.98 * (x12. + x1a2. + x1aa2.) /
                     (x12. +
                      x1a2. / s2_cjs +
                      x1aa2./(s2_cjs * s3_cjs))

vt_cjs_NEWMod <- 0.98 * (x12. + x1a2. + x1aa2.) /
                     (x12. +
                      x1a2. / s2_cjs +
                      x1aa2./(s2s3_cjsMod))

#Smolt populations--------------------------------------------------------------
#LGR equivalent transported smolts
tlgr_cjs     <- x12.
tlgs_cjs     <- x102.   / s2_cjs
tlmn_cjs     <- x1002.  / (s2_cjs * s3_cjs)
tmcn_cjs     <- x10002. / vc_mcn
tlgs2_cjs    <- x1a2.   / s2_cjs
tlmn2_cjs    <- x1aa2.  / (s2_cjs * s3_cjs)
tmcn2_cjs    <- x1aaa2. / vc_mcn

tlmn_cjsMod     <- x1002.  / (s2s3_cjsMod)
tlmn2_cjsMod    <- x1aa2.  / (s2s3_cjsMod)

#in-river removals summary
#<<!!ROUNDED AS PER OLD BOOTSTRAP PROGRAM!!>>
delta0  <- roundT(delta5_0 / vc_mcn + delta6_0 / vc_jda + delta7_0 / vc_cjs)
delta1  <- roundT(delta5_1. / vc_mcn + delta6_1. / vc_jda + delta7_1. / vc_cjs)

#popula_cjs can be calculated 2 ways as:
# 1) R1*s1_cjs
#       Or
# 2) m2 / p2
#
#Since, m2/p2     = m2 / [m2/(m2+(z2*(R2/lr2)))]
#	                = m2+(z2*(R2/lr2))
#and
#       R1*s1_cjs = [(m2+(z2*(R2/lr2)))/R1] * R1
#		              = m2+(z2*(R2/lr2))
#The old bootstrap program used version 1:
#Line 2179: CRT_popula_cjs = CRT_S1_CJS * R1
#Line 2180: popula_cjs     = CRT_popula_cjs
#Line 2180: popula_cjs     = CRT_popula_cjs
#So, sticking with version 1 here for consistancy

#strangley, this smolt population is not rounded in the old bootstrap programs
#doing similarly here to match
popula_cjs  <- R1 * s1_cjs
popula_cjs. <- R1.* s1_cjs

#popula_cjs  <- m2  / p2   # having some problems getting this one to work today Apr2014
#popula_cjs. <- m2. / p2

###<<!!These Smolts ROUNDED AS PER OLD BOOTSTRAP PROGRAM!!>>
#t0 smolts at LGR
t0_cjs      <- roundT(tlgr_cjs + tlgs_cjs  + tlmn_cjs)
t0_cjsMod   <- roundT(tlgr_cjs + tlgs_cjs  + tlmn_cjsMod)

#t1 smolts at LGR
t1_cjs      <- roundT(tlgr_cjs + tlgs2_cjs + tlmn2_cjs)
t1_cjsMod   <- roundT(tlgr_cjs + tlgs2_cjs + tlmn2_cjsMod)

#c0 smolts at LGR
c0_cjs      <- roundT(popula_cjs - m12 - m13 / s2_cjs - m14 /(s2_cjs * s3_cjs) - delta0)
c0_cjsMod   <- roundT(popula_cjs - m12 - m13 / s2_cjs - m14 /(s2s3_cjsMod)     - delta0)


#c1 smolts at LGR [[old]]
# c1_cjs      <- roundT(m12. - delta2. + (m13. - delta3.) / s2_cjs + (m14. - delta4.) /(s2_cjs * s3_cjs) - delta1)
# c1_cjsMod   <- roundT(m12. - delta2. + (m13. - delta3.) / s2_cjs + (m14. - delta4.) /(s2s3_cjsMod    ) - delta1)

#c1 smolts at LGR [[new]]
c1_cjsNEW     <- roundT(R1. * s1_cjs * (p2 + (1-p2)*p3 + (1-p2)*(1-p3)*p4) -
    (delta2. + (delta3. / s2_cjs) + (delta4. / (s2_cjs * s3_cjs)) + delta1))
c1_cjsNEWMod  <- roundT(R1. * s1_cjs * (p2 + (1-p2)*p3 + (1-p2)*(1-p3)*p4) -
		(delta2. + (delta3. / s2_cjs) + (delta4. / (s2s3_cjsMod    )) + delta1))

#SAR, TIR, and D----------------------------------------------------------------
#a.k.a. C0 SAR
sar_c0_cjs   <- c0adults  / c0_cjs
sar_c0_cjsMod<- c0adults  / c0_cjsMod

#a.k.a. C1 SAR [[new version]]
sar_c1_cjs   <- c1adults.  / c1_cjsNEW
sar_c1_cjsMod<- c1adults.  / c1_cjsNEWMod

#a.k.a. T0 SAR
sart0cjsad   <- roundT(t0adults.  / t0_cjs   , digits = 6)
sart0cjsadMod<- roundT(t0adults.  / t0_cjsMod, digits = 6)

#a.k.a. Tx SAR
sart1cjsad   <- roundT(t1adults.  / t1_cjs   , digits = 6)
sart1cjsadMod<- roundT(t1adults.  / t1_cjsMod, digits = 6)

#a.k.a. TIR Tx / C0
t1_c0_cjsu   <- sart1cjsad / sar_c0_cjs
t1_c0_cjsuMod<- sart1cjsad / sar_c0_cjsMod

#a.k.a. D Tx & C0
d_tm1_cjsu   <- t1_c0_cjsu * (vc_cjs / vt_cjs_NEW)
d_tm1_cjsuMod<- t1_c0_cjsu * (vc_cjs / vt_cjs_NEWMod)

#overall SARs by location and with and without jacks
sar_tws_cr <- totaladult./(s1_cjs * R1.)
sar_tws_cr_gj <- totaladult_gj./(s1_cjs * R1.)
sar_tws_cr_b <-  totaladult_b./(s1_cjs * R1.)
sar_tws_cr_bj <- totaladult_bj./(s1_cjs * R1.)

#NEW Proportions----------------------------------------------------------------
E.c0          <-  popula_cjs. * (1-p2) * (1-p3) * (1-p4)
E.c1          <-  popula_cjs. * (p2 + (1 - p2) * p3
                                    + (1 - p2) * (1 - p3) * p4)-
                       (delta2. +
                        delta3. /  s2_cjs +
                        delta4. / (s2_cjs * s3_cjs))


pr_trans_new  <-  t1_cjs    / (E.c0 + t1_cjs + E.c1)
pr_c0_new     <-  E.c0      / (E.c0 + t1_cjs + E.c1)
pr_c1_new     <-  E.c1      / (E.c0 + t1_cjs + E.c1)

#dam transport SARs-------------------------------------------------------------
sarLGR. <- lgradults. /  x12.
sarLGS. <- lgsadults. / (x1a2.)
sarLMN. <- lmnadults. / (x1aa2.)

#not done yet . . .
##with delayed transportation, estimate with T0 + T1 =or=> Tx--------------------
#SAR_Tlgr  <-  get.CIs(dat$tlgradults/dat$x12)
#Tlgradults<-  dat$tlgradults[1]
#Tlgrsmolts<-  dat$x12[1]
#
#SAR_Tlgs  <- get.CIs(dat$lgsadults2/(dat$x1a2/dat$s2_cjs))
#Tlgsadults<- dat$lgsadults2[1]
#Tlgssmolts<- (dat$x1a2/dat$s2_cjs)[1]
#
#SAR_Tlmn  <- get.CIs(dat$lmnadults2/(dat$x1aa2/(dat$s2_cjs * dat$s3_cjs)))
#Tlmnadults<- dat$lmnadults2[1]
#Tlmnsmolts<- (dat$x1aa2/(dat$s2_cjs * dat$s3_cjs))[1]
#
#
##without preassignment or delayed transportation--------------------------------
#SAR_Tlgr  <-  get.CIs(dat$tlgradults/dat$x12)
#Tlgradults<-  dat$tlgradults[1]
#Tlgrsmolts<-  dat$x12[1]
#
#SAR_Tlgs  <-  get.CIs(dat$tlgsadults/(dat$x102/dat$s2_cjs))
#Tlgsadults<-  dat$tlgsadults[1]
#Tlgssmolts<-  (dat$x102/dat$s2_cjs)[1]
#
#SAR_Tlmn  <-  get.CIs(dat$tlmnadults/(dat$x1002/(dat$s2_cjs * dat$s3_cjs)))
#Tlmnadults<-  dat$tlmnadults[1]
#Tlmnsmolts<-  (dat$x1002/(dat$s2_cjs * dat$s3_cjs))[1]

#load up answers----------------------------------------------------------------
releaseCRT <- fill.CIs(R1)
# releaseT   <- fill.CIs(t$r1)
# releaseR   <- fill.CIs(r$r1)
# BSreaches  <- fill.CIs(crt$numreaches)
doCalcsreaches <- fill.CIs(reaches)

#
releaseT   <- get.CIs(R1.)
S1 <- get.CIs(s1_cjs)
S2.3.4 <- get.CIs(s2_cjs * s3_cjs * s4_cjs)
S2 <- get.CIs(s2_cjs)
S3 <- get.CIs(s3_cjs)
S4 <- get.CIs(s4_cjs)
S5.6  <- get.CIs(s5_cjs * s6_cjs)
S5 <- get.CIs(s5_cjs)
S6 <- get.CIs(s6_cjs)
S2.3 <- get.CIs(s2_cjs * s3_cjs)

SR_0expan <- get.CIs(vc_cjs_0expan)
SR_1expan <- get.CIs(vc_cjs_1expan)
SR_2expan <- get.CIs(vc_cjs_2expan)
SR_3expan <- get.CIs(vc_cjs_3expan)

p2 <- get.CIs(p2)
p3 <- get.CIs(p3)
p4 <- get.CIs(p4)
p5 <- get.CIs(p5)
p6 <- get.CIs(p6)
p7 <- get.CIs(p7)

x12   <- get.CIs(x12.)
x1a2  <- get.CIs(x1a2.)
x1aa2 <- get.CIs(x1aa2.)

delta2 <- get.CIs(delta2.)
delta3 <- get.CIs(delta3.)
delta4 <- get.CIs(delta4.)    # got an issue here (what issue?)

SR   <- get.CIs(vc_cjs)
C0   <- get.CIs(sar_c0_cjs* 100, exci='y',
  exactci, x=c0adults[1], n=c0_cjs[1])
C1   <- get.CIs(sar_c1_cjs* 100, exci='y',
  exactci, x=c1adults.[1], n=c1_cjsNEW[1])
T0   <- get.CIs(sart0cjsad* 100, exci='y',
  exactci, x=t0adults.[1], n=t0_cjs[1])
Tx   <- get.CIs(sart1cjsad* 100, exci='y',
  exactci, x=t1adults.[1], n=t1_cjs[1])
TIR  <- get.CIs(t1_c0_cjsu)
D    <- get.CIs(d_tm1_cjsu)

C0Modified   <- get.CIs(sar_c0_cjsMod* 100, exci='y',
  exactci, x=c0adults[1], n=c0_cjsMod[1])
C1Modified   <- get.CIs(sar_c1_cjsMod* 100, exci='y',
  exactci, x=c1adults.[1], n=c1_cjsNEWMod[1])
T0Modified   <- get.CIs(sart0cjsadMod* 100, exci='y',
  exactci, x=t0adults.[1], n=t0_cjsMod[1])
TxModified   <- get.CIs(sart1cjsadMod* 100, exci='y',
  exactci, x=t1adults.[1], n=t1_cjsMod[1])
TIRModified  <- get.CIs(t1_c0_cjsuMod)
DModified    <- get.CIs(d_tm1_cjsuMod)

Pr_T         <- get.CIs(pr_trans_new)
Pr_C0        <- get.CIs(pr_c0_new)
Pr_C1        <- get.CIs(pr_c1_new)

overallSAR   <- get.CIs(sar_tws_cr* 100, exci='y',
  exactci, x=totaladult.[1], n=popula_cjs.[1])
overallSAR_gj   <- get.CIs(sar_tws_cr_gj* 100, exci='y',
  exactci, x=totaladult_gj.[1], n=popula_cjs.[1])
overallSAR_b   <- get.CIs(sar_tws_cr_b* 100, exci='y',
  exactci, x=totaladult_b.[1], n=popula_cjs.[1])
overallSAR_bj   <- get.CIs(sar_tws_cr_bj* 100, exci='y',
  exactci, x=totaladult_bj.[1], n=popula_cjs.[1])
LGRpop       <- get.CIs(popula_cjs.)

sarLGR <- get.CIs(sarLGR.* 100, exci='y',
  exactci, x=lgradults.[1], n=x12.[1])
sarLGS <- get.CIs(sarLGS.* 100, exci='y',
  exactci, x=lgsadults.[1], n=x1a2.[1])
sarLMN <- get.CIs(sarLMN.* 100, exci='y',
  exactci, x=lmnadults.[1], n=x1aa2.[1])

adultLGR <- get.CIs(lgradults.)
adultLGS <- get.CIs(lgsadults.)
adultLMN <- get.CIs(lmnadults.)

MCNtransLGRequival <- get.CIs(tmcn2_cjs)

T0pop         <- get.CIs(t0_cjs)
T0popModified <- get.CIs(t0_cjsMod)
Txpop         <- get.CIs(t1_cjs)
TxpopModified <- get.CIs(t1_cjsMod)
C0pop         <- get.CIs(c0_cjs)
C0popModified <- get.CIs(c0_cjsMod)

C0adults  <- get.CIs(c0adults)
C1adults  <- get.CIs(c1adults.)
T0adults  <- get.CIs(t0adults.)
Txadults  <- get.CIs(t1adults.)

adults_g  <- get.CIs(totaladult.)
adults_gj  <- get.CIs(totaladult_gj.)
adults_b  <- get.CIs(totaladult_b.)
adults_bj  <- get.CIs(totaladult_bj.)

C1pop         <- get.CIs(c1_cjsNEW)
C1popModified <- get.CIs(c1_cjsNEWMod)

#wrangle format and output------------------------------------------------------
#choose parameters:
parm <- c(
  # general stuff for evaluation; has various reach survival versions
  'doCalcsreaches', #'BSreaches',
  'releaseCRT', 'releaseT', #'releaseR',
  'SR_0expan', 'SR_1expan', 'SR_2expan', 'SR_3expan',
  'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'SR', 'S2.3',
  'S2.3.4', 'S5.6', ##adding for LGR-MCN S2.3.4, and S5.6 for MCN-BON

  'p2', 'p3', 'p4', 'p5', 'p6', 'p7',

  # overall SARs [from T group]
  'overallSAR', 'overallSAR_gj', 'overallSAR_b', 'overallSAR_bj',

  # for SARs by route of passage or component SARs and proportion for each
  'TIR', 'D', 'C0', 'C1', 'Tx',

  # adults
  'C0adults', 'C1adults', 'Txadults',
  'adults_g', 'adults_gj', 'adults_b', 'adults_bj',

  # smolt populations
  'LGRpop', 'C0pop', 'C1pop', 'Txpop',

  # with modified s2*s3 truncated at 1.00
  'TIRModified', 'DModified',
  'C0Modified', 'C1Modified', 'TxModified',
  'C0popModified', 'C1popModified', 'TxpopModified',

  # for the proportion transport appendix
  'delta2', 'delta3', 'delta4',
  'Pr_T', 'Pr_C0', 'Pr_C1',

  # for the transport SARS appendix
  'x12', 'x1a2', 'x1aa2',
  'sarLGR', 'adultLGR',
  'sarLGS', 'adultLGS',
  'sarLMN', 'adultLMN',

  # LGR equivalents of McNary Transports [T group]
  'MCNtransLGRequival'
  )

ansDF <- get(parm[1])
for(i in parm[-1]) {ansDF <- rbind(ansDF, get(i))}
ans <- cbind(parm, ansDF)


# option to save output to working directory and sql server
if (makefile == 'y' | makefile == 'Y') {
  # save bootystrap output to sql server
  rand_str <- function(n) {
    do.call(paste0, replicate(7, sample(c(letters,letters,0:9), n, TRUE), FALSE))
  }
  channel <- odbcDriverConnect("case=nochange;
                                Description=CSSOUTPUT;
                                DRIVER=SQL Server;
                                SERVER=PITTAG_2016SQL;
                                UID=sa;
                                PWD=frznool;
                                WSID=CUTTHROAT;
                                DATABASE=CSSOUTPUT;
                                Network=DBMSSOCN")
  sqlSave(channel, data.frame(crt), tablename=paste0('CRT_', target, '_bootylator_', rand_str(1), format(Sys.time(), '%m%d%Y')) )
  odbcCloseAll()

  # write doCalcs output in csv file to working directory
  nm <- paste("doCalcs_v3 ", target, format(Sys.time(), '%Y-%m-%d %H%M%S'),
    " ",reaches," reaches",".csv", sep = "")
  # ans$css_group<-css_group
  write.table(ans, paste(nm, sep = "\\"),
    col.names = T, row.names = F, sep = ',', quote = F)

  # also send to clipboard
  # write.table(ans,
  #             'clipboard',
  #             col.names = T, row.names = F, sep = ',',
  #             quote = F)

  # save to sql server CSS report
  channel2 <- odbcDriverConnect("case=nochange;
                                Description=CSSREPORT;
                                DRIVER=SQL Server;
                                SERVER=PITTAG_2016SQL;
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
  inthesqltable<- sqlQuery(channel2, paste("select css_group, migr_year, parm, thecount = count(initial) from BOOTSTRAP_RESULTS where css_group = ", "'", css_group, "'", " and migr_year = ", "'", migr_yr, "'", " and parm = 'overallSAR' group by css_group, migr_year, parm", sep=""))

  # EITHER UPDATE OR SAVE DEPENDING ON ANSWER
  if(length(inthesqltable$thecount) != '0') {
    sqlUpdate(channel2, data.frame(ans),tablename="BOOTSTRAP_RESULTS", index = c('css_group', 'migr_year', 'parm'), verbose = FALSE, test = FALSE, nastring = NULL,fast = TRUE)
  } else {
      sqlSave(channel2,data.frame(ans),tablename="BOOTSTRAP_RESULTS",safer=TRUE,append=TRUE)
  }

  odbcCloseAll()
}

# finally print in R
# print(format(ans, scientific = F))
rownames(ans)<- parm
return(ans)
}
#------------------------------------------------------------------------------
