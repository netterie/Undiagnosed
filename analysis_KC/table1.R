
########################
# table1.R
# 
# Purpose: reconcile counts reported in Table 1
# of the paper and later tabulations
# 
# 
# Dependencies:
# data-cleaning_JKB.R
#
# History: code from email exchange 'Effect of 289 cases (fwd)'
########################


######################################################################
# SETUP and DATA CLEANING
######################################################################

# TEMPORARY: SOURCE FUNCTIONS
source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/other.R')

# Load libraries, data and data-cleaning file
# Eventually this function should return the cleaned dataset,
# but data-cleaning has the name hardcoded as msm and I'm not
# going to generalize that right now
setup_hivbackcalc(workd='/Users/jeanette/Dropbox/School/PhD/HIV_WA',
                  datafile='data/MSM_KingCounty_rev.csv',
                  source_these='analysis/data-cleaning_JKB.R',
                  loadlib=FALSE)

#str(msm)

######################################################################
# INVESTIGATE TABULATIONS - Martina & Ian's code
######################################################################

# Martina's

ehars <- msm$lagHIV_LastNegEhars
his <- msm$lagHIV_HISNeg
ps <- msm$lagHIV_LastNegPCRSrpt
#cor(cbind(ehars, his, ps), use="pairwise.complete.obs")
#pairs(cbind(ehars,his, ps), ylim=c(0,5000), xlim=c(0,5000))
table(is.na(his[is.na(ehars)])) # 382 obs
table(ps[is.na(ehars)] <  his[is.na(ehars)]) # 64 obs from PS rest
table(ehars=is.na(ehars), his=is.na(his), ps=is.na(ps))

# Ian's
#total n
nrow(msm)

#number with both HIS and PS
with(msm,table(!is.na(lagHIV_HISNeg) , !is.na(lagHIV_LastNegPCRSrpt)))

#number of HIS used
tmp <- with(msm[!is.na(msm$lagHIV_HISNeg) & is.na(msm$lagHIV_LastNegEhars),],
                 lagHIV_HISNeg <= lagHIV_LastNegPCRSrpt)
tmp[is.na(tmp)] <- TRUE
sum(tmp)

#number of PS used
tmp <- with(msm[!is.na(msm$lagHIV_LastNegPCRSrpt) & is.na(msm$lagHIV_LastNegEhars),],
            lagHIV_HISNeg > lagHIV_LastNegPCRSrpt)
tmp[is.na(tmp)] <- TRUE
sum(tmp)

#number where PS == HIS (table assumes that HIS takes presidence)
tmp <- with(msm[!is.na(msm$lagHIV_HISNeg) & is.na(msm$lagHIV_LastNegEhars),],
            lagHIV_HISNeg == lagHIV_LastNegPCRSrpt)
tmp[is.na(tmp)] <- TRUE
sum(tmp)

670+80+382

#never tested
sum(!is.na(msm$everTested) & !msm$everTested)


######################################################################
# INVESTIGATE TABULATIONS - my code
######################################################################

# Ever tested info (Y=yes, N=no, R=refused, U=unknown)
with(msm, table(eharsHIS=ungtst, pcrs=EverTestedpcrs, useNA='ifany'))

lastNeg <- data.frame(ehars=msm$lagHIV_LastNegEhars, 
                      his=msm$lagHIV_HISNeg, 
                      ps=msm$lagHIV_LastNegPCRSrpt)

# For me, the key to getting unconfused was to make comparisons
# by constructing variables that do not have NA values
lastNeg <- within(lastNeg, { 
                  ehars_noMiss <- !is.na(ehars)
                  his_noMiss <- !is.na(his)
                  ps_noMiss <- !is.na(ps)
                  his_and_ps <- his_noMiss & ps_noMiss 
                  his_and_ehars <- his_noMiss & ehars_noMiss 
                  ps_and_ehars <- ps_noMiss & ehars_noMiss
                  his_and_ps_no_ehars <- his_noMiss & ps_noMiss & !ehars_noMiss
                  his_and_ehars_no_ps <- his_noMiss & ehars_noMiss & !ps_noMiss
                  ps_and_ehars_no_his <- ps_noMiss & ehars_noMiss & !his_noMiss
                  his_no_ehars <- his_noMiss & !ehars_noMiss
                  his_no_ps <- his_noMiss & !ps_noMiss
                  ps_no_his <- ps_noMiss & !his_noMiss
                  ps_no_ehars <- ps_noMiss & !ehars_noMiss
                  ehars_only <- ehars_noMiss & !his_noMiss & !ps_noMiss
                  his_only <- his_noMiss & !ehars_noMiss & !ps_noMiss
                  ps_only <- ps_noMiss & !ehars_noMiss & !his_noMiss
                  his_lessOReq_ps <- his<=ps & his_and_ps_no_ehars
                  his_greaterthan_ps <- his>ps & his_and_ps_no_ehars
                  atleast_onesource <- ehars_noMiss | his_noMiss | ps_noMiss
                  nosource <- !atleast_onesource
                  allsources <- ehars_noMiss & his_noMiss & ps_noMiss
                  use_his_not_ps <- his_lessOReq_ps | his_no_ps
                  use_ps_not_his <- his_greaterthan_ps | ps_no_his
                  use_his_final <- his_no_ehars & use_his_not_ps
                  use_ps_final <- ps_no_ehars & use_ps_not_his
                  })

                     
sums=data.frame(colSums(lastNeg))































