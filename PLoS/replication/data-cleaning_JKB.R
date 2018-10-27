########################
# data-cleaning_JKB.R
# 
# Purpose: 
# Takes in MSM_KingCounty_rev.csv and turns it into 
# a data frame, msm, with additional variables
# 
# Dependencies:
# MSM_KingCounty_rev.csv
#
# History: 
# This file is a commented version of Ian's data-cleaning.R
########################

#setwd('/Users/jeanette/Dropbox/School/PhD/HIV_WA')
#msm <- read.csv("data/MSM_KingCounty_rev.csv",na.string="",stringsAsFactor=FALSE)

# Year, quarter, and quarter-year of Dx (diagnosis)
msm$yearDx <- as.numeric(substring(msm$HIVdxQTR,0,4))
msm$quarterDx <- as.numeric(substring(msm$HIVdxQTR,5,5))
msm$timeDx <- msm$yearDx + (msm$quarterDx-1)/4
# AIDS at Dx - if missing, assumed to be false 
msm$aidsAtDx <- msm$HIVdxQTR == msm$AIDSdxQTR
msm$aidsAtDx[is.na(msm$aidsAtDx)] <- FALSE
# Year, quarter, and quarter-year of AIDS (diagnosis)
msm$yearAids <- as.numeric(substring(msm$AIDSdxQTR,0,4))
msm$quarterAids <- as.numeric(substring(msm$AIDSdxQTR,5,5))
msm$timeAids <- msm$yearAids + (msm$quarterAids-1)/4

# Was subject ever tested previously, combining eHARS-HIS and PCRS hierarchically
# Y or N or NA (no response, unknowns and refusals)
msm$everTested <- msm$ungtst
msm$everTested <- ifelse(msm$everTested=="Y" | is.na(msm$EverTestedpcrs),
                         msm$everTested, msm$EverTestedpcrs)
msm$everTested[msm$everTested=="U"] <- NA
msm$everTested[msm$everTested=="R"] <- NA

# Create last negative variable
# If count in days for any of the systems is <0, it's missing
msm$lagHIV_LastNegEhars[msm$lagHIV_LastNegEhars<0] <- NA
msm$lagHIV_HISNeg[msm$lagHIV_HISNeg<0] <- NA
msm$lagHIV_LastNegPCRSrpt[msm$lagHIV_LastNegPCRSrpt<0] <- NA
# Record the most recent negative test from all sources
msm$lastNegMin <- with(msm,pmin(lagHIV_LastNegEhars, lagHIV_HISNeg, 
                             lagHIV_LastNegPCRSrpt,na.rm=TRUE))
# In a new variable, prioritize EHARS if it's non-missing
msm$lastNeg <- with(msm,ifelse(is.na(lagHIV_LastNegEhars),lastNegMin,lagHIV_LastNegEhars))
msm$lastNeg <- msm$lastNeg/365
msm$lastNegMin <- msm$lastNegMin/365

# everTested - set to Y those who have a date of last test
# even though their everTested was N or NA (39 cases), 
# and turn into a logical 
table(msm$everTested, !is.na(msm$lastNeg), useNA='ifany')
msm$everTested[!is.na(msm$lastNeg)] <- "Y"
msm$everTested <- msm$everTested=="Y"

# Set at-risk period (infPeriod): for those with testing history, 
# time since last negative test, and for others, 
# time since age 16 (sexual debut assumption). In both 
# cases, at-risk time is capped at aidsUB (~18 yrs). 
aidsUB <- qweibull(.95,shape=2.516,scale=1/0.086) #17.98418
msm$infPeriod <- ifelse(!is.na(msm$everTested) & !msm$everTested, # Not NA and never tested
                      pmin(msm$hiv_age_yrs-16,aidsUB), msm$lastNeg)
msm$infPeriod <- pmin(msm$infPeriod,aidsUB)

msm$infPeriod_imputeNA <- ifelse(is.na(msm$infPeriod),
                            pmin(msm$hiv_age_yrs-16,aidsUB), 
                            msm$infPeriod)
msm$infPeriod_imputeNA <- pmin(msm$infPeriod_imputeNA,aidsUB)

