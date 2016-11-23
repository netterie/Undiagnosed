msm <- read.csv("data/MSM_KingCounty_rev.csv",na.string="",stringsAsFactor=FALSE)

msm$yearDx <- as.numeric(substring(msm$HIVdxQTR,0,4))
msm$quarterDx <- as.numeric(substring(msm$HIVdxQTR,5,5))
msm$timeDx <- msm$yearDx + (msm$quarterDx-1)/4
msm$aidsAtDx <- msm$HIVdxQTR == msm$AIDSdxQTR
msm$aidsAtDx[is.na(msm$aidsAtDx)] <- FALSE
msm$yearAids <- as.numeric(substring(msm$AIDSdxQTR,0,4))
msm$quarterAids <- as.numeric(substring(msm$AIDSdxQTR,5,5))
msm$timeAids <- msm$yearAids + (msm$quarterAids-1)/4

#was subject ever tested previously
msm$everTested <- msm$ungtst
msm$everTested <- ifelse(msm$everTested=="Y" | is.na(msm$EverTestedpcrs),
                         msm$everTested, msm$EverTestedpcrs)
msm$everTested[msm$everTested=="U"] <- NA
msm$everTested[msm$everTested=="R"] <- NA

#create last negative variable
msm$lagHIV_LastNegEhars[msm$lagHIV_LastNegEhars<0] <- NA
msm$lagHIV_HISNeg[msm$lagHIV_HISNeg<0] <- NA
msm$lagHIV_LastNegPCRSrpt[msm$lagHIV_LastNegPCRSrpt<0] <- NA
msm$lastNegMin <- with(msm,pmin(lagHIV_LastNegEhars, lagHIV_HISNeg, 
                             lagHIV_LastNegPCRSrpt,na.rm=TRUE))
msm$lastNeg <- with(msm,ifelse(is.na(lagHIV_LastNegEhars),lastNegMin,lagHIV_LastNegEhars))
msm$lastNeg <- msm$lastNeg/365
msm$lastNegMin <- msm$lastNegMin/365

msm$everTested[!is.na(msm$lastNeg)] <- "Y"
msm$everTested <- msm$everTested=="Y"

aidsUB <- qweibull(.95,shape=2.516,scale=1/0.086)
msm$infPeriod <- ifelse(!is.na(msm$everTested) & !msm$everTested,
                      pmin(msm$hiv_age_yrs-16,aidsUB), msm$lastNeg)
msm$infPeriod <- pmin(msm$infPeriod,aidsUB)

