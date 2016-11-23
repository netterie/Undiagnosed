########################
# format_data.R
# 
# Purpose: 
# Explores and formats WA state testing history
# data
#
# Dependencies:
# None
#
# History: 
# Based off of Ian's data-cleaning.R
#
# Notes:
# Meant to be called by other files that have loaded the data
# Uses read_chunks() notation (## ---- tag ----) so that
# it can be used in reports as well as sourced by other
# R files. See format_data_report.Rnw
########################


#############################################################
# SETUP
#############################################################

# Note: I realized that the standalone indicator is actually
# not necessary as long as the report that read_chunk()s this
# file does not run this initial chunk, but instead has its
# own code for setup. read_chunk() is pretty brilliant! But, 
# just to remind myself to be careful in general, 
# I am setting standalone as FALSE.

standalone <- FALSE
if (standalone) {
    rm(list=ls())
    workd <- '/Users/jeanette/Dropbox/School/PhD/HIV_WA'
    #filename <- 'WA_BACKCALC_DATA_v2.csv'
    filename <- 'wa_backcalc_data_201506.csv'
    dataf <- read.csv(file.path(workd,'data', filename), 
                      stringsAsFactors=FALSE, na.strings="")
    library(plyr)
}

## ---- printSummaries ----
# I'm setting an indicator for printing summaries because
# they're printing in reports even when echo=FALSE...
printSummaries <- FALSE

## ---- overview ----
#############################################################
# OVERVIEW: N and expected variables
#############################################################

Nobs0 <- nrow(dataf)
colnames(dataf) <- tolower(colnames(dataf))
variables_in_codebook <- tolower(c('hdx_age',
                           'NEW_RACE',
                           'NEW_MODE',
                           'HDX_YR_QTR',
                           'HDX_DT_FLAG',
                           'ADX_YR_QTR',
                           'ADX_DT_FLAG',
                           'TTH_EVER_NEG',
                           'LAG_LNEG_HDX_DT',
                           'TTH_LNEG_DT_FLAG',
                           'TTH_PREV_POS',
                           'LAG_PPOS_HDX_DT',
                           'TTH_PPOS_DT_FLAG',
                           'METH_USE',
                           'CD4_DAYS',
                           'hst',
                           'FirstCD4cnt',
                           'VL_DAYS',
                           'FirstVL'))
not_in_dataf <- variables_in_codebook[!variables_in_codebook%in%colnames(dataf)]
not_in_codebook <- colnames(dataf)[!colnames(dataf)%in%variables_in_codebook]

# Factors and labels
races <- c('White', 'Black', 'Hisp', 'Asian', 'NHoPI', 
           'AI/AN', 'Multi', 'Unknown')
modes <- c('MSM', 'IDU', 'MSM/IDU', 'Transfus', 'Hemo', 'Hetero', 
           'Ped', 'F Pres Hetero', 'NIR')
dataf <- transform(dataf,
                   new_race=factor(new_race, levels=as.character(1:8),
                                   labels=races),
                   new_mode=factor(new_mode, levels=as.character(1:9),
                                   labels=modes))

## ---- summarize ----
#############################################################
# SUMMARIZE: Summarize or tabulate variables
#############################################################
if (printSummaries) {
    for (x in 1:ncol(dataf)) {
       cat('\n\n\nVARIABLE',x,':',colnames(dataf)[x],'\n     ')
       var = dataf[,x]
       if (length(unique(var))>10) {
           if (is.numeric(var)) {
               sum <- summary(var) 
               nmiss <- sum["NA's"]
           } else {
               nmiss <- table(var)[is.na(names(table(var)))]
               sum <- ''
           }
       } else {
           sum <- table(var, useNA='always')
           nmiss <- sum[is.na(names(sum))]
       }
       print(sum)
       nmiss <- as.numeric(nmiss)
       cat('\n     Percent missing:')
       print(round(100*nmiss/nrow(dataf),2))
    }
}


## ---- split_yr_qtr ----
#############################################################
# SPLIT COMBINED YR-QTR VARIABLE
#############################################################
# Year, quarter, and quarter-year of Dx (diagnosis)
dataf$yearDx <- as.numeric(substring(dataf$hdx_yr_qtr,0,4))
dataf$quarterDx <- as.numeric(substring(dataf$hdx_yr_qtr,6,6))
dataf$timeDx <- dataf$yearDx + (dataf$quarterDx-1)/4
# AIDS at Dx - if missing, assumed to be false 
dataf$aidsAtDx <- dataf$hdx_yr_qtr == dataf$adx_yr_qtr
dataf$aidsAtDx[is.na(dataf$aidsAtDx)] <- FALSE
# Year, quarter, and quarter-year of AIDS (diagnosis)
dataf$yearAids <- as.numeric(substring(dataf$adx_yr_qtr,0,4))
dataf$quarterAids <- as.numeric(substring(dataf$adx_yr_qtr,6,6))
dataf$timeAids <- dataf$yearAids + (dataf$quarterAids-1)/4

## ---- subset_initial ----
#############################################################
# SUBSET THE DATA - INITIAL RESTRICTIONS
#############################################################
if (!'year_min'%in%ls()) year_min <- 2005
if (!'year_max'%in%ls()) year_max <- 2013

# Year min and max for this run
c(year_min, year_max)

# Non-sequential look
table(hst_included=dataf$hst=='WA', useNA='ifany')
table(yearDx_included=dataf$yearDx>=year_min & dataf$yearDx<=year_max, 
      useNA='ifany')
table(yearDx_missing=is.na(dataf$hdx_yr_qtr))
table(age_missing_and_missing_lastNeg=(is.na(dataf$hdx_age) & 
                                       is.na(dataf$lag_lneg_hdx_dt)))
# Sequential look
(hst_included <- table(hst_included=dataf$hst=='WA', useNA='ifany'))
dataf <- subset(dataf, hst=='WA')
(yearDx_included <- table(yearDx_included=(dataf$yearDx>=year_min & dataf$yearDx<=year_max), useNA='ifany'))
dataf <- subset(dataf, yearDx>=year_min & yearDx<=year_max)
(age_included <- table(age_and_lastNeg_present=!(is.na(dataf$hdx_age) & 
                                                 is.na(dataf$lag_lneg_hdx_dt))))
dataf <- subset(dataf, !(is.na(hdx_age) & is.na(lag_lneg_hdx_dt)))
(Nobs1 <- nrow(dataf))

## ---- impute_qtr ----
#############################################################
# IMPUTE A QUARTER IF ONLY YEAR IS KNOWN
#############################################################
impute_qtr <- !is.na(dataf$yearDx) & is.na(dataf$quarterDx)
set.seed(98103)
dataf$quarterDx[impute_qtr] <- sample(4, size=sum(impute_qtr), 
                                      replace=TRUE)
dataf$timeDx <- dataf$yearDx + (dataf$quarterDx-1)/4
summary(dataf$timeDx, digits=6)

time_min <- min(dataf$timeDx)
time_max <- max(dataf$timeDx)

# Time min and max for this run
c(time_min, time_max)

## ---- collapse ----
#############################################################
# COLLAPSE RACE AND MODE OF DIAGNOSIS
#############################################################

race_levels <- c('White', 'Black', 'Hisp', 'Asian', 'Native', 'Multi')
mode_levels <- c('MSM', 'Hetero', 'Blood/Needle')
dataf <- within(dataf, {
                race <- as.character(new_race)
                race[race=='AI/AN' | race == 'NHoPI'] <- 'Native'
                race <- factor(race,
                                labels=race_levels,
                                levels=race_levels)
                mode <- as.character(new_mode)
                mode[mode=='MSM/IDU'] <- 'MSM'
                mode[mode=='F Pres Hetero' | mode=='NIR'] <- 'Hetero'
                mode[mode=='IDU'|mode=='Transfus'|mode=='Hemo'|
                      mode=='Ped'] <- 'Blood/Needle'
                mode <- factor(mode,
                                levels=mode_levels,
                                labels=mode_levels)
                mode2 <- factor(ifelse(mode=='MSM', 'MSM', 'non-MSM'))
                })


## ---- everHadNegTest ----
#############################################################
# CREATE everHadNegTest
#############################################################
# Define everHadNegTest based on tth_ever_neg
# 2015 data update: this variable was coded numerically, so I have
# added that option in. 
dataf <- transform(dataf, 
                  everHadNegTest=ifelse(tth_ever_neg=='Y' | tth_ever_neg==1, TRUE, 
                                        ifelse(tth_ever_neg=='N' | tth_ever_neg==2, FALSE, NA)))
with(dataf,table(everHadNegTest, tth_ever_neg, useNA='always'))

# Now cross-check it with the lag_lneg_hdx_dt, which actually has the 
# time since last negative test
(checkEver <- with(dataf,table(everHadNegTest, 
                               TID_NA=is.na(lag_lneg_hdx_dt), useNA='always')))

# Look at actual lag_lneg_hdx_dt values by everHadNegTest
ddply(dataf, .(everHadNegTest), function(x) c(summary(x$lag_lneg_hdx_dt)))


## ---- fix_everHadNegTest_toTRUE ----
toTRUE1 <- !dataf$everHadNegTest & !is.na(dataf$lag_lneg_hdx_dt)
toTRUE2 <- is.na(dataf$everHadNegTest) & !is.na(dataf$lag_lneg_hdx_dt)
dataf$everHadNegTest[toTRUE1] <- TRUE
dataf$everHadNegTest[toTRUE2] <- TRUE

## ---- fix_everHadNegTest_toFALSE ----
toFALSE <- dataf$everHadNegTest & is.na(dataf$lag_lneg_hdx_dt)
dataf$everHadNegTest[toFALSE] <- FALSE

## ---- check_everHadNegTest ----
(checkEver <- with(dataf,table(everHadNegTest, 
                               TID_NA=is.na(lag_lneg_hdx_dt), useNA='always')))

## ---- infPeriod ----
#############################################################
# CREATE infPeriod and then look at it
#############################################################

#### TEMPORARY: 
#dataf$age=35

aidsUB <- qweibull(.95,shape=2.516,scale=1/0.086) #17.98418
dataf <- within(dataf,{
                  lastNeg_yrs=lag_lneg_hdx_dt/365
                  infPeriod=ifelse(everHadNegTest,
                                   pmin(lastNeg_yrs, aidsUB),
                                   ifelse(!everHadNegTest,
                                          pmin(hdx_age-16, aidsUB),
                                          NA))
                  earliestInf=hdx_age-infPeriod
                  })

## ---- infPeriod_investigate_neg ----
# Number of cases who got a negative infPeriod
(neginfPeriod <- sum(dataf$infPeriod<0,na.rm=TRUE))
# Diagnoses at or under age 16 by everHadNegTest
(a1 <- table(atunder16=dataf$hdx_age<=16, 
             everHadNegTest=dataf$everHadNegTest, useNA='ifany'))
# Diagnoses at or under age 16 by year, 2005-2013
table(atunder16count=subset(dataf, yearDx>=year_min & yearDx<=year_max)$hdx_age<=16, 
      year=subset(dataf, yearDx>=year_min & yearDx<=year_max)$yearDx, useNA='ifany')
# Now just under 16, excluding hdx_age=16
# Diagnoses under age 16 by everHadNegTest
(a2 <- table(under16=dataf$hdx_age<16, 
             everHadNegTest=dataf$everHadNegTest, useNA='ifany'))
# Diagnoses under age 16 by year
table(under16count=subset(dataf, yearDx>=year_min & yearDx>=year_max)$hdx_age<16, 
      year=subset(dataf, yearDx>=year_min & yearDx>=year_max)$yearDx, useNA='ifany')
# Among those diagnosed at or under 16: everHadNegTest by mode
table(everHadNegTest=subset(dataf,hdx_age<=16)$everHadNegTest,
      mode=subset(dataf,hdx_age<=16)$new_mode, useNA='ifany')

## ---- infPeriod_cap_neg
(young_included <- with(dataf,
                       table(over16_or_atunder16_with_obs_infPeriod= 
                             (hdx_age>16 | 
                             !(hdx_age<=16 & (!everHadNegTest |
                                             is.na(everHadNegTest)))))))
dataf <- subset(dataf, !(hdx_age<=16 & (!everHadNegTest | 
                                        is.na(everHadNegTest))))
(Nobs2 <- nrow(dataf))
summary(dataf$infPeriod, digits=3)

## ---- infPeriod_cap_tested ----
# We did cap some people whose TID's were >aidsUB
(check_cap1 <- with(subset(dataf, everHadNegTest), 
                    table(original_over_aidsUB=lastNeg_yrs>aidsUB, 
                          infPeriod_over_aidsUB=infPeriod>aidsUB, 
                          useNA='ifany')))

## ---- infPeriod_cap_nottested ----
(check_cap2 <- with(subset(dataf, !everHadNegTest), 
                    table(original_over_aidsUB=lastNeg_yrs>aidsUB, 
                          infPeriod_over_aidsUB=infPeriod>aidsUB, 
                          useNA='ifany')))

## ---- infPeriod_cap_NAtested ----
(check_cap3 <- with(subset(dataf, is.na(everHadNegTest)), 
                    table(original_over_aidsUB=lastNeg_yrs>aidsUB, 
                          infPeriod_over_aidsUB=infPeriod>aidsUB, 
                          useNA='ifany')))

## ---- fix16cases ----
# This is no longer needed because these 16 cases were fixed in the 
# 'fix_everHadNegTest_toTRUE' section
# CAREFUL because running this after the "infPeriod_investigate_zero" section below
# will overwrite the changes to those zeroes
these_cases <- with(dataf, is.na(everHadNegTest) & !is.na(lastNeg_yrs))
sum(these_cases)
dataf$everHadNegTest[these_cases] <- TRUE 
dataf$infPeriod[these_cases] <- dataf$lastNeg_yrs[these_cases]
with(subset(dataf, is.na(everHadNegTest)), table(original_over_aidsUB=lastNeg_yrs>aidsUB,
                 infPeriod_over_aidsUB=infPeriod>aidsUB,
                 useNA='ifany'))

## ---- infPeriod_investigate_zero ----
# Cases who still have a zero infPeriod - treat like missing
zeroinf <- dataf$infPeriod==0 & !is.na(dataf$infPeriod)
(table(dataf$everHadNegTest[zeroinf], useNA='ifany'))
# Change their everHadNeg flag to NA and their infPeriod to NA,
# since TID=0 does not make sense
dataf$everHadNegTest[zeroinf] <- NA
dataf$infPeriod[zeroinf] <- NA

## ---- imputeinfPeriod ----
#############################################################
# CREATE infPeriod_imputeNA and then look at it
#############################################################
dataf <- within(dataf,{ 
                infPeriod_imputeNA=ifelse(is.na(everHadNegTest),
                                          pmin(hdx_age-16, aidsUB),
                                          infPeriod)
                  })

with(subset(dataf, is.na(everHadNegTest)), summary(infPeriod))
with(subset(dataf, is.na(everHadNegTest)), summary(infPeriod_imputeNA))


## ---- imputeinfPeriod_cap_neg ----
summary(dataf$infPeriod_imputeNA, digits=3)


## ---- final_summary ----
nrow(dataf)

if (printSummaries) {
    for (var in c('hdx_age', 'timeDx', 'everHadNegTest', 
                  'lastNeg_yrs', 'infPeriod')) {
        cat('\nVARIABLE:', var, '\n')
        print(summary(dataf[,var]))
    }
}


## ---- agegroups ----
#############################################################
# DEFINE AGE GROUPS
#############################################################
dataf <- transform(dataf,
                   agecat5=cut(hdx_age,
                              breaks=c(0,seq(20,70,by=5),85),
                              include.lowest=TRUE,
                              right=TRUE,
                              labels=c('<=20',
                                       '21-25',
                                       '26-30',
                                       '31-35',
                                       '36-40',
                                       '41-45',
                                       '46-50',
                                       '51-55',
                                       '56-60',
                                       '61-65',
                                       '66-70',
                                       '71-85')))

