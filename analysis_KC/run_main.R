
########################
# run_main.R
# 
# Purpose: run the analysis using the new wrapper R functions
# in HIVBackCalc/R/other, package version 1.1
# 
# 
# Dependencies:
# data-cleaning_JKB.R
# packages HIVBackCalc, reshape2, ggplot2, scales, Hmisc, plyr
#
# History: based on run_JKB.R
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
                  source_these='analysis/data-cleaning_JKB.R')

#str(msm)

######################################################################
# RUN BACK-CALCULATION
######################################################################


# Run both cases: first the base case, then the upper bound

all_noimpute <- runBackCalc(TID=msm$infPeriod, 
                   impute=FALSE,
                   age=msm$hiv_age_yrs,
                   diagnosedCounts=table(msm$timeDx), 
                   upperBound=FALSE, 
                   runBoth=TRUE,
                   intervalLength=0.25, 
                   printProgress=FALSE) 

all_impute <- runBackCalc(TID=msm$infPeriod, 
                   impute=TRUE,
                   age=msm$hiv_age_yrs,
                   diagnosedCounts=table(msm$timeDx), 
                   upperBound=FALSE, 
                   runBoth=TRUE,
                   intervalLength=0.25, 
                   printProgress=FALSE) 

######################################################################
# RUN BACK-CALCULATION
######################################################################

summaries_noimpute <- summarize_runBackCalc(results=all_noimpute,
                                   diagnosedCounts=table(msm$timeDx),
                                   times=seq(2006,2012.75,by=0.25))

summaries_impute <- summarize_runBackCalc(results=all_impute,
                                   diagnosedCounts=table(msm$timeDx),
                                   times=seq(2006,2012.75,by=0.25))



######################################################################
# SAVE RESULTS
######################################################################

# Stats
stats = data.frame(imputed=c(rep('Yes',nrow(summaries_impute[[1]])),
                             rep('No',nrow(summaries_noimpute[[1]]))),
                   rbind(summaries_impute[[1]],
                         summaries_noimpute[[1]]))

write.csv(stats, file='analysis/results/run_main.csv',
          row.names=FALSE)

