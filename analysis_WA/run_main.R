
########################
# run_main.R
# 
# Purpose: run the analysis using the new wrapper R functions
# in HIVBackCalc/R/other, package version 1.1
# 
# 
# Dependencies:
# format_data.R
# packages HIVBackCalc, rootSolve, reshape2, ggplot2, scales, Hmisc, plyr
#
# History: based on run_JKB.R
########################

## ---- setup ----
######################################################################
# SETUP and DATA CLEANING
######################################################################
#rm(list=ls())
# TEMPORARY: SOURCE FUNCTIONS
source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/other.R')

# Load libraries, data and data-cleaning file
# Eventually this function should return the cleaned dataset,
# but data-cleaning has the name hardcoded as msm and I'm not
# going to generalize that right now
setup_hivbackcalc(workd='/Users/jeanette/Dropbox/School/PhD/HIV_WA',
                  datafile='data/WA_BACKCALC_DATA_v2.csv',
                  source_these='analysis_WA/format_data.R')

## ---- subset ----
######################################################################
# OPTIONAL SUBSET TO A SUBGROUP
######################################################################

if ('subset_before_run'%in%ls()) {
    if (u=='non-MSM') {
        dataf <- dataf[dataf[,'mode']!='MSM',]
    } else {
        dataf <- dataf[dataf[,svar]==u,]
    }
    cat('\nSubgroup:', svar, '=', u, '\n')
    cat('\n       sample size =', nrow(dataf))
}

## ---- qtrCounts ----
######################################################################
# DEFINE DIAGNOSED COUNTS PER QUARTER
######################################################################

# This step is important because just using obsCounts, a tabulation
# of timeDx, will not include quarters with 0 counts. Those are
# effectively ignored, incorrectly treating the surrounding quarters
# as if they are adjacent
allTimes <- seq(time_min, time_max, by=0.25)
obsCounts <- table(dataf$timeDx)
allCounts <- structure(rep(0,length(allTimes)),
                       class='table',
                       names=allTimes)
allCounts[names(allCounts)%in%names(obsCounts)] <- obsCounts

######################################################################
# RUN BACK-CALCULATION
######################################################################

# Flag for working
allruns_worked <- FALSE

# Run both base case and upper bound for no imputing then imputing

## ---- run_noimpute ----
cat('\nBackcalculating with no/minimal imputation...\n')
all_noimpute <- runBackCalc(TID=dataf$infPeriod, 
                   impute=FALSE,
                   age=dataf$hdx_age,
                   diagnosedCounts=allCounts,
                   upperBound=FALSE, 
                   runBoth=TRUE,
                   intervalLength=0.25, 
                   printProgress=FALSE) 

## ---- run_impute ----
cat('\nBackcalculating with full imputation...\n')
all_impute <- runBackCalc(TID=dataf$infPeriod, 
                   impute=TRUE,
                   age=dataf$hdx_age,
                   diagnosedCounts=allCounts,
                   upperBound=FALSE, 
                   runBoth=TRUE,
                   intervalLength=0.25, 
                   printProgress=FALSE) 

# flag for working
allruns_worked <- TRUE

## ---- summarize ----
######################################################################
# SUMMARIZE
######################################################################

# Stats
stats_worked <- FALSE

# Plots
summaries_noimpute <- summarize_runBackCalc(results=all_noimpute,
                                   diagnosedCounts=allCounts,
                                   times=allTimes)

summaries_impute <- summarize_runBackCalc(results=all_impute,
                                   diagnosedCounts=allCounts,
                                   times=allTimes)

summaries_both <- summarize_runBackCalc_combined(
                            results=list(noimpute=all_noimpute, 
                                         impute=all_impute), 
                            diagnosedCounts=allCounts,
                            times=allTimes)

stats = data.frame(imputed=c(rep('Yes',nrow(summaries_impute[['stats']])),
                             rep('No',nrow(summaries_noimpute[['stats']]))),
                   rbind(summaries_impute[['stats']],
                         summaries_noimpute[['stats']]))
statsAll = data.frame(imputed=c(rep('Yes',
                                    nrow(summaries_impute[['statsAll']])), 
                                rep('No',
                                 nrow(summaries_noimpute[['statsAll']])),
                                rep('Yes',
                                    nrow(summaries_impute[['statsUndiag']])), 
                                rep('No',
                                 nrow(summaries_noimpute[['statsUndiag']]))),
                   rbind(summaries_impute[['statsAll']],
                         summaries_noimpute[['statsAll']], 
                         summaries_impute[['statsUndiag']],
                         summaries_noimpute[['statsUndiag']]))

stats_worked <- TRUE


## ---- tidy ----
######################################################################
# SAVE RESULTS
######################################################################


## ---- save ----
if (!'subset_before_run'%in%ls()) {
    write.csv(stats, file=file.path(workd, 'analysis_WA', 'results', 'run_main_summary.csv'), row.names=FALSE)
    write.csv(statsAll, file=file.path(workd, 'analysis_WA', 'results', 'run_main.csv'), row.names=FALSE)
}

