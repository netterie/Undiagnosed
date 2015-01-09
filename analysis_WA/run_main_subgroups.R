
########################
# run_main_subgroups.R
# 
# Purpose: run the analysis using the new wrapper R functions
# in HIVBackCalc/R/other, package version 1.1, FOR SUBGROUPS
# 
# 
# Dependencies:
# format_data.R
# packages HIVBackCalc, reshape2, ggplot2, scales, Hmisc, plyr
#
# History: based on run_JKB.R
#
########################

######################################################################
# SETUP and DATA CLEANING
######################################################################
rm(list=ls())
debug <- FALSE
# TEMPORARY: SOURCE FUNCTIONS
source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/other.R')

# Load libraries, data and data-cleaning file
# Eventually this function should return the cleaned dataset,
# but data-cleaning has the name hardcoded as msm and I'm not
# going to generalize that right now
setup_hivbackcalc(workd='/Users/jeanette/Dropbox/School/PhD/HIV_WA',
                  datafile='data/WA_BACKCALC_DATA_v2.csv',
                  source_these='analysis_WA/format_data.R')

if (debug) {
    warning('In debugging mode')
    source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/model.R')
}
warning('You probably want to set printSummaries=FALSE in format_data.R')
######################################################################
# DEFINE SUBGROUPS
######################################################################

# Signal to run_main.R that we will do a subgroup analysis
subset_before_run <- TRUE

# Define subgroup variables
subgroup_vars <- c('mode', 'race')


######################################################################
# RUN
######################################################################
warning('This does not save the quarterly estimates, only the summary over quarters')

# Save full sample
datafull <- dataf

# For storing results
statsfull <- list(NULL)
warningsu <- list(NULL)
fullres <- list(NULL)

for (svar in subgroup_vars) {
    cat('\nSUBGROUP VARIABLE BREAKDOWN:')
    print(table(datafull[,svar]))
    cat('\n\n')

    svar_levels <- unique(datafull[,svar])
    these_levels <- svar_levels

    if (is.factor(datafull[,svar])) {
        charlevels=as.character(svar_levels)
        if (svar=='mode') charlevels=c('non-MSM', charlevels)
        these_levels <- charlevels
    }

    for (u in these_levels) {

        # Run subgroup analysis
        tryCatch.W.E(source(file.path(workd, 'analysis_WA', 'run_main.R')))
        
        # Store results
        if (is.factor(datafull[,svar])) u <- charlevels[which(svar_levels==u)]
        if (allruns_worked & stats_worked) {
            pl <- summaries_both + ggtitle(u)
            fullres[[paste0(svar,'=',u)]] <- pl             
            usave <- gsub('/', '', u)
            ggsave(file=file.path(workd, 'analysis_WA', 'results', 
                                  paste0(usave,'.pdf')), plot=pl)
            stats <- transform(stats,
                               Subgroup=paste0(svar,'=',u))
            statsAll <- transform(statsAll,
                                  Subgroup=paste0(svar,'=',u))
            write.csv(stats, 
                      file=file.path(workd, 'analysis_WA', 'results', 
                                     paste0(usave,'summary.csv')), 
                      row.names=FALSE)
            write.csv(statsAll, 
                      file=file.path(workd, 'analysis_WA', 'results', 
                                     paste0(usave,'.csv')), 
                      row.names=FALSE)
        } else if (allruns_worked & !stats_worked) {
            fullres[[paste0(svar,'=',u)]] <- list(noimpute=all_noimpute,
                                                  impute=all_impute)
            stats <- NULL
        } else {
            fullres[[paste0(svar,'=',u)]] <- NULL
            stats <- NULL
        }
        statsfull[[paste0(svar,'=',u)]] <- stats
                                          
        # Reset dataf
        dataf <- datafull
    }
}

######################################################################
# FORMAT AND STORE RESULTS
######################################################################

isnull_statsfull  <- sapply(statsfull, is.null)
statsfull_f <- lapply(statsfull[!isnull_statsfull], format_stats)
statsfull_d <- do.call(rbind, statsfull_f)
statsfull_d  <- statsfull_d[,c(8,1:7)]
statsfull_d$Subgroup <- gsub('[a-z]+=','',statsfull_d$Subgroup)
colnames(statsfull_d)[2] <- 'Quantity'
write.csv(statsfull_d, 
          file=file.path(workd, 'analysis_WA', 'results', 'run_main_subgroups.csv'), 
          row.names=FALSE)


