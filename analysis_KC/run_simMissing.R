
########################
# run_simMissing.R
# 
# Purpose: using the run structure from run_main.R, run the KC MSM analysis
# for different levels of missing data and summarize results
# 
# 
# Dependencies:
# data-cleaning_JKB.R
# packages HIVBackCalc, reshape2, ggplot2, scales, Hmisc, plyr
#
# History: based on run_main.R
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
# CREATE A 100% NONMISSING infPeriod WITH YEARLY DIAGNOSES
# EQUALING THOSE IN THE ORIGINAL MSM DATASET
######################################################################

set.seed(98103)
# Create a dataset that has no missing infPeriod but has the 
# same number of diagnoses per year as the original msm, 1522
sampleN <- table(msm$yearDx)
theserows <- lapply(2006:2012,
                    function(t) {
                        samplefrom=rownames(msm)[msm$yearDx==t & 
                                                 !is.na(msm$infPeriod)]
                        rows=sample(samplefrom, 
                                    size=sampleN[which(names(sampleN)==as.character(t))], 
                                    replace=TRUE)
                        return(rows)
                    })
rows <- unlist(theserows)
msmNoMiss <- msm[rows,]

# Check result
sum(is.na(msmNoMiss$infPeriod))
table(msmNoMiss$yearDx)
table(msm$yearDx)

######################################################################
# CREATE A DATAFRAME CONTAINING infPeriod (TID) VECTORS WITH
# VARYING LEVELS OF MISSINGNESS
######################################################################

# Define levels of missing, in percent
missLev <- seq(0,100,by=10)

# Create matrix of infPeriods with those levels of missing
infPeriods <- sapply(missLev,
                     function(m) {
                         infPeriod <- msmNoMiss$infPeriod
                         changeToMiss <- sample.int(nrow(msm), 
                                                    size=(m/100)*nrow(msm))
                         infPeriod[changeToMiss] <- NA
                         cat('Percent missing is', sum(is.na(infPeriod))/length(infPeriod),'\n')
                         return(infPeriod)
                     })

# Check result
round(100*colSums(is.na(infPeriods))/nrow(msmNoMiss),2)


######################################################################
# RUN BACKCALC, LOOPING THROUGH COLUMNS OF infPeriods
######################################################################

results <- vector(mode='list', length=ncol(infPeriods))

for (i in 1:ncol(infPeriods)) {

    cat('\n\nLEVEL OF MISSINGNESS IS', missLev[i], "%\n")
    ######################################################################
    # RUN BACK-CALCULATION
    ######################################################################

    # Run both cases: first the base case, then the upper bound

    all <- runBackCalc(TID=infPeriods[,i],
                       impute=TRUE,
                       age=msm$hiv_age_yrs,
                       diagnosedCounts=table(msmNoMiss$timeDx), 
                       upperBound=FALSE, 
                       runBoth=TRUE,
                       intervalLength=0.25, 
                       printProgress=FALSE) 


    ######################################################################
    # SUMMARIZE runBackCalc
    ######################################################################

    summaries <- summarize_runBackCalc(results=all,
                                       diagnosedCounts=table(msmNoMiss$timeDx),
                                       times=seq(2006,2012.75,by=0.25))

    results[[i]] <- data.frame(missPerc=missLev[i],summaries[[1]])
}

######################################################################
# SUMMARIZE RESULTS
######################################################################

results.df <- do.call(rbind,results)

figp_plot = ggplot(results.df, aes(x=as.factor(`var`), y=`Median`, fill=as.factor(missPerc))) +
  geom_bar(position=position_dodge(), stat='identity') +
  geom_errorbar(aes(ymin=`Min.`, ymax=`Max.`), width=0.2, 
                position=position_dodge(0.9)) +
  scale_x_discrete(name="") + 
  labs(title="") +
  scale_fill_manual(name="Percent Missing",values=colors()[1:11]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file='analysis/results/run_simMissing.pdf', width=6, height=5)
figp_plot
dev.off()


