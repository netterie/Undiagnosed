########################
# describe_data.R
# 
# Purpose: 
# Formal descriptive statistics for the analysis
# sample determined by format_data.R and described
# in format_data_report.pdf
#
# Dependencies:
# format_data.R
# packages ggplot2, reshape2, MASS, survival, gridExtra
#
# History: 
# Based off of Ian's descriptives.R
#
# Notes:
# Meant to be called by other files that have loaded the data
# Uses read_chunks() notation (## ---- tag ----) so that
# it can be used in reports as well as sourced by other
# R files. See describe_data_report.Rnw
########################

## ---- diagnoses ----
#############################################################
# DIAGNOSIS COUNTS PER QUARTER
#############################################################

qtrDx <- table(dataf$timeDx)
qtrDx
sum(qtrDx)
qtrs <- as.numeric(names(qtrDx))

plot(qtrs, as.numeric(qtrDx), cex=(qtrDx/100),
     xlab='Quarter of diagnosis',
     ylab='Number of diagnosied cases',
     main='Diagnoses per quarter')

plot_diag <- plot_qtrDx(dataf)
plot_diag_subgroup <- plot_qtrDx(dataf,c(Race='race',Mode='mode'))

## ---- age ----
#############################################################
# AGE AT DIAGNOSIS
#############################################################

hist(dataf$hdx_age, main='Age at diagnosis')

## ---- everHadNegTest_subgroups ----
#############################################################
# everHadNegTest
#############################################################

variables <- c(`Age Group`='agecat5', 
               `Race/Ethnicity`='race', 
               `Mode of Transmission`='mode')

(everHadNegTest_subgrouptab <- tabulate_everHadNegTest(dataf,
                                                      variables,
                                                      supercolumn=TRUE))

(everHadNegTest_racebydx <- tabulate_everHadNegTest(dataf,
                                                    list(c('mode', 'race'))))

## ---- everHadNegTest ----
#############################################################
# everHadNegTest
#############################################################

summarize_infPeriod(dataf$infPeriod,
                    bygroup=dataf$everHadNegTest)

## ---- time ----
#############################################################
# TIME
#############################################################

summarize_infPeriod(dataf$infPeriod,
                    bygroup=dataf$yearDx)

ggplot(aes(y=infPeriod,x=factor(yearDx)),data=dataf) + 
  geom_jitter(alpha=.25) + 
  geom_boxplot(color="darkred",fill=NA,outlier.size=0)

## ---- time_everHadNegTest ----
(everHadNegTest_time <- tabulate_everHadNegTest(dataf,'yearDx'))
                                                    
(everHadNegTest_modebytime <- tabulate_everHadNegTest(dataf,
                                                    list(c('mode', 'yearDx'))))
(everHadNegTest_mode2bytime <- tabulate_everHadNegTest(dataf,
                                                    list(c('mode2', 'yearDx'))))
(everHadNegTest_racebytime <- tabulate_everHadNegTest(dataf,
                                                    list(c('race', 'yearDx'))))
plot_time <- plot_everHadNegTest(everHadNegTest_time)
plot_racetime <- plot_everHadNegTest(everHadNegTest_racebytime, 
                                     panel='race')
plot_modetime <- plot_everHadNegTest(everHadNegTest_modebytime, 
                                     panel='mode')
plot_mode2time <- plot_everHadNegTest(everHadNegTest_mode2bytime, 
                                     panel='mode2')
  
## ---- race ----
#############################################################
# RACE
#############################################################

summarize_infPeriod(dataf$infPeriod,
                    bygroup=dataf$new_race)

ggplot(aes(y=infPeriod,x=factor(new_race)),data=dataf) + 
  geom_jitter(alpha=.25) + 
  geom_boxplot(color="darkred",fill=NA,outlier.size=0)

## ---- mode2 ----
#############################################################
# MODE2
#############################################################

summarize_infPeriod(dataf$infPeriod,
                    bygroup=dataf$mode2)

# MSM vs non-MSM counts and % of all diagnoses
table(dataf$mode2)
table(dataf$mode2)/nrow(dataf)
# Missing TID (is.na(infPeriod)) by mode2
table(dataf$mode2, is.na(dataf$infPeriod))
(missTIDmode2 <- table(dataf$mode2, is.na(dataf$infPeriod)))
# % of mode2 cases comprising non-missing TID (is.na=FALSE) and missing (is.na=TRUE) 
missTIDmode2/colSums(missTIDmode2)

ggplot(aes(y=infPeriod,x=factor(mode2)),data=dataf) + 
  geom_jitter(alpha=.25) + 
  geom_boxplot(color="darkred",fill=NA,outlier.size=0)

## ---- density ----
#############################################################
# DENSITY, INFECTION OCCURS AT LAST NEG TEST (UPPER BOUND)
#############################################################

kmdat <- subset(dataf, !is.na(infPeriod))
par(mfrow=c(2,2))
hist(kmdat$infPeriod, freq=TRUE, breaks=100,
     xlab='Time from last negative test to diagnosis',
     main='Histogram and density of TID\nExcluding everHadNegTest=NA')
plot(density(kmdat$infPeriod), main='Density\nExcluding everHadNegTest=NA')
hist(dataf$infPeriod_imputeNA, freq=TRUE, breaks=100,
     xlab='Time from last negative test to diagnosis',
     main='Histogram and density of TID\nImputing everHadNegTest=NA')
plot(density(dataf$infPeriod_imputeNA, na.rm=TRUE),
     main='Density\nImputing everHadNegTest=NA')

## ---- km ----
#############################################################
# KAPLAN-MEIER, INFECTION OCCURS AT LAST NEG TEST (UPPER BOUND)
#############################################################
library(survival)

kmdat <- transform(kmdat, status=1)
dataf <- transform(dataf, status=1)
kmfit <- survfit(Surv(infPeriod,status)~1, data=kmdat)
kmfit2 <- survfit(Surv(infPeriod_imputeNA,status)~1, data=dataf)


## ---- kmplot ----
par(mfrow=c(1,2))
plot(kmfit, xlab='Time from last negative test to diagnosis\nexcluding everHadNegTest=NA', ylab='Proportion undiagnosed', main='TID excluding everHadNegTest=NA\nKM of upper bound')
plot(kmfit2, xlab='Time from last negative test to diagnosis\nImputing everHadNegTest=NA', ylab='Proportion undiagnosed', main='TID imputing everHadNegTest=NA\nKM of upper bound')

## ---- smoothsurv ----
#############################################################
# SMOOTHED SURVIVAL CURVE, BASE AND UPPER BOUND
#############################################################

library(gridExtra)
fig1plot <- fig1(dataf$infPeriod) + ggtitle('excluding everHadNegTest=NA')
fig2plot_imputeNA <- fig1(dataf$infPeriod_imputeNA) + 
    ggtitle('imputing everHadNegTest=NA')
grid.arrange(fig1plot, fig2plot_imputeNA, nrow=1)

## ---- smoothsurv_combined ----
#############################################################
# SMOOTHED SURVIVAL CURVE, BASE AND WORST CASE (OBS) AND
# WORST CASE (MISS)
#############################################################

newfig1 <- fig1combined(dataf, legendposition='right')

racefig1 <- fig1combined(dataf, panel='race')
modefig1 <- fig1combined(dataf, panel='mode')
