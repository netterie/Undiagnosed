source("log-probs.R")
source("data-cleaning.R")
#library(mi)

#*********************************************
# WA data
source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/internal_fxns.R')

# Change year min and max
year_min <- 2005
year_max <- 2014

# Load libraries, data and data-cleaning file
setup_hivbackcalc(workd='/Users/jeanette/Dropbox/School/PhD/HIV_WA',
                  datafile='data/wa_backcalc_data_201602.csv',
                  source_these='analysis_WA/format_data.R',
                  package_updated=TRUE,
                  packagefile='HIVBackCalc/R/internal_fxns.R')

dat <- dataf[c('firstcd4cnt', 'cd4_days', 'infPeriod')]
#*********************************************



sero <- 0:(365*20 - 1)
m <- length(sero)

#*********************************************
#dat <- msm[c("CD4result","CD4days","infPeriod")]
#*********************************************
dat$infPeriod <- dat$infPeriod * 365
names(dat) <- c("cd4","cd4lag","interval")
dat$cd4[dat$cd4lag > 60] <- NA

#' analysis omitting missing
dat1 <- na.omit(dat)
nrow(dat)
nrow(dat1)
n <- nrow(dat1)

probDist <- rep(0, m)
probDistUnif <- rep(0, m) 
for(i in 1:n){
  row <- dat1[i,,drop=FALSE]
  logProb <- rep(0,m)
  logProb[sero > row$interval] <- -Inf
  probDistUnif <- probDistUnif + exp(logProb) / sum(exp(logProb))
  logProb <- logProb + logProbCd4(row$cd4, sero + row$cd4lag)
  prob <- exp(logProb - max(logProb))
  prob <- prob / sum(prob)
  probDist <- probDist + prob
}
probDist <- probDist / n
probDistUnif <- probDistUnif / n

plot(sero/365, 
     1 - cumsum(probDistUnif),
     type='l',
     ylab="Undiagnosed Fraction",
     xlab="Time Since Infection",
     main="Base case: black , cd4 corrected: red")
points(sero/365, 1 - cumsum(probDist),type='l',col='red')


#' expected tid
sum(sero * probDist)

#'expected tid base case
sum(sero * probDistUnif)
