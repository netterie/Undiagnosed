########################
# run_JKB.R
# 
# Purpose: Generate results presented in the HIV backcalculation paper
# 
# 
# Dependencies:
# data-cleaning_JKB.R
# packages HIVBackCalc, reshape2, ggplot2, scales, Hmisc
#
# History: 
# This file is a commented version of Ian's run.R
########################

setwd('/Users/jeanette/Dropbox/School/PhD/HIV_WA')

# Setup
source("analysis/data-cleaning_JKB.R")
library(HIVBackCalc)
library(reshape2)
library(ggplot2)
library(scales)
library(Hmisc)

######################################################################
# PART 1: STATS BY DATA SOURCE (TABLE 1) 
######################################################################

#N
n <- dim(msm)[1]
n

# HIS - 1080
# Total n, less those with missing lagHIV_HISNeg OF those with everTested=NA or TRUE
# This is confusing, because you get the same thing with 
# n-sum(is.na(msm$lagHIV_HISNeg))
n-sum(is.na(msm[is.na(msm$everTested) || msm$everTested,]$lagHIV_HISNeg))

# EHARS - 382
# Again, equivalent to  n-sum(is.na(msm$lagHIV_LastNegEhars))
n-sum(is.na(msm[is.na(msm$everTested) || msm$everTested,]$lagHIV_LastNegEhars))

# PCRS - 476
# Again, equivalent to n-sum(is.na(msm$lagHIV_LastNegPCRSrpt))
n-sum(is.na(msm[is.na(msm$everTested) || msm$everTested,]$lagHIV_LastNegPCRSrpt))

# Correlations
cor(msm[c("lagHIV_HISNeg","lagHIV_LastNegEhars","lagHIV_LastNegPCRSrpt")],use="pair")

# Never tested - 101
# More easily seen using table(msm$everTested, useNA='ifany')
sum(!is.na(msm$everTested) & !msm$everTested)

# Combined - 1233 have a date of last negative test
table(is.na(msm$infPeriod))

######################################################################
# PART 2: TREND IN infPeriod
######################################################################

# Stability of time to dx
by(msm$infPeriod,msm$yearDx,function(a)mean(a,na.rm=T))
mean(msm$infPeriod,na.rm=T)
oneway.test(infPeriod ~ yearDx,data=msm)

######################################################################
# PART 3: STATS BY RACE, INCL TREND IN infPERIOD
######################################################################

# Race
table(msm$racel)
msm$race2 <- as.character(msm$racel)
msm$race2[msm$race2 %in% c("4PI","5AmInd","6Multi")] <- "4OtherMulti"
msm$race2 <- factor(msm$race2)
levels(msm$race2) <- c("Caucasian","African American","Hispanic","Asian","Other/Multi")
table(msm$race2)
by(msm,msm$race2,function(x)mean(x$infPeriod,na.rm=T))

# Violin plots - this is the part that needs Hmisc
pdf("plots/race.pdf",width=6,height=3)
p <- ggplot(aes(y=infPeriod,x=race2),data=msm) + 
  geom_violin(adjust=3,fill="grey83",color=NA) + 
  stat_summary(color="darkred",fun.data = "mean_cl_boot") +
  ylab("Years Since Last Negative") +
  xlab("") +
  theme_bw() +
  coord_flip()
plot(p)
dev.off()

# I think this code was abandoned but not deleted
race <- ""
#race <- "Caucasian"
#race <- "African American"
#race <- "Hispanic"
#race <- "Asian"
#race <- "Other"
if(race!=""){
  tmp <- race
  if(race == "Other"){
    tmp <- "Other/Multi"
  }
  msm <- msm[msm$race2==tmp,]
}

######################################################################
# FIGURE 1: CUMULATIVE TIME FROM LAST NEG TO INFECTION (1-EMPIRICAL CDF) 
######################################################################

# Sort infPeriod times - exclude zeroes, who never had a negative test
ti <- sort(msm$infPeriod[!is.na(msm$infPeriod) & msm$infPeriod > 0])
nUB <- length(ti) #1217
# Function to compute 1-cdf
# This reflects the UB assumption because it assumes that
# time of infection=time of last negative test
qUB <- function(u) {
  uind <- sum(ti<=u)/nUB
  if (is.na(uind)) 
    return(0)
  uind
}
# Function to compute at each infPeriod,
# the ratio of the # of diagnoses occurring at/after the time
# to total observed diagnoses...so, # implied diagnoses per 1 time unit
# What was the intention for eta?
pi <- function(i, eta, ti = ti) {
  sapply(i, function(ii) {
    # Get those infPeriod times that are >= ii, a unique infPeriod time
    ints <- ti[ti >= ii]
    # Average # diagnoses/unit time for times >= ii
    # Is this what assumes infection is uniformly distributed in an infPeriod?
    sum(1/ints)/length(ti)
  })
}
# Unique infPeriod times 
uti <- unique(ti) #652, of #1217 possible
# Multiply implied diagnoses/unit time by actual time step (vs 1 unit time)
# diff(c(0,uti)) - for each of 0 and the unique infPeriods, this generates a 1st order lag
p <- pi(uti, , ti) * diff(c(0, uti))
# Cumulative distribution of p
cs <- cumsum(p)
# Function to evaluate cs at selected times (u), 
# returning 0 for presumably only u=0 
qi <- function(u) { 
  # Indicators for the 1st time in uti (unique infPeriod times) >= u
  uind <- rev(which(uti <= u))[1]
  if (is.na(uind)) 
    return(0)
  # Return the cumulative distribution of p for that indicator
  cs[uind]
}
s <- seq(from=0,to=20,length.out=500)
# Base case 1-cdf
est <- sapply(s,qi)
# Upper bound 1-cdf
ub <- sapply(s,qUB)

# Plot
d1 <- rbind(data.frame(var="Base Case  ",value=est,Time=s),
            data.frame(var="Upper Bound  ",value=ub,Time=s))
pdf(paste0("plots/dist",race,".pdf"),width=5,height=3)
p <- ggplot(d1) + geom_line(aes(x=Time,y=1-value,color=var)) + 
  scale_color_hue(name="") +
  theme_bw() + 
  ylab("Undiagnosed Fraction") +
  xlab("Time Since Infection") +
  scale_x_continuous(expand=c(0,.2)) +
  theme(legend.position="bottom")
print(p)
dev.off()

########START HERE

######################################################################
# UNDIAGNOSED INCIDENCE (INFECTION) COUNTS, BASE CASE
######################################################################

# Set time step to quarter-years
intLength <- .25

# Define the discrete time probability distribution function of 
# time from infection to diagnosis, using the base case 
# assumption that infection is uniformly distributed between
# time of last negative test and time of diagnosis
pid <- estimateProbDist(infPeriod=msm$infPeriod,intLength=intLength)

# Set y = 100 NA's + number of diagnoses per quarter-year to indicate
# that we want to backcalculate incidence for 100 time steps prior to 
# our data
y <- c(rep(NA,100),table(msm$timeDx))

# Backcalculate quarterly incidence for all
mod <- estimateIncidence(y,pid,gamma=.1,verbose=TRUE,tol=10^-4)

# Plot estimated incidence counts over the period for which we observed 
# diagnoses, 2006-2012.75
#plot(mod,time=c(2006,2012.75))

# Estimate undiagnosed incidence given incidence, diagnosis, and pid function
undiag <- estimateUndiagnosed(mod)
obs <- !is.na(y)
time <- c(2006,2012.75)
time <- seq(from=time[1],to=time[2],length.out=sum(obs))
#plot(time,undiag[obs],ylim=c(300,400),type="l")

# Non-zero infPeriods
ti <- with(msm,sort(infPeriod[!is.na(infPeriod) & infPeriod>0]))

######################################################################
# NAIVE METHOD? SEE IAN'S analysis.Rnw/pdf
######################################################################

# Average diagnoses/year * average time to infection in years * # quarters
# Based on analysis.Rnw/pdf, I think he means 12/3 not 12/4
#mean(y,na.rm=TRUE) * mean(ti/2)*12/4
mean(y,na.rm=TRUE) * mean(ti/2)*12/3

######################################################################
# UNDIAGNOSED INCIDENCE (INFECTION) COUNTS, UPPER BOUND
######################################################################

# Alternative to estimateProbDist, the base case function,
# that uses the upper bound assumption
empirProbDist <- function(infPeriod,intLength=1){
  ti <- sort(infPeriod[!is.na(infPeriod) & infPeriod > 0])
  n <- length(ti)
  qi <- function(u) {
    uind <- sum(ti<=u)/n
    if (is.na(uind)) 
      return(0)
    uind
  }
  pidCalc <- function(i) {
    sapply(i, function(ii) {
      qi((ii + 1) * intLength) - qi(ii * intLength)
    })
  }
  m <- max(ti/intLength) + 1
  pidProbs <- pidCalc(0:m)
  pid <- function(i) {
    ifelse(i > m, 0, pidProbs[i + 1])
  }
  pid
}

# Set y = 100 NA's + number of diagnoses per quarter-year to indicate
# that we want to backcalculate incidence for 100 time steps prior to 
# our data
y <- c(rep(NA,100),table(msm$timeDx))

# Define the discrete time probability distribution function of 
# time from infection to diagnosis, using the upper bound
# assumption that infection is occurs at the
# time of last negative test and time of diagnosis
pidUB <- empirProbDist(infPeriod=msm$infPeriod,intLength=intLength)

# Backcalculate quarterly incidence for all
modUB <- estimateIncidence(y,pidUB,gamma=.1,verbose=TRUE,tol=10^-4)
#plot(modUB,time=c(2006,2012.75))

# Estimate undiagnosed incidence given incidence, diagnosis, and pid function
undiagUB <- estimateUndiagnosed(modUB)
obs <- !is.na(y)
time <- c(2006,2012.75)
time <- seq(from=time[1],to=time[2],length.out=sum(obs))
#plot(time,undiagUB[obs],ylim=c(300,800),type="l")


######################################################################
# RESULTS SUMMARIES
######################################################################

d <- rbind(
  data.frame(time=time,var="# Diagnosed  ",value=y[obs]),
  data.frame(time=time,var="Incidence (Base Case)  ",value=mod$lambda[obs]),
  data.frame(time=time,var="Incidence (Upper Bound)  ",value=modUB$lambda[obs])
)
by(d$value,d$var,summary)

pdf(paste0("plots/inc",race,".pdf"),height=5)
p <- ggplot(d,aes(x=time,y=value,linetype=var))  +   
  geom_line(aes(alpha=var)) +
  geom_point(aes(color=var)) + 
  theme_bw() + 
  scale_alpha_manual(values=c(.5,1,1),name="") + 
  scale_color_hue(name="") + 
  scale_linetype_manual(name="",values=c(3,1,2)) + 
  xlab("Time") + ylab("Counts") + ylim(c(0,75)) +
  theme(legend.position="bottom")
print(p)
dev.off()

d1 <- rbind(
  data.frame(time=time,var="Base Case  ",value=undiag[obs]),
  data.frame(time=time,var="Upper Bound  ",value=undiagUB[obs])
)
cols <- hue_pal()(3)[-1]
pdf(paste("plots/undiag",race,".pdf"),height=5)
p <- ggplot(d1,aes(x=time,y=value,linetype=var))  +   
  geom_line() +
  geom_point(aes(color=var)) + 
  theme_bw() + 
  scale_color_manual(name="",values=cols) + 
  scale_linetype(name="") + 
  xlab("Time") + ylab("# Undiagnosed HIV+") + ylim(c(0,750)) +
  theme(legend.position="bottom")
print(p)
dev.off()

#Count ranges
by(d1$value,d1$var,summary)


######################################################################
# UNDIAGNOSED INCIDENCE (INFECTION) COUNTS, CONSTANT INCIDENCE
######################################################################

# Average diagnoses/year = average incidence
incidence <- mean(y[obs])

# Base case
k <- 10000
m <- max(ti)
s <- seq(from=0,to=m,length.out=k)
l <- length(ti)
v <- sum(sapply(s,function(x) (sum(1-qi(x)))))
v * 4 * incidence / (k/m)

# Upper bound
v <- sum(sapply(s,function(x) (sum(1-qUB(x)))))
v * 4 * incidence / (k/m)



