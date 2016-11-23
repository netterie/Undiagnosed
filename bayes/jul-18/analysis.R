source("log-probs.R")
source("data-cleaning.R")
#library(mi)


sero <- 0:(365*20 - 1)
m <- length(sero)

dat <- msm[c("CD4result","CD4days","infPeriod")]
dat$infPeriod <- dat$infPeriod * 365
names(dat) <- c("cd4","cd4lag","interval")
dat$cd4[dat$cd4lag > 60] <- NA

#' analysis omitting missing
dat1 <- na.omit(dat)
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
