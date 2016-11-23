source("log-probs.R")
library(ggplot2)
dat <- read.csv("odn-n\ Sero.csv")
names(dat) <- c("sero","odn","a","b")

qplot(pmax(0,sero),odn,data=dat) + geom_smooth()
qplot(log(pmax(0,sero)+50),sqrt(odn),data=dat) + geom_smooth()

#' transform to linear relationship
lsero <- log(pmax(0,dat$sero)+50)
sqodn <- sqrt(dat$odn)

#' linear model for relationship
mod <- lm(sqodn ~ lsero)
summary(mod)
resid <- mod$residuals
qplot(lsero,abs(resid))

#' cut point for variance estimation
#' variances differ below and above cut point
lcut <- lsero > 4.5
qplot(resid[lcut])
qplot(resid[!lcut])

stdlow <- sd(resid[!lcut])
stdhigh <- sd(resid[lcut])

#' model for the probability positive given sero
pdat <- read.csv("bed_figure7.csv")
names(pdat) <- c("sero","prob_pos")
logitp <-log(pdat$prob_pos / (1-pdat$prob_pos))
sero <- pdat$sero
mod1 <- lm(logitp ~ sero + I(sero^2))
mod1$coefficients

#' 5 year descrimination of bed recent
plot(exp(logProbBedRecent(TRUE, (0:(365*5)))),type='l')
#' 20 year descrimination of bed recent
plot(exp(logProbBedRecent(TRUE, (0:(365*20)))),type='l')

#' plots showing bed od-n test discrimination at od-n = .5, .75, 1, 1.5, 2, 3
plot(exp(logProbOdn(.5, (0:2000))),type='l')
plot(exp(logProbOdn(.75, (0:2000))),type='l')
plot(exp(logProbOdn(1, (0:2000))),type='l')
plot(exp(logProbOdn(1.5, (0:2000))),type='l')
plot(exp(logProbOdn(2, (0:2000))),type='l')
plot(exp(logProbOdn(3, (0:2000))),type='l')



#' 20 year discrimination if aids=FALSE
plot(exp(logProbAids(FALSE, (0:(365*20)))),type='l')
#' 20 year discrimination if aids=TRUE
plot(exp(logProbAids(TRUE, (0:(365*20)))),type='l')



#' 20 year discrimination for cd4 count at cd4=0, 200, 350, 500, 700
plot(exp(logProbCd4(0,0:(365*20))),type='l')
plot(exp(logProbCd4(200,0:(365*20))),type='l')
plot(exp(logProbCd4(350,0:(365*20))),type='l')
plot(exp(logProbCd4(500,0:(365*20))),type='l')
plot(exp(logProbCd4(700,0:(365*20))),type='l')


#' #Probability distributions for different cases
#' Computes the p(tid | odn, aids) using naive bayes
naiveBayes <- function(odn, aids, lastNegHIV){
  times <- 0:(365*20)
  nt <- length(times)
  logProb <- rep(0,nt)
  
  logProb[times > lastNegHIV] <- -Inf
  logProb <- logProb + logProbAids(aids, times)
  logProb <- logProb + logProbOdn(odn, times)
  prob <- exp(logProb - max(logProb))
  prob <- prob / sum(prob)
  prob
}
#' non-aids, never tested, lowish odn. Note the hgher concentration toward
#' recent infection
plot(naiveBayes(odn=1.5, aids=FALSE, lastNegHIV=365*20),type='l')

#' last negative 2 years ago, comes in with aids and high odn. note the high concentration
#' toward the end of the 2 year interval
plot(naiveBayes(odn=3, aids=TRUE, lastNegHIV=365*2),type='l')


#' very small odn, but with aids and never tested. Possible misclassification of aids
#' due to primary infection. note the high concentration on recent dates
plot(naiveBayes(odn=.25, aids=TRUE, lastNegHIV=365*20),type='l')


#' #Calculation of TID distribution across a sample
tidDat <- data.frame(odn=c(1.5,2,.25, 2, 1, 3, 1), aids=c(F,T,T,F, F, F, F), ln=365*c(20,2,20,4, 10, 8, 1))
individDists <- apply(tidDat, 1, function(x) naiveBayes(x[1],x[2],x[3]))
tid <- rowMeans(individDists)
plot(tid,type='l')



