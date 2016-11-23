#' Log probability of an BED od-n given date since sero conversion in days 
logProbOdn <- function(odn, sero, aids=NULL, ...){
  if(!is.null(aids))
    probAids <- aids
  else
    probAids <- pweibull(sero / 365, shape = 2.5, scale = 1/0.086)
  sero <- pmin(800, sero) # Don't project out of validated range
  coef <- structure(c(-1.40644598916766, 0.458516848702799), .Names = c("(Intercept)", "lsero"))
  stdlow <- 0.06160604
  stdhigh <- 0.2168296
  stdaids <- .4 # set so that p(odn < 1 | sero >= 800) = .043 . The missclassification rate for AIDS patients
  lsero <- log(pmax(0, sero)+50)
  sqodn <- sqrt(odn)
  
  #standard deviation is stdlow below 4.3 and stdhigh above 4.6, and averaged in between
  std <- ifelse(lsero < 4.6 & lsero > 4.3, 
                stdlow* (1 - (lsero-4.3)/.3) + stdhigh * (lsero-4.3)/.3,
                ifelse(lsero <= 4.3, stdlow, stdhigh)
  )
  resid <- sqodn - coef[1] - coef[2]*lsero
  #dnorm(resid, mean=0, sd=std)
  ifelse(sero==800, 
         log(dnorm(resid, mean=0, sd=stdaids) * probAids + dnorm(resid, mean=0, sd=std) * (1-probAids)),
         dnorm(resid, mean=0, sd=std,log=TRUE)
  )
}

#' log probability of positive from a bed recency test
logProbBedRecent <- function(bedRecent, sero, aids=NULL, ...){
  if(length(bedRecent) == 1)
    bedRecent <- rep(bedRecent, length(sero))
  if(!is.null(aids))
    probAids <- aids
  else
    probAids <- pweibull(sero / 365, shape = 2.5, scale = 1/0.086)
  
  # based on figure 7 of Parekh et. al.
  lin <- 6.036492e+00 - pmin(500,sero)* 4.596835e-02 + pmin(500,sero)^2 * 4.444895e-05
  prob <- exp(lin) / (1 + exp(lin))
  
  # older hiv infections have a .0169 false positive rate
  prob <- pmax(.0169, prob)
  
  # aids patients have a .0438 false positive rate
  prob <- prob * (1 - probAids) + .0438 * probAids
  
  ifelse(bedRecent, log(prob), log(1-prob))
}


#' log probability of having AIDS at time of diagnosis based on wiebull
#' Assumes individuals become diagnosed once they progress to AIDS
#' 
#' If subject has aids at dx, the probability that they have progressed is
#' dweibull(time, 2.5, 1/0.086)
#' 
#' If the subject does not have aids, then the probability that they
#' have not progressed at any time prior to the current time is
#' 1- pweibull(time, 2.5, 1/0.086)
logProbAids <- function(aids, sero, ...){
  if(length(aids) == 1)
    aids <- rep(aids, length(sero))
  sero <- pmax(0, sero) / 365
  sh <- 2.5
  sc <- 1/0.086
  ifelse(aids, dweibull(sero,shape = sh,scale = sc, log=TRUE),
         pweibull(sero,shape = sh,scale = sc, log.p=TRUE, lower.tail = FALSE))
}



#' probability of having a CD-4 count given sero conversion interval
logProbCd4 <- function(cd4, sero, ...){
  
  #' From s. lodi private communication
  sero <- sero/365
  sigmaA <- 6.12947
  sigmaB <- 1.515038
  sigmaError <- 2.788371
  corAB <- -.441032
  intercept <- 23.81368
  slope <- -1.219816
  
  sqcd4 <- sqrt(cd4)
  
  # Mixed model Probability. Normal distribution truncated at 0.
  mean <- intercept + slope*sero
  sd <- sqrt(sigmaA^2 + sero^2 * sigmaB^2 + sigmaError^2 + 2 * sero * corAB * sigmaA * sigmaB)
  dnorm(sqcd4, mean=mean, sd=sd, log=TRUE) - pnorm(mean/sd,log=TRUE)
}




