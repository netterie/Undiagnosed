source("data-cleaning.R")
library(ggplot2)
library(reshape2)
library(MASS)

by(msm,msm$yearDx,function(x)mean(x$infPeriod,na.rm=T))
oneway.test(infPeriod~yearDx,data=msm)

#blacks and asians have less frequent tests
by(msm,msm$racel,function(x)mean(x$infPeriod,na.rm=T))
table(msm$racel)
oneway.test(infPeriod~racel,data=msm)
ggplot(aes(y=infPeriod,x=factor(racel)),data=msm) + 
  geom_jitter(alpha=.25) + 
  geom_boxplot(color="darkred",fill=NA,outlier.size=0)


by(msm,paste(msm$racel,msm$yearDx),
   FUN=function(x)c(mean(x$infPeriod,na.rm=T),sum(!is.na(x$infPeriod)))
)


#blacks + asians more likely to never have been tested
x <- xtabs(~ racel + everTested,data=msm)
x / rowSums(x)
chisq.test(x,simulate=T)


#exponential fit
lastTest <- na.omit(msm$infPeriod)
dst <- dexp
fit <- fitdistr(lastTest, "exponential")
xFit <- seq(from = 0, to = max(lastTest), length.out = 500)
yFit <- do.call(dst, c(as.list(fit$estimate), list(x = xFit)))
ggplot() + geom_density(aes(x = lastTest)) + geom_line(aes(x = xFit, y = yFit), 
                                                       color = "red")


#exponential fit (non-aids at dx)
lastTest <- na.omit(msm$infPeriod[!msm$aidsAtDx])
dst <- dexp
fit <- fitdistr(lastTest, "exponential")
xFit <- seq(from = 0, to = max(lastTest), length.out = 500)
yFit <- do.call(dst, c(as.list(fit$estimate), list(x = xFit)))
ggplot() + geom_density(aes(x = lastTest)) + geom_line(aes(x = xFit, y = yFit), 
                                                       color = "red")




#
# Upper bounds on AIDS incubation periods evaluated using those who are diagnosed with AIDS at HIV
# diagnosis
#
msm1 <- msm[msm$aidsAtDx,]
msm1$infPeriodHe <- pmin(msm1$lagHIV_LastNegEhars,msm1$lagHIV_HISNeg,na.rm=TRUE) / 365
msm1$infPeriodHe <- with(msm1,ifelse(is.na(ungtst) & ungtst=="N",hiv_age_yrs-15,infPeriodHe))
msm1$infPeriodUB <- with(msm1,ifelse(is.na(infPeriod),hiv_age_yrs-15,infPeriodHe))
mean(msm1$hiv_age_yrs)
msm1 <- msm1[c("infPeriodHe","infPeriodUB","lagHIV_LastNegEhars","lagHIV_LastNegPCRSrpt","lagHIV_HISNeg")]
msm1[3:5] <- msm1[c(3:5)]/365
plot(msm1,xlim=c(0,30),ylim=c(0,30))
cor(msm1,use="pair")

nrow(msm1)
summary(msm1)
m <- melt(msm1[-c(2,4)])
ggplot(aes(x=value,color=variable),data=m) + 
  stat_function(fun = dweibull, colour = "grey50", arg = list(shape=2.516,scale=1/0.086), size=2) + 
  annotate("text", x = 17, y = .092, label = "Wiebull(shape=2.516,scale=11.628)",color="grey50") +
  geom_density(adjust=2.25,na.rm=TRUE) + 
  xlim(c(0,30)) + 
  theme_bw() + 
  scale_color_discrete("Report",labels=c("Combined","EHARS","HIS")) + 
  xlab("Years since last negative HIV test") + 
  ylab("Density")

s <- 1:300 / 10
a <- data.frame(a=ecdf(msm1[[1]])(s),b=ecdf(msm1[[3]])(s),
                c=ecdf(msm1[[4]])(s))
m1 <- melt(a)
ggplot(aes(x=rep(s,3),y=value,color=variable),data=m1) + 
  stat_function(fun = pweibull, colour = "grey50", arg = list(shape=2.516,scale=1/0.086), size=2) + 
  annotate("text", x = 17, y = .092, label = "Wiebull(shape=2.516,scale=11.628)",color="grey50") +
  geom_step() + 
  xlim(c(0,30)) + 
  theme_bw() + 
  scale_color_discrete("Report",labels=c("Combined","EHARS","HIS")) + 
  xlab("Years since last negative HIV test") + 
  ylab("Proportion without AIDS")

medConf <- function(x){
  x <- na.omit(x)
  ul <- sort(x)[qbinom(c(.025,.975), length(x), 0.5)]
  c(median = median(x),lower=ul[1],upper=ul[2])
}
sapply(msm1,medConf)


#
# Use pcrs and create hard upper bound assuming all missings have never been tested before.
#
msm1 <- msm[msm$aidsAtDx & msm$timeDx>=2010,]
msm1$infPeriodHe <- msm1$infPeriod
msm1$infPeriodUB <- with(msm1,ifelse(is.na(infPeriod),hiv_age_yrs-15,infPeriodHe))
mean(msm1$hiv_age_yrs)
msm1 <- msm1[c("infPeriodHe","infPeriodUB","lagHIV_LastNegEhars","lagHIV_LastNegPCRSrpt","lagHIV_HISNeg")]
msm1[3:5] <- msm1[c(3:5)]/365
plot(msm1,xlim=c(0,30),ylim=c(0,30))
cor(msm1,use="pair")

nrow(msm1)
summary(msm1)
m <- melt(msm1[-2])
ggplot(aes(x=value,color=variable),data=m) + 
  stat_function(fun = dweibull, colour = "grey50", arg = list(shape=2.516,scale=1/0.086), size=2) + 
  annotate("text", x = 17, y = .092, label = "Wiebull(shape=2.516,scale=11.628)",color="grey50") +
  geom_density(adjust=2.25,na.rm=TRUE) + 
  xlim(c(0,30)) + 
  theme_bw() + 
  scale_color_discrete("Report",labels=c("Combined","EHARS","PCRS","HIS")) + 
  xlab("Years since last negative HIV test") + 
  ylab("Density")

s <- 1:300 / 10
a <- data.frame(a=ecdf(msm1[[1]])(s),b=ecdf(msm1[[3]])(s),
                c=ecdf(msm1[[4]])(s),d=ecdf(msm1[[5]])(s))
m1 <- melt(a)
ggplot(aes(x=rep(s,4),y=value,color=variable),data=m1) + 
  stat_function(fun = pweibull, colour = "grey50", arg = list(shape=2.516,scale=1/0.086), size=2) + 
  annotate("text", x = 17, y = .15, label = "Wiebull(shape=2.516,scale=11.628)",color="grey50") +
  geom_step() + 
  xlim(c(0,30)) + 
  theme_bw() + 
  scale_color_discrete("Report",labels=c("Combined","EHARS","PCRS","HIS")) + 
  xlab("Years since last negative HIV test") + 
  ylab("Proportion without AIDS")

sapply(msm1,medConf)



