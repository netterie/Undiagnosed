source("data-cleaning.R")
library(HIVBackCalc)
library(reshape2)
library(ggplot2)
library(scales)
#N
n <- dim(msm)[1]
n

#HIS
n-sum(is.na(msm[is.na(msm$everTested) || msm$everTested,]$lagHIV_HISNeg))

#EHARS
n-sum(is.na(msm[is.na(msm$everTested) || msm$everTested,]$lagHIV_LastNegEhars))


#PCRS
n-sum(is.na(msm[is.na(msm$everTested) || msm$everTested,]$lagHIV_LastNegPCRSrpt))

#cors
cor(msm[c("lagHIV_HISNeg","lagHIV_LastNegEhars","lagHIV_LastNegPCRSrpt")],use="pair")

#never tested
sum(!is.na(msm$everTested) & !msm$everTested)

#combined
table(is.na(msm$infPeriod))

#stability of time to dx
by(msm$infPeriod,msm$yearDx,function(a)mean(a,na.rm=T))
mean(msm$infPeriod,na.rm=T)
oneway.test(infPeriod ~ yearDx,data=msm)

#Race
table(msm$racel)
msm$race2 <- as.character(msm$racel)
msm$race2[msm$race2 %in% c("4PI","5AmInd","6Multi")] <- "4OtherMulti"
msm$race2 <- factor(msm$race2)
levels(msm$race2) <- c("Caucasian","African American","Hispanic","Asian","Other/Multi")
table(msm$race2)
by(msm,msm$race2,function(x)mean(x$infPeriod,na.rm=T))

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

#Plot distributions
ti <- sort(msm$infPeriod[!is.na(msm$infPeriod) & msm$infPeriod > 0])
nUB <- length(ti)
qUB <- function(u) {
  uind <- sum(ti<=u)/nUB
  if (is.na(uind)) 
    return(0)
  uind
}
pi <- function(i, eta, ti = ti) {
  sapply(i, function(ii) {
    ints <- ti[ti >= ii]
    sum(1/ints)/length(ti)
  })
}
uti <- unique(ti)
p <- pi(uti, , ti) * diff(c(0, uti))
cs <- cumsum(p)
qi <- function(u) {
  uind <- rev(which(uti <= u))[1]
  if (is.na(uind)) 
    return(0)
  cs[uind]
}
s <- seq(from=0,to=20,length.out=500)
est <- sapply(s,qi)
ub <- sapply(s,qUB)
d1 <- rbind(data.frame(var="Base Case  ",value=est,Time=s),
            data.frame(var="Upper Bound  ",value=ub,Time=s)
)
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



#estimated
intLength <- .25
y <- c(rep(NA,100),table(msm$timeDx))
pid <- estimateProbDist(infPeriod=msm$infPeriod,intLength=intLength)
mod <- estimateIncidence(y,pid,gamma=.1,verbose=TRUE,tol=10^-4)
#plot(mod,time=c(2006,2012.75))

undiag <- estimateUndiagnosed(mod)
obs <- !is.na(y)
time <- c(2006,2012.75)
time <- seq(from=time[1],to=time[2],length.out=sum(obs))
#plot(time,undiag[obs],ylim=c(300,400),type="l")

ti <- with(msm,sort(infPeriod[!is.na(infPeriod) & infPeriod>0]))
mean(y,na.rm=TRUE) * mean(ti/2)*12/4


#upper bound
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
y <- c(rep(NA,100),table(msm$timeDx))
pidUB <- empirProbDist(infPeriod=msm$infPeriod,intLength=intLength)
modUB <- estimateIncidence(y,pidUB,gamma=.1,verbose=TRUE,tol=10^-4)
#plot(modUB,time=c(2006,2012.75))

undiagUB <- estimateUndiagnosed(modUB)
obs <- !is.na(y)
time <- c(2006,2012.75)
time <- seq(from=time[1],to=time[2],length.out=sum(obs))
#plot(time,undiagUB[obs],ylim=c(300,800),type="l")

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


#constant incidence models
incidence <- mean(y[obs])

k <- 10000
m <- max(ti)
s <- seq(from=0,to=m,length.out=k)
l <- length(ti)
v <- sum(sapply(s,function(x) (sum(1-qi(x)))))
v * 4 * incidence / (k/m)


v <- sum(sapply(s,function(x) (sum(1-qUB(x)))))
v * 4 * incidence / (k/m)



