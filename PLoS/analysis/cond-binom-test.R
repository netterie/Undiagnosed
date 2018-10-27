library(sampling)

p <- c(.1,.05,.2)
n <- c(10,50,30)
z <- 10
a <- rbinom(1000000,size=n[1],p=p[1])
b <- rbinom(1000000,size=n[2],p=p[2])
c <- rbinom(1000000,size=n[3],p=p[3])
d <- cbind(a,b,c)
s <- a+b+c
d1 <- d[s==z,]
q <- p/(1-p)
pr <- q/sum(q)
pr*n


cn <-matrix(NA,nrow=10000,ncol=3)
# for(ind in 1:10000){
#   q <- p/(1-p)
#   vals <- c(rep(q[1],n[1]),
#             rep(q[2],n[2]),
#             rep(q[3],n[3]))
#   #i <- sampford(vals,z)
#   #i <- pps.sampling(vals,z,method="sampford")$sample
#   i <- as.logical(UPMEsfromq(UPMEqfromw(vals,z)))
#   cn[ind,] <- c(sum(vals[i]==q[1]),
#                 sum(vals[i]==q[2]),
#                 sum(vals[i]==q[3]))
#   
# }

for(ind in 1:5000){
  cn[ind,] <- genConditionalBinom(n,p,z)
}
colMeans(cn)
colMeans(d1)

round(prop.table(table(cn[,2])),4)
round(prop.table(table(d1[,2])),4)

pr <- table(d1[,1])
format(r/sum(r),sci=F)
format(dbinom(0:z,size=z,p=p[1]),sci=F)
