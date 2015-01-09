library(rootSolve)

#' Estimates the discrete time probability distribution function for
#' the time between infection and diagnosis
#' @param infPeriod a series of times from last HIV test to diagnosis
#' @param intLength the length of the discrete time intervals
estimateProbDist <- function(infPeriod,intLength=1){
  ti <- sort(infPeriod[!is.na(infPeriod) & infPeriod>0])
  #continuous density of time between infection and diagnosis
  pi <- function(i,eta,ti=ti){
    sapply(i,function(ii){
      ints <- ti[ti>=ii]
      sum(1/ints)/length(ti)
    })
  }
  uti <- unique(ti)
  p<-pi(uti,,ti) * diff(c(0,uti))
  cs <- cumsum(p)
  #cdf of density
  qi <- function(u){
    uind <- rev(which(uti<=u))[1]
    if(is.na(uind))
      return(0)
    cs[uind]
  }
  
  #the descrete probability that a diagnosis is made i time units after infection
  pidCalc <- function(i){
    sapply(i,function(ii){
      qi((ii+1)*intLength) - qi(ii*intLength)
    })
  }

  #calc discrete time prob
  m <- max(ti/intLength) + 1
  pidProbs <- pidCalc(0:m)
  pid <- function(i){
    ifelse(i>m,0,pidProbs[i+1])
  }
  pid
}


#' EM update step
#' @param y y
#' @param pid pid
#' @param lambda lambda
#' @param gamma gamma
meanEmUpdate <- function(y,pid,lambda,gamma=0){
  T <- length(y)
  obs <- !is.na(y)
  a <- b <- c <- lamNew <- rep(NA,T)
  for(k in 1:length(lambda)){
    s <- 0:(T-k)
    b[k] <- sum(pid(s))
    no <- !obs[s+k]
    if(any(no))
      a[k] <- sum(pid(s[no])) / b[k]
    else
      a[k] <- 0
    c[k] <- 0
    for(d in s){
      if(obs[k+d]){
        c[k] <- c[k] + y[k+d]*pid(d) / sum(lambda[1:(k+d)]*pid(k+d-(1:(k+d))))
      }
    }
  }
  if(gamma > 0){
    f <- function(ll){
      (1/ll) * (a*b+c) * lambda - b - 2 * gamma * c(0, ll[2:T] - ll[-T]) - 
        2 * gamma * c(ll[1:(T-1)] - ll[-1] ,0)
    }
    j <- function(ll){
      diag <- (-1/ll^2)* (a*b+c)*lambda - 4 * gamma
      diag[1] <- diag[1] + 2 * gamma
      diag[T] <- diag[T] + 2 * gamma
      off <- rep(0,T)
      off[2:(T-1)] <- 2*gamma
      rbind(-off,diag,off)
    }
    #browser()
    lamNew <- multiroot(f=f,start=lambda,positive=TRUE,jacfunc=j,jactype="bandint")$root
  }else{
    lamNew <- lambda * (a + c / b)
  }
  lamNew
}

#' estimates incedence via back calculation
#' @param y diagnosis counts
#' @param pid a function giving the probability of diagnosis given time units from infection
#' @param gamma smoothing penulty
#' @param tol The tolerance for the EM algorithm
#' @param verbose level of verbosity
estimateIncidence <- function(y,pid,gamma=0,tol=10^-5,verbose=FALSE){
  lambda <- rep(mean(y,na.rm=TRUE),length(y))
  ll <- lambda
  dev <- Inf
  while(dev>tol){
    lambda <- meanEmUpdate(y,pid,lambda,gamma)
    dev <- sum((ll-lambda)^2/ll)
    if(verbose){
      cat("lambda: ",paste(round(lambda,1),collapse=" "),"\n",
          "parameter change: ",dev,"\n",sep="")
    }
    ll <- lambda
  }
  mod <- list(lambda=lambda,y=y,pid=pid,gamma=gamma,tol=tol)
  class(mod) <- "backproj"
  mod
}

#' print function
#' @param x object of type backproj
#' @param ... passed to cat
print.backproj <- function(x,...) {
  cat("Back Projection Incidence Model\n","Estimated incidence: ",x$lambda,...)
}

#' plots the model
#' @param x the incedence model
#' @param time a vector of length two giving the start and end times of the observed period
#' @param showDiagCounts should diagnosis counts be plotted
#' @param ... passed to plot
plot.backproj <- function(x,time,showDiagCounts=TRUE, ...){
  obs <- !is.na(x$y)
  if(missing(time)) 
    time <- 1:length(x$y[obs])
  else
    time <- seq(from=time[1],to=time[2],length.out=sum(obs))
  plot(time,x$lambda[obs],ylim=c(0,100),type="l",main="Estimated Incidence",ylab="Count",...)
  if(showDiagCounts)
    points(time,x$y[obs],col="red")
}

#' Estimates the number of undiagnosed in the population
#' @param mod the model
#' @nExt the number of time units beyond the end of the observed period to use in calculation
estimateUndiagnosed <- function(mod,nExt=500){
  pid <- mod$pid
  lambda <- mod$lambda
  y <- mod$y
  nExt <- 500
  obs <- c(!is.na(y),rep(FALSE,nExt))
  T <- length(y)
  P <- T + nExt
  n <- matrix(NA,nrow=T,ncol=P)
  for(i in 1:T){
    for(j in i:P){
      if(!obs[j]){
        n[i,j] <- lambda[i]*pid(j-i)
      }else{
        n[i,j] <- y[j] * lambda[i]*pid(j-i) / sum( lambda[1:j]*pid(j-(1:j)) )
      }
    }
  }
  undiag <- rep(NA,T)
  for(i in 1:T) undiag[i] <- sum(n[1:(i),(i+1):P])
  undiag
}


#fit model with EM
# tb <- table(msm$timeDx)
# lambda <- rep(mean(y,na.rm=TRUE),length(y))
# ll <- lambda
# for(i in 1:400){
#   print(lambda <- meanEmUpdate(y,lambda,.1))
#   print(sum(ll-lambda)^2)
#   ll <- lambda
# }

#view incendece




# genConditionalBinom <- function(n,p,y){
#   q <- p/(1-p)
#   vals <- unlist(sapply(1:length(n),function(x) rep(p[x],n[x])))
#   i <- as.logical(UPMEsfromq(UPMEqfromw(vals,y)))
#   j<-1
#   cn <- cumsum(n)
#   res <- rep(NA,length(n))
#   for(k in 1:length(n)){
#     if(n[k]==0){
#       res[k] <- 0
#       next
#     }
#     res[k] <- sum(i[j:cn[k]])
#     j <- cn[k]+1
#   }
#   res
# }
# 
# genZ <- function(x,y){
#   z <- list()
#   p <- list()
#   n <- length(x)
#   tmp <- pid(0:n)
#   for(i in 1:n){
#     z[[i]] <- rep(NA,n+1-i)
#   }
#   p[[1]] <- tmp
#   avail <- x[1]
#   for(j in 1:n){
#     ps <- sapply(p,function(x)x[1])
#     if(is.na(y[j])){
#       zj <- rbinom(rep(1:length(ps)),size=avail,p=ps)
#     }else{
#       zj <- genConditionalBinom(n=avail,p=ps,y=y[j])
#     }
#     for(i in 1:length(ps)){
#       z[[i]][j] <- zj[i]
#       avail[i] <- avail[i] - zj[i]
#       p[[i]] <- p[[1]][-1]
#       p[[i]] <- p[[i]]/sum(p[[i]])
#     }
#     p[[j+1]] <- tmp
#     avail[j+1] <- x[j+1]
#   }
#   do.call(rbind,z)
# }
# 
# x <- rpois(20,40)
# y <- c(rep(NA,200),rpois(10,40))
# debug(genZ)
# z <- genZ(x,y)


# xform <- x ~ 1
# # poisson log probabilities for incendence
# lpx <- function(x, beta, xform=xform){
#   t <- 1:length(x)
#   df <- data.frame(t)
#   mm <- model.matrix(xform,df)
#   lambda <- exp(mm %*% beta) 
#   dpois(x,lambda,log=TRUE)
# }



# #log probability of z_{i,\dot} | x_i
# lpzi <- function(z,xi,pis){
#   dmultinom(z,size=xi,log=TRUE)
# }
# 
# #joint log probability
# lp <- function(z,x,beta,eta,xform=xform,ti=ti){
#   l <- lpx(x,beta,xform)
#   for(i in 1:length(x)){
#     l <- l + lpzi(z[[i]],x[i],
#   }
# }




#hist(ti,freq=F,breaks=100)
#points((0:50)/4,pid(0:50),type="l",col="red")

# meanEmUpdateOld <- function(y,lambda,gamma=0){
#   T <- length(y)
#   obs <- !is.na(y)
#   lamNew <- rep(NA,length(lambda))
#   for(k in 1:length(lambda)){
#     s <- 0:(T-k)
#     b <- sum(pid(s))
#     no <- !obs[s+k]
#     if(any(no))
#       a <- sum(pid(s[no])) / b
#     else
#       a <- 0
#     c <- 0
#     for(d in s){
#       if(obs[k+d]){
#         c <- c + y[k+d]*pid(d) / sum(lambda[1:(k+d)]*pid(k+d-(1:(k+d))))
#       }
#     }
#     if(gamma > 0){
#       r1 <- - 4 * gamma * ((k!=1) + (k!=T))
#       if(k!=1 && k!=T) { 
#         r2 <- 2 * gamma * (lambda[k] + lambda[k-1] + lambda[k] + lambda[k+1]) - b
#       }else if(k==1){
#         r2 <- 2 * gamma * (lambda[k] + lambda[k+1]) - b
#       }else{
#         r2 <- 2 * gamma * (lambda[k] + lambda[k-1]) - b
#       }
#       r3 <- lambda[k] * (a*b + c)
#       roots <- polyroot(c(r3,r2,r1))
#       if(any(Im(roots)>10^-8)){
#         browser()
#         stop("non-real roots")
#       }
#       #browser()
#       lamNew[k] <- max(Re(roots))
#     }else{
#       lamNew[k] <- lambda[k] * (a + c / b)
#     }
#   }
#   lamNew
# }


