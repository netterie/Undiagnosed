########################
# internal_fxns.R
# 
# Purpose: 
# Store functions used internally for HIV backcalculation, and
# provide a place for developing functions that can eventually
# be moved into the next package release
#
# Dependencies:
# None
#
# History: 
#
# Notes:
# Begun 6/29/15 by JKB
#
# assignInNamespace('runBackCalc', runBackCalc, ns='HIVBackCalc')
########################

# Switch on/off sourcing of functions in progress
source_inProgress_Aug16=TRUE
source_inProgress_May16=FALSE
source_inProgress_Oct15=FALSE

######################################################################
# THE FOLLOWING FUNCTIONS WERE ADDED HERE STARTING 8/9/16
######################################################################
if (source_inProgress_Aug16) {

# ADDED A PANELING OPTION
plot.results <- function(x, panel=NULL) {

    if (is.null(panel)) {
        d <- x$resultsAll 
    } else {
        d <- x
        d$Group <- d[,panel]
    }

    p <- ggplot(d,aes(x=time,y=value, linetype=var))  +   
      geom_line(aes(color=var), size=0.5) +
      geom_point(aes(color=var, shape=var), size=2) + 
      scale_alpha_manual(values=c(.5,1,1),name="") + 
      scale_color_manual(name="", values=c("gray3", "blue", "orange2")) + 
      scale_linetype_manual(name="",values=c(6,3,3)) + 
      scale_shape_manual(name="", values=c(3,16,16)) +
      xlab("Time") + ylab("Counts") + 
      geom_blank(aes(x=2008,y=0)) +
      scale_y_continuous(expand=c(.15,0)) + 
      theme_bw() + 
      theme(text = element_text(size = 10)) +
      theme(legend.position="bottom",axis.text.x=element_text(angle=90))
   
    if (!is.null(panel)) {
        p <- p+facet_grid(group~Group, scales='free_y')
    } else p <- p+ facet_grid(group~.,scales="free_y") 

    return(p)
}
assignInNamespace('plot.results', plot.results, ns='HIVBackCalc')

# FIXED THE PANEL OPTION
plotTestHist <- function(testhist, panel=NULL) {

    if (is.null(panel)) vars <- 'yearDx' else vars <- list(c('yearDx', panel))
    tabTime <- tabTestHist(testhist, vars)

    keep.vars <- c(panel, 'yearDx', grep('Percent ', colnames(tabTime),
                                       value=TRUE))
    these.ids <- c(panel, 'yearDx')
    tabTime <- tabTime[,keep.vars]
    colnames(tabTime) <- gsub('Percent ','',colnames(tabTime))
    tabTime <- melt(tabTime, id.vars=these.ids)
    if (!is.null(panel)) tabTime$Group <- tabTime[,panel]

    p <- ggplot(tabTime,aes(x=yearDx,y=value,group=variable))  +   
    geom_line(aes(color=variable)) +
    geom_point(aes(color=variable)) + 
    theme_bw()+
    theme(legend.position='bottom',axis.text.x=element_text(angle=90)) + 
    scale_color_hue(name="Ever had negative test?") + 
    scale_x_continuous(breaks=seq(min(tabTime$yearDx),max(tabTime$yearDx),by=2))+
    xlab("Time") + ylab("Percent") 

    if (!is.null(panel)) p <- p + facet_grid(.~Group)
    return(p)
}
assignInNamespace('plotTestHist', plotTestHist, ns='HIVBackCalc')

combineResults <- function(x) {

  # Times with observed diagnoses
  allTimes <- as.numeric(names(x[[1]][[1]]$y))
  obsTimes <- !is.na(allTimes)
  times <- allTimes[obsTimes]
  x$times <- times
  
  # Diagnoses
  x$diagnoses <- x[[1]][[1]]$y[obsTimes]
  
  # Incidence and Undiagnosed organized by case
  cases <- names(x)[!names(x)%in%c('times', 'diagnoses')]
  for (c in cases) {
    incidence <- x[[c]][[1]]$lambda[obsTimes]
    undiagnosed <- x[[c]][[2]][obsTimes]
    x[[c]] <- list(incidence=incidence, undiagnosed=undiagnosed)
  }
  
  # Data frame with all results
  x$resultsAll <- data.frame(time=times, 
                             group='Diagnoses and Incidence', 
                             var='# Diagnosed',
                             value=x$diagnoses)

  cat('Is it here???\n')
  for (c in cases) {
    x$resultsAll  <- rbind(x$resultsAll,
                           data.frame(time=times, 
                                      group='Diagnoses and Incidence', 
                                      var=c,
                                      value=x[[c]]$incidence),
                           data.frame(time=times,
                                      group='Undiagnosed Cases',
                                      var=c,
                                      value=x[[c]]$undiagnosed))
  }
  cat('Is this too late???\n')
  
  # Data frame with summarized results
  x$resultsSummary <- ddply(x$resultsAll, .(var, group), function(x) c(summary(x$value)))
  x$resultsSummary <- within(x$resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(x$resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  x$resultsSummaryYear <- ddply(transform(x$resultsAll, Year=floor(time)), 
                            .(var, group, Year), function(x) c(summary(x$value)))
  x$resultsSummaryYear <- within(x$resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(x$resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  class(x) <- append(class(x), 'results')
  return(x)
}
assignInNamespace('combineResults', combineResults, ns='HIVBackCalc')

runBackCalc = function(testhist, intLength, cases=NULL, prev=NULL,...) {
  
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  
  # Estimate TIDs
  TIDs <- estimateTID(testhist$infPeriod, 
                      intLength,
                      cases,
                      infPeriodOrig=testhist$infPeriod,
                      ...)
  
  # Diagnoses
  diagCounts = tabulateDiagnoses(testhist, intLength)
  
  # Initialize incidence and undiagnosed count lists
  incidence <- vector(mode='list', length=length(cases))
  names(incidence) <- cases
  undiagnosed <- incidence
  
  # Estimate incidence and undiagnosed
  for (c in cases) {
    cat('\nEstimating case', c, '...\n')
    incidence[[c]] = estimateIncidence(y=diagCounts,
                                      pid=TIDs[[c]]$pdffxn,
                                      gamma=0.1,
                                      verbose=FALSE)
    undiagnosed[[c]] <- estimateUndiagnosed(incidence[[c]])
  }

  # Compile results - this code allows there to be
  # any number of cases
  toCombine <- vector(mode='list',length=length(cases))
  names(toCombine) <- names(cases)
  for (c in cases) {
      toCombine[[which(cases==c)]] <- list(incidence[[c]], undiagnosed[[c]])
  }
  results <- combineResults(toCombine)

  # True prevalence
  if (!is.null(prev)) trueprev <- calcTruePrev(results, prev) else trueprev <- NULL
  
  return(list(TIDs=TIDs, results=results, trueprev=trueprev, N=nrow(testhist)))
}
assignInNamespace('runBackCalc', runBackCalc, ns='HIVBackCalc')

#' Optional wrapper function to run and compile results for subgroups
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param subvar Name of the variable within the testhist data frame that defines
#'        subgroups within which to run the backcalculation
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional frame with 1st column 'Year' and a 2nd column with PLWHA
#'        for the population represented in the testhist object
#' @param save  Optional file path to save compiled true prevalence results

runSubgroups = function(testhist, subvar, intLength, cases=NULL, 
                        prev=NULL, save=NULL, ...) {
  
  # Subvar
  if (is.numeric(testhist[,subvar])) {
    warning('Subgroup variable will be coerced to character')
    testhist[,subvar] <- as.character(testhist[,subvar])
  }
  subgroups <- unique(testhist[,subvar])
  numsub <- length(subgroups)
  
  # Cases
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  # Prevalence
  if (!is.null(prev)) {
      # Check that if prevalence is given, it is given for all the 
      # subgroups. If subvar = 'stageGroup', just make fake prev data
      # because we don't need the subgroup results, just the total-weighted
      if (sum(subgroups%in%colnames(prev))!=numsub) stop('In runSubGroups, 
                              prevalence data are insufficient')
  }
  
  # Prepare to store results for each subgroup
  subResults <- vector('list', length=(numsub+1))
  names(subResults) <- c(as.character(subgroups), 'Total-stratified')
  
  # Loop through subgroups
  for (s in subgroups) {
    
    cat('\nSUBGROUP: ', s, '\n')

    # FOR CD4 CASE (make more clever): Identify if median windows have been passed
    if ('medWindowsVar'%in%names(list(...))) {
        medWindforCD4 <- testhist[testhist[,subvar]==s,list(...)[['medWindowsVar']]]
    } else medWindforCD4 <- NULL

    # Run the backcalculation for this subgroup, selecting
    # the correct prevalence column if applicable
    if (!is.null(prev)) {
      subPrev <- prev[, c('Year', as.character(s))]
    } else subPrev <- NULL
    subResults[[s]] <- runBackCalc(testhist[testhist[,subvar]==s,], 
                                   intLength,
                                   cases,
                                   prev=subPrev,
                                   medWindows=medWindforCD4,
                                   ...)
    
      # Add a subgroup identifier to the compiled results
      for (r in c('resultsAll', 'resultsSummary', 'resultsSummaryYear')) {
        subResults[[s]]$results[[r]] = data.frame(Subgroup = s,
                                                  subResults[[s]]$results[[r]],
                                                  check.names=FALSE)
      }
  }

  # Extract the results in order to get subgroup-stratified totals
  resultsAllList <- lapply(lapply(subResults, `[[`, 'results'), `[[`, 'resultsAll')

  # Standardize the times common across the groups - some groups may have diagnoses
  # in years or quarters earlier or later than others. Because we're taking averages
  # of incident/undiagnosed cases across time periods rather than sums, let's just
  # remove the extra time periods. Otherwise we would have to impute something reasonable,
  # since we do have to sum across subgroups - can't just put in a zero.

  times <- lapply(subgroups, function(x) subResults[[x]]$results$resultsAll$time)
  mintime <- max(sapply(times, min))
  maxtime <- min(sapply(times, max))
  keeptimes <- seq(mintime, maxtime, by=intLength)

  # Vectors of results for just the keeptimes
  resultsAllList <- lapply(resultsAllList, function(x) x[x$time%in%keeptimes,]$value)

  # Back to extracting results
  resultsAll <- cbind(subResults[[1]]$results$resultsAll[,c('time', 'group', 'var')],
                      do.call(cbind, resultsAllList))
  resultsAll$value <- apply(as.matrix(resultsAll[,as.character(subgroups)]),1,sum)
  
  # Summarize subgroup-stratified totals
  # Data frame with summarized results
  resultsSummary <- ddply(resultsAll, .(var, group), function(x) c(summary(x$value)))
  resultsSummary <- within(resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  resultsSummaryYear <- ddply(transform(resultsAll, Year=floor(time)), 
                                .(var, group, Year), function(x) c(summary(x$value)))
  resultsSummaryYear <- within(resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Create the results list of class 'results'
    resultsL=list(resultsAll=resultsAll,
                      resultsSummary=resultsSummary,
                      resultsSummaryYear=resultsSummaryYear)
    class(resultsL) <- append(class(resultsL), 'results')

  # Save in subResults[['Total-stratified']]$results object
  subResults[['Total-stratified']] <- list(results=resultsL)

  if (!is.null(prev)) {
    
    # Calculate total-stratified true prevalence
    subResults[['Total-stratified']]$trueprev <- 
      calcTruePrev(subResults[['Total-stratified']]$results,
                   prev=data.frame(Year=prev$Year,
                                   Total=apply(as.matrix(prev[,as.character(subgroups)]),
                                               1,sum)))
  }
  
  if (!is.null(save)) {
    trueprev <- do.call(rbind, lapply(names(subResults), 
                               function(x) {
                                 data.frame(Subgroup=x,
                                            subResults[[x]]$trueprev,
                                            check.names=FALSE)
                                }))
    write.csv(trueprev[, !colnames(trueprev) %in% c('1st Qu.', 'Median', '3rd Qu.')],
              file=save,
              row.names=FALSE)
  }
  
  return(subResults)
}
assignInNamespace('runSubgroups', runSubgroups, ns='HIVBackCalc')

TID_cdfContinuous <- function(infPeriod,case,intLength,
                              survivor=FALSE, 
                              maxWindow=17.98, 
                              ...) {

    # Grab optional arguments
    optArgs <- list(...)
    # Remove missing infPeriods, sort, and warn if there are zeroes
    infPeriod <- sort(infPeriod[!is.na(infPeriod)]) 
    if(sum(infPeriod==0)!=0) stop('Please resolve infPeriod=0 cases before proceeding, e.g. 
                                  set them to infPeriod=NA (missing) if you cannot identify
                                  a valid infPeriod.')
    n <- length(infPeriod)

    # Define the function that distributes probability uniformly
    # Function that gets the heights of the probability
    # bars, but not the widths
    pi_uniform <- function(uniqueInfs,infPeriod){
        sapply(uniqueInfs,function(eachInf){
          # Select only larger windows because the 
          # smaller ones have stopped contributing probability 
          # of infection
          largerInfs <- infPeriod[infPeriod>=eachInf]
          # Sum the probabilities (1/x_i) across people and then
          # divide by the number of people. Because of fractional
          # windows, this sum can be greater than 1
          sum(1/largerInfs)/length(infPeriod)
        })
    }

    # Define the CDF of time from infection to diagnosis in 
    # the function 'qi'
    switch(case,

         'upper_bound' = {

             # Simply the empirical CDF of the infPeriods
              qi <- function(u) {
                uind <- sum(infPeriod<=u)/n
                if (is.na(uind)) 
                  return(0)
                uind
              }
         },

         'cd4_case' = {

              if (is.null(optArgs$medWindows)) stop('Problem in CD4 Case')
              keepThese <- !is.na(optArgs$infPeriodOrig) & 
                  optArgs$infPeriodOrig!=0
              infPeriod <- optArgs$infPeriodOrig[keepThese]
              medWindows <- optArgs$medWindows[keepThese]
              cdfCD4 <- function(t, infP=infPeriod, 
                                 m=medWindows) {
                  # Heights of probability bars for each person, pre-median
                  heightPre <- 0.5/(medWindows)
                  # Post-median heights
                  heightPost <- 0.5/(infP-medWindows)
                  # Multiply by width (time of eval) for each person's probability
                  # contribution, but the max for each person is 0.5 for each
                  # of the pre- and post-median periods
                  probPre <- pmin(heightPre*t,0.5)
                  probPost <- pmin(pmax(heightPost*(t-medWindows),0),0.5)
                  # Sum and scale by n
                  cdf <- sum((probPre+probPost)/length(infP))
                  return(cdf)
              }
              # Evaluate at the unique infPeriods
              # medWindows <- infPeriod/2
              # cs_cd4 <- sapply(unique(infPeriod), cdf)

              # CDF function
              qi <- function(u) return(cdfCD4(u))
        },
         'base_case_alt' = {

              # A Base Case replica but using a different computational
              # approach
              cdf <- function(t, infP=infPeriod) {
                  # Heights of probability bars for each person
                  height <- 1/infP
                  # Multiply by width (time of eval) for each person's probability
                  # contribution, but the max for each person is 1.
                  prob <- pmin(height*t,1)
                  # Sum and scale by n
                  cdf <- sum(prob/length(infP))
                  return(cdf)
              }
              # Evaluate at the unique infPeriods
              # cs <- sapply(unique(infPeriod), cdf)

              # CDF function
              qi <- function(u) return(cdf(u))
        },
         'base_case' = {

              # Have to multiply by the diff because 
              # each person's area of 1, (1/x)*x, has to 
              # be apportioned out over the different sections of 
              # their infPeriod (x), and the sections are defined
              # by the unique infs that come before x.
              uniqueInf <- unique(infPeriod)
              p<-pi_uniform(uniqueInf,infPeriod) * diff(c(0,uniqueInf))
              cs <- cumsum(p)

              # CDF of density
              qi <- function(u){
                  # cs is a vector where the index (1:length(infPeriod)) 
                  # indicates the cumulative probability for the 
                  # unique infPeriod corresponding to that index. So if u is
                  # the time of eval, uind is the corresponding index in
                  # uniqueInf
                uind <- rev(which(uniqueInf<=u))[1]
                # If uind is NA, that means we're outside our data min - 
                # if there are no uniqueInf<=u, then the uind phrase
                # above returns NA
                if(is.na(uind))
                  return(0)
                cs[uind]
              }
        }


    ) # end switch
    if (survivor) {
        Sx <- function(u) { 1-qi(u) }
        return(Sx)
    } else return(qi)
} # end TID_cdfContinuous
assignInNamespace('TID_cdfContinuous', TID_cdfContinuous, ns='HIVBackCalc')

#' Evaluates the discrete PDF of time from infection to diagnosis in reporting intervals.
#'  
#' Calls the TID_cdfContinuous function to estimate CDF of the TID, then
#' evaluates the discrete PDF in time intervals that match reporting period
#' units (intLength).
#'  
#' @param infPeriod A vector of continuous times from last HIV test to diagnosis
#'          for a population
#' @param case One of "base_case" or "upper_bound", indicating the 
#'          assumption to apply for when infection occurred within infPeriod
#' @param intLength A single number indicating the length in years of discrete 
#'          time intervals by which HIV diagnoses are recorded. The default of 
#'          0.25 represents a quarter-year.
#' @param ... Additional parameters to pass to TID_cdfContinuous
#'  
#' @return A function that takes integer arguments 0 and higher. Will return 0 
#'          for integers greater than floor(max(infPeriod/intLength)) + 1
TID_pdf <- function(infPeriod,case,intLength,...) {

    # Remove missing infPeriods 
    infPeriod <- sort(infPeriod[!is.na(infPeriod)]) 
    n <- length(infPeriod)

    # Define the CDF of time from infection to diagnosis in 
    # the function 'qi'
    qi <- TID_cdfContinuous(infPeriod,case,intLength,...)

    # Use the CDF to define the discrete probability of infection 
    # during interval i to i+1
    pidCalc <- function(i){
    sapply(i,function(ii){
      qi((ii+1)*intLength) - qi(ii*intLength)
    })
    }

    # Calculate the discrete PDF over the m observed intervals, and set
    # probability to zero for longer intervals
    m <- max(infPeriod/intLength) + 1
    pidProbs <- pidCalc(0:m)
    pid <- function(i){
    ifelse(i>m,0,pidProbs[i+1])
    }

    # Return the PDF function
    return(pid)
}

#' Creates a TID (time from infection to diagnosis) object
#'  
#' Calls the TID_pdf() function to estimate the TID from testing history
#' data for each assumption/case. Evaluates the function for the observed
#' time span in the data and returns both the probability and cumulative
#' density distributions.
#'  
#' @param infPeriod The vector containing infection periods, or times 
#'          (in years) between last negative test and diagnosis
#' @param intLength The interval length by which diagnoses are reported
#'          (also in years, 1=1 year)
#' @param cases Cases to estimate; default is c('base_case', 'upper_bound')
#' @param ... Additional parameters to pass to TID_cdfContinuous
#'  
#' @return A nested list. The first tier indicates the assumption used
#'          to estimate the TID. The second tier contains 3 elements:
#'          "pdffxn", the PDF function, and "pdf" and "cdf" which 
#'          indicate the respective distributions
estimateTID <- function(infPeriod, intLength, cases=NULL,...) {

    # Default cases
    if (is.null(cases)) cases <- c('base_case', 'upper_bound')

    # TID object
    TIDobject <- vector(length=length(cases), mode='list')
    names(TIDobject) <- cases

    # Intervals observed in infPeriod
    maxInt <- max(infPeriod/intLength, na.rm=TRUE)+1

    # Populate with TID functions and actual distributions
    for (c in cases) {
        TIDobject[[c]] <- vector(mode='list', length=3)
        names(TIDobject[[c]]) <- c('pdffxn', 'pdf', 'cdf')
        # PDF function
        TIDobject[[c]]$pdffxn <- TID_pdf(infPeriod, c, intLength,...)
        # PDF
        TIDobject[[c]]$pdf <- sapply(0:maxInt, TIDobject[[c]]$pdffxn)
        # CDF
        TIDobject[[c]]$cdf <- cumsum(TIDobject[[c]]$pdf)
    }

    class(TIDobject) <- append(class(TIDobject), 'TID')
    return(TIDobject)
}

} # end source_inProgress_Aug16

######################################################################
# THE FOLLOWING FUNCTIONS WERE ADDED HERE STARTING 5/12/16 
######################################################################
if (source_inProgress_May16) {
######################################################################
# runSubgroups
######################################################################

#' Optional wrapper function to run and compile results for subgroups
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param subvar Name of the variable within the testhist data frame that defines
#'        subgroups within which to run the backcalculation
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional frame with 1st column 'Year' and a 2nd column with PLWHA
#'        for the population represented in the testhist object
#' @param save  Optional file path to save compiled true prevalence results

runSubgroups = function(testhist, subvar, intLength, cases=NULL, 
                        prev=NULL, save=NULL) {
  
  # Subvar
  if (is.numeric(testhist[,subvar])) {
    warning('Subgroup variable will be coerced to character')
    testhist[,subvar] <- as.character(testhist[,subvar])
  }
  subgroups <- unique(testhist[,subvar])
  numsub <- length(subgroups)
  
  cat('Here1')
  # Cases
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  cat('Here2')
  # Prevalence
  if (!is.null(prev)) {
      # Check that if prevalence is given, it is given for all the 
      # subgroups. If subvar = 'stageGroup', just make fake prev data
      # because we don't need the subgroup results, just the total-weighted
      if (sum(subgroups%in%colnames(prev))!=numsub) stop('In runSubGroups, 
                              prevalence data are insufficient')
  }
  
  # Prepare to store results for each subgroup
  subResults <- vector('list', length=(numsub+1))
  names(subResults) <- c(as.character(subgroups), 'Total-stratified')
  
  cat('Here3')
  # Loop through subgroups
  for (s in subgroups) {
    
    cat('\nSUBGROUP: ', s, '\n')
    # Run the backcalculation for this subgroup, selecting
    # the correct prevalence column if applicable
    if (!is.null(prev)) {
      subPrev <- prev[, c('Year', as.character(s))]
    } else subPrev <- NULL
    subResults[[s]] <- runBackCalc(testhist[testhist[,subvar]==s,], 
                                   intLength,
                                   cases,
                                   prev=subPrev)
    
      # Add a subgroup identifier to the compiled results
      for (r in c('resultsAll', 'resultsSummary', 'resultsSummaryYear')) {
        subResults[[s]]$results[[r]] = data.frame(Subgroup = s,
                                                  subResults[[s]]$results[[r]],
                                                  check.names=FALSE)
      }
  }

  cat('Here4')
  # Extract the results in order to get subgroup-stratified totals
  resultsAllList <- lapply(lapply(subResults, `[[`, 'results'), `[[`, 'resultsAll')

  # Standardize the times common across the groups - some groups may have diagnoses
  # in years or quarters earlier or later than others. Because we're taking averages
  # of incident/undiagnosed cases across time periods rather than sums, let's just
  # remove the extra time periods. Otherwise we would have to impute something reasonable,
  # since we do have to sum across subgroups - can't just put in a zero.

  times <- lapply(subgroups, function(x) subResults[[x]]$results$resultsAll$time)
  mintime <- max(sapply(times, min))
  maxtime <- min(sapply(times, max))
  keeptimes <- seq(mintime, maxtime, by=intLength)

  cat('Here5')
  # Vectors of results for just the keeptimes
  resultsAllList <- lapply(resultsAllList, function(x) x[x$time%in%keeptimes,]$value)

  # Back to extracting results
  resultsAll <- cbind(subResults[[1]]$results$resultsAll[,c('time', 'group', 'var')],
                      do.call(cbind, resultsAllList))
  resultsAll$value <- apply(as.matrix(resultsAll[,as.character(subgroups)]),1,sum)
  
  # Summarize subgroup-stratified totals
  # Data frame with summarized results
  resultsSummary <- ddply(resultsAll, .(var, group), function(x) c(summary(x$value)))
  resultsSummary <- within(resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  resultsSummaryYear <- ddply(transform(resultsAll, Year=floor(time)), 
                                .(var, group, Year), function(x) c(summary(x$value)))
  resultsSummaryYear <- within(resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Create the results list of class 'results'
    resultsL=list(resultsAll=resultsAll,
                      resultsSummary=resultsSummary,
                      resultsSummaryYear=resultsSummaryYear)
    class(resultsL) <- append(class(resultsL), 'results')

  # Save in subResults[['Total-stratified']]$results object
  subResults[['Total-stratified']] <- list(results=resultsL)

  if (!is.null(prev)) {
    
    # Calculate total-stratified true prevalence
    subResults[['Total-stratified']]$trueprev <- 
      calcTruePrev(subResults[['Total-stratified']]$results,
                   prev=data.frame(Year=prev$Year,
                                   Total=apply(as.matrix(prev[,as.character(subgroups)]),
                                               1,sum)))
  }
  
  if (!is.null(save)) {
    trueprev <- do.call(rbind, lapply(names(subResults), 
                               function(x) {
                                 data.frame(Subgroup=x,
                                            subResults[[x]]$trueprev,
                                            check.names=FALSE)
                                }))
    write.csv(trueprev[, !colnames(trueprev) %in% c('1st Qu.', 'Median', '3rd Qu.')],
              file=save,
              row.names=FALSE)
  }
  
  return(subResults)
}
assignInNamespace('runSubgroups', runSubgroups, ns='HIVBackCalc')
######################################################################
# runBackCalc
######################################################################

#' Optional wrapper function to run all the backcalculation steps
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional data frame with 1st column 'Year' and a 2nd column with
#'        PLWHA for the population represented in the testhist object
runBackCalc = function(testhist, intLength, cases=NULL, prev=NULL) {
  
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  
  # Estimate TIDs
  TIDs <- estimateTID(testhist$infPeriod, 
                      intLength,
                      cases)
  
  # Diagnoses
  diagCounts = tabulateDiagnoses(testhist, intLength)
  
  # Initialize incidence and undiagnosed count lists
  incidence <- vector(mode='list', length=length(cases))
  names(incidence) <- cases
  undiagnosed <- incidence
  
  # Estimate incidence and undiagnosed
  for (c in cases) {
    cat('\nEstimating case', c, '...\n')
    incidence[[c]] = estimateIncidence(y=diagCounts,
                                      pid=TIDs[[c]]$pdffxn,
                                      gamma=0.1,
                                      verbose=FALSE)
    undiagnosed[[c]] <- estimateUndiagnosed(incidence[[c]])
  }

  # Compile results - this code allows there to be
  # any number of cases
  toCombine <- vector(mode='list',length=length(cases))
  names(toCombine) <- names(cases)
  for (c in cases) {
      toCombine[[which(cases==c)]] <- list(incidence[[c]], undiagnosed[[c]])
  }
  results <- combineResults(toCombine)

  # True prevalence
  if (!is.null(prev)) trueprev <- calcTruePrev(results, prev) else trueprev <- NULL
  
  return(list(TIDs=TIDs, results=results, trueprev=trueprev, N=nrow(testhist)))
}
assignInNamespace('runBackCalc', runBackCalc, ns='HIVBackCalc')

}


######################################################################
# setup_hivbackcalc
######################################################################

#' Setup function
#' @param workd Working directory full file path
#' @param datafile Path to data, either full or from the workd
#' @param source_these Vector of paths to other files to source, either full or from the workd
#' @param loadlib Logical indicating whether to load libries
#' @param msm Logical indicating whether this is the KC MSM analysis where the data frame was called msm
#' @param na Character indicating whether NAs are "NA" or spaces ("")
#' @param package_updated DEPRECATED
#' @param load_package Logical indicating whether the HIVBackCalc package 
#'        should be loaded (TRUE). If you want to have some package functions but overwrite others,
#'        specify TRUE here and a file in the 'packagefile' param
#' @param packagefile Within the workd, path to file with newer package functions
setup_hivbackcalc = function(workd, datafile, source_these, loadlib=TRUE,
                             msm=FALSE, na="", load_package=TRUE,
                             package_updated=NULL, packagefile=NULL) {
  
  if (loadlib) {
    cat('Loading libraries...\n')
    # Load libraries
    # Eventually, make these dependencies of the HIBBackCalc package?
    library(reshape2)
    library(plyr)
    library(ggplot2)
    library(scales)
    library(Hmisc)
  }

  if (!is.null(package_updated)) stop('package_updated is deprecated; use load_package instead')
  if (load_package) library(HIVBackCalc)
  if (!is.null(packagefile)) source_these = c(source_these, packagefile)
  
  # Working directory
  workd <<- workd
  
  # Load data
  cat('Loading data and storing it in object msm or dataf...\n')
  if (msm) {
    msm <<- read.csv(file.path(workd,datafile),na.string=na,stringsAsFactor=FALSE) 
  } else {
    dataf <<- read.csv(file.path(workd,datafile),na.string=na,stringsAsFactor=FALSE) 
  }
  
  # Source files
  for (f in source_these) {
    cat('Sourcing', f, '...\n')
    source(file.path(workd,f))
  }
  
}

######################################################################
# THE FOLLOWING FUNCTIONS WERE ADDED HERE STARTING 12/29/15 TO 
# Oct15 BRANCH OF THE PACKAGE
######################################################################

if (source_inProgress_Oct15) {

TID_cdfContinuous <- function(infPeriod,case,intLength,
                              survivor=FALSE,
                              stageGroup=NULL) {

    # Address stageGroup. This is for the extended method that 
    # uses BED/dual diagnosis information. Groups 
    # that don't have modified code will be re-routed to
    # the original base case
    if (case=='base_case_withStage') {
        if (stageGroup%in%c('BEDmDD+', 'BEDmDD-', 
                            'BED+DD+', 'BEDm+DD+')) case <- 'base_case'
    }

    # Remove missing infPeriods, sort, and warn if there are zeroes
    infPeriod <- sort(infPeriod[!is.na(infPeriod)]) 
    if(sum(infPeriod==0)!=0) stop('Please resolve infPeriod=0 cases before proceeding, e.g. 
                                  set them to infPeriod=NA (missing) if you cannot identify
                                  a valid infPeriod.')
    n <- length(infPeriod)

    # Define the function that distributes probability uniformly

    # Function that gets the heights of the probability
    # bars, but not the widths
    pi_uniform <- function(uniqueInfs,infPeriod){
        sapply(uniqueInfs,function(eachInf){
          # Select only larger windows because the 
          # smaller ones have stopped contributing probability 
          # of infection
          largerInfs <- infPeriod[infPeriod>=eachInf]
          # Sum the probabilities (1/x_i) across people and then
          # divide by the number of people. Because of fractional
          # windows, this sum can be greater than 1
          sum(1/largerInfs)/length(infPeriod)
        })
    }
    # Modified function that gets the heights of the probability
    # bars as well as the widths
    pi_BEDneg <- function(uniqueInfs,infPeriod,BEDw){
        sapply(uniqueInfs,function(eachInf){
                # The windows that are still going at eachInf
                largerInfs <- infPeriod[infPeriod>=eachInf]
                # Compare infs to BEDw
                largerThanBED <- largerInfs>BEDw
                # Width of probability bar
                eachInfIndex <- which(uniqueInfs==eachInf)
                width <- ifelse(eachInfIndex>1,
                                uniqueInfs[eachInfIndex]-
                                uniqueInfs[eachInfIndex-1],
                                eachInf)
                if (eachInf<=BEDw) {
                    # Height of probability bar
                    toSum <- 1/largerInfs[!largerThanBED]
                } else {
                    # Height of probability bar
                    toSum <- 1/(largerInfs[largerThanBED]-BEDw) 
                    # Which uniqueInf is closest and larger than the BEDw
                    diffs <- uniqueInfs-BEDw
                    afterBEDw <- which(diffs==min(diffs[diffs>0]))
                    # Modify the width for this uniqueInf
                    if (eachInfIndex==afterBEDw) width <- diffs[eachInfIndex]
                }
                
                # Get height*width and scale by # ppl
                return(sum(toSum)*width/length(infPeriod))
            })
    }

    # Define the CDF of time from infection to diagnosis in 
    # the function 'qi'
    switch(case,

         'upper_bound' = {

             # Simply the empirical CDF of the infPeriods
              qi <- function(u) {
                uind <- sum(infPeriod<=u)/n
                if (is.na(uind)) 
                  return(0)
                uind
              }
         },

         'base_case' = {

              # Have to multiply by the diff because 
              # each person's area of 1, (1/x)*x, has to 
              # be apportioned out over the different sections of 
              # their infPeriod (x), and the sections are defined
              # by the unique infs that come before x.
              uniqueInf <- unique(infPeriod)
              p<-pi_uniform(uniqueInf,infPeriod) * diff(c(0,uniqueInf))
              cs <- cumsum(p)

              # CDF of density
              qi <- function(u){
                  # cs is a vector where the index (1:length(infPeriod)) 
                  # indicates the cumulative probability for the 
                  # unique infPeriod corresponding to that index. So if u is
                  # the infPeriod, uind is the corresponding index in
                  # uniqueInf
                uind <- rev(which(uniqueInf<=u))[1]
                # If uind is NA, that means we're outside our data min - 
                # if there are no uniqueInf<=u, then the uind phrase
                # above returns NA
                if(is.na(uind))
                  return(0)
                cs[uind]
              }
        },

         'base_case_withStage' = {

              if (is.null(stageGroup)) stop('In TID_cdfContinuous, specify stageGroup')

              # Define BED window
              BEDw <- 162/365

              if (stageGroup=='BED+DD-' | stageGroup=='BED+') {

                  # Create a modified infPeriod
                  infPeriodMod <- sort(ifelse(infPeriod>=BEDw, BEDw, infPeriod))

                  uniqueInf <- unique(infPeriodMod)
                  p<-pi_uniform(uniqueInf,infPeriodMod) * diff(c(0,uniqueInf))
              } else if (stageGroup=='BED-DD+') {

                  # For now, we are not using the AIDS incubation distribution - 
                  # ultimately this section should be different from BED-DD-
                  uniqueInf <- unique(infPeriod)
                  p<-pi_BEDneg(uniqueInf,infPeriod,BEDw) 
              } else if (stageGroup=='BED-DD-' | stageGroup=='BED-') {

                  uniqueInf <- unique(infPeriod)
                  p<-pi_BEDneg(uniqueInf,infPeriod,BEDw) 
              } else {
                  stop('In TID_cdfContinuous, not coded yet')
              }

              cs <- cumsum(p)
              # Warn if probability error
              if (abs(sum(p)-1)>.000000001) stop(paste('In TID_cdfContinuous, sum of p!=1 for group', stageGroup))

              # CDF of density
              qi <- function(u){
                uind <- rev(which(uniqueInf<=u))[1]
                if(is.na(uind))
                  return(0)
                cs[uind]
              }
        },

         'upper_bound_withStage' = {

              if (is.null(stageGroup)) stop('In TID_cdfContinuous, specify stageGroup')
              BEDw <- 162/365

              if (stageGroup=='BED+DD-' | stageGroup=='BED+') {

                  # Create a modified infPeriod
                  infPeriodMod <- sort(ifelse(infPeriod>=BEDw, BEDw, infPeriod))

              } else infPeriodMod <- infPeriod

             # Simply the empirical CDF of the infPeriods
              qi <- function(u) {
                uind <- sum(infPeriodMod<=u)/n
                if (is.na(uind)) 
                  return(0)
                uind
              }
        }

    ) # end switch
    if (survivor) {
        Sx <- function(u) { 1-qi(u) }
        return(Sx)
    } else return(qi)
} # end TID_cdfContinuous
assignInNamespace('TID_cdfContinuous', TID_cdfContinuous, ns='HIVBackCalc')

#' Evaluates the discrete PDF of time from infection to diagnosis in reporting intervals.
#'  
#' Calls the TID_cdfContinuous function to estimate CDF of the TID, then
#' evaluates the discrete PDF in time intervals that match reporting period
#' units (intLength).
#'  
#' @param infPeriod A vector of continuous times from last HIV test to diagnosis
#'          for a population
#' @param case One of "base_case" or "upper_bound", indicating the 
#'          assumption to apply for when infection occurred within infPeriod
#' @param intLength A single number indicating the length in years of discrete 
#'          time intervals by which HIV diagnoses are recorded. The default of 
#'          0.25 represents a quarter-year.
#'  
#' @return A function that takes integer arguments 0 and higher. Will return 0 
#'          for integers greater than floor(max(infPeriod/intLength)) + 1
TID_pdf <- function(infPeriod,case,intLength,
                    stageGroup=NULL) {

    if (!is.null(stageGroup)) {
        # Remove missing infPeriods or, for BED+, impute BEDw
        if (grepl('withStage', case) & (stageGroup=='BED+DD-' | stageGroup=='BED+')) {
            BEDw <- 162/365
            infPeriod <- ifelse(is.na(infPeriod), BEDw, infPeriod)
        }
    }
    infPeriod <- sort(infPeriod[!is.na(infPeriod)]) 
    n <- length(infPeriod)

    # Define the CDF of time from infection to diagnosis in 
    # the function 'qi'
    qi <- TID_cdfContinuous(infPeriod,case,intLength,
                            survivor=FALSE,
                            stageGroup)

    # Use the CDF to define the discrete probability of infection 
    # during interval i to i+1
    pidCalc <- function(i){
    sapply(i,function(ii){
      qi((ii+1)*intLength) - qi(ii*intLength)
    })
    }

    # Calculate the discrete PDF over the m observed intervals, and set
    # probability to zero for longer intervals
    m <- max(infPeriod/intLength) + 1
    pidProbs <- pidCalc(0:m)
    pid <- function(i){
    ifelse(i>m,0,pidProbs[i+1])
    }

    # Return the PDF function
    return(pid)
}
assignInNamespace('TID_pdf', TID_pdf, ns='HIVBackCalc')

#' Creates a TID (time from infection to diagnosis) object
#'  
#' Calls the TID_pdf() function to estimate the TID from testing history
#' data for each assumption/case. Evaluates the function for the observed
#' time span in the data and returns both the probability and cumulative
#' density distributions.
#'  
#' @param infPeriod The vector containing infection periods, or times 
#'          (in years) between last negative test and diagnosis
#' @param intLength The interval length by which diagnoses are reported
#'          (also in years, 1=1 year)
#'  
#' @return A nested list. The first tier indicates the assumption used
#'          to estimate the TID. The second tier contains 3 elements:
#'          "pdffxn", the PDF function, and "pdf" and "cdf" which 
#'          indicate the respective distributions
estimateTID <- function(infPeriod, intLength, stageGroup=NULL, cases=NULL) {

    # Default cases
    if (is.null(cases)) cases <- c('base_case', 'upper_bound')

    # TID object
    TIDobject <- vector(length=length(cases), mode='list')
    names(TIDobject) <- cases

    # Intervals observed in infPeriod
    maxInt <- max(infPeriod/intLength, na.rm=TRUE)+1

    # Populate with TID functions and actual distributions
    for (c in cases) {
        TIDobject[[c]] <- vector(mode='list', length=3)
        names(TIDobject[[c]]) <- c('pdffxn', 'pdf', 'cdf')
        # PDF function
        TIDobject[[c]]$pdffxn <- TID_pdf(infPeriod, c, intLength,
                                         stageGroup)
        # PDF
        TIDobject[[c]]$pdf <- sapply(0:maxInt, TIDobject[[c]]$pdffxn)
        # CDF
        TIDobject[[c]]$cdf <- cumsum(TIDobject[[c]]$pdf)
    }

    class(TIDobject) <- append(class(TIDobject), 'TID')
    return(TIDobject)
}
assignInNamespace('estimateTID', estimateTID, ns='HIVBackCalc')

######################################################################
# runSubgroups
######################################################################

#' Optional wrapper function to run and compile results for subgroups
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param subvar Name of the variable within the testhist data frame that defines
#'        subgroups within which to run the backcalculation
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional frame with 1st column 'Year' and a 2nd column with PLWHA
#'        for the population represented in the testhist object
#' @param save  Optional file path to save compiled true prevalence results
runSubgroups = function(testhist, subvar, intLength, cases=NULL, 
                        prev=NULL, save=NULL) {
  
  # Subvar
  if (is.numeric(testhist[,subvar])) {
    warning('Subgroup variable will be coerced to character')
    testhist[,subvar] <- as.character(testhist[,subvar])
  }
  subgroups <- unique(testhist[,subvar])
  numsub <- length(subgroups)
  
  # Cases
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  # Prevalence
  if (!is.null(prev)) {
      # Check that if prevalence is given, it is given for all the 
      # subgroups. If subvar = 'stageGroup', just make fake prev data
      # because we don't need the subgroup results, just the total-weighted
      if (subvar!='stageGroup') {
          if (sum(subgroups%in%colnames(prev))!=numsub) stop('In runSubGroups, 
                                  prevalence data are insufficient')
      } else {
          prev <- cbind(prev,replicate(numsub,prev$Total/numsub))
          colnames(prev)[(ncol(prev)-numsub+1):ncol(prev)] <- 
              as.character(subgroups)
      }
  }
  
  # Prepare to store results for each subgroup
  subResults <- vector('list', length=(numsub+1))
  names(subResults) <- c(as.character(subgroups), 'Total-stratified')
  
  # Loop through subgroups
  for (s in subgroups) {
    
    cat('\nSUBGROUP: ', s, '\n')
    # Run the backcalculation for this subgroup, selecting
    # the correct prevalence column if applicable
    if (!is.null(prev)) {
      subPrev <- prev[, c('Year', as.character(s))]
    } else subPrev <- NULL
    subResults[[s]] <- runBackCalc(testhist[testhist[,subvar]==s,], 
                                   intLength,
                                   cases,
                                   prev=subPrev)
    
      # Add a subgroup identifier to the compiled results
      for (r in c('resultsAll', 'resultsSummary', 'resultsSummaryYear')) {
        #cat('\n    Saving results for',s,'\n')
        subResults[[s]]$results[[r]] = data.frame(Subgroup = s,
                                                  subResults[[s]]$results[[r]],
                                                  check.names=FALSE)
      }
  }

  # Extract the results in order to get subgroup-stratified totals
  resultsAllList <- lapply(lapply(subResults, `[[`, 'results'), `[[`, 'resultsAll')

  # Standardize the times common across the groups - some groups may have diagnoses
  # in years or quarters earlier or later than others. Because we're taking averages
  # of incident/undiagnosed cases across time periods rather than sums, let's just
  # remove the extra time periods. Otherwise we would have to impute something reasonable,
  # since we do have to sum across subgroups - can't just put in a zero.

  times <- lapply(subgroups, function(x) subResults[[x]]$results$resultsAll$time)
  mintime <- max(sapply(times, min))
  maxtime <- min(sapply(times, max))
  keeptimes <- seq(mintime, maxtime, by=intLength)

  # Vectors of results for just the keeptimes
  resultsAllList <- lapply(resultsAllList, function(x) x[x$time%in%keeptimes,]$value)

  # Back to extracting results
  resultsAll <- cbind(subResults[[1]]$results$resultsAll[,c('time', 'group', 'var')],
                      do.call(cbind, resultsAllList))
  resultsAll$value <- apply(as.matrix(resultsAll[,as.character(subgroups)]),1,sum)
  
  # Summarize subgroup-stratified totals
  # Data frame with summarized results
  resultsSummary <- ddply(resultsAll, .(var, group), function(x) c(summary(x$value)))
  resultsSummary <- within(resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  resultsSummaryYear <- ddply(transform(resultsAll, Year=floor(time)), 
                                .(var, group, Year), function(x) c(summary(x$value)))
  resultsSummaryYear <- within(resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Save in subResults[['Total-stratified']]$results object
#  subResults[['Total-stratified']] <- 
#    list(results=list(resultsAll=resultsAll,
#                      resultsSummary=resultsSummary,
#                      resultsSummaryYear=resultsSummaryYear))
  
  # Create the results list of class 'results'
    resultsL=list(resultsAll=resultsAll,
                      resultsSummary=resultsSummary,
                      resultsSummaryYear=resultsSummaryYear)
    class(resultsL) <- append(class(resultsL), 'results')

  # Save in subResults[['Total-stratified']]$results object
  subResults[['Total-stratified']] <- list(results=resultsL)

  if (!is.null(prev)) {
    
    # Calculate total-stratified true prevalence
    subResults[['Total-stratified']]$trueprev <- 
      calcTruePrev(subResults[['Total-stratified']]$results,
                   prev=data.frame(Year=prev$Year,
                                   Total=apply(as.matrix(prev[,as.character(subgroups)]),
                                               1,sum)))
  }
  
  if (!is.null(save)) {
    trueprev <- do.call(rbind, lapply(names(subResults), 
                               function(x) {
                                 data.frame(Subgroup=x,
                                            subResults[[x]]$trueprev,
                                            check.names=FALSE)
                                }))
    write.csv(trueprev[, !colnames(trueprev) %in% c('1st Qu.', 'Median', '3rd Qu.')],
              file=save,
              row.names=FALSE)
  }
  
  return(subResults)
}
assignInNamespace('runSubgroups', runSubgroups, ns='HIVBackCalc')

######################################################################
# runBackCalc
######################################################################

#' Optional wrapper function to run all the backcalculation steps
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional data frame with 1st column 'Year' and a 2nd column with
#'        PLWHA for the population represented in the testhist object
runBackCalc = function(testhist, intLength, cases=NULL, prev=NULL) {
  
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  # Check if this is stage of infection 
  if ('base_case_withStage'%in%cases | 'upper_bound_withStage'%in%cases) {
      stageGroup <- unique(testhist$stageGroup)
      if (length(stageGroup)!=1) stop('In runBackCalc, error with stage groups')
  } else stageGroup <- NULL
  
  # Estimate TIDs
  TIDs <- estimateTID(testhist$infPeriod, 
                      intLength,
                      stageGroup,
                      cases)
  
  # Diagnoses
  diagCounts = tabulateDiagnoses(testhist, intLength)
  
  # Initialize incidence and undiagnosed count lists
  incidence <- vector(mode='list', length=length(cases))
  names(incidence) <- cases
  undiagnosed <- incidence
  
  # Estimate incidence and undiagnosed
  for (c in cases) {
    cat('\nEstimating case', c, '...\n')
    incidence[[c]] = estimateIncidence(y=diagCounts,
                                      pid=TIDs[[c]]$pdffxn,
                                      gamma=0.1,
                                      verbose=FALSE)
    undiagnosed[[c]] <- estimateUndiagnosed(incidence[[c]])
  }

  # Compile results
  toCombine <- vector(mode='list',length=length(cases))
  names(toCombine) <- names(cases)
  for (c in cases) {
      toCombine[[which(cases==c)]] <- list(incidence[[c]], undiagnosed[[c]])
  }
  results <- combineResults(toCombine)


#  results <- combineResults(list(`Base Case`=list(incidence[['base_case']],
#                                                  undiagnosed[['base_case']]),
#                                 `Upper Bound`=list(incidence[['upper_bound']],
#                                                    undiagnosed[['upper_bound']])))

  # True prevalence
  if (!is.null(prev)) trueprev <- calcTruePrev(results, prev) else trueprev <- NULL

#  warning('TO REMOVE: hardcoded writing to disk in runBackCalc')
#  write.csv(trueprev,
#            paste0('~/Dropbox/School/PhD/HIV_WA/analysis_WA/2016_StageOfInfection_temp_',
#                   stageGroup, '.csv')
#  )

  
  return(list(TIDs=TIDs, results=results, trueprev=trueprev, N=nrow(testhist)))
}
assignInNamespace('runBackCalc', runBackCalc, ns='HIVBackCalc')

######################################################################
# FIXED PLOT.RESULTS
######################################################################

#' Plot estimates incidence and undiagnosed cases
#'  
#' Overlays the backcalculated incidence on diagnoses in one panel and
#' undiagnosed counts in another. Cases are indicated by colors.
#'  
#' @param x List of class "results", the output of combineResults()
#'  
#' @return Paneled ggplot2 object 
plot.results <- function(x) {

    d <- x$resultsAll
    xmintime <- min(d$time)

    ncases <- length(unique(d$var))
    alpha.linetype.shape.rep <- ncases-1
    linecolors <- c('gray3', 'blue', 'orange2')
   if (ncases>3) linecolors <- c(linecolors, 'green4', 'tomato1', 
                                  'green4', 'chocolate4', 
                                  'firebrick4', 'steelblue')

    alphaval <- c(0.5, rep(1, alpha.linetype.shape.rep))
    ltval <- c(6, rep(3, alpha.linetype.shape.rep))
    shapeval <- c(3, rep(16, alpha.linetype.shape.rep))

    p <- ggplot(d,aes(x=time,y=value, linetype=var))  +   
      geom_line(aes(color=var), size=0.5) +
      geom_point(aes(color=var, shape=var), size=2) + 
      scale_alpha_manual(values=alphaval,name="") + 
      scale_color_manual(name="", values=linecolors) + 
      scale_linetype_manual(name="",values=ltval) + 
      scale_shape_manual(name="", values=shapeval) +
      facet_grid(group~.,scales="free_y") +
      xlab("Time") + ylab("Counts") + 
      geom_blank(aes(x=min(time),y=0)) +
      scale_y_continuous(expand=c(.15,0)) + 
      theme_bw() + 
      theme(text = element_text(size = 10)) +
      theme(legend.position="bottom",axis.text.x=element_text(angle=90))

    return(p)
}
assignInNamespace('plot.results', plot.results, ns='HIVBackCalc')

######################################################################
# TABULATE everHadNegTest in percents
######################################################################
#' Tabulate responses to 'Have you ever had a negative test?'
#'
#' Tabulates the everHadNegTest variable by any number of stratification/
#' subgroup variables
#'
#' @param testhist Data frame of class 'testinghistories' with variable
#'        everHadNegTest 
#' @param variables Character vector of stratification variable names. 
#'        Variables must exist in the testhist data frame. Use a list
#'        of a character vector to get cross-tabs.
#' @param supercolumn Set to TRUE to include a pretty column indicating
#'        stratification variables
#' @param fullsample_row Set to TRUE to have the 1st row be the results
#'        for the entire sample
tabTestHist <- function(testhist, variables, supercolumn=FALSE,
                                    fullsample_row=FALSE) {

  vars <- list(NULL)
  for (v in 1:length(variables)) {
    
    # Tabulate everHadNegTest for this subgroup
    tab <- ddply(testhist, variables[[v]], function(x, TN=nrow(testhist)) {
      n <- nrow(x)
      c(N=n,
        `Column Percent`=round(100*n/TN,0),
        `Percent Yes`=round(100*sum(x$everHadNegTest, na.rm=TRUE)/n,0),
        `Percent No`=round(100*sum(!x$everHadNegTest, na.rm=TRUE)/n,0),
        `Percent Missing`=round(100*sum(is.na(x$everHadNegTest))/n,0))
    })
    # Add pretty column with subgroup name or remove repeating labels
    if (supercolumn) {
        if (!is.list(variables)) {
            colnames(tab)[1] <- 'Subgroup'
            tab <- data.frame(Characteristic=rep('',nrow(tab)),
                            tab,
                            check.names=FALSE,
                            stringsAsFactors=FALSE)
            tab$Characteristic[1] <- names(variables)[v]
        } else {
            nVarsToFix <- which(colnames(tab)=='N')-2
            for (i in 1:nVarsToFix) {
                if (is.factor(tab[,i])) tab[,i] <- as.character(tab[,i])
                ulab <- unique(tab[,i])
                # Changepoint
                w <- min(which(tab[,i]==ulab[2]))
                numDel <- w-2
                reps <- nrow(tab)/(w-1)
                vecDel <- c(rep(c(FALSE, rep(TRUE,numDel)),reps))
                tab[,i][vecDel] <- ''
            }
        }
    }
    
    vars[[v]] <- tab
  }
  
  # Compile results into one table
  vars <- do.call(rbind, vars)
  
  # Now add a row for the full sample
  if (fullsample_row) {
    fulleverHadNegTest <- round(100*table(testhist$everHadNegTest, 
                                          useNA='ifany')/nrow(testhist),0)
    fullrow <- data.frame(vars[1,])
    thiscol <- which(colnames(fullrow)=='N')
    fullrow[,c((ncol(fullrow)-4):ncol(fullrow))] <-
      c(nrow(testhist), 
        100, 
        fulleverHadNegTest['TRUE'], 
        fulleverHadNegTest['FALSE'],
        100-fulleverHadNegTest['TRUE']-fulleverHadNegTest['FALSE'])
    fullrow[,1:(thiscol-1)] <- c('Total', rep('', thiscol-2))
    colnames(fullrow) <- colnames(vars)
    vars <- rbind(vars, fullrow)
  }
  
  return(vars)
}
assignInNamespace('tabTestHist', tabTestHist, ns='HIVBackCalc')

} # end source_inProgress

######################################################################
# THE FOLLOWING CODE IS BRACKETED OUT BECAUSE IT IS INCLUDED IN THE
# Oct15 BRANCH OF THE PACKAGE - except runSubgroups?
######################################################################



if (1==0) {

#############################################################
# PLOT DIAGNOSIS COUNTS PER QUARTER
#############################################################
#' Plots diagnosis counts over time
#'
#' Plot of diagnoses vs time, with option to panel by subgroups
#'
#' @param testhist Data frame of class 'testinghistories' with variable
#'        timeDx 
#' @param panel Name of variable in testhist by which to panel the plot
plotDiagnoses  <- function(testhist, panel=NULL) {
  
  if (is.null(names(panel))) names(panel) <- panel
  variables <- c('timeDx', panel[1])
  
  counts <- ddply(testhist, variables, function(x) nrow(x))
  
  if (!is.null(panel)) {
    counts$group <- counts[,names(panel)[1]] 
    legendpos <- 'bottom'
  } else {
    counts$group <- 'All'
    legendpos <- 'none'
  }
  
  if (!is.null(panel)) {
    if (length(panel)==2) {
      variables2 <- c('timeDx', panel[2])
      counts2 <- ddply(testhist, variables2, function(x) nrow(x))
      counts2$group <- counts2[,names(panel)[2]]
      counts$biggroup <- names(panel)[1]
      counts2$biggroup <- names(panel)[2]
      counts <- rbind(counts2[,c('timeDx','V1','group', 'biggroup')], 
                         counts[,c('timeDx','V1','group', 'biggroup')])
    }
  }
  
  p <- ggplot(counts,aes(x=timeDx,y=V1,group=group))  +   
    geom_line(aes(color=group)) +
    geom_point(aes(color=group)) + 
    theme_bw()+
    scale_color_hue(name="") +
    theme(legend.position=legendpos,axis.text.x=element_text(angle=90)) + 
    scale_x_continuous(breaks=seq(min(testhist$yearDx),max(testhist$yearDx),by=1))+
    xlab("Time") + ylab("Diagnoses") 
  
  if (!is.null(panel)) {
    if (length(panel)==2) p <- p +facet_grid(.~biggroup) 
  }
  
  return(p)
}

######################################################################
# PLOT everHadNegTest over time in percents
######################################################################
#' Plot responses to 'Have you ever had a negative test?' over time
#'
#' Plots the everHadNegTest results over time
#'
#' @param testhist Data frame of class 'testinghistories' with variables
#'        everHadNegTest and yearDx
#' @param panel Optional variable by which to panel the plot
plotTestHist <- function(testhist, panel=NULL) {

    tabTime <- tabTestHist(testhist, 'yearDx')

    keep.vars <- c(panel, 'yearDx', grep('Percent ', colnames(tabTime),
                                       value=TRUE))
    these.ids <- c(panel, 'yearDx')
    tabTime <- tabTime[,keep.vars]
    colnames(tabTime) <- gsub('Percent ','',colnames(tabTime))
    tabTime <- melt(tabTime, id.vars=these.ids)
    if (!is.null(panel)) tabTime$Group <- tabTime[,panel]

    p <- ggplot(tabTime,aes(x=yearDx,y=value,group=variable))  +   
    geom_line(aes(color=variable)) +
    geom_point(aes(color=variable)) + 
    theme_bw()+
    theme(legend.position='bottom',axis.text.x=element_text(angle=90)) + 
    scale_color_hue(name="Ever had negative test?") + 
    scale_x_continuous(breaks=seq(min(tabTime$yearDx),max(tabTime$yearDx),by=2))+
    xlab("Time") + ylab("Percent") 

    if (!is.null(panel)) p <- p + facet_grid(.~Group)
    return(p)
}

######################################################################
# runBackCalc
######################################################################

#' Optional wrapper function to run all the backcalculation steps
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. 
#' @param prev  Optional data frame with 1st column 'Year' and a 2nd column with
#'        PLWHA for the population represented in the testhist object
runBackCalc = function(testhist, intLength, cases=NULL, prev=NULL) {
  
  if (!is.null(cases)) stop("In runBackCalc: not coded yet")
  
  # Estimate TIDs
  TIDs <- estimateTID(testhist$infPeriod, intLength=diagInterval)
  cases <- names(TIDs)
  
  # Diagnoses
  diagCounts = tabulateDiagnoses(testhist, intLength=diagInterval)
  
  # Initialize incidence and undiagnosed count lists
  incidence <- vector(mode='list', length=length(cases))
  names(incidence) <- cases
  undiagnosed <- incidence
  
  # Estimate incidence and undiagnosed
  for (c in cases) {
    cat('\nEstimating case', c, '...\n')
    incidence[[c]] = estimateIncidence(y=diagCounts,
                                      pid=TIDs[[c]]$pdffxn,
                                      gamma=0.1,
                                      verbose=FALSE)
    undiagnosed[[c]] <- estimateUndiagnosed(incidence[[c]])
  }
  
  # Compile results
  results <- combineResults(list(`Base Case`=list(incidence[['base_case']],
                                                  undiagnosed[['base_case']]),
                                 `Upper Bound`=list(incidence[['upper_bound']],
                                                    undiagnosed[['upper_bound']])))
  
  # True prevalence
  if (!is.null(prev)) trueprev <- calcTruePrev(results, prev)
  
  return(list(TIDs=TIDs, results=results, trueprev=trueprev, N=nrow(testhist)))
}

######################################################################
# runSubgroups
######################################################################

#' Optional wrapper function to run and compile results for subgroups
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @subvar Name of the variable within the testhist data frame that defines
#'        subgroups within which to run the backcalculation
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. 
#' @prev  Optional frame with 1st column 'Year' and a 2nd column with PLWHA
#'        for the population represented in the testhist object
#' @save  Optional file path to save compiled true prevalence results
runSubgroups = function(testhist, subvar, intLength, cases=NULL, 
                        prev=NULL, save=NULL) {
  
    print(colnames(testhist))
  if (!is.null(cases)) stop("In runSubgroups: not coded yet")
  if (is.numeric(testhist[,subvar])) {
    warning('Subgroup variable will be coerced to character')
    testhist[,subvar] <- as.character(testhist[,subvar])
  }
  
  # Check that if prevalence is given, it is given for all the 
  # subgroups
  subgroups <- unique(testhist[,subvar])
  numsub <- length(subgroups)
  if (sum(subgroups%in%colnames(prev))!=numsub) stop('In runSubGroups, 
                          prevalence data are insufficient')
  
  # Prepare to store results for each subgroup
  subResults <- vector('list', length=(numsub+1))
  names(subResults) <- c(as.character(subgroups), 'Total-stratified')
  
  # Loop through subgroups
  for (s in subgroups) {
    
    cat('\nSUBGROUP: ', s, '\n')
    # Run the backcalculation for this subgroup, selecting
    # the correct prevalence column if applicable
    if (!is.null(prev)) {
      subPrev <- prev[, c('Year', as.character(s))]
    } else subPrev <- NULL
    subResults[[s]] <- runBackCalc(testhist[testhist[,subvar]==s,], 
                                   intLength=diagInterval, 
                                   prev=subPrev)
    
      # Add a subgroup identifier to the compiled results
      for (r in c('resultsAll', 'resultsSummary', 'resultsSummaryYear')) {
        subResults[[s]]$results[[r]] = data.frame(Subgroup = s,
                                                  subResults[[s]]$results[[r]],
                                                  check.names=FALSE)
      }
  }
  
  # Extract the results in order to get subgroup-stratified totals
  resultsAllList <- lapply(subResults, function(x) x$results$resultsAll$value)
  resultsAll <- cbind(subResults[[1]]$results$resultsAll[,c('time', 'group', 'var')],
                      do.call(cbind, resultsAllList))
  resultsAll$value <- apply(resultsAll[,as.character(subgroups)],1,sum)
  
  # Summarize subgroup-stratified totals
  # Data frame with summarized results
  resultsSummary <- ddply(resultsAll, .(var, group), function(x) c(summary(x$value)))
  resultsSummary <- within(resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  resultsSummaryYear <- ddply(transform(resultsAll, Year=floor(time)), 
                                .(var, group, Year), function(x) c(summary(x$value)))
  resultsSummaryYear <- within(resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Create the results list of class 'results'
    resultsL=list(resultsAll=resultsAll,
                      resultsSummary=resultsSummary,
                      resultsSummaryYear=resultsSummaryYear)
    class(resultsL) <- append(class(x), 'results')

  # Save in subResults[['Total-stratified']]$results object
  subResults[['Total-stratified']] <- list(results=resultsL)
  
  if (!is.null(prev)) {
    
    # Calculate total-stratified true prevalence
    subResults[['Total-stratified']]$trueprev <- 
      calcTruePrev(subResults[['Total-stratified']]$results,
                   prev=data.frame(Year=prev$Year,
                                   Total=apply(trueprev_data[,as.character(subgroups)],
                                               1,sum)))
  }
  
  if (!is.null(save)) {
    trueprev <- do.call(rbind, lapply(names(subResults), 
                               function(x) {
                                 data.frame(Subgroup=x,
                                            subResults[[x]]$trueprev,
                                            check.names=FALSE)
                                }))
    write.csv(trueprev[, !colnames(trueprev) %in% c('1st Qu.', 'Median', '3rd Qu.')],
              file=save,
              row.names=FALSE)
  }
  
  return(subResults)
}

######################################################################
# calcTruePrev
######################################################################

#' Function to combine yearly PLWHA prevalence with undiagnosed estimates to provide
#' true prevalence estimates
#' @param x Object of class 'results', the output of combineResults
#' @param prev Data frame with 1st column 'Year' and a 2nd column with PLWH
#'        for the population represented in the results
calcTruePrev = function(x, prev) {
  
  # Fix the column name of the prevalence estimate as 'PLWHA'
  colnames(prev)[2] <- 'PLWHA'
  
  # Subset to Undiagnosed estimates and merge on prevalence
  undiag <- subset(x$resultsSummaryYear, Estimate=='Undiagnosed Cases')
  undiag <- merge(undiag, prev, all.y=TRUE, by='Year')
  
  # Now that we've reduced estimates to those years for which we 
  # were given prevalence, lose the PLWHA column and extract the 
  # estimates and PLWHA as a matrix separately
  estMatrix <- undiag[,4:9]
  prevMatrix <- replicate(ncol(estMatrix), undiag$PLWHA)
  undiag$PLWHA <- NULL
  
  # Compute true prevalence and % undiagnosed
  trueprevMatrix <- estMatrix + prevMatrix
  undiagFracMatrix <- 100*(estMatrix/trueprevMatrix)
  
  # Turn those back into data frames
  trueprev <- undiag
  trueprev$Estimate <- 'True Prevalence'
  trueprev[,4:9] <- trueprevMatrix
  undiagFrac <- undiag
  undiagFrac$Estimate <- 'Undiagnosed Fraction (%)'
  undiagFrac[,4:9] <- undiagFracMatrix
  
  # Prepare prevalence to be able to insert it into the trueprev table
  previnsert = transform(prev, diag='PLWHA', est='PLWHA', min=NA, q1=NA, med=NA, q3=NA, max=NA)
  previnsert = previnsert[,c('Year', 'diag', 'est', 'min', 'q1', 'med', 
                             'PLWHA', 'q3', 'max')]
  colnames(previnsert) <- colnames(undiag)
  
  # Combine and sort
  allResults <- rbind(previnsert, undiag, trueprev, undiagFrac)
  allResults[,4:9] <- round(allResults[,4:9], 1)
  allResults <- allResults[order(allResults$Year, allResults[,2]),]
  
  return(allResults)
}

######################################################################
# combineResults
######################################################################

#' Create an object of class "results" that contains all results
#'  
#' Incidence and undiagnosed estimates for different cases are initially
#' saved in separate objects. This function combines them into one 
#' object of class "results" to facilitate presenting results. It also
#' summarizes results across all time periods and by year.
#'  
#' @param x List with two tiers: the first tier identifies the cases.
#'          Each case is a list of 2: the first element is the "backproj"
#'          object returned by estimateIncidence(), and the second is
#'          the vector of undiagnosed counts returned by estimateUndiagnosed()
#'  
#' @return List object of class "results" 
combineResults <- function(x) {
  # Times with observed diagnoses
  allTimes <- as.numeric(names(x[[1]][[1]]$y))
  obsTimes <- !is.na(allTimes)
  times <- allTimes[obsTimes]
  x$times <- times
  
  # Diagnoses
  x$diagnoses <- x[[1]][[1]]$y[obsTimes]
  
  # Incidence and Undiagnosed organized by case
  cases <- names(x)[!names(x)%in%c('times', 'diagnoses')]
  for (c in cases) {
    incidence <- x[[c]][[1]]$lambda[obsTimes]
    undiagnosed <- x[[c]][[2]][obsTimes]
    x[[c]] <- list(incidence=incidence, undiagnosed=undiagnosed)
  }
  
  # Data frame with all results
  x$resultsAll <- data.frame(time=times, 
                             group='Diagnoses and Incidence', 
                             var='# Diagnosed',
                             value=x$diagnoses)
  for (c in cases) {
    x$resultsAll  <- rbind(x$resultsAll,
                           data.frame(time=times, 
                                      group='Diagnoses and Incidence', 
                                      var=c,
                                      value=x[[c]]$incidence),
                           data.frame(time=times,
                                      group='Undiagnosed Cases',
                                      var=c,
                                      value=x[[c]]$undiagnosed))
  }
  
  # Data frame with summarized results
  x$resultsSummary <- ddply(x$resultsAll, .(var, group), function(x) c(summary(x$value)))
  x$resultsSummary <- within(x$resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(x$resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  x$resultsSummaryYear <- ddply(transform(x$resultsAll, Year=floor(time)), 
                            .(var, group, Year), function(x) c(summary(x$value)))
  x$resultsSummaryYear <- within(x$resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(x$resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  class(x) <- append(class(x), 'results')
  return(x)
}

#############################################################
# format_eHARS
#############################################################

#' Format eHARS person view data for use with the HIVBackCalc package
#'  
#' Creates a data frame of class 'testinghistories' with formatted and
#' calculated variables. Additionally returns a record of formatting
#' decisions for inconsistent entries
#'  
#' @param rawdata Unformatted data frame or file path to CSV file
#' @param assumptionNo Choice of assumption for those reporting never 
#'          having had a negative test. Default is 'age16', which imputes
#'          time from infection to diagnosis as min(age-16, 18 yrs). 'age16mid'
#'          instead uses the midpoint between age and 16, i.e. 
#'          min((age-16)/2, 18 yrs). 
#' @return Data frame of class "testinghistories" 
format_eHARS <- function(rawdata, assumptionNo='age16') {
  
    require(plyr)

    # To delete
    #rawdata  <- '/Users/jeanette/Dropbox/School/PhD/HIV_WA/data/eHARS_testsubset.csv'
    
    # Read in data
    if (!is.data.frame(rawdata)) {
        rawdata  <- read.csv(rawdata,
                             na.strings='',
                             stringsAsFactors=FALSE)
    }

    # Helper function to flag records for cleaning
    recordFlag  <- function(df, logical, message) {
        # Set NA's in logical to FALSE
        logical[is.na(logical)] <- FALSE
        # Create flag variable if it doesn't exist
        if (!'flag'%in%colnames(df)) df$flag=''
        # Record flag and tidy up
        df$flag[logical] <- paste0(df$flag[logical], '; ', message) 
        df$flag <- gsub('^;', '', df$flag)
        # Add to assumptions table
        assumptionsN[assumptionsN$Issue==message,'N'] <<- sum(logical) 
        return(df)
    }

    # Construct empty assumptions table
    assumptionsN <- data.frame(Issue=c('Missing month',
                                    'Missing day',
                                    'Illogical last negative',
                                    'everHadNegTest inconsistent with infPeriod',
                                    'infPeriod capped at aidsUB',
                                    'Age <=16 and no infPeriod'),
                               Assumption=c('Month (diagnosis or last neg test) assumed July for computing infPeriod; diagnosis quarter randomly imputed',
                                            'Day (diagnosis or last neg test) assumed 15th for computing infPeriod',
                                            'Last negative date overwritten as missing because recorded as at or after diagnosis',
                                            'everHadNegTest flag altered to match presence/absence of last neg test date',
                                            'infPeriod capped at 18 years',
                                            'Removed from dataset because <=16 yrs and no recorded infPeriod'),
                               N=NA)

    #############################################################
    # OVERVIEW: N, expected variables and labels for factors
    #############################################################
    dataf <- rawdata
    colnames(dataf) <- tolower(colnames(dataf))
    variables_expected <- c('rsh_state_cd',
                          'aids_dx_dt',
                          'cd4_first_hiv_dt',
                          'cd4_first_hiv_type',
                          'cd4_first_hiv_value',
                          'hiv_aids_age_yrs',
                          'hiv_dx_dt',
                          'race',
                          'screen_last_neg_dt',
                          'trans_categ',
                          'vl_first_det_dt',
                          'vl_first_det_value',
                          'birth_sex',
                          'stateno',
                          'tth_last_neg_dt',
                          'tth_first_pos_dt',
                          'tth_ever_neg')
    not_in_dataf <- variables_expected[!variables_expected%in%colnames(dataf)]
    if (length(not_in_dataf)!=0) stop('Some eHARS variables are missing: \n', 
                                      paste(not_in_dataf,
                                            collapse='\n'))
    
    # Factors and labels
    races <- c('Hispanic',
                'American Indian/Alaska Native',
                'Asian',
                'Black',
                'Native Hawaiian/Pacific Islander',
                'White',
                'Legacy Asian/Pacific Islander',
                'Multi-race',
                'Unknown')
    modes <- c('Adult MSM', 
                'Adult IDU',
                'Adult MSM & IDU',
                'Adult received clotting factor',
                'Adult heterosexual contact',
                'Adult received transfusion/transplant',
                'Perinatal exposure, HIV diagnosed at age 13 years or older',
                'Adult with other confirmed risk',
                'Adult with no identified risk (NIR)',
                'Adult with no reported risk (NRR)',
                'Child received clotting factor',
                'Perinatal exposure',
                'Child received transfusion/transplant',
                'Child with other confirmed risk',
                'Child with no identified risk (NIR)',
                'Child with no reported risk (NRR)',
                'Risk factors selected with no age at diagnosis')
    dataf <- transform(dataf,
                     new_race=as.character(factor(race, levels=1:9,
                                     labels=races)),
                     new_mode=as.character(factor(trans_categ, 
                                                  levels=c(1:13,18:20,99), 
                                                  labels=modes)),
                       stringsAsFactors=FALSE)
    # Some renaming
    dataf <- rename(dataf,c('hiv_aids_age_yrs'='hdx_age',
                              'birth_sex'='sex'))
    

    #############################################################
    # SUMMARIZE: Summarize or tabulate variables and store in
    #            a table
    #############################################################
    varsummaries <- vector('list', length=ncol(dataf))
        names(varsummaries) <- colnames(dataf)
        for (x in 1:ncol(dataf)) {
        varname <- colnames(dataf)[x]
        var = dataf[,x]
        if (length(unique(var))>25) {
          if (is.numeric(var)) {
            sum <- summary(var)            
            if ("NA's"%in%names(sum)) {
                nmiss <- sum["NA's"]
            } else nmiss <- 0
            sum <- sum[c('Min.', 'Mean', 'Max.')]
          } else {
            nmiss <- table(var)[is.na(names(table(var)))]
            sum <- ''
          }
        } else {
          sum <- table(var, useNA='always')
          nmiss <- sum[is.na(names(sum))]
        }
        nmiss <- as.numeric(nmiss)
        result <- c(sum, round(100*nmiss/nrow(dataf),2))
        names(result)[length(result)] <- 'Percent Miss'
        result <- result[!is.na(names(result))]
        varsummaries[[x]] <- data.frame(Variable=c(varname, 
                                                   rep('', length(result)-1)),
                                        Values=names(result),
                                        N=result,
                                        row.names=NULL)
    }
    varsummaries <- data.frame(do.call('rbind', varsummaries),
                               row.names=NULL)


    #############################################################
    # COLLAPSE RACE AND MODE OF DIAGNOSIS
    #############################################################

    collapsed_race <- c('White', 'Black', 'Hispanic', 'Asian', 
                        'Native', 'Multi/Other')
    collapsed_mode <- c('MSM', 'Hetero', 'Blood/Needle/Other')
    dataf <- within(dataf, {
        race6 <- as.character(new_race)
        race6[race6 %in% c("American Indian/Alaska Native", 
                           "Native Hawaiian/Pacific Islander")] <- 'Native'
        race6[race6 %in% "Legacy Asian/Pacific Islander"] <- 'Asian'
        race6[race6 %in% c("Multi-race","Unknown")] <- 'Multi/Other'
        mode3 <- as.character(new_mode)
        mode3[mode3 %in% c('Adult MSM','Adult MSM & IDU')] <- 'MSM'
        mode3[mode3 %in% c('Adult heterosexual contact',
                           'Adult with no identified risk (NIR)')] <- 'Hetero'
        mode3[!mode3 %in% c('MSM', 'Hetero')] <- 'Blood/Needle/Other'
#        race6 <- factor(race6,
#                       labels=collapsed_race,
#                       levels=collapsed_race)
#        mode3 <- factor(mode3,
#                       levels=collapsed_mode,
#                       labels=collapsed_mode)
#        mode2 <- factor(ifelse(mode3 %in% 'MSM', 'MSM', 'non-MSM'))
         mode2 <- ifelse(mode3 %in% 'MSM', 'MSM', 'non-MSM')
        # FOR NOW: make the main mode=mode2
        mode <- mode2
    })

    #############################################################
    # DEFINE AGE GROUPS
    #############################################################
    dataf <- transform(dataf,
                     agecat5=cut(hdx_age,
                                 breaks=c(0,seq(20,70,by=5),85),
                                 include.lowest=TRUE,
                                 right=TRUE,
                                 labels=c('<=20',
                                          '21-25',
                                          '26-30',
                                          '31-35',
                                          '36-40',
                                          '41-45',
                                          '46-50',
                                          '51-55',
                                          '56-60',
                                          '61-65',
                                          '66-70',
                                          '71-85')))

    #############################################################
    # FORMAT DATES AND CREATE INFPERIOD
    #############################################################

    # Helper function to work with dates
    # For each date, need a fake date if month and/or day are missing 
    # plus an imputed quarter if month is missing

    get_dates <- function(timevar) {
        year <- suppressWarnings(as.numeric(substr(timevar,1,4)))
        month <- substr(timevar,5,6)
        day <- substr(timevar,7,8)
        missing_month <- !is.na(year) & month=='..'
        missing_day <- !is.na(year) & day=='..' & !month=='..'
        # Create a year-quarter variable, imputing a quarter if necessary
        set.seed(98103)
        yrqtr <- year + 
            suppressWarnings(ifelse(missing_month, sample(c(0,0.25,0.5,0.75)), 
                                    floor(as.numeric(month)/4)*0.25))
        # Create an  _imputed date for calculating inter-test intervals
        # 15th of the month if only month is known; July 1 if only year known
        day <- ifelse(missing_month, '01', ifelse(missing_day, '15', day))
        month <- ifelse(missing_month, '07', month)
        dateChar <- apply(cbind(year,month,day),1,paste,collapse='-')
        dateChar[dateChar=='NA-NA-NA'] <- ''
        dateImp <- as.Date(dateChar,"%Y-%m-%d")
        return(list(dateImp=dateImp, 
                    year=year,
                    yrqtr=yrqtr,
                    missMonth=missing_month, 
                    missDay=missing_day))
    }

    # Diagnosis date
    dxDate <- get_dates(dataf$hiv_dx_dt)
    # Last negative test date
    negDate <- get_dates(dataf$tth_last_neg_dt)
    dataf <- transform(dataf,
                       yearDx=dxDate$year,
                       timeDx=dxDate$yrqtr,
                       infPeriod=as.numeric(dxDate$dateImp-
                                            negDate$dateImp)/365,
                       stringsAsFactors=FALSE)
    # Record assumptions
    dataf <- recordFlag(dataf, dxDate$missMonth | negDate$missMonth,
                        'Missing month')
    dataf <- recordFlag(dataf, dxDate$missDay | negDate$missDay,
                        'Missing day')

    # Illogical last negative
    dataf <- recordFlag(dataf, dataf$infPeriod<=0,
                        'Illogical last negative')
    dataf <- within(dataf, {
                    infPeriod[infPeriod<=0] <- NA
                       })


    #############################################################
    # CREATE everHadNegTest after creating infPeriod
    #############################################################
    # Define everHadNegTest based on tth_ever_neg
    dataf <- transform(dataf, 
                     everHadNegTest=ifelse(tth_ever_neg=='Y', TRUE, 
                                           ifelse(tth_ever_neg=='N', FALSE, NA)))
    #with(dataf,table(everHadNegTest, tth_ever_neg, useNA='always'))

    # Look at actual infPeriod values by everHadNegTest
    #ddply(dataf, .(everHadNegTest), function(x) c(summary(x$infPeriod)))

    ## ---- fix_everHadNegTest_toTRUE ----
    toTRUE1 <- !dataf$everHadNegTest & !is.na(dataf$infPeriod)
    toTRUE2 <- is.na(dataf$everHadNegTest) & !is.na(dataf$infPeriod)
    dataf$everHadNegTest[toTRUE1] <- TRUE
    dataf$everHadNegTest[toTRUE2] <- TRUE

    ## ---- fix_everHadNegTest_toFALSE ----
    toFALSE <- dataf$everHadNegTest & is.na(dataf$infPeriod)
    dataf$everHadNegTest[toFALSE] <- FALSE

    # Record assumptions and flag
    dataf <- recordFlag(dataf, toTRUE1 | toTRUE2 | toFALSE,
                     'everHadNegTest inconsistent with infPeriod')
              
    ## ---- check_everHadNegTest ----
    #checkEver <- with(dataf,table(everHadNegTest, 
    #                             TID_NA=is.na(infPeriod), useNA='always')))

    #############################################################
    # EDIT infPeriod
    #############################################################

    # Cap at AIDS upper bound of ~18 years
    aidsUB <- qweibull(.95,shape=2.516,scale=1/0.086) #17.98418
    infTemp <- dataf$infPeriod
    
    fixNo <- function(dataf, assumptionNo) {
        switch(assumptionNo,
               age16 = {
                    dataf <- transform(dataf,
                                       infPeriod=ifelse(everHadNegTest, 
                                                        pmin(infPeriod, aidsUB), 
                                                        ifelse(!everHadNegTest, 
                                                               pmin(hdx_age-16, 
                                                                    aidsUB), 
                                                               NA)))
               },
               age16mid = {
                    dataf <- transform(dataf,
                                       infPeriod=ifelse(everHadNegTest, 
                                                        pmin(infPeriod, aidsUB), 
                                                        ifelse(!everHadNegTest, 
                                                               pmin((hdx_age-16)/2, 
                                                                    aidsUB), 
                                                               NA)))
               }
               )
        return(dataf)
    }
    dataf <- fixNo(dataf, assumptionNo)

    dataf <- recordFlag(dataf, infTemp!=dataf$infPeriod,
                        'infPeriod capped at aidsUB')

    # Remove cases who are too young for the impute-infPeriod assumption
    dataf <- recordFlag(dataf, 
                        dataf$hdx_age<=16 & (!dataf$everHadNegTest | 
                                        is.na(dataf$everHadNegTest)),
                        message='Age <=16 and no infPeriod')
    dataf <- subset(dataf, !(hdx_age<=16 & (!everHadNegTest | 
                                            is.na(everHadNegTest))))

    #############################################################
    # CREATE infPeriod_imputeNA
    #############################################################
    dataf <- within(dataf,{ 
        infPeriod_imputeNA=ifelse(is.na(everHadNegTest),
                                  pmin(hdx_age-16, aidsUB),
                                  infPeriod)
    })


    class(dataf) <- append(class(dataf), 'testinghistories')
    return(list(data=dataf,
                assumptions=assumptionsN,
                rawVarSum=varsummaries))
}

} # end if(1==0)

