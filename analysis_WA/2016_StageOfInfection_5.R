#############################################################
# SETUP
#############################################################
rm(list=ls())
# TEMPORARY: SOURCE FUNCTIONS
source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/internal_fxns.R')

# Change year min and max
year_min <- 2005
year_max <- 2014

# Load libraries, data and data-cleaning file
# Eventually this function should return the cleaned dataset,
# but data-cleaning has the name hardcoded as msm and I'm not
# going to generalize that right now
setup_hivbackcalc(workd='/Users/jeanette/Dropbox/School/PhD/HIV_WA',
                  datafile='data/wa_backcalc_data_201506.csv',
                  source_these='analysis_WA/format_data.R',
                  package_updated=TRUE,
                  packagefile='HIVBackCalc/R/internal_fxns.R')

library(xtable)
library(gridExtra)
#############################################################
# KNITR
#############################################################

library(knitr)
knit_hooks$set(inline = function(x) {
                    prettyNum(round(x,2), big.mark=",")
                  })
# set global chunk options
opts_chunk$set(fig.path='figure/minimal-', fig.align='center', 
               fig.show='hold', concordance=TRUE, external=TRUE,
               tidy=TRUE, size='footnotesize', verbose=TRUE,
               purl=TRUE)
options(formatR.arrow=TRUE,width=80,digits=7)

read_chunk(file.path(workd,'analysis_WA/run_main.R'))
#read_chunk(file.path(workd,'analysis_WA/describe_data.R'))

#############################################################
# ADDITIONAL FORMATTING FOR EXPLORING BED + DUAL DIAGNOSES
#############################################################
BEDlabels <- c('Long-standing infection', 'Recent infection',
            'Other/unknown')
BED.DD <- c(sapply(c('BED+', 'BED-', 'BEDm'), 
                   FUN=function(x) {paste0(x, c('DD+', 'DD-'))}))
           
dataf <- within(dataf, {
                est_infect_period=factor(est_infect_period,
                                         levels=1:3,
                                         labels=BEDlabels)
                BED=factor(est_infect_period,
                           levels=BEDlabels,
                           labels=c('-','+','Miss'))
                BED2=paste('BED', BED)
                BED2fac=factor(BED2, levels=c('BED +', 'BED -', 'BED Miss'),
                               labels=c('BED +', 'BED -', 'BED Miss'))
                AidsMinusHIVDx=timeAids-timeDx
                dualDiag=ifelse(AidsMinusHIVDx<1,TRUE,FALSE)
                dualDiag[is.na(timeAids)] <- FALSE
                dualDiagChar=ifelse(dualDiag, 'AIDS w/in 1yr of HIV',
                                    'no AIDS or not w/in 1yr of HIV')
                DD2=ifelse(dualDiag, 'DD +', 'DD -')
                DD2fac=factor(DD2, levels=c('DD +', 'DD -'),
                              labels=c('DD +', 'DD -'))
                stageGroup=NA
                stageGroup[BED=='+' & dualDiag] <- BED.DD[1]
                stageGroup[BED=='+' & !dualDiag] <- BED.DD[2]
                stageGroup[BED=='-' & dualDiag] <- BED.DD[3]
                stageGroup[BED=='-' & !dualDiag] <- BED.DD[4]
                stageGroup[BED=='Miss' & dualDiag] <- BED.DD[5]
                stageGroup[BED=='Miss' & !dualDiag] <- BED.DD[6]
                stageGroup=factor(stageGroup,
                                levels=BED.DD, labels=BED.DD)
            })

# Function to get the empirical cdf of LNT windows and compare 
# to the BED window

LNTcdf <- function(df, thisgroup='all', BEDw=162/365) {
    # Subset?
    if (thisgroup!='all') df <- subset(df, stageGroup==thisgroup)
    # Empirical cdf (the next line removes NAs)
    cdf <- ecdf(df$infPeriod)
    infs <- unique(sort(df$infPeriod))
    # Plot the survivor fxn
    cdfdf <- data.frame(infPeriod=infs,
                        cdf=cdf(infs),
                        surv=1-cdf(infs))

    BEDpoint <- data.frame(x=BEDw,y=1-cdf(BEDw))
    BEDlabel <- paste0(round(100*BEDpoint$y),'% have \nwindows>BED window')
    g <- ggplot(cdfdf, aes(x=infPeriod, y=surv)) + 
        geom_line() + 
        theme_bw() + 
        scale_x_continuous(name='Infection window length') + 
        scale_y_continuous(name='S(x)') +
        geom_text(data=BEDpoint, aes(x=x+4,y=y),label=BEDlabel) +
        geom_point(data=BEDpoint, aes(x=x,y=y), size=3) + 
        ggtitle(paste0('S(x) of infection windows\nfor ', thisgroup, ' cases'))
    
    return(list(g=g, BEDsurv=round(100*BEDpoint$y)))
}



# Proportions: total, and conditional on rows or columns
percTab <- function(proptab) { round(100*proptab) }

ctab <- with(dataf, table(BEDResult=BED2fac, dualDiagChar, useNA='ifany'))
ptab <- percTab(prop.table(ctab))

tabMain <- function(df) {
    # Crosstab
    ctab <- with(df, table(BEDResult=BED2fac, dualDiagChar, useNA='ifany'))

    # Not using all, right now
    ptab <- percTab(prop.table(ctab))
    ptab1 <- percTab(prop.table(ctab, margin=1))
    ptab2 <- percTab(prop.table(ctab, margin=2))

    # Margins
    m1 <- table(df$dualDiagChar)
    m2 <- table(df$BED2fac)
    pm1 <- percTab(prop.table(m1))
    pm2 <- percTab(prop.table(m2))

    # Combine N and row percents and margins
    cptab1 <- data.frame(ctab[,1], ptab1[,1], ctab[,2], ptab1[,2], 
                         c(m2), c(pm2),
                         stringsAsFactors=FALSE)
    cptab1 <- rbind(cptab1, c(m1[1], pm1[1], m1[2], pm1[2], sum(m2), 100))
    rownames(cptab1)[4] <- 'Total'
    colnames(cptab1) <- c('N', 'Row %', 'N', 'Row %', 'N', 'Col %')

    return(cptab1)
}

cptab <- tabMain(dataf)
cptab2014 <- tabMain(subset(dataf,yearDx==2014))
cptabMSM <- tabMain(subset(dataf,mode2=='MSM'))

# EverHadNegTest by the groups
variables <- 
    c(`BED Result`='BED2fac')
variables2 <- 
    list(c(`Transmission Mode`='mode2',`BED Result`='BED2fac'))
testhistTab <- tabTestHist(dataf,
                                                      variables,
                                                      supercolumn=TRUE,
                                                      fullsample_row=TRUE)
testhistTab2 <- tabTestHist(dataf,
                                                      variables2,
                                                      supercolumn=TRUE,
                                                      fullsample_row=TRUE)
# rownames(testhistTab) <- testhistTab$stageGroup
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline & \\multicolumn{2}{c}{DD+} & \\multicolumn{2}{c}{DD-} & \\multicolumn{2}{|c}{Total} \\\\'
print(xtable(cptab,
             caption='Cross-tabulation of BED and dual diagnoses, 2005-2014. Row percents show the percent of dual diagnoses (DD) within each BED group or the total sample (last row)',
             label='tab:cptab',
             align=c('l', rep('r',4), '|r', 'r'),
             digits=0),
      table.placement='!h',
      caption.placement='top',
      add.to.row=addtorow,
      include.rownames=TRUE,
      hline.after=c(-1,0,3,4),
      size='small'
      )
print(xtable(cptab2014,
             caption='Cross-tabulation of BED and dual diagnoses, 2014 only',
             label='tab:cptab2014',
             align=c('l', rep('r',4), '|r', 'r'),
             digits=0),
      table.placement='!h',
      caption.placement='top',
      add.to.row=addtorow,
      include.rownames=TRUE,
      hline.after=c(-1,0,3,4),
      size='small'
      )
if (1==0) {
print(xtable(ptab,
             caption='Cross-tabulation of BED and dual diagnoses, \\%',
             label='tab:ptab',
             digits=0),
      table.placement='!h',
      caption.placement='top',
      include.rownames=TRUE,
      size='small'
      )
print(xtable(ptab1,
             caption='Dual diagnoses \\%, conditional on BED',
             label='tab:ptab1',
             digits=0),
      table.placement='!h',
      caption.placement='top',
      include.rownames=TRUE,
      size='small'
      )
print(xtable(ptab2,
             caption='BED \\%, conditional on dual diagnosis',
             label='tab:ptab2',
             digits=0),
      table.placement='!h',
      caption.placement='top',
      include.rownames=TRUE,
      size='small'
      )
}

    stageTab <- ddply(dataf, c('yearDx'), function(x, TN=nrow(dataf)) {
      n <- nrow(x)
      c(N=n,
        `Column Percent`=round(100*n/TN,0),
        `Percent BED+`=round(100*sum(x$BED2fac=='BED +', na.rm=TRUE)/n,0),
        `Percent BED-`=round(100*sum(!x$BED2fac=='BED -', na.rm=TRUE)/n,0),
        `Percent BED Miss`=round(100*sum(x$BED2fac=='BED Miss')/n,0),
        `Percent DD+`=round(100*sum(x$DD2=='DD +')/n,0))
    })

    sg <- ggplot(stageTab, aes(yearDx, `Percent BED Miss`)) + 
        stat_smooth(method='loess') + geom_point() +
        theme_bw() + scale_x_continuous(name='Year')

    plot(sg)

    dg <- ggplot(stageTab, aes(yearDx, `Percent DD+`)) + 
        stat_smooth(method='loess') + geom_point() +
        theme_bw() + scale_x_continuous(name='Year')

    plot(dg)
# Marginal
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline &&&& \\multicolumn{3}{c}{Ever Had a Negative Test}  \\\\'
print(xtable(testhistTab,
             caption='Testing history responses by BED status. Column \\% sums to 100 across all rows. Availability of testing history data within each subgroup level is shown as row percents',
             label='tab:sample',
             digits=0),
      table.placement='h',
      caption.placement='top',
      add.to.row=addtorow,
      include.rownames=FALSE,
      size='small',
      sanitize.text.function=function(str) {
          gsub('(\\.)*Percent(\\.)*', ' \\\\% ',str);
      })
addtorow$command <- '\\hline &&&& \\multicolumn{3}{c}{Ever Had a Negative Test}  \\\\'
print(xtable(testhistTab2,
             caption='Testing history responses by BED status, separately for MSM and non-MSM. Column \\% sums to 100 across all rows. Availability of testing history data within each subgroup level is shown as row percents',
             label='tab:sampleMSM',
             digits=0),
      table.placement='h',
      caption.placement='top',
      add.to.row=addtorow,
      include.rownames=FALSE,
      size='small',
      sanitize.text.function=function(str) {
          gsub('(\\.)*Percent(\\.)*', ' \\\\% ',str);
      })
    dataf <- within(dataf, {
                    LNT=factor(as.character(everHadNegTest),
                               levels=c('TRUE', 'FALSE', NA),
                               labels=c('LNT', 'no LNT'))
          })
    datafBED <- subset(dataf, BED2fac=='BED +')
    datafBEDlong <- rbind(datafBED, datafBED)
    datafBEDlong <- transform(datafBEDlong,
                              WindowType=c(rep('Original', nrow(datafBED)),
                                             rep('BED-Modified', nrow(datafBED))),
                              infWindow=c(datafBED$infPeriod, 
                                          datafBED$infPeriodBED))

    v <- ggplot(subset(dataf,!is.na(LNT)), 
                aes(factor(LNT, levels=rev(levels(LNT))), infPeriod))
    v <- v + geom_violin(width=1,  aes(fill=LNT), alpha=0.5) + 
        scale_y_log10(breaks=c(0.25,0.5,1,2,5,18), name='Infection Window, Years') + 
        geom_boxplot(width=0.05, outlier.size=0, 
                     aes(fill=LNT)) + 
        coord_flip() + 
        theme_bw() + 
#        facet_grid(BED2fac~.) + 
        scale_x_discrete(name='', labels=NULL, breaks=NULL) + 
        scale_fill_discrete(guide=FALSE)
#    b <- ggplot(subset(dataf,!is.na(LNT)), aes(LNT))
    dataf$all = 'Testing History'
    b <- ggplot(subset(dataf, !is.na(LNT)), aes(all))
#    b <- b + geom_bar(aes(y = (..count..)/sum(..count..), fill=LNT),
    b <- b + geom_bar(aes(fill=LNT), position='fill') + 
        coord_flip() +
#        facet_grid(BED2fac~.) +
#        scale_fill_discrete(guide=FALSE) + 
        theme_bw() + 
        scale_y_continuous(name='', expand=c(0,0)) +
        scale_x_discrete(name='', labels=NULL, breaks=NULL)
    b <- b + theme(legend.position='top', legend.title=element_text('')) + labs(fill='')
    grid.arrange(b,v, ncol=1, heights=c(1.5,4))

#    ggplot(datafBEDlong, aes(y=Window, x=Type)) + geom_boxplot()
    grid.arrange(b+facet_grid(.~mode2),v+facet_grid(.~mode2),
                 heights=c(1.5,4))
dataf <- within(dataf, {
                ChangeBED=(BED=='+' & (infPeriod>162/365 | is.na(infPeriod)))
      })
cBEDdf <- table(everHadNegTest=dataf$everHadNegTest, 
                changeBED=dataf$ChangeBED, useNA='ifany')
cBEDdf <- rbind(cBEDdf, colSums(cBEDdf))
rownames(cBEDdf) <- c('LNT', 'No LNT', 'Missing TH', 'Total')
colnames(cBEDdf) <- c('No Change', 'Changed by BED+')
cBEDdfp <- cbind(round(100*prop.table(as.matrix(cBEDdf),1),1),Total=rowSums(cBEDdf))
cBEDdfall <- cbind(cBEDdf[,'No Change'], cBEDdfp[,'No Change'],
                   cBEDdf[,'Changed by BED+'], cBEDdfp[,'Changed by BED+'],
                   cBEDdfp[,'Total'])
colnames(cBEDdfall) <- c(rep(c('N', 'Row %'),2), 'Total')

# Print it
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline & \\multicolumn{2}{c}{No Change} & \\multicolumn{2}{c}{Changed by BED+} & \\\\'
print(xtable(cBEDdfall,
             caption='Impact of BED+ results on the entire sample, by testing history (TH) status. Cases have no change if they are not BED+ or they are BED+ but with a reported LNT that is shorter than the BED window of 162 days. Cases are impacted by BED+ results if they are BED+ with an LNT farther back than 162 days, they report no LNT, or they have missing testing history.',
             label='tab:changesummary',
             digits=0),
      caption.placement='top',
      table.placement='ht',
      add.to.row=addtorow,
      hline.after=c(-1,0,3,4),
      size='small',
      include.rownames=TRUE)
    dataf <- transform(dataf,
                       infPeriodBED=ifelse(BED2fac=='BED +' & (is.na(infPeriod) |
                                           infPeriod>(162/365)), 
                                           162/365,
                                           infPeriod))

    sumTable <- rbind(summary(subset(dataf, BED2fac=='BED +')$infPeriod), 
                      c(summary(subset(dataf, BED2fac=='BED +')$infPeriodBED), 0))
    rownames(sumTable) <- c('Original BED+ infection windows', 'After modifying windows for BED+')

    print(xtable(sumTable,
                 caption='Among the BED+, summary of infection windows pre- and post- modification by BED+, in years',
                 label='tab:sumTable',
                 digits=3),
          caption.placement='top',
          table.placement='ht',
          size='small',
          include.rownames=TRUE)
    datafBED <- subset(dataf, BED2fac=='BED +')
    BEDpTID <- estimateTID(datafBED$infPeriod,
                           intLength=0.25, 
                           stageGroup='BED+',
                           cases=c('base_case','base_case_withStage'))
    BEDpTIDu <- estimateTID(datafBED$infPeriod,
                           intLength=0.25, 
                           stageGroup='BED+',
                           cases=c('upper_bound','upper_bound_withStage'))
    Sx_tab <- cbind(summary(BEDpTID, times=c(0,0.25,0.5, 1), intLength=0.25)[,c(1,3,5)],
                    summary(BEDpTIDu, times=c(0,0.25,0.5, 1), intLength=0.25)[,c(3,5)])
    colnames(Sx_tab) <- c('Time', 'Original Base Case', 'Modified Base Case', 
                          'Original Upper Bound', 'Modified Upper Bound')
    print(xtable(Sx_tab,
                 caption='Among the BED+, the fraction remaining undiagnosed at quarter-year time steps, using the different cases. Time indicates the left/lower bound of the time interval',
                 label='tab:Sx_tab',
                 digits=3),
          caption.placement='top',
          table.placement='ht',
          size='small',
          include.rownames=FALSE)
    BED.DD4 <- c('BED+', 'BED-', 'BEDmDD+', 'BEDmDD-')
               
    dataf2 <- within(dataf, {
                    stageGroup <- as.character(stageGroup)
                    stageGroup[BED=='+'] <- BED.DD4[1]
                    stageGroup[BED=='-'] <- BED.DD4[2]
                    stageGroup[BED=='Miss' & dualDiag] <- BED.DD4[3]
                    stageGroup[BED=='Miss' & !dualDiag] <- BED.DD4[4]
                    stageGroup=factor(stageGroup,
                                    levels=BED.DD4, labels=BED.DD4)
                })

    table(dataf2$stageGroup, dataf2$BED)
    table(dataf2$stageGroup, dataf2$dualDiag)
diagInterval = 0.25
TIDs <- estimateTID(dataf$infPeriod, intLength=diagInterval)

# Original 6
theseGroups <- levels(dataf$stageGroup)
TIDsExt <- vector(mode='list',length=length(theseGroups))
for (i in 1:length(theseGroups)) {
    group  <- theseGroups[i]
    TIDsExt[[i]] <- plot(estimateTID(subset(dataf,stageGroup==group)$infPeriod,
                           intLength=diagInterval, 
                           stageGroup=group,
                           cases=c('base_case','base_case_withStage')),
                         intLength=0.25) + 
                    ggtitle(group) + 
                    theme(text=element_text(size=8))
                          
}
# New 4
theseGroups2 <- levels(dataf2$stageGroup)
TIDsExt2 <- vector(mode='list',length=length(theseGroups2))
for (i in 1:length(theseGroups2)) {
    group  <- theseGroups2[i]
    TIDsExt2[[i]] <- plot(estimateTID(subset(dataf2,stageGroup==group)$infPeriod,
                           intLength=diagInterval, 
                           stageGroup=group,
                           cases=c('base_case','base_case_withStage')),
                         intLength=0.25) + 
                    ggtitle(group) + 
                    theme(text=element_text(size=8))
}
# UB                          
group <- 'BED+'
    TIDsU <- plot(estimateTID(subset(dataf2,stageGroup==group)$infPeriod,
                           intLength=diagInterval, 
                           stageGroup=group,
                           cases=c('upper_bound','upper_bound_withStage')),
                         intLength=0.25) + 
                    ggtitle('BED+, Upper Bound') + 
                    theme(text=element_text(size=8))
#plot(baseVbase, intLength=diagInterval, 
#     cases = c('Original BC', 'Extended BC'))
#plot(TIDs, intLength=diagInterval, 
#     cases = c('Base Case', 'Upper Bound'))
grid.arrange(TIDsExt2[[1]],TIDsExt2[[2]],TIDsExt2[[3]],TIDsExt2[[4]],
             nrow=2)

TIDsU
  diagCounts = tabulateDiagnoses(dataf, intLength=diagInterval)
  incidenceBase = estimateIncidence(y=diagCounts,
                                    pid=TIDs[['base_case']]$pdffxn,
                                    gamma=0.1,
                                    verbose=FALSE)
# Here, do incidenceExt
  incidenceUpper = estimateIncidence(y=diagCounts,
                                    pid=TIDs[['upper_bound']]$pdffxn,
                                    gamma=0.1,
                                    verbose=FALSE)
  undiagnosedBase <- estimateUndiagnosed(incidenceBase)
  undiagnosedUpper <- estimateUndiagnosed(incidenceUpper)
# Here, do undiagnosedExt. Then,
# Here, combine results from the two base cases
#  results <- combineResults(list(`Original BC`=list(incidenceBase, 
#                                                    undiagnosedBase),
#                                 `Extended BC`=list(incidenceBaseExt,
#                                                    undiagnosedBaseExt))
  results <- combineResults(list(`Base Case`=list(incidenceBase,
                                              undiagnosedBase),
                               `Upper Bound`=list(incidenceUpper,
                                                undiagnosedUpper)))
plot(results)
print(xtable(results$resultsSummary,
             caption='Observed diagnoses and estimated quarterly incidence and undiagnosed counts over 2005-2014 in WA state',
             label='tab:res_main',
             digits=0),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)
  trueprev_data = read.csv(file.path(workd,'data/Reported_prevalence_2010-2014.csv'),
                           na.string="",stringsAsFactor=FALSE, check.names=FALSE)
  trueprev <- calcTruePrev(results, subset(trueprev_data, select=c('Year', 'Total')))
print(xtable(subset(trueprev, Year==2014),
             caption='Estimated true prevalence and the undiagnosed fraction in WA state, limited to just 2014',
             label='tab:res_trueprev',
             digits=1),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)

  these_cases <- c('base_case', 'base_case_withStage', 
                   'upper_bound', 'upper_bound_withStage')
  names(these_cases) <- c('Base Case', 'Base Case using Stage',
                          'Upper Bound', 'Upper Bound using Stage')
  subgroups <- runSubgroups(dataf2,
                            subvar='stageGroup',
                            intLength=diagInterval, 
                            cases=these_cases,
                            prev=trueprev_data,
                            save=file.path(workd, 'analysis_WA/results/2016_trueprev_StageOfInfection.csv'))

  totRes <- subgroups[['Total-stratified']]$results
  totTP <- subgroups[['Total-stratified']]$trueprev

  # Summary of summaries
  sumtable <- subset(totRes$resultsSummary[order(totRes$resultsSummary$Estimate),],
                     Estimate=='Undiagnosed Cases')
  sumtable2 <- subset(totRes$resultsSummaryYear[order(totRes$resultsSummaryYear$Estimate),],
                      Estimate=='Undiagnosed Cases' & Year==2014, select=-Year)
  sumtable$Year <- '2005-2014'
  sumtable2$Year <- '2014'
  sumtable <- rbind(sumtable, sumtable2)[,c('Year', colnames(sumtable)[-ncol(sumtable)])]
  colnames(sumtable)[2] <- 'Case'
  m <- subset(melt(sumtable), variable=='Mean', select=-Estimate)
  m <- transform(m,
                 variable=NULL,
                 Stage=rep(c('Without Stage', 'With Stage'),4),
                 Case=gsub(' using Stage', '', Case))
  mwide <- dcast(m, Year+Case~Stage, value.var='value')
  mwide$Difference <- mwide[,3]-mwide[,4]
  mwide <- transform(mwide, `Percent Change`=round(100*Difference/`Without Stage`),
                     check.names=FALSE)
plot(totRes)
print(xtable(sumtable,
             caption='Observed diagnoses and estimated quarterly incidence and undiagnosed counts over 2005-2014 and just 2014 in WA state, using stage-subgroup strata',
             label='tab:res_main2',
             digits=0),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)
print(xtable(mwide,
             caption='Impact of using BED result to modify the TID on mean undiagnosed estimates',
             label='tab:stagediff',
             digits=1),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)
print(xtable(subset(totTP, Year==2014)[order(subset(totTP, Year==2014)$Estimate),],
             caption='Estimated true prevalence and the undiagnosed fraction for 2014 in WA state, using stage-subgroup strata',
             label='tab:res_trueprev2',
             digits=1),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)
