
standalone <- TRUE
runEstimation <- FALSE

if (standalone) {
    # rm(list=ls())
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
                      datafile='data/wa_backcalc_data_201602.csv',
                      source_these='analysis_WA/format_data.R',
                      load_package=TRUE,
                      packagefile='HIVBackCalc/R/internal_fxns.R')

    library(xtable)
    library(gridExtra)
    library(plyr)
    library(reshape2)
    library(ggplot2)
}

#############################################################
# For explaining CD4 Case methods - plotting functions etc

## ---- plotProbModel
plot.BC <- function(maxWindow,
                    xmax=18, 
                    ymax=NULL,
                    extra.ticks=NULL, t='',
                    shading=NULL,
                    ...){
  x <- seq(0,maxWindow,by=0.1)
  y <- rep(1/maxWindow,length.out=length(x))
  if (is.null(ymax)) ymax=2*max(y)
  plot(0:xmax, seq(0,ymax,length.out=xmax+1), type='n', 
       xlab='Infection window (years prior to dx)', 
       ylab='Hazard of infection', 
       bty='l', xaxt='n', yaxt='n',
       main=t)
  if (!is.null(extra.ticks)) xticks=c(0,extra.ticks,xmax) else xticks=c(0,xmax)
  axis(side=1, at=xticks)
  if (!is.null(shading)) {
    xcord=c(0,0,shading, shading)
    ycord=c(0,y[1],y[1],0)
    polygon(xcord,ycord,col='gray',border=FALSE)
    label=paste0(round(100*shading/maxWindow),"%")
    text(x=shading/2, y=1/(maxWindow*2), labels=label)
  }
  lines(x, y, type='l',xaxt='n', yaxt='n', col='red',lwd=3)
}

plot.CD4 <- function(...){
  args <- list(...)
  med <- args$median
  plot.BC(...)
  x1=c(0,med)
  x2=c(med,args$maxWindow)
  y1=rep(0.5/med,2)
  y2=rep(0.5/diff(x2),2)
  if ((args$shading.median)) {
    xcord=c(0,0,med, med)
    ycord=c(0,y1[1],y1[1],0)
    polygon(xcord,ycord,col='gray',border=FALSE)
    label="50%"
    text(x=med/2, y=args$ymax/10, labels=label)
  }
  lines(c(x1,x2),c(y1,y2),col='blue',lwd=3)
  x <- seq(0,args$maxWindow,by=0.1)
  y <- rep(1/args$maxWindow,length.out=length(x))
  lines(x, y, type='l',xaxt='n', yaxt='n', col='red',lwd=3)
}


#############################################################
# Data transformations for descriptives

## ---- transformations
BEDlabels <- c('Long-standing infection', 'Recent infection',
               'Other/unknown')
cd4breaks <- c(0,200,350,500,2000)
windowbreaks <- c(0,3,8,16,18)

dataf <- transform(dataf, 
                   everHad=ifelse(is.na(everHadNegTest), 'Missing', 
                                  ifelse(everHadNegTest, 'Tested', 
                                         'Never Tested')))
dataf <- transform(dataf, 
                   everHad=factor(everHad, levels=c('Tested', 
                                                    'Never Tested', 
                                                    'Missing'), 
                                  label=c('Tested', 'Never Tested', 
                                          'Missing')),
                   infPeriodLength=cut(infPeriod, breaks=c(0,1,2,5,18))
                   )

dataf <- within(dataf, {
  infType=as.character(infPeriodLength)
  infType=ifelse(is.na(infType), infType, paste(infType, 'yrs'))
  infType[!is.na(infPeriod) & infPeriod>17.9 & everHad=="Never Tested"] <- 
    "18 yrs (NT)"
  # Non-missing testing history
  hasTestHist <- !is.na(everHadNegTest)
  # CD4 measured within 30d
  cd4within30 <- hasTestHist & !is.na(cd4_days) & cd4_days<=30 &
    !is.na(firstcd4cnt)
  # Another CD4 within 30 that doesn't limit to non-missing TH
  cd4within302=ifelse(!is.na(firstcd4cnt) & !is.na(cd4_days) & 
                        cd4_days<=30, TRUE, FALSE)
  # Categories
  cd4cat=cut(firstcd4cnt, breaks=cd4breaks,
             include.lowest=TRUE, right=FALSE)
  origcd4cat=cd4cat
  cd4cat=as.character(cd4cat)
  cd4cat[!cd4within30]<-NA
  cd4cat[is.na(cd4cat)]<-'No CD4 within 30d'
  # Windows
  infBreaks=cut(infPeriod, breaks=windowbreaks,
                include.lowest=TRUE, right=FALSE)
  # Windows with 'yrs'
  infBreaksYrs=ifelse(is.na(infBreaks), infBreaks,
                      paste(infBreaks,'yrs'))
  # BED
  est_infect_period=factor(est_infect_period,
                           levels=1:3,
                           labels=BEDlabels)
  BED=factor(est_infect_period,
             levels=BEDlabels,
             labels=c('-','+','Miss'))
  BED2=paste('BED', BED)
  BED2fac=factor(BED2, levels=c('BED +', 'BED -', 'BED Miss'),
                 labels=c('BED +', 'BED -', 'BED Miss'))
  everHad2=ifelse(everHad=='Missing', 'TH Missing', 'TH Non-Missing')
  
})

#############################################################
# Testing history distributions (includes plots)

## ---- everHadDistr
everHadDistr<- with(dataf, table(mode2, everHad))
everHadDistr2 <- with(dataf, table(everHadNegTest))
wd <- data.frame(t(everHadDistr))
wd <- rbind(wd,
            data.frame(mode2='Total', 
                       ddply(wd, .(everHad), 
                             summarise, Freq=sum(Freq))))

# Figure advice: http://stackoverflow.com/questions/22231124/how-to-draw-stacked-bars-in-ggplot2-that-show-percentages-based-on-group

# Calculate the percentages
wd = ddply(wd, .(mode2),
           transform,
           Percent = round(100*Freq/sum(Freq)))

# Format the labels and calculate their positions
wd = ddply(wd, .(mode2), transform, pos = (cumsum(Freq) - 0.5 * Freq))
wd$label = paste0(sprintf("%.0f", wd$Percent), "%")

# Plot
everHadPlot <- ggplot(subset(wd, mode2=='Total'), aes(x = factor(mode2), y = Freq, fill = everHad)) +
  geom_bar(stat = "identity", width = .7) +
  geom_text(aes(y = pos, label = label), size = 4) +
  coord_flip() + theme(text = element_text(size=15)) + 
  scale_x_discrete(name='', label='', breaks='') + 
  scale_y_continuous('Number of cases') +
  geom_label(aes(y=pos, fill = everHad, label=everHad), colour = "white", fontface = "bold", vjust=-1) + theme_bw() + guides(fill=FALSE)

everHadPlotMSM <- 
  ggplot(subset(wd, mode2!='Total'), 
         aes(x = factor(mode2,levels=c('non-MSM', 'MSM'),
                        labels=c('non-MSM', 'MSM')), 
             y = Freq, fill = everHad)) +
  geom_bar(stat = "identity", width = .7) +
  geom_text(aes(y = pos, label = label), size = 4) +
  coord_flip() + theme(text = element_text(size=15)) + 
  scale_x_discrete(name='') + 
  scale_y_continuous('Number of cases') + theme_bw() +
  theme(legend.title=element_blank())

#############################################################
# Window  distributions (includes plots)

## ---- windowDistr
windowDistr <- with(dataf, table(mode2, infType))
windowDistr2 <- with(dataf, table(infType))
wd <- data.frame(t(windowDistr))
wd <- rbind(wd,
            data.frame(mode2='Total', 
                       ddply(wd, .(infType), 
                             summarise, Freq=sum(Freq))))

# Calculate the percentages
wd = ddply(wd, .(mode2),
           transform,
           Percent = round(100*Freq/sum(Freq)))

# Format the labels and calculate their positions
wd = ddply(wd, .(mode2), transform, pos = (cumsum(Freq) - 0.5 * Freq))
wd$label = paste0(sprintf("%.0f", wd$Percent), "%")

# Plot
wdPlot <- ggplot(subset(wd, mode2=='Total'), aes(x = factor(mode2), y = Freq, fill = infType)) +
  geom_bar(stat = "identity", width = .7) +
  geom_text(aes(y = pos, label = label), size = 4) +
  coord_flip() + theme(text = element_text(size=15)) + 
  scale_x_discrete(name='') + 
  scale_y_continuous('Number of cases') +
  geom_label(aes(y=pos, fill=infType, label=infType), 
             colour = "white", fontface = "bold", vjust=-1.5) + 
  theme_bw() + guides(fill=FALSE)

# Plot with a box around 0-2 years
wdBox1 <- data.frame(ymin=0,
                    ymax=sum(subset(wd, mode2=='Total')$Freq[1:2]),
                    xmin=-Inf, xmax=Inf)
wdPlotBox1 <- wdPlot + 
  geom_rect(data=wdBox1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color="grey20",
            alpha=0.2,
            inherit.aes = FALSE)

# Plot with a box around 5-18 years
wdBox <- data.frame(ymin=sum(subset(wd, mode2=='Total')$Freq[1:3]),
                    ymax=sum(subset(wd, mode2=='Total')$Freq),
                    xmin=-Inf, xmax=Inf)
wdPlotBox <- wdPlot + 
  geom_rect(data=wdBox, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
          color="grey20",
          alpha=0.2,
          inherit.aes = FALSE)

wdPlotMSM <- ggplot(subset(wd, mode2!='Total'), 
                    aes(x = factor(mode2,levels=c('non-MSM', 'MSM'),
                                   labels=c('non-MSM', 'MSM')), 
                        y = Freq, fill = infType)) +
#                    aes(x = factor(mode2), 
 #                       y = Freq, fill = infType)) +
  geom_bar(stat = "identity", width = .7) +
  geom_text(aes(y = pos, label = label), size = 2.75) +
  coord_flip() + theme(text = element_text(size=15)) + 
  scale_x_discrete(name='') + 
  scale_y_continuous('Number of cases') +
  labs(fill='') +
  theme_bw()

#############################################################
# Concurrent diagnoses (includes plots)

## ---- concurrentDx
dataf <- transform(dataf, aidsAtDx2=ifelse(aidsAtDx, 'Concurrent Dx', 'No concurrent Dx'))
aidsDistr <- with(dataf, table(mode2, aidsAtDx2))
aidsDistr2 <- with(dataf, table(aidsAtDx2))
wd <- data.frame(t(aidsDistr))
wd <- rbind(wd,
            data.frame(mode2='Total', 
                       ddply(wd, .(aidsAtDx2), 
                             summarise, Freq=sum(Freq))))
# Calculate the percentages
wd = ddply(wd, .(mode2),
           transform,
           Percent = round(100*Freq/sum(Freq)))

# Format the labels and calculate their positions
wd = ddply(wd, .(mode2), transform, pos = (cumsum(Freq) - 0.5 * Freq))
wd$label = paste0(sprintf("%.0f", wd$Percent), "%")

# Plot
concurrentDxPlot <- ggplot(subset(wd, mode2!='Total'), 
                           aes(x = factor(mode2,levels=c('non-MSM', 'MSM'),
                                          labels=c('non-MSM', 'MSM')), 
                           #aes(x = factor(mode2), 
                               y = Freq, fill = aidsAtDx2)) +
  geom_bar(stat = "identity", width = .7) +
  geom_text(aes(y = pos, label = label), size = 2.75) +
  coord_flip() + theme(text = element_text(size=15)) + 
  scale_x_discrete(name='') + 
  scale_y_continuous('Number of cases') +
  labs(fill='') +
  theme_bw() +
  scale_fill_brewer(palette="Spectral")

#############################################################
# Mean age (includes plot)

## ---- meanAge
meanAge <- ddply(subset(dataf, everHad=='Never Tested'), .(mode2), summarise, meanAge=mean(hdx_age))
# Plot
meanAgePlot <- ggplot(subset(dataf, everHad=='Never Tested'), 
                      aes(hdx_age, fill = mode2)) +
  geom_density(aes(x=hdx_age, group=mode2, 
                   colour=mode2, fill=mode2), fill=NA) +
  theme_bw() +
  scale_x_continuous(name='Age at diagnosis') + 
  scale_y_continuous(name='') +
  labs(fill='', colour='')



#############################################################
# CD4 info

## ---- cd4not30

df28 <- subset(dataf,!is.na(everHadNegTest) & !cd4within30)
df28 <- transform(df28, cd4_days_cat=cut(cd4_days,breaks=c(0,30,60,3500)))
cd4_days_cat_lev <- c(levels(df28$cd4_days_cat)[2:3],'No CD4')
df28 <- within(df28,{
  cd4_days_cat=as.character(cd4_days_cat)
  cd4_days_cat[is.na(cd4_days_cat)] <-'No CD4'
  cd4_days_cat=factor(cd4_days_cat,
                      levels=cd4_days_cat_lev,
                      labels=cd4_days_cat_lev)
})

cd4not30Plot <- ggplot(df28, aes(x=cd4_days_cat)) + 
    geom_bar(aes(fill=origcd4cat)) + theme_bw() +
    scale_x_discrete(name="# days after dx when CD4 was measured") +
    scale_y_continuous(name='# of cases') + labs(fill='')

# A second attempt, since the 1st graph does not communicate
# the main message...
df28 <- transform(df28, highCD4=ifelse(firstcd4cnt>=350, 
                                       'CD4>=350', 
                                       'CD4<350'))
hicd4distr <- with(df28, table(cd4_days_cat, highCD4, useNA='ifany'))
hd <- data.frame(t(hicd4distr))
# Calculate the percentages
hd = ddply(hd, .(cd4_days_cat),
           transform,
           Percent = round(100*Freq/sum(Freq)))

# Format the labels and calculate their positions
hd = ddply(hd, .(cd4_days_cat), transform, pos = (cumsum(Freq) - 0.5 * Freq))
hd$label = paste0(sprintf("%.0f", hd$Percent), "%")

# Plot
cd4not30Plot2 <- ggplot(subset(hd,Percent!=0),
                        aes(x = factor(cd4_days_cat), 
                                y = Freq, fill = highCD4)) +
  geom_bar(stat = "identity", width = .7) +
  geom_text(aes(y = pos, label = label), size = 2.75) +
  theme(text = element_text(size=15)) + 
  scale_x_discrete(name='') + 
  scale_y_continuous('Number of cases') +
  labs(fill='') +
  theme_bw() 
  #scale_fill_manual(values=c("#b2df8a", "#33a02c", "#000000"), name='')
  #scale_fill_brewer(palette="Set1")

## ---- cd4scatter
cd4scatter <- ggplot(subset(dataf, cd4within30), aes(x=infPeriod, y=firstcd4cnt)) + 
  geom_point() + facet_grid(everHad~mode2) + 
  theme_bw()
cd4scatter2 <- cd4scatter+ geom_smooth() 
cd4scatter2 <- cd4scatter+ geom_rug(col=rgb(.5,0,0,alpha=.2)) +
  scale_x_continuous(name='Infection window length') +
  scale_y_continuous(name='First CD4 Count')
# I added a rug plot to show densities. See here for marginal histograms. The
# tricky part would be the panels, I think.
# http://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2

cd4dens <- ggplot(subset(dataf, !is.na(cd4_days) & cd4_days<=30 &
                           !is.na(firstcd4cnt) & everHad!='Missing')) + 
  geom_density(aes(x=firstcd4cnt, group=everHad, 
                   colour=everHad, fill=everHad), fill=NA) +
  facet_grid(.~mode2) +
  scale_x_continuous(name='First CD4 count',limits=c(0,2000)) + 
  labs(colour='') +
  theme_bw() +
  scale_colour_brewer(palette="Set1")

# Focus on never-testers and shading below or above CD4=500
tempdf <- subset(dataf, !is.na(cd4_days) & cd4_days<=30 &
                   !is.na(firstcd4cnt) & everHad!='Missing')
gg <- ddply(tempdf, .(mode2, everHad), 
            .fun=function(x) {
            xdens=density(x$firstcd4cnt)$x
            ydens=density(x$firstcd4cnt)$y
            return(data.frame(xdens,ydens))
        })

cd4densUnder500 <- cd4dens + geom_ribbon(data=subset(gg,everHad=="Never Tested" & xdens<500),
            aes(x=xdens,ymax=ydens),ymin=0,alpha=0.3) 
cd4densOver500 <- cd4dens + geom_ribbon(data=subset(gg,everHad=="Never Tested" & xdens>500),
                                         aes(x=xdens,ymax=ydens),ymin=0,alpha=0.3)



#############################################################
# Small investigation of missing at randomness - includes BED

## ---- missAtRand
missAtRand <- ddply(dataf, .(everHad2, cd4within302), summarise, 
                    percBEDpos = 100*sum(BED2fac=='BED +')/length(BED2fac),
                    meanCD4 = mean(firstcd4cnt, na.rm=TRUE),
                    medCD4 = median(firstcd4cnt, na.rm=TRUE),
                    percConcurrent = 100*sum(aidsAtDx2=='Concurrent Dx')/
                      length(aidsAtDx2)
)

# Tidy up
missAtRand <- subset(missAtRand, cd4within302==TRUE)
rownames(missAtRand) <- gsub('TH', '', missAtRand$everHad2)
missLong <- t(missAtRand[,3:6])
rownames(missLong) <- c('Percent BED +', 'Mean CD4', 'Median CD4', 
                        'Percent with Concurrent Dx')
missLong <- missLong[,2:1]

#############################################################
# CD4 versus infPeriod crosstab

## ---- cd4vinf
cd4vinf <- with(subset(dataf, everHad!='Missing'), 
                table(infBreaks, cd4cat))

#http://stackoverflow.com/questions/20673584/visualizing-crosstab-tables-with-a-plot-in-r
# Create the data frame
df <- melt(as.data.frame(cd4vinf))
df <- transform(df, prop=round(100*value/sum(value)))
df$label <- paste0(sprintf("%.0f",
                           df$prop), 
                   "%")
df$dcolor <- 'No Impact'

infl <- levels(df$infBreaks)
cd4l <- levels(df$cd4cat)

df$dcolor[df$infBreaks==infl[4] & df$cd4cat%in%cd4l[2:4]] <- 'Impact'
df$dcolor[df$infBreaks==infl[3] & df$cd4cat%in%cd4l[3:4]] <- 'Impact'
df$dcolor[df$infBreaks==infl[2] & df$cd4cat%in%cd4l[4]] <- 'Impact'

df <- transform(df, cd4cat=factor(cd4cat, 
                                  levels=c(cd4l[c(5,1:4)]),
                                  labels=c(cd4l[c(5,1:4)])))

#Plot the Data
cd4vinfPlot <- ggplot(df, aes(infBreaks, cd4cat)) + 
  geom_point(aes(size = value, colour=dcolor)) + 
  theme_bw() + xlab("Window Length") + ylab("CD4 Bin") +
  scale_colour_brewer(palette="Set2")
cd4vinfPlot <- cd4vinfPlot + scale_size_continuous(range=c(10,30)) + geom_text(aes(label = label)) +
  guides(size=FALSE) + labs(color='')

# Now by MSM

df2 <- ddply(subset(dataf,everHad!='Missing'), 
             .(mode2, infBreaks, cd4cat), summarise, 
             value=length(cd4cat))

#http://stackoverflow.com/questions/20673584/visualizing-crosstab-tables-with-a-plot-in-r
# Create the data frame
df2$prop=NA
df2 <- within(df2, {
  prop[mode2=='MSM'] <- round(100*value[mode2=='MSM']/sum(value[mode2=='MSM']))
  prop[mode2=='non-MSM'] <- round(100*value[mode2=='non-MSM']/sum(value[mode2=='non-MSM']))
})

df2$label <- paste0(sprintf("%.0f",
                            df2$prop), 
                    "%")
df2$dcolor <- 'No Impact'

#infl <- levels(df2$infBreaks)
#cd4l <- levels(df2$cd4cat)

df2$dcolor[df2$infBreaks==infl[4] & df2$cd4cat%in%cd4l[2:4]] <- 'Impact'
df2$dcolor[df2$infBreaks==infl[3] & df2$cd4cat%in%cd4l[3:4]] <- 'Impact'
df2$dcolor[df2$infBreaks==infl[2] & df2$cd4cat%in%cd4l[4]] <- 'Impact'

df2 <- transform(df2, cd4cat=factor(cd4cat, 
                                    levels=levels(df$cd4cat),
                                    labels=levels(df$cd4cat))
)

#Plot the Data
cd4vinfPlotMSM <- ggplot(df2, aes(infBreaks, cd4cat)) + geom_point(aes(size = value, colour=dcolor)) + theme_bw() + xlab("Window Length") + ylab("CD4 Bin")
cd4vinfPlotMSM <- cd4vinfPlotMSM + scale_size_continuous(range=c(10,30)) + geom_text(aes(label = label)) +
  guides(size=FALSE) + labs(color='') + facet_grid(.~mode2) +
  scale_colour_brewer(palette="Set2")


#############################################################
# Proof of equivalence between Base Case and Base Case Alt (Continuous)

## ---- bcAltTID
bcVbcalt <- estimateTID(dataf$infPeriod,
                       intLength=0.25, 
                       cases=c('base_case','base_case_alt'))
Sx_tab <- summary(bcVbcalt, times=c(0,0.25,0.5, 1,5,18), intLength=0.25)[,c(1,3,5)]
colnames(Sx_tab) <- c('Time', 'Original Base Case', 'Alternate Base Case')
print(xtable(Sx_tab,
             caption='Base Case TIDs using different computational approaches',
             label='tab:Sx_bcAlt',
             digits=3),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)


## ---- bcAltTIDPlot
plot(bcVbcalt, intLength=0.25)


#############################################################
# Proof of equivalence between Fake CD4 Case and Base Case Continuous

## ---- cd4caseTID
cd4fake <- estimateTID(dataf$infPeriod,
                       intLength=0.25, 
                       cases=c('base_case_alt','cd4_case'),
                       medWindows=dataf$infPeriod/2,
                       infPeriodOrig=dataf$infPeriod)
Sx_tab <- summary(cd4fake, times=c(0,0.25,0.5, 1,5,18), 
                  intLength=0.25)[,c(1,3,5)]
colnames(Sx_tab) <- c('Time', 'Alternative Base Case', 'Fake CD4 Case')
print(xtable(Sx_tab,
             caption='Base Case versus Fake CD4 Case TIDs',
             label='tab:Sx_cd4fake',
             digits=3),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)


## ---- cd4caseTIDPlot
plot(cd4fake, intLength=0.25)


#############################################################
# Setting up real CD4-based medians

## ---- cd4meds
# Define our literature-based median times to infection by CD4 bin
(cd4meds <- data.frame(cd4lower=c(500,350,200),
                      cd4upper=c(2000, 500, 350),
                      medWindow=c(1.5, 4, 8))
)

#*********************
# Define who should get a CD4-based median
cd4breaks <- c(0,200,350,500,2000)
windowbreaks <- c(0,3,8,16,18)

dataf <- within(dataf, {
            # Non-missing testing history
            hasTestHist <- !is.na(everHadNegTest)
            # CD4 measured within 30d
            cd4within30 <- hasTestHist & !is.na(cd4_days) & cd4_days<=30 &
                !is.na(firstcd4cnt)
            # Categories
            cd4cat=cut(firstcd4cnt, breaks=cd4breaks,
                   include.lowest=TRUE, right=FALSE)
            # Windows
            infBreaks=cut(infPeriod, breaks=windowbreaks,
                          include.lowest=TRUE, right=FALSE)
            # Windows with 'yrs'
            infBreaksYrs=ifelse(is.na(infBreaks), infBreaks,
                                paste(infBreaks,'yrs'))
          })
with(dataf, table(hasTestHist))
with(dataf, table(cd4within30))

#*********************
# Assign medians

# Start with 1/2 of infPeriod, which is just the Base Case.
# Update to CD4-based median if indicated by infPeriod (infection window)
# Define our literature-based median times to infection by CD4 bin
cd4meds <- data.frame(cd4lower=c(500,350,200),
                      cd4upper=c(2000, 500, 350),
                      medWindow=c(1.5, 4, 8))

#*********************
# Assign medians

# Start with 1/2 of infPeriod, which is just the Base Case.
# Update to CD4-based median if indicated by infPeriod (infection window)
dataf <- transform(dataf, medWindows=infPeriod/2, impacted=0)

for (i in 1:nrow(cd4meds)) {
  dataf <- transform(dataf, temp=cd4within30 &
                                   firstcd4cnt>=cd4meds[i,'cd4lower'] &
                                   firstcd4cnt<cd4meds[i, 'cd4upper'] &
                                   infPeriod>=2*cd4meds[i, 'medWindow'])
    dataf <- transform(dataf, impacted=ifelse(temp==1,1,impacted))
    dataf <- within(dataf, {
                        medWindows[hasTestHist & cd4within30 &
                                   firstcd4cnt>=cd4meds[i,'cd4lower'] &
                                   firstcd4cnt<cd4meds[i, 'cd4upper'] &
                                   infPeriod>=2*cd4meds[i, 'medWindow']] <- 
                                   cd4meds[i,'medWindow']
              })
}

# Was expecting 296 cases impacted; need to find the 6
with(dataf, sum(medWindows!=infPeriod/2, na.rm=TRUE))
# Ok
with(dataf, table(mode2, impacted))
with(dataf, table(mode2, impacted)/rowSums(table(mode2,impacted)))
# Now look among the 3016 with testing history
with(subset(dataf,!is.na(everHadNegTest)), table(mode2, impacted))
with(subset(dataf,!is.na(everHadNegTest)), table(mode2, impacted)/rowSums(table(mode2,impacted)))

# Show old and new median windows AMONG the 3016 contributing to testing histories
ddply(subset(dataf,!is.na(everHadNegTest)), .(mode2,cd4cat), summarise, 
      N_impacted=sum(impacted),
      avgOldMedian=round(mean(infPeriod/2, na.rm=TRUE),1),
      avgNewMedian=round(mean(medWindows, na.rm=TRUE),1),
      Difference=avgOldMedian-avgNewMedian)
                    

#############################################################
# Median differences plot (Continuous)

## ---- medDiff

tabMedDiff <- ddply(subset(dataf,everHad!='Missing'),
                    .(cd4cat), summarise, 
                    N_CD4bin=length(impacted),
                    N_impacted=sum(impacted, na.rm=TRUE),
                    Perc_CD4bin=100*sum(impacted)/length(impacted),
                    avgOldMedian=round(mean(infPeriod/2, na.rm=TRUE),1),
                    avgNewMedian=round(mean(medWindows, na.rm=TRUE),1),
                    Difference=avgOldMedian-avgNewMedian,
                    PercChange=round(100*(avgOldMedian-avgNewMedian)/avgOldMedian,1))
tabMedDiff <- transform(tabMedDiff, Prop_CD4bin=N_CD4bin/sum(N_CD4bin))
colnames(tabMedDiff) <- c('CD4 Bin', 
                          'CD4 Bin Size',
                          'Number Impacted', '% of CD4 Bin Impacted',
                          'Base Case Median', 
                          'CD4 Case Median', 'Difference', 'Percent Change',
                          'Proportion in CD4 Bin')
tabMedDiff <- tabMedDiff[1:4,]

t = melt(tabMedDiff, id.vars=c('CD4 Bin', 'CD4 Bin Size'))
t1 = subset(t, variable%in%c('Base Case Median', 'CD4 Case Median'))
          #  &`CD4 Bin`!='No CD4 within 30d')

#Plot the Data
medDiffPlot <- ggplot(t1, aes(`CD4 Bin`, y=value, fill=variable)) + 
  geom_bar(position='dodge', stat='identity') + theme_bw() + labs(fill='') +
  scale_fill_manual(values=c('blue', 'orange2'), name='') +
  scale_y_continuous(name='Median time since infection (years)')

# Not stratified by CD4, but by MSM
medwindDF <- data.frame(rbind(transform(subset(dataf,everHad!='Missing'),
                                   medWindows=infPeriod/2,
                                   case='Base Case Median'
                                   ),
                             transform(subset(dataf,everHad!='Missing'),
                                         case='CD4 Case Median'
                                         ))
                         )
medDiffPlotMarg <- ggplot(medwindDF, 
                      aes(mode2, y=medWindows, fill=case)) + 
  geom_bar(position='dodge', stat='summary', fun.y='mean') + theme_bw() + labs(fill='') +
  scale_fill_manual(values=c('blue', 'orange2'), name='') +
  scale_y_continuous(name='Median time since \ninfection (years)') +
  scale_x_discrete(name='') +
  theme(legend.position='bottom')

medwindStats <- ddply(medwindDF, .(mode2, case), 
                      summarise, AvgMed=mean(medWindows))
medwindStats <- dcast(medwindStats, mode2~case)
medwindStats <- transform(medwindStats, 
                          `Absolute Difference` = `CD4 Case Median`-`Base Case Median`,
                          `Percent Change` = -100*(1-(`CD4 Case Median`/`Base Case Median`)),
                          check.names=FALSE)
colnames(medwindStats)[1] <- ''

# By CD4 and MSM
tabMedDiffMSM <- ddply(subset(dataf,everHad!='Missing'), 
                       .(mode2, cd4cat), summarise, 
                       N_CD4bin=length(impacted),
                       N_impacted=sum(impacted, na.rm=TRUE),
                       Perc_CD4bin=100*sum(impacted)/length(impacted),
                       avgOldMedian=round(mean(infPeriod/2, na.rm=TRUE),1),
                       avgNewMedian=round(mean(medWindows, na.rm=TRUE),1),
                       Difference=avgOldMedian-avgNewMedian,
                       PercChange=round(100*(avgOldMedian-avgNewMedian)/avgOldMedian,1))
colnames(tabMedDiffMSM) <- c('Mode', 'CD4 Bin', 
                             'CD4 Bin Size',
                             'Number Impacted', '% of CD4 Bin Impacted',
                             'Base Case Median', 
                             'CD4 Case Median', 'Difference', 
                             'Percent Change'
                             )
tabMedDiffMSM <- tabMedDiffMSM[tabMedDiffMSM[,'CD4 Bin']!='No CD4 within 30d',]

t = melt(tabMedDiffMSM)
t1 = subset(t, variable%in%c('Base Case Median', 'CD4 Case Median'))
       #     &`CD4 Bin`!='No CD4 within 30d')

#Plot the Data
medDiffPlotMSM <- ggplot(t1, aes(`CD4 Bin`, y=value, fill=variable)) + 
  geom_bar(position='dodge', stat='identity') + theme_bw() +
  facet_grid(.~Mode) + labs(fill='') +
  theme(legend.position='bottom')

## ---- medDiff2

medDiffs <- transform(subset(dataf, !(is.na(everHadNegTest))),
                      origMed=infPeriod/2,
                      medDiff=(infPeriod/2)-medWindows)

medDiff2tab <- ddply(subset(medDiffs,medDiff>0), .(mode2), summarise,
                     p25=round(quantile(medDiff, probs=0.25),2),
                     mean=round(mean(medDiff),2),
                     med=round(median(medDiff),2),
                     p75=round(quantile(medDiff, probs=0.75),2))
medDiff2tab <- transform(medDiff2tab, 
                         lab=paste(mode2, med, sep='\n'),
                         xmed=med+c(0.3,-0.3),
                         xmean=mean+c(0.3,-0.3)
                         )
print(medDiff2tab, digits=3)

medDiff2Plot <- ggplot(subset(medDiffs,medDiff>0), aes(medDiff, group=mode2)) + 
  geom_histogram(aes(fill=mode2), position='dodge') +
  theme_bw() + labs(fill='') + 
  scale_x_continuous(name='Difference (years)') +
  scale_y_continuous(name='') +
  geom_vline(data=medDiff2tab,
             aes(xintercept=med, colour=mode2)) +
  geom_text(data=medDiff2tab, aes(x=xmed, y=10, label=med,colour=mode2), 
            size=5, vjust = -1.5,hjust=0.5) +
  guides(colour=FALSE) 

medDiff3Plot <- ggplot(subset(medDiffs,medDiff>0), aes(medDiff, group=mode2)) + 
  geom_density(aes(x=medDiff, group=mode2, 
                   colour=mode2, fill=mode2), fill=NA) +
  scale_x_continuous(name='Difference (years)') + 
  labs(colour='') +
  theme_bw() +
  geom_vline(data=medDiff2tab,
             aes(xintercept=med, colour=mode2)) +
  geom_text(data=medDiff2tab, aes(x=xmed, y=0.2, label=med,colour=mode2), 
            size=5, vjust = -0.5,hjust=0.5) +
  guides(fill=FALSE)

medDiff4Plot <- ggplot(subset(medDiffs,medDiff>0), aes(medDiff, group=mode2)) + 
  geom_density(aes(x=medDiff, group=mode2, 
                   colour=mode2, fill=mode2), fill=NA) +
  scale_x_continuous(name='Difference (years)') + 
  labs(colour='') +
  theme_bw() +
  geom_vline(data=medDiff2tab,
             aes(xintercept=mean, colour=mode2)) +
  geom_text(data=medDiff2tab, aes(x=xmean, y=0.2, label=mean,colour=mode2), 
            size=5, vjust = -0.5,hjust=0.5) +
  guides(fill=FALSE) 


#############################################################
# Estimate TIDs


## ---- cd4caseTIDReal
cd4real <- estimateTID(dataf$infPeriod,
                       intLength=0.25, 
                       cases=c('base_case_alt','cd4_case'),
                       medWindows=dataf$medWindows,
                       infPeriodOrig=dataf$infPeriod)
cd4real.MSM <- estimateTID(subset(dataf, mode2=='MSM')$infPeriod,
                       intLength=0.25, 
                       cases=c('base_case_alt','cd4_case'),
                       medWindows=subset(dataf, mode2=='MSM')$medWindows,
                       infPeriodOrig=subset(dataf, mode2=='MSM')$infPeriod)
cd4real.nonMSM <- estimateTID(subset(dataf, mode2!='MSM')$infPeriod,
                       intLength=0.25, 
                       cases=c('base_case_alt','cd4_case'),
                       medWindows=subset(dataf, mode2!='MSM')$medWindows,
                       infPeriodOrig=subset(dataf, mode2!='MSM')$infPeriod)

Sx_tab <- summary(cd4real, times=c(0,0.25,0.5, 1,5,18), 
                  intLength=0.25)[,c(1,3,5)]
Sx_tab2 <- summary(cd4real.MSM, times=c(0,0.25,0.5, 1,5,18), 
                  intLength=0.25)[,c(1,3,5)]
Sx_tab3 <- summary(cd4real.nonMSM, times=c(0,0.25,0.5, 1,5,18), 
                  intLength=0.25)[,c(1,3,5)]
Sx_tab <- data.frame(Pop=c('All', rep('', nrow(Sx_tab)-1),
                           'MSM', rep('', nrow(Sx_tab)-1),
                           'non-MSM', rep('', nrow(Sx_tab)-1)),
                     rbind(Sx_tab, Sx_tab2, Sx_tab3))
colnames(Sx_tab) <- c('Population', 'Time', 
                      'Alternative Base Case', 'CD4 Case')
print(xtable(Sx_tab,
             caption='Base Case versus CD4 Case TIDs',
             label='tab:cd4real_tab',
             digits=3),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)


## ---- cd4caseTIDRealPlotMSM
plot(cd4real.MSM, 0.25)

## ---- cd4caseTIDRealPlotnonMSM
plot(cd4real.nonMSM, 0.25)


## ---- cd4caseTIDpdfs
pdfMSM <- summary(cd4real.MSM, times=c(0,0.25,0.5, 1,5,18), 
                  intLength=0.25)[,c(1,2,4)]
pdfnonMSM <- summary(cd4real.nonMSM, times=c(0,0.25,0.5, 1,5,18), 
                  intLength=0.25)[,c(1,2,4)]
pdfs <- data.frame(Pop=c('MSM', rep('', nrow(pdfMSM)-1),
                           'non-MSM', rep('', nrow(pdfMSM)-1)),
                     rbind(pdfMSM, pdfnonMSM))
colnames(pdfs)[3:4] <- c('bc_pdf', 'cd4case_pdf')
pdfs <- transform(pdfs, diff=cd4case_pdf-bc_pdf, ratio=bc_pdf/cd4case_pdf)

print(xtable(pdfs,
             caption='Base Case versus CD4 Case PDFs',
             label='tab:cd4real_pdf',
             digits=3),
      caption.placement='top',
      table.placement='ht',
      size='small',
      include.rownames=FALSE)

#############################################################
# Estimate TIDs - focusing on the PDFs

#############################################################
# Investigate probability reassigned and impact on  TIDs

## ---- percProbReassigned

# Compute BC probability assigned within the median window:
# just 1/infPeriod times the medWindow. Then compare that to 
# 0.5, which is how much the CD4 Case assigs within the median window

dataf <- within(dataf, {
                  medProbBC=(medWindows)*(1/infPeriod)
                  probReassigned=0.5-medProbBC
      })

summary(dataf$probReassigned)

ddply(subset(dataf,!is.na(everHadNegTest)), .(mode2), summarise, 
      totalReassigned=sum(probReassigned, na.rm=TRUE),
      propReassigned=sum(probReassigned)/length(probReassigned))

# Look separately among impacted cases
ddply(subset(dataf,!is.na(everHadNegTest) & impacted==1), .(mode2), summarise, 
      totalReassigned=sum(probReassigned),
      propReassigned=sum(probReassigned)/length(probReassigned))


## ---- meanTID
timeStep <- 0.25
yearTimes <- seq(0,18,by=timeStep)

# Get full TID curves
msmSx <- summary(cd4real.MSM, times=yearTimes, intLength=0.25)
msmSx$mode='MSM'
nonmsmSx <- summary(cd4real.nonMSM, times=yearTimes, intLength=0.25)
nonmsmSx$mode='non-MSM'
FullSx <- rbind(msmSx, nonmsmSx)

# Multiply to get discrete AUC
means <- ddply(FullSx, .(mode), summarise,
      bc_auc=0.25*sum(`base_case_alt S(x)`),
      cd4_auc=0.25*sum(`cd4_case S(x)`))
means <- transform(means, ratio=cd4_auc/bc_auc, diff=bc_auc-cd4_auc)
print(means, digits=2)

## ---- medianTID

## ---- start if (runEstimation)
if (runEstimation) {

#############################################################
# Prepare for estimation

## ---- true_prevalence

# Read in true prevalence
  trueprev_data = read.csv(file.path(workd,'data/Reported_prevalence_2010-2014.csv'),
                           na.string="",
                           stringsAsFactor=FALSE, 
                           check.names=FALSE)


#############################################################
# Estimate undiagnosed cases and incidence

## ---- run_subgroups

  these_cases <- c('base_case_alt', 'cd4_case')
  names(these_cases) <- c('Base Case', 'CD4 Case')
  subgroups <- runSubgroups(dataf,
                            subvar='mode2',
                            intLength=0.25,
                            cases=these_cases,
                            medWindowsVar='medWindows',
                            prev=trueprev_data,
                            save=file.path(workd, 'analysis_WA/results/2016_trueprev_CD4Case.csv'))

#############################################################
# Process undiagnosed results

## ---- process_subgroups

# Function to extract desired comparative results
  compareUndx <- function(x, subgroups, name='') {
      totRes <- subgroups[[x]]$results

      # Summary of summaries
      sumtable <- subset(totRes$resultsSummary[order(totRes$resultsSummary$Estimate),],
                         Estimate=='Undiagnosed Cases')
      sumtable2 <- subset(totRes$resultsSummaryYear[order(totRes$resultsSummaryYear$Estimate),],
                          Estimate=='Undiagnosed Cases')
      sumtable$Year <- '2005-2014'
      sumtable <- rbind(sumtable, sumtable2)[,c('Year', colnames(sumtable)[-ncol(sumtable)])]
      colnames(sumtable)[which(colnames(sumtable)=='Diagnoses/Case')] <- 'Case'
      m <- subset(melt(sumtable), variable=='Mean', select=-Estimate)
      mwide <- dcast(m, Year~Case, value.var='value')
      mwide$Difference <- mwide[,2]-mwide[,3]
      mwide <- transform(mwide, `Percent Change`=round(100*Difference/`Base Case`),
                         check.names=FALSE)
      return(data.frame(Group=name,mwide,check.names=FALSE))
  }

# Combined comparison of undiagnosed estimates
compareAll <- rbind(compareUndx('Total-stratified', subgroups, 'Total'),
                    compareUndx('MSM', subgroups, 'MSM'),
                    compareUndx('non-MSM', subgroups, 'non-MSM'))

# Compare undiagnosed fractions
compareFrac  <- function(x, subgroups, name='') {
    totTP <- subgroups[[x]]$trueprev
    totTP <- rename(totTP, c('Diagnoses/Case'='Case'))
    totTP <- subset(totTP, Estimate==unique(Estimate)[2] | 
                    Estimate==unique(Estimate)[4],
                select=c('Year', 'Case', 'Estimate', 'Mean'))
    mwide <- dcast(totTP, Year+Estimate~Case, value.var='Mean')
    mwide$Difference <- mwide[,3]-mwide[,4]
    mwide <- transform(mwide, `Percent Change`=round(100*Difference/`Base Case`),
                         check.names=FALSE)
    return(data.frame(Group=name,mwide,check.names=FALSE))
}

compareAll <- rbind(compareFrac('Total-stratified', subgroups, 'Total'),
                    compareFrac('MSM', subgroups, 'MSM'),
                    compareFrac('non-MSM', subgroups, 'non-MSM'))

#############################################################
# Process incidence

## ---- incidence
resForInc <- rbind(subgroups$MSM$results$resultsAll,
                    subgroups[['non-MSM']]$results$resultsAll)


## ---- end if (runEstimation)
}



