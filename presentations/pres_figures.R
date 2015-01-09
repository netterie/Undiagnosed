########################
# pres_figures.R
# 
# Purpose: 
# Generate graphs specifically needed for presentations
# See headers for which presentations
#
# Dependencies:
# format_data.R
# packages ggplot2, reshape2, MASS, survival, gridExtra
#
# History: 
# Based off of analysis_WA/describe_data.R
#
# Notes:
########################

#############################################################
# SETUP
#############################################################
source('/Users/jeanette/Dropbox/School/PhD/HIV_WA/HIVBackCalc/R/other.R')

# Load libraries, data and data-cleaning file
# Eventually this function should return the cleaned dataset,
# but data-cleaning has the name hardcoded as msm and I'm not
# going to generalize that right now
setup_hivbackcalc(workd='/Users/jeanette/Dropbox/School/PhD/HIV_WA',
                  datafile='data/MSM_KingCounty_rev.csv',
                  source_these='analysis_KC/data-cleaning_JKB.R',
                  msm=TRUE)

# Also load the survival package
library(survival)

#############################################################
# PRES: UNDIAGNOSED FRACTION UPDATE 2014-Oct
#############################################################

# OVERLAYED HISTOGRAMS/DENSITIES

histdat <- rbind(data.frame(History='Yes',
                            infPeriod=subset(msm, 
                                             everTested==TRUE)$infPeriod),
                 data.frame(History='Yes+No',
                            infPeriod=subset(msm, 
                                             !is.na(infPeriod))$infPeriod),
                 data.frame(History='Yes+No+Miss',
                            infPeriod=msm$infPeriod_imputeNA))
histdat <- transform(histdat, 
                     History=factor(History, 
                                    levels=c('Yes+No+Miss', 'Yes+No', 'Yes'), 
                                    labels=c('Yes+No+Miss', 'Yes+No', 'Yes')))

hist1 <- ggplot(subset(histdat, History!='Yes+No+Miss'),
                aes(infPeriod, color=History)) +
         geom_density(alpha=0.8, size=2) +
         scale_color_manual(values=c('#66c2a5', '#fc8d62')) +
         theme_bw()

hist1 <- ggplot(subset(histdat, History!='Yes+No+Miss'),
                aes(infPeriod, fill=History, color=History)) +
         geom_density(alpha=0.2, size=1.5) +
         theme_bw() +
         scale_color_manual(values=c('#66c2a5', '#fc8d62')) +
         scale_fill_manual(values=c('#66c2a5', '#fc8d62')) +

hist1 <- ggplot(subset(histdat, History!='Yes+No+Miss'),
                aes(infPeriod, fill=History)) +
         geom_histogram(position='identity') +
         theme_bw() +
         scale_fill_manual(values=c('#66c2a5', '#fc8d62'))

hist2 <- ggplot(histdat,
                aes(infPeriod, color=History, fill=History)) +
         geom_density(alpha=0.2, size=1.5) +
         scale_color_manual(values=c('#66c2a5', '#fc8d62', '#8da0cb')) +
         scale_fill_manual(values=c('#66c2a5', '#fc8d62', '#8da0cb')) +
         theme_bw()
        

# SEQUENTIAL OVERLAYED TIDS

fig1basedat <- fig1(msm$infPeriod, returnFig=FALSE)
fig1base <- fig1(msm$infPeriod, d1=subset(fig1basedat, var=='Base Case  '), 
                 returnFig=TRUE) + ggtitle('Base Case') + aes(size=1.5)
fig1upper <- fig1(msm$infPeriod) + ggtitle('Worst Case (Obs)') + aes(size=1.5)

# Color scheme edit
ggplotColours <- function(n=6, h=c(0, 360) +15){
      if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
fig1worstmiss <- fig1combined(msm) + aes(size=1.5) +
                 scale_color_manual(values=c(ggplotColours(n=2), '#8da0cb'))

# WA SUBGROUP RESULTS 

sg <- read.csv('/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results/run_main_subgroups.csv')

sg <- sg[grepl('Undiagnosed', sg$Quantity),]
sg <- within(sg,{
             Case=as.character(Quantity)
             Case[grep('Base Case', Case)] <- 'Base Case'
             Case[grep('Worst Case (Obs)', Case)] <- 'Worst Case (Obs)'
             Case[grep('Worst Case (Miss)', Case)] <- 'Worst Case (Miss)'
                 })

sg1 <- subset(sg, Subgroup=='MSM' | Subgroup=='Hetero')
sg2 <- subset(sg, Subgroup=='White' | Subgroup=='Black' | Subgroup=='Hisp')

sgf <- rbind(data.frame(sg1, Group='MSM v Hetero'),
             data.frame(sg2, Group='Race'))

pl <- ggplot(sgf,
             aes(x=Subgroup, y=Median,
                 colour=Case,
                 shape=Case,
                 size=2))+
      scale_x_discrete(name='') +
      scale_y_continuous(name='') + 
      geom_point() + 
      theme_bw() + 
      facet_grid(.~Group)

pl1 <- ggplot(sg1,
             aes(x=Subgroup, y=Median,
                 colour=Case,
                 shape=Case,
                 size=2))+
      scale_x_discrete(name='') +
      scale_y_continuous(name='') + 
      geom_point() + 
      theme_bw()
pl2 <- ggplot(sg2,
             aes(x=Subgroup, y=Median,
                 colour=Case,
                 shape=Case,
                 size=2))+
      scale_x_discrete(name='') +
      scale_y_continuous(name='') + 
      geom_point() + 
      theme_bw()
