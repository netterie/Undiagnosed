
########################
# true_prevalence.R
# 
# Purpose: combines results from run_main.R with reported
# prevalences to produce true prevalence estimates
# 
# Dependencies:
# results in the style of run_main.csv
#
# History: none
#
########################

            
######################################################################
# RUN ALL
######################################################################

true_prevalence(thiscategory='Total',
                thisgroup='Total',
                thispop='WA',
                results_file=file.path('analysis_WA', 
                                       'results',
                                       'run_main.csv'),
                prevalence_file=file.path('data',
                                          'Reported_Prevalence.csv'))

true_prevalence(thiscategory='Mode-consolidated',
                thisgroup='MSM',
                thispop='WA',
                results_file=file.path('analysis_WA', 
                                       'results',
                                       'MSM.csv'),
                prevalence_file=file.path('data',
                                          'Reported_Prevalence.csv'))

true_prevalence(thiscategory='Mode-consolidated',
                thisgroup='MSM',
                thispop='WA',
                results_file=file.path('analysis_WA', 
                                       'results',
                                       'MSM.csv'),
                prevalence_file=file.path('data',
                                          'Reported_Prevalence.csv'),
                use_imputed_results='Yes',
                otherID='-imputed')

true_prevalence(thiscategory='Mode-consolidated',
                thisgroup='non-MSM',
                thispop='WA',
                results_file=file.path('analysis_WA', 
                                       'results',
                                       'non-MSM.csv'),
                prevalence_file=file.path('data',
                                          'Reported_Prevalence.csv'))

true_prevalence(thiscategory='Mode-consolidated',
                thisgroup='non-MSM',
                thispop='WA',
                results_file=file.path('analysis_WA', 
                                       'results',
                                       'non-MSM.csv'),
                prevalence_file=file.path('data',
                                          'Reported_Prevalence.csv'),
                use_imputed_results='Yes',
                otherID='-imputed')

true_prevalence(thiscategory='Mode-consolidated',
                thisgroup=c('MSM', 'non-MSM'),
                thispop='WA',
                results_file=c(file.path('analysis_WA', 
                                       'results',
                                       'non-MSM.csv'),
                               file.path('analysis_WA',
                                         'results',
                                         'MSM.csv')),
                prevalence_file=file.path('data',
                                          'Reported_Prevalence.csv'),
                aggregate='Total-weighted',
                returntp=TRUE)

######################################################################
# GRAPH ALL
######################################################################


trueprev_fig_WA <- trueprev_figure(
                        files=c('true_prevalence_Mode-consolidated_MSM.csv',
                                'true_prevalence_Mode-consolidated_non-MSM.csv'),
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2014-07-14_WA',
                        group='Group')

trueprev_fig_MSM <- trueprev_figure(
                        files='true_prevalence_from_compare_MSM_report.csv',
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2014-07-14_MSM',
                        group='Population')

trueprev_fig_MSM_Undiag <- trueprev_figure(
                        files='true_prevalence_from_compare_MSM_report.csv',
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2014-07-14_MSM_Undiag',
                        this_subset='Estimated Undiagnosed',
                        group='Population')

trueprev_fig_MSM_PLWHA <- trueprev_figure(
                        files='true_prevalence_from_compare_MSM_report.csv',
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2014-07-14_MSM_PLWHA',
                        this_subset='Estimated PLWHA',
                        group='Population')

trueprev_fig_MSM_PercUndiag <- trueprev_figure(
                        files='true_prevalence_from_compare_MSM_report.csv',
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2014-07-14_MSM_PercUndiag',
                        this_subset='Estimated % Undiagnosed',
                        group='Population')

trueprev_fig_All <- trueprev_figure(
                        files='true_prevalence_Total_Total.csv',
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2014-07-14_All',
                        group=NULL)
             
trueprev_fig_AllandWeighted <- trueprev_figure(
                        files=c('true_prevalence_Total_Total.csv',
                                'true_prevalence_Total-weighted_MSM+non-MSM.csv'),
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='total_v_msm+non',
                        group='Category')
trueprev_fig_4cases_MSMvNon <- trueprev_figure_4cases(
                        files=c('true_prevalence_Mode-consolidated_MSM.csv',
                                'true_prevalence_Mode-consolidated_non-MSM.csv',
                        'true_prevalence_Mode-consolidated_MSM-imputed.csv',
                                'true_prevalence_Mode-consolidated_non-MSM-imputed.csv'),
                        fileIDs=c('Used Missing', 'No', 'No', 'Yes', 'Yes'),
                        dir_files='/Users/jeanette/Dropbox/School/PhD/HIV_WA/analysis_WA/results',
                        suffix='2013-allcases-MSMvnon-MSM',
                        years=c(2013),
                        group='Group')
             


