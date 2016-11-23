#############################################################
# run_undx.R
#
# use read.chunk to run in a report
# expects the following to be defined:
# - undx_repo (path to Undiagnosed repository)
# - prev_path (path to prevalence file within undx_repo)
# - results_path (path to results within undx_repo)
# - these_cases (a named vector of allowed cases)
# - diagInterval (0.25, 0.5 or 1)
# - run_estimation (logical indicating whether to run the back-calc)
# Eventually this should be fleshed out and turned into
# a run_analysis function, I think
#############################################################


#############################################################
# Prepare for estimation

## ---- true_prevalence

# Read in true prevalence
trueprev_data = read.csv(file.path(undx_repo, prev_path),
                         na.string="",
                         stringsAsFactor=FALSE, 
                         check.names=FALSE)


#############################################################
# Estimate undiagnosed cases and incidence

## ---- run_subgroups
if (run_estimation) {
    subgroups <- runSubgroups(dataf,
                              subvar='mode2',
                              intLength=diagInterval,
                              cases=estimation_cases,
                              prev=trueprev_data,
                              save=file.path(undx_repo,results_path))
}
