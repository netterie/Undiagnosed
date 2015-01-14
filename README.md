
## Estimating Undiagnosed HIV using Testing Histories

*Updated 1/13/15*

### Coding To-Do's

* Embed non-proprietary datasets to use with tutorials, e.g. simulated KC data
* Consider renaming some variables and functions, e.g. the "infPeriod" variable would be clearer as "tid" and the "fig1" function would be clearer as "tidSurvivorFxn" or something. But, this would require a lot of work to make all the previous code compatible.

### Overview of folders

Folder | Description | Notes
------ | ----------- | -----
HIVBackCalc | Folder for the R package | Has not yet been compiled to include any code in "other.R" 
analysis_ian | Ian's original files |   
analysis_KC | Jeanette's replication of KC analysis |   
analysis_WA | Jeanette's extension of method to WA state data | Also contains some comparisons of the original KC analysis to a KC analysis only using the WA dataset
literature | Some relevant papers and data |   
presentations | Presentations by our group |   

### Code structure and workflow

_These are notes to expand upon and organize_

* Two modes: package updated = TRUE (just load HIVBackCalc package) vs FALSE (need to source other.R)
* .R files have knitr hooks
* .Rnw files execute code from .R files within a report

### Replicating analyses

_Steps to replicate each analysis that has been done so far_

* KC MSM
* WA state - all and subgroups
* KC MSM using KC database vs KC MSM using WA state database

### Index to functions

_Estimation Functions_

Name | Functionality | Demo files | Customized aspects | Comments
---- | ------------- | ---------- | ------------------ | --------------
estimateProbDist | Returns a PDF of the TID, using the Base Case assumption | | None known | Output is referred to as "pid" in later code
empirProbDist | Returns discrete time PDF of the TID, using the Upper Bound assumption | | None known | Output is referred to as "pid" in later code. I think this is a custom-coding of a standard empirical PDF, but I haven't checked it. This functon is not in the 1.0 version of the package
meanEMupdate | EM update step | | None known | 
estimateIncidence | Backcalculates incidence | | None known | 
estimateUndiagnosed | From estimated incidence and TID, estimates undiagnosed | | None known | Interpret results in light of the time step by which diagnoses were entered (e.g., per quarter-year)
print.backproj | Prints incidence results to screen | | None known | 

_Formatting Functions_

_Plotting Functions_

Name | Functionality | Demo files | Customized aspects | Comments
---- | ------------- | ---------- | ------------------ | --------------
plot.backproj | Plots backcalculated incidence with diagnoses overlayed | | None known | Is this used? I think all plots are now ggplot2