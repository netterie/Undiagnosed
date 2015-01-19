
## Estimating Undiagnosed HIV using Testing Histories

*Updated 1/19/15*

### Coding To-Do's

####More urgent
* Split files in HIVBackCalc/R by function, for helping developers. Use nice section headings within those .R files 
* Embed non-proprietary datasets to use with tutorials, e.g. simulated KC data
* Create vignette 
* Release minimal package to accompany PLoS One paper
* Update Shiny to use this package

####Less urgent
* Map out code terminology - in the "Index to functions" below
* Improve terminology
    + "infPeriod" variable -> "tid" 
    + "fig1" function -> "tidSurvivorFxn" ?
    + "impute" -> includeMiss vs excludeMiss ?
    + Learn about search and replace across multiple files
    + Test code and make sure it works with updated terminology
* Turn scripts (e.g., format_data.R) into functions that can be embedded in the package and essentially serve as a wrapper for running a sequence of functions?
* Get Sam's input on 
    + Function-ifying scripts
    + When to use Markdown and HTMl vs Sweave and PDF formats
    + Pros/cons of eventually migrating to using S2, S3 or S4 classes
* Add two modes to setup_hivbackcalc: package updated = TRUE (just load HIVBackCalc package) vs FALSE (need to source the package files and source the sim'd data) ?


### Overview of folders

Folder | Description | Notes
------ | ----------- | -----
HIVBackCalc | Folder for the R package | Has not yet been compiled to include any code in "other.R" 
analysis_ian | Ian's original files |   
analysis_KC | Jeanette's replication of KC analysis |   
analysis_WA | Jeanette's extension of method to WA state data | Also contains some comparisons of the original KC analysis to a KC analysis only using the WA dataset
literature | Some relevant papers and data |   
presentations | Presentations by our group |   

### Code structure, workflow and replication for the analysis_KC folder

Structure is fairly obvious: the data formatting file is sourced by the run file, so the run file is self-contained.

Source "run_JKB.R" or knit "replication_JKB.Rmd". You should get the results reported in the submitted paper. Both these files are basically Ian's code from run.R, but with section headers and comments to clarify the code.

### Code structure and workflow for the analysis_WA folder

####Guide to .R code

R file | Description | Notes 
------ | ----------- | -----
format_data.R | Formats the WA state data for use with the method | Look in .pdf report for explanations/details
describe_data.R | Does some EDA on the formatted WA state data | This was very preliminary and probably could be greatly improved 
run_main.R | Runs the analysis for a single group, whether it is a subset of the data or the full data | If subsetting, assumes that object "subset_before_run" has been set to TRUE beforehand and the subgroup name is stored in the object "u"
run_main_subgroups.R | Loops through subgroups and runs the analysis for each one

_Note on knitr hooks and code chunks_ 

The .R files in this folder contain knitr "hooks" that define code "chunks". These "chunks" are read in at the beginning of reports and then executed at subsequent points in the report. 

The advantage of this approach is that it makes it easy to to have R code that you can execute _by section_ throughout a report, as opposed to executing it all at once by sourcing the whole .R file. But it's really convenient that the code is also fully contained and maintainable through a standalone .R file that _can_ be sourced all together, when desired. 

The process is:

1. Insert "hooks" into R code to define "chunks" of code. Hook syntax is
```{r}
## ---- hookname ----
```
2. At the beginning of your .Rnw or .Rmd report file, after you load the knitr library, use the read_chunk() function to read in the R code, e.g. 
```{r}
read_chunk('path_to_file.R')
```
3. Throughout the report, execute code chunks using the hook name, e.g.
```{r}
<<hookname, echo=FALSE>>=
@
```

####Guide to analyses stored in this folder

Analysis | Files are named... | Notes
-------- | ------------------ | -----
WA state analysis of undiagnosed counts | HIVBackCalc_full_report | At the end, this report runs the analysis for a bunch of subgroups. That part takes a while.
WA state analysis of undiagnosed fraction | true_prevalence_report | Requires that undiagnosed counts have already been obtained--just reads in the results file plus a denominator file
Comparison of KC MSM undiagnosed counts using data provided by KC versus data provided by WA state | compare_MSM_report | Relies on code within the report to set up and perform the comparison analysis

####Guide to file types

Extension | Purpose
--------- | -------
.R | R code for that (part of) analysis. Expect knitr "hooks" that identify "chunks" of code
.Rnw | R-noweb file, an alternative to .Rmd, that uses LaTeX instead of Markdown for the writing of text | After loading knitr library, run and compile using "knit2pdf('filename.Rnw')"
.tex | Created by knit2pdf() - this is a standalone TeX file that can be compiled into a pdf
.pdf | PDF report


### Replicating analyses in analysis_WA 

#### Suggested procedure

1. Download HIVBackCalc_1.01.tar.gz (Unix) or HIVBackCalc_1.01.zip (Windows) from this repository, and install the package using the following instructions:
```{r}
# Open R in/set the R workding directory to the folder containing the downloaded zip file and type:
# (switch in .tar.gz for .zip if on Unix)
install.packages('HIVBackCalc_1.01.zip', repos=NULL)
library(HIVBackCalc)
```
2. Identify the .Rnw file for the appropriate analysis
3. In the setup_hivbackcalc() function call(s) near the top of the report, change the working directory to reflect yours.
4. Load the knitr library in R
5. Use the purl() function from the knitr library to run only the R code without generating the report. I think that if you use the tangle=TRUE option, you will get an .R file of all the code in the report. You could then execute code from that .R file.

### Index to functions

####Estimation Functions

Name | Functionality | Demo files | Customized aspects | Comments
---- | ------------- | ---------- | ------------------ | --------------
estimateProbDist | Returns a PDF of the TID, using the Base Case assumption | | None known | Output is referred to as "pid" in later code
empirProbDist | Returns discrete time PDF of the TID, using the Upper Bound assumption | | None known | Output is referred to as "pid" in later code. I think this is a custom-coding of a standard empirical PDF, but I haven't checked it. This functon is not in the 1.0 version of the package
meanEMupdate | EM update step | | None known | 
estimateIncidence | Backcalculates incidence | | None known | 
estimateUndiagnosed | From estimated incidence and TID, estimates undiagnosed | | None known | Interpret results in light of the time step by which diagnoses were entered (e.g., per quarter-year)
print.backproj | Prints incidence results to screen | | None known | 

####Formatting Functions

####Plotting Functions

Name | Functionality | Demo files | Customized aspects | Comments
---- | ------------- | ---------- | ------------------ | --------------
plot.backproj | Plots backcalculated incidence with diagnoses overlayed | | None known | Is this used? I think all plots are now ggplot2
