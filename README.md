# multivariate_phenotype

This repository contains the R scripts needed to reproduce the results reported 
in the manuscript 'Illuminating the mammalian genome with multivariate phenotype analysis'. 

## Installation

To run these scripts, you will need R version 3.6.3 or later, widely available on 
Unix-like, Windows and Mac families of operating systems. If you have R version 4.0.0
or above on Windows, you will also need to install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). 
# The demo below has been tested on macOS 10.15 Catalina, Windows 10 and CentOS Linux 7. 
To get started, first [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
this repository onto your local machine. Next, open an R console, and install the 
[renv](https://rstudio.github.io/renv/index.html) R package if you don't have it 
already (e.g. via `install.packages("renv")`). Then, run the following, 
changing `path_to_dir` to the path of your local version of this repository,
```
path_to_dir <- "path/to/multivariate_phenotype"
setwd(path_to_dir)
renv::activate(path_to_dir)
renv::restore(path_to_dir)
```
to install the required packages for the scripts. 

## Demo

The repository includes a subset of the data for the lower-tier local authority 
of Adur and its corresponding PHE region of the South East. To load the data and
the `prevdebiasr` package:

```
base::load("data/example.RData")
library(prevdebiasr)
```

This loads two data frames: `southeast_df` containing weekly Pillar 1+2 testing 
data and REACT study data for the whole of the South East, and `adur_df` 
containing weekly Pillar 1+2 data for Adur, from the end of May 2020 to the
beginning of August 2021. The columns `nt` and `Nt` contain the number of
of positive and total Pillar 1+2 tests respectively, while the columns `nr` and 
`Nr` contain the corresponding number of REACT tests. The column `M` contains 
the population of the relevant region.  
  
To generate estimates for &delta; for the South East, i.e. the log odds ratio 
of being tested in the infected versus the non-infected subpopulations: 

```
control_debias <- get_control_parameters(alpha_testing = 3e-4)
delta_df <- specify_delta_prior(southeast_df, control_debias)
# 2 minutes, 28 seconds
```

The `get_control_parameters()` function is a utility to set the hyperparameters 
of the statistical model. See the manuscript and `?get_control_parameters` for 
more details. The output `delta_df` is a data frame with two columns: 
`delta_prior_mean` and `delta_prior_sd` contain the mean and standard deviation 
of a moment-matched Gaussian approximation for &delta;, to be used downstream 
in estimating debiased prevalence.

```
head(delta_df)
#   delta_prior_mean delta_prior_sd
# 1         2.888957      0.3226672
# 2         2.881859      0.3045518
# 3         2.875330      0.2829014
# 4         2.869367      0.2568006
# 5         2.895197      0.2599046
# 6         2.878135      0.2641643
```


To obtain estimates of debiased prevalence for Adur:
```
adur_df <- cbind(adur_df, delta_df)
adur_prevalence <- local_prevalence(adur_df, control_debias)
# 1 minute, 33 seconds
```
The output `adur_prevalence` is a list containing three matrices: `log_post`, 
`norm_post` and `log_lik` corresponding to the log posterior, the normalised
posterior, and the log likelihood for debiased prevalence. Each matrix is of 
dimension 'number of weeks' x 'number of prevalence bins'. The prevalence bins 
used can be found in `control_debias$bin.d`.

```
adur_prevalence$norm_post[1:3, 30:32]
#             45           48           51
# 1 0.0188419722 0.0231530560 0.0266946627
# 2 0.0002964074 0.0002060873 0.0001500834
# 3 0.0483394905 0.0429635767 0.0402508229
```


## Scripts to reproduce manuscript results

To subdirectory `scripts` contains all the code needed to reproduce the results 
in the manuscript. 

`00_wrapper.R`  fetches the required data.  
 


### Notes

