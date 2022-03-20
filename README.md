# multivariate_phenotype

This repository contains the R scripts needed to reproduce the results reported 
in the manuscript 'Illuminating the mammalian genome with multivariate phenotype analysis'. 

## Installation

To run these scripts, you will need R version 3.6.3 or later, widely available on 
Unix-like, Windows and Mac families of operating systems. If you have R version 4.0.0
or above on Windows, you will also need to install 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). 

To get started, first [clone](https://git-scm.com/book/en/v2/Git-Basics-Getting-a-Git-Repository)
this repository onto your local machine in folder "<your_path>/multivariate_phenotype"

Next, open an R console, and install the 
[renv](https://rstudio.github.io/renv/index.html) R package if you don't have it 
already (e.g. via `install.packages("renv")`). Then, run the following, 
changing `path_to_dir` to the path of your local version of this repository to install the required packages for the scripts
```
path_to_dir <- "<your_path>/multivariate_phenotype"
setwd(path_to_dir)
renv::activate(path_to_dir)
renv::restore(path_to_dir)
```
 

## Data

Accompanying data provided in "data/Data_all.RDS" which can be loaded in R via
```
Data_all <- readRDS("data/Data_all.RDS")
```

## Demo run

You can run a demonstration of the software on a smaller data set via
```
source("scripts/01_model_fitting_wrapper.R")
```

## Full analysis


