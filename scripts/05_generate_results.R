rm(list = ls())
try({
  path_to_dir <- "C:/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/github_multivariate_phenotype/multivariate_phenotype"
  setwd(path_to_dir)
  renv::activate(path_to_dir)
  renv::restore(path_to_dir)
})
getwd()

Data <- "impc"

##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control contains parameter settings
control <- get_control_parameters_mv()

###################################
# Load data
Data_all <- readRDS(control$Data_all_file)
linemap <- Data_all$impc$linemap
phmap <- Data_all$impc$phmap

resl.comp <- readRDS(file = control$file.resl.comp)
compl <- readRDS(file = control$file.compl)
objl <- readRDS(file = control$file.objl)
resll <- readRDS(file = control$file.resll)

source("scripts/05a_power_comparisons.R")
