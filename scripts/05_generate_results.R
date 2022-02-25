rm(list = ls())
try({
  path_to_dir <- "C:/Users/nicho/Documents/GitHub/multivariate_phenotype"
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
restabl <- readRDS(file = file.path(control$global_res_dir, paste0("restabl_comb.RDS")))



# source("scripts/05a_heatmap_of_hits.R")
# source("scripts/05b_power_comparisons.R")
# source("scripts/05c_reference_lines.R")
source("scripts/05d_gene_ontology.R")
