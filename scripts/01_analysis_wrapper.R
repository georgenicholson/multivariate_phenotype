# rm(list = ls())

path_to_dir <- "C:/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/github_multivariate_phenotype/multivariate_phenotype"
renv::activate(path_to_dir)
renv::restore(path_to_dir)
arguments <- commandArgs()

if("--args" %in% arguments){
  argnam.in <- data.frame(nam = c("run_type", "scen", "subsamseed"), 
                          coersion.fn = c("as.character", "as.numeric", "as.numeric"), 
                          stringsAsFactors = F)
  for(i in 1:nrow(argnam.in))
    assign(argnam.in$nam[i], eval(call(argnam.in$coersion.fn[i], arguments[grep("--args", arguments) + i])))
} else {
  run_type <- c("demo", "main", "benchmark", "test_benchmark")[1]
  scen <- 1
  subsamseed <- 1
}

##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control to contain parameters
control <- get_control_parameters_mv()

##########################################
# Create directory structure
for(dirc in c(control$output_dir, control$methods_comp_dir, control$global_res_dir, control$data_dir, control$train_test_samples_dir))
  dir.create(dirc, showWarnings = F)

##########################################
# Load data
Data_all <- readRDS(file = control$Data_all_file)

##########################################
# Get table of analyses
analysis_table <- create_table_of_analyses(control = control, check_status = T, run_type = run_type)

Data <- analysis_table$Data[scen]
Meth <- analysis_table$Meth[scen]
N <- analysis_table$N[scen]
P <- analysis_table$P[scen]
nSig <- analysis_table$nSig[scen]
Data <- analysis_table$Data[scen]
XDmeth <- Meth

Y_raw <- Data_all[[Data]]$Y_raw
Y_zeroed <- Data_all[[Data]]$Y_zeroed
S_raw <- Data_all[[Data]]$S_raw
S_zeroed <- Data_all[[Data]]$S_zeroed
Y.eml <- S.eml <- list()
dat.typev <- switch(Data, 
                    impc = c("raw", "zero"), 
                    eqtl = "raw")
Y.eml$raw <- Data_all[[Data]]$Y_raw
Y.eml$zero <- Data_all[[Data]]$Y_zeroed
S.eml$raw <- Data_all[[Data]]$S_raw
S.eml$zero <- Data_all[[Data]]$S_zeroed

###################################
# Choose phen subset if 
ph_all <- colnames(Y_raw)
if (run_type %in% c("demo", "test_benchmark")) {
  phens_ok_to_use <- ph_all[colMeans(is.na(Y_raw[, ph_all])) < .8]
  set.seed(1)
  phens_to_use <- sample(x = phens_ok_to_use, size = P)
} else {
  phens_to_use <- ph_all
}

###################################
# Create and Load train-test samples information 
train_test_list <- get_train_test_split(control = control, N = N, P = P, Data = Data)
train_test_list_file_curr <- file.path("output", "train_test_splits", paste0(Data, "_N_", N, "_P_", P, ".RDS"))
if (!file.exists(train_test_list_file_curr)) {
  saveRDS(object = train_test_list, file = train_test_list_file_curr)
}
    
###################################
# Assign samples
for (var_assign_subsample in c("sams_for_lik_cross_val", 
                               "sams_for_model_testing", 
                               "sams_for_model_training", 
                               "sams_for_cor_est",
                               "sams_for_strong_cov_est")) {
  subsam_mat_curr <- train_test_list[[var_assign_subsample]]
  assign(x = var_assign_subsample, value = rownames(subsam_mat_curr)[subsam_mat_curr[, subsamseed]])
}

##############################################
# Estimate R
if(Data == "impc"){
  R.init <- cor((Y_zeroed / S_zeroed)[sams_for_cor_est, phens_to_use], use = "p", meth = "p")
  R.init[is.na(R.init)] <- 0
  diag(R.init) <- 1
  if(qr(R.init)$rank < P)
    R.init <- (1 - control$rank_deficient_R_eps) * R.init + control$rank_deficient_R_eps * diag(rep(1, P))
}
if(Data == "eqtl"){
  data.in <- list(Bhat = Y_zeroed[sams_for_cor_est, phens_to_use], Shat = S_zeroed[sams_for_cor_est, phens_to_use])
  R.init <- mashr::estimate_null_correlation_simple(data = data.in, z_thresh = 2)
}
dimnames(R.init) <- list(phens_to_use, phens_to_use)


##############################################
# Run MVphen
if (grepl("MVphen", Meth)) {
  source(file = "scripts/01a_fit_MVphen.R")
}


##################################################
# Fit models using XD and/or MASH
if(Meth %in% c("mash", "XD")){
  source(file = "scripts/01b_fit_XD_and_mash.R")
}

