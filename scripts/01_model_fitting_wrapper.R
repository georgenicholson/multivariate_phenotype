
arguments <- commandArgs()
if("--args" %in% arguments){
  argnam.in <- data.frame(nam = c("run_type", "scen", "subsamseed"), 
                          coersion.fn = c("as.character", "as.numeric", "as.numeric"), 
                          stringsAsFactors = F)
  for(i in 1:nrow(argnam.in))
    assign(argnam.in$nam[i], eval(call(argnam.in$coersion.fn[i], arguments[grep("--args", arguments) + i])))
} else {
  run_type <- c("demo", "main", "benchmark", "test_benchmark")[2]
  scen <- 1
  if(run_type == "main") scen <- 1
  subsamseed <- 1
}
print(Sys.info())

library(foreach)
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
for (dirc in c(control$output_dir, control$methods_comp_dir, control$global_res_dir, control$data_dir, control$train_test_samples_dir)) {
  dir.create(dirc, recursive = TRUE)
}
print(getwd())
print(paste("data dir exits:", dir.exists("data")))

##########################################
# Download data
if (!file.exists(control$Data_all_file)) {
  source("scripts/00_download_data.R")
}

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
MVphen_K <- analysis_table$MVphen_K[scen]
n_subsamples <- analysis_table$n_subsamples[scen]
XDmeth <- Meth

cl <- parallel::makeCluster(25)
doParallel::registerDoParallel(cl = cl)
# parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
foreach(subsamseed = 1:analysis_table$n_subsamples, .errorhandling = "pass", .verbose = TRUE) %dopar% {
  
  file_list <- get_file_list(control = control, 
                             file_core_name = analysis_table[scen, "file_core_name"], 
                             subsamseed = subsamseed)
  
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
  # Create and Load train-test samples information 
  force_overwrite_train_test_splits <- FALSE
  if (!file.exists(train_test_list_file_curr) | force_overwrite_train_test_splits) {
    train_test_list <- get_train_test_split(control = control, Data_all = Data_all, N = N, P = P, Data = Data, n_subsamples = n_subsamples)
    train_test_list_file_curr <- file.path("output", "train_test_splits", paste0(Data, "_N_", N, "_P_", P, ".RDS"))
    saveRDS(object = train_test_list, file = train_test_list_file_curr)
  } else {
    train_test_list <- readRDS(file = train_test_list_file_curr)
  }
  phens_to_use <- train_test_list$phens_to_use
  
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
  # Initilize Sig list
  Sigl.em.init <- initialize_Sig_list(control = control, 
                                      Y = Y_zeroed[sams_for_model_training, phens_to_use], 
                                      nSig = nSig, 
                                      random_init = grepl("rand", Meth))
  
  
  ##############################################
  # Run MVphen
  if (grepl("MVphen", Meth)) {
    force_rerun_em <- FALSE
    if(force_rerun_em | !file.exists(file_list$emout.file.namc)){
      print("Running EM algorithm")
      emout.mix <- EM_algo_mixture_multi_Sig(control = control, 
                                             Y.em = Y_raw[sams_for_model_training, phens_to_use], 
                                             S.em = S_raw[sams_for_model_training, phens_to_use], 
                                             MVphen_K = MVphen_K,
                                             Sigl.em.init = Sigl.em.init, 
                                             R.em.init = R.init)
      print("Calculating posterior means")
      resl.store <- list()
      for(dat.type in dat.typev){
        out.post.mix <- em.update.function(Y.em = Y.eml[[dat.type]][sams_for_model_testing, phens_to_use], 
                                           S.em = S.eml[[dat.type]][sams_for_model_testing, phens_to_use],
                                           Sigl = emout.mix$Sigl, 
                                           R = emout.mix$R, 
                                           omegaseq = emout.mix$omegaseq,
                                           pimat = emout.mix$pi, 
                                           meth = "post.mn")
        resl.store[[dat.type]] <- list(mn = out.post.mix$mnmat[sams_for_model_testing, phens_to_use], 
                                       sd = out.post.mix$sdmat[sams_for_model_testing, phens_to_use],
                                       loglikv = out.post.mix$loglikv[sams_for_model_testing],
                                       lfsr = out.post.mix$lfsrmat[sams_for_model_testing, phens_to_use], 
                                       Sigl = lapply(emout.mix$Sigl, function(M) M[phens_to_use, phens_to_use]), 
                                       Sig.mn = emout.mix$Sigmn[phens_to_use, phens_to_use],
                                       Ksig = emout.mix$Ksig, R = emout.mix$R[phens_to_use, phens_to_use], 
                                       pimat = emout.mix$pi, omegaseq = emout.mix$omegaseq)
      }
      saveRDS(object = emout.mix, file = file_list$emout.file.namc)
      saveRDS(object = resl.store, file = file_list$res.store.namc)
    } else {
      print("Loading previously run EM output")
      resl.store <- readRDS(file = file_list$res.store.namc)
      emout.mix <- readRDS(file = file_list$emout.file.namc)
    }
    
    if(run_type == "main" & Data == "impc"){
      print("Calculating leave-one-procedure-out predictions")
      emout.mix <- readRDS(file = file_list$emout.file.namc)
      procun <- unique(Data_all$impc$phmap$procnam)
      matout <- matrix(NA, length(sams_for_model_testing), length(phens_to_use), dimnames = list(sams_for_model_testing, phens_to_use))
      loocv.store <- list(mn = matout, sd = matout)
      for(procc in procun){
        ph.leave.out <- Data_all$impc$phmap$ph[Data_all$impc$phmap$procnam == procc]
        Y_raw.left.out <- Data_all$impc$Y_raw
        Y_raw.left.out[, ph.leave.out] <- NA
        S_raw.left.out <- Data_all$impc$S_raw
        S_raw.left.out[, ph.leave.out] <- NA
        post.mn.loocv <- em.update.function(Y.em = Y_raw.left.out[sams_for_model_testing, phens_to_use], 
                                            S.em = S_raw.left.out[sams_for_model_testing, phens_to_use],
                                            Sigl = emout.mix$Sigl, 
                                            R = emout.mix$R, 
                                            omegaseq = emout.mix$omegaseq,
                                            pimat = emout.mix$pi, 
                                            meth = "post.mn")
        loocv.store$mn[sams_for_model_testing, ph.leave.out] <- post.mn.loocv$mnmat[sams_for_model_testing, ph.leave.out]
        loocv.store$sd[sams_for_model_testing, ph.leave.out] <- post.mn.loocv$sdmat[sams_for_model_testing, ph.leave.out]
      }
      saveRDS(loocv.store, file = file_list$loocv.res.store.namc)
    }
    # source(file = "scripts/01a_fit_MVphen.R")
  }
  
  ##################################################
  # Fit models using XD and/or MASH
  if(Meth %in% c("mash", "XD")){
    source(file = "scripts/01b_fit_XD_and_mash.R")
  }
}  
parallel::stopCluster(cl)
