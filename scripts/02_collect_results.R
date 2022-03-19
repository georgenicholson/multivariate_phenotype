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

##########################################
# Get table of analyses
# run_type <- c("demo", "main", "benchmark", "test_benchmark")[3]
analysis_table <- create_table_of_analyses(control = control, check_status = T, run_type = "benchmark")

resl.comp <- compl <- Sigll <- Ksigl <- pimatl <- Sighatl <- Rl <- omegaseql <- samll <- objl <- list()
resl.comp$uv <- list(mn = Data_all$impc$Y_raw[Data_all$impc$sam_names, Data_all$impc$meas_names], 
                     sd = Data_all$impc$S_raw[Data_all$impc$sam_names, Data_all$impc$meas_names])

for (scen in 1:nrow(analysis_table)) {
  for (var_name in c("N", "P", "Data", "Meth", "nSig", "MVphen_K", "n_subsamples", "file_core_name")) {
    assign(var_name, analysis_table[scen, var_name])
  }
  sam_names <- Data_all[[Data]]$sam_names
  N_all <- Data_all[[Data]]$N_all
  # train_test_list_file_curr <- file.path("output", "train_test_splits", paste0(Data, "_N_", N, "_P_", P, ".RDS"))
  # train_test_list <- readRDS(file = train_test_list_file_curr)
  train_test_list <- get_train_test_split(control = control, Data_all = Data_all, N = N, P = P, Data = Data, n_subsamples = n_subsamples)
  
  meas_names <- train_test_list$phens_to_use
  P_all <- length(meas_names)
  
  if (is.na(MVphen_K)) {
    MVphen_K <- control$nfac
  }
  facnam <- paste0("fac_", 1:MVphen_K)
  array1 <- array(NA, dim = c(N_all, P_all, n_subsamples), dimnames = list(sam_names, meas_names, 1:n_subsamples))
  array2 <- array(NA, dim = c(N_all, n_subsamples), dimnames = list(sam_names, 1:n_subsamples))
  array.fac <- array(NA, dim = c(N_all, MVphen_K, n_subsamples), dimnames = list(sam_names, facnam, 1:n_subsamples))
  crossval.llv <- rep(NA, n_subsamples)
  shared.formatl <- list(mnarr = array1, sdarr = array1, loocv.mnarr = array1, loocv.sdarr = array1, 
                         lfsrarr = array1, llmat = array2, llmat.raw = array2, llmat.zero = array2)
  if(grepl("MVphen", Meth)) {
    namc <- paste(Data, Meth, "nSig", nSig, "K", MVphen_K, sep = "_")
  } else {
    namc <- paste(Data, Meth, "nSig", nSig, sep = "_")
  }
  compl[[namc]] <- shared.formatl
  objl[[namc]] <- list()
  pimatl[[namc]] <- Ksigl[[namc]] <- Sigll[[namc]] <- Sighatl[[namc]] <- Rl[[namc]] <- omegaseql[[namc]] <- samll[[namc]] <- list()
  for (subsamseed in 1:n_subsamples) {
    objl[[namc]][[subsamseed]] <- list()
    file_list <- get_file_list(control = control, file_core_name = analysis_table[scen, "file_core_name"], subsamseed = subsamseed)
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
    
    samll[[namc]][[subsamseed]] <- objl[[namc]][[subsamseed]]$saml <- 
      list(sams.for.cor.est = sams_for_cor_est, 
           sams.for.model.fitting = sams_for_model_training, 
           sams.for.testing = sams_for_model_testing, 
           sams.for.lik.cross.val = sams_for_lik_cross_val)
    res.store <- NULL
    if(grepl("MVphen", Meth)) {
      emloaded <- resloaded <- loocv.loaded <- fac.loaded <- FALSE
      if (file.exists(file_list$emout.file.namc)) {
        emout.mix <- readRDS(file = file_list$emout.file.namc)
        emloaded <- TRUE
      }
      if (file.exists(file_list$res.store.namc)) {
        resl.store <- readRDS(file = file_list$res.store.namc)
        resloaded <- TRUE
        res.store <- resl.store$raw
      }
      if (file.exists(file_list$loocv.res.store.namc)) {
        loocv.store <- readRDS(file = file_list$loocv.res.store.namc)
        loocv.loaded <- TRUE
      }
      if (file.exists(file_list$fac.res.store.namc)) {
        loocv.store <- readRDS(file = file_list$fac.res.store.namc)
        fac.loaded <- TRUE
      }
      if (!emloaded | !resloaded) {
        warning(paste0(namc, " subsamseed ", subsamseed, " not loaded"))
      }
      if (namc == control$mv_meth_nam_use) {#Meth == "MVphen" & nSig == 1 & Data == "impc"){
        if (loocv.loaded) {
          compl[[namc]]$loocv.mnarr[match(rownames(loocv.store$mn), sam_names), 
                                    match(colnames(loocv.store$mn), meas_names), subsamseed] <- loocv.store$mn
          compl[[namc]]$loocv.sdarr[match(rownames(loocv.store$sd), sam_names), 
                                    match(colnames(loocv.store$sd), meas_names), subsamseed] <- loocv.store$sd
        }
      }
    }
    
    if(Meth %in% c("mash", "XD")) {
      res.store <- NULL
      file_to_load <- switch(Meth,
                             mash = file_list$mash.resl.file.namc,
                             XD = file_list$XD.resl.file.namc)
      mash_or_XD_loaded <- FALSE
      if (file.exists(file_to_load)) {
        resl.store <- readRDS(file = file_to_load)
        datatypec <- switch(Data, eqtl = "raw", impc = "zero")
        if (datatypec  %in% names(resl.store)) {
          res.store <- resl.store[[datatypec]]
          mash_or_XD_loaded <- TRUE
        } 
      }
      if (!mash_or_XD_loaded) {
        warning(paste0(file_to_load, " not loaded"))
      }
    }
    if (!is.null(res.store)) {
      Sigll[[namc]][[subsamseed]] <- objl[[namc]][[subsamseed]]$Sigl <- 
        lapply(res.store$Sigl, function(M) M[match(meas_names, rownames(M)), match(meas_names, colnames(M))])
      Rl[[namc]][[subsamseed]] <- objl[[namc]][[subsamseed]]$R <- res.store$R[match(meas_names, rownames(res.store$R)), match(meas_names, colnames(res.store$R))]
      pimatl[[namc]][[subsamseed]] <- objl[[namc]][[subsamseed]]$pimat <- res.store$pimat
      omegaseql[[namc]][[subsamseed]] <- objl[[namc]][[subsamseed]]$omegaseq <- res.store$omegaseq
      Sighatl[[namc]][[subsamseed]] <- objl[[namc]][[subsamseed]]$Sighat <- 0
      for(sc in 1:length(Sigll[[namc]][[subsamseed]])){
        Sighatl[[namc]][[subsamseed]] <- Sighatl[[namc]][[subsamseed]] + sum(pimatl[[namc]][[subsamseed]][, sc] * omegaseql[[namc]][[subsamseed]]) *  Sigll[[namc]][[subsamseed]][[sc]]
        objl[[namc]][[subsamseed]]$Sighat <- objl[[namc]][[subsamseed]]$Sighat + sum(pimatl[[namc]][[subsamseed]][, sc] * omegaseql[[namc]][[subsamseed]]) *  Sigll[[namc]][[subsamseed]][[sc]]
      }
      compl[[namc]]$mnarr[match(rownames(res.store$mn), sam_names), match(colnames(res.store$mn), meas_names), subsamseed] <- res.store$mn
      compl[[namc]]$sdarr[match(rownames(res.store$sd), sam_names), match(colnames(res.store$sd), meas_names), subsamseed] <- res.store$sd
      compl[[namc]]$lfsrarr[match(rownames(res.store$lfsr), sam_names), match(colnames(res.store$lfsr), meas_names), subsamseed] <- res.store$lfsr
      compl[[namc]]$llmat[match(names(res.store$loglikv), sam_names), subsamseed] <- res.store$loglikv
      compl[[namc]]$llmat.raw[match(names(resl.store$raw$loglikv), sam_names), subsamseed] <- resl.store$raw$loglikv
      compl[[namc]]$llmat.zero[match(names(resl.store$zero$loglikv), sam_names), subsamseed] <- resl.store$zero$loglikv
    }
  }
  
  ########################################################
  # Calculate combined Sig and R
  suppressWarnings(lmat.norm <- exp(compl[[namc]]$llmat - apply(compl[[namc]]$llmat, 1, function(v) max(v, na.rm = T))))
  pmix <- compl[[namc]]$pmix <- lmat.norm / rowSums(lmat.norm, na.rm = T)
  Sigmix.p <- colSums(pmix, na.rm = T) / sum(pmix, na.rm = T)
  compl[[namc]]$Sig.comb <- compl[[namc]]$R.comb <- 0
  for(subsamseed in 1:n_subsamples){
    if (Sigmix.p[subsamseed] > 0) {
      compl[[namc]]$Sig.comb <- compl[[namc]]$Sig.comb + Sighatl[[namc]][[subsamseed]][meas_names, meas_names] * Sigmix.p[subsamseed]
      compl[[namc]]$R.comb <- compl[[namc]]$R.comb + Rl[[namc]][[subsamseed]][meas_names, meas_names] * Sigmix.p[subsamseed]
    }
  }
  resl.comp[[namc]]$Sig.comb <- compl[[namc]]$Sig.comb
  resl.comp[[namc]]$Sigcor.comb <- t(resl.comp[[namc]]$Sig.comb / sqrt(diag(resl.comp[[namc]]$Sig.comb))) / sqrt(diag(resl.comp[[namc]]$Sig.comb))
  eigc <- eigen(resl.comp[[namc]]$Sigcor.comb)
  resl.comp[[namc]]$facs.varimax <- varimax(eigc$vectors[, 1:MVphen_K])$loadings
  resl.comp[[namc]]$facs.promax <- promax(eigc$vectors[, 1:MVphen_K])$loadings
  dimnames(resl.comp[[namc]]$facs.varimax) <- dimnames(resl.comp[[namc]]$facs.promax) <- list(meas_names, facnam)
  wt.mnarr <- apply(sweep(compl[[namc]]$mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  wt.varr <- apply(sweep(compl[[namc]]$sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  wt.mnsqarr <- apply(sweep(compl[[namc]]$mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  wt.lfsrarr <- apply(sweep(compl[[namc]]$lfsrarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  compl[[namc]]$mn.comb <- wt.mnarr
  compl[[namc]]$sd.comb <- sqrt(wt.varr + wt.mnsqarr - wt.mnarr^2)
  compl[[namc]]$lfsr.comb <- wt.lfsrarr
  if(Meth == "MVphen" & nSig == 1 & Data == "impc"){
    loocv.wt.mnarr <- apply(sweep(compl[[namc]]$loocv.mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
    loocv.wt.varr <- apply(sweep(compl[[namc]]$loocv.sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
    loocv.wt.mnsqarr <- apply(sweep(compl[[namc]]$loocv.mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
    compl[[namc]]$loocv.mn.comb <- loocv.wt.mnarr
    compl[[namc]]$loocv.sd.comb <- sqrt(loocv.wt.varr + loocv.wt.mnsqarr - loocv.wt.mnarr^2)
  }
  resl.comp[[namc]] <- c(resl.comp[[namc]], list(mn = compl[[namc]]$mn.comb[sam_names, meas_names], 
                                                 sd = compl[[namc]]$sd.comb[sam_names, meas_names], 
                                                 loocv.mn = compl[[namc]]$loocv.mn.comb[sam_names, meas_names], 
                                                 loocv.sd = compl[[namc]]$loocv.sd.comb[sam_names, meas_names],
                                                 lfsr = compl[[namc]]$lfsr.comb[sam_names, meas_names]))
}


saveRDS(resl.comp, file = control$file.resl.comp)
saveRDS(compl, file = control$file.compl)
saveRDS(objl, file = control$file.objl)


