rm(list = ls())
path_to_dir <- "C:/Users/nicho/Dropbox/GitHub_Dropbox/multivariate_phenotype/scripts"
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
  run_type <- c("demo", "main", "benchmark", "test_benchmark")[2]
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
# Download and organise data



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

dat.typev <- switch(Data, impc = c("raw", "zero"), eqtl = "raw")
Y.eml <- S.eml <- list()
Y_raw <- Data_all[[Data]]$Y_raw
Y_zeroed <- Data_all[[Data]]$Y_zeroed
S_raw <- Data_all[[Data]]$S_raw
S_zeroed <- Data_all[[Data]]$S_zeroed
Y.eml$raw <- Data_all[[Data]]$Y_raw
Y.eml$zero <- Data_all[[Data]]$Y_zeroed
S.eml$raw <- Data_all[[Data]]$S_raw
S.eml$zero <- Data_all[[Data]]$S_zeroed

###################################
# Choose phen subset
ph_all <- colnames(Y_raw)
if (run_type %in% c("demo", "test_benchmark")) {
  phens_ok_to_use <- ph_all[colMeans(is.na(Y_raw[, ph_all])) < .8]
  phens_to_use <- sample(x = phens_ok_to_use, size = P)
} else {
  phens_to_use <- ph_all
}

###################################
# Create and Load train-test samples information 
train_test_list <- get_train_test_split(control = control, N = N, P = P, Data_all = Data_all, Data = Data, n_subsamples = 50)
train_test_list_file_curr <- file.path("output", "train_test_splits", paste0(Data, "_N_", N, "_P_", P, ".RDS"))
if (!file.exists(train_test_list_file_curr)) {
  saveRDS(object = train_test_list, file = train_test_list_file_curr)
}
# train_test_list <- readRDS(file = file.path("output", "train_test_splits", paste0(Data, "_N_", N, "_P_", P))
    
    
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
#Initialise R
if(Data == "impc"){
  R.init <- cor((Y_zeroed / S_zeroed)[sams_for_cor_est, phens_to_use], use = "p", meth = "p")
  R.init[is.na(R.init)] <- 0
  diag(R.init) <- 1
  if(qr(R.init)$rank < P)
    R.init <- (1 - control$rank_deficient_R_eps) * R.init + control$rank_deficient_R_eps * diag(rep(1, P))
}
if(Data == "eqtl"){
  # snps.est.cor <- snpmap.sub$snp[snpmap.sub$task == "random.train.est.cor"]
  data.in <- list(Bhat = Y_zeroed[sams_for_cor_est, phens_to_use], Shat = S_zeroed[sams_for_cor_est, phens_to_use])
  R.init <- mashr::estimate_null_correlation_simple(data = data.in, z_thresh = 2)
  # Sig.init <- cov(Y_zeroed[sams_for_model_training, phens_to_use], use = "p", meth = "p")
}
dimnames(R.init) <- list(phens_to_use, phens_to_use)

##############################################
#Initialise list of cov. mat's Sigl.em.init
if (!grepl("rand", Meth)) {
  if(nSig > 1){
    require(mclust)
    mcl.out <- Mclust(data = Y_zeroed[sams_for_model_training, phens_to_use], G = nSig)
    cluster.membership <- apply(mcl.out$z, 1, which.max)
  } else {
    cluster.membership <- rep(1, length(sams_for_model_training))
  }
  names(cluster.membership) <- sams_for_model_training
  Sigl.em.init <- list()
  for(j in 1:nSig){
    sams.in <- names(cluster.membership)[which(cluster.membership == j)]
    Sig.init <- cov(Y_zeroed[sams.in, phens_to_use, drop = F], use = "p", meth = "p")
    Sig.init[is.na(Sig.init)] <- 0
    diag(Sig.init)[is.na(diag(Sig.init))] <- 1
    if(qr(Sig.init)$rank < P)
      Sig.init <- (1 - control$rank_deficient_R_eps) * Sig.init + control$rank_deficient_R_eps * diag(rep(1, P))
    dim(Sig.init)
    dimnames(Sig.init) <- list(phens_to_use, phens_to_use)
    Sigl.em.init[[j]] <- Sig.init
  }
} else {
  Sigl.em.init <- list()
  for(j in 1:nSig){
    Sig.temp <- solve(rWishart(1, P, diag(rep(1, P)))[, , 1])
    Sigl.em.init[[j]] <- t(Sig.temp / sqrt(diag(Sig.temp))) / sqrt(diag(Sig.temp))
  }
}

##############################################
# Run MVphen
if (grepl("MVphen", Meth)) {
  variables_in_filename_use <- control$variables_in_filename
  file.base <- paste0(control$methods_comp_dir, "/", 
                      paste(paste(variables_in_filename_use, sapply(variables_in_filename_use, get), sep = "_"), collapse = "_"))
  
  emout.file.namc <- paste0(file.base, "_emout.RDS")
  res.store.namc <- paste0(file.base, "_res.RDS")
  resl.store.namc <- paste0(file.base, "_resl.RDS")
  fac.res.store.namc <- paste0(file.base, "_facres.RDS")
  loocv.res.store.namc <- paste0(file.base, "_loocv_res.RDS")
  rerun.em <- T
  if(rerun.em | !file.exists(emout.file.namc)){
    print("Running EM algorithm")
    emout.mix <- EM_algo_mixture_multi_Sig(control = control, 
                                           Y.em = Y_raw[sams_for_model_training, phens_to_use], 
                                           S.em = S_raw[sams_for_model_training, phens_to_use], 
                                           MVphen_K = MVphen_K,
                                           Sigl.em.init = Sigl.em.init, 
                                           R.em.init = R.init)
    saveRDS(object = emout.mix, file = emout.file.namc)
  } else {
    print("Loading previously run EM output")
    emout.mix <- readRDS(file = emout.file.namc)
  }
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
  saveRDS(object = resl.store, file = res.store.namc)
  calc.loocv <- F
  if(run_type == "main"){
    print("Calculating leave-one-procedure-out predictions")
    emout.mix <- readRDS(file = emout.file.namc)
    procun <- unique(phmap$procnam)
    matout <- matrix(NA, length(sams_for_model_testing), length(phens_to_use), dimnames = list(sams_for_model_testing, phens_to_use))
    loocv.store <- list(mn = matout, sd = matout)
    for(procc in procun){
      ph.leave.out <- phmap$ph[phmap$procnam == procc]
      Y_raw.left.out <- Y_raw
      Y_raw.left.out[, ph.leave.out] <- NA
      S_raw.left.out <- S_raw
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
    save(loocv.store, file = loocv.res.store.namc)
  }
}


##################################################
# Fit models using Extreme Deconvolution and/or MASH
if(Meth %in% c("mash", "XD")){
  # si <- include.singletons
  # bo <- run.XD
  XDmeth <- Meth
  file.base.mash <- paste0(control$methods_comp_dir, "/", paste(paste(control$variables_in_filename, sapply(control$variables_in_filename, get), sep = "_"), collapse = "_"))
  file.base.XD <- paste0(control$methods_comp_dir, "/", paste(paste(control$variables_in_filename, sapply(control$variables_in_filename, get), sep = "_"), collapse = "_"))
  mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RDS")
  mash.raw.results.file.namc <- paste0(file.base.mash, "_mash_raw_results.RDS")
  XD.output.file.namc <- paste0(file.base.XD, "_XD_output.RDS")
  XD.resl.file.namc <- paste0(file.base.XD, "_XD_resl.RDS")
  library(mashr)
  normU <- T
  use.pt.mass <- T
  ############################################################
  # Get big effects covariance matrices for MASH
  if(Data == "impc"){
    #Choose strongest effects from UV analysis
    uv.t.mat <- (Y_zeroed / S_zeroed)[sams_for_model_training, ]
    uv.t.mat <- uv.t.mat[order(-apply(abs(uv.t.mat), 1, function(v) max(v, na.rm = T))), ]
    max.abs.t <- apply(abs(uv.t.mat), 1, function(v) max(v, na.rm = T))
    t.th <- 4
    min.sams.for.Ztil <- 100
    lines.strong.use <- names(max.abs.t)[max.abs.t > t.th]
    if(length(lines.strong.use) < min.sams.for.Ztil)
      lines.strong.use <- names(max.abs.t)[1:min.sams.for.Ztil]
    big.eff.use <- lines.strong.use
    Ztil <- scale((Y_zeroed / S_zeroed)[big.eff.use, phens_to_use], scale = F)
    if(any(colMeans((Ztil == 0)) == 1)) #If phen's completely missing among big effects then use all samples
      Ztil <- scale((Y_zeroed / S_zeroed)[sams_for_model_training, phens_to_use], scale = F)
  }
  if(Data == "eqtl"){
    # table(snpmap.sub$task)
    snps.strong.use <- sams_for_strong_cov_est
    # snps.strong.use <- snpmap.sub$snp[snpmap.sub$task == "strong.est.cov"]
    big.eff.use <- snps.strong.use
    Ztil <- scale((Y_zeroed / S_zeroed)[big.eff.use, phens_to_use], scale = F)
  }
  Sigl.init.bigeff.mash <- list()
  Sigl.init.bigeff.mash[[1]] <- cov(Ztil)
  svd.Ztil <- svd(Ztil)
  Sigl.init.bigeff.mash[[2]] <- 0
  K1 <- min(3, P)
  K2 <- min(5, P)
  for(i in 1:K1)
    Sigl.init.bigeff.mash[[2]] <- Sigl.init.bigeff.mash[[2]] + svd.Ztil$d[i] * svd.Ztil$v[, i] %*% t(svd.Ztil$v[, i])
  Sigl.init.bigeff.mash[[3]] <- 0
  for(i in 1:K2)
    Sigl.init.bigeff.mash[[3]] <- Sigl.init.bigeff.mash[[3]] + svd.Ztil$d[i] * svd.Ztil$v[, i] %*% t(svd.Ztil$v[, i])
  names(Sigl.init.bigeff.mash) <- paste0("U.", 1:3)
  
  mashdata.for.model.fitting <- mash_set_data(Bhat = Y_zeroed[sams_for_model_training, phens_to_use], 
                                              Shat = S_zeroed[sams_for_model_training, phens_to_use], V = R.init)
  
  
  
  mashdata.big.effects.for.XD <- mash_set_data(Bhat = Y_zeroed[big.eff.use, phens_to_use], 
                                               Shat = S_zeroed[big.eff.use, phens_to_use], V = R.init)
  
  ##############################################
  # Run Extreme Deconvolution (XD)
  if(Meth == "XD"){
    if(XDmeth == "mash"){
      mashdata.XD <- mashdata.big.effects.for.XD
      Ul.XD.init <- Sigl.init.bigeff.mash
    }
    if(XDmeth == "XD"){
      mashdata.XD <- mashdata.for.model.fitting
      Ul.XD.init <- Sigl.em.init
    }
    print("Running Extreme Deconvolution")
    replace.XD <- T
    if(!file.exists(XD.output.file.namc) | replace.XD){
      XD.out <- mashr:::bovy_wrapper(data = mashdata.XD, Ulist_init = Ul.XD.init, tol = control$XD_conv_tol)
      Sigl.XD <- XD.out$Ulist <- lapply(XD.out$Ulist, function(M){ dimnames(M) <- list(phens_to_use, phens_to_use); M})
      pimat.XD <- t(XD.out$pi)
      XD_result <- list(XD.out = XD.out, mashdata.XD = mashdata.XD, Sigl.XD = Sigl.XD, pimat.XD = pimat.XD)
      save(XD_result, file = XD.output.file.namc)
      print(XD.output.file.namc)
    } else {
      load(file = XD.output.file.namc)
    }
    Sigl.XD <- lapply(XD_result$Sigl.XD, function(M) M[phens_to_use, phens_to_use])
    omegaseq.XD <- 1
    resl.store <- list()
    for(dat.type in dat.typev){
      out.post.mix.XD <- em.update.function(Y.em = Y.eml[[dat.type]][sams_for_model_testing, phens_to_use], S.em = S.eml[[dat.type]][sams_for_model_testing, phens_to_use],
                                            Sigl = Sigl.XD, R = R.init, omegaseq = omegaseq.XD,
                                            pimat = pimat.XD, meth = "post.mn")
      resl.store[[dat.type]] <- list(mn = out.post.mix.XD$mnmat[sams_for_model_testing, phens_to_use], sd = out.post.mix.XD$sdmat[sams_for_model_testing, phens_to_use],
                                     loglikv = out.post.mix.XD$loglikv[sams_for_model_testing],
                                     lfsr = out.post.mix.XD$lfsrmat[sams_for_model_testing, phens_to_use], Sigl = Sigl.XD, Sig.mn = NULL,
                                     Ksig = NULL, R = R.init, pimat = pimat.XD, omegaseq = omegaseq.XD)
    }
    save(resl.store, file = XD.resl.file.namc)
  }
  
  ##############################################
  # Run MASH
  if(Meth == "mash"){
    # load(file = file.objl)
    # Sigl.for.mash <- objl$impc_em.fit_nSig_1[[seed]]$Sigl[1]
    # names(Sigl.for.mash) <- "eb"
    
    replace.XD <- F
    if(!file.exists(XD.output.file.namc) | replace.XD){
      print("Running Extreme Deconvolution")
      mashdata.XD <- mashdata.big.effects.for.XD
      Ul.XD.init <- Sigl.init.bigeff.mash
      XD.out <- mashr:::bovy_wrapper(data = mashdata.XD, Ulist_init = Ul.XD.init, tol = control$XD_conv_tol_for_mash)
      Sigl.XD <- XD.out$Ulist <- lapply(XD.out$Ulist, function(M){ dimnames(M) <- list(phens_to_use, phens_to_use); M})
      pimat.XD <- t(XD.out$pi)
      XD_result <- list(XD.out = XD.out, mashdata.XD = mashdata.XD, Sigl.XD = Sigl.XD, pimat.XD = pimat.XD)
      save(XD_result, file = XD.output.file.namc)
      save(XD.out, mashdata.XD, Sigl.XD, pimat.XD, file = XD.output.file.namc)
      print(XD.output.file.namc)
    } else {
      print("Loading Extreme Deconvolution results")
      load(file = XD.output.file.namc)
    }
    Ul.XD.use <- XD.out$Ulist
    if(K2 > 1){
      sfa <- varimax(svd.Ztil$v[, 1:K2])
      sfa.F <- t(svd.Ztil$v[, 1:K2] %*% sfa$rotmat)
      sfa.L <- (svd.Ztil$u[, 1:K2, drop = F] %*% diag(svd.Ztil$d[1:K2], nrow = K2, ncol = K2)) %*% sfa$rotmat
    } else {
      sfa.F <- t(svd.Ztil$v[, 1:K2])
      sfa.L <- (svd.Ztil$u[, 1:K2, drop = F] %*% diag(svd.Ztil$d[1:K2], nrow = K2, ncol = K2))
    }
    Ul.rank1 <- list()
    for(i in 1:K2)
      Ul.rank1[[i]] <- t(sfa.L[, i] %*% sfa.F[i, , drop = F]) %*% sfa.L[, i] %*% sfa.F[i, , drop = F] / nrow(Ztil)
    Ul.data <- c(Ul.XD.use, Ul.rank1)
    names(Ul.data) <- paste0("U.", 1:(3 + K2))
    cov.meth <- c("identity", "equal_effects", "simple_het")
    cov.meth <- c(cov.meth, "singletons")
    Ul.canon <- cov_canonical(mashdata.for.model.fitting, cov_methods = cov.meth)
    # Ul.eb <- list(Sig)
    mashmethv <- list(c("data", "canon"), c("data"), c("eb"), c("canon"), c("data", "canon", "eb"))[1]#[c(1, 2, 3, 4)]
    if(!"resl" %in% ls())
      resl <- list()
    for(priorc in c("nullbiased", "uniform")[1]){#null.bias in c(T, F)){#priorc <- "nullbiased"#
    for(mashmeth in mashmethv){#mashmeth <- mashmethv[[1]]#
        if("eb" %in% mashmeth && (!"Sigl.for.mash" %in% ls() | is.null(Sigl.for.mash)))
          next
        Ul.all <- list()
        if("data" %in% mashmeth)
          Ul.all <- c(Ul.all, Ul.data)
        if("canon" %in% mashmeth)
          Ul.all <- c(Ul.all, Ul.canon)
        if("eb" %in% mashmeth)
          Ul.all <- c(Ul.all, Sigl.for.mash[1])#emout.mix$Sigl)
        res.mash.fitted.model <- mash(data = mashdata.for.model.fitting, Ulist = Ul.all, outputlevel = 2, usepointmass = use.pt.mass,
                                      prior = priorc, normalizeU = normU, add.mem.profile = F)#, grid = sqrt(omegaseq))
        if("eb" %in% mashmeth)
          eb.model <- res.mash.fitted.model
        mashdata.all.testing <- mash_set_data(Bhat = Y_zeroed[sams_for_model_testing, phens_to_use], Shat = S_zeroed[sams_for_model_testing, phens_to_use], alpha = 0, V = R.init)
        res.mash.all.testing <- mash(mashdata.all.testing, g = res.mash.fitted.model$fitted_g, fixg = T, outputlevel = 2, add.mem.profile = F)
        dimnames(res.mash.all.testing$vloglik) <- list(sams_for_model_testing, "llik")
        mashnam <- paste0("mash_", paste(mashmeth, collapse = "+"), "_prior_", priorc)
        resl[[mashnam]] <- list(mn = res.mash.all.testing$result$PosteriorMean[sams_for_model_testing, phens_to_use],
                                sd = res.mash.all.testing$result$PosteriorSD[sams_for_model_testing, phens_to_use],
                                loglikv = res.mash.all.testing$vloglik,
                                lfdr = res.mash.all.testing$result$lfdr[sams_for_model_testing, phens_to_use],
                                lfsr = res.mash.all.testing$result$lfsr[sams_for_model_testing, phens_to_use],
                                loglik = sum(res.mash.all.testing$vloglik[sams_for_lik_cross_val, ]))
        names(resl[[mashnam]]$loglikv) <- sams_for_model_testing
      }
    }
    phnam.mash <- colnames(res.mash.fitted.model$result$PosteriorMean)
    mash.omegaseq <- res.mash.fitted.model$fitted_g$grid^2
    mash.n.om <- length(mash.omegaseq)
    null.Sigl <- list(diag(rep(0, P)))
    names(null.Sigl) <- "null"
    mash.Sigl <- lapply(c(null.Sigl, res.mash.fitted.model$fitted_g$Ulist), function(M){ dimnames(M) <- list(phens_to_use, phens_to_use); M})
    mash.nSig <- length(mash.Sigl)
    mash.pimat.t <- matrix(0, mash.nSig, mash.n.om, 
                         dimnames = list(c("null", names(res.mash.fitted.model$fitted_g$Ulist)), 
                                                                    as.character(mash.omegaseq)))
    mash.piv <- res.mash.fitted.model$fitted_g$pi
    mash.pimat.t[2:mash.nSig, ] <- mash.piv[2:length(mash.piv)]
    mash.pimat.t[1, 1] <- mash.piv[1]
    mash.pimat.t.use <- mash.pimat.t[, , drop = F]
    mash.pimat.t.use <- mash.pimat.t.use / sum(mash.pimat.t.use)
    resl.store <- resl[names(resl) != "uv"]
    for(dat.type in dat.typev){#dat.type <- "raw"#
      out.post.mix.mash <- em.update.function(Y.em = Y.eml[[dat.type]][sams_for_model_testing, phens_to_use], 
                                                S.em = S.eml[[dat.type]][sams_for_model_testing, phens_to_use],
                                                Sigl = mash.Sigl, R = R.init,
                                                omegaseq = mash.omegaseq, prior.in.obj = F,
                                                pimat = t(mash.pimat.t.use), meth = "post.mn")
      
      resl.store[[dat.type]] <- list(mn = out.post.mix.mash$mnmat[sams_for_model_testing, phens_to_use], sd = out.post.mix.mash$sdmat[sams_for_model_testing, phens_to_use],
                                     loglikv = out.post.mix.mash$loglikv[sams_for_model_testing],
                                     lfsr = out.post.mix.mash$lfsrmat[sams_for_model_testing, phens_to_use], Sigl = mash.Sigl, Sig.mn = NULL,
                                     Ksig = NULL, R = R.init, pimat = t(mash.pimat.t.use), omegaseq = mash.omegaseq)
    }
    save(resl.store, file = mash.resl.file.namc)
    mash_results <- list(res.mash.all.testing = res.mash.all.testing, res.mash.fitted.model = res.mash.fitted.model)
    save(mash_results, file = mash.raw.results.file.namc)
  }
}
  # }
# }

