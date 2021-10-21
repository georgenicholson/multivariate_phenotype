



require(mashr)
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
  snps.strong.use <- sams_for_strong_cov_est
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
  if(!file.exists(file_list$XD.output.file.namc) | replace.XD){
    XD.out <- mashr:::bovy_wrapper(data = mashdata.XD, Ulist_init = Ul.XD.init, tol = control$XD_conv_tol)
    Sigl.XD <- XD.out$Ulist <- lapply(XD.out$Ulist, function(M){ dimnames(M) <- list(phens_to_use, phens_to_use); M})
    pimat.XD <- t(XD.out$pi)
    XD_result <- list(XD.out = XD.out, mashdata.XD = mashdata.XD, Sigl.XD = Sigl.XD, pimat.XD = pimat.XD)
    saveRDS(XD_result, file = file_list$XD.output.file.namc)
  } else {
    XD_result <- readRDS(file = file_list$XD.output.file.namc)
  }
  Sigl.XD <- lapply(XD_result$Sigl.XD, function(M) M[phens_to_use, phens_to_use])
  omegaseq.XD <- 1
  resl.store <- list()
  for(dat.type in dat.typev){
    out.post.mix.XD <- em.update.function(Y.em = Y.eml[[dat.type]][sams_for_model_testing, phens_to_use], 
                                          S.em = S.eml[[dat.type]][sams_for_model_testing, phens_to_use],
                                          Sigl = Sigl.XD, 
                                          R = R.init, 
                                          omegaseq = omegaseq.XD,
                                          pimat = pimat.XD, 
                                          meth = "post.mn")
    resl.store[[dat.type]] <- list(mn = out.post.mix.XD$mnmat[sams_for_model_testing, phens_to_use], 
                                   sd = out.post.mix.XD$sdmat[sams_for_model_testing, phens_to_use],
                                   loglikv = out.post.mix.XD$loglikv[sams_for_model_testing],
                                   lfsr = out.post.mix.XD$lfsrmat[sams_for_model_testing, phens_to_use], 
                                   Sigl = Sigl.XD, 
                                   Sig.mn = NULL,
                                   Ksig = NULL, 
                                   R = R.init, 
                                   pimat = pimat.XD, 
                                   omegaseq = omegaseq.XD)
  }
  saveRDS(resl.store, file = file_list$XD.resl.file.namc)
}

##############################################
# Run MASH
if(Meth == "mash"){
  replace.XD <- F
  if(!file.exists(file_list$XD.output.file.namc) | replace.XD){
    print("Running Extreme Deconvolution")
    mashdata.XD <- mashdata.big.effects.for.XD
    Ul.XD.init <- Sigl.init.bigeff.mash
    XD.out <- mashr:::bovy_wrapper(data = mashdata.XD, Ulist_init = Ul.XD.init, tol = control$XD_conv_tol_for_mash)
    Sigl.XD <- XD.out$Ulist <- lapply(XD.out$Ulist, function(M){ dimnames(M) <- list(phens_to_use, phens_to_use); M})
    pimat.XD <- t(XD.out$pi)
    XD_result <- list(XD.out = XD.out, mashdata.XD = mashdata.XD, Sigl.XD = Sigl.XD, pimat.XD = pimat.XD)
    saveRDS(XD_result, file = file_list$XD.output.file.namc)
  } else {
    print("Loading Extreme Deconvolution results")
    print(file_list$XD.output.file.namc)
    print(Data)
    XD_result <- readRDS(file = file_list$XD.output.file.namc)
  }
  Ul.XD.use <- XD_result$XD.out$Ulist
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
  resl <- list()
  mashmeth <- c("data", "canon")
  Ul.all <- c(Ul.data, Ul.canon)
  res.mash.fitted.model <- mash(data = mashdata.for.model.fitting, 
                                Ulist = Ul.all, 
                                outputlevel = 2, 
                                add.mem.profile = F)
  mashdata.all.testing <- mash_set_data(Bhat = Y_zeroed[sams_for_model_testing, phens_to_use], 
                                        Shat = S_zeroed[sams_for_model_testing, phens_to_use], 
                                        alpha = 0, 
                                        V = R.init)
  res.mash.all.testing <- mash(mashdata.all.testing, 
                               g = res.mash.fitted.model$fitted_g, 
                               fixg = T, 
                               outputlevel = 2, 
                               add.mem.profile = F)
  dimnames(res.mash.all.testing$vloglik) <- list(sams_for_model_testing, "llik")
  mashnam <- paste0("mash_", paste(mashmeth, collapse = "+"))
  resl[[mashnam]] <- list(mn = res.mash.all.testing$result$PosteriorMean[sams_for_model_testing, phens_to_use],
                          sd = res.mash.all.testing$result$PosteriorSD[sams_for_model_testing, phens_to_use],
                          loglikv = res.mash.all.testing$vloglik,
                          lfdr = res.mash.all.testing$result$lfdr[sams_for_model_testing, phens_to_use],
                          lfsr = res.mash.all.testing$result$lfsr[sams_for_model_testing, phens_to_use],
                          loglik = sum(res.mash.all.testing$vloglik[sams_for_lik_cross_val, ]))
  names(resl[[mashnam]]$loglikv) <- sams_for_model_testing
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
    
    resl.store[[dat.type]] <- list(mn = out.post.mix.mash$mnmat[sams_for_model_testing, phens_to_use], 
                                   sd = out.post.mix.mash$sdmat[sams_for_model_testing, phens_to_use],
                                   loglikv = out.post.mix.mash$loglikv[sams_for_model_testing],
                                   lfsr = out.post.mix.mash$lfsrmat[sams_for_model_testing, phens_to_use], 
                                   Sigl = mash.Sigl, 
                                   Sig.mn = NULL,
                                   Ksig = NULL, 
                                   R = R.init, 
                                   pimat = t(mash.pimat.t.use), 
                                   omegaseq = mash.omegaseq)
  }
  saveRDS(resl.store, file = file_list$mash.resl.file.namc)
  mash_results <- list(res.mash.all.testing = res.mash.all.testing, res.mash.fitted.model = res.mash.fitted.model)
  saveRDS(mash_results, file = file_list$mash.raw.results.file.namc)
}
