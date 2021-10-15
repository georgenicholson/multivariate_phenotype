variables_in_filename_use <- control$variables_in_filename
file.base <- paste0(control$methods_comp_dir, "/", 
                    paste(paste(variables_in_filename_use, sapply(variables_in_filename_use, get), sep = "_"), collapse = "_"))

emout.file.namc <- paste0(file.base, "_emout.RDS")
res.store.namc <- paste0(file.base, "_res.RDS")
resl.store.namc <- paste0(file.base, "_resl.RDS")
fac.res.store.namc <- paste0(file.base, "_facres.RDS")
loocv.res.store.namc <- paste0(file.base, "_loocv_res.RDS")

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