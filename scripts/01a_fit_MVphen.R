
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