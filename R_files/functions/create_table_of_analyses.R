
create_table_of_analyses <- function(check_status = T) {
  ########################################
  # Variables to be included in filename
  var.in.name <- c("N", "P", "nSig", "EMtol", "EMfm", "EMbic", "EMwish", "EMK", "EMKup", "seed", "Data")
  var.in.name.mash <- c("N", "P", "si", "bo", "EDtol", "seed", "Data")
  var.in.name.ed <- c("N", "P", "EDmeth", "EDtol", "nSig", "seed", "Data")
  
  ####################################################
  # General parameters to pass to run_EM.R
  #############################################################
  data_sets_to_analyse <- c("impc", "eqtl")[1:2]
  methods_to_apply <- c("em.fit", "ed", "mash", "em.fit.rand", "em.fit.N.500", 
                 "em.test", "em.loocv", "em.fac", "em.all.test", "em.all")[1:5]
  N_values_to_run <- NA# c(100, 200, 500, 1000, 2000, 5000)[3]
  P_values_to_run <- NA# c(10, 20, 40, 60, 100, 148)[6]
  perform_full_data_analysis <- T
  n_Sig_mix_to_run <- 1:2
  n_seed_to_run <- 10
  EM_dim_fac_space_to_run <- c(15, 20, 30, 40)
  EM_model_type_to_run <- c("fa", "pca")[1]
  EM_dim_fac_space_to_use_for_main_analysis <- 20
  
  EM_BIC_penalty <- 0
  EM_use_wishart_prior <- F
  EM_update_fac_space_dim <- F

  ###############################################
  # Our EM specific parameters
  EM_conv_tolerance_to_run <- c(1e-4, 1e-5, 5e-6)[1]

  ###############################
  # MASH specific parameters
  ED_conv_tol_for_mash <- 1e-5
  mash_include_ED_components <- 1
  mash_include_singletons <- 1
  
  ########################################
  # Extreme Deconvolution specific  parameters
  ED_analyses_to_run <- c("mash", "justED")[2]
  ED_conv_tol <- c(1e-3, 1e-4, 1e-5, 1e-6)[2]
  
  
  
  #############################
  # Create runtab
  runtab <- expand.grid(N = N_values_to_run, P = P_values_to_run, da = data_sets_to_analyse, 
                        me = methods_to_apply, fu = perform_full_data_analysis, 
                        nSig = n_Sig_mix_to_run, EDmeth = ED_analyses_to_run, 
                        EDtol = ED_conv_tol, EMtol = EM_conv_tolerance_to_run, EMfm = EM_model_type_to_run, 
                        EMK = EM_dim_fac_space_to_run, 
                        EMbic = EM_BIC_penalty, EMwish = EM_use_wishart_prior, 
                        EMKup = EM_update_fac_space_dim, bo = mash_include_ED_components, 
                        si = mash_include_singletons, stringsAsFactors = F)
  runtab[which(!grepl("em", runtab$me)), c("EMtol", "EMfm", "EMbic",  "EMK", "EMKup")] <- NA
  
  #############################
  # Filter runtab
  #############################
  runtab[!runtab$me %in% c("ed", "mash"), c("EDtol", "EDmeth")] <- NA
  runtab[runtab$me == "mash", "EDmeth"] <- "mash"
  runtab[runtab$me == "mash", "EDtol"] <- ED_conv_tol_for_mash
  runtab[runtab$me == "ed" & runtab$EDmeth == "mash", "EDtol"] <- ED_conv_tol_for_mash
  runtab[runtab$me == "ed" & runtab$EDmeth == "mash", c("nSig")] <- 1
  runtab[runtab$me == "mash", c("nSig")] <- 1
  runtab[!(runtab$me == "em.fit" & runtab$nSig == 1), "fu"] <- T
  runtab[runtab$me == "em.fit" & runtab$nSig == 1, "EMtol"] <- 1e-4#5e-6
  runtab[runtab$fu & runtab$da == "impc", "P"] <- 148
  runtab[runtab$fu & runtab$da == "impc", "N"] <- 2000
  runtab[runtab$fu & runtab$da == "eqtl", "P"] <- 44
  runtab[runtab$fu & runtab$da == "eqtl", "N"] <- 5000
  runtab <- unique(runtab)
  runtab$ntimes <- n_seed_to_run
  runtab[runtab$me == "em.fit.N.500", "N"] <- 500
  runtab <- runtab[!(runtab$me == "em.fit.N.500" & (runtab$da == "eqtl" | runtab$nSig > 1)), ]
  runtab <- runtab[!(runtab$me == "em.fit.rand" & (runtab$da == "eqtl")), ]
  runtab <- runtab[!(runtab$me %in% c("em.fit.rand", "em.fit.N.500") & runtab$EMK != EM_dim_fac_space_to_use_for_main_analysis), ]
  runtab <- runtab[!(runtab$me == "em.fit" & runtab$nSig > 1 & runtab$EMK != EM_dim_fac_space_to_use_for_main_analysis), ]
  runtab <- runtab[!(runtab$me == "em.fit.rand" & runtab$nSig > 1), ]
  
  ###############################################################
  #Specify memory requirements for each type of job
  #################################################################
  runtab$mem <- 2300
  runtab[runtab$da %in% c("impc", "eqtl") & runtab$me == "em.fit", "mem"] <- 2000
  runtab[runtab$da  == "impc" & runtab$me == "mash" & runtab$fu, "mem"] <- 6000
  runtab[runtab$da == "eqtl" & runtab$me == "mash" & runtab$fu, "mem"] <- 4600
  
  runtab <- unique(runtab)
  runtab <- runtab[order(runtab$da, match(runtab$me, methods_to_apply), runtab$nSig), ]
  runtab$rand <- ifelse(grepl("rand", runtab$me), T, F)
  runtab$loocv <- ifelse(runtab$me == "em.fit" & runtab$nSig == 1 & runtab$N == 2000, T, F)
  rownames(runtab) <- 1:nrow(runtab)
  runtab$Data <- runtab$da
  
  #############################################
  # Check when each scenario was last run
  ###########################################
  if(check_status) {
    fcheckfields <- c("emout.ok", "emout.st", "emout.en", "res.ok", "res.st", "res.en", "loocv.ok", 
                       "loocv.st", "loocv.en", "fac.ok", "fac.st", "fac.en")
    runtab[, fcheckfields] <- NA
    runtab[, c("file.base", "meth")] <- NA
    for(scen in 1:nrow(runtab)){
      print(scen)
      for(j in 1:ncol(runtab))
        assign(colnames(runtab)[j], runtab[scen, j], pos = sys.frame())
      Data <- da
      if(grepl("em.fit", me)){
        var.in.name.use <- var.in.name
        if(rand) {
          var.in.name.use <- c(var.in.name.use, "rand")
        }
        meth <- "eb"
      }
      if(me == "mash"){
        var.in.name.use <- var.in.name.mash
        meth <- "mash"
      }
      if(me == "ed"){
        var.in.name.use <- var.in.name.ed
        meth <- "XD"
      }
      
      file.base <- c()
      for(seed in 1:n_seed_to_run) {
        file.base <- c(file.base, paste0(paste(paste(var.in.name.use, 
                                           sapply(var.in.name.use, function(x) get(x, envir = environment())), sep = "_"), collapse = "_")))
      }
      print("made it")
      loocv.res.store.namc <- paste0(meth.comp.output.dir, "/", file.base, "_loocv_res.RData")
      fac.res.store.namc <- paste0(meth.comp.output.dir, "/", file.base, "_facres.RData")
      res.store.namc <- paste0(meth.comp.output.dir, "/", file.base, switch(meth, eb = "_res.RData", mash = "_mash_resl.RData", XD = "_bovy_resl.RData"))
      emout.file.namc <- paste0(meth.comp.output.dir, "/", file.base, switch(meth, eb = "_emout.RData", mash = NA, XD = NA))
      runtab[scen, fcheckfields] <- list(all(file.exists(emout.file.namc)),
                                         format(min(file.info(emout.file.namc)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(emout.file.namc)[, "ctime"]), format = "%d-%b"),
                                         all(file.exists(res.store.namc)),
                                         format(min(file.info(res.store.namc)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(res.store.namc)[, "ctime"]), format = "%d-%b"),
                                         all(file.exists(loocv.res.store.namc)),
                                         format(min(file.info(loocv.res.store.namc)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(loocv.res.store.namc)[, "ctime"]), format = "%d-%b"),
                                         all(file.exists(fac.res.store.namc)),
                                         format(min(file.info(fac.res.store.namc)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(fac.res.store.namc)[, "ctime"]), format = "%d-%b"))
      runtab$meth[scen] <- meth
      runtab$file.base[scen] <- gsub("seed\\_1", "seed\\_XXX", file.base[1])
    }
  }

  return(runtab)  
}




# runtab <- expand.grid(N = N_values_to_run, P = P_values_to_run, bo = bovyseq, si = singleseq, da = data_sets_to_analyse, 
#                       me = methods_to_apply, ma = max.test.seq, fu = perform_full_data_analysis, 
#                       ss = sexspecific.seq, nSig = n_Sig_mix_to_run, EDmeth = ED_analyses_to_run, 
#                       EDtol = ED_conv_tol, EMtol = EM_conv_tolerance_to_run, EMfm = EM_model_type_to_run, 
#                       EMbic = EMbic.seq, EMwish = EMwish.seq, EMmash = EMmash.seq, 
#                       EMK = EM_dim_fac_space_to_run, EMKup = EMKup.seq, stringsAsFactors = F)
# sex_specific_analysis <- F
# max.test.seq <- 20000 # ??
# EMmash.seq <- F
# EMbic.seq <- c(0, .5, 1, 2, 3)[1]
# EMwish.seq <- F#c(0, .25, .5)
# K.use <- 20
# EMKup.seq <- F
# bovyseq <- 1
# singleseq <- 1
