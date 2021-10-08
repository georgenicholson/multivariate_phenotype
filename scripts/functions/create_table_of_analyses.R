
create_table_of_analyses <- function(control, check_status = T, run_type = c("demo", "main", "benchmark", "test_benchmark")[2]) {
  
  if(run_type == "demo") {
    runtab <- data.frame(N = 100, P = 10, Data = "impc", 
                          Meth = c("MVphen", "XD", "mash"), 
                          nSig = 1, 
                          MVphen_K = 3, 
                          stringsAsFactors = F)
  }  
  
  if(run_type == "main") {
    runtab <- data.frame(N = 2000, P = 148, Data = "impc", 
                         Meth = "MVphen", 
                         nSig = 1, 
                         MVphen_K = 20, 
                         stringsAsFactors = F)
  }  
  
  
  if (run_type %in% c("benchmark", "test_benchmark")) {
    runtab <- expand.grid(N = NA, P = NA, Data = c("impc", "eqtl"), 
                          Meth = c("MVphen", "XD", "mash", "MVphen_rand", "MVphen_N_500"), 
                          nSig = 1:2, 
                          MVphen_K = c(15, 20, 30, 40), 
                          stringsAsFactors = F)
    runtab[which(!grepl("MVphen", runtab$Meth)), "MVphen_K"] <- NA
  
    #############################
    # Filter runtab
    #############################
    runtab[runtab$Meth == "mash", c("nSig")] <- 1
    runtab[runtab$Data == "impc", "P"] <- control$default_parameters$impc$P
    runtab[runtab$Data == "impc", "N"] <- control$default_parameters$impc$N
    runtab[runtab$Data == "eqtl", "P"] <- control$default_parameters$eqtl$P
    runtab[runtab$Data == "eqtl", "N"] <- control$default_parameters$eqtl$N
    runtab$n_subsamples <- control$n_subsamples
    runtab[runtab$Meth == "MVphen_N_500", "N"] <- 500
    runtab <- runtab[which(!(runtab$Meth == "MVphen_N_500" & (runtab$Data == "eqtl" | runtab$nSig > 1 | runtab$MVphen_K != 20))), ]
    runtab <- runtab[which(!(runtab$Meth == "MVphen_rand" & (runtab$Data == "eqtl" | runtab$MVphen_K != 20))), ]
    runtab <- runtab[which(!(runtab$Meth == "MVphen_rand" & runtab$nSig > 1)), ]
    
    ###############################################################
    #Specify memory requirements for each type of job
    #################################################################
    runtab$mem <- 2300
    runtab[runtab$Data %in% c("impc", "eqtl") & runtab$Meth == "MVphen", "mem"] <- 2000
    runtab[runtab$Data  == "impc" & runtab$Meth == "mash", "mem"] <- 6000
    runtab[runtab$Data == "eqtl" & runtab$Meth == "mash", "mem"] <- 4600
    runtab <- runtab[order(runtab$Data, runtab$Meth, runtab$nSig), ]
    runtab$rand <- ifelse(grepl("rand", runtab$Meth), T, F)
    runtab$loocv <- ifelse(runtab$Meth == "MVphen" & runtab$nSig == 1 & runtab$N == 2000, T, F)
    rownames(runtab) <- 1:nrow(runtab)
    if (run_type == "test_benchmark") {
      runtab$n_subsamples <- 1  
      runtab$N <- 200  
      runtab$P <- 20
      runtab$MVphen_K <- 5
    }
    runtab <- unique(runtab)
  }  
  runtab
  
  #############################################
  # Check when each scenario was last run
  ###########################################
  if(check_status) {
    fcheckfields <- c("emout_ok", "emout_st", "emout_en", "res_ok", "res_st", "res_en", "loocv_ok", 
                       "loocv_st", "loocv_en", "fac_ok", "fac_st", "fac_en")
    runtab[, fcheckfields] <- NA
    runtab[, c("file_core_name", "meth")] <- NA
    for(scen in 1:nrow(runtab)){
      # print(scen)
      for(j in 1:ncol(runtab))
        assign(colnames(runtab)[j], runtab[scen, j], pos = sys.frame(which = 0))
      Data <- da
      if(grepl("MVphen", Meth)){
          variables_in_filename_use <- control$variables_in_filename_MVphen
        if(rand) {
          variables_in_filename_use <- c(variables_in_filename_use, "rand")
        }
      }
      if(Meth == "mash"){
        variables_in_filename_use <- control$variables_in_filename_mash
      }
      if(Meth == "XD"){
        variables_in_filename_use <- control$variables_in_filename_XD
      }
      
      file_core_name <- c()
      for(seed in 1:n_seed_to_run) {
        file_core_name <- c(file_core_name, paste0(paste(paste(variables_in_filename_use, 
                                           sapply(variables_in_filename_use, function(x) get(x, envir = environment())), sep = "_"), collapse = "_")))
      }
      # print("made it")
      loocv_results_filename <- paste0(control$methods_comp_dir, "/", file_core_name, "_loocv_res.RData")
      factor_results_filename <- paste0(control$methods_comp_dir, "/", file_core_name, "_facres.RData")
      basic_results_filename <- paste0(control$methods_comp_dir, "/", file_core_name, 
                               switch(Meth, eb = "_res.RData", mash = "_mash_resl.RData", XD = "_bovy_resl.RData"))
      emout_file_name <- paste0(control$methods_comp_dir, "/", file_core_name, 
                                switch(Meth, eb = "_emout.RData", mash = NA, XD = NA))
      runtab[scen, fcheckfields] <- list(all(file.exists(emout_file_name)),
                                         format(min(file.info(emout_file_name)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(emout_file_name)[, "ctime"]), format = "%d-%b"),
                                         all(file.exists(basic_results_filename)),
                                         format(min(file.info(basic_results_filename)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(basic_results_filename)[, "ctime"]), format = "%d-%b"),
                                         all(file.exists(loocv_results_filename)),
                                         format(min(file.info(loocv_results_filename)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(loocv_results_filename)[, "ctime"]), format = "%d-%b"),
                                         all(file.exists(factor_results_filename)),
                                         format(min(file.info(factor_results_filename)[, "ctime"]), format = "%d-%b"),
                                         format(max(file.info(factor_results_filename)[, "ctime"]), format = "%d-%b"))
      runtab$file_core_name[scen] <- gsub("seed\\_1", "seed\\_XXX", file_core_name[1])
    }
  }
  rownames(runtab) <- 1:nrow(runtab)
  return(runtab)  
}


