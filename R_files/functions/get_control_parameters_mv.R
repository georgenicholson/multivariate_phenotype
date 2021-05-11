
get_control_parameters_mv <- function(control = list()) {
  
  defaults <- list()
  ################################################
  # Parameters for simulations
  defaults$Nseq <- c(100, 200, 500, 1000, 2000, 5000)
  defaults$Pseq <- c(10, 20, 40, 60, 100)
  
  ######################################
  #Parameters for FDR control
  defaults$fdr.th <- .05
  defaults$fdr.th.seq <- c(seq(.1, 10, by = .1), Inf)
  
  #####################################################################
  #Parameters for EM algorithm for MAP estimation
  defaults$nam.truemut <- "trueMut"
  defaults$nam.negcon <- "negCon"
  defaults$mv.meth.nam.use <- "impc_em.fit_nSig_1_fm_fa_bic_0_EMK_20"
  defaults$estimation.meth <- c("map", "mcmc")[1]
  
  defaults$prior.sd.on.unobserved.thetahat <- 10
  defaults$nfac <- 20
  defaults$facnam <- paste("f", 1:nfac, sep = ".")
  
  defaults$var.in.name <- c("N", "P", "nSig", "EMtol", "EMfm", "EMbic", "EMwish", "EMK", "EMKup", "seed", "Data")
  defaults$var.in.name.mash <- c("N", "P", "si", "bo", "EDtol", "seed", "Data")
  defaults$var.in.name.ed <- c("N", "P", "EDmeth", "EDtol", "nSig", "seed", "Data")
  defaults$default.parameters <- list()
  defaults$default.parameters$impc <- list(N = 2000, P = 148,
                                  var.in.name = var.in.name,
                                  var.in.name.mash = var.in.name.mash,
                                  var.in.name.ed = var.in.name.ed)
  defaults$default.parameters$eqtl <- list(N = 5000, P = 44,
                                  var.in.name = var.in.name,
                                  var.in.name.mash = var.in.name.mash,
                                  var.in.name.ed = var.in.name.ed)
  # if(!"Data" %in% ls())
  #   Data <- "impc"
  # if(!exists("full.analysis"))
  #   full.analysis <- T
  
  # defaults$uv.results.Y.S <- paste(uv.results.dir, "/uv_res_Y_S_phmap_linemap_", date.data.arrived, ".RData", sep = "")
  # file.copy(uv.results.Y.S, "raw_files/impc_data_in.RData")
  # getwd()
  
  defaults$uv.results.Y.S <- "raw_files/impc_data_in.RData"
  defaults$procord.file <- paste(base.dir, "/RData_files/processed_data/procord.RData", sep = "")
  
  #####################################################################
  # Directories for updated EM and MASH
  mv.results.dir <- paste0(base.dir, "/RData_files/mv_results")
  file.runtab <- paste0(mv.results.dir, "/runtab.RData")
  global.res.dir <- paste0(mv.results.dir, "/global_results")
  file.compl <- paste0(global.res.dir, "/global_compl.RData")
  file.resl.comp <- paste0(global.res.dir, "/global_reslcomp.RData")
  file.resl.comp.fac <- paste0(global.res.dir, "/global_reslcompfac.RData")
  file.objl <- paste0(global.res.dir, "/global_objl.RData")
  file.glob.res <- paste0(global.res.dir, "/global_eb_results.RData")
  file.glob.loadings <- paste0(global.res.dir, "/global_eb_loadings.RData")
  file.go.results <- paste0(global.res.dir, "/global_go_results.RData")
  file.ebi.impc.results <- paste0(processed.data.dir, "/preprocessed_ebi_impc_results_v11.RData")
  # meth.comp.output.dir <- paste0(mv.results.dir, "/meth_comp")
  meth.comp.output.dir <- paste0(mv.results.dir, "/methods_comparison")
  for(dirc in c(mv.results.dir, meth.comp.output.dir, global.res.dir))
    dir.create(dirc, showWarnings = F)
  
  # sub.data.sets.dir <- paste0(base.dir, "/RData_files/sub_datasets", sep = "")
  sub.data.sets.dir <- paste0(base.dir, "/RData_files/sub_datasets_v2", sep = "")
  dir.create(sub.data.sets.dir, showWarnings = F)
  
  control_out <- control
  for(par_input in names(defaults)) {
    if(!par_input %in% names(control)) {
      control_out[[par_input]] <- defaults[[par_input]]
    }
  }
  return(control_out)
}

