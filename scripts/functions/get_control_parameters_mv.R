
get_control_parameters_mv <- function(control = list()) {
  
  defaults <- list()
  
  ################################################
  # Parameters for simulations
  defaults$Nseq <- 500#c(100, 200, 500)#, 1000, 2000, 5000)
  defaults$Pseq <- c(10, 20)#, 40, 60, 100)
  defaults$n_subsamples_benchmark <- 10
  defaults$n_subsamples_main <- 50
  
  ######################################
  #Parameters for FDR control
  defaults$fdr.th <- .05
  defaults$fdr.th.seq <- c(seq(.1, 10, by = .1), Inf)
  
  #####################################################################
  #Parameters for EM algorithm for MAP estimation
  defaults$nam.truemut <- "trueMut"
  defaults$nam.negcon <- "negCon"
  # defaults$mv.meth.nam.use <- "impc_em.fit_nSig_1_fm_fa_bic_0_EMK_20"
  # defaults$estimation.meth <- c("map", "mcmc")[1]
  
  defaults$prior.sd.on.unobserved.thetahat <- 10
  defaults$nfac <- 20
  defaults$facnam <- paste("fac", 1:defaults$nfac, sep = "_")
  defaults$omegaseq <- exp((-10):10 * log(2))
  defaults$rank_deficient_R_eps <- .05
  
  
  defaults$mv_meth_nam_use <- "impc_MVphen_nSig_1_K_20"
  
  
  ########################################
  # Variables to be included in filename
  # defaults$variables_in_filename_MVphen <- c("Data", "Meth", "N", "P", "subsamseed", "nSig", "MVphen_K")
  # defaults$variables_in_filename_mash <- c("Data", "Meth", "N", "P", "subsamseed")
  # defaults$variables_in_filename_XD <- c("Data", "Meth", "N", "P", "subsamseed", "nSig", "XDmeth")
  defaults$variables_in_filename <- c("Data", "Meth", "N", "P", "subsamseed", "nSig", "MVphen_K", "XDmeth")
  # defaults$var.in.name <- c("N", "P", "nSig", "EMtol", "EMfm", "EMbic", "EMwish", "EMK", "EMKup", "subsamseed", "Data")
  # defaults$var.in.name.mash <- c("N", "P", "si", "bo", "EDtol", "subsamseed", "Data")
  # defaults$var.in.name.ed <- c("N", "P", "EDmeth", "EDtol", "nSig", "subsamseed", "Data")
  defaults$default_parameters <- defaults$default_parameters <- list()
  defaults$default_parameters$impc <- list(N = 2000, P = 148, N_demo = 200, P_demo = 20,
                                           variables_in_filename_MVphen = defaults$variables_in_filename_MVphen,
                                           variables_in_filename_mash = defaults$variables_in_filename_mash,
                                           variables_in_filename_XD = defaults$variables_in_filename_XD)
  defaults$default_parameters$eqtl <- list(N = 5000, P = 44, N_demo = 200, P_demo = 20,
                                           variables_in_filename_MVphen = defaults$variables_in_filename_MVphen,
                                           variables_in_filename_mash = defaults$variables_in_filename_mash,
                                           variables_in_filename_XD = defaults$variables_in_filename_XD)
  
  ###############################################
  # Our MVphen specific parameters
  defaults$MVphen_conv_tolerance <- 1e-4
  defaults$MVphen_conv_tolerance_default <- 1e-4
  
  ###############################
  # MASH specific parameters
  defaults$XD_conv_tol_for_mash <- 1e-5
  
  
  ########################################
  # Extreme Deconvolution specific  parameters
  defaults$XD_conv_tol <- 1e-4
  
  # if(!"Data" %in% ls())
  #   Data <- "impc"
  # if(!exists("full.analysis"))
  #   full.analysis <- T
  
  # defaults$uv.results.Y.S <- paste(uv.results.dir, "/uv_res_Y_S_phmap_linemap_", date.data.arrived, ".RData", sep = "")
  # file.copy(uv.results.Y.S, "raw_files/impc_data_in.RData")
  # getwd()
  
  ##########################################
  # Set RNG
  RNGkind(kind = "Mersenne-Twister", normal.kind =  "Inversion", sample.kind = "Rejection")
  
  #########################################
  # Heatmap colours
  # defaults$heat_col_palette <- c(rgb(0, 0, 1, alpha = seq(1, 0, len = 500)), rgb(1, 0, 0, alpha = seq(0, 1, len = 500)))
  # defaults$na_col <- grey(.5)
  defaults$heat_col_palette <- rainbow(n = 1000, start = 0, end = .66)[1000:1]
  
  #####################################################################
  # Directory structure
  defaults$data_dir <- "data"
  defaults$output_dir <- "output"
  defaults$global_res_dir <- paste0(defaults$output_dir, "/global_results")
  defaults$methods_comp_dir <- paste0(defaults$output_dir, "/methods_comparison")
  defaults$train_test_samples_dir <- paste0(defaults$output_dir, "/train_test_splits")
  # output_dir <- paste0(base.dir, "/RData_files/mv_results")
  defaults$figure_dir <- "figures"
  # defaults$text_numbers_dir <- "text_numbers"
  defaults$dropbox_figure_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Illuminating the mammalian genome with multivariate phenotype analysis/revision_figures"
  defaults$dropbox_text_numbers_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Illuminating the mammalian genome with multivariate phenotype analysis/revision_text_numbers"
  defaults$dropbox_table_dir <- "C:/Users/nicho/Dropbox/Apps/Overleaf/Illuminating the mammalian genome with multivariate phenotype analysis/revision_tables"
  
  
  #####################################################################
  # File names
  # defaults$file.runtab <- paste0(defaults$output_dir, "/runtab.RData")
  defaults$file.compl <- paste0(defaults$global_res_dir, "/global_compl.RDS")
  defaults$file.resl.comp <- paste0(defaults$global_res_dir, "/global_reslcomp.RDS")
  defaults$file.resll <- paste0(defaults$global_res_dir, "/resll_comb.RDS")
  defaults$file.resl.comp.fac <- paste0(defaults$global_res_dir, "/global_reslcompfac.RDS")
  defaults$file.objl <- paste0(defaults$global_res_dir, "/global_objl.RDS")
  defaults$file.resimp <- paste0(defaults$global_res_dir, "/resimp_comb.RDS")
  defaults$file_raw_factor_results <- paste0(defaults$global_res_dir, "/raw_factor_results.RDS")
  defaults$file_raw_factor_results_parallel_output <- paste0(defaults$global_res_dir, "/raw_factor_results_parallel.RDS")
  defaults$Data_all_file <- file.path(defaults$data_dir, "Data_all.RDS") 
  # defaults$file.glob.res <- paste0(defaults$global_res_dir, "/global_eb_results.RData")
  # defaults$file.glob.loadings <- paste0(defaults$global_res_dir, "/global_eb_loadings.RData")
  defaults$file.go.results <- paste0(defaults$global_res_dir, "/global_go_results.RData")
  # defaults$file.ebi.impc.results <- paste0(processed.data.dir, "/preprocessed_ebi_impc_results_v11.RData")
  # meth.comp.output.dir <- paste0(output_dir, "/meth_comp")
  
  control_out <- control
  for(par_input in names(defaults)) {
    if(!par_input %in% names(control)) {
      control_out[[par_input]] <- defaults[[par_input]]
    }
  }
  return(control_out)
}

