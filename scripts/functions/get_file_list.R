get_file_list <- function (control, file_core_name, subsamseed) {
  file_core_name_w_seed <- gsub("XXX", subsamseed, file_core_name)
  file_list <- list()
  file_list$emout.file.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_emout.RDS"))
  file_list$res.store.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_res.RDS"))
  file_list$resl.store.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_resl.RDS"))
  file_list$fac.res.store.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_facres.RDS"))
  file_list$loocv.res.store.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_loocv_res.RDS"))
  file_list$mash.resl.file.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_mash_resl.RDS"))
  file_list$mash.raw.results.file.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_mash_raw_results.RDS"))
  file_list$XD.output.file.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_XD_output.RDS"))
  file_list$XD.resl.file.namc <- file.path(control$methods_comp_dir, paste0(file_core_name_w_seed, "_XD_resl.RDS"))
  file_list
}
