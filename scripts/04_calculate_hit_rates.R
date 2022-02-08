rm(list = ls())
try({
  path_to_dir <- "C:/Users/nicho/Documents/GitHub/multivariate_phenotype"
  setwd(path_to_dir)
  # renv::activate(path_to_dir)
  # renv::restore(path_to_dir)
})
getwd()

Data <- "impc"

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
linemap <- Data_all$impc$linemap

resl.comp <- readRDS(file = control$file.resl.comp)
compl <- readRDS(file = control$file.compl)
objl <- readRDS(file = control$file.objl)
resl.comp.fac <- readRDS(file = control$file_raw_factor_results)

str(resl.comp.fac)

resl.err.rates.comb <- c(resl.comp[grepl(Data, names(resl.comp)) | names(resl.comp) == "uv"],
                          list(varimax = resl.comp.fac[[control$mv_meth_nam_use]]))
split.use <- 1
resl.err.rates.one.split <- c(list(uv = resl.comp$uv), lapply(compl[grepl(Data, names(compl))], 
                                   function(x) list(mn = x$mnarr[, , split.use], sd = x$sdarr[, , split.use], lfsr = x$lfsrarr[, , split.use])))



# str(resl.err.rates.comb[[1]])
# mean(resl.err.rates.comb$varimax$mn == 0)
# 
# tc_mash <- resl.err.rates.comb$impc_mash_nSig_1$mn / resl.err.rates.comb$impc_mash_nSig_1$sd
# tc_mvphen <- resl.err.rates.comb$impc_MVphen_nSig_1_K_20$mn / resl.err.rates.comb$impc_MVphen_nSig_1_K_20$sd
# tc_XD <- resl.err.rates.comb$impc_XD_nSig_1$mn / resl.err.rates.comb$impc_XD_nSig_1$sd
# str(tc)
# par(mfrow = c(2, 1))
# boxplot(tc[1:1000, ])
# boxplot(tc_mvphen[1:1000, ])
# 
# 
# negcons <- Data_all$impc$linemap$geno[Data_all$impc$linemap$line.type == control$nam.negcon]
# par(mfrow = c(3, 1))
# boxplot(tc[negcons[1:2000], ])
# abline(h = c(-2, 2))
# boxplot(tc_mvphen[negcons[1:2000], ])
# abline(h = c(-2, 2))
# boxplot(tc_XD[negcons[1:2000], ])
# abline(h = c(-2, 2))

##########################################
# Get table of analyses
analysis_table <- create_table_of_analyses(control = control, check_status = T, run_type = "benchmark")

file_list <- get_file_list(control = control, file_core_name = analysis_table[12, "file_core_name"], subsamseed = 1)

mash_results <- readRDS(file = file_list$mash.raw.results.file.namc)
str(mash_results)

for(err.data.type in c("comb", "single")[1]){
  resl.err.rates <- switch(err.data.type, 
                           comb = resl.err.rates.comb, 
                           single = resl.err.rates.one.split)
  
  out.perm <- err.rate.control(control = control, 
                               resl = resl.err.rates,
                               err.rate.meth = "perm", 
                               test.stat = "z", 
                               linemap = Data_all$impc$linemap, 
                               reflinemap = Data_all$impc$reflinemap, 
                               phmap = Data_all$impc$phmap, 
                               cenmap = Data_all$impc$cenmap,
                               Yhat = Data_all$impc$Y_raw,
                               control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                               err.thresh = .05, 
                               p.complete.null.true = 1, 
                               p.test.null.true = 1)
  print(out.perm$restab)
  out.perm.lfsr <- err.rate.control(control = control, 
                               resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],
                               err.rate.meth = "perm", 
                               test.stat = "lfsr", 
                               linemap = Data_all$impc$linemap, 
                               reflinemap = Data_all$impc$reflinemap, 
                               phmap = Data_all$impc$phmap, 
                               cenmap = Data_all$impc$cenmap,
                               Yhat = Data_all$impc$Y_raw,
                               control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                               err.thresh = .05, 
                               p.complete.null.true = 1, 
                               p.test.null.true = 1)
  
  print(out.perm.lfsr$restab)
  out.lfsr <- err.rate.control(control = control, 
                                    resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],
                                    err.rate.meth = "lfsr", 
                                    test.stat = "lfsr", 
                                    linemap = Data_all$impc$linemap, 
                                    reflinemap = Data_all$impc$reflinemap, 
                                    phmap = Data_all$impc$phmap, 
                                    cenmap = Data_all$impc$cenmap,
                                    Yhat = Data_all$impc$Y_raw,
                                    control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                    err.thresh = .05, 
                                    p.complete.null.true = 1, 
                                    p.test.null.true = 1)
  print(out.lfsr$restab)
  restabl <- list(perm = out.perm$restab, perm.lfsr = out.perm.lfsr$restab, lfsr = out.lfsr$restab)
  resll <- list(perm = out.perm$resl, perm.lfsr = out.perm.lfsr$resl, lfsr = out.lfsr$resl)
  resimpl <- list(perm = out.perm$resimp, perm.lfsr = out.perm.lfsr$resimp, lfsr = out.lfsr$resimp)
  
  saveRDS(resll, file = file.path(control$global_res_dir, paste0("resll_", err.data.type, ".RDS")))
  saveRDS(resimpl, file = file.path(control$global_res_dir, paste0("resimpl_", err.data.type, ".RDS")))
  saveRDS(restabl, file = file.path(control$global_res_dir, paste0("restabl_", err.data.type, ".RDS")))
}  

resimpl <- readRDS(file = file.path(control$global_res_dir, "resimpl_comb.RDS"))
resimp <- resimpl$perm
names(resimp)
resimp$eb.signsig <- resimp[, paste0(control$mv_meth_nam_use, ".perm.signsig")]
resimp$eb.sig <- abs(resimp$eb.signsig)
resimp$eb.t <- resimp[, paste0(control$mv_meth_nam_use, ".t")]
resimp$eb.th.final <- resimp[, paste0(control$mv_meth_nam_use, ".th.final")]
resimp <- resimp[resimp$line.type == control$nam.truemut, ]
resimp$uv.sig <- abs(resimp$uv.perm.signsig)
resimp$uv.signsig <- resimp$uv.perm.signsig
saveRDS(resimp, file = paste0(control$global_res_dir, "/resimp_comb.RDS"), version = 2)  


























