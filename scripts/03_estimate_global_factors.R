rm(list = ls())
try({
  path_to_dir <- "C:/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/github_multivariate_phenotype/multivariate_phenotype"
  setwd(path_to_dir)
  renv::activate(path_to_dir)
  renv::restore(path_to_dir)
})
getwd()

##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control contains parameter settings
control <- get_control_parameters_mv()

resl.comp <- readRDS(file = control$file.resl.comp)
compl <- readRDS(file = control$file.compl)
objl <- readRDS(file = control$file.objl)

##########################################
# Load data
Data_all <- readRDS(file = control$Data_all_file)

# ##########################################
# # Get table of analyses
# analysis_table <- create_table_of_analyses(control = control, run_type = "main")

resl.comp.fac <- list()
Data <- "impc"
namc <-  control$mv_meth_nam_use
fac.meth <- "varimax"
facnam <- paste0("fac_", 1:control$nfac)
n_subsamples <- dim(compl[[namc]]$mnarr)[3]
ncore <- min(20, n_subsamples)
require(doParallel)
if(!"clust" %in% ls())
  clust <- parallel::makeCluster(rep("localhost", ncore), type = "SOCK")
doParallel::registerDoParallel(clust)
fac.res.store <- foreach::foreach(subsamseed = 1:n_subsamples, .verbose = T) %dopar% {
  sams.for.testing <- objl[[namc]][[subsamseed]]$saml$sams.for.testing
  meas.names <- rownames(objl[[namc]][[subsamseed]]$Sigl[[1]])
  facs <- switch(fac.meth, 
                 varimax = resl.comp[[namc]]$facs.varimax, 
                 promax = resl.comp[[namc]]$facs.promax)
  fac.out.post.mix <- em.update.function(Y.em = Data_all[[Data]]$Y_raw[sams.for.testing, meas.names], 
                                         S.em = Data_all[[Data]]$S_raw[sams.for.testing, meas.names],
                                         Sigl = objl[[namc]][[subsamseed]]$Sigl, 
                                         R = objl[[namc]][[subsamseed]]$R, 
                                         omegaseq = objl[[namc]][[subsamseed]]$omegaseq,
                                         pimat = objl[[namc]][[subsamseed]]$pimat, 
                                         meth = "post.mn.fac", 
                                         loadings = facs[meas.names, ])
  out <- list(mn = fac.out.post.mix$mnmat[sams.for.testing, facnam], 
              sd = fac.out.post.mix$sdmat[sams.for.testing, facnam],
              loglikv = fac.out.post.mix$loglikv[sams.for.testing],
              lfsr = fac.out.post.mix$lfsrmat[sams.for.testing, facnam], 
              loadings = facs)
  return(out)
}
# saveRDS(object = fac.res.store, file = )

suppressWarnings(lmat.norm <- exp(compl[[namc]]$llmat - apply(compl[[namc]]$llmat, 1, function(v) max(v, na.rm = T))))
pmix <- lmat.norm / rowSums(lmat.norm, na.rm = T)
fac.mnarr <- fac.sdarr <- fac.lfsrarr <- fac.wt.mnarr <- fac.wt.varr <- fac.wt.mnsqarr <- fac.wt.lfsrarr <- 
  array(NA, dim = c(Data_all[[Data]]$N_all, control$nfac, n_subsamples), dimnames = list(Data_all[[Data]]$sam_names, facnam, 1:n_subsamples))
for(subsamseed in 1:n_subsamples){
  fac.mnarr[objl[[namc]][[subsamseed]]$sams.for.testing, , subsamseed] <- fac.res.store[[subsamseed]]$mn
  fac.sdarr[objl[[namc]][[subsamseed]]$sams.for.testing, , subsamseed] <- fac.res.store[[subsamseed]]$sd
  fac.lfsrarr[objl[[namc]][[subsamseed]]$sams.for.testing, , subsamseed] <- fac.res.store[[subsamseed]]$lfsr
}
fac.wt.mnarr <- apply(sweep(fac.mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
fac.wt.varr <- apply(sweep(fac.sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
fac.wt.mnsqarr <- apply(sweep(fac.mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
fac.wt.lfsrarr <- apply(sweep(fac.lfsrarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
resl.comp.fac[[namc]] <- list(mn = fac.wt.mnarr, sd = sqrt(fac.wt.varr + fac.wt.mnsqarr - fac.wt.mnarr^2), lfsr = fac.wt.lfsrarr)
saveRDS(resl.comp.fac, file = control$file_raw_factor_results)

parallel::stopCluster(clust)

str(resl.comp.fac, m = 2)
