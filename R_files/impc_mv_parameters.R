send.to.latex <- function(){
  require("KeyboardSimulator")
  keybd.press('ctrl+alt+u')
} 

linux <- Sys.info()["sysname"] == "Linux"
# platdir <- ifelse(linux, "/mnt/x", "X:")
if(Sys.info()["nodename"] == "TW")
  platdir <- paste0(ifelse(linux, "/mnt/c", "C:"), "/Users/nicho/Documents/bauer_sync")
if(Sys.info()["nodename"] != "TW")
  platdir <- ifelse(linux, "/mnt/x", "X:")
#########################################################
#Load packages for R
# if(Sys.info()["nodename"] %in% c("NICHOLSO-LAPTOP", "BAUER") & !linux)
#   .libPaths(paste0(platdir, "/R_libraries"))
libs.need <- c("doParallel", "Rcpp", "Matrix", "XML", "MASS", "httr")#, "mclust")
# libs.install <- libs.need[!libs.need %in% library()$results[, 1]]
# for(libc in libs.install)
#   install.packages(libc)
for(libc in libs.need)
  library(libc, character.only = T)


#########################################################
#Directory structure
base.dir <- paste0(platdir, "/projects/impc_mv_analysis")
R.file.dir <- paste0(base.dir, "/R_files", sep = "")
bash.file.dir <- paste0(base.dir, "/bash_files", sep = "")
dir.create(bash.file.dir, showWarnings = F)
#########################################################
#Dropbox directories
dropdir <- ifelse(linux, "/mnt/c/Users/*/Dropbox", "C:/Users/nicho/Dropbox")
chris.pres.dropbox <- paste0(dropdir, "/chris_isba_slides_2018_june/figures")
paper.figures.dropbox <- paste0(dropdir, "/impc_mv_paper/tex_files/figures")
text.numbers.dropbox <- paste0(dropdir, "/impc_mv_paper/tex_files/text_numbers")
chris.pres.numbers.dropbox <- paste0(dropdir, "/chris_isba_slides_2018_june/numbers")
revision.paper.figures.dropbox <- paste0(dropdir, "/impc_mv_paper/tex_files/revision_figures")
revision.text.numbers.dropbox <- paste0(dropdir, "/impc_mv_paper/tex_files/revision_text_numbers")
revision.paper.figures <- paste0(base.dir, "/PLOS_bio_revision_tex_files/revision_figures")
revision.text.numbers <- paste0(base.dir, "/PLOS_bio_revision_tex_files/revision_text_numbers")
dir.create(revision.paper.figures.dropbox, showWarnings = F)
dir.create(revision.text.numbers.dropbox, showWarnings = F)
dir.create(revision.paper.figures, showWarnings = F)
dir.create(revision.text.numbers, showWarnings = F)


#########################################################
#Definitions for raw data import
options(stringsAsFactors = F)
raw.data.in.dir <- paste(base.dir, "/data_in", sep = "")
date.data.arrived <- c("2018-03-26", "2019-04-23")[1]
if(date.data.arrived == "2018-03-26")
  raw.data.file <- paste(raw.data.in.dir, "/impcDataForGeorge_26-3-18.csv", sep = "")
if(date.data.arrived == "2019-04-23")
  raw.data.file <- paste(raw.data.in.dir, "/impcDataForGeorge_23-4-19.csv", sep = "")

#########################################################
#Definitions for initially processed data
processed.data.dir <- paste(base.dir, "/RData_files/processed_data", sep = "")
data.in.filename <- paste(processed.data.dir, "/", date.data.arrived, "_impc_data_in.RData", sep = "")
negative.controls.filename <- paste(processed.data.dir, "/", date.data.arrived, "_negative_controls.RData", sep = "")

#####################################################
#Directories for storing UV results
uv.results.dir <- paste(base.dir, "/RData_files/uv_results", sep = "")
dir.num <- paste(uv.results.dir, "/numerical_", date.data.arrived, sep = "")
dir.cat <- paste(uv.results.dir, "/categorical_", date.data.arrived, sep = "")
resmat.dir.num <- paste(dir.num, "/processed_resmats", sep = "")
raw.results.dir.num <- paste(dir.num, "/raw_results", sep = "")
raw.results.errors.dir.num <- paste(dir.num, "/raw_results_errors", sep = "")
# raw.results.sexspecific.dir.num <- paste(dir.num, "/raw_results_sexspecific", sep = "")
# raw.results.errors.sexspecific.dir.num <- paste(dir.num, "/raw_results_sexspecific_errors", sep = "")
raw.results.sexspecific.dir.num <- paste(dir.num, "/raw_results_sexgeno_inter", sep = "")
raw.results.errors.sexspecific.dir.num <- paste(dir.num, "/raw_results_sexgeno_inter_errors", sep = "")
uv.resmat.file <- paste(uv.results.dir, "/univariate_resmat_", date.data.arrived, ".RData", sep = "")
sexspecific.uv.resmat.file <- paste(uv.results.dir, "/sexspecific_univariate_resmat_", date.data.arrived, ".RData", sep = "")
uv.resmat.file.post.qc <- paste(uv.results.dir, "/post_qc_univariate_resmat_", date.data.arrived, ".RData", sep = "")
uv.resmat.file.pre.qc <- paste(uv.results.dir, "/pre_qc_univariate_resmat_", date.data.arrived, ".RData", sep = "")
uv.resmat.post.qc.with.sig.thresholds.file <- paste(uv.results.dir, "/uv_resmat_post_qc_with_sig_thresholds_", date.data.arrived, ".RData", sep = "") 
uv.resmat.pre.qc.with.sig.thresholds.file <- paste(uv.results.dir, "/uv_resmat_pre_qc_with_sig_thresholds_", date.data.arrived, ".RData", sep = "") 
uv.results.Y.S <- paste(uv.results.dir, "/uv_res_Y_S_phmap_linemap_", date.data.arrived, ".RData", sep = "")
for(dirc in c(uv.results.dir, dir.num, dir.cat, resmat.dir.num, raw.results.dir.num))
  dir.create(dirc, showWarnings = F)

#####################################################################
#Parameters for UV analysis
trans = T
fit.spline = T
burnin = 50000
nits = 550000
thin = 1000


################################################
# Parameters for simulations
Nseq <- c(100, 200, 500, 1000, 2000, 5000)
Pseq <- c(10, 20, 40, 60, 100)
# sub.data.sets.dir <- paste0(base.dir, "/RData_files/sub_datasets", sep = "")
sub.data.sets.dir <- paste0(base.dir, "/RData_files/sub_datasets_v2", sep = "")
dir.create(sub.data.sets.dir, showWarnings = F)

######################################
#Parameters for FDR control
fdr.th <- .05
fdr.th.seq <- c(seq(.1, 10, by = .1), Inf)

####################################################################  

#Heuristics for QC of UV analysis
#Examine UV heatmap outputted by plot_annotation_heatmaps.R
#(cen, proc) paris to remove, due to obvious temporal effects
nline.min <- 500
omitl <- list(list(cen = c("Bcm"), proc = c("Auditory Brain Stem Response"), omit = T),
              list(cen = c("Gmc"), proc = c("Eye Morphology"), omit = T),
              list(cen = c("Gmc"), proc = c("Echo"), omit = T),
              list(cen = c("Gmc"), proc = c("FACS"), omit = T),
              list(cen = c("Ics"), proc = c("Hematology"), omit = T),
              list(cen = c("Ics"), proc = c("Open Field"), omit = T),
              list(cen = c("Ics"), proc = c("Acoustic Startle and Pre-pulse Inhibition (PPI)"), omit = F),
              list(cen = c("J"), proc = c("Clinical Chemistry"), omit = T),
              list(cen = c("J"), proc = c("FACS"), omit = T),
              list(cen = c("J"), proc = c("Hematology"), omit = T),
              list(cen = c("J"), proc = c("Heart Weight"), omit = T),
              list(cen = c("J"), proc = c("Electroconvulsive Threshold Testing"), omit = T),
              list(cen = c("Tcp"), proc = c("FACS"), omit = T),
              list(cen = c("Tcp"), proc = c("Combined SHIRPA and Dysmorphology"), omit = T),
              list(cen = c("Tcp"), proc = c("Acoustic Startle and Pre-pulse Inhibition (PPI)"), omit = T),
              list(cen = c("Tcp"), proc = c("Hematology"), omit = T),
              list(cen = c("Tcp"), proc = c("Open Field"), omit = T),
              list(cen = c("Wtsi"), proc = c("Grip Strength"), omit = T))


#####################################################################
#Parameters for EM algorithm for MAP estimation
methnam <- c("MVphen", "MV")[2]
nam.truemut <- "trueMut"
nam.negcon <- "negCon"
mv.meth.nam.use <- "impc_em.fit_nSig_1_fm_fa_bic_0_EMK_20"
estimation.meth <- c("map", "mcmc")[1]
uv.scaling <- c("simple", "scaledInvWish")[1]
runEM <- F
scale.uv.pre <- F
init.R.at <- c("prev", "posDefMode", "all")[2]
# prior.sd.on.unobserved.thetahat <- 1
prior.sd.on.unobserved.thetahat <- 10
varpri <- c("wish", "jeff", "no.prior")[1]
nparallel <- 5
# prop.geno.sub <- c(0.8, 1)
prop.geno.sub <- c(0.2, 1)
pro.geno.swap.at <- 100
# pro.geno.swap.at <- 0
em.nits <- 600
R.optim.itmax <- 5
fac.specify <- c("corr", "nfac")[2]
if(fac.specify == "corr")
  corr.explained <- .8
if(fac.specify == "nfac")
  nfac <- 20
facnam <- paste("f", 1:nfac, sep = ".")

n.seed.run.subsams <- 10
var.in.name.old <- c("N", "P", "K", "nSig", "seed", "Data")
var.in.name <- c("N", "P", "nSig", "EMtol", "EMfm", "EMbic", "EMwish", "EMK", "EMKup", "seed", "Data")
var.in.name.mash <- c("N", "P", "si", "bo", "EDtol", "seed", "Data")
var.in.name.ed <- c("N", "P", "EDmeth", "EDtol", "nSig", "seed", "Data")
default.parameters <- list()
# default.parameters$impc <- list(N = 2000, P = 148, K = 148, nSig = 1,
#                                 var.in.name = c("N", "P", "K", "nSig", "seed", "Data", "sexspecific"),
#                                 var.in.name.mash = c("N", "P", "si", "bo", "seed", "Data", "sexspecific"))
default.parameters$impc <- list(N = 2000, P = 148,
                                var.in.name = var.in.name,
                                var.in.name.mash = var.in.name.mash,
                                var.in.name.ed = var.in.name.ed)
default.parameters$eqtl <- list(N = 5000, P = 44,
                                var.in.name = var.in.name,
                                var.in.name.mash = var.in.name.mash,
                                var.in.name.ed = var.in.name.ed)






if(!"Data" %in% ls())
  Data <- "impc"
if(!exists("full.analysis"))
  full.analysis <- T
# if(full.analysis){
#   for(i in 1:length(default.parameters[[Data]]))
#     assign(names(default.parameters[[Data]][i]), default.parameters[[Data]][[i]])
# }

em.finfo <- paste("final_scaleUVpre_", scale.uv.pre,
                  "_initR_", init.R.at, "_swapat_", pro.geno.swap.at, sep = "")
# em.finfo <- paste("revision_scaleUVpre_", scale.uv.pre,
#                   "_initR_", init.R.at, "_swapat_", pro.geno.swap.at, sep = "")

emdir <- switch(date.data.arrived, '2018-03-26' = paste(base.dir, "/RData_files/em_output", sep = ""),
                '2019-04-23'  = paste(base.dir, "/RData_files/em_output_data_arrived_2019-04-23", sep = ""))
em.parallel.proc.dir <- paste(emdir, "/objects_from_parallel_proc", sep = "")
for(dirc in c(emdir, em.parallel.proc.dir))
  dir.create(dirc, showWarnings = F)
em.curated.results.file <- paste(emdir, "/curated_em_output_", em.finfo, ".RData", sep = "")
em.resimpWithSignificanceFile <- paste(emdir, "/", estimation.meth, "_resimp_with_sig_", em.finfo, ".RData", sep = "")
em.resimp.with.sig.thresholds.file <- paste(emdir, "/", estimation.meth, "_resimp_with_sig_thresholds_", em.finfo, ".RData", sep = "")
if(estimation.meth == "map"){
  curated.results.file <- em.curated.results.file
  resimpWithSignificanceFile <- em.resimpWithSignificanceFile
  resimp.with.sig.thresholds.file <- em.resimp.with.sig.thresholds.file
}

paper.plot.dir <- paste(base.dir, "/plots/impc_mv_paper_plots", sep = "")
fig.obj.dir <- paste(base.dir, "/RData_files/fig_tab_related_objects", sep = "")
procord.file <- paste(base.dir, "/RData_files/processed_data/procord.RData", sep = "")

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

#######################################
#Import centre details 
cenmap = read.csv(paste(raw.data.in.dir, "/centreIds.csv", sep = ""))
colnames(cenmap) = c("cen", "nam")
cenForPerm = c("Gmc", "H", "Ics", "J", "Tcp", "Rbrc", "Wtsi", "Bcm", "Ucd")
cenmap$forPerm = F
cenmap[cenmap$nam %in% cenForPerm, "forPerm"] = T
cenun = cenmap[cenmap$forPerm, "cen"]
cennamun = cenmap[cenmap$forPerm, "nam"]
cen.omit <- "15"


################################################
#Load various data
load(file = paste(processed.data.dir, "/impc_analysis_plan.RData", sep = "")) #pout is object containing phenotype descriptions
pout <- pout[pout$analysis.current %in% c("numerical"), ] #Only working with quantitative data currently

################################################
#Define various functions
#Function for estimating heterogeneity
calc.i2 = function(Y, S){
  if(!all(is.na(Y) == is.na(S)))
    stop("Missing values do not match across Y and S")
  dfv = rowSums(is.finite(Y)) - 1
  W = 1 / S^2
  muhat = rowSums(Y * W, na.rm = T) / rowSums(W, na.rm = T)
  Q = rowSums(W * (Y - muhat)^2, na.rm = T)
  I = (Q - dfv) / Q
  I[I < 0] = 0
  I[dfv < 1] = NA
  Q[dfv < 1] = NA
  return(list(Q = Q, I = I))
}


#Function for estimating FDR from a concordance table
fdr.est.tab <- function(tab){
  require(Hmisc)
  namv <- seq(-1, 1, by = 1)
  tab <- tab[match(namv, rownames(tab)), match(namv, colnames(tab))]
  dimnames(tab) <- list(namv, namv)
  tab[is.na(tab)] <- 0
  x <- tab["1", "-1"] + tab["-1", "1"] 
  n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
  ci <- binconf(x = x, n = n)
  ci[ci > 2 / 3] <- 2 / 3 # 2/3 is maximum q compatible with model
  fdr.est <- (4 - sqrt(16 - 24 * ci)) / (2 * 3)
  return(list(est = fdr.est[1], ci = fdr.est[2:3], n = n, x = x))
}


#Function for estimating Fsr from a concordance table
fdr.est.tab <- function(tab){
  require(Hmisc)
  namv <- seq(-1, 1, by = 1)
  tab <- tab[match(namv, rownames(tab)), match(namv, colnames(tab))]
  dimnames(tab) <- list(namv, namv)
  tab[is.na(tab)] <- 0
  x <- tab["1", "-1"] + tab["-1", "1"] 
  n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
  ci <- binconf(x = x, n = n)
  ci[ci > 2 / 3] <- 2 / 3 # 2/3 is maximum q compatible with model
  fsr.est <- .5 * (1 - sqrt(1 - 2 * ci))
  return(list(est = fsr.est[1], ci = fsr.est[2:3], n = n, x = x))
}


# fdr.est.tab(tab)
# fsr.est.tab(tab)
# 
# fdr.est.tab <- function(tab){
#   require(Hmisc)
#   namv <- seq(-1, 1, by = 1)
#   tab <- tab[match(namv, rownames(tab)), match(namv, colnames(tab))]
#   dimnames(tab) <- list(namv, namv)
#   tab[is.na(tab)] <- 0
#   x <- tab["1", "-1"] + tab["-1", "1"]
#   n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
#   if(n == 0){
#     warning("Table must have a non-zero entry in at least one corner for FDR estimation; returning FDR = NA")
#     return(list(est = NA, ci = NA, n = n, x = x))
#   }
#   k <- (sum(rowSums(tab[c("-1", "1"), , drop = F])) * sum(colSums(tab[, c("-1", "1"), drop = F]))) / (sum(tab) * n)
#   ci <- binconf(x = x, n = n)
#   ci[ci > 2 / 3] <- 2 / 3 # 2/3 is maximum q compatible with model
#   fdr.est <- (2 - sqrt(4 - 2 * ci * (4 - k))) / (4 - k)
#   return(list(est = fdr.est[1], ci = fdr.est[2:3], n = n, x = x))
# }
# 
# resl <- out.lfsr$resl
# par(mfrow = c(2, 2))
# nptseq <- c(2, 5, 10, 20) * 100
# npt <- nptseq[1]#for(npt in nptseq){
#   # npt <- 200
#   methc <- c("eb", "mash")[2]
#   tabin <- switch(methc, mash = resl$impc_mash_nSig_1$ref.lines$tabl$post,
#             eb = resl$impc_em.fit_nSig_1$ref.lines$tabl$post)
#   tab <- tabin
#   # tab <- (tabin + tabin[3:1, 3:1])/2
#   # fdr.est.tab <- function(tab){
#     require(Hmisc)
#     namv <- as.character(seq(-1, 1, by = 1))
#     tab <- tab[match(namv, rownames(tab)), match(namv, colnames(tab))]
#     dimnames(tab) <- list(namv, namv)
#     tab[is.na(tab)] <- 0
#     x <- tab["1", "-1"] + tab["-1", "1"]
#     n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
#     if(n == 0){
#       warning("Table must have a non-zero entry in at least one corner for FDR estimation; returning FDR = NA")
#       return(list(est = NA, ci = NA, n = n, x = x))
#     }
# 
#     p0 <- .75
#     k <- (sum(rowSums(tab[c("-1", "1"), , drop = F])) * sum(colSums(tab[, c("-1", "1"), drop = F]))) / (sum(tab) * n) / p0
#     vecsgn <- c("-1", "1")
#     cmat <- matrix(NA, 2, 2, dimnames = list(vecsgn, vecsgn))
#     fsrseq <- fdrseq <- 10^seq(-3, 0, length.out = npt + 1)#seq(0, .25, length.out = npt + 1)
#     fsrseq <- fdrseq <- seq(0, .25, length.out = npt + 1)
#     # fsrseq <- seq(0, .15, length.out = npt + 1)
#     loglikmat <- 0
#     # a="-1"
#     # b="1"
#     # prob.eq <- outer(fdrseq, fsrseq, function(fdr, fsr) 0.25 * k * fdr^2 + 0.5 * (1 - k * fdr^2) / (1 - fdr)^2 * (fsr^2 + (1 - fsr - fdr)^2))
#     # prob.neq <- outer(fdrseq, fsrseq, function(fdr, fsr) 0.25 * k * fdr^2 + 0.5 * (1 - k * fdr^2) / (1 - fdr)^2 * (fsr^2 + (1 - fsr - fdr)^2))
#     probl <- list()
#     for(a in vecsgn){
#       for(b in vecsgn){
#         if(tab[a, b] > 0){
#           # if(a == "-1" & b == "-1")
#           #   probc <- outer(fdrseq, fsrseq, function(fdr, fsr) 0.25 * k * fdr^2 + (1 - k * fdr^2) / (1 - fdr)^2 * (fsr^2))
#           # if(a == "1" & b == "1")
#           #   probc <- outer(fdrseq, fsrseq, function(fdr, fsr) 0.25 * k * fdr^2 + (1 - k * fdr^2) / (1 - fdr)^2 * (1 - fsr - fdr)^2)
#           if(a == b)
#             probc <- outer(fdrseq, fsrseq, function(fdr, fsr) 0.25 * k * fdr^2 + 0.5 * (1 - k * fdr^2) / (1 - fdr)^2 * (fsr^2 + (1 - fsr - fdr)^2))
#           if(a != b)
#             probc <- outer(fdrseq, fsrseq, function(fdr, fsr) 0.25 * k * fdr^2 + 0.5 * (1 - k * fdr^2) / (1 - fdr)^2 * 2 * fsr * (1 - fsr - fdr))
#           probc[outer(fdrseq, fsrseq, function(fdr, fsr) return(fdr + fsr > 1))] <- NA
#           # probc[outer(fdrseq, fsrseq, function(fdr, fsr) return(fdr + fsr > 1))] <- NA
#           # probc[outer(fdrseq, fsrseq, function(fdr, fsr) return(fsr > fdr))] <- NA
#           probl[[paste(a, b)]] <- probc
#           # print(c(a, b))
#           # print(round(probc, 2))
#           # probc[probc > 1] <- 0
#           # cmat[a, b] <- k - sum(rowSums(tab[a, , drop = F])) * sum(colSums(tab[, b, drop = F])) / (tab[a, b] * n)
#           # loglikseq <- loglikseq + tab[a, b] * log(pmax(0, .25 * (1 + fdrseq^2 * cmat[a, b])))
#           loglikmat <- loglikmat + tab[a, b] * log(probc)
#         }
#       }
#     }
#     dimnames(loglikmat) <- list(fdrseq, fsrseq)
#     pmat <- 0
#     for(i in 1:length(probl))
#       pmat <- pmat + probl[[i]]
#     all((pmat[is.finite(pmat)] - 1) < 1e-10)
#     # pmat
#     atv <- 1:(npt + 1) 
#     graphics.off()
#     image(x = atv, y = atv, z = exp(loglikmat), col = rain, xaxt = "n", yaxt = "n", xlab = "Fdr", ylab = "Fsr")
#     axis(1, at = atv, labels = round(fdrseq, 3), las = 2)
#     axis(2, at = atv, labels = round(fsrseq, 3), las = 2)
#     tab
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     max(loglikmat, na.rm = T)
#     # loglikmat
#     maxind <- which.max(c(loglikmat))
# 
#     best.fdr <- fdrseq[row(loglikmat)[maxind]]
#     best.fsr <- fsrseq[col(loglikmat)[maxind]]
#     print(best.fdr)
#     print(best.fsr)
#     loglikmat[maxind]
#     sapply(probl, function(M) M[maxind])
#     tab
# # }
# 
# 
#   sumseq <- outer(fdrseq, fsrseq, function(fdr, fsr) return(fdr + fsr))
#   # image(sumseq, col = rain)
#   sumun <- sort(unique(c(sumseq)))
#   liksum <- sapply(sumun, function(sumc) mean(loglikmat[which(sumseq == sumc)], na.rm = T))
#   liksum <- sapply(sumun, function(sumc) sum(loglikmat[which(sumseq == sumc)], na.rm = T))
#   plot(sumun, liksum)
#   post.dist <- exp(loglikmat - max(loglikmat, na.rm = T))
#   post.dist[is.na(post.dist)] <- 0
#   post.dist <- post.dist / sum(post.dist)
#   post.dist.sum <- sapply(sumun, function(sumc) sum(post.dist[which(sumseq == sumc)], na.rm = T))
#   plot(sumun, post.dist.sum)
#   sum(post.dist.sum)
#   postmn <- sapply(sumun, function(sumc) mean(loglikmat[which(sumseq == sumc)], na.rm = T))
#   plot(sumun, postmn, ty = "l")
#   fdr.prof <- exp(loglikmat[, as.character(best.fsr)])
#   fdr.prof <- fdr.prof / sum(fdr.prof)
#   fsr.prof <- exp(loglikmat[as.character(best.fdr), ])
#   fsr.prof <- fsr.prof / sum(fsr.prof)
#   plot(fdrseq, fdr.prof, ylim = c(0, max(fdr.prof)), type = "l")
#   plot(fsrseq, fsr.prof, ylim = c(0, max(fsr.prof)), type = "l")
#   dimnames(pmat) <- list(fdrseq, fsrseq)
#   probl <- lapply(probl, function(M){ dimnames(M) <- list(fdrseq, fsrseq); M})
#   fdrlook <- best.fdr
#   fsrlook <- best.fsr
#   sapply(probl, function(M) M[as.character(fdrlook), as.character(fsrlook)])
# 
#   sqrt(x / n * 4 / k)
# 

#   
#   # probc
#   postseq <- exp(loglikseq - max(loglikseq, na.rm = T))
#   # postseq[is.na(postseq)] <- 0
#   postseq <- postseq / sum(postseq)
#   ci.al <- .05
#   plot(postseq)
#   ci <- fdrseq[c(findInterval(ci.al / 2, cumsum(postseq)), findInterval(1 - ci.al / 2, cumsum(postseq)))]
#   # ci <- binconf(x = x, n = n)
#   # ci[ci > 2 / 3] <- 2 / 3 # 2/3 is maximum q compatible with model
#   fdr.est <- postseq[findInterval(.5, cumsum(postseq))]
#   return(list(est = fdr.est[1], ci = fdr.est[2:3], n = n, x = x))
# }

# 
# fdr.est.tab <- function(tab){
#   require(Hmisc)
#   namv <- seq(-1, 1, by = 1)
#   tab <- tab[match(namv, rownames(tab)), match(namv, colnames(tab))]
#   dimnames(tab) <- list(namv, namv)
#   tab[is.na(tab)] <- 0
#   x <- tab["1", "-1"] + tab["-1", "1"]
#   n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
#   if(n == 0){
#     warning("Table must have a non-zero entry in at least one corner for FDR estimation; returning FDR = NA")
#     return(list(est = NA, ci = NA, n = n, x = x))
#   }
#   k <- (sum(rowSums(tab[c("-1", "1"), , drop = F])) * sum(colSums(tab[, c("-1", "1"), drop = F]))) / (sum(tab) * n)
#   vecsgn <- c("-1", "1")
#   cmat <- matrix(NA, 2, 2, dimnames = list(vecsgn, vecsgn))
#   npt <- 1000
#   fdrseq <- seq(0, 1, length.out = npt)
#   priseq <- eq()
#   loglikseq <- 0
#   a="-1"
#   b="1"
#   for(a in vecsgn){
#     for(b in vecsgn){
#       if(tab[a, b] > 0){
#         cmat[a, b] <- k - sum(rowSums(tab[a, , drop = F])) * sum(colSums(tab[, b, drop = F])) / (tab[a, b] * n)
#         loglikseq <- loglikseq + tab[a, b] * log(pmax(0, .25 * (1 + fdrseq^2 * cmat[a, b])))
#       }
#     }
#   }
#   postseq <- exp(loglikseq - max(loglikseq, na.rm = T))
#   # postseq[is.na(postseq)] <- 0
#   postseq <- postseq / sum(postseq)
#   ci.al <- .05
#   plot(postseq)
#   ci <- fdrseq[c(findInterval(ci.al / 2, cumsum(postseq)), findInterval(1 - ci.al / 2, cumsum(postseq)))]
#   # ci <- binconf(x = x, n = n)
#   # ci[ci > 2 / 3] <- 2 / 3 # 2/3 is maximum q compatible with model
#   fdr.est <- postseq[findInterval(.5, cumsum(postseq))]
#   return(list(est = fdr.est[1], ci = fdr.est[2:3], n = n, x = x))
# }
# 
# ciseq <- seq(0, 1, len = 1000)
# matplot(ciseq, cbind((4 - sqrt(16 - 24 * ciseq)) / (2 * 3), (4 + sqrt(16 - 24 * ciseq)) / (2 * 3)), ty = "l")
# 
# #Function for estimating FDR from a concordance table
# fdr.est.tab <- function(tab){
#   require(Hmisc)
#   x <- tab["1", "-1"] + tab["-1", "1"] 
#   n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
#   ci <- binconf(x = x, n = n)
#   fdr.est <- (4 - sqrt(16 - 24 * ci)) / (2 * 3)
#   return(list(est = fdr.est[1], ci = fdr.est[2:3]))
# }
# 

# Pallette for plotting colours
rain = rainbow(n = 1000, start = 0, end = .66)[1000:1]

####################################################
#Function to check memory usage
whos <- function(){ 
  objv <- eval(quote(sapply(ls(), function(x){object.size(get(x))})), 
               envir = parent.frame())
  print(objv)
  datout <- data.frame(obj = names(objv), size = objv / 1e9);
  datout <- datout[order(datout$size), ]
  datout <- rbind(datout, data.frame(obj = "total", size = sum(datout$size)))
  rownames(datout) <- NULL
  datout$size <- round(datout$size, 2)
  return(datout)
}


estimate_null_correlation_simple <- function (data, z_thresh = 2, est_cor = TRUE){
  z = data$Bhat/data$Shat
  max_absz = apply(abs(z), 1, max)
  nullish = which(max_absz < z_thresh)
  if (length(nullish) < ncol(z)) {
    stop("not enough null data to estimate null correlation")
  }
  nullish_z = z[nullish, ]
  Vhat = cov(nullish_z)
  if (est_cor) {
    Vhat = cov2cor(Vhat)
  }
  return(Vhat)
}


# #####################################################################
# #Parameters for full MCMC on effect estimates approach
# full.data = T
# centre.specific.R = F
# run.mcmc <- F
# save.th.from.mcmc <- F
# new.run <- T
# if(new.run){
#   mv.nits = 250000
#   thin = 100
#   mv.burnin = 50000 #will collect all iterations anyway
#   mcmc.finfo = paste("saveTh_", save.th.from.mcmc, "_nIts_", mv.nits, "_thin_", thin, sep = "")
#   keepinds <- (ceiling(mv.burnin / thin) + 1):floor(mv.nits / thin)
# } else {
#   mv.nits = 25000
#   mv.burnin = 5000 #will collect all iterations anyway
#   mcmc.finfo = paste("saveTh_", save.th.from.mcmc, "_nIts_", mv.nits, sep = "")
#   keepinds <- (mv.burnin + 1):mv.nits
# }
# nkeep <- length(keepinds)
# prop.line.update <- .02
# # big.sd = 1e4 #for synthetic large SE for /hat{beta} at unmeasured phenotypes 
# prior.sd.on.mu <- 10
# mcmcdir <- paste(base.dir, "/RData_files/gibbs_output", sep = "")
# mcmcfile <- paste(mcmcdir, "/", mcmc.finfo, ".RData", sep = "")
# max.nsam = 800
# ncen.sub = 7
# 
# dir.create(mcmcdir, showWarnings = F)
# mcmc.curated.results.file <- paste(mcmcdir, "/curated_mcmc_output_", mcmc.finfo, "_nKeep_", nkeep, ".RData", sep = "")
# mcmc.resimpWithSignificanceFile <- paste(mcmcdir, "/", estimation.meth, "_resimp_with_sig_", mcmc.finfo, ".RData", sep = "") 
# mcmc.resimp.with.sig.thresholds.file <- paste(mcmcdir, "/", estimation.meth, "_resimp_with_sig_thresholds_", mcmc.finfo, ".RData", sep = "") 
# if(estimation.meth == "mcmc"){
#   curated.results.file <- mcmc.curated.results.file
#   resimpWithSignificanceFile <- mcmc.resimpWithSignificanceFile
#   resimp.with.sig.thresholds.file <- mcmc.resimp.with.sig.thresholds.file
# }
# 
# dirnam = "basic.lm"
# out.dir = paste(base.dir, "/RData_files/mcmc_output/", dirnam, sep = "")

# calc.i2 = function(Y, S){
#   if(!all(is.na(Y) == is.na(S)))
#     stop("Missing values do not match across Y and S")
#   dfv = rowSums(is.finite(Y)) - 1
#   table(dfv)
#   W = 1 / S^2
#   muhat = rowSums(Y * W) / rowSums(W)
#   Q = rowSums(W * (Y - muhat)^2)
#   I = (Q - dfv) / Q
#   I[I < 0] = 0
#   range(I)
#   return(list(Q = Q, I = I))
# }





# #Methods comparison
# scale.data <- c(T, F)[1]
# old.R.sampling <- c(T, F)[2]
# perform.uv.qc <- c(T, F)[1]
# for.hugh <- as.data.frame(t(sapply(omitl[sapply(omitl, function(x) x$omit)], unlist))[, 1:2])
# write.csv(for.hugh, file = "X:/projects/impc_mv_analysis/data_out/omitted_cen_proc_combinations.csv")
# em.finfo <- paste("scaleUVpre_", scale.uv.pre, "_initR_", init.R.at,
#                   "_missingPriorSd_", prior.sd.on.unobserved.thetahat, sep = "")
# em.finfo <- paste("uvAnalysisScaling_", uv.scaling, "_scaleUVpre_", scale.uv.pre,
#                   "_initR_", init.R.at, "_swapat_", pro.geno.swap.at, sep = "")
# em.finfo <- paste("scaleUVpre_", scale.uv.pre, "_initR_", init.R.at, "_swapat_", pro.geno.swap.at, sep = "")
#sigth = c(5e-3, 1e-3, 5e-4, 1e-4, 10^c(seq(-5, -10, by = -1), -20), 0)
#centre.specific.hyperpars = F
#hetcol = gray(.35)
#homcol = gray(.8)
#cenmap <- rbind(cenmap, data.frame(cen = 15, nam = "15", forPerm = F))
#cenun = cenun[cenun != "15"]#Only one line at this centre, and don't have name of cen


# proc.omit = "Challenge Whole Body Plethysmography"#c("FACS")
# cor.meth = c("sp", "p")[1]
# scaling = c("residsd", "none")[2]
# cenc = "all"



