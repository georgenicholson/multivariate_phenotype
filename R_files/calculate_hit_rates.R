rm(list = ls())
Data <- "impc"
xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
load(file = file.runtab)
load(file = file.resl.comp)
load(file = file.compl)
load(file = uv.results.Y.S)
load(file = file.resl.comp.fac)
for(i in 1:length(default.parameters[[Data]]))
  assign(names(default.parameters[[Data]][i]), default.parameters[[Data]][[i]])
# table(linemap$line.type)
# linemap.for.err.rate <- linemap
# linemap.for.err.rate$line.type <- ifelse(linemap$line.type == "trueMut", "trueMutTes", "negConTes")
resl.err.rates.comb <- c(resl.comp[grepl(Data, names(resl.comp)) | names(resl.comp) == "uv"],
                         list(varimax = resl.comp.fac[[mv.meth.nam.use]]))
split.use <- 1
resl.err.rates.one.split <- c(list(uv = resl.comp$uv), lapply(compl[grepl(Data, names(compl))], 
                                   function(x) list(mn = x$mnarr[, , split.use], sd = x$sdarr[, , split.use], lfsr = x$lfsrarr[, , split.use])))
for(err.data.type in c("comb", "single")){
  resl.err.rates <- switch(err.data.type, comb = resl.err.rates.comb, single = resl.err.rates.one.split)
  out.perm <- err.rate.control(resl = resl.err.rates,
                               err.rate.meth = "perm", sep.imp.thresh = F,
                               test.stat = "z", linemap = linemap, phmap = phmap, Yhat = Yhat,
                               use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                               err.thresh = .05, centre.specific.thresh = F, p.complete.null.true = 1, p.test.null.true = 1)
  print(out.perm$restab)
  out.perm.lfsr <- err.rate.control(resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],
                                    err.rate.meth = "perm", sep.imp.thresh = F,
                                    test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
                                    use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                    err.thresh = .05, centre.specific.thresh = F)
  print(out.perm.lfsr$restab)
  out.lfsr <- err.rate.control(resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],#[names(resl.comp) %in% c("eb", "ed", "mash")],
                               err.rate.meth = "lfsr", sep.imp.thresh = F,
                               test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
                               use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                               err.thresh = .05, centre.specific.thresh = F)
  print(out.lfsr$restab)
  restabl <- list(perm = out.perm$restab, perm.lfsr = out.perm.lfsr$restab, lfsr = out.lfsr$restab)
  resll <- list(perm = out.perm$resl, perm.lfsr = out.perm.lfsr$resl, lfsr = out.lfsr$resl)
  resimpl <- list(perm = out.perm$resimp, perm.lfsr = out.perm.lfsr$resimp, lfsr = out.lfsr$resimp)
  save(resll, file = paste0(global.res.dir, "/resll_", err.data.type, ".RData"), version = 2)  
  save(resimpl, file = paste0(global.res.dir, "/resimpl_", err.data.type, ".RData"), version = 2)  
  save(restabl, file = paste0(global.res.dir, "/restabl_", err.data.type, ".RData"), version = 2)  
}  

load(file = paste0(global.res.dir, "/resimpl_comb.RData"))
resimp <- resimpl$perm
resimp$eb.signsig <- resimp[, paste0(mv.meth.nam.use, ".perm.signsig")]
resimp$eb.sig <- abs(resimp$eb.signsig)
resimp$eb.t <- resimp[, paste0(mv.meth.nam.use, ".t")]
resimp$eb.th.final <- resimp[, paste0(mv.meth.nam.use, ".th.final")]
resimp <- resimp[resimp$line.type == nam.truemut, ]
resimp$uv.sig <- abs(resimp$uv.perm.signsig)
resimp$uv.signsig <- resimp$uv.perm.signsig
save(resimp, file = paste0(global.res.dir, "/resimp_comb.RData"), version = 2)  




























# resl.seed <- list()
# seed <- 1
# dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
# npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
# file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
# suppressWarnings(load(file.in))
# resl.seed$uv <- list(mn = Yhat[, ph.use], sd = smat[, ph.use])
# for(namc in names(compl)){
#   resl.seed[[namc]] <- list(mn = compl[[namc]]$mnarr[, ph.use, seed], sd = compl[[namc]]$sdarr[, ph.use, seed], 
#                             lfsr = compl[[namc]]$lfsrarr[, ph.use, seed])
# }
# out.perm.seed <- err.rate.control(resl = resl.seed[!names(resl.seed) %in% c("uv.ss")],
#                                   err.rate.meth = "perm", sep.imp.thresh = F,
#                                   test.stat = "z", linemap = linemap.for.err.rate, phmap = phmap, Yhat = Yhat,
#                                   use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                                   err.thresh = .05, centre.specific.thresh = F)
# out.perm.seed$restab
# out.perm.lfsr.seed <- err.rate.control(resl = resl.seed[!names(resl.seed) %in% c("uv", "uv.ss")],
#                                        err.rate.meth = "perm", sep.imp.thresh = F,
#                                        test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                                        use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                                        err.thresh = .05, centre.specific.thresh = F)
# out.perm.lfsr.seed$restab
# out.lfsr.seed <- err.rate.control(resl = resl.seed[!names(resl.seed) %in% c("uv", "uv.ss")],#[names(resl.comp) %in% c("eb", "ed", "mash")],
#                                   err.rate.meth = "lfsr", sep.imp.thresh = F,
#                                   test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                                   use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                                   err.thresh = .05, centre.specific.thresh = F)
# out.lfsr.seed$restab
# 

# 
# seed <- 1
# dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
# npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
# file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
# suppressWarnings(load(file.in))
# var.in.namec <- var.in.name.old
# file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.namec, sapply(var.in.namec, get), sep = "_"), collapse = "_"))
# emout.file.namc <- paste0(file.base, "_emout.RData")
# load(file = emout.file.namc)
# sub.sams.for.testing <- sams.for.lik.cross.val[1:100]
# out.post.mix <- em.update.function(Y.em = Yhat[sub.sams.for.testing, ph.use], S.em = smat[sub.sams.for.testing, ph.use],
#                                    Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
#                                    pimat = emout.mix$pi, meth = "post.mn")

# addtorow$pos <- list()
# addtorow$pos <- list(0)
# addtorow$command <- c("& \\multicolumn{9}{c}{test} \\\\\n")
# 
# 
# tabout <- print(xtable(restaball.out, label = "tab:hitrates", align = rep("r", ncol(restaball.out) + 1), add.to.row = addtorow,
#                        caption = "Hit rates and error rate comparison across methods"),
#                 caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
# 
# 
# tabout <- print(xtable(restaball.out, label = "tab:hitrates", add.to.row = addtorow))
# 
# ,
#                        caption = "Hit rates and error rate comparison across methods"),
#                 caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
# 
# 
# for(j in 1:nrow(colkeepmap))
#   tabout <- gsub(colkeepmap$keep[j], colkeepmap$keep.as[j], tabout)
# tabout
# 
# 
# Grade3 <- c("A","B","B","A","B","C","C","D","A","B",
#             "C","C","C","D","B","B","D","C","C","D")
# Grade6 <- c("A","A","A","B","B","B","B","B","C","C",
#             "A","C","C","C","D","D","D","D","D","D")
# Cohort <- table(Grade3, Grade6)
# Cohort
# 
# addtorow <- list()
# addtorow$pos <- list(0)
# addtorow$command <- c("& \\multicolumn{4}{c}{Grade 6} \\\\\n")
# print(xtable(Cohort), add.to.row = addtorow, include.colnames = FALSE)
# addtorow <- list()
# # #assign a position argument to addtorow
# # #rws are the row indexes for the row to be colored, 
# # #0 is the row index for longtable argument
# addtorow$pos <- list(1)
# # addtorow$pos <- list()
# # addtorow$pos <- list(0)
# addtorow$command <- c("&& & \\multicolumn{2}{c}{test}&\\multicolumn{2}{c}{test2} \\\\\n")
# 
# #   rws, #positions for first commands(highlighting rows)
# #   0    #position for second command (longtable argument)
# # ))
# 
# 
# #############################################################################
# #Compare hit rates on split by split


# grep("\\begin\\{tabular\\}\\{rrrrrr\\}\n", tabout.test)
# tabout.test.with.header <- gsub("\\\\begin\\{tabular\\}\\{rrrrrr\\}\n", 
#               "\\\\begin\\{tabular\\}\\{rrrrrr\\}\n & &\\\\multicolumn\\{2\\}\\{c\\}\\{Hit rate \\(\\\\%\\) 
#               when data\\}&\\\\multicolumn\\{2\\}\\{c\\}\\{Estimated FDR \\(\\\\%\\)\\} \\\\\\\\\n", tabout.test)
# nSig <- 1
# K <- 74
# si <- 1
# bo <- 1
# EDtol <- 1e-5
# EDmeth <- "justED"

#############################################################################
#Compare hit rates on combined data
# seed <- 4
# dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
# npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
# file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
# suppressWarnings(load(file.in))
# var.in.namec <- var.in.name.old
# file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.namec, sapply(var.in.namec, get), sep = "_"), collapse = "_"))
# emout.file.namc <- paste0(file.base, "_emout.RData")
# load(file = emout.file.namc)
# sub.sams.for.testing <- sams.for.lik.cross.val[1:100]
# sams.for.testing <- dimnames(compl$impc_eb_1$mnarr)[[1]][rowSums(!is.na(compl$impc_eb_1$mnarr[, ph.use, seed])) > 0]
# 
# cv.sams.curr <- intersect(linemap$geno[linemap$line.type == "trueMut"], names(compl$impc_eb_1$llmat[, seed])[!is.na(compl$impc_eb_1$llmat[, seed])])
# 
# out.post.mix.eb.obj <- em.update.function(Y.em = Y.em[cv.sams.use, ph.use], S.em = S.em[cv.sams.use, ph.use],
#                                           Sigl = lapply(emout.mix$Sigl, function(M) M[ph.use, ph.use]), R = emout.mix$R[ph.use, ph.use],
#                                           omegaseq = emout.mix$omegaseq, prior.in.obj = F,
#                                           pimat = em.pi, meth = "just.obj")

# resl = resl.seed[!names(resl.seed) %in% c("uv.ss")] 
# err.rate.meth = "perm"; sep.imp.thresh = F;
# test.stat = "z"; linemap = linemap.for.err.rate; phmap = phmap; Yhat = Yhat;
# use.upper.fp.est = F; control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1];
# err.thresh = .05; centre.specific.thresh = F

# par(mfrow = c(2, 1))
# hist(pvall[[8]])
# hist(pvall[[2]])
# pvall <- lapply(resl.comp, function(x){pvout <- x$mn; pvout[] <- pnorm(-abs(x$mn / x$sd)) * 2; pvout})
# qvall <- lapply(pvall, function(M){ qvout <- M; qvout[] <- p.adjust(c(M), method = "BY")})
# sapply(qvall, function(M) mean(M < .05, na.rm = T))
