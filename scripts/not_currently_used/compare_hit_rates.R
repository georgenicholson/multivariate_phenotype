rm(list = ls())
Data <- "impc"
xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
for(i in 1:length(default.parameters[[Data]]))
  assign(names(default.parameters[[Data]][i]), default.parameters[[Data]][[i]])

recalculate.hitrates <- F
if(recalculate.hitrates){
  load(file = file.resl.comp)
  load(file = file.compl)
  load(file = uv.results.Y.S)
  linemap.for.err.rate <- linemap
  linemap.for.err.rate$line.type <- ifelse(linemap$line.type == "trueMut", "trueMutTes", "negConTes")
  # linemap <- linemap[linemap$geno %in% rownames(resl.comp[[1]]$mn), ]
  resl.err.rates.comb <- resl.comp[grepl(Data, names(resl.comp)) | names(resl.comp) == "uv"]
  split.use <- 1
  resl.err.rates.one.split <- c(list(uv = resl.comp$uv), lapply(compl[grepl(Data, names(compl)) | names(compl) == "uv"], 
                                     function(x) list(mn = x$mnarr[, , split.use], sd = x$sdarr[, , split.use], lfsr = x$lfsrarr[, , split.use])))
  for(err.data.type in c("comb", "single")){
    resl.err.rates <- switch(err.data.type, comb = resl.err.rates.comb, single = resl.err.rates.one.split)
    resl = resl.err.rates
    err.rate.meth = "perm"; sep.imp.thresh = F;
    test.stat = "z"; linemap = linemap.for.err.rate; phmap = phmap; Yhat = Yhat;
    use.upper.fp.est = F; control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1];
    err.thresh = .05; centre.specific.thresh = F
    p.complete.null.true = 1; p.test.null.true = 1
    out.perm <- err.rate.control(resl = resl.err.rates,
                                 err.rate.meth = "perm", sep.imp.thresh = F,
                                 test.stat = "z", linemap = linemap.for.err.rate, phmap = phmap, Yhat = Yhat,
                                 use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                 err.thresh = .05, centre.specific.thresh = F, p.complete.null.true = 1, p.test.null.true = 1)
    print(out.perm$restab)
    out.perm.lfsr <- err.rate.control(resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],
                                      err.rate.meth = "perm", sep.imp.thresh = F,
                                      test.stat = "lfsr", linemap = linemap.for.err.rate, phmap = phmap, Yhat = Yhat,
                                      use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                      err.thresh = .05, centre.specific.thresh = F)
    print(out.perm.lfsr$restab)
    resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")]#[names(resl.comp) %in% c("eb", "ed", "mash")],
    err.rate.meth = "lfsr"; sep.imp.thresh = F;
    test.stat = "lfsr"; linemap = linemap.for.err.rate; phmap = phmap; Yhat = Yhat;
    use.upper.fp.est = F; control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1];
    err.thresh = .05; centre.specific.thresh = F
    out.lfsr <- err.rate.control(resl = resl.err.rates[!names(resl.err.rates) %in% c("uv", "uv.ss")],#[names(resl.comp) %in% c("eb", "ed", "mash")],
                                 err.rate.meth = "lfsr", sep.imp.thresh = F,
                                 test.stat = "lfsr", linemap = linemap.for.err.rate, phmap = phmap, Yhat = Yhat,
                                 use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                 err.thresh = .05, centre.specific.thresh = F)
    print(out.lfsr$restab)
    restabl <- list(perm = out.perm$restab, perm.lfsr = out.perm.lfsr$restab, lfsr = out.lfsr$restab)
    resimpl <- list(perm = out.perm$resimp, perm.lfsr = out.perm.lfsr$resimp, lfsr = out.lfsr$resimp)
    save(resimpl, file = paste0(global.res.dir, "/resimpl_", err.data.type, ".RData"), version = 2)  
    save(restabl, file = paste0(global.res.dir, "/restabl_", err.data.type, ".RData"), version = 2)  
  }  
}




load(file = uv.results.Y.S)
truemuts <- linemap$geno[linemap$line.type == "trueMut"]
negcons <- linemap$geno[linemap$line.type == "negCon"]
llmean.splitmat.true <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[truemuts, ], na.rm = T))
llmean.splitmat.null <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[negcons, ], na.rm = T))
# llmean.splitmat.zero <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.zero[truemuts, ], na.rm = T))
true.mns <- colMeans(llmean.splitmat.true, na.rm = T)
true.ses <- apply(llmean.splitmat.true, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))
null.mns <- colMeans(llmean.splitmat.null, na.rm = T)
null.ses <- apply(llmean.splitmat.null, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))


for(err.data.type in c("comb", "single")){#err.data.type <- "comb"#
  load(file = paste0(global.res.dir, "/restabl_", err.data.type, ".RData"))  
  formatfn <- function(v)  formatC(v * 100, digits = 1, format = "f")
  format.cvlik.fn <- function(v) formatC(v, digits = 1, format = "f")
  caption <- switch(err.data.type, comb = "IMPC data: hit rate and error rate comparison across methods", 
                    single = "IMPC data: hit rate and error rate comparison across methods (single CV split)")
  colkeep1 <- c("meth", "err.rate.meth", "test.stat", 
                "mvhitimp", "hit.rate.imp.ci.l", "hit.rate.imp.ci.u", 
                "mvhitnonimp", "hit.rate.nonimp.ci.l", "hit.rate.nonimp.ci.u",
                "line.fdr.est", "line.fdr.ci.l", "line.fdr.ci.u", 
                "fdr.est", "fdr.ci.l", "fdr.ci.u", "fdr.est.imp", "fdr.est.nonimp", 
                "ref.lines.post", "ref.lines.post.l", "ref.lines.post.u")
  restaball <- rbind(restabl$perm[, colkeep1], restabl$perm.lfsr[, colkeep1], restabl$lfsr[, colkeep1])
  restaball <- restaball[!(grepl("N", restaball$meth) | grepl("rand", restaball$meth)), ]
  names(restaball)
  restaball$hit.imp.ci <- paste0(formatfn(restaball$mvhitimp), " (", formatfn(restaball$hit.rate.imp.ci.l), "-", formatfn(restaball$hit.rate.imp.ci.u), ")")
  restaball$hit.nonimp.ci <- paste0(formatfn(restaball$mvhitnonimp), " (", formatfn(restaball$hit.rate.nonimp.ci.l), "-", formatfn(restaball$hit.rate.nonimp.ci.u), ")")
  restaball$line.fdr.ci <- paste0(formatfn(restaball$line.fdr.est), " (", formatfn(restaball$line.fdr.ci.l), "-", formatfn(restaball$line.fdr.ci.u), ")")
  restaball$fdr.ci <- paste0(formatfn(restaball$fdr.est), " (", formatfn(restaball$fdr.ci.l), "-", formatfn(restaball$fdr.ci.u), ")")
  restaball$fsr.ci <- paste0(formatfn(restaball$ref.lines.post), " (", formatfn(restaball$ref.lines.post.l), "-", formatfn(restaball$ref.lines.post.u), ")")
  restaball$Method <- NA
  restaball$facmod <- ifelse(grepl("_fa_", restaball$meth), "FA", "PCA")
  restaball$Method[grepl("em.fit", restaball$meth)] <- methnam
  # restaball$Method[grepl("em.fit", restaball$meth)] <- paste0(methnam, " (", restaball$facmod[grepl("em.fit", restaball$meth)], ")")
  restaball$Method[grepl("mash", restaball$meth)] <- "mash"
  restaball$Method[grepl("ed", restaball$meth)] <- "XD"
  restaball$Method[grepl("uv", restaball$meth)] <- "UV"
  restaball$mvhitimp[grepl("uv", restaball$meth)] <- NA
  restaball.out <- restaball
  restaball$Test <- paste0(ifelse(restaball$err.rate.meth == "perm", 1, 2), ifelse(restaball$test.stat == "z", "A", "B"))
  restaball$S <- sapply(strsplit(restaball$meth, spl = "_"), function(x) x[4])
  restaball$S[restaball$Method == "mash"] <- P + 10
  restaball$S[restaball$Method == "UV"] <- ""
  restaball$K <- sapply(strsplit(restaball$meth, spl = "_"), function(x) x[10])
  restaball$fm <- sapply(strsplit(restaball$meth, spl = "_"), function(x) x[6])
  restaball <- restaball[which(restaball$fm != "pca" | is.na(restaball$fm)), ]
  restaball$Control <- ifelse(restaball$err.rate.meth == "perm", "$\\mathrm{Fdr}_{\\mathrm{complete}}\\leq 5\\%$", "$lfsr \\leq 5\\%$")
  restaball$Statistic <- ifelse(restaball$test.stat == "z", "$z$", "$lfsr$")
  restaball$cvlik <- true.mns[restaball$meth]
  restaball$cvlik.l <- true.mns[restaball$meth] - 2 * true.ses[restaball$meth]
  restaball$cvlik.u <- true.mns[restaball$meth] + 2 * true.ses[restaball$meth]
  restaball$cvlik.ci <- paste0(format.cvlik.fn(restaball$cvlik), " (", format.cvlik.fn(restaball$cvlik.l), ",", format.cvlik.fn(restaball$cvlik.u), ")")
  restaball <- restaball[order(match(restaball$Method, c("UV", "XD", "mash", methnam))), ]
  rates.format <- c("mvhitimp", "mvhitnonimp", "line.fdr.est", "fdr.est", "fdr.est.imp", "fdr.est.nonimp")
  for(j in rates.format){
    if(is.numeric(restaball[, j]))
      restaball[, j] <- formatC(restaball[, j] * 100, digits = 1, format = "f")
  }
  restaball
  #############################################################################
  #Output power and type I error table
  library(xtable)
  colkeepmap <- data.frame(keep = c("meth", "err.rate.meth", "test.stat", "hit.imp.ci", "hit.nonimp.ci", "line.fdr.est", "line.fdr.ci",
                                    "fdr.est", "fdr.ci", "ref.lines.post", "fsr.ci", "fdr.est.imp", "fdr.est.nonimp", "Method", "Test", 
                                    "S", "K", "Control", "Statistic", "cvlik.ci"),
                           keep.as = c("meth", "err.rate.meth", "test.stat", "missing", "measured", 
                                       "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{complete}}$", "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{complete}}$",
                                       "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{single}}$", "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{single}}$",
                                       "$\\widehat{\\mathrm{Fsr}}_{\\mathrm{replicate}}$", "$\\widehat{\\mathrm{Fsr}}_{\\mathrm{replicate}}$",
                                       "fdr.est.imp", "fdr.est.nonimp", "Method", 
                                       "Test", "S", "K", "Control", "Statistic", "CV Log Likelihood"))
  colkeep <- c("Test", "Control", "Statistic", "Method", "S", "K", "hit.nonimp.ci", "hit.imp.ci", "line.fdr.ci",
               "fdr.ci", "fsr.ci")
  restaball.out <- restaball[, colkeep]
  colnames(restaball.out) <- colkeepmap[match(colkeep, colkeepmap$keep), "keep.as"]
  restaball.out
  
  tabout <- print(xtable(restaball.out, label = "tab:hitrates", align = rep("r", ncol(restaball.out) + 1),
                         caption = "Hit rates and error rate comparison across methods"),
                  caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
  cat(tabout, file = paste(revision.text.numbers, "/hitrates_combined.txt", sep = ""))
  
  
  hitList <- split(restaball.out[, 4:ncol(restaball.out)], f = restaball.out$Test)
  attr(hitList, "subheadings") <- paste0(names(hitList), ". Controlling ", 
                                         restaball.out$Control[match(names(hitList), restaball.out$Test)],
                                         " using ", restaball.out$Statistic[match(names(hitList), restaball.out$Test)],  " statistic ")
  xList <- xtableList(hitList, label = paste0("tab:hitrates.", err.data.type), 
                      align = gsub(" llllll", "llllll|", paste0(" ", paste0(rep("l", ncol(hitList[[1]]) + 1), collapse = ""))),
                      caption = caption)
  tabout.test <- print(xList, table.placement = "h!", colnames.format = "multiple", add.to.row = addtorow,
                                  caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
  tabout.test.with.header <- gsub("\n\\\\hline\nMe", 
                                  paste0("\n\\\\hline\n& & &\\\\multicolumn\\{2\\}\\{c\\}\\{Hit rate  in \\\\% when data\\}&",
                                  "\\\\multicolumn\\{3\\}\\{|c\\}\\{Estimated error rate in \\\\% (95\\\\% CI) \\}\\\\\\\\\n",
                                  # "\\\\cline\\{3\\-4\\}\\\\cline\\{5\\-7\\}\\\nMe"), tabout.test)
                                  "\\\nMe"), tabout.test)
  tabout.test.with.header
  cat(tabout.test.with.header, file = paste(revision.text.numbers, "/hitrates_", err.data.type, ".txt", sep = ""))

  
  colkeep.lik <- c("Method", "S", "K", "cvlik.ci")#, "hit.nonimp.ci", "hit.imp.ci")
  restaball.lik <- unique(restaball[, colkeep.lik])
  restaball.lik <- restaball.lik[restaball.lik$Method != "UV", ]
  colnames(restaball.lik) <- colkeepmap[match(colkeep.lik, colkeepmap$keep), "keep.as"]
  restaball.lik
  
  tabout.lik <- print(xtable(restaball.lik, label = "tab:cvlik_impc", align = rep("l", ncol(restaball.lik) + 1),
                         caption = "Comparison of cross-validated log likelihood across MV methods"),
                  caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
  cat(tabout.lik, file = paste(revision.text.numbers, "/cvlik_table.txt", sep = ""))
  
}


































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
