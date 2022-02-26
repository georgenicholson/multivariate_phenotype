# rm(list = ls())
data <- "impc"
source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")
# nfac <- 24
# vc.type <- c("vari", "pro")[1]
# load(file = paste0(meth.comp.output.dir, "/ebmix_results_nf_", nfac, "_vt_", vc.type, ".RData"))
file.glob.res <- paste0(global.res.dir, "/global_eb_results.RData")
load(file = file.glob.res)
load(file = uv.results.Y.S)
power.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/power_plots_temp"

# load(em.curated.results.file)
# load(file = resimp.with.sig.thresholds.file)
load(file = procord.file)
resimp <- resimp[resimp$line.type == "trueMut" & resimp$ph %in% phord, ]
genuse <- unique(resimp$geno)
phuse <- unique(resimp$ph)
true.use <- unique(resimp[resimp$line.type == "trueMut", "geno"])


Siguse <- resl.out$eb$Sig.comb
Ruse <- resl.out$eb$R
# 
# S = smat[genuse, phuse]
# Y <- Yhat[genuse, phuse]
# S[is.na(S)] <- prior.sd.on.unobserved.thetahat
# Y[is.na(Y)] <- 0
# Rinv = solve(Ruse)
# Siguse.inv = solve(Siguse)
# ebmn = ebsd = ebmn.proc = ebsd.proc = ebt.proc = matrix(NA, length(genuse), length(phuse), dimnames = list(genuse, phuse))
# procun = unique(pout[match(phuse, pout$ph), "procnam"])
# library(doParallel)
# if(!"clust" %in% ls())
#   clust = makePSOCKcluster(names = 15)
# registerDoParallel(clust)
# genlook = unique(resimp[resimp$line.type == "trueMut", "geno"])#sample(resimp[resimp$line.type == "trueMut", "geno"], 1000)
# t1 <- Sys.time()
# outl = foreach(genc = genlook) %dopar% {#genc = genlook#
#   #print(genc)
#   ph.meas = na.omit(colnames(S)[S[genc, ] != prior.sd.on.unobserved.thetahat])
#   proc.meas = unique(pout[match(ph.meas, pout$ph), "procnam"])
#   for(proc.loc in proc.meas){
#     Sc = S[genc, phuse]
#     yc <- Y[genc, ]
#     ph.loc = pout[pout$procnam == proc.loc, "ph"]
#     ph.loc = ph.loc[ph.loc %in% ph.meas]
#     Sc[ph.loc] = prior.sd.on.unobserved.thetahat
#     yc[ph.loc] <- 0
#     SRS.inv = t(t(Rinv[phuse, phuse]) / Sc[phuse]) / Sc[phuse]#t(t(Ruse) * S[genc, ]) * S[genc, ]
#     varc = solve(Siguse.inv[phuse, phuse] + SRS.inv[phuse, phuse])
#     ebsd.proc[genc, ph.loc] = sqrt(diag(varc))[ph.loc]
#     ebmn.proc[genc, ph.loc] = (varc %*% SRS.inv[phuse, phuse] %*% yc[phuse])[ph.loc, ]
#   }
#   out = list(sd = ebsd.proc[genc, ], mn = ebmn.proc[genc, ])
#   return(out)
# }
# t2 <- Sys.time()
# t2 - t1
# 
# outl[[1]]
# genc <- genlook[1]
# names(outl) = genlook
# for(genc in genlook){
#   ebmn.proc[genc, phuse] = outl[[genc]]$mn[phuse]
#   ebsd.proc[genc, phuse] = outl[[genc]]$sd[phuse]
# }
# 
# ######################################
# #Calibrate the output
# sd.scale.fac <- (ebl$sd / ebl.raw$sd)[genuse, phuse]
# ebsd.proc.scaled <- ebsd.proc[genuse, phuse] * sd.scale.fac
# ebt.proc[genuse, phuse] = (ebmn.proc / ebsd.proc.scaled)[genuse, phuse]
# resimp$ph_geno = paste(resimp$ph, resimp$geno, sep = "_")
# inds = match(paste(colnames(ebt.proc)[col(ebt.proc)], rownames(ebt.proc)[row(ebt.proc)], sep = "_"), resimp$ph_geno)
# notnas <- which(!is.na(inds))
# resimp[inds[notnas], "eb.mn.loo.proc"] = c(ebmn.proc)[notnas]
# resimp[inds[notnas], "eb.sd.loo.proc"] = c(ebsd.proc)[notnas]
# resimp[inds[notnas], "eb.t.loo.proc"] = c(ebt.proc)[notnas]


str(resl.out$eb$loocv.mn)

resimp$eb.mn.loo.proc <- resl.out$eb$loocv.mn[cbind(resimp$geno, resimp$ph)]
resimp$eb.sd.loo.proc <- resl.out$eb$loocv.sd[cbind(resimp$geno, resimp$ph)]
resimp$eb.t.loo.proc <- resimp$eb.mn.loo.proc / resimp$eb.sd.loo.proc


######################################
#Numbers for text
tab.loo.uv = table(c(sign(resimp$uv.t) * (abs(resimp$uv.t) > resimp$uv.th.final)),
      c(sign(resimp$eb.t.loo.proc) * (abs(resimp$eb.t.loo.proc) > resimp$eb.th.final)))
tab.loo.mv = table(c(sign(resimp$eb.t) * (abs(resimp$eb.t) > resimp$eb.th.final)),
                   c(sign(resimp$eb.t.loo.proc) * (abs(resimp$eb.t.loo.proc) > resimp$eb.th.final)))
loomv.uv.fdr <- fdr.est.tab(tab.loo.uv)
fdr.est.tab(tab.loo.mv)


loomv.uv.ci.numv <- formatC(unlist(loomv.uv.fdr) * 100, format = "f", digits = 1)
loomv.uv.ci <- paste(loomv.uv.ci.numv[1], "\\% (95\\% CI: ", loomv.uv.ci.numv[2], "\\% - ", loomv.uv.ci.numv[3], "\\%)", sep = "")
save.text <- c("loomv.uv.ci")
for(numc in save.text)
  write.table(eval(as.name(numc)), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)


###############################################
#Leave-one-out imputed vs UV scatterplot
fnamc <- "looEB_uv_scatter.jpg"
jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 12, 7, units = "in", res = 500)
par(mfrow = c(1, 2), oma = c(3, 3, 3, 3))
for(compc in c("uv", "eb")){
  tc <- switch(compc, uv = resimp$uv.t, eb = resimp$eb.t)
  thc <- switch(compc, uv = resimp$uv.th.final, eb = resimp$eb.th.final)
  labc <- switch(compc, uv = "UV", eb = "MV")
  loo.tc <- resimp$eb.t.loo.proc
  loo.thc <- resimp$eb.th.final
  ntab <- table(sign(tc) * (abs(tc) > thc), sign(loo.tc) * (abs(loo.tc) > loo.thc))
  ptab <- ntab / sum(!is.na(tc) & !is.na(loo.tc))
  limn <- 3
  cexlab <- 1.25
  legcex <- .8
  plot(tc / thc, loo.tc / loo.thc, pch = 1, cex = .2, 
       xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
  abline(0, 1)
  mtext(side = 2, line = 3, text = expression(paste("LOO-MV ", italic(tilde(z)))), cex = cexlab)
  mtext(side = 1, line = 3, text = substitute(paste(a, " ", italic(tilde(z))), list(a = labc)), cex = cexlab)
  abline(h = c(-1, 1), col = 2, lty = 2)
  abline(v = c(-1, 1), col = 2, lty = 2)
  atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
  for(i in 1:3){
    for(j in 1:3){
      # labc <- paste(formatC(100 * ptab[i, j], format = "f", digits = 1), "% (", prettyNum(ntab[i, j], big.mark = ","), ")", sep = "")
      labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
      legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .85), xjust = .5, cex = legcex)
    }
  }
  mtext(side = 3, line = 2, at = -limn, text = ifelse(compc == "uv", "(a)", "(b)"), cex = 1.3)
  if(compc == "uv")
    fdr.res <- fdr.est.tab(tab.loo.uv)
  if(compc == "eb")
    fdr.res <- fdr.est.tab(tab.loo.mv)
  numv <- formatC(c(fdr.res[[1]], fdr.res[[2]]) * 100, format = "f", digits = 1)
  fdrlab <- paste("FDR ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
  mtext(side = 3, text = fdrlab, line = .5, cex = 1)
}
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)


###############################################
#Plot estimated power by procedure
resimp$procnam <- phmap[match(resimp$ph, phmap$ph), "procnam"]
procun <- unique(resimp$procnam)
powtab <- data.frame(procnam = procun, uv = NA, eb = NA, loo.eb = NA, n.uv = NA, n.eb.non = NA, n.loo.eb = NA)
rownames(powtab) <- procun
resimp$loo.eb.sig <- abs(resimp$eb.t.loo.proc) > resimp$eb.th.final
for(procc in procun){
  resc <- resimp[!resimp$imputed & resimp$procnam == procc & resimp$line.type == "trueMut", ]
  powtab[procc, "x.uv"] <- length(unique(resc[which(resc$uv.perm.signsig != 0), "geno"]))
  powtab[procc, "x.eb.non"] <- length(unique(resc[which(resc$eb.perm.signsig != 0), "geno"]))
  # powtab[procc, "x.loo.eb"] <- length(unique(resc[resc$loo.eb.sig, "geno"]))
  powtab[procc, c("n.uv", "n.eb.non")] <- rep(length(unique(resc$geno)), 2)
  resc.imp <- resimp[resimp$imputed & resimp$procnam == procc & resimp$line.type == "trueMut", ]
  powtab[procc, "x.eb.imp"] <- length(unique(resc.imp[which(resc.imp$eb.perm.signsig != 0), "geno"]))
  powtab[procc, "n.eb.imp"] <- length(unique(resc.imp$geno))
  # powtab[procc, "x.loo.eb"] <- length(unique(resc.imp[which(resc.imp$loo.eb.sig), "geno"]))
  # powtab[procc, "n.loo.eb"] <- length(unique(resc.imp$geno))
  # resc.bo <- resimp[resimp$procnam == procc & resimp$line.type == "trueMut", ]
  # powtab[procc, "x.loo.eb"] <- length(unique(resc.bo[which(resc.bo$loo.eb.sig), "geno"]))
  # powtab[procc, "n.loo.eb"] <- length(unique(resc.bo$geno))
  # powtab[procc, "x.uv"] <- sum(resc[, "uv.sig"])
  # powtab[procc, "x.eb"] <- sum(resc[, "eb.sig"])
  # powtab[procc, "x.loo.eb"] <- sum(resc[, "loo.eb.sig"])
  # powtab[procc, "n"] <- nrow(resc)
}
uv.ci <- binconf(x = powtab$x.uv, n = powtab$n.uv)
powtab[, c("uv.est", "uv.l", "uv.u")] <- uv.ci
eb.non.ci <- binconf(x = powtab$x.eb.non, n = powtab$n.eb.non)
powtab[, c("eb.non.est", "eb.non.l", "eb.non.u")] <- eb.non.ci
# loo.eb.ci <- binconf(x = powtab$x.loo.eb, n = powtab$n.loo.eb)
# powtab[, c("loo.eb.est", "loo.eb.l", "loo.eb.u")] <- loo.eb.ci
eb.imp.ci <- binconf(x = powtab$x.eb.imp, n = powtab$n.eb.imp)
powtab[, c("eb.imp.est", "eb.imp.l", "eb.imp.u")] <- eb.imp.ci
powtab <- powtab[order(powtab$uv.est), ]
# ydum <- powtab[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
#                    "eb.imp.est", "eb.imp.l", "eb.imp.u", "loo.eb.est", "loo.eb.l", "loo.eb.u")]
ydum <- powtab[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
                   "eb.imp.est", "eb.imp.l", "eb.imp.u")]
nproc <- length(procun)
xpl <- 1:nproc
colv <- c(uv = "black", eb.non = "red", eb.imp = "blue", eb.loocv = "green")
colv
cexlab = 1.3
fnamc <- "power_MV_looMV_UV.jpg"
jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 12, 12, units = "in", res = 500)
par(mar = c(25, 8, 4, 4))
matplot(x = 1:length(procun), y = ydum, ty = "n", xaxt = "n", xlab = "", ylab = "", xlim = c(.5, nproc + .5), ylim = c(0, max(ydum)),
        xaxs = "i", yaxs = "i", las = 2, cex.axis = cexlab)
eps = .2
# points(x = rep(xpl, 4) + rep(c(-.5, -.16, .16, .5) * eps, each = nproc), 
#        y = c(powtab$uv.est, powtab$eb.non.est, powtab$eb.imp.est, powtab$loo.eb.est), col = rep(colv, each = nproc),
#        pch = 19)
points(x = rep(xpl, 3) + rep(c(-.5, -.16, .16) * eps, each = nproc), 
       y = c(powtab$uv.est, powtab$eb.non.est, powtab$eb.imp.est), col = rep(colv[c("uv", "eb.non", "eb.imp")], each = nproc), pch = 19)
for(i in 1:nproc){
  lines(x = rep(i - .5 * eps, 2), y = unlist(powtab[i, c("uv.l", "uv.u")]), col = colv["uv"])
  lines(x = rep(i - .16 * eps, 2), y = unlist(powtab[i, c("eb.non.l", "eb.non.u")]), col = colv["eb.non"])
  lines(x = rep(i + .16 * eps, 2), y = unlist(powtab[i, c("eb.imp.l", "eb.imp.u")]), col = colv["eb.imp"])
  # lines(x = rep(i + .5 * eps, 2), y = unlist(powtab[i, c("loo.eb.l", "loo.eb.u")]), col = colv["eb.loocv"])
}
lines(x = 1:nproc - .5 * eps, y = powtab$uv.est, col = colv["uv"])
lines(x = 1:nproc - .16 * eps, y = powtab$eb.non.est, col = colv["eb.non"])
lines(x = 1:nproc + .16 * eps, y = powtab$eb.imp.est, col = colv["eb.imp"])
# lines(x = 1:nproc + .5 * eps, y = powtab$loo.eb.est, col = colv["eb.loocv"])
abline(v = .5 + 0:nproc)
axis(side = 1, at = 1:nproc, labels = powtab$procnam, las = 2, cex.axis = cexlab)
mtext(side = 2, line = 4, text = "Proportion of KO lines with one or more perturbations identified", cex = cexlab)
# legend(x = "topleft", legend = c("UV model", "MV model (non-imputed)", "MV model (LOO-imputed)"), col = colv, lty = 1,
#        cex = cexlab)
# legend(x = "topleft", legend = c("MV model (non-imputed)", "MV model (imputed)", "MV model (LOO-imputed)", "UV model"), 
#       col = colv[c("eb.non", "eb.imp", "eb.loocv", "uv")], lty = 1, cex = cexlab)
legend(x = "topleft", legend = c("MV model (non-imputed)", "MV model (imputed)", "UV model"), 
       col = colv[c("eb.non", "eb.imp", "uv")], lty = 1, cex = cexlab)
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)




###############################################
#Plot estimated power by phenotype
phenun <- unique(resimp$ph)
powtab.ph <- data.frame(ph = phenun, uv = NA, eb = NA, loo.eb = NA, n.uv = NA, n.eb.non = NA, n.loo.eb = NA)
rownames(powtab.ph) <- phenun
resimp$loo.eb.sig <- abs(resimp$eb.t.loo.proc) > resimp$eb.th.final
for(phenc in phenun){
  resc <- resimp[!is.na(resimp$uv.t) & resimp$ph == phenc & resimp$line.type == "trueMut", ]
  powtab.ph[phenc, "x.uv"] <- length(unique(resc[which(resc$uv.perm.signsig != 0), "geno"]))
  powtab.ph[phenc, "x.eb.non"] <- length(unique(resc[which(resc$eb.perm.signsig != 0), "geno"]))
  powtab.ph[phenc, "x.loo.eb"] <- length(unique(resc[resc$loo.eb.sig, "geno"]))
  powtab.ph[phenc, c("n.uv", "n.eb.non", "n.loo.eb")] <- rep(length(unique(resc$geno)), 3)
  resc.imp <- resimp[is.na(resimp$uv.t) & resimp$ph == phenc & resimp$line.type == "trueMut", ]
  powtab.ph[phenc, "x.eb.imp"] <- length(unique(resc.imp[which(resc.imp$eb.perm.signsig != 0), "geno"]))
  powtab.ph[phenc, "n.eb.imp"] <- length(unique(resc.imp$geno))
  powtab.ph[phenc, "x.loo.eb"] <- length(unique(resc.imp[resc.imp$loo.eb.sig, "geno"]))
  powtab.ph[phenc, "n.loo.eb"] <- length(unique(resc.imp$geno))
  # resc.bo <- resimp[resimp$ph == phenc & resimp$line.type == "trueMut", ]
  # powtab.ph[phenc, "x.loo.eb"] <- length(unique(resc.bo[resc.bo$loo.eb.sig, "geno"]))
  # powtab.ph[phenc, "n.loo.eb"] <- length(unique(resc.bo$geno))
}
uv.ci <- binconf(x = powtab.ph$x.uv, n = powtab.ph$n.uv)
powtab.ph[, c("uv.est", "uv.l", "uv.u")] <- uv.ci
eb.non.ci <- binconf(x = powtab.ph$x.eb.non, n = powtab.ph$n.eb.non)
powtab.ph[, c("eb.non.est", "eb.non.l", "eb.non.u")] <- eb.non.ci
loo.eb.ci <- binconf(x = powtab.ph$x.loo.eb, n = powtab.ph$n.loo.eb)
powtab.ph[, c("loo.eb.est", "loo.eb.l", "loo.eb.u")] <- loo.eb.ci
eb.imp.ci <- binconf(x = powtab.ph$x.eb.imp, n = powtab.ph$n.eb.imp)
powtab.ph[, c("eb.imp.est", "eb.imp.l", "eb.imp.u")] <- eb.imp.ci
powtab.ph <- powtab.ph[match(phord, powtab.ph$ph), ]
powtab.ph$procnam <- pout[match(powtab.ph$ph, pout$ph), "procnam"]
# powtab.ph <- powtab.ph[order(match(powtab.ph$procnam, procord), powtab.ph$uv.est), ]
ydum <- powtab.ph[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
                   "eb.imp.est", "eb.imp.l", "eb.imp.u", "loo.eb.est", "loo.eb.l", "loo.eb.u")]
nphen <- length(phenun)
xpl <- 1:nphen
cexlab = 1.3
fnamc <- "power_MV_looMV_UV_by_phenotype.jpg"
jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 20, 12, units = "in", res = 1000)
par(mar = c(25, 8, 4, 4))
matplot(x = xpl, y = ydum, ty = "n", xaxt = "n", xlab = "", ylab = "", xlim = c(.5, nphen + .5), ylim = c(0, max(ydum)), 
        xaxs = "i", yaxs = "i", las = 2, cex.axis = cexlab)
abline(v = match(procun, powtab.ph$procnam) - .5, col = "black")
eps = .4
points(x = rep(xpl, 3) + rep(c(-.5, -.0, .5) * eps, each = nphen), 
       y = c(powtab.ph$uv.est, powtab.ph$eb.non.est, powtab.ph$eb.imp.est), col = rep(colv[c("uv", "eb.non", "eb.imp")], each = nphen),
       pch = 19)
for(i in 1:nphen){
  lines(x = rep(i - .5 * eps, 2), y = unlist(powtab.ph[i, c("uv.l", "uv.u")]), col = colv["uv"])
  lines(x = rep(i - .0 * eps, 2), y = unlist(powtab.ph[i, c("eb.non.l", "eb.non.u")]), col = colv["eb.non"])
  # lines(x = rep(i - .0 * eps, 2), y = unlist(powtab.ph[i, c("eb.non.l", "eb.non.u")]), col = colv[2])
  lines(x = rep(i + .5 * eps, 2), y = unlist(powtab.ph[i, c("eb.imp.l", "eb.imp.l")]), col = colv["eb.imp"])
}
lines(x = 1:nphen - .5 * eps, y = powtab.ph$uv.est, col = colv["uv"])
lines(x = 1:nphen - .0 * eps, y = powtab.ph$eb.non.est, col = colv["eb.non"])
lines(x = 1:nphen + .5 * eps, y = powtab.ph$eb.imp.est, col = colv["eb.imp"])
ats <- sapply(procord, function(procc) mean(which(powtab.ph$procnam == procc)))
axis(side = 1, at = ats, labels = procord, las = 2, cex.axis = cexlab)
mtext(side = 2, line = 4, text = "Proportion of KO lines annotated", cex = cexlab)
# legend(x = "topright", legend = c("UV", "MV (non-imputed)", "LOO-MV"), col = colv, lty = 1,
#        cex = cexlab)
legend(x = "topleft", legend = c("MV model (non-imputed)", "MV model (imputed)", "UV model"), col = colv[c("eb.non", "eb.imp", "uv")], lty = 1,
       cex = cexlab, lwd = 2)
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)


#























































# 
# ebt = (ebmn / ebsd)[genuse, phuse]
# inds = match(paste(colnames(ebt)[col(ebt)], rownames(ebt)[row(ebt)], sep = "_"), resimp$ph_geno)
# resimp[inds, "eb.mean.mn.loo.ph"] = c(ebmn)
# resimp[inds, "eb.mean.sd.loo.ph"] = c(ebsd)
# resimp[inds, "eb.mean.t.loo.ph"] = c(ebt)
# table(c(sign(resimp$uv.t) * (abs(resimp$uv.t) > resimp$th.uv.phcen.interaction.fdr)),
#       c(sign(resimp$eb.mean.t.loo.ph) * (abs(resimp$eb.mean.t.loo.ph) > resimp$th.eb.phcen.interaction.fdr)))






# for(loc in ph.meas){
#   Sc = S
#   Sc[genc, loc] = big.sd
#   SRS.inv = t(t(Rinv[phuse, phuse]) / Sc[genc, phuse]) / Sc[genc, phuse]#t(t(Ruse) * S[genc, ]) * S[genc, ]
#   varc = solve(Siguse.inv[phuse, phuse] + SRS.inv[phuse, phuse])
#   ebsd[genc, loc] = sqrt(diag(varc))[loc]
#   ebmn[genc, loc] = (varc %*% (Siguse.inv[phuse, phuse] %*% mu.st[phuse] + SRS.inv[phuse, phuse] %*% Y[genc, phuse]))[loc, ]
# }
# fnamc <- "power_MV_looMV_UV.jpg"
# jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 12, 12, units = "in", res = 500)
# par(mar = c(25, 8, 4, 4))
# matplot(x = 1:length(procun), y = ydum, ty = "n", xaxt = "n", xlab = "", ylab = "", xlim = c(.5, nproc + .5), 
#         xaxs = "i", las = 2, cex.axis = cexlab)
# eps = .2
# points(x = rep(xpl, 4) + rep(c(-.75, -.25, .25, .75) * eps, each = nproc), 
#        y = c(powtab.ph$uv.est, powtab$eb.non.est, powtab$eb.imp.est, powtab$loo.eb.est), col = rep(colv, each = nproc),
#        pch = 19)
# for(i in 1:nproc){
#   lines(x = rep(i - .75 * eps, 2), y = unlist(powtab[i, c("uv.l", "uv.u")]), col = colv[1])
#   lines(x = rep(i - .25 * eps, 2), y = unlist(powtab[i, c("eb.non.l", "eb.non.u")]), col = colv[2])
#   lines(x = rep(i + .25 * eps, 2), y = unlist(powtab[i, c("eb.imp.l", "eb.imp.u")]), col = colv[3])
#   lines(x = rep(i + .75 * eps, 2), y = unlist(powtab[i, c("loo.eb.l", "loo.eb.u")]), col = colv[4])
# }
# lines(x = 1:nproc - .75 * eps, y = powtab$uv.est, col = colv[1])
# lines(x = 1:nproc - .25 * eps, y = powtab$eb.non.est, col = colv[2])
# lines(x = 1:nproc + .25 * eps, y = powtab$eb.imp.est, col = colv[3])
# lines(x = 1:nproc + .75 * eps, y = powtab$loo.eb.est, col = colv[4])
# abline(v = .5 + 0:nproc)
# axis(side = 1, at = 1:nproc, labels = powtab$procnam, las = 2, cex.axis = cexlab)
# mtext(side = 2, line = 4, text = "Proportion of KO lines with one or more perturbations identified", cex = cexlab)
# legend(x = "topleft", legend = c("UV model", "MV model (non-imputed)", "MV model (imputed)", "MV model (LOO-imputed)"), col = colv, lty = 1,
#        cex = cexlab)
# dev.off()
# file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
#           to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
# 


# ###############################################
# #Plot estimated power by procedure
# resimp$procnam <- phmap[match(resimp$ph, phmap$ph), "procnam"]
# procun <- unique(resimp$procnam)
# powtab <- data.frame(procnam = procun, uv = NA, eb = NA, loo.eb = NA, n.uv = NA, n.eb.non = NA, n.loo.eb = NA)
# rownames(powtab) <- procun
# resimp$loo.eb.sig <- abs(resimp$eb.t.loo.proc) > resimp$eb.th.final
# for(procc in procun){
#   resc <- resimp[!resimp$imputed & resimp$procnam == procc & resimp$line.type == "trueMut", ]
#   powtab[procc, "x.uv"] <- length(unique(resc[which(resc$uv.perm.signsig != 0), "geno"]))
#   powtab[procc, "x.eb.non"] <- length(unique(resc[which(resc$eb.perm.signsig != 0), "geno"]))
#   powtab[procc, "x.loo.eb"] <- length(unique(resc[resc$loo.eb.sig, "geno"]))
#   powtab[procc, c("n.uv", "n.eb.non", "n.loo.eb")] <- rep(length(unique(resc$geno)), 3)
#   resc.imp <- resimp[resimp$imputed & resimp$procnam == procc & resimp$line.type == "trueMut", ]
#   powtab[procc, "x.eb.imp"] <- length(unique(resc.imp[which(resc.imp$eb.perm.signsig != 0), "geno"]))
#   powtab[procc, "n.eb.imp"] <- length(unique(resc.imp$geno))
#   # powtab[procc, "x.uv"] <- sum(resc[, "uv.sig"])
#   # powtab[procc, "x.eb"] <- sum(resc[, "eb.sig"])
#   # powtab[procc, "x.loo.eb"] <- sum(resc[, "loo.eb.sig"])
#   # powtab[procc, "n"] <- nrow(resc)
# }
# uv.ci <- binconf(x = powtab$x.uv, n = powtab$n.uv)
# powtab[, c("uv.est", "uv.l", "uv.u")] <- uv.ci
# eb.non.ci <- binconf(x = powtab$x.eb.non, n = powtab$n.eb.non)
# powtab[, c("eb.non.est", "eb.non.l", "eb.non.u")] <- eb.non.ci
# loo.eb.ci <- binconf(x = powtab$x.loo.eb, n = powtab$n.loo.eb)
# powtab[, c("loo.eb.est", "loo.eb.l", "loo.eb.u")] <- loo.eb.ci
# eb.imp.ci <- binconf(x = powtab$x.eb.imp, n = powtab$n.eb.imp)
# powtab[, c("eb.imp.est", "eb.imp.l", "eb.imp.u")] <- eb.imp.ci
# powtab <- powtab[order(powtab$uv.est), ]
# ydum <- powtab[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
#                    "eb.imp.est", "eb.imp.l", "eb.imp.u", "loo.eb.est", "loo.eb.l", "loo.eb.u")]
# nproc <- length(procun)
# xpl <- 1:nproc
# colv <- 1:3
# cexlab = 1.3
# fnamc <- "power_MV_looMV_UV.jpg"
# jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 12, 12, units = "in", res = 500)
# par(mar = c(25, 8, 4, 4))
# matplot(x = 1:length(procun), y = ydum, ty = "n", xaxt = "n", xlab = "", ylab = "", xlim = c(.5, nproc + .5), ylim = c(0, max(ydum)),
#         xaxs = "i", yaxs = "i", las = 2, cex.axis = cexlab)
# eps = .2
# points(x = rep(xpl, 3) + rep(c(-.5, -.0, .5) * eps, each = nproc), 
#        y = c(powtab$uv.est, powtab$eb.non.est, powtab$loo.eb.est), col = rep(colv, each = nproc),
#        pch = 19)
# for(i in 1:nproc){
#   lines(x = rep(i - .5 * eps, 2), y = unlist(powtab[i, c("uv.l", "uv.u")]), col = colv[1])
#   lines(x = rep(i - .0 * eps, 2), y = unlist(powtab[i, c("eb.non.l", "eb.non.u")]), col = colv[2])
#   lines(x = rep(i + .5 * eps, 2), y = unlist(powtab[i, c("loo.eb.l", "loo.eb.u")]), col = colv[3])
# }
# lines(x = 1:nproc - .5 * eps, y = powtab$uv.est, col = colv[1])
# lines(x = 1:nproc - .0 * eps, y = powtab$eb.non.est, col = colv[2])
# lines(x = 1:nproc + .5 * eps, y = powtab$loo.eb.est, col = colv[3])
# abline(v = .5 + 0:nproc)
# axis(side = 1, at = 1:nproc, labels = powtab$procnam, las = 2, cex.axis = cexlab)
# mtext(side = 2, line = 4, text = "Proportion of KO lines with one or more perturbations identified", cex = cexlab)
# # legend(x = "topleft", legend = c("UV model", "MV model (non-imputed)", "MV model (LOO-imputed)"), col = colv, lty = 1,
# #        cex = cexlab)
# legend(x = "topleft", legend = c("MV model (non-imputed)", "MV model (LOO-imputed)", "UV model"), col = colv[c(2:3, 1)], lty = 1,
#        cex = cexlab)
# dev.off()
# file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
#           to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
# 
# 
# #
#
# plot(powtab$loo.eb.est, powtab$eb.imp.est)
# abline(0, 1)
# 
# graphics.off()
# plot(resimp$eb.t, resimp$eb.t.loo.proc)
# mean(abs(resimp$eb.t.loo.proc) > abs(resimp$eb.t), na.rm = T)
# mean(abs(resimp$eb.t.loo.proc) > abs(resimp$eb.t) & abs(resimp$eb.t.loo.proc) > resimp$eb.th.final & abs(resimp$eb.t) < resimp$eb.th.final, na.rm = T)
# 

