Data <- "impc"
source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")
source(paste(R.file.dir, "/impc_mv_paper_plots_code/plot_power_comparisons_revision.R", sep = ""))
# source(paste(R.file.dir, "//impc_mv_paper_plots_code/post_eb_assess_imputation_revision.R", sep = ""))
# source(paste(R.file.dir, "/impc_plots/plot_ref_lines.R", sep = ""))

# load(file = em.curated.results.file)
# load(file = resimp.with.sig.thresholds.file)
# power.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/power_plots"
# dir.create(power.plot.dir, showWarnings = F)
# phen.un <- unique(resimp$ph)
# resimp <- resimp[resimp$line.type == nam.truemut, ]
# trueline.un <- unique(resimp[resimp$line.type == nam.truemut, "geno"])
# n.uv <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & !is.na(resimp$uv.t), "uv.sig"], na.rm = T))
# n.eb.non <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & !is.na(resimp$uv.t), "eb.sig"], na.rm = T))
# n.eb.imp <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & is.na(resimp$uv.t), "eb.sig"], na.rm = T))


###############################################
#UV/MV global line-by-line, phen-by-phen scatterplot
resimp1 <- resimp[resimp$line.type == nam.truemut, ]
# resimp1 <- resimp[resimp$line.type == nam.truemut, ]
resimp1 <- resimp1[order(resimp1$geno), ]
trueline.un <- unique(resimp1$geno)
froms <- match(trueline.un, resimp1$geno)
tos <- c(froms[2:length(froms)] - 1, nrow(resimp1))
# resimp1$uv.sig <- resimp1$uv.signsig != 0
# resimp1$eb.sig <- resimp1$eb.perm.signsig != 0
signifl <- Map(function(from, to) resimp1[from:to, c("uv.sig", "eb.sig")], from = froms, to = tos)
line.n.uv <- sapply(signifl, function(x) sum(x$uv.sig, na.rm = T))
line.n.eb.non <- sapply(signifl, function(x) sum(x$eb.sig[!is.na(x$uv.sig)], na.rm = T))
line.n.eb.imp <- sapply(signifl, function(x) sum(x$eb.sig[is.na(x$uv.sig)], na.rm = T))

wideps <- .5
heieps <- .5
graphics.off()
fnamc <- "combined_power_plot.jpg"
jpeg(file = paste(power.plot.dir, "/", fnamc, sep = ""), 12, 15, units = "in", res = 500)
layout(matrix(c(1:2, 3, 3), 2, 2, byrow = T))
par(oma = c(20, 4, 0, 0))
# par(mfrow = c(1, 2), mar = c(6, 6, 4, 4), oma = c(0, 0, 0, 0))
plty <- "non"
xc <- n.uv
yc <- switch(plty, non = n.eb.non,
             imp = n.eb.non + n.eb.imp)
lims <- range(c(xc, yc))#range(c(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp)))
repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
circmult1 <- 5
circmult2 <- 1 / 3
cexax <- 1.2
cexax2 <- 1.2
cexax3 <- 1.6
axnumcex <- 1.4
symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) * circmult1, inches = F, xlim = lims, ylim = lims, las = 1,
        xlab = "", ylab = "", cex.axis = axnumcex)
linc <- 3.5
linc2 <- 2
mtext(side = 1, line = linc, text = "UV analysis", cex = cexax)
mtext(side = 2, line = linc, text = "MV analysis", cex = cexax)
mtext(side = 3, line = linc2, text = "Annotations per phenotype", cex = cexax2)
mtext(side = 3, line = linc2, at = 0, text = "(a)", cex = cexax3)
abline(0, 1)

xc <- line.n.uv
yc <- switch(plty, non = line.n.eb.non,
             imp = line.n.eb.non + line.n.eb.imp)
lims <- range(c(xc, yc))#range(c(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp)))
repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) * circmult2, inches = F, xlim = lims, ylim = lims, las = 1,
        xlab = "", ylab = "", cex.axis = axnumcex)
abline(0, 1)
mtext(side = 1, line = linc, text = "UV analysis", cex = cexax)
mtext(side = 2, line = linc, text = "MV analysis", cex = cexax)
mtext(side = 3, line = linc2, text = "Annotations per KO gene", cex = cexax2)
mtext(side = 3, line = linc2 / 3, text = "(point area proportional to number)", cex = cexax2 * .75)
mtext(side = 3, line = linc2, at = 0, text = "(b)", cex = cexax3)
# dev.off()
# file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
#           to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)


###############################################
#Plot estimated power by procedure
resimp$procnam <- phmap[match(resimp$ph, phmap$ph), "procnam"]
procun <- unique(resimp$procnam)
powtab <- data.frame(procnam = procun, uv = NA, eb = NA, loo.eb = NA, n.uv = NA, n.eb.non = NA, n.loo.eb = NA)
rownames(powtab) <- procun
# resimp$loo.eb.sig <- abs(resimp$eb.t.loo.proc) > resimp$eb.th.final
for(procc in procun){
  resc <- resimp[!resimp$imputed & resimp$procnam == procc & resimp$line.type == nam.truemut, ]
  powtab[procc, "x.uv"] <- length(unique(resc[which(resc$uv.signsig != 0), "geno"]))
  powtab[procc, "x.eb.non"] <- length(unique(resc[which(resc$eb.signsig != 0), "geno"]))
  # powtab[procc, "x.loo.eb"] <- length(unique(resc[resc$loo.eb.sig, "geno"]))
  powtab[procc, c("n.uv", "n.eb.non")] <- rep(length(unique(resc$geno)), 2)
  resc.imp <- resimp[resimp$imputed & resimp$procnam == procc & resimp$line.type == nam.truemut, ]
  powtab[procc, "x.eb.imp"] <- length(unique(resc.imp[which(resc.imp$eb.signsig != 0), "geno"]))
  powtab[procc, "n.eb.imp"] <- length(unique(resc.imp$geno))
  # powtab[procc, "x.loo.eb"] <- length(unique(resc.imp[which(resc.imp$loo.eb.sig), "geno"]))
  # powtab[procc, "n.loo.eb"] <- length(unique(resc.imp$geno))
  # resc.bo <- resimp[resimp$procnam == procc & resimp$line.type == nam.truemut, ]
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
ydum <- powtab[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
                   "eb.imp.est", "eb.imp.l", "eb.imp.u")]
# ydum <- powtab[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
#                    "eb.imp.est", "eb.imp.l", "eb.imp.u", "loo.eb.est", "loo.eb.l", "loo.eb.u")]
nproc <- length(procun)
xpl <- 1:nproc
colv <- c(uv = "black", eb.non = "red", eb.imp = "blue", eb.loocv = "green")
colv
cexlab = 1.3
# par(mar = c(25, 8, 4, 4))
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
mtext(side = 3, line = linc2, at = 0, text = "(c)", cex = cexax3)
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)



diff.text <- ydum[, "eb.non.est"] - ydum[, "uv.est"]
names(diff.text) <- rownames(ydum)
for(rowc in rownames(ydum))
dir.save <- paste(revision.text.numbers, "/proc_pow_up", sep = "")
dir.create(dir.save, showWarnings = F)
  #Export numbers to text
for(numc in save.prop)
for(rowc in rownames(ydum))
    write.table(formatC(100 * diff.text[rowc], digits = 0, format = "f"), file = paste(dir.save, "/", gsub(" ", "_", gsub("/", "", rowc)), ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)


# 
# ####################################################
# #Reference lines concordance scatter
# fnamc <- "ref_lines_z_scatter.jpg"
# jpeg(paste(reflines.plot.dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
# par(mfrow = c(2, 2), oma = c(3, 2, 1, 1), mar = c(3, 3, 3, 3))
# namv <- c("pre", "postcomp", "postimp", "postimpcomp")
# for(namc in namv){
#   ntab <- table(matl.sig[[namc]][, 1], matl.sig[[namc]][, 2])
#   ptab <- ntab / sum(ntab)
#   limn <- 3
#   legcex <- .65
#   plot(matl[[namc]][, 1], matl[[namc]][, 2], pch = 1, cex = .2, 
#        xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
#   abline(0, 1)
#   labc1 <- switch(namc, pre = "UV replicate 1",
#                   postcomp = "MV non-imputed replicate 1",
#                   postimp = "MV imputed replicate 1",
#                   postimpcomp = "MV imputed replicate 1")
#   labc2 <- switch(namc, pre = "UV replicate 2, ",
#                   postcomp = "MV non-imputed replicate 2",
#                   postimp = "MV imputed replicate 2",
#                   postimpcomp = "MV non-imputed replicate 2")
#   mtext(side = 1, line = 2.5, text = substitute(paste(labc1, ", ", italic(tilde(z))), list(labc1 = labc1)))
#   mtext(side = 2, line = 2.5, text = substitute(paste(labc2, ", ", italic(tilde(z))), list(labc2 = labc2)))
#   abline(h = c(-1, 1), col = 2, lty = 2)
#   abline(v = c(-1, 1), col = 2, lty = 2)
#   atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
#   for(i in 1:3){
#     for(j in 1:3){
#       # labc <- paste(formatC(100 * ptab[i, j], format = "f", digits = 1), "% (", prettyNum(ntab[i, j], big.mark = ","), ")", sep = "")
#       labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
#       legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .75), xjust = .5, cex = legcex)
#     }
#   }
#   mtext(side = 3, text = paste("(", letters[match(namc, namv)], ")", sep = ""), at = -limn, line = 1, cex = 1.3)
#   numv <- formatC(c(fdr.ests[[namc]][[1]], fdr.ests[[namc]][[2]]) * 100, format = "f", digits = 1)
#   fdrlab <- paste("Discordance-implied FDR ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
#   mtext(side = 3, text = fdrlab, line = .5, cex = .75)
# }
# dev.off()
# file.copy(from = paste(reflines.plot.dir, "/", fnamc, sep = ""),
#           to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
# 
