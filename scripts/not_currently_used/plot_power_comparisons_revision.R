rm(list = ls())
Data <- "impc"
source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")
load(file = uv.results.Y.S)
# load(file = paste0(global.res.dir, "/resimpl_comb.RData"))  
load(file = paste0(global.res.dir, "/resimp_comb.RData"))
phen.un <- unique(resimp$ph)
phen.un <- phen.un[!phen.un %in% paste0("f.", 1:nfac)]
resimp <- resimp[resimp$ph %in% phen.un, ]

power.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/power_plots_temp"
dir.create(power.plot.dir, showWarnings = F)
trueline.un <- unique(resimp[resimp$line.type == nam.truemut, "geno"])
n.uv <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & !is.na(resimp$uv.t), "uv.sig"], na.rm = T))
n.eb.non <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & !is.na(resimp$uv.t), "eb.sig"], na.rm = T))
n.eb.imp <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & is.na(resimp$uv.t), "eb.sig"], na.rm = T))

###############################################
#First set of numbers for text
n.ph.measured <- length(phen.un)
n.ko.line.measured <- length(trueline.un)
prop.miss <- mean(is.na(Yhat[trueline.un, ]))
prop.sig.eb.non <- mean(resimp$eb.sig[!resimp$imputed], na.rm = T)
n.sig.eb.non <- sum(resimp$eb.sig[!resimp$imputed], na.rm = T)
prop.sig.eb.imp <- mean(resimp$eb.sig[resimp$imputed], na.rm = T)
n.sig.eb.imp <- sum(resimp$eb.sig[resimp$imputed], na.rm = T)
prop.sig.uv <- mean(resimp$uv.sig[!resimp$imputed], na.rm = T)
n.sig.uv <- sum(resimp$uv.sig[!resimp$imputed], na.rm = T)
prop.sig.eb <- ((n.sig.eb.non + n.sig.eb.imp) / nrow(resimp))
fold.increase <- (n.sig.eb.non + n.sig.eb.imp) / n.sig.uv
nonimp.fold.increase <- prop.sig.eb.non / prop.sig.uv
imp.fold.increase <- prop.sig.eb.imp / prop.sig.uv
n.sig.eb.tot <- n.sig.eb.non + n.sig.eb.imp
prop.sig.eb.tot <- mean(resimp$eb.sig, na.rm = T)
mean.s.measured.phens <- mean(smat, na.rm = T)
sd.y.measured.phens <- sd(Yhat, na.rm = T)
for(dir.save in c(revision.text.numbers)){
  save.prop <- c("prop.sig.eb.non", "prop.sig.eb.imp", "prop.sig.uv", "prop.sig.eb.tot", "prop.miss")
  for(numc in save.prop)
    write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
                                                  col.names = F, row.names = F, quote = F)
  save.num <- c("n.sig.eb.non", "n.sig.eb.imp", "n.sig.uv", "n.sig.eb.tot", "n.ko.line.measured", "n.ph.measured")
  for(numc in save.num)
    write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
                col.names = F, row.names = F, quote = F)
  save.num2 <- c("mean.s.measured.phens", "sd.y.measured.phens")
  for(numc in save.num2)
    write.table(formatC(eval(as.name(numc)), format = "f", digits = 2), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
                col.names = F, row.names = F, quote = F)
  save.fold <- c("nonimp.fold.increase", "imp.fold.increase", "fold.increase")
  for(numc in save.fold)
    write.table(formatC(eval(as.name(numc)), digits = 1, format = "f"), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
                col.names = F, row.names = F, quote = F)
}

sum(is.na(resimp$uv.th.final) & !is.na(resimp$uv.t))
resimp[(is.na(resimp$uv.sig) & !resimp$imputed), ]


###############################################
# UV/MV numbers for text
ntab <- table(resimp$uv.signsig, resimp$eb.signsig)
ptab <- table(resimp$uv.signsig, resimp$eb.signsig) / sum(!is.na(resimp$uv.signsig))
uv.z <- resimp$uv.t / resimp$uv.th.final
eb.z <- resimp$eb.t / resimp$eb.th.final
save.tab <- c("uv_mv_prop_table")
for(numc in save.tab)
  write.table(formatC(100 * ptab, format = "f", digits = 2),
              file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
uv.sign <- sign(uv.z) * (abs(uv.z) > 1)
eb.sign <- sign(eb.z) * (abs(eb.z) > 1)

tabuvmv <- table(resimp$uv.signsig, resimp$eb.signsig)
fdrci <- fdr.est.tab(tabuvmv)
fdrc <- fdrci
loomv.uv.ci.numv <- formatC(unlist(fdrc) * 100, format = "f", digits = 1)
loomv.uv.ci <- paste(loomv.uv.ci.numv[1], "\\% (95\\% CI: ", loomv.uv.ci.numv[2], "\\% - ", loomv.uv.ci.numv[3], "\\%)", sep = "")
save.text <- "uvmv.fdr.ci"
for(numc in save.text)
  write.table(loomv.uv.ci, file = paste(revision.text.numbers, "/", save.text, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


###############################################
#UV/MV global scatterplot

fnamc <- "uv_mv_scatter.jpg"
jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
limn <- 3
cexax <- 1.3
cexax2 <- 1.5
plot(uv.z, eb.z, pch = 1, cex = .2, cex.axis = cexax,
     xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
abline(0, 1)
mtext(side = 2, line = 3, text = expression(paste("MV ", italic(tilde(z)))), cex = cexax2)
mtext(side = 1, line = 3, text = expression(paste("UV ", italic(tilde(z)))), cex = cexax2)
abline(h = c(-1, 1), col = 2, lty = 2)
abline(v = c(-1, 1), col = 2, lty = 2)
atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
for(i in 1:3){
  for(j in 1:3){
    # labc <- paste(formatC(100 * ptab[i, j], format = "f", digits = 1), "% (", prettyNum(ntab[i, j], big.mark = ","), ")", sep = "")
    labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
    legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .85), xjust = .5)
  }
}
numv <- formatC(c(fdrci[[1]], fdrci[[2]]) * 100, format = "f", digits = 1)
fdrlab <- paste("Discordance-implied FDR ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
mtext(side = 3, text = fdrlab, line = .5, cex = cexax2)
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)


###############################################
#Second set of numbers for text

uv.mv.agree.n <- ntab["1", "1"] + ntab["-1", "-1"]
uv.mv.agree.p <- ptab["1", "1"] + ptab["-1", "-1"]
uv.mv.agree.p.of.uv <- uv.mv.agree.n / n.sig.uv
uv.not.mv.n <- ntab["1", "0"] + ntab["-1", "0"]
uv.not.mv.p <- ptab["1", "0"] + ptab["-1", "0"]
mv.not.uv.n <- ntab["0", "1"] + ntab["0", "-1"]
mv.not.uv.p <- ptab["0", "1"] + ptab["0", "-1"]
uv.mv.disagree.n <- ntab["-1", "1"] + ntab["1", "-1"]
uv.mv.disagree.p <- ptab["-1", "1"] + ptab["1", "-1"]
uv.mv.disagree.p.of.uv <- uv.mv.disagree.n / n.sig.uv
save.num <- c("uv.mv.agree.n", "uv.not.mv.n", "mv.not.uv.n", "uv.mv.disagree.n")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.prop <- c("uv.mv.agree.p", "uv.mv.agree.p.of.uv", "uv.not.mv.p", "mv.not.uv.p", "uv.mv.disagree.p", "uv.mv.disagree.p.of.uv")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


###############################################
#UV/MV global line-by-line, phen-by-phen scatterplot
resimp1 <- resimp[resimp$line.type == nam.truemut, ]
resimp1 <- resimp1[order(resimp1$geno), ]
trueline.un <- unique(resimp1$geno)
froms <- match(trueline.un, resimp1$geno)
tos <- c(froms[2:length(froms)] - 1, nrow(resimp1))


signifl <- Map(function(from, to) resimp1[from:to, c("uv.sig", "eb.sig")], from = froms, to = tos)
line.n.uv <- sapply(signifl, function(x) sum(x$uv.sig, na.rm = T))
line.n.eb.non <- sapply(signifl, function(x) sum(x$eb.sig[!is.na(x$uv.sig)], na.rm = T))
line.n.eb.imp <- sapply(signifl, function(x) sum(x$eb.sig[is.na(x$uv.sig)], na.rm = T))
fnamc <- "line_by_line_uv_mv_comp.jpg"
jpeg(file = paste(power.plot.dir, "/", fnamc, sep = ""), 11.5, 6, units = "in", res = 500)
par(mfrow = c(1, 2), mar = c(6, 6, 4, 4), oma = c(0, 0, 0, 0))
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
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)



#################################################
#Numbers for line-by-line comparison
phen.mv.greater.p <- mean(n.eb.non > n.uv)
phen.uv.greater.p <- mean(n.eb.non < n.uv)
phen.uv.equal.p <- mean(n.eb.non == n.uv)
phen.mv.greater.n <- sum(n.eb.non > n.uv)
phen.uv.greater.n <- sum(n.eb.non < n.uv)
phen.uv.equal.n <- sum(n.eb.non == n.uv)
phen.average.increase.n <- mean(n.eb.non - n.uv)
line.mv.greater.p <- mean(line.n.eb.non > line.n.uv)
line.uv.greater.p <- mean(line.n.eb.non < line.n.uv)
line.uv.equal.p <- mean(line.n.eb.non == line.n.uv)
line.mv.greater.n <- sum(line.n.eb.non > line.n.uv)
line.uv.greater.n <- sum(line.n.eb.non < line.n.uv)
line.uv.equal.n <- sum(line.n.eb.non == line.n.uv)
line.average.increase.n <- mean(line.n.eb.non - line.n.uv)

plot(n.uv, n.eb.non / n.uv)

save.num <- c("phen.mv.greater.n", "phen.uv.greater.n", "phen.uv.equal.n",
              "line.mv.greater.n", "line.uv.greater.n", "line.uv.equal.n")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.num1 <- c("phen.average.increase.n", "line.average.increase.n")
for(numc in save.num1)
  write.table(formatC(eval(as.name(numc)), digits = 1, format = "f"), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.prop1 <- c("phen.mv.greater.p", "phen.uv.greater.p", "phen.uv.equal.p")
for(numc in save.prop1)
  write.table(formatC(100 * eval(as.name(numc)), digits = 0, format = "f"), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.prop2 <- c("line.mv.greater.p", "line.uv.greater.p", "line.uv.equal.p")
for(numc in save.prop2)
  write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)

#####################################
#Plot of discordance vs FDR
qseq <- seq(0, .3, len = 100000)
fdrseq1 <- (4 - sqrt(4^2 - 4 * 3 * 2 * qseq)) / (2 * 3)
fdrseq2 <- (4 + sqrt(4^2 - 4 * 3 * 2 * qseq)) / (2 * 3)
fnamc <- "discordance_fdr.jpg"
# pdf(paste(power.plot.dir, "/mv_vs_uv_power.pdf", sep = ""), 12, 12)
jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 6, 6, units = "in", res = 500)
par(mfrow = c(1, 1))
par(mar = c(6, 6, 2, 2))
plot(qseq, fdrseq1, ty = "l", ylab = "", xlab = "", las = 1, xaxs = "i", yaxs = "i")
# abline(0, 1, lty = 3)
mtext(side = 1, text = "P(two methods discordant | both methods annotate)", line = 4)
mtext(side = 2, text = "Discordance-implied FDR", line = 4)
dev.off()
file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)






######################################################
######################################################
# Combined power plot
######################################################
######################################################

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










more.power.plots <- F

if(more.power.plots){
  fnamc <- "line_by_line_uv_mv_comp.jpg"
  # pdf(paste(power.plot.dir, "/mv_vs_uv_power.pdf", sep = ""), 12, 12)
  jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), oma = c(2, 2, 2, 2))
  labs <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  for(plty in c("non", "imp")){
    xc <- log(1 + n.uv)
    yc <- switch(plty, non = log(1 + n.eb.non),
                 imp = log(1 + n.eb.non + n.eb.imp))
    lims <- c(0, log(1000))#range(c(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp)))
    repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
    symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) / 30, inches = F, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, las = 2,
            xlab = "UV analysis, number of annotations per phenotype (log scale)", ylab = "MV analysis, number of annotations per phenotype (log scale)", cex = .5)
    # plot(xc, yc, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, cex = .5,
    #      xlab = "UV analysis, number of annotations per phenotype (log scale)", ylab = "MV analysis, number of annotations per phenotype (log scale)")
    abline(0, 1)
    for(ax in 1:2) 
      axis(side = ax, at = log(1 + labs), labels = labs)
  }
  labs <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  for(plty in c("non", "imp")){
    xc <- log(1 + line.n.uv)
    yc <- switch(plty, non = log(1 + line.n.eb.non),
                 imp = log(1 + line.n.eb.non + line.n.eb.imp))
  #  lims <- range(c(0, log(1 + line.n.uv), log(1 + line.n.eb.non + line.n.eb.imp)))
    repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
    symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) / 30, inches = F, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, las = 2, 
            xlab = "UV analysis, number of annotations per line (log scale)", ylab = "MV analysis, number of annotations per line (log scale)", cex = .5)
    #plot(xc, yc, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, xlab = "Number of ")
    abline(0, 1)
    for(ax in 1:2) 
      axis(side = ax, at = log(1 + labs), labels = labs)
  }
  dev.off()
  file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
            to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
  
  
  
  
  
  labs <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  for(plty in c("non", "imp")){
    xc <- log(1 + n.uv)
    yc <- switch(plty, non = log(1 + n.eb.non),
                 imp = log(1 + n.eb.non + n.eb.imp))
    lims <- c(0, log(1000))#range(c(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp)))
    repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
    symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) / 30, inches = F, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, las = 2,
            xlab = "UV analysis, number of annotations per phenotype (log scale)", ylab = "MV analysis, number of annotations per phenotype (log scale)", cex = .5)
    # plot(xc, yc, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, cex = .5,
    #      xlab = "UV analysis, number of annotations per phenotype (log scale)", ylab = "MV analysis, number of annotations per phenotype (log scale)")
    abline(0, 1)
    for(ax in 1:2) 
      axis(side = ax, at = log(1 + labs), labels = labs)
  }
  labs <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  for(plty in c("non", "imp")){
    xc <- log(1 + line.n.uv)
    yc <- switch(plty, non = log(1 + line.n.eb.non),
                 imp = log(1 + line.n.eb.non + line.n.eb.imp))
    #  lims <- range(c(0, log(1 + line.n.uv), log(1 + line.n.eb.non + line.n.eb.imp)))
    repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
    symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) / 30, inches = F, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, las = 2, 
            xlab = "UV analysis, number of annotations per line (log scale)", ylab = "MV analysis, number of annotations per line (log scale)", cex = .5)
    #plot(xc, yc, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, xlab = "Number of ")
    abline(0, 1)
    for(ax in 1:2) 
      axis(side = ax, at = log(1 + labs), labels = labs)
  }
  dev.off()
  
  
  
  # 
  # 
  # 
  # plot(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp), xaxt = "n", yaxt = "n")
  # for(ax in 1:2) 
  #   axis(side = ax, at = log(1 + labs), labels = labs)
  # 
  # 
  # 
  # 
  # grp.max <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
  # grp.min <- c(0, 1, 2, grp.max[3:(length(grp.max) - 1)] + 1)
  # grpnam <- apply(cbind(grp.min, grp.max), 1, function(v) paste(unique(v), collapse = "-"))
  # grpnamtab <- data.frame(grpnam = grpnam, grp.min = grp.min, grp.max = grp.max)
  # grpdat <- data.frame(geno = trueline.un, n.uv = line.n.uv, n.eb.non = line.n.eb.non, n.eb.imp = line.n.eb.imp,
  #                      n.eb.both = line.n.eb.non + line.n.eb.imp)
  # grpdat$uv.grp.num <- apply(outer(grpdat$n.uv, grp.max, "<=") * outer(grpdat$n.uv, grp.min, ">="), 1, function(v) match(T, v))
  # grpdat$eb.non.grp.num <- apply(outer(grpdat$n.eb.non, grp.max, "<=") * outer(grpdat$n.eb.non, grp.min, ">="), 1, function(v) match(T, v))
  # 
  # image(log(1 + table(grpdat$uv.grp.num, grpdat$eb.non.grp.num)))
  # image()
  # 
  # hist(line.n.uv)
  # hist(line.n.eb.non)
  
  
  pdf("X:/projects/impc_mv_analysis/plots/mv_uv_scatter.pdf", 12, 12)
  lc = 5
  colv <- rep(1, nrow(resimp))
  colv[abs(resimp$uv.t) > resimp$uv.th.final & abs(resimp$eb.t) < resimp$eb.th.final] <- 2
  colv[abs(resimp$uv.t) < resimp$uv.th.final & abs(resimp$eb.t) > resimp$eb.th.final] <- 3
  colv[abs(resimp$uv.t) > resimp$uv.th.final & abs(resimp$eb.t) > resimp$eb.th.final] <- 4
  plot(resimp$uv.t[order(colv)], resimp$eb.t[order(colv)], xlim = c(-1, 1) * lc, ylim = c(-1, 1) * lc, xlab = "UV", 
       ylab = "uv", pch = ".", col = colv[order(colv)], cex = 1.5)
  abline(0, 1)
  #abline(h = c(-1, 1) * mean(resimp$eb.th.final, na.rm = T), v = c(-1, 1) * mean(resimp$uv.th.final, na.rm = T))
  dev.off()
  
  for(curpl in c("null", "true")){
    fnamc <- paste("meth_", estimation.meth, "_", curpl, "_uv_mv_comparison.jpeg", sep = "")
    jpeg(paste(chris.pres.dropbox, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
    par(mfrow = c(2, 2))
    resimpc <- switch(curpl,
                      null = resimp[resimp$line.type == nam.negcon & !is.na(resimp$uv.t), ],
                      true = resimp[resimp$line.type == nam.truemut & !is.na(resimp$uv.t), ])
    plot(resimpc$uv.mn, resimpc$eb.mn)
    abline(0, 1)
    plot(resimpc$uv.sd, resimpc$eb.sd)
    abline(0, 1)
    plot(resimpc$uv.t, resimpc$eb.t)
    abline(0, 1)
    dev.off()
  }
  
  smoothScatter(resimp$uv.t / resimp$uv.th.final, resimp$eb.t / resimp$eb.th.final, pch = ".", 
                xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, nbin = 1000)
  abline(h = c(-1, 1), col = 2)
  abline(v = c(-1, 1), col = 2)
  
  
  
  print(table((abs(resimp$uv.t) > resimp$uv.th.final), 
              (abs(resimp$eb.t) > resimp$eb.th.final)) / sum(!is.na(resimp$uv.t)))
  numtab <- table((abs(resimp$uv.t) > resimp$uv.th.final), 
                  (abs(resimp$eb.t) > resimp$eb.th.final))
  mean(resimp$eb.sig[is.na(resimp$uv.t)], na.rm = T)
  # table((abs(resimp$uv.t) > resimp$uv.th.final), 
  #       (abs(resimp$mv.t) > resimp$mv.th.final))
  # table((abs(resimp$mv.t) > resimp$mv.th.final), 
  #       (abs(resimp$eb.t) > resimp$eb.th.final))
  sum(resimp$eb.sig[is.na(resimp$uv.t)], na.rm = T)
  
  mean(resimp$eb.sig[is.na(resimp$uv.t)], na.rm = T)
  fpr.true.imp <- mean(resimp$eb.sig[is.na(resimp$uv.t) & resimp$line.type == nam.truemut], na.rm = T)
  fpr.neg.imp <- mean(resimp$eb.sig[is.na(resimp$uv.t) & resimp$line.type == nam.negcon], na.rm = T)
  fpr.neg.imp / fpr.true.imp
  fpr.true.non <- mean(resimp$eb.sig[!is.na(resimp$uv.t) & resimp$line.type == nam.truemut], na.rm = T)
  fpr.neg.non <- mean(resimp$eb.sig[!is.na(resimp$uv.t) & resimp$line.type == nam.negcon], na.rm = T)
  fpr.neg.non / fpr.true.non
  
  
  
  
  #t-statistic QQ-plot for Chris' talk
  fnamc <- paste("t_statistic_qqplot_mv_uv.jpeg", sep = "")
  jpeg(paste(chris.pres.dropbox, "/", fnamc, sep = ""), 7, 7, units = "in", res = 500)
  par(mfrow = c(1, 1), mar = c(5, 5, 5, 5))
  cexax <- 2
  linc <- 3.5
  cexax1 <- 1.4
  colv <- rep(1:2, times = c(length(qq.mv$x), length(qq.uv$x)))
  plot(c(qq.mv$x, qq.uv$x), c(qq.mv$y, qq.uv$y), ty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  axis(side = 1, cex.axis = cexax1, las = 1)
  axis(side = 2, cex.axis = cexax1, las = 1)
  points(c(qq.mv$x, qq.uv$x), c(qq.mv$y, qq.uv$y), cex = .6, col = colv)
  mtext(side = 1, text = "Permutation null t quantile", cex = cexax, line = linc)
  mtext(side = 2, text = "KO gene t quantile", cex = cexax, line = linc)
  legend(x = "topleft", legend = c("MV", "UV"), pch = 1, col = 1:2, cex = 1.4)
  abline(0, 1)
  dev.off()
  
  
  cdf.mv.qval <- ecdf(p.adjust(p.mv, meth = "BH"))
  cdf.uv.qval <- ecdf(p.adjust(p.uv, meth = "BH"))
  cdf.mv <- ecdf(p.mv)
  cdf.uv <- ecdf(p.uv)
  pseq <- 10^seq(-6, -2, len = 1000)
  matplot(pseq, cbind(cdf.mv(pseq), cdf.uv(pseq)), ty = "l", log = "x", lty = 1)
  matplot(pseq, cbind(cdf.mv.qval(pseq), cdf.uv.qval(pseq)), ty = "l", log = "x", lty = 1)
  
  
  
  tseq <- 10^seq(0, 1, len = 1000)
  matplot(tseq, 1 - cbind(ecdf(abs(true.t.mv))(tseq), ecdf(abs(true.t.uv))(tseq)), ty = "l", log = "", lty = 1)
  
  
  
  # fnamc <- paste("perm_pvalue_histograms_densities.jpeg", sep = "")
  # jpeg(paste(chris.pres.dropbox, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
  par(mfrow = c(2, 2))
  hist(log10(p.mv), main = "MV Permutation p-values", xlab = "p")
  hist(log10(p.uv), main = "UV Permutation p-values", xlab = "p")
  plot(0, ty = "l", xlim = c(0, 1), ylim = c(0, max(c(d.mv$y, d.uv$y))), xaxs = "i", xlab = "p", ylab = "Density")
  lines(x = d.mv$x, y = d.mv$y)
  lines(x = d.uv$x, y = d.uv$y, col = "red")
  legend(x = "topright", legend = c("MV perm. p-values", "UV perm. p-values"), lty = 1, col = 1:2)
  # dev.off()
  
  
  par(mfrow = c(2, 2))
  qq.mv <- qqplot(null.t.mv, true.t.mv, plot.it = FALSE)
  qq.uv <- qqplot(null.t.uv, true.t.uv, plot.it = FALSE)
  
  # dc <- d[d$ph == "IMPC_DXA_010_001", ]
  # mean(table(dc$geno))
  
  
  namc <- "isba_pres_mv_vs_uv_power.jpeg"
  jpeg(file = paste(power.plot.dir, "/", namc, sep = ""), 11.5, 6, units = "in", res = 500)
  # pdf(paste(power.plot.dir, "/", namc, sep = ""), 12, 6)
  par(mfrow = c(1, 2), mar = c(6, 6, 4, 4), oma = c(0, 0, 0, 0))
  plty <- "non"
  xc <- n.uv
  yc <- switch(plty, non = n.eb.non,
               imp = n.eb.non + n.eb.imp)
  lims <- range(c(xc, yc))#range(c(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp)))
  repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
  circmult1 <- 5
  circmult2 <- 1 / 3
  cexax <- 1.4
  cexax2 <- 1.7
  axnumcex <- 1.4
  symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) * circmult1, inches = F, xlim = lims, ylim = lims, las = 1,
          xlab = "", ylab = "", cex.axis = axnumcex)
  linc <- 3.5
  linc2 <- 2
  mtext(side = 1, line = linc, text = "UV analysis", cex = cexax)
  mtext(side = 2, line = linc, text = "MV analysis", cex = cexax)
  mtext(side = 3, line = linc2, text = "Annotations per phenotype", cex = cexax2)
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
  dev.off()
  file.copy(from = paste(power.plot.dir, "/", namc, sep = ""),
            to = paste(chris.pres.dropbox, "/", namc, sep = ""), overwrite = TRUE)
}





# pow.uv <- sapply(phen.un, function(ph) mean(resimp[resimp$ph == ph, "uv.sig"][!is.na(resimp$uv.t)], na.rm = T))
# pow.eb.non <- sapply(phen.un, function(ph) mean(resimp[resimp$ph == ph, "eb.sig"][!is.na(resimp$uv.t)], na.rm = T))
# pow.eb.imp <- sapply(phen.un, function(ph) mean(resimp[resimp$ph == ph, "eb.sig"][is.na(resimp$uv.t)], na.rm = T))
