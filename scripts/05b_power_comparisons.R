resimp <- readRDS(file = control$file.resimp)
phen.un <- unique(resimp$ph)
phen.un <- phen.un[!phen.un %in% paste0("fac_", 1:control$nfac)]
resimp <- resimp[resimp$ph %in% phen.un, ]

trueline.un <- unique(resimp[resimp$line.type == control$nam.truemut, "geno"])
n.uv <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & !is.na(resimp$uv.t), "uv.sig"], na.rm = T))
n.eb.non <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & !is.na(resimp$uv.t), "eb.sig"], na.rm = T))
n.eb.imp <- sapply(phen.un, function(ph) sum(resimp[resimp$ph == ph & is.na(resimp$uv.t), "eb.sig"], na.rm = T))

###############################################
#First set of numbers for text
n.ph.measured <- length(phen.un)
n.ko.line.measured <- length(trueline.un)
prop.miss <- mean(is.na(Data_all$impc$Y_raw[trueline.un, ]))
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
mean.s.measured.phens <- mean(Data_all$impc$S_raw, na.rm = T)
sd.y.measured.phens <- sd(Data_all$impc$Y_raw, na.rm = T)
for(dir.save in c(control$dropbox_text_numbers_dir)){
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
              file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
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
  write.table(loomv.uv.ci, file = paste(control$dropbox_text_numbers_dir, "/", save.text, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


###############################################
#UV/MV global scatterplot
fnamc <- "uv_mv_scatter.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
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
    labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
    legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .85), xjust = .5)
  }
}
numv <- formatC(c(fdrci[[1]], fdrci[[2]]) * 100, format = "f", digits = 1)
fdrlab <- paste("Discordance-implied FDR ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
# mtext(side = 3, text = fdrlab, line = .5, cex = cexax2)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

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
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.prop <- c("uv.mv.agree.p", "uv.mv.agree.p.of.uv", "uv.not.mv.p", "mv.not.uv.p", "uv.mv.disagree.p", "uv.mv.disagree.p.of.uv")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)

###############################################
#UV/MV global line-by-line, phen-by-phen scatterplot
resimp1 <- resimp[resimp$line.type == control$nam.truemut, ]
resimp1 <- resimp1[order(resimp1$geno), ]
trueline.un <- unique(resimp1$geno)
froms <- match(trueline.un, resimp1$geno)
tos <- c(froms[2:length(froms)] - 1, nrow(resimp1))

signifl <- Map(function(from, to) resimp1[from:to, c("uv.sig", "eb.sig")], from = froms, to = tos)
line.n.uv <- sapply(signifl, function(x) sum(x$uv.sig, na.rm = T))
line.n.eb.non <- sapply(signifl, function(x) sum(x$eb.sig[!is.na(x$uv.sig)], na.rm = T))
line.n.eb.imp <- sapply(signifl, function(x) sum(x$eb.sig[is.na(x$uv.sig)], na.rm = T))
fnamc <- "line_by_line_uv_mv_comp.jpg"
jpeg(file = paste(control$figure_dir, "/", fnamc, sep = ""), 11.5, 6, units = "in", res = 500)
par(mfrow = c(1, 2), mar = c(6, 6, 4, 4), oma = c(0, 0, 0, 0))
plty <- "non"
xc <- n.uv
yc <- switch(plty, non = n.eb.non,
             imp = n.eb.non + n.eb.imp)
lims <- range(c(xc, yc))
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
mtext(side = 3, line = linc2, text = "Hits per phenotype", cex = cexax2)
mtext(side = 3, line = linc2, at = 0, text = "(a)", cex = cexax3)
abline(0, 1)

xc <- line.n.uv
yc <- switch(plty, non = line.n.eb.non,
             imp = line.n.eb.non + line.n.eb.imp)
lims <- range(c(xc, yc))
repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) * circmult2, inches = F, xlim = lims, ylim = lims, las = 1,
        xlab = "", ylab = "", cex.axis = axnumcex)
abline(0, 1)
mtext(side = 1, line = linc, text = "UV analysis", cex = cexax)
mtext(side = 2, line = linc, text = "MV analysis", cex = cexax)
mtext(side = 3, line = linc2, text = "Hits per KO gene", cex = cexax2)
mtext(side = 3, line = linc2 / 3, text = "(point area proportional to number)", cex = cexax2 * .75)
mtext(side = 3, line = linc2, at = 0, text = "(b)", cex = cexax3)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

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
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.num1 <- c("phen.average.increase.n", "line.average.increase.n")
for(numc in save.num1)
  write.table(formatC(eval(as.name(numc)), digits = 1, format = "f"), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.prop1 <- c("phen.mv.greater.p", "phen.uv.greater.p", "phen.uv.equal.p")
for(numc in save.prop1)
  write.table(formatC(100 * eval(as.name(numc)), digits = 0, format = "f"), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
save.prop2 <- c("line.mv.greater.p", "line.uv.greater.p", "line.uv.equal.p")
for(numc in save.prop2)
  write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)

#####################################
#Plot of discordance vs FDR
qseq <- seq(0, .3, len = 100000)
fdrseq1 <- (4 - sqrt(4^2 - 4 * 3 * 2 * qseq)) / (2 * 3)
fdrseq2 <- (4 + sqrt(4^2 - 4 * 3 * 2 * qseq)) / (2 * 3)
fnamc <- "discordance_fdr.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 6, 6, units = "in", res = 500)
par(mfrow = c(1, 1))
par(mar = c(6, 6, 2, 2))
plot(qseq, fdrseq1, ty = "l", ylab = "", xlab = "", las = 1, xaxs = "i", yaxs = "i")
mtext(side = 1, text = "P(two methods discordant | both methods annotate)", line = 4)
mtext(side = 2, text = "Discordance-implied FDR", line = 4)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

######################################################
######################################################
# Combined power plot
######################################################
######################################################

###############################################
#UV/MV global line-by-line, phen-by-phen scatterplot
resimp1 <- resimp[which(resimp$line.type == control$nam.truemut), ]
resimp1 <- resimp1[order(resimp1$geno), ]
trueline.un <- unique(resimp1$geno)
froms <- match(trueline.un, resimp1$geno)
tos <- c(froms[2:length(froms)] - 1, nrow(resimp1))
signifl <- Map(function(from, to) resimp1[from:to, c("uv.sig", "eb.sig")], from = froms, to = tos)
line.n.uv <- sapply(signifl, function(x) sum(x$uv.sig, na.rm = T))
line.n.eb.non <- sapply(signifl, function(x) sum(x$eb.sig[!is.na(x$uv.sig)], na.rm = T))
line.n.eb.imp <- sapply(signifl, function(x) sum(x$eb.sig[is.na(x$uv.sig)], na.rm = T))

wideps <- .5
heieps <- .5
graphics.off()
fnamc <- "combined_power_plot.jpg"
jpeg(file = paste(control$figure_dir, "/", fnamc, sep = ""), 12, 15, units = "in", res = 500)
layout(matrix(c(1:2, 3, 3), 2, 2, byrow = T))
par(oma = c(20, 4, 0, 0))
plty <- "non"
xc <- n.uv
yc <- switch(plty, non = n.eb.non,
             imp = n.eb.non + n.eb.imp)
lims <- range(c(xc, yc))
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
mtext(side = 3, line = linc2, text = "Hits per phenotype", cex = cexax2)
mtext(side = 3, line = linc2, at = 0, text = "(a)", cex = cexax3)
abline(0, 1)

xc <- line.n.uv
yc <- switch(plty, non = line.n.eb.non,
             imp = line.n.eb.non + line.n.eb.imp)
lims <- range(c(xc, yc))
repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) * circmult2, inches = F, xlim = lims, ylim = lims, las = 1,
        xlab = "", ylab = "", cex.axis = axnumcex)
abline(0, 1)
mtext(side = 1, line = linc, text = "UV analysis", cex = cexax)
mtext(side = 2, line = linc, text = "MV analysis", cex = cexax)
mtext(side = 3, line = linc2, text = "Hits per KO gene", cex = cexax2)
mtext(side = 3, line = linc2 / 3, text = "(point area proportional to number)", cex = cexax2 * .75)
mtext(side = 3, line = linc2, at = 0, text = "(b)", cex = cexax3)

###############################################
#Plot estimated power by procedure
resimp$procnam <- phmap[match(resimp$ph, phmap$ph), "procnam"]
procun <- unique(resimp$procnam)
powtab <- data.frame(procnam = procun, uv = NA, eb = NA, loo.eb = NA, n.uv = NA, n.eb.non = NA, n.loo.eb = NA)
rownames(powtab) <- procun
for(procc in procun){
  resc <- resimp[!resimp$imputed & resimp$procnam == procc & resimp$line.type == control$nam.truemut, ]
  powtab[procc, "x.uv"] <- length(unique(resc[which(resc$uv.signsig != 0), "geno"]))
  powtab[procc, "x.eb.non"] <- length(unique(resc[which(resc$eb.signsig != 0), "geno"]))
  powtab[procc, c("n.uv", "n.eb.non")] <- rep(length(unique(resc$geno)), 2)
  resc.imp <- resimp[resimp$imputed & resimp$procnam == procc & resimp$line.type == control$nam.truemut, ]
  powtab[procc, "x.eb.imp"] <- length(unique(resc.imp[which(resc.imp$eb.signsig != 0), "geno"]))
  powtab[procc, "n.eb.imp"] <- length(unique(resc.imp$geno))
}
uv.ci <- binconf(x = powtab$x.uv, n = powtab$n.uv)
powtab[, c("uv.est", "uv.l", "uv.u")] <- uv.ci
eb.non.ci <- binconf(x = powtab$x.eb.non, n = powtab$n.eb.non)
powtab[, c("eb.non.est", "eb.non.l", "eb.non.u")] <- eb.non.ci
eb.imp.ci <- binconf(x = powtab$x.eb.imp, n = powtab$n.eb.imp)
powtab[, c("eb.imp.est", "eb.imp.l", "eb.imp.u")] <- eb.imp.ci
powtab <- powtab[order(powtab$uv.est), ]
ydum <- powtab[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
                   "eb.imp.est", "eb.imp.l", "eb.imp.u")]
nproc <- length(procun)
xpl <- 1:nproc
colv <- c(uv = "black", eb.non = "red", eb.imp = "blue", eb.loocv = "green")
colv
cexlab = 1.3
matplot(x = 1:length(procun), y = ydum, ty = "n", xaxt = "n", xlab = "", ylab = "", xlim = c(.5, nproc + .5), ylim = c(0, max(ydum)),
        xaxs = "i", yaxs = "i", las = 2, cex.axis = cexlab)
eps = .2
points(x = rep(xpl, 3) + rep(c(-.5, -.16, .16) * eps, each = nproc), 
       y = c(powtab$uv.est, powtab$eb.non.est, powtab$eb.imp.est), col = rep(colv[c("uv", "eb.non", "eb.imp")], each = nproc), pch = 19)
for(i in 1:nproc){
  lines(x = rep(i - .5 * eps, 2), y = unlist(powtab[i, c("uv.l", "uv.u")]), col = colv["uv"])
  lines(x = rep(i - .16 * eps, 2), y = unlist(powtab[i, c("eb.non.l", "eb.non.u")]), col = colv["eb.non"])
  lines(x = rep(i + .16 * eps, 2), y = unlist(powtab[i, c("eb.imp.l", "eb.imp.u")]), col = colv["eb.imp"])
}
lines(x = 1:nproc - .5 * eps, y = powtab$uv.est, col = colv["uv"])
lines(x = 1:nproc - .16 * eps, y = powtab$eb.non.est, col = colv["eb.non"])
lines(x = 1:nproc + .16 * eps, y = powtab$eb.imp.est, col = colv["eb.imp"])
abline(v = .5 + 0:nproc)
axis(side = 1, at = 1:nproc, labels = powtab$procnam, las = 2, cex.axis = cexlab)
mtext(side = 2, line = 4, text = "Proportion of KO lines with one or more hits", cex = cexlab)
legend(x = "topleft", legend = c("MV model (observed data)", "MV model (missing data)", "UV model"), 
       col = colv[c("eb.non", "eb.imp", "uv")], lty = 1, cex = cexlab)
mtext(side = 3, line = linc2, at = 0, text = "(c)", cex = cexax3)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

diff.text <- ydum[, "eb.non.est"] - ydum[, "uv.est"]
names(diff.text) <- rownames(ydum)
for(rowc in rownames(ydum))
  dir.save <- paste(control$dropbox_text_numbers_dir, "/proc_pow_up", sep = "")
dir.create(dir.save, showWarnings = F)
#Export numbers to text
for(numc in save.prop)
  for(rowc in rownames(ydum))
    write.table(formatC(100 * diff.text[rowc], digits = 0, format = "f"), file = paste(dir.save, "/", gsub(" ", "_", gsub("/", "", rowc)), ".txt", sep = ""), 
                col.names = F, row.names = F, quote = F)


###############################################
#Plot estimated power by phenotype
##########################################
# Calc power and missingness by phenotype
phenun <- unique(resimp$ph)
powtab.ph <- data.frame(ph = phenun, uv = NA, eb = NA, loo.eb = NA, n.uv = NA, n.eb.non = NA, n.loo.eb = NA)
rownames(powtab.ph) <- phenun
miss_df <- data.frame(ph = phenun, missing_prop = NA)
for (phenc in phenun) {
  resc <- resimp[!is.na(resimp$uv.t) & resimp$ph == phenc & resimp$line.type == "trueMut", ]
  powtab.ph[phenc, "x.uv"] <- length(unique(resc[which(resc$uv.signsig != 0), "geno"]))
  powtab.ph[phenc, "x.eb.non"] <- length(unique(resc[which(resc$eb.signsig != 0), "geno"]))
  powtab.ph[phenc, c("n.uv", "n.eb.non", "n.loo.eb")] <- rep(length(unique(resc$geno)), 3)
  resc.imp <- resimp[resimp$imputed & resimp$ph == phenc & resimp$line.type == "trueMut", ]
  powtab.ph[phenc, "x.eb.imp"] <- length(unique(resc.imp[which(resc.imp$eb.signsig != 0), "geno"]))
  powtab.ph[phenc, "n.eb.imp"] <- length(unique(resc.imp$geno))
  miss_df[match(phenc, miss_df$ph), "missing_prop"] <- nrow(resc) / (nrow(resc) + nrow(resc.imp))
}
uv.ci <- binconf(x = powtab.ph$x.uv, n = powtab.ph$n.uv)
powtab.ph[, c("uv.est", "uv.l", "uv.u")] <- uv.ci
eb.non.ci <- binconf(x = powtab.ph$x.eb.non, n = powtab.ph$n.eb.non)
powtab.ph[, c("eb.non.est", "eb.non.l", "eb.non.u")] <- eb.non.ci
eb.imp.ci <- binconf(x = powtab.ph$x.eb.imp, n = powtab.ph$n.eb.imp)
powtab.ph[, c("eb.imp.est", "eb.imp.l", "eb.imp.u")] <- eb.imp.ci
powtab.ph <- powtab.ph[match(Data_all$impc$phord, powtab.ph$ph), ]
powtab.ph$procnam <- Data_all$impc$phmap[match(powtab.ph$ph, Data_all$impc$phmap$ph), "procnam"]

##########################################
# Plot the data
ydum <- powtab.ph[, c("uv.est", "uv.l", "uv.u", "eb.non.est", "eb.non.l", "eb.non.u", 
                      "eb.imp.est", "eb.imp.l", "eb.imp.u")]
nphen <- length(phenun)
xpl <- 1:nphen
cexlab = 1
fnamc <- "power_by_phenotype.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 20, 12, units = "in", res = 500)
# pdf(paste(control$figure_dir, "/", fnamc, sep = ""), 20, 12)#, units = "in", res = 1000)
hei_prop <- .7
layout(mat = matrix(1:2, 2, 1), heights = c(1 - hei_prop, hei_prop))
par(oma = c(25, 8, 2, 6), mar = c(0, 0, 0, 0))

plot(x = xpl, y = 100 * miss_df$missing_prop, ty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(.5, nphen + .5), ylim = c(0, 100), 
        xaxs = "i", yaxs = "i", las = 2, cex.axis = cexlab)
axis(side = 4, las = 1)
points(x = xpl, y = 100 * miss_df$missing_prop, pch = 19)
lines(x = xpl, y = 100 * miss_df$missing_prop, pch = 19)
abline(v = match(procun, powtab.ph$procnam) - .5, col = "black")
mtext(side = 4, line = 2.5, text = "% of data missing", cex = cexlab)

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
  lines(x = rep(i + .5 * eps, 2), y = unlist(powtab.ph[i, c("eb.imp.l", "eb.imp.l")]), col = colv["eb.imp"])
}
lines(x = 1:nphen - .5 * eps, y = powtab.ph$uv.est, col = colv["uv"])
lines(x = 1:nphen - .0 * eps, y = powtab.ph$eb.non.est, col = colv["eb.non"])
lines(x = 1:nphen + .5 * eps, y = powtab.ph$eb.imp.est, col = colv["eb.imp"])
ats <- sapply(Data_all$impc$procord, function(procc) mean(which(powtab.ph$procnam == procc)))
axis(side = 1, at = ats, labels = Data_all$impc$procord, las = 2, cex.axis = cexlab)
mtext(side = 2, line = 4, text = "Proportion of KO lines annotated", cex = cexlab)
par(xpd = NA)
legend(x = "topleft",  legend = c("MV model (observed data)", "MV model (missing data)", "UV model"), 
       col = colv[c("eb.non", "eb.imp", "uv")], lty = 1,
       cex = cexlab, lwd = 2)
par(xpd = F)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

