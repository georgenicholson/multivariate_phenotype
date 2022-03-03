resl.comp.fac <- readRDS(file = control$file_raw_factor_results)
resimp <- readRDS(file = control$file.resimp)
# resl <- resl.comp
resl.comp$eb <- resl.comp[[control$mv_meth_nam_use]]
sig <- resl.comp$eb$Sig.comb
fac.meth <- c("varimax", "promax")[1]
facs <- resl.comp$eb$facs.varimax
fl <- list()
fl$t <- resl.comp.fac[[control$mv_meth_nam_use]]$mn / resl.comp.fac[[control$mv_meth_nam_use]]$sd
fl$th <- fl$t
fl$th[] <- resimp[, paste0(fac.meth, ".th.final")][1]
if(length(unique(resimp[, paste0(fac.meth, ".th.final")])) > 1)
  stop("Multiple factor sig thresh, but code written for single thresh")

sigcor <- t(t(sig / sqrt(diag(sig))) / sqrt(diag(sig)))
eval <- eigen(sigcor)$values
graphics.off()
cumexpl <- cumsum(eval) / sum(eval)
corr.explained <- (cumsum(eval) / sum(eval))[control$nfac]
tabc <- data.frame(n = 1:ncol(sig), cump = cumexpl)
nfac <- control$nfac
facnam <- paste0("fac_", 1:control$nfac)

propexp <- tabc[match(nfac, tabc$n), "cump"]

####################################################
#plot cumulative correlation explained
#######################################
graphics.off()
fnamc <- "cumulative_correlation_explained.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 6, 6, units = "in", res = 1000)
par(mar = c(5, 5, 2, 2))
plot(c(0, 1:ncol(sig)), c(0, cumexpl), ty = "l", xlab = "Number of eigenvectors", 
     ylab = "", ylim = c(0, 1), xaxs = "i", yaxs = "i", las = 1)
lines(x = c(nfac, nfac, 0), y = c(0, propexp, propexp), lty = 3)
# axis(side = 2, at = propexp, labels = formatC(propexp, format = "f", digits = 3), las = 2, col = 1)
mtext(side = 2, text = "Cumulative proportion of correlation explained", line = 4)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

save.num <- c("nfac")
for(numc in save.num)
  write.table(nfac, file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
save.prop <- c("corr.explained")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 0, format = "f"), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)

##############################################################################
#Calculate odds ratio pairwise between factors for use in text and figure
############################################################
true.use <- unique(resimp[resimp$line.type == "trueMut", "geno"])
sigmat.in <- t(t((sign(fl$t) * (abs(fl$t) > fl$th))[true.use, facnam]))
fac.load.sign.adjust <- sign(colMeans(sigmat.in))
fac.load.sign.adjust[fac.load.sign.adjust == 0] <- 1
sigmat <- t(t(sigmat.in) * fac.load.sign.adjust)
sigmat <- sigmat[, order(-colMeans(abs(sigmat)))]
facord <- colnames(sigmat)
propup <- colMeans(sigmat == 1)
propdo <- colMeans(sigmat == -1)
propall <- colMeans(sigmat != 0)

all.dir <- opp.dir <- same.dir <- ormat <- orpmat <- matrix(NA, nfac, nfac)
for(i in 1:(nfac - 1)){
  for(j in (i + 1):nfac){
    all.dir[i, j] <- all.dir[j, i] <- (mean(sigmat[, i] == 1 & sigmat[, j] == 1) + mean(sigmat[, i] == -1 & sigmat[, j] == -1)
                                       + mean(sigmat[, i] == 1 & sigmat[, j] == -1) + mean(sigmat[, i] == -1 & sigmat[, j] == 1)) /
      (propup[i] * propup[j] + propdo[i] * propdo[j] + propup[i] * propdo[j] + propdo[i] * propup[j])
  }
}

for(i in 1:(nfac - 1)){
  for(j in (i + 1):nfac){
    n00 <- sum(sigmat[, i] == 0 & sigmat[, j] == 0)
    n11 <- sum(sigmat[, i] != 0 & sigmat[, j] != 0)
    n01 <- sum(sigmat[, i] == 0 & sigmat[, j] != 0)
    n10 <- sum(sigmat[, i] != 0 & sigmat[, j] == 0)
    ormat[i, j] <- ormat[j, i] <- n11 * n00 / n01 / n10
    matc <- matrix(c(n00, n01, n10, n11), 2, 2)
    orpmat[i, j] <- orpmat[j, i] <- fisher.test(matc)$p.value
  }
}

n.fac.pr.test <- .5 * nfac * (nfac - 1)
orqmat <- p.adjust(orpmat[lower.tri(orpmat)], method = "BY")# * n.fac.pr.test
fac.sig.th <- .05
n.fac.pr.sig <- sum(orqmat < fac.sig.th, na.rm = T)
# n.fac.pr.sig <- sum(orqmat < fac.sig.th, na.rm = T) / 2

#order phenotypes within procedure for factor plot
phord.fac <- Data_all$impc$phord
for(procc in Data_all$impc$procord){
  phunc <- as.matrix(Data_all$impc$phmap[Data_all$impc$phmap$procnam == procc, "ph"])
  phunc <- phunc[phunc %in% Data_all$impc$phord]
  if(length(phunc) > 2){
    phunc.ord <- phunc[hclust(dist(abs(facs[phunc, ]), meth = "manhattan"))$order]
    Data_all$impc$phord[Data_all$impc$phord %in% phunc] <- phunc.ord
  }
}

###################################################################
#Export numbers for text, particularly the number of pairwise tests and the number significant
###################################################################
max.prop.annot <- max(propall)
min.prop.annot <- min(propall)
mean.prop.up <- mean(propup / propall)
save.num <- c("max.prop.annot", "min.prop.annot", "mean.prop.up")
for(numc in save.num)
  write.table(formatC(eval(as.name(numc)) * 100, format = "f", digits = 1), 
              file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
save.num <- c("n.fac.pr.test", "n.fac.pr.sig")
for(numc in save.num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 0), 
              file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)



fac.maxload <- apply(facs[phord.fac, facord], 2, function(v) max(abs(v)))
corfac <- cor(sigmat)
nloadlook <- 10
zpl <- t(facs[phord.fac, facord]) * fac.load.sign.adjust[facord] / fac.maxload[facord]
phmat <- apply(zpl, 1, function(v) colnames(zpl)[order(-abs(v))[1:nloadlook]])
toploadmat <- apply(zpl, 1, function(v) v[order(-abs(v))[1:nloadlook]])
nammat <- phmat
nammat[] <- Data_all$impc$phmap[match(phmat, Data_all$impc$phmap$ph), "nam"]
nph <- length(phord.fac)
tab.fac.interp <- data.frame(phen.name = c(nammat), loading = round(c(toploadmat), 2),
                             factor.num = rep(1:nfac, each = nloadlook))
facmap <- data.frame(original = facnam, reordered = facord)
tab.fac.interp$factor_num_original <- rep(facord, each = nloadlook)#facmap[match(paste0("fac_", tab.fac.interp$factor.num), facmap$reordered), "original"]
write.csv(tab.fac.interp, file = "output/factor_interp_to_be_completed.csv")

fac_annot_tab <- read.csv(file = "output/factor_interp_completed.csv", stringsAsFactors = FALSE)
facannot <- unique(fac_annot_tab$factor.annotation)
names(facannot) <- paste0("ord", 1:nfac)


######################################################
#Plot factor interpretation panel figure for paper
################################################
library("KeyboardSimulator")
jpegc <- T
devwid <- 9
devhei <- 11
if(jpegc){
  fnamc <- "factor_interpretation_plot.jpg"
  jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), devwid, 11, units = "in", res = 1000)
} else {
  graphics.off()
  windows(devwid, 11, xpos = 1250, ypos = 300)
}
wc <- .6
hc <- .8
laymat <- matrix(c(1, 2, 4, 3), 2, 2)
layl <- list(p1 = c(2, 2), p2 = c(3, 2), p3 = c(3, 3))
npl <- length(layl)
laymat <- matrix(npl + 1, 4, 4)
for(i in 1:length(layl))
  laymat[t(layl[[i]])] <- i
laymat
vp <- c(.025, .6, .2, .2)
hp <- c(.3, .25, .2, .3)
layout(laymat, hei = vp, wid = hp)
oma.use <- rep(0, 4)
par(mar = rep(.5, 4), oma = oma.use, xpd = F)
propsig <- colMeans(sigmat[, facord])
barmat <- 100 * t(cbind(propup, propdo))
cexax <- .9
cexax2 <- .9
cexax3 <- .95
cexax4 <- .8
cexlet <- 1.15
cexleg <- 1
cexax.labrt <- .6
cexax.lablt <- 1
cex.axis <- .9
cexax.mtext.big <- .8

###############################
#Panel (a)
#######################
image(x = 1:nfac, y = 1:nph, z = zpl, col = control$heat_col_palette, yaxt = "n", ylab = "", cex.axis = cexax, xaxt = "n")
ypl1 <- grconvertY(c(0, 1), from = "nfc", "ndc")
mtext(side = 3, at = -2, line = .5, text = "(a)", cex = cexlet)
abline(v = 1:nfac - .5)
procv <- Data_all$impc$phmap[match(phord.fac, Data_all$impc$phmap$ph), "procnam"]
procats <- sapply(Data_all$impc$procord, function(procc) mean(which(procv == procc)))
axis(side = 2, labels = Data_all$impc$procord, at = procats, las = 2, cex.axis = cexax.lablt)
phnam <- Data_all$impc$phmap[match(phord.fac, Data_all$impc$phmap$ph), "nam"]
odds <- seq(1, nph, by = 2)
evens <- seq(2, nph, by = 2)
ats <- 1:length(phord.fac)
phnamsh <- phnam
lineout <- 17
axis(side = 4, labels = NA, at = ats[odds], las = 2, cex.axis = cexax2)
axis(side = 4, labels = NA, at = ats[evens], las = 2, cex.axis = cexax2, tcl = -(lineout - .5), col.ticks = "grey")
mtext(side = 4, line = 1, text = phnamsh[odds], at = ats[odds], las = 2, cex = cexax.labrt, adj = 0)
mtext(side = 4, line = lineout, text = phnamsh[evens], at = ats[evens], las = 2, cex = cexax.labrt, adj = 0)
abline(h = match(Data_all$impc$procord, procv) - .5, lwd = 2)

###############################
#Panel (b)
#######################
orth <- 10
ormatpl <- ormat
ormatpl[ormatpl > orth] <- orth
image(x = 1:nfac, y = 1:nfac, log(ormatpl), zlim = c(-1, 1) * log(orth),
      xaxt = "n", yaxt = "n", col = control$heat_col_palette, xlab = "", ylab = "")
ypl2 <- grconvertY(c(0, 1), from = "nfc", "ndc")
mtext(side = 1, at = -3, line = 1, text = "(b)", cex = cexlet)

for(j in 1:nrow(zpl)){
  phvc <- colnames(zpl)[order(-abs(zpl[j, ]))[1:10]]
  phnamvc <- Data_all$impc$phmap[match(phvc, Data_all$impc$phmap$ph), "nam"]
  datc <- data.frame(nam = phnamvc, val = zpl[j, phvc])
}
axis(side = 2, labels = facannot, at = 1:nfac, las = 2, cex.axis = cexax3)
axis(side = 1, labels = facannot, at = 1:nfac, las = 2, cex.axis = cexax3)

###############################
#Panel (c)
#######################
barplot(barmat, col = c("red", "blue"), xaxs = "i", horiz = TRUE, yaxs = "i",
        xaxt = "s", xlab = "", ylab = "", las = 1, cex.axis = cex.axis, yaxt = "n", las = 2)
mtext(side = 1, at = max(barmat) * 1.3, line = 1, text = "(c)", cex = cexlet)
mtext(side = 1, text = "% lines annotated", cex = cexax.mtext.big, line = 2)
# legend(x = "topright", legend = c("Positive", "Negative"), title = "Perturbation in factor score", pch = 22, pt.bg = c("red", "blue"), col = 1,
#        cex = cexleg)

barwid <- .015
barhei <- .05
barx2 <- barx1 <- .095
bary1 <- mean(ypl1)
bary2 <- mean(ypl2)

ax.mult <- .75
lin.mult <- .65

###############################
#Panel (a) scale
#######################
par(fig = c(barx1, barx1 + barwid, bary1 - barhei, bary1 + barhei), new = T)
cexax <- 1.1
linec <- 4
par(mar = c(0, 0, 0, 0))
image(z = t(as.matrix(1:1000)), y = seq(-1, 1, len = 1000), x = 1, col = control$heat_col_palette, xaxt = "n", yaxt = "n",
      ylab = "")
axis(side = 2, las = 2, cex.axis = cexax4, labels = c("-1.0", -0.5, 0.0, 0.5, "1.0"), at = seq(-1, 1, by = .5), las = 1)
mtext(side = 2, line = linec, text = 'Loadings', cex = cexax.mtext.big, las = 0)
mtext(side = 2, line = linec * lin.mult, text = "in panel (a)", cex = cexax.mtext.big * ax.mult, las = 0)

###############################
#Panel (b) scale
#######################
par(fig = c(barx2, barx2 + barwid, bary2 - barhei, bary2 + barhei), new = T)
par(mar = c(0, 0, 0, 0))
range(log(ormatpl))
image(z = t(as.matrix(1:1000)), y = seq(-1 * log(orth), 1 * log(orth), len = 1000), x = 1, col = control$heat_col_palette, xaxt = "n", yaxt = "n",
      ylab = "")
if(orth == 10){
  orat <- c(.1, .2, .5, 1, 2, 5, 10)
  oratlab <- c('< 0.1', .2, .5, 1, 2, 5, '> 10')
}
if(orth == 50){
  orat <- c(.02, .05, .1, .2, .5, 1, 2, 5, 10, 20, 50)
  oratlab <- c('< 0.02', .05, .1, .2, .5, 1, 2, 5, 10, 20, '> 50')
}
if(orth == 20){
  orat <- c(.05, .1, .2, .5, 1, 2, 5, 10, 20)
  oratlab <- c('< 0.05', .1, .2, .5, 1, 2, 5, 10, '> 20')
}

axis(side = 2, las = 2, cex.axis = cexax4, labels = oratlab, at = log(orat), las = 1)
mtext(side = 2, line = linec, text = "Odds Ratio", cex = cexax.mtext.big, las = 0)
mtext(side = 2, line = linec * lin.mult, text = "in panel (b)", cex = cexax.mtext.big * ax.mult, las = 0)
# 
if(jpegc){
  dev.off()
    file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
              to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
}
# 
# image(x = 1:nfac, y = 1:nfac, log(ormatpl), zlim = c(-1, 1) * log(orth),
#       xaxt = "n", yaxt = "n", col = control$heat_col_palette, xlab = "", ylab = "")
# keybd.press('ctrl+alt+u')









