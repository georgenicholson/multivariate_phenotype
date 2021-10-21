# rm(list = ls())
# Data <- "impc"
# full.analysis <- T
# source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")
# load(file = file.resl.comp)
# load(file = paste0(global.res.dir, "/resimp_comb.RData"))
# load(file = uv.results.Y.S)
# reflines.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/ref_lines_temp"
# dir.create(reflines.plot.dir, showWarnings = F)

resl.comp.fac <- readRDS(file = control$file_raw_factor_results)
resimp <- readRDS(file = control$file.resimp)


# file.glob.res <- paste0(global.res.dir, "/global_eb_results.RData")
# load(file = file.glob.res)
resl <- resl.comp
resl$eb <- resl[[control$mv_meth_nam_use]]
# str(resl.comp, m = 1)

# load(file = file.resl.comp.fac)
# load(file = procord.file)
# factor.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/factor_plots_temp"
# dir.create(factor.plot.dir, showWarnings = F)


# load(file = file.glob.res)
# nfac <- 24
# vc.type <- c("vari", "pro")[1]
# load(file = paste0(meth.comp.output.dir, "/ebmix_results_nf_", nfac, "_vt_", vc.type, ".RData"))
# load(file = paste0(meth.comp.output.dir, "/ebmix_results.RData"))
# load(file = em.curated.results.file)
# load(file = resimp.with.sig.thresholds.file)

# load(file = file.glob.loadings) #facs
sig <- resl$eb$Sig.comb
fac.meth <- c("varimax", "promax")[1]
# facs <- switch(fac.meth, varimax = facs.varimax, promax = facs.promax)

facs <- resl$eb$facs.varimax

# load(file = fac.res.store.namc)
# facs <- fac.res.store[[fac.meth]]$loadings
# facs <- resl.out$ebf$loadings.ord
fl <- list()
# resl.comp.fac
fl$t <- resl.comp.fac[[control$mv_meth_nam_use]]$mn / resl.comp.fac[[control$mv_meth_nam_use]]$sd
# fl$t <- resl.out[[fac.meth]]$mn / resl.out[[fac.meth]]$sd
fl$th <- fl$t

unique(resimp$ph)
str(resimp)

fl$th[] <- resimp[, paste0(fac.meth, ".th.final")][1]
if(length(unique(resimp[, paste0(fac.meth, ".th.final")])) > 1)
  stop("Multiple factor sig thresh, but code written for single thresh")
# facdat <- resl$f$facdat

sigcor <- t(t(sig / sqrt(diag(sig))) / sqrt(diag(sig)))
eval <- eigen(sigcor)$values
graphics.off()
cumexpl <- cumsum(eval) / sum(eval)
corr.explained <- (cumsum(eval) / sum(eval))[nfac]
tabc <- data.frame(n = 1:ncol(sig), cump = cumexpl)
nfac <- ncol(facs)
propexp <- tabc[match(nfac, tabc$n), "cump"]

# Cormn <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
# eigc <- eigen(Cormn)
# corr.explained <- (cumsum(eigc$values) / sum(eigc$values))[nfac]

####################################################
#plot cumulative correlation explained
#######################################
graphics.off()
fnamc <- "cumulative_correlation_explained.jpg"
jpeg(paste(factor.plot.dir, "/", fnamc, sep = ""), 6, 6, units = "in", res = 1000)
par(mar = c(5, 5, 2, 2))
plot(c(0, 1:ncol(sig)), c(0, cumexpl), ty = "l", xlab = "Number of eigenvectors", 
     ylab = "", ylim = c(0, 1), xaxs = "i", yaxs = "i", las = 1)
lines(x = c(nfac, nfac, 0), y = c(0, propexp, propexp), lty = 3)
# axis(side = 2, at = propexp, labels = formatC(propexp, format = "f", digits = 3), las = 2, col = 1)
mtext(side = 2, text = "Cumulative proportion of correlation explained", line = 4)
dev.off()
file.copy(from = paste(factor.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)



save.num <- c("nfac")
for(numc in save.num)
  write.table(nfac, file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
save.prop <- c("corr.explained")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 0, format = "f"), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
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
orqmat <- orpmat * n.fac.pr.test
fac.sig.th <- .05
n.fac.pr.sig <- sum(orqmat < fac.sig.th, na.rm = T) / 2

#order phenotypes within procedure for factor plot
phord.fac <- phord
for(procc in procord){
  phunc <- as.matrix(pout[pout$procnam == procc, "ph"])
  phunc <- phunc[phunc %in% phord]
  if(length(phunc) > 2){
    phunc.ord <- phunc[hclust(dist(abs(facs[phunc, ]), meth = "manhattan"))$order]
    phord.fac[phord %in% phunc] <- phunc.ord
  }
}
# ormat.clust <- -log(ormat)
# diag(ormat.clust) <- 0
# facord.clust <- paste("f.", 1:nfac, sep = "")[hclust(as.dist(-ormat.clust))$order]
lean.ph <- phmap[phmap$nam == "Lean mass", "ph"]
table(resl.out$eb$signsig[, lean.ph])
table(resl.out$eb.fac$signsig[, "f.1"])
table(sigmat.in[, "f.1"])
table(sign(fl$t[, "f.1"]))
str(phmap)

# 
# write.csv(tab.fac.interp, file = "X:/projects/impc_mv_analysis/data_out/factor_interpretation_table.csv", row.names = F)
# loadmat <- zpl
# save(loadmat, file = "X:/projects/impc_mv_analysis/data_out/loadmat.RData", row.names = F)


###################################################################
#Export numbers for text, particularly the number of pairwise tests and the number significant
###################################################################
max.prop.annot <- max(propall)
min.prop.annot <- min(propall)
mean.prop.up <- mean(propup / propall)
save.num <- c("max.prop.annot", "min.prop.annot", "mean.prop.up")
for(numc in save.num)
  write.table(formatC(eval(as.name(numc)) * 100, format = "f", digits = 1), 
              file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
save.num <- c("n.fac.pr.test", "n.fac.pr.sig")
for(numc in save.num)
  write.table(formatC(eval(as.name(numc)), format = "f", digits = 0), 
              file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)




plot(facs[, "f.1"])
facs[lean.ph, "f.1"]
table(sigmat.in[, "f.1"] )

fac.maxload <- apply(facs[phord.fac, facord], 2, function(v) max(abs(v)))
corfac <- cor(sigmat)
nloadlook <- 10
zpl <- t(facs[phord.fac, facord]) * fac.load.sign.adjust[facord] / fac.maxload[facord]
phmat <- apply(zpl, 1, function(v) colnames(zpl)[order(-abs(v))[1:nloadlook]])
toploadmat <- apply(zpl, 1, function(v) v[order(-abs(v))[1:nloadlook]])
nammat <- phmat
nammat[] <- pout[match(phmat, pout$ph), "nam"]
nph <- length(phord.fac)
tab.fac.interp <- data.frame(phen.name = c(nammat), loading = round(c(toploadmat), 2),
                             factor.num = rep(1:nfac, each = nloadlook))
facannot <- 1:nfac
# facannot <- c("Body size (-)", "Fat mass/insulin (-)", "Cardiopulmonary function (-)", "Bone mineral density (-)",
#               "Cardiac function (-)", "Red blood cell count (-)", "Activity/exploration 1 (+)", "Grip strength (-)", 
#               "Sensorimotor gating (-)", "Liver function (-)", "Blood cholestrol (-)", "Heart rate (+)", 
#               "Balance/coordination (-)", "Hemoglobin (-)", "Deafness (+)", "Activity/exploration 2 (+)", 
#               "Neutrophil & monocyte differential (+)", "White blood cell count (-)", "Kidney function (-)", "Sleep time (+)")  
facannot <- c(ord1 = "Body size (-)", ord2 = "Liver function (-)", ord3 = "???Lymphocyte differential (-)", ord4 = "Blood cholestrol (-)", 
              ord5 = "???Balance/coordination/Sleep (-)", ord6 = "??? Cardiopulmonary function (-)",
              ord7 = "Activity/exploration 1 (+)", ord8 = "Bone mineral density (-)", ord9 = "Sensorimotor gating (-)",
              ord10 = "Deafness (+)", ord11 = "Grip strength (-)", ord12 = "Heart rate (+)", 
              ord13 = "Activity/exploration 2 (+)", ord14 = "Red blood cell count (-)", ord15 = "Respiratory exchange ratio (-)",
              ord16 = "Hemoglobin (-)", ord17 = "White blood cell count (+)", ord18 = "Eosinophil differential (-)",
              ord19 = "??? Anaemia (+)", ord20 = "??? Sleep bout SD (-)")



# "Fat mass/insulin (-)", "Bone mineral density (-)",
# "Cardiac function (-)", ,                 
# "Balance/coordination (-)",  
# "Neutrophil & monocyte differential (+)", "Kidney function (-)", "Sleep time (+)")
tab.fac.interp$factor.annotation <- rep(facannot, each = nloadlook)
tab.fac.interp[]
# save(facord, fac.load.sign.adjust, file = paste(fig.obj.dir, "/facord.RData", sep = ""))

str(zpl[1, ])
# write.csv(tab.fac.interp, file = "X:/projects/impc_mv_analysis/data_out/factor_interpretation_table.csv", row.names = F)


######################################################
#Plot factor interpretation panel figure for paper
################################################
library("KeyboardSimulator")
jpegc <- T
devwid <- 9
devhei <- 11
if(jpegc){
  fnamc <- "factor_interpretation_plot.jpg"
  jpeg(paste(factor.plot.dir, "/", fnamc, sep = ""), devwid, 11, units = "in", res = 1000)
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
image(x = 1:nfac, y = 1:nph, z = zpl, col = rain, yaxt = "n", ylab = "", cex.axis = cexax, xaxt = "n")
ypl1 <- grconvertY(c(0, 1), from = "nfc", "ndc")
mtext(side = 3, at = -2, line = .5, text = "(a)", cex = cexlet)
abline(v = 1:nfac - .5)
procv <- pout[match(phord.fac, pout$ph), "procnam"]
procats <- sapply(procord, function(procc) mean(which(procv == procc)))
axis(side = 2, labels = procord, at = procats, las = 2, cex.axis = cexax.lablt)
phnam <- pout[match(phord.fac, pout$ph), "nam"]
odds <- seq(1, nph, by = 2)
evens <- seq(2, nph, by = 2)
ats <- 1:length(phord.fac)
phnamsh <- phnam
datswap <- data.frame(old = c("Forelimb and hindlimb grip strength measurement mean", 
                              "Forelimb grip strength normalised against body weight",
                              "Forelimb and hindlimb grip strength normalised against body weight"),
                      new = c("Forelimb and hindlimb strength grip mean", 
                              "Forelimb grip strength (/weight)",
                              "Forelimb and hindlimb grip strength (/weight)"))
phnamsh[match(datswap$old, phnamsh)] <- datswap$new
lineout <- 17
axis(side = 4, labels = NA, at = ats[odds], las = 2, cex.axis = cexax2)
axis(side = 4, labels = NA, at = ats[evens], las = 2, cex.axis = cexax2, tcl = -(lineout - .5), col.ticks = "grey")
mtext(side = 4, line = 1, text = phnamsh[odds], at = ats[odds], las = 2, cex = cexax.labrt, adj = 0)
mtext(side = 4, line = lineout, text = phnamsh[evens], at = ats[evens], las = 2, cex = cexax.labrt, adj = 0)
abline(h = match(procord, procv) - .5, lwd = 2)
# ?mtext

###############################
#Panel (b)
#######################
orth <- 10
ormatpl <- ormat
ormatpl[ormatpl > orth] <- orth
image(x = 1:nfac, y = 1:nfac, log(ormatpl), zlim = c(-1, 1) * log(orth),
      xaxt = "n", yaxt = "n", col = rain, xlab = "", ylab = "")
ypl2 <- grconvertY(c(0, 1), from = "nfc", "ndc")
mtext(side = 1, at = -3, line = 1, text = "(b)", cex = cexlet)

for(j in 1:nrow(zpl)){
  phvc <- colnames(zpl)[order(-abs(zpl[j, ]))[1:10]]
  phnamvc <- pout[match(phvc, pout$ph), "nam"]
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
image(z = t(as.matrix(1:1000)), y = seq(-1, 1, len = 1000), x = 1, col = rain, xaxt = "n", yaxt = "n",
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
image(z = t(as.matrix(1:1000)), y = seq(-1 * log(orth), 1 * log(orth), len = 1000), x = 1, col = rain, xaxt = "n", yaxt = "n",
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
    file.copy(from = paste(factor.plot.dir, "/", fnamc, sep = ""),
              to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)
}
# 
# image(x = 1:nfac, y = 1:nfac, log(ormatpl), zlim = c(-1, 1) * log(orth),
#       xaxt = "n", yaxt = "n", col = rain, xlab = "", ylab = "")
# keybd.press('ctrl+alt+u')





























# facannot <- c("Decreased body size", "Decreased fat mass/insulin, increased anxiety", "Reduced balance/coordination", 
#               "Decreased bone mineral content", "Decreased liver function", "Decreased liver function", 
#               "Decreased grip strength", "Increased activity/exploration 1", "Hearing loss 1",
#               "Increased activity/exploration 2", "Reduced heart function", "Reduced RBC count & increased cell volume",
#               "Increased heart rate", "Reduced hemoglobin", "Hearing loss 2", 
#               "Decreased kidney function", "Increased neutrophil & monocyte differential", "Increased white blood cell count",
#               "Decreased blood cholestrol", "Increased sleep time")



# 
# all.dir <- opp.dir <- same.dir <- matrix(NA, nfac, nfac)
# for(i in 1:(nfac - 1)){
#   for(j in (i + 1):nfac){
#     all.dir[i, j] <- same.dir[i, j] <- (mean(sigmat[, i] == 1 & sigmat[, j] == 1) + mean(sigmat[, i] == -1 & sigmat[, j] == -1)) /
#         (propup[i] * propup[j] + propdo[i] * propdo[j])
#     all.dir[j, i] <- opp.dir[i, j] <- (mean(sigmat[, i] == 1 & sigmat[, j] == -1) + mean(sigmat[, i] == -1 & sigmat[, j] == 1)) /
#       (propup[i] * propdo[j] + propdo[i] * propup[j])
#     # print(c(same.dir, opp.dir))
#   }
# }
# all.dir <- opp.dir <- same.dir <- matrix(NA, nfac, nfac)
# for(i in 1:(nfac - 1)){
#   for(j in (i + 1):nfac){
#     all.dir[i, j] <- all.dir[j, i] <- (mean(sigmat[, i] == 1 & sigmat[, j] == 1) + mean(sigmat[, i] == -1 & sigmat[, j] == -1)
#                                        + mean(sigmat[, i] == 1 & sigmat[, j] == -1) + mean(sigmat[, i] == -1 & sigmat[, j] == 1)) /
#       (propup[i] * propup[j] + propdo[i] * propdo[j] + propup[i] * propdo[j] + propdo[i] * propup[j])
#   }
# }
# 
# image(log(same.dir), col = rain)
# image(log(opp.dir), col = rain)
# image(log(all.dir), col = rain, zlim = c(-1, 1) * log(max(abs(all.dir), na.rm = T)))
# 
# 
# max(all.dir, na.rm = T)
# 
# 
# 
# 
# graphics.off()
# fnamc <- "factor_interpretation_plot.jpg"
# jpeg(paste(factor.plot.dir, "/", fnamc, sep = ""), 13, 7, units = "in", res = 1000)
# wc <- .7
# hc1 <- .15
# hc2 <- .75
# layout(matrix(3:1, 1, 3), wid = c(wc, 1 - wc), hei = c(hc2, 1 - hc2))
# par(mar = c(20, .5, 20, .5), oma = c(0, 15, 0, 6))
# propsig <- colMeans(sigmat[, facord])
# barmat <- t(cbind(propup, propdo))
# cexax <- .75
# # plot(propsig, xaxt = "n", xlab = "", ylab = "", yaxt = "n", ylim = c(0, max(propsig)))
# image(y = 1:nfac, x = 1:nph, z = t(zpl), col = rain, yaxt = "n", ylab = "", cex.axis = cexax, xlab = "", xaxt = "n")
# axis(side = 3, cex.axis = cexax)
# procv <- pout[match(phord, pout$ph), "procnam"]
# procats <- sapply(procord, function(procc) mean(which(procv == procc)))
# axis(side = 1, labels = procord, at = procats, las = 2, cex.axis = .9)
# phnam <- pout[match(phord, pout$ph), "nam"]
# axis(side = 3, labels = phnam, at = 1:length(phord), las = 2, cex.axis = .6)
# abline(v = match(procord, procv) - .5, lwd = 2)
# image(y = 1:nfac, x = 1:nfac, corfac, zlim = c(-1, 1), xaxt = "n", yaxt = "n", col = rain, xlab = "", ylab = "")
# facannot <- c("Decreased body size 1", "Decreased body size 2", "Decreased body size 3", "Decreased liver function",
#               "Decreased fat mass and insulin", "Decreased liver function", "Decreased kidney/liver function", 
#               "Decreased bone mineral content", "Decreased pre-pulse inhibition", "Decreased blood protein", "Hearing loss", 
#               "Decreased blood cholestrol", "Decreased red blood cell volume", "Increased sleep time",
#               "Increased activity", "Increased monocyte/neutrophil content", "Decreased haemaglobin levels",
#               "Increased exploratory behaviour", "Decreased grip strength", "Abnormal ECG",
#               "Increased white blood cell count", "Decreased cardiac function")
# axis(side = 2, labels = facannot, at = 1:nfac, las = 2, cex.axis = 1)
# # axis(side = 1, labels = facannot, at = 1:nfac, las = 2, cex.axis = 1)
# dev.off()
# file.copy(from = paste(factor.plot.dir, "/", fnamc, sep = ""),
#           to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
# 
# 
# 
# jpeg(paste(power.plot.dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
# par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), oma = c(2, 2, 2, 2))
# labs <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
# for(plty in c("non", "imp")){
#   xc <- log(1 + n.uv)
#   yc <- switch(plty, non = log(1 + n.eb.non),
#                imp = log(1 + n.eb.non + n.eb.imp))
#   lims <- c(0, log(1000))#range(c(log(1 + n.uv), log(1 + n.eb.non + n.eb.imp)))
#   repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
#   symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) / 30, inches = F, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, las = 2,
#           xlab = "UV analysis, number of annotations per phenotype (log scale)", ylab = "MV analysis, number of annotations per phenotype (log scale)", cex = .5)
#   # plot(xc, yc, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, cex = .5,
#   #      xlab = "UV analysis, number of annotations per phenotype (log scale)", ylab = "MV analysis, number of annotations per phenotype (log scale)")
#   abline(0, 1)
#   for(ax in 1:2) 
#     axis(side = ax, at = log(1 + labs), labels = labs)
# }
# labs <- c(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
# for(plty in c("non", "imp")){
#   xc <- log(1 + line.n.uv)
#   yc <- switch(plty, non = log(1 + line.n.eb.non),
#                imp = log(1 + line.n.eb.non + line.n.eb.imp))
#   #  lims <- range(c(0, log(1 + line.n.uv), log(1 + line.n.eb.non + line.n.eb.imp)))
#   repdat <- aggregate(data.frame(xc, yc), by = list(xc, yc), length)
#   symbols(repdat$Group.1, repdat$Group.2, circles = sqrt(repdat$xc) / 30, inches = F, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, las = 2, 
#           xlab = "UV analysis, number of annotations per line (log scale)", ylab = "MV analysis, number of annotations per line (log scale)", cex = .5)
#   #plot(xc, yc, xaxt = "n", yaxt = "n", xlim = lims, ylim = lims, xlab = "Number of ")
#   abline(0, 1)
#   for(ax in 1:2) 
#     axis(side = ax, at = log(1 + labs), labels = labs)
# }
# dev.off()
# file.copy(from = paste(power.plot.dir, "/", fnamc, sep = ""),
#           to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
# 
# 
# 
# 
# str(fl)
# mean(is.na(abs(fl$t) > fl$th))
# mean(abs(fl$t) > fl$th, na.rm = T)
# sort(colMeans(sigmat, na.rm = T))
# str(sigmat)
# table(rowSums(sigmat, na.rm = T))
# par(mfrow = c(1, 2))
# hist(colMeans(sigmat, na.rm = T))
# barplot(table(rowSums(sigmat, na.rm = T)))
# 
# facannot <- c("Decreased body size 1", "Decreased body size 2", "Decreased body size 3", "Decreased liver function",
#               "Decreased fat mass and insulin", "Decreased liver function", "Decreased kidney/liver function",
#               "Decreased bone mineral content", "Decreased pre-pulse inhibition", "Decreased blood protein", "Hearing loss",
#               "Decreased blood cholestrol", "Decreased red blood cell volume", "Increased sleep time",
#               "Increased activity", "Increased monocyte/neutrophil content", "Decreased haemaglobin levels",
#               "Increased exploratory behaviour", "Decreased grip strength", "Abnormal ECG",
#               "Increased white blood cell count", "Decreased cardiac function")

# 
#  facdir <- c(-1, -1, 1, 
#              -1, -1, -1, 
#              -1, 1, -1, 
#              -1, 1, -1, )
# facannot <- c("Body size (-)", "Fat mass/insulin (-), Anxiety (+)", "Balance/coordination (-)", 
#               "Bone mineral content (-)", "Liver function (-)", "Liver function (-)", 
#               "Grip strength (-)", "Activity/exploration 1 (+)", "Pre-pulse inhibition (-)",
#               "Activity/exploration 2 (+)", "Heart function (-)", "Red blood cell count (-)",
#               "Heart rate (+)", "Hemoglobin (-)", "ABR threshold (+)", 
#               "Kidney function (-)", "Neutrophil & monocyte differential (+)", "White blood cell count (+)",
#               "Blood cholestrol (-)", "Sleep time (+)")
# facannot <- c("Body size (-)", "Fat mass/insulin (-)", "Bone mineral density (-)", "Liver function (-)", 
#               "Heart function (-)", "Balance/coordination (-)", "Grip strength (-)", "Respiratory exchange ratio (-)",
#               "Heart rate (+)", "Activity/exploration 1 (+)", "Pre-pulse inhibition (-)", "Activity/exploration 2 (+)", 
#               "Blood cholestrol (-)", "ABR threshold (+)", "Red blood cell count (-)", "Hemoglobin (-)", 
#               "Kidney function (-)", "White blood cell count (+)", "Neutrophil & monocyte differential (+)", "Sleep time (+)")  

# sigmatcor <- cor(sigmat)
# 
# sigmatcor <- cor(fl$t, me = "sp")
# prec <- solve(sigmatcor)
# diag(prec) <- NA
# image(prec, col = rain)
# 
# install.packages("glasso")
# library(glasso)
# library(help = glasso)
# library(GGMselect)
# library(help = GGMselect)
# graphout <- selectFast(X = abs(sigmat), family = "LA")
# plot(graphout)
# sum(is.na(fl$mn))
# 
# gl1 <- glasso(s = cor(fl$t[, facord], me = "sp"), rho = .01)
# invpl <- gl1$wi
# diag(invpl) <- NA
# image(invpl, col = rain)
# for(i in 1:(nfac - 1)){
#   for(j in (i + 1):nfac){
#     all.dir[i, j] <- same.dir[i, j] <- (mean(sigmat[, i] == 1 & sigmat[, j] == 1) + mean(sigmat[, i] == -1 & sigmat[, j] == -1)) /
#       (propup[i] * propup[j] + propdo[i] * propdo[j])
#     all.dir[j, i] <- opp.dir[i, j] <- (mean(sigmat[, i] == 1 & sigmat[, j] == -1) + mean(sigmat[, i] == -1 & sigmat[, j] == 1)) /
#       (propup[i] * propdo[j] + propdo[i] * propup[j])
#     # print(c(same.dir, opp.dir))
#   }
# }
# graphics.off()
# fnamc <- "factor_interpretation_plot.jpg"
# # jpeg(paste(factor.plot.dir, "/", fnamc, sep = ""), 7, 11, units = "in", res = 1000)
# wc <- .6
# hc <- .8
# layout(matrix(c(1, 2, 4, 3), 2, 2), hei = c(hc, 1 - hc), wid = c(wc, 1 - wc))
# par(mar = c(0, .5, 2, .5), oma = c(14, 14, 1.5, 11))
# propsig <- colMeans(sigmat[, facord])
# barmat <- 100 * t(cbind(propup, propdo))
# cexax <- .75
# cexax2 <- .45
# cexax3 <- .675
# cexlet <- 1.15
# cexleg <- .6
# # plot(propsig, xaxt = "n", xlab = "", ylab = "", yaxt = "n", ylim = c(0, max(propsig)))
# image(x = 1:nfac, y = 1:nph, z = zpl, col = rain, yaxt = "n", ylab = "", cex.axis = cexax, xaxt = "n")
# mtext(side = 3, at = -2, line = .5, text = "(a)", cex = cexlet)
# abline(v = 1:nfac - .5)
# # axis(side = 3, cex.axis = cexax)
# procv <- pout[match(phord, pout$ph), "procnam"]
# procats <- sapply(procord, function(procc) mean(which(procv == procc)))
# axis(side = 2, labels = procord, at = procats, las = 2, cex.axis = cexax)
# phnam <- pout[match(phord, pout$ph), "nam"]
# axis(side = 4, labels = phnam, at = 1:length(phord), las = 2, cex.axis = cexax2)
# abline(h = match(procord, procv) - .5, lwd = 2)
# par(mar = c(0, .5, .5, .5))
# # image(x = 1:nfac, y = 1:nfac, corfac, zlim = c(-1, 1), xaxt = "n", yaxt = "n", col = rain, xlab = "", ylab = "")
# # image(x = 1:nfac, y = 1:nfac, log(all.dir), zlim = c(-1, 1) * log(max(abs(all.dir), na.rm = T)),
# #       xaxt = "n", yaxt = "n", col = rain, xlab = "", ylab = "")
# orth <- 8
# ormatpl <- ormat
# ormatpl[ormatpl > orth] <- orth
# image(x = 1:nfac, y = 1:nfac, log(ormatpl), zlim = c(-1, 1) * log(orth),
#       xaxt = "n", yaxt = "n", col = rain, xlab = "", ylab = "")
# mtext(side = 1, at = -3, line = 1, text = "(b)", cex = cexlet)
# # facannot <- rep("", nfac)
# 
# for(j in 1:nrow(zpl)){
#   phvc <- colnames(zpl)[order(-abs(zpl[j, ]))[1:10]]
#   phnamvc <- pout[match(phvc, pout$ph), "nam"]
#   datc <- data.frame(nam = phnamvc, val = zpl[j, phvc])
#   print(j)
#   print(datc)
# }
# axis(side = 2, labels = facannot, at = 1:nfac, las = 2, cex.axis = cexax3)
# axis(side = 1, labels = facannot, at = 1:nfac, las = 2, cex.axis = cexax3)
# barplot(barmat, col = c("red", "blue"), xaxs = "i", horiz = TRUE, yaxs = "i",
#         xaxt = "s", xlab = "", ylab = "", las = 1, cex.axis = cexax3, yaxt = "n", las = 2)
# mtext(side = 1, at = max(barmat) * 1.3, line = 1, text = "(c)", cex = cexlet)
# mtext(side = 1, text = "% lines annotated", cex = cexax2, line = 2)
# legend(x = "topright", legend = c("+", "-"), title = "Effect sign", pch = 22, pt.bg = c("red", "blue"), col = 1,
#        cex = cexleg)
# peps <- .2
# par(fig = c(wc * peps, wc * (1 - peps), .98, .99), new = T)
# cexax <- 1.1
# par(mar = c(0, 0, 0, 0))
# image(z = as.matrix(1:1000), x = seq(-1, 1, len = 1000), y = 1, col = rain, xaxt = "n", yaxt = "n", 
#       ylab = "")
# axis(side = 3, las = 2, cex.axis = cexax2, labels = c("-1.0", -0.5, 0.0, 0.5, "1.0"), at = seq(-1, 1, by = .5), las = 0)
# mtext(side = 3, text = expression(italic(tilde(z))), line = 2, cex = cexax, las = 0)
# dev.off()
# file.copy(from = paste(factor.plot.dir, "/", fnamc, sep = ""),
#           to = paste(paper.figures.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
# 
# 
# 
# 
# 
# sort(orpmat, dec = T)[1:10]
# sort(orqmat, dec = T)[1:10]


# 
# all.dir <- matrix(NA, nfac, nfac)
# for(i in 1:(nfac - 1)){
#   for(j in (i + 1):nfac){
#         all.dir[i, j] <- all.dir[j, i] <- 
#           (mean(sigmat[, i] == 1 & sigmat[, j] == 1) + mean(sigmat[, i] == -1 & sigmat[, j] == -1)) /
#           (propup[i] * propup[j] + propdo[i] * propdo[j])
#   }
# }
# all.dir <- matrix(NA, nfac, nfac)
# for(i in 1:(nfac - 1)){
#   for(j in (i + 1):nfac){
#     all.dir[i, j] <- (mean(sigmat[, i] == 1 & sigmat[, j] == 1)) /
#       (propup[i] * propup[j])
#     all.dir[j, i] <- (mean(sigmat[, i] == -1 & sigmat[, j] == -1)) /
#                                              (propdo[i] * propdo[j])
#   }
# }


# 
# write.csv(tab.fac.interp, file = "X:/projects/impc_mv_analysis/data_out/factor_interpretation_table.csv", row.names = F)
# loadmat <- zpl
# save(loadmat, file = "X:/projects/impc_mv_analysis/data_out/loadmat.RData", row.names = F)
