data <- "impc"
source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")
corr.mat.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/correlation_matrices_revision"
dir.create(corr.mat.plot.dir, showWarnings = F)
annot.heatmaps.plot.dir <- "X:/projects/impc_mv_analysis/plots/impc_mv_paper_plots/annotation_heatmaps_revision"
dir.create(annot.heatmaps.plot.dir, showWarnings = F)

file.glob.res <- paste0(global.res.dir, "/global_eb_results.RData")
load(uv.results.Y.S)
load(file = file.glob.res)
# load(file = em.curated.results.file)
load(file = paste(processed.data.dir, "/", date.data.arrived, "_impc_data_in.RData", sep = ""))
d <- d[order(d$day), ]
geno.un <- unique(resimp$geno)
ngeno <- length(geno.un)
first.day <- d[match(geno.un, d$geno), "day"]
names(first.day) <- geno.un
n.line.plot <- 500

# load(file = resimp.with.sig.thresholds.file)
resimp$geno.first.day <- first.day[match(resimp$geno, names(first.day))]
resimp$procnam <- phmap[match(resimp$ph, phmap$ph), "procnam"]
  
#############################################################
#create and save ordering of procedures for plotting, so that labels don't overlap
###################################
procall <- unique(resimp$procnam)
phenall <- unique(resimp$ph)
proctab <- table(pout[match(phenall, pout$ph), "procnam"])
proctabord <- sort(proctab)
procord <- c()
smallbig <- 0
while(length(proctabord) > 0){
  if(smallbig == 0){
    procord <- c(procord, names(proctabord)[1])
    proctabord <- proctabord[-1]
  }
  if(smallbig == 1){
    procord <- c(procord, names(proctabord)[length(proctabord)])
    proctabord <- proctabord[-length(proctabord)]
  }
  smallbig <- 1 - smallbig
}
inds.swap <- which(procord %in% c("Intraperitoneal glucose tolerance test (IPGTT)", "Heart Weight"))
procord[inds.swap] <- procord[rev(inds.swap)]
procord[]
phord <- unique(resimp[order(match(resimp$procnam, procord), resimp$ph), "ph"])
phord <- phord[!phord %in% facnam]
save(procord, phord, file = procord.file)

set.seed(123)
gensub <- c(sample(unique(resimp$geno[resimp$line.type == "trueMut"]), n.line.plot),
            sample(unique(resimp$geno[resimp$line.type == "negCon"]), n.line.plot))
gensub <- gensub[order(resimp[match(gensub, resimp$geno), "cen"])]
resimp <- resimp[resimp$geno %in% gensub & resimp$ph %in% phord, ]
geno.un <- unique(resimp$geno)
ngeno <- length(geno.un)

# resimp <- resimp[sample(1:nrow(resimp)), ]
resimp <- resimp[order(match(resimp$procnam, procord), resimp$cen), ]
# resimp <- resimp[order(match(resimp$procnam, procord), resimp$cen, resimp$geno.first.day, resimp$geno), ]
# resimp <- resimp[order(match(resimp$procnam, procord), resimp$cen, match(resimp$geno, gensub)), ]
phen.un <- phord#unique(resimp$ph)
nphen <- length(phen.un)
t.uv.pl <- t.eb.pl <- matrix(NA, ngeno, nphen, dimnames = list(geno.un, phen.un))
t.uv.pl[cbind(resimp$geno, resimp$ph)] <- resimp$uv.t / resimp$uv.th.final
t.eb.pl[cbind(resimp$geno, resimp$ph)] <- resimp$eb.t / resimp$eb.th.final

# str(t.uv.pl[cbind(resimp$geno, resimp$ph)] )
# str(t.uv.pl)
# mean(resimp$geno %in% rownames(t.uv.pl))
# mean(resimp$ph %in% colnames(t.uv.pl))
# length(unique(resimp$ph))


##########################################################
#Plot heatmaps of t-stats with UV and MV on same page
for(justSig in c(F, T)){
  zth <- 1
  cexax <- .9
  for(incl in c("negcons", "truekos")){
    fnamc <- paste("estmeth_", estimation.meth, "_EMits_", em.nits, "_annotation_heatmap_justSig_", 
                   justSig, "_compuvmv_", incl, ".jpg", sep = "")
    plnamc <- paste(annot.heatmaps.plot.dir, "/", fnamc, sep = "")
    # pdf(plnamc, 12, 9)
    jpeg(paste(annot.heatmaps.plot.dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
    par(oma = c(4, 4, 4, 4))
    layout(mat = matrix(c(1, 1, 2, 2, 3, 5, 4, 6), 4, 2), widths = c(.95, .05), heights = rep(.25, 4))
    par(mar = c(2, 12, 2, 1))#, oma = c(4, 1, 4, 1))
    for(meth in c("uv", "eb")){
      if(meth == "uv")
        tpl = t.uv.pl
      if(meth == "eb")
        tpl = t.eb.pl
      if(incl == "negcons")
        tpl = tpl[resimp[match(rownames(tpl), resimp$geno), "line.type"] == "negCon", ]
      if(incl == "truekos")
        tpl = tpl[resimp[match(rownames(tpl), resimp$geno), "line.type"] == "trueMut", ]
      tpl = tpl[, colnames(tpl) %in% resimp$ph]
      tpl = tpl[order(match(rownames(tpl), resimp$geno)), order(match(colnames(tpl), resimp$ph))]
      # tpl = tpl[order(match(rownames(tpl), gensub)), order(match(colnames(tpl), resimp$ph))]
      tpl[which(abs(tpl) > zth)] = (sign(tpl) * zth)[which(abs(tpl) > zth)]
      if(justSig)
        tpl[abs(tpl) != zth] <- NA
      if(meth == "uv")
        tpna = is.na(tpl)
      image(x = 1:nrow(tpl), y = 1:ncol(tpl), z = tpl, col = rain, zlim = c(-1, 1) * zth, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
      # tit <- paste(ifelse(meth == "uv", "UV analysis on ", "MV analysis on "), 
      #              ifelse(incl == "negcons", "negative controls", "KO lines"), sep = "")
      tit <- ifelse(meth == "uv", "UV Model", "MV Model")
      mtext(side = 3, line = 0.5, text = tit)
      mtext(side = 3, line = .5, text = ifelse(meth == "uv", "(a)", "(b)"), cex = 1.25, at = 1)
      procv = resimp[match(colnames(tpl), resimp$ph), "procnam"]
      procun = unique(procv)
      ats = sapply(procun, function(x) mean(which(procv == x)))
      linats = sapply(procun, function(x) match(x, procv) - .5)
      axis(side = 2, at = ats, labels = procun, las = 2, cex.axis = .7)
      abline(h = linats)
      cenv = cenmap[match(resimp[match(rownames(tpl), resimp$geno), "cen"], cenmap$cen), "nam"]
      cenats = sapply(cennamun, function(x) mean(which(cenv == x)))
      cenlinats = sapply(cennamun, function(x) match(x, cenv) - .5)
      axis(side = 1, at = cenats, labels = cennamun, las = 2, cex.axis = .7)
      if(meth == "eb")
        mtext(side = 1, text = "Phenotyping laboratory", line = 3, cex = .7, las = 1)
      mtext(side = 2, text = "Procedure", line = 11, cex = .7)
      abline(v = cenlinats)
      # par(xpd = NA)
      # legend(x = nrow(tpl) * 1.025, y = ncol(tpl) * .975, yjust = 1, xjust = 0, col = rain[c(1000, 500, 1)], 
      #        legend = c(3, 0, -3), pch = 15, title = "t", cex = .9)
      # par(xpd = F)
    }
    par(mar = c(2, 1, 2, 1), xpd = F)
    image(z = t(as.matrix(1:1000)), x = 1, y = seq(-1, 1, len = 1000), col = rain, xaxt = "n", yaxt = "n", 
          ylab = "")
    axis(side = 4, las = 2, cex.axis = cexax, labels = c("< -1.0", -0.5, 0.0, 0.5, "> 1.0"), at = seq(-1, 1, by = .5))
    mtext(side = 4, text = expression(italic(tilde(z))), line = 3, cex = cexax, las = 2)
    dev.off()
    # file.copy(from = paste(annot.heatmaps.plot.dir, "/", fnamc, sep = ""),
    #           to = paste(chris.pres.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
    if(incl == "truekos" & meth == "eb" & justSig == FALSE)
      file.copy(from = paste(annot.heatmaps.plot.dir, "/", fnamc, sep = ""),
              to = paste(revision.paper.figures, "/paper_annotation_heatmaps_figure.jpg", sep = ""), overwrite = TRUE)
  }
}


##########################################################
#Plot heatmaps of t-stats
for(justSig in c(F, T)){
  zth <- 1
  cexout <- 1.5
  cexax <- .9
  for(meth in c("uv", "eb")){
    for(incl in c("negcons", "truekos")){
      fnamc <- paste("estmeth_", estimation.meth, "_EMits_", em.nits, "_annotation_heatmap_justSig_", 
                     justSig, "_", meth, "_", incl, ".pdf", sep = "")
      plnamc <- paste(annot.heatmaps.plot.dir, "/", fnamc, sep = "")
      pdf(plnamc, 12, 8)
      par(oma = c(4, 4, 4, 4))
      layout(mat = matrix(c(1, 1, 2, 3), 2, 2), widths = c(.95, .05), heights = c(.25, .75))
      par(mar = c(6, 19, 1, 1))#, oma = c(4, 1, 4, 1))
      # par(mfrow = c(1, 1), mar = c(6, 19, 3, 4))
      if(meth == "uv")
        tpl = t.uv.pl
      if(meth == "eb")
        tpl = t.eb.pl
      if(incl == "negcons")
        tpl = tpl[resimp[match(rownames(tpl), resimp$geno), "line.type"] == "negCon", ]
      if(incl == "truekos")
        tpl = tpl[resimp[match(rownames(tpl), resimp$geno), "line.type"] == "trueMut", ]
      # tpl = tpl[!resimp[match(rownames(tpl), resimp$geno), "negcon"], ]
      tpl = tpl[, colnames(tpl) %in% resimp$ph]
      tpl = tpl[order(match(rownames(tpl), resimp$geno)), order(match(colnames(tpl), resimp$ph))]
      tpl[which(abs(tpl) > zth)] = (sign(tpl) * zth)[which(abs(tpl) > zth)]
      if(justSig)
        tpl[abs(tpl) != zth] <- NA
      if(meth == "uv")
        tpna = is.na(tpl)
      image(x = 1:nrow(tpl), y = 1:ncol(tpl), z = tpl, col = rain, zlim = c(-1, 1) * zth, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
      tit <- paste(ifelse(meth == "uv", "UV analysis on ", "MV analysis on "), 
                   ifelse(incl == "negcons", "negative controls", "KO lines"), sep = "")
      # mtext(side = 3, line = 0, text = tit)
      procv = resimp[match(colnames(tpl), resimp$ph), "procnam"]
      procun = unique(procv)
      ats = sapply(procun, function(x) mean(which(procv == x)))
      linats = sapply(procun, function(x) match(x, procv) - .5)
      axis(side = 2, at = ats, labels = procun, las = 2, cex.axis = cexax)
      abline(h = linats)
      cenv = cenmap[match(resimp[match(rownames(tpl), resimp$geno), "cen"], cenmap$cen), "nam"]
      cenats = sapply(cennamun, function(x) mean(which(cenv == x)))
      cenlinats = sapply(cennamun, function(x) match(x, cenv) - .5)
      axis(side = 1, at = cenats, labels = cennamun, las = 2, cex.axis = cexax)
      abline(v = cenlinats)
      mtext(side = 1, line = 3.3, text = "KO genes (grouped by phenotyping centre)", cex = cexout)
      mtext(side = 2, line = 17, text = "Phenotypes (grouped by procedure)", cex = cexout)
      # par(xpd = NA)
      # legend(x = nrow(tpl) * 1.025, y = ncol(tpl) * .975, yjust = 1, xjust = 0, col = rain[c(1000, 500, 1)], 
      #        legend = c(3, 0, -3), pch = 15, title = "t")
      # par(xpd = F)
      par(mar = c(1, 1, 1, 1), xpd = F)
      image(z = t(as.matrix(1:1000)), x = 1, y = seq(-3, 3, len = 1000), col = rain, xaxt = "n", yaxt = "n", ylab = "t")
      axis(side = 4, las = 2, cex.axis = cexax)
      mtext(side = 4, text = expression(italic(t)), line = 3, cex = cexax, las = 2)
      dev.off()
      # file.copy(from = paste(annot.heatmaps.plot.dir, "/", fnamc, sep = ""),
      #           to = paste(chris.pres.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
    }
  }
}
  

##########################################################
#Plot correlation matrix
# keepinds <- (mv.burnin + 1):mv.nits
# Sig.mn <- (Reduce('+', Sigl[keepinds]) / nkeep)[phen.un, phen.un]
# R.mn <- (Reduce('+', Rl.collect[keepinds]) / nkeep)[phen.un, phen.un]

Sig.mn <- resl.out$eb$Sig.comb
R.mn <- resl.out$eb$R
SigCormn <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
str(SigCormn)
for(plty in c("Sig", "R")){
  if(plty == "Sig"){
    fnamc <- paste("estmeth_", estimation.meth, "_EMits_", em.nits, "_Sig_correl_matrix.pdf", sep = "")
    pdf(paste(corr.mat.plot.dir, "/", fnamc, sep = ""), 12, 12)
    zpl <- SigCormn[phord, phord]
  }
  if(plty == "R"){
    fnamc <- paste("estmeth_", estimation.meth, "_EMits_", em.nits, "_R_matrix.pdf", sep = "")
    pdf(paste(corr.mat.plot.dir, "/", fnamc, sep = ""), 12, 12)
    zpl <- R.mn[phord, phord]
  }
  par(oma = c(4, 4, 4, 4))
  layout(mat = matrix(c(1, 1, 2, 3), 2, 2), widths = c(.95, .05), heights = c(.25, .75))
  par(mar = c(22, 22, 1, 1), xpd = F)
  cexax <- 1.15
  image(x = 1:nrow(zpl), y = 1:ncol(zpl), z = zpl, zlim = c(-1, 1), col = rain, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  str(zpl)
  procv <- resimp[match(colnames(zpl), resimp$ph), "procnam"]
  procun <- unique(procv)
  ats <- sapply(procun, function(x) mean(which(procv == x)))
  linats <- sapply(procun, function(x) match(x, procv) - .5)
  axis(side = 2, at = ats, labels = procun, las = 2, cex.axis = cexax)
  axis(side = 1, at = ats, labels = procun, las = 2, cex.axis = cexax)
  abline(h = linats, v = linats)
  par(mar = c(1, 1, 1, 1), xpd = F)
  image(z = t(as.matrix(1:1000)), x = 1, y = seq(-1, 1, len = 1000), col = rain, xaxt = "n", yaxt = "n", ylab = "t")
  axis(side = 4, las = 2, cex.axis = cexax)
  # mtext(side = 4, text = "r", line = 3, cex = cexax, las = 2)
  # 
  # c(x1 = nrow(zpl) * 1.025, y1 = ncol(zpl) * .9, x2 = nrow(zpl) * 1.05, y2 = ncol(zpl) * .975)
  # par(new = T, usr = c(nrow(zpl) * 1.025, ncol(zpl) * .9, nrow(zpl) * 1.05, ncol(zpl) * .975))
  # image(z = as.matrix(1:1000), x = 1:1000, y = 1, col = rain)
  # legend(x = nrow(zpl) * 1.025, y = ncol(zpl) * .975, yjust = 1, xjust = 0, col = rain[c(1000, 500, 1)], 
  #        legend = c(1, 0, -1), pch = 15, title = "cor", cex = 1.2)
  # par(xpd = NA)
  # legend(x = nrow(zpl) * 1.025, y = ncol(zpl) * .975, yjust = 1, xjust = 0, col = rain[c(1000, 500, 1)], 
  #        legend = c(1, 0, -1), pch = 15, title = "cor", cex = 1.2)
  # par(xpd = F)
  dev.off()
  # file.copy(from = paste(corr.mat.plot.dir, "/", fnamc, sep = ""),
  #           to = paste(chris.pres.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)
}


##########################################################
#Plot Sig and R on same plot
fnamc <- paste("estmeth_", estimation.meth, "_EMits_", em.nits, "_Sig_and_R_heatmaps.jpg", sep = "")
jpeg(paste(corr.mat.plot.dir, "/", fnamc, sep = ""), 12, 8, units = "in", res = 500)
par(mfrow = c(1, 2))
layout(mat = matrix(c(1, 1, 2, 2, 3, 4), 2, 3), widths = c(.47, .47, .06), heights = c(.25, .75))
par(oma = c(4, 21, 4, 4))
zpl <- SigCormn[phen.un, phen.un]
# SigCormn <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
for(plty in c("Sig", "R")){
  if(plty == "Sig"){
    zpl <- SigCormn[phord, rev(phord)]
    par(mar = c(22, 1, 1, 1), xpd = F)
  }
  if(plty == "R"){
    zpl <- R.mn[phord, rev(phord)]
    par(mar = c(22, 1, 1, 1), xpd = F)
  }
  cexax <- 1.15
  image(x = 1:nrow(zpl), y = 1:ncol(zpl), z = zpl, zlim = c(-1, 1), col = rain, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  str(zpl)
  for(sidec in 1:2){
    if(sidec == 2)
      procv <- resimp[match(colnames(zpl), resimp$ph), "procnam"]
    if(sidec == 1)
      procv <- resimp[match(rownames(zpl), resimp$ph), "procnam"]
    
    procun <- unique(procv)
    ats <- sapply(procun, function(x) mean(which(procv == x)))
    linats <- sapply(procun, function(x) match(x, procv) - .5)
    if(sidec == 2){
      if(plty == "Sig")
        axis(side = 2, at = ats, labels = procun, las = 2, cex.axis = cexax)
      abline(h = linats)
    }
    if(sidec == 1){
      axis(side = 1, at = ats, labels = procun, las = 2, cex.axis = cexax)
      abline(v = linats)
    }
  }
  mtext(side = 3, line = .5, text = ifelse(plty == "Sig", "(a)", "(b)"), cex = 1.25, at = 1)
  if(plty == "Sig")
    mtext(side = 3, line = 1, text = expression(paste("Knockout-induced correlation, ", italic(hat(Sigma))), sep = ""), cex = 1.3)
  if(plty == "R")
    mtext(side = 3, line = 1, text = expression(paste("Experimental correlation, ", italic(hat(R))), sep = ""), cex = 1.3)
}
par(mar = c(1, 1, 1, 1), xpd = F)
image(z = t(as.matrix(1:1000)), x = 1, y = seq(-1, 1, len = 1000), col = rain, xaxt = "n", yaxt = "n", ylab = "t")
axis(side = 4, las = 2, cex.axis = cexax)
dev.off()
file.copy(from = paste(corr.mat.plot.dir, "/", fnamc, sep = ""),
          to = paste(revision.paper.figures, "/paper_correlation_heatmaps_figure.jpg", sep = ""), overwrite = TRUE)
# file.copy(from = paste(corr.mat.plot.dir, "/", fnamc, sep = ""),
#           to = paste(chris.pres.dropbox, "/", fnamc, sep = ""), overwrite = TRUE)


# 
# pout[grep("IMM", pout$ph), ]
# 
# pout[grep("Eye", pout$procnam), ]
# 
# 
# pl.look <- pout[grep("Pleth", pout$procnam), "ph"]
# plc <- pl.look[6]
# par(mfrow = c(2, 2))
# hist(resimp[resimp$ph == plc & resimp$imputed, "eb.mn"])
# hist(resimp[resimp$ph == plc & !resimp$imputed, "eb.mn"])
# hist(resimp[resimp$ph == plc & resimp$imputed, "eb.t"])
# hist(resimp[resimp$ph == plc & !resimp$imputed, "eb.t"])
# 
# 
# hist(resimp[grepl("Pleth", resimp$procnam) & !resimp$imputed, "eb.mn"])
# 
# 
# ?grepl
# 


# #for(plot.type in c("cen.fpr", "cen.fdr", "phcen.interaction.fpr", "phcen.interaction.fdr")){
# # t.uv.pl = uvl$t / uvl[[plot.type]]$thmat# resimp[match(rownames(uvl$t), resimp$geno), "th.uv"]
# # t.mv.pl = ebl$t / mvl[[plot.type]]$thmat
# 
# 
# 
#   
#       if(meth == "eb"){
#         #par(new = T)
#         #image(x = 1:nrow(tpl), y = 1:ncol(tpl), z = tpna, col = c(grey(.9, .2), grey(.8, 0)), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#       }
#       # for(j in 1:nrow(tpl)){
#       #   for(k in 1:ncol(tpl)){
#       #     if(any(na.omit(abs(tpl[j, k])) == 1)){
#       #       lines(x = j + c(-.5, .5), y = k + c(-.5, .5))
#       #       lines(x = j + c(-.5, .5), y = k + c(.5, -.5))
#       #     }
#       #   }
#       # }
#     }
#   }
#   # dev.off()
# }
# if(estimation.meth == "mcmc"){
#   load(mcmcfile)
#   load(file = mcmc.curated.results.file)
#   mumat <- sapply(mul, function(v) v)
#   mu.tmn <- sort(rowMeans(mumat) / apply(mumat, 1, function(v) sd(v) / sqrt(ncol(mumat))))
#   sort(sapply(procun, function(prc) mean(mu.tmn[pout[pout$procnam == prc, "ph"]], na.rm = T)))
#   keepinds.plot <- keepinds#4500:5000
#   nkeep.plot <- length(keepinds.plot)
#   Sig.mn <- (Reduce('+', Sigl[keepinds.plot]) / nkeep.plot)
#   R.mn <- (Reduce('+', Rl.collect[keepinds.plot]) / nkeep.plot)
# }
# if(estimation.meth == "map"){
  