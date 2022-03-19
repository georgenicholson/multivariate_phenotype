resimp <- readRDS(file = control$file.resimp)
n.line.plot <- 500
resimp$procnam <- phmap[match(resimp$ph, phmap$ph), "procnam"]

set.seed(123)
gensub <- c(sample(unique(resimp$geno), n.line.plot))
gensub <- gensub[order(resimp[match(gensub, resimp$geno), "cen"])]
resimp <- resimp[resimp$geno %in% gensub & resimp$ph %in% Data_all$impc$phord, ]
geno.un <- unique(resimp$geno)
ngeno <- length(geno.un)
resimp <- resimp[order(match(resimp$procnam, Data_all$impc$procord), resimp$cen), ]
phen.un <- Data_all$impc$phord
nphen <- length(phen.un)
t.uv.pl <- t.eb.pl <- matrix(NA, ngeno, nphen, dimnames = list(geno.un, phen.un))
t.uv.pl[cbind(resimp$geno, resimp$ph)] <- resimp$uv.t / resimp$uv.th.final
t.eb.pl[cbind(resimp$geno, resimp$ph)] <- resimp$eb.t / resimp$eb.th.final

##########################################################
#Plot heatmaps of t-stats with UV and MV on same page
justSig <- F
zth <- 1
cexax <- .9
incl <- "truekos"#for(incl in c("negcons", "truekos")[2]){
plot_file_name <- "paper_annotation_heatmaps_figure.jpg"

jpeg(file.path(control$figure_dir, plot_file_name), width = 9, height = 9, units = "in", res = 500)
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
  tpl[which(abs(tpl) > zth)] = (sign(tpl) * zth)[which(abs(tpl) > zth)]
  if(justSig)
    tpl[abs(tpl) != zth] <- NA
  if(meth == "uv")
    tpna = is.na(tpl)
  image(x = 1:nrow(tpl), y = 1:ncol(tpl), z = tpl, col = control$heat_col_palette, zlim = c(-1, 1) * zth, 
        xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  tit <- ifelse(meth == "uv", "UV Model", "MV Model")
  mtext(side = 3, line = 0.5, text = tit)
  mtext(side = 3, line = .5, text = ifelse(meth == "uv", "(a)", "(b)"), cex = 1.25, at = 1)
  procv = resimp[match(colnames(tpl), resimp$ph), "procnam"]
  procun = unique(procv)
  ats = sapply(procun, function(x) mean(which(procv == x)))
  linats = sapply(procun, function(x) match(x, procv) - .5)
  axis(side = 2, at = ats, labels = procun, las = 2, cex.axis = .7)
  abline(h = linats)
  cenv = Data_all$impc$cenmap[match(resimp[match(rownames(tpl), resimp$geno), "cen"], Data_all$impc$cenmap$cen), "nam"]
  cenats = sapply(Data_all$impc$cenmap$nam, function(x) mean(which(cenv == x)))
  cenlinats = sapply(Data_all$impc$cenmap$nam, function(x) match(x, cenv) - .5)
  axis(side = 1, at = cenats, labels = Data_all$impc$cenmap$nam, las = 2, cex.axis = .7)
  if(meth == "eb")
    mtext(side = 1, text = "Phenotyping laboratory", line = 3, cex = .7, las = 1)
  mtext(side = 2, text = "Procedure", line = 11, cex = .7)
  abline(v = cenlinats)
}
par(mar = c(2, 1, 2, 1), xpd = F)
image(z = t(as.matrix(1:1000)), x = 1, y = seq(-1, 1, len = 1000), col = control$heat_col_palette, xaxt = "n", yaxt = "n", 
      ylab = "")
axis(side = 4, las = 2, cex.axis = cexax, labels = c("< -1.0", -0.5, 0.0, 0.5, "> 1.0"), at = seq(-1, 1, by = .5))
mtext(side = 4, text = expression(italic(tilde(z))), line = 3, cex = cexax, las = 2)
dev.off()
file.copy(from = file.path(control$figure_dir, "paper_annotation_heatmaps_figure.jpg"),
          to = file.path(control$dropbox_figure_dir, "paper_annotation_heatmaps_figure.jpg"), overwrite = TRUE)

##########################################################
#Plot Sig and R on same plot
Sig.mn <- compl[[control$mv_meth_nam_use]]$Sig.comb
R.mn <- compl[[control$mv_meth_nam_use]]$R.comb
SigCormn <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
fnamc <- "paper_correlation_heatmaps_figure.jpg"#paste("estmeth_", estimation.meth, "_EMits_", em.nits, "_Sig_and_R_heatmaps.jpg", sep = "")
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), width = 12, height = 8, units = "in", res = 500)
par(mfrow = c(1, 2))
layout(mat = matrix(c(1, 1, 2, 2, 3, 4), 2, 3), widths = c(.47, .47, .06), heights = c(.25, .75))
par(oma = c(4, 21, 4, 4))
zpl <- SigCormn[phen.un, phen.un]
# SigCormn <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
for(plty in c("Sig", "R")){
  if(plty == "Sig"){
    zpl <- SigCormn[Data_all$impc$phord, rev(Data_all$impc$phord)]
    par(mar = c(22, 1, 1, 1), xpd = F)
  }
  if(plty == "R"){
    zpl <- R.mn[Data_all$impc$phord, rev(Data_all$impc$phord)]
    par(mar = c(22, 1, 1, 1), xpd = F)
  }
  cexax <- 1.15
  image(x = 1:nrow(zpl), y = 1:ncol(zpl), z = zpl, zlim = c(-1, 1), col = control$heat_col_palette, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
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
image(z = t(as.matrix(1:1000)), x = 1, y = seq(-1, 1, len = 1000), col = control$heat_col_palette, xaxt = "n", yaxt = "n", ylab = "t")
axis(side = 4, las = 2, cex.axis = cexax)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
