resimp <- readRDS(file = control$file.resimp)
####################################################
#Reference lines numbers for text
resimp0.imp <- resimp[resimp$line.type == control$control$nam.negcon & resimp$imputed, ]
resimp1.imp <- resimp[resimp$line.type == control$control$nam.truemut & resimp$imputed, ]
glob.fdr.imp <- mean(abs(resimp0.imp$eb.t) > resimp0.imp$eb.th.final) / 
  mean(abs(resimp1.imp$eb.t) > resimp1.imp$eb.th.final)
resimp0.non <- resimp[resimp$line.type == control$nam.negcon & !resimp$imputed, ]
resimp1.non <- resimp[resimp$line.type == control$nam.truemut & !resimp$imputed, ]
glob.fdr.non <- mean(abs(resimp0.non$eb.t) > resimp0.non$eb.th.final) / 
  mean(abs(resimp1.non$eb.t) > resimp1.non$eb.th.final)

matl <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$matl
matl$pre <- resll$perm$uv$ref.lines$matl$post
str(matl)
matl.sig <- lapply(matl, function(M){M[abs(M) < 1] <- 0; sign(M)})
matl.sig.ord <- lapply(matl.sig, function(M)t(apply(M, 1, function(v) sort(v, decreasing = T))))
tabl <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$tabl

fdr.ests <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$fdr.ests
fdr.ests$pre <- resll$perm$uv$ref.lines$fdr.ests$post


#FDR ests and CIs
for(i in 1:length(fdr.ests)){
  fdrc <- fdr.ests[[i]]
  loomv.uv.ci.numv <- formatC(unlist(fdrc) * 100, format = "f", digits = 1)
  loomv.uv.ci <- paste(loomv.uv.ci.numv[1], "\\% (95\\% CI: ", loomv.uv.ci.numv[2], "\\% - ", loomv.uv.ci.numv[3], "\\%)", sep = "")
  save.text <- paste("reflines.fdr.ci.", names(fdr.ests)[i], sep = "")
  for(numc in save.text)
    write.table(loomv.uv.ci, file = paste(control$dropbox_text_numbers_dir, "/", save.text, ".txt", sep = ""),
                col.names = F, row.names = F, quote = F)
}
#global concordance/discordance
ref.lines.concordant.n <- tabl$post["-1", "-1"] + tabl$post["1", "1"]
ref.lines.discordant.n <- tabl$post["-1", "1"] + tabl$post["1", "-1"]
save.num <- c("ref.lines.concordant.n", "ref.lines.discordant.n")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


####################################################
#Reference lines concordance scatter
fnamc <- "ref_lines_z_scatter.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
par(mfrow = c(2, 2), oma = c(3, 2, 1, 1), mar = c(3, 3, 3, 3))
namv <- c("pre", "postcomp", "postimp", "postimpcomp")
for(namc in namv){
  ntab <- table(matl.sig[[namc]][, 1], matl.sig[[namc]][, 2])
  ptab <- ntab / sum(ntab)
  limn <- 3
  legcex <- .65
  cexax <- .9
  plot(matl[[namc]][, 1], matl[[namc]][, 2], pch = 1, cex = .2, 
       xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
  abline(0, 1)
  labc1 <- switch(namc, pre = "UV replicate 1 (measured)",
                  postcomp = "MV replicate 1 (measured)",
                  postimp = "MV replicate 1 (missing)",
                  postimpcomp = "MV replicate 1 (missing)")
  labc2 <- switch(namc, pre = "UV replicate 2 (measured), ",
                  postcomp = "MV replicate 2 (measured)",
                  postimp = "MV replicate 2 (missing)",
                  postimpcomp = "MV replicate 2 (measured)")
  mtext(side = 1, line = 2.5, text = substitute(paste(labc1, ", ", italic(tilde(z))), list(labc1 = labc1)), cex = cexax)
  mtext(side = 2, line = 2.5, text = substitute(paste(labc2, ", ", italic(tilde(z))), list(labc2 = labc2)), cex = cexax)
  abline(h = c(-1, 1), col = 2, lty = 2)
  abline(v = c(-1, 1), col = 2, lty = 2)
  atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
  xl <- list(c(1, 1, 3, 3), c(1, 1, 3, 3) * -1, c(1, 1, 3, 3), c(1, 1, 3, 3) * -1)
  yl <- list(c(1, 3, 3, 1), c(1, 3, 3, 1) * -1, c(1, 3, 3, 1) * -1, c(1, 3, 3, 1))
  for(quad in 1:4){
    if(quad %in% 1:2)
      colc <- rgb(0, 0, 1, .35)
    if(quad %in% 3:4)
      colc <- rgb(1, 0, 0, .35)
    # adjustcolor("grey", .35)
    polygon(x = xl[[quad]], y = yl[[quad]], col = colc)
  }
  for(i in 1:3){
    for(j in 1:3){
      # labc <- paste(formatC(100 * ptab[i, j], format = "f", digits = 1), "% (", prettyNum(ntab[i, j], big.mark = ","), ")", sep = "")
      labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
      legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .75), xjust = .5, cex = legcex)
    }
  }
  mtext(side = 3, text = paste("(", letters[match(namc, namv)], ")", sep = ""), at = -limn - .25, line = 1, cex = 1.2)
  numv <- formatC(c(fdr.ests[[namc]][[1]], fdr.ests[[namc]][[2]]) * 100, format = "f", digits = 1)
  methlab <- switch(namc, pre = "UV (measured)", postcomp = "MV (measured)", postimp = "MV (missing)", postimpcomp = "MV (meas. vs miss.)")
  fdrlab <- paste0(methlab, " Fsr = ", numv[1], "% (", numv[2], "% - ", numv[3], "%)")
  mtext(side = 3, text = fdrlab, line = .5, cex = .9)
  # mtext(side = 3, text = fdrlab, line = 1.5, cex = .75)
}
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)



#######################################
#Heatmap of reference lines results
graphics.off()
tl.pre <- resll$perm$uv$ref.lines$tl
thl.pre <- resll$perm$uv$ref.lines$thl
tl.post <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$tl
thl.post <- resll$perm[[control$mv_meth_nam_use]]$ref.lines$thl

str(tl.pre)
nlin <- length(tl.pre)
nlinrep <- sapply(tl.pre, ncol)
zth = 1
names(tl.pre)
big.mar <- 1.5
small.mar <- .5
fnamc <- "ref_lines_heatmap.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 12, 18, units = "in", res = 600)
layout(matrix(1:(2 * nlin), 1), width = rep(nlinrep, each = 2) + big.mar)
par(oma = c(6, 22, 0, 20))
martop <- 10
for(i in 1:length(tl.pre)){
  for(ty in c("pre", "pos")){
    cenv <- resimp[match(colnames(tl.pre[[i]]), resimp$geno), "cenlong"]
    nlinrepc <- nlinrep[i]
    if(ty == "pre"){
      plc = (tl.pre[[i]] / thl.pre[[i]])[Data_all$impc$phord, ]
      par(mar = c(1, big.mar / 2, martop, small.mar / 2))
    }
    if(ty == "pos"){
      plc = (tl.post[[i]] / thl.post[[i]])[Data_all$impc$phord, ]
      par(mar = c(1, small.mar / 2, martop, big.mar / 2))
    }
    plc[which(abs(plc) > zth)] = zth * sign(plc[which(abs(plc) > zth)])
    image(x = 1:ncol(plc), y = 1:nrow(plc), z = t(plc), col = control$heat_col_palette, zlim = c(-1, 1) * zth, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for(j in 1:nrow(plc)){
      for(k in 1:ncol(plc)){
        if(any(na.omit(abs(plc[j, k])) == 1)){
          lines(x = k + c(-.5, .5), y = j + c(-.5, .5))
          lines(x = k + c(.5, -.5), y = j + c(-.5, .5))
        }
      }
    }
    procv <- Data_all$impc$phmap[match(Data_all$impc$phord, Data_all$impc$phmap$ph), "procnam"]
    if(i == 1 & ty == "pre"){
      procats <- sapply(Data_all$impc$procord, function(procc) mean(which(procv == procc)))
      axis(side = 2, labels = Data_all$impc$procord, at = procats, las = 2, cex.axis = 1.2)
    }
    if(i == nlin & ty == "pos"){
      phnam <- Data_all$impc$phmap[match(Data_all$impc$phord, Data_all$impc$phmap$ph), "nam"]
      axis(side = 4, labels = phnam, at = 1:length(Data_all$impc$phord), las = 2, cex.axis = 1)
    }
    abline(h = match(Data_all$impc$procord, procv) - .5, lwd = 2)
    polygon(x = rep(c(0, nlinrepc), each = 2) + .5, y = c(0, nrow(plc))[c(1, 2, 2, 1)] + .5, lwd = 2)
    mtext(side = 3, text = ifelse(ty == "pre", "UV", "MV"), cex = .8)
    axis(side = 1, labels = cenv, at = 1:ncol(plc), cex.axis = .9, las = 2)
    splc <- strsplit(names(tl.pre)[i], spl = "_")[[1]]
    gnam1 <- quote(substitute(italic(a), list(a = splc)))
    gnam2 <- quote(substitute(paste("(", b, ")", sep = ""), list(b = ifelse(splc[2] == 1, "hom", "het"))))
    if(ty == "pre"){
      mtext(side = 3, at = nlinrepc + 1, line = 2.5, text = eval(gnam1), cex = 1)
      mtext(side = 3, at = nlinrepc + 1, line = 1, text = eval(gnam2), cex = 1)
    }
  }
}
par(fig = c(.3, .7, .96, .97), new = T)
cexax <- 1.1
par(mar = c(0, 0, 0, 0))
image(z = as.matrix(1:1000), x = seq(-1, 1, len = 1000), y = 1, col = control$heat_col_palette, xaxt = "n", yaxt = "n", 
      ylab = "")
axis(side = 3, las = 2, cex.axis = cexax, labels = c("< -1.0", -0.5, 0.0, 0.5, "> 1.0"), at = seq(-1, 1, by = .5), las = 0)
mtext(side = 3, text = expression(italic(tilde(z))), line = 2, cex = cexax, las = 0)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)



#

hom.z <- (resll$perm[[control$mv_meth_nam_use]]$het.hom$tl$het.hom / resll$perm[[control$mv_meth_nam_use]]$het.hom$thl$het.hom)[, 1]
het.z <- (resll$perm[[control$mv_meth_nam_use]]$het.hom$tl$het.hom / resll$perm[[control$mv_meth_nam_use]]$het.hom$thl$het.hom)[, 2]
hom.z.sign <- sign(hom.z) * (abs(hom.z) > 1)
het.z.sign <- sign(het.z) * (abs(het.z) > 1)
tabhethom <- table(hom.z.sign, het.z.sign)

###############################################
#Het hom FDR ests and CIs
tabhethom <- table(hom.z.sign, het.z.sign)
fdrci <- fdr.est.tab(tabhethom)
fdrc <- fdrci
loomv.uv.ci.numv <- formatC(unlist(fdrc) * 100, format = "f", digits = 1)
loomv.uv.ci <- paste(loomv.uv.ci.numv[1], "\\% (95\\% CI: ", loomv.uv.ci.numv[2], "\\% - ", loomv.uv.ci.numv[3], "\\%)", sep = "")
save.text <- "hethom.fdr.ci"
for(numc in save.text)
  write.table(loomv.uv.ci, file = paste(control$dropbox_text_numbers_dir, "/", save.text, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)



fnamc <- "hethom_scatter.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 9, 9, units = "in", res = 500)
par(mfrow = c(1, 1), oma = c(3, 2, 1, 1), mar = c(3, 3, 3, 3))
ntab <- t(tabhethom)# table(het.z.sign, hom.z.sign)
ptab <- ntab / sum(ntab)
limn <- 3
legcex <- .65
plot(het.z, hom.z, pch = 1, cex = .2, 
     xlim = c(-1, 1) * limn, ylim = c(-1, 1) * limn, ylab = "", xlab = "", xaxs = "i", yaxs = "i")
abline(0, 1)
labc1 <- "Heterozygote MV "
labc2 <- "Homozygote MV"
mtext(side = 1, line = 2.5, text = substitute(paste(labc1, " ", italic(tilde(z))), list(labc1 = labc1)))
mtext(side = 2, line = 2.5, text = substitute(paste(labc2, " ", italic(tilde(z))), list(labc2 = labc2)))
abline(h = c(-1, 1), col = 2, lty = 2)
abline(v = c(-1, 1), col = 2, lty = 2)
atv <- c((-limn - 1) / 2, 0, (limn + 1) / 2)
for(i in 1:3){
  for(j in 1:3){
    # labc <- paste(formatC(100 * ptab[i, j], format = "f", digits = 1), "% (", prettyNum(ntab[i, j], big.mark = ","), ")", sep = "")
    labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
    legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .75), xjust = .5, cex = legcex)
  }
}
# mtext(side = 3, text = paste("(", letters[match(namc, namv)], ")", sep = ""), at = -limn, line = 1, cex = 1.3)
numv <- formatC(c(fdrci[[1]], fdrci[[2]]) * 100, format = "f", digits = 1)
fdrlab <- paste("Estimated Fsr = ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
mtext(side = 3, text = fdrlab, line = .5, cex = 1)
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)


#global concordance/discordance
hethom.concordant.n <- tabhethom["-1", "-1"] + tabhethom["1", "1"]
hethom.discordant.n <- tabhethom["-1", "1"] + tabhethom["1", "-1"]
save.num <- c("hethom.concordant.n", "hethom.discordant.n")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)


#het/hom hit rate
hom.hitrate <- mean(hom.z.sign != 0, na.rm = T)
het.hitrate <- mean(het.z.sign != 0, na.rm = T)
save.num1 <- c("hom.hitrate", "het.hitrate")
for(numc in save.num1)
  write.table(formatC(eval(as.name(numc)) * 100, format = "f", digits = 1), 
              file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)




extra.plots <- F
if(extra.plots){
  # resimp.hom$poscon = resimp.het$poscon = sign(resimp.hom$sumsig) == sign(resimp.het$sumsig) & resimp.het$sumsig != 0
  # hethom.poscon = data.frame(ph = resimp.hom[which(resimp.hom$poscon), "ph"], genosh = resimp.hom[which(resimp.hom$poscon), "genosh"])
  # for(j in 1:nrow(hethom.poscon)){
  #   resimp[resimp$ph == hethom.poscon[j, "ph"] & resimp$genosh == hethom.poscon[j, "genosh"], "poscon"] = T  
  # }
  
  
  
  ###################################################
  #Positive controls analysis
  #Pull out ref lines
  refline.poscon = data.frame()
  for(genc in names(tl.post)){
    sigmat.post = sign(tl.post[[genc]]) * (abs(tl.post[[genc]]) > thl.post[[genc]])
    sigmat.pre = sign(tl.pre[[genc]]) * (abs(tl.pre[[genc]]) > thl.pre[[genc]])
    sigmat.pre[is.na(sigmat.pre)] = 0
    sigmat.comb = abs(sigmat.post) + abs(sigmat.pre)
    sigmat.comb[] = pmin(1, sigmat.comb)
    allsamesign = apply(cbind(sigmat.pre, sigmat.post), 1, function(v) all(sign(v) <= 0) | all(sign(v) >= 0))
    atleasttwosig = rowSums(sigmat.comb) > 1
    phsig = names(allsamesign)[allsamesign & atleasttwosig]
    if(length(phsig) == 0)
      next
    add = data.frame(ph = phsig, genename = genc)
    refline.poscon =  rbind(refline.poscon, add)
  }
  resimp$poscon = F
  for(j in 1:nrow(refline.poscon)){
    zygc = strsplit(refline.poscon[j, "genename"], spl = "_")[[1]][2]
    gennamec = strsplit(refline.poscon[j, "genename"], spl = "_")[[1]][1]
    gennums = ref[ref$gene_symbol == gennamec, "genotype_id"]
    phc = refline.poscon[j, "ph"]
    resimp[which(resimp$ph == phc & resimp$geno %in% paste(gennums, zygc, sep = "_")), "poscon"] = T  
  }
  sum(resimp$poscon)
  
  resimp.pos = resimp[resimp$poscon & !resimp$negcon, ]
  print(table(sign(resimp.pos$uv.t) * (abs(resimp.pos$uv.t) > resimp.pos$uv.th.final),
              sign(resimp.pos$eb.t) * (abs(resimp.pos$eb.t) > resimp.pos$eb.th.final)))
  
  
  
  #########################################
  #Examine meta-analysis output
  mean(eb.mean <= pre, na.rm = T)
  mean(post <= pre, na.rm = T)
  mean(post.comp <= pre, na.rm = T)
  plot(pre, post.comp)
  abline(0, 1)
  # plot(pre, eb.old)
  # abline(0, 1)
  ith = 0.55
  mean(pre <= ith, na.rm = T)
  mean(post <= ith, na.rm = T)
  mean(post.comp <= ith, na.rm = T)
  mean(eb.mean <= ith, na.rm = T)
  #mean(f <= ith, na.rm = T)
  
  #boxplot(resimp$mv.mn ~ resimp$cen)
  cenmnmad = sapply(unique(resimp$cen), function(cenc) mad(resimp[resimp$cen == cenc & resimp$line.type == control$nam.truemut, "mv.mn"], na.rm = T))
  censdmad = sapply(unique(resimp$cen), function(cenc) mad(resimp[resimp$cen == cenc & resimp$line.type == control$nam.truemut, "mv.sd"], na.rm = T))
  cenmnmad
  censdmad
  
  
  #hist(sort(na.omit(rowMeans(post < pre, na.rm = T))))
  # library(tools)
  # install.packages(package_dependencies("Hmisc")$Hmisc)
  library(Hmisc)
  mean(post[!is.na(pre)] <= ith, na.rm = T)
  mvnonimplog = na.omit(post[which(!is.na(pre))] <= ith)
  mvimplog = na.omit(post[which(is.na(pre))] <= ith)
  eb.mean.nonimplog = na.omit(eb.mean[which(!is.na(pre))] <= ith)
  eb.mean.mvimplog = na.omit(eb.mean[which(is.na(pre))] <= ith)
  uvlog = na.omit(pre[which(!is.na(pre))] <= ith)
  flog = na.omit(c(f <= ith))
  binconf(x = sum(uvlog), n = length(mvnonimplog))
  binconf(x = sum(mvnonimplog), n = length(mvnonimplog))
  binconf(x = sum(mvimplog), n = length(mvimplog))
  binconf(x = sum(eb.mean.nonimplog), n = length(eb.mean.nonimplog))
  binconf(x = sum(eb.mean.mvimplog), n = length(eb.mean.mvimplog))
  save.tab <- c("uv_mv_prop_table")
  for(numc in save.tab)
    write.table(formatC(100 * ptab[i, j], format = "f", digits = 2), 
                file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
                col.names = F, row.names = F, quote = F)
  
  ##  
  
  str(postpaircomp)
  str(postpairimp)
  
  
  postpair[1:10, ]
  postpaircomp[1:10, ]
  postpairimp[1:10, ]
  
  sum(!is.na(postpaircomp)) + sum(!is.na(postpairimp)) + sum(!is.na(postpairimpcomp))
  sum(!is.na(postpair))
  
  
  str(postpaircomp)
  str(postpairimp)
  
  prepair <- t(apply(prepair, 1, function(v) sort(v, decreasing = T)))
  tab.pre <- table(prepair[, 1], prepair[, 2])
  postpaircomp <- t(apply(postpaircomp, 1, function(v) sort(v, decreasing = T)))
  tab.post.comp <- table(postpaircomp[, 1], postpaircomp[, 2])
  postpairimp <- t(apply(postpairimp, 1, function(v) sort(v, decreasing = T)))
  tab.post.imp <- table(postpairimp[, 1], postpairimp[, 2])
  postpairimpcomp <- t(apply(postpairimpcomp, 1, function(v) sort(v, decreasing = T)))
  tab.post.impcomp <- table(postpairimpcomp[, 1], postpairimpcomp[, 2])
  postpair <- t(apply(postpair, 1, function(v) sort(v, decreasing = T)))
  tab.post <- table(postpair[, 1], postpair[, 2])
  print(glob.fdr.imp)
  print(glob.fdr.non)
  
  
  table(resimp$line.type)
  
  # write.csv(tab.pre, file = "X:/projects/impc_mv_analysis/tables/uv_ref_lines.csv")
  # write.csv(tab.post.comp, file = "X:/projects/impc_mv_analysis/tables/mv_comp_uv_ref_lines.csv")
  # write.csv(tab.post.imp, file = "X:/projects/impc_mv_analysis/tables/mv_imp_only_ref_lines.csv")
  # write.csv(tab.post, file = "X:/projects/impc_mv_analysis/tables/mv_ref_lines.csv")
  # n.disagree.post.for.this.geno = sum(rowSums(tabpost == 1, na.rm = T) > 0 & rowSums(tabpost == -1, na.rm = T) > 0)
  # n.disagree.post = n.disagree.post + n.disagree.post.for.this.geno
  # n.disagree.pre.for.this.geno = sum(rowSums(tabpre == 1, na.rm = T) > 0 & rowSums(tabpre == -1, na.rm = T) > 0)
  # n.disagree.pre = n.disagree.pre + n.disagree.pre.for.this.geno
  # tabpost = tabpostcomp = tabpostimp = sign(tl.post[[namc]]) * (abs(tl.post[[namc]]) > thl.post[[namc]])
  # tabpre = sign(tl.pre[[namc]]) * (abs(tl.pre[[namc]]) > thl.pre[[namc]])
  
  #propn of pairwise comparisons disagreeing
  (tab.post["-1", "1"] + tab.post["1", "-1"]) / (sum(tab.post))# - tab.post["0", "0"])
  (tab.post.comp["-1", "1"] + tab.post.comp["1", "-1"]) / (sum(tab.post.comp))# - tab.post.comp["0", "0"])
  (tab.pre["-1", "1"] + tab.pre["1", "-1"]) / (sum(tab.pre))# - tab.pre["0", "0"])
  sum(tab.post.comp)
  sum(tab.pre)
  #compare distribution of t-stats for imputed vs non-imputed
  summary(resimp$mv.t[!resimp$imputed & resimp$line.type == control$nam.negcon])
  summary(resimp$mv.t[resimp$imputed & resimp$line.type == control$nam.negcon])
  summary(resimp$mv.t[!resimp$imputed & resimp$line.type == control$nam.truemut])
  summary(resimp$mv.t[resimp$imputed & resimp$line.type == control$nam.truemut])
  
  
  table(substr(rownames(postpair[(postpair[, 1] == 1 & postpair[, 2] == -1) | (postpair[, 1] == -1 & postpair[, 2] == 1), ]), 1, 8))
  
  #############################
  #Extra Ref lines plot
  
  mnl.pre[[1]][1:10, ]
  mnl.post[[1]][1:10, ]
  
  cenun = na.omit(cenmap[match(unique(resimp$cen), cenmap$cen), "nam"])
  badphens = names(sort(rowSums(pre > .25, na.rm = T), decreasing = T))
  colnames(pre) = colnames(post) = names(mnl.pre)
  pldir = "X:/projects/impc_mv_analysis/plots/ref_lines"
  pdf(file = paste(pldir, "/ref_by_phen_ordered_by_quality.pdf", sep = ""), 12, 8)
  for(phc in badphens){#phc = rownames(mnl.pre[[1]])[1]#
    par(mfrow = c(2, 1), mar = c(4, 4, 4, 8))
    for(phase in c("pre", "post")){
      mnv = sev = linev = genev = c()
      if(phase == "pre"){
        mnlc = mnl.pre
        selc = sel.pre
        imat = pre
      }
      if(phase == "post"){
        mnlc = mnl.post
        selc = sel.post
        imat = post
      }
      for(genec in names(mnlc)){
        mnv = c(mnv, mnlc[[genec]][phc, ])
        sev = c(sev, selc[[genec]][phc, ])
        linev = c(linev, colnames(mnlc[[genec]]))
        genev = c(genev, rep(genec, ncol(mnlc[[genec]])))
      }
      refdat = data.frame(gene = genev, mn = mnv, se = sev, geno = linev)
      refdat$cen = resimp[match(refdat$geno, resimp$geno), "cen"]
      refdat$cenlong = cenmap[match(refdat$cen, cenmap$cen), "nam"]
      refdat = refdat[rowSums(is.na(refdat)) == 0, ]
      if(nrow(refdat) > 0){
        genun = unique(refdat$gene)
        refdat = refdat[order(refdat$gene, refdat$cenlong), ]
        colun = control$heat_col_palette[floor((1:length(cenun)) / length(cenun) * 1000)]
        names(colun) = cenun
        colv = colun[refdat$cenlong]
        ylimc = range(c(refdat$mn + 2 * refdat$se, refdat$mn - 2 * refdat$se))
        plot(1:nrow(refdat), refdat$mn, xaxt = "n", xlab = "", col = colv, pch = 19, cex = 1.3, 
             ylim = ylimc, xlim = c(0.5, nrow(refdat) + 0.5), ylab = "Effect estimate and 95% CI", xaxs = "i")
        abline(h = 0, lty = 3)
        abline(v = which(diff(match(refdat$gene, genun)) != 0) + .5)
        for(i in 1:nrow(refdat))
          lines(x = rep(i, 2), y = refdat$mn[i] + c(-2, 2) * refdat$se[i], col = colv[i])
        mtext(round(imat[phc, genun], 2), at = sapply(genun, function(genc) mean(which(refdat$gene == genc))))
        mtext(side = 1, genun, at = sapply(genun, function(genc) mean(which(refdat$gene == genc))))
        par(xpd = NA)
        legend(x = nrow(refdat) + .5, y = mean(ylimc), legend = cenun, col = colun, pch = 19, yjust = .5, xjust = 0)
        par(xpd = F)
        mtext(outer = T, paste(phc, Data_all$impc$phmap[match(phc, Data_all$impc$phmap$ph), "nam"]), side = 3, line = -1.5)
      }
    }
  }
  dev.off()
  
}  

# 
# par(mfrow = c(3, 3))
# for(cenc in cenmap$nam){
#   try({
#     tlook <- resimp[resimp$cenlong == cenc & resimp$line.type == c("negCon", control$nam.truemut)[1], "uv.t"]
#     hist(tlook, main = cenc)
#     print(cenc)
#     print(mean(abs(tlook) > 5, na.rm = T))
#     
#   })
# }
# 
# table(paste(resimp$cen, resimp$line.type))
# sum(abs(tlook) > 6, na.rm = T)
# # source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")
# # nfac <- 24
# # vc.type <- c("vari", "pro")[1]
# # load(file = paste0(meth.comp.output.dir, "/ebmix_results_nf_", nfac, "_vt_", vc.type, ".RData"))
# # load(file = uv.results.Y.S)
# # str(resimp)
# # load(file = resimp.with.sig.thresholds.file)
# # resimp$geno.sh = sapply(strsplit(resimp$geno, spl = "_"), function(x) x[1])
# 
# 
# # test.stat <- "z"
# # sep.imp.thresh <- F
# # calibrate <- F
# # err.rate.meth <- "perm"
# # Sig <- emout.mix$Sig
# # R <- emout.mix$R
# # piv <- emout.mix$pi
# 
# # source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
# # err.rate.out <- err.rate.control(resl = resl, err.rate.meth = err.rate.meth, calibrate = calibrate, sep.imp.thresh = sep.imp.thresh, 
# #                                  test.stat = test.stat, linemap = linemap, phmap = phmap, Yhat = Yhat)
# # err.rate.out$restab
# # resl <- err.rate.out$resl
# # resimp <- err.rate.out$resimp
# # resimp0 <- resimp[resimp$line.type == "negConCheck", ]
# # resimp1 <- resimp[resimp$line.type == control$nam.truemut, ]
# # sd.null.non <- sd.null.non.uv <- sd.null.imp <- sig.null.imp <- sig.null.non <- sig.null.non.uv <- 
# #   sig.true.imp <- sig.true.non <- sig.true.non.uv <- c()
# # for(cenc in cennamun){
# #   sd.null.imp[cenc] <- mad(resimp0[resimp0$cenlong == cenc & resimp0$imputed, "eb.mix.t"], na.rm = T)
# #   sd.null.non[cenc] <- mad(resimp0[resimp0$cenlong == cenc & !resimp0$imputed, "eb.mix.t"], na.rm = T)
# #   sd.null.non.uv[cenc] <- mad(resimp0[resimp0$cenlong == cenc & !resimp0$imputed, "uv.t"], na.rm = T)
# #   sig.null.imp[cenc] <- mean(abs(resimp0[resimp0$cenlong == cenc & resimp0$imputed, "eb.mix.perm.signsig"]), na.rm = T)
# #   sig.null.non[cenc] <- mean(abs(resimp0[resimp0$cenlong == cenc & !resimp0$imputed, "eb.mix.perm.signsig"]), na.rm = T)
# #   sig.null.non.uv[cenc] <- mean(abs(resimp1[resimp1$cenlong == cenc & !resimp1$imputed, "uv.perm.signsig"]), na.rm = T)
# #   sig.true.imp[cenc] <- mean(abs(resimp1[resimp1$cenlong == cenc & resimp1$imputed, "eb.mix.perm.signsig"]), na.rm = T)
# #   sig.true.non[cenc] <- mean(abs(resimp1[resimp1$cenlong == cenc & !resimp1$imputed, "eb.mix.perm.signsig"]), na.rm = T)
# #   sig.true.non.uv[cenc] <- mean(abs(resimp1[resimp1$cenlong == cenc & !resimp1$imputed, "uv.perm.signsig"]), na.rm = T)
# #   # mean(is.na(resimp0[resimp$cenlong == cenc, "eb.mix.t"]))
# # }
# # lookmat <- cbind(sd.null.non.uv, sd.null.non, sd.null.imp, sig.null.non.uv, sig.null.non, sig.null.imp, sig.true.non.uv, sig.true.non, sig.true.imp)
# # round(lookmat, 2)
# 
# postmat <- resl$eb$ref.lines$matl$post
# which.discordant <- which((postmat[, 1] > 1 & postmat[, 2] < -1) | (postmat[, 1] < -1 & postmat[, 2] > 1))
# str(postmat)
# tabc <- table(rownames(postmat)[which.discordant])
# sum(tabc)
# sort(tabc)
# 
# phen.un <- unique(resimp$ph)
# # phen.un <- phen.un[!phen.un %in% paste0("f.", 1:nfac)]
# resimp <- resimp[resimp$ph %in% phen.un, ]
# resimp <- resimp[resimp$line.type == control$nam.truemut, ]
# resimp$uv.sig <- abs(resimp$uv.perm.signsig)
# resimp$eb.sig <- abs(resimp$eb.perm.signsig)
# resimp$eb.t <- resimp$eb.t
# resimp$eb.th.final <- resimp$eb.th.final
# resimp$eb.signsig <- resimp$eb.perm.signsig
# resimp$uv.signsig <- resimp$uv.perm.signsig
# 
# 
# 
# 
# ph.check <- c("IMPC_ACS_035_001", "JAX_SLW_014_001", "IMPC_GRS_011_001")[3]
# pout[pout$ph == ph.check, ]
# rescheck <- resimp[resimp$ph == ph.check, ]
# tabl <- list()
# for(cenc in cenun){
#   namc <- cenmap[match(cenc, cenmap$cen), "nam"]
#   tabc <- table(rescheck[rescheck$cen == cenc, "eb.perm.signsig"])
#   # tabl[[namc]] <- round(100 * tabc / sum(tabc))
#   tabl[[namc]] <- tabc
# }
# tabl
# str
# table(resche)
# str(resimp0)
# str(resimp)

# 
# 
# ########################################################
# #Extract reference lines results
# proc.omit = ""
# ph.imputed = pout[pout$procnam == proc.omit, "ph"]
# ref = read.csv("W:/projects/major/impc_stats_analysis/data_in/impc_data/reference_line_genotypeIds.csv")
# linun = unique(ref$gene_symbol)
# phun = unique(resimp$ph)
# resimp$zyg = sapply(strsplit(resimp$geno, spl = "_"), function(x) x[2])
# tl.post = thl.post = tl.pre = thl.pre = tl.f = mnl.pre = sel.pre = mnl.post = sel.post = list()
# pre = post = post.comp = f = eb.old = eb.initial = eb.mean = NULL
# for(linc in linun){#linc = "Prkab1"#
#   for(zyg in 0:1){#zyg = 0#
#     namc = paste(linc, zyg, sep = "_")
#     gen.shv = ref[ref$gene_symbol == linc, "genotype_id"]
#     # cen.gen <- resimp[match(paste(gen.shv, zyg, sep = "_"), resimp$geno), "cenlong"]
#     resc = resimp[resimp$zyg == zyg & resimp$line.type == control$nam.truemut, ]
#     genv = na.omit(resc[match(gen.shv, resc$geno.sh), "geno"])
#     mn.post = se.post = t.post = th.post = mn.f = se.f = mn.pre = se.pre = t.pre = th.pre = NULL
#     if(length(genv) <= 1)
#       next
#     for(genc in genv){
#       resimpc = resimp[resimp$geno == genc, ]
#       mn.post = cbind(mn.post, resimpc[match(phun, resimpc$ph), "eb.mn"])
#       th.post = cbind(th.post, resimpc[match(phun, resimpc$ph), "eb.th.final"])
#       se.post = cbind(se.post, resimpc[match(phun, resimpc$ph), "eb.sd"])
#       t.post = cbind(t.post, resimpc[match(phun, resimpc$ph), "eb.t"])
#       mn.pre = cbind(mn.pre, resimpc[match(phun, resimpc$ph), "uv.mn"])
#       th.pre = cbind(th.pre, resimpc[match(phun, resimpc$ph), "uv.th.final"])
#       se.pre = cbind(se.pre, resimpc[match(phun, resimpc$ph), "uv.sd"])
#       t.pre = cbind(t.pre, resimpc[match(phun, resimpc$ph), "uv.t"])
#     }
#     rownames(mn.post) <- rownames(se.post) <- rownames(t.post) <- rownames(th.post) <- 
#       rownames(mn.pre) <- rownames(se.pre) <- rownames(t.pre) <- rownames(th.pre) <- phun
#     colnames(mn.post) <- colnames(se.post) <- colnames(t.post) <- colnames(th.post) <- 
#       colnames(mn.pre) <- colnames(se.pre) <- colnames(t.pre) <- colnames(th.pre) <- genv
#     mn.post.comp = mn.post.imp = mn.post
#     se.post.comp = se.post.imp = se.post
#     mn.post.comp[is.na(mn.pre)] <- se.post.comp[is.na(se.pre)] <- 
#       mn.post.imp[!is.na(mn.pre)] <- se.post.imp[!is.na(se.pre)] <- NA
#     ipost = calc.i2(mn.post, se.post)
#     post = cbind(post, ipost$I)
#     ipost.comp = calc.i2(mn.post.comp, se.post.comp)
#     post.comp = cbind(post.comp, ipost.comp$I)
#     ipre = calc.i2(mn.pre, se.pre)
#     pre = cbind(pre, ipre$I)
#     ieb.mean = calc.i2(mn.post, se.post)
#     eb.mean = cbind(eb.mean, ieb.mean$I)
#     tl.post[[namc]] = t.post
#     thl.post[[namc]] = th.post
#     tl.pre[[namc]] = t.pre
#     thl.pre[[namc]] = th.pre
#     mnl.pre[[namc]] = mn.pre
#     sel.pre[[namc]] = se.pre
#     mnl.post[[namc]] = mn.post
#     sel.post[[namc]] = se.post
#   }
# }
# 
# ##################################################################
# #Analyse reference lines results for hit concordance
# n.disagree.post = n.disagree.pre = 0 
# matl <- list()
# postpair = postpaircomp = postpairimp = postpairimpcomp = prepair = NULL
# for(namc in names(tl.post)){
#   tabpost = tabpostcomp = tabpostimp = tl.post[[namc]] / thl.post[[namc]]
#   tabpre = tl.pre[[namc]] / thl.pre[[namc]]
#   tabpostcomp[is.na(tabpre)] = NA
#   tabpostimp[!is.na(tabpre)] = NA
#   for(i in 1:(ncol(tabpost) - 1)){
#     for(j in (i + 1):ncol(tabpost)){
#       matl$post <- rbind(matl$post, tabpost[, c(i, j)])
#       matl$postcomp = rbind(matl$postcomp, tabpostcomp[, c(i, j)])
#       matl$postimp = rbind(matl$postimp, tabpostimp[, c(i, j)])
#       matl$postimpcomp <- rbind(matl$postimpcomp, cbind(tabpostimp[, i], tabpostcomp[, j]), 
#                                cbind(tabpostimp[, j], tabpostcomp[, i]))
#       matl$pre = rbind(matl$pre, tabpre[, c(i, j)])
#     }
#   }
# }
# matl <- lapply(matl, function(M) M[rowSums(is.na(M)) == 0, ])
# matl.sig <- lapply(matl, function(M){M[abs(M) < 1] <- 0; sign(M)})
# matl.sig.ord <- lapply(matl.sig, function(M)t(apply(M, 1, function(v) sort(v, decreasing = T))))
# tabl <- lapply(matl.sig.ord, function(M) table(M[, 1], M[, 2]))
# nfac <- 24
# vc.type <- c("vari", "pro")[1]
# load(file = paste0(meth.comp.output.dir, "/ebmix_results_nf_", nfac, "_vt_", vc.type, ".RData"))
# load(file = file.resl.comp)
