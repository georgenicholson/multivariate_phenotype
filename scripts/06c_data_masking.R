resimp <- readRDS(file = control$file.resimp)
resimp <- resimp[!grepl("fac", resimp$ph), ]

resimp$eb.mn.loo.proc <- compl[[control$mv_meth_nam_use]]$loocv.mn.comb[cbind(resimp$geno, resimp$ph)]
resimp$eb.sd.loo.proc <- compl[[control$mv_meth_nam_use]]$loocv.sd.comb[cbind(resimp$geno, resimp$ph)]
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
  write.table(eval(as.name(numc)), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)


###############################################
#Leave-one-out imputed vs UV scatterplot
fnamc <- "looEB_uv_scatter.jpg"
jpeg(paste(control$figure_dir, "/", fnamc, sep = ""), 7, 7, units = "in", res = 500)
par(mfrow = c(1, 1), oma = c(3, 3, 3, 3))
for(compc in c("uv", "eb")[1]){
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
      labc <- paste(prettyNum(ntab[i, j], big.mark = ","), " (", formatC(100 * ptab[i, j], format = "f", digits = 1), "%)", sep = "")
      legend(x = atv[i], y = atv[j], legend = labc, text.col = 1, bg = adjustcolor("white", .85), xjust = .5, cex = legcex)
    }
  }
  # mtext(side = 3, line = 2, at = -limn, text = ifelse(compc == "uv", "(a)", "(b)"), cex = 1.3)
  if(compc == "uv")
    fdr.res <- fdr.est.tab(tab.loo.uv)
  if(compc == "eb")
    fdr.res <- fdr.est.tab(tab.loo.mv)
  numv <- formatC(c(fdr.res[[1]], fdr.res[[2]]) * 100, format = "f", digits = 1)
  fdrlab <- paste("Estimated Fsr = ", numv[1], "% (", numv[2], "% - ", numv[3], "%)", sep = "")
  mtext(side = 3, text = fdrlab, line = .5, cex = 1)
}
dev.off()
file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
          to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)

