


#####################################################################
# Calculate KL divergence between split models and combined model
kl1 <- kl2 <- c()
for(seed in 1:n.seed.run.subsams){
  if(!is.null(Sigll[[seed]][[1]])){
    kl1[seed] <- .5 * (sum(diag(solve(Sig.comb) %*% Sigll[[seed]][[1]])) - P + 
                         determinant(Sig.comb, logarithm = T)$modulus - determinant(Sigll[[seed]][[1]], logarithm = T)$modulus)
    kl2[seed] <- .5 * (sum(diag(solve(Sigll[[seed]][[1]]) %*% Sig.comb)) - P + 
                         determinant(Sigll[[seed]][[1]], logarithm = T)$modulus - determinant(Sig.comb, logarithm = T)$modulus)
  }
}
plot(kl1, kl2)
seed.largest.kl <- which.max(kl1 + kl2)


#####################################################
# Checking concordance across splits when using LFSR as error rate control
compl <- lapply(compl, function(x){ x$signsigarr.lfsr <- as.numeric(x$lfsrarr < .025) * sign(x$mnarr); return(x) })
str(compl, m = 2)
table(c(compl$mash$signsigarr.lfsr[, , 1]), c(compl$mash$signsigarr.lfsr[, , 2]))
table(c(compl$eb$signsigarr.lfsr[, , 1]), c(compl$eb$signsigarr.lfsr[, , 2]))


######################################################
# Beginning to plot factor comparison for reviewer sensitivity analysis  
par(mfrow = c(1, 2))#layout(matrix(c(1, 1, 1, 2, 3, 4), 3, 2))
par(oma = c(4, 20, 4, 4))
for(j in 1:2){
  if(j == 1)
    Sig.mn <- Sig.comb
  if(j == 2)
    Sig.mn <- Sigll[[seed.largest.kl]][[1]]
  corout <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
  nfac.pl <- min(floor(P / 4), nfac)
  loadmat <- varimax(svd(corout)$v[, 1:nfac.pl])$loadings
  loadmat <- sweep(loadmat, 2, apply(loadmat, 2, function(v) v[which.max(abs(v))]), '/')
  loadmat <- loadmat[, order(colMeans(abs(loadmat) * 1:nrow(loadmat)))]
  image(x = 1:nfac.pl, y = 1:P, t(loadmat), yaxt = "n", ylab = "", zlim = c(-1, 1), col = rain, xaxt = "n", xlab = "")
  if(j == 1)
    axis(side = 2, labels = pout[match(ph.use, pout$ph), "nam"], at = 1:P, las = 2, cex.axis = .81)
}

# ll.comb <- sapply(compl, function(x) rowMeans(x$llmat, na.rm = T))
# all.true.kos <- linemap$geno[linemap$old.line.type == "trueMut"]
# plot(ll.comb[all.true.kos, 1], ll.comb[all.true.kos, 2])
# abline(0, 1)
# geno.look <- all.true.kos[which.max(abs(ll.comb[all.true.kos, 1] - ll.comb[all.true.kos, 2]))]
restab <- data.frame()
resl.store <- list()
use.upper.fp.est <- T
calibrate <- F
sep.imp.thresh <- F
err.rate.methv <- c("perm", "lfsr")
control.level <- c("line.fdr", "line.fwer", "phcen.fdr")[1]
resl.out <- list()
for(err.rate.meth in err.rate.methv){
  if(err.rate.meth == "perm"){
    test.statv <- c("z", "lfsr")
    if(control.level == "line.fdr"){
      centre.specific.threshv <- c(T, F)
    } else {
      centre.specific.threshv <- F
    }
  }
  if(err.rate.meth == "lfsr")
    test.statv <- "lfsr"
  for(test.stat in test.statv){
    for(centre.specific.thresh in centre.specific.threshv){
      linemap$line.type <- ifelse(linemap$old.line.type == "trueMut", "trueMutTes", "negConTes")
      source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
      out <- err.rate.control(resl = resl.comp, err.rate.meth = err.rate.meth,
                              calibrate = calibrate, sep.imp.thresh = sep.imp.thresh, test.stat = test.stat,
                              linemap = linemap, phmap = phmap, Yhat = Yhat,
                              use.upper.fp.est = use.upper.fp.est, control.level = control.level,
                              err.thresh = .05, centre.specific.thresh = centre.specific.thresh)
      restab <- rbind(restab, out$restab)
    }
  }
}
restab

#################################################
# looking at concordance of methods across splits, based on permutation based FDR control
resl.comp.largest.kl <- resl.comp
file.base.largest.kl <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name, tab.all.seed[seed.largest.kl, var.in.name], sep = "_"), collapse = "_"))
resl.file.largest.kl <- paste0(file.base.largest.kl, "_resl.RData")
load(file = resl.file.largest.kl)
resl.comp.largest.kl$eb.furthest <- resl.store[[1]][c("mn", "sd", "lfsr")]
mash.file.base.largest.kl <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, tab.all.seed[seed.largest.kl, var.in.name.mash], sep = "_"), collapse = "_"))
mash.resl.file.largest.kl <- paste0(mash.file.base.largest.kl, "_mash_resl.RData")
load(file = mash.resl.file.largest.kl)
resl.comp.largest.kl$mash.furthest <- resl.store[[1]][c("mn", "sd", "lfsr")]
lines.largest.kl <- rownames(resl.comp.largest.kl$eb.furthest$mn)
resl.comp.largest.kl <- rapply(resl.comp.largest.kl, function(x) return(x[lines.largest.kl, ph.use]), how = "replace")
str(resl.comp.largest.kl,m=2)
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
out <- err.rate.control(resl = resl.comp.largest.kl, err.rate.meth = "perm",
                        calibrate = calibrate, sep.imp.thresh = sep.imp.thresh, test.stat = "z",
                        linemap = linemap[match(lines.largest.kl, linemap$geno), ], phmap = phmap, Yhat = Yhat,
                        use.upper.fp.est = use.upper.fp.est, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                        err.thresh = .05)
out$restab
table(out$resimp$eb.perm.signsig, out$resimp$eb.furthest.perm.signsig)
table(out$resimp$mash.perm.signsig, out$resimp$mash.furthest.perm.signsig)



uv.ash.out <- ash(betahat = na.omit(c(resl$uv$mn)), sebetahat = na.omit(c(resl$uv$sd)))
resl <- lapply(resl, function(x) { x$z.pval <- 2 * pnorm(-abs(x$mn / x$sd)); return(x) })

str(resl,m=2)


z.comb <- (mn.comb / sd.comb)[lines.truemut, ]
p.comb <- pnorm(-abs(z.comb)) * 2
q.comb <- p.adjust(p.comb, meth = "BY")
resl$eb$fdr.bh <- mn.comb


str(signsigarr)
sum.neg.eff <- apply(signsigarr == -1, 1:2, function(v) sum(v, na.rm = T))
sum.pos.eff <- apply(signsigarr == 1, 1:2, function(v) sum(v, na.rm = T))
table(c(sum.neg.eff > 0 & sum.pos.eff > 0))
sum(c(sum.neg.eff > 0 | sum.pos.eff > 0))

