
linemap.use <- Data_all$impc$linemap
linemap.use$line.type <- ifelse(linemap$line.type == "trueMut", "trueMutTes", "negConTes")
linemap.use <- linemap.use[linemap.use$geno %in% rownames(resl.comp$eb$mn), ]

# objl[[control$mv_meth_nam_use]][[1]]$Sigl[[1]]

# Sig.comb <- resl.comp$impc_eb_1$Sig.comb
# str(Sigll, m = 1)

#####################################################################
# Calculate KL divergence between split models and combined model for Factor Sensitivity Analysis
kl1 <- kl2 <- c()
Sig.comb <- resl.comp[[control$mv_meth_nam_use]]$Sig.comb
P <- control$default_parameters$impc$P
for(seed in 1:control$n_subsamples_main){
  Sig_for_curr_data_split <- objl[[control$mv_meth_nam_use]][[seed]]$Sigl[[1]]
  if(!is.null(Sig_for_curr_data_split)){
    kl1[seed] <- .5 * (sum(diag(solve(Sig.comb) %*% Sig_for_curr_data_split)) - P + 
                         determinant(Sig.comb, logarithm = T)$modulus - determinant(Sig_for_curr_data_split, logarithm = T)$modulus)
    kl2[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split) %*% Sig.comb)) - P + 
                         determinant(Sig_for_curr_data_split, logarithm = T)$modulus - determinant(Sig.comb, logarithm = T)$modulus)
  }
}



names(resl.comp)
plot(kl1, kl2)



seed.largest.kl <- which.min(kl1 + kl2)
Sig_main <- Sig.comb
Sig_compare <- objl[[control$mv_meth_nam_use]][[seed.largest.kl]]$Sigl[[1]]
Sig_corr_compare <- t(Sig_compare / sqrt(diag(Sig_compare))) / sqrt(diag(Sig_compare))
eigc_compare <- eigen(Sig_corr_compare)
Sig_corr_main <- t(Sig_main / sqrt(diag(Sig_main))) / sqrt(diag(Sig_main))
eigc_main <- eigen(Sig_corr_main)
loadings_main <- varimax(eigc_main$vectors[, 1:control$nfac])$loadings
loadings_compare <- varimax(eigc_compare$vectors[, 1:control$nfac])$loadings
loadings_main <- sweep(loadings_main, 2, apply(loadings_main, 2, function(v) v[which.max(abs(v))]), "/")
loadings_compare <- sweep(loadings_compare, 2, apply(loadings_compare, 2, function(v) v[which.max(abs(v))]), "/")
# apply(loadings_main, 2, max)

t(loadings_main) %*% loadings_main
rownames(loadings_compare) <- rownames(Sig_corr_compare)
rownames(loadings_main) <- rownames(Sig_corr_main)
loadings_compare <- loadings_compare[rownames(loadings_main), ]
# reorder compare loadings to match main
prox_mat <- abs(t(loadings_main) %*% loadings_compare)
map_to <- rep(0, control$nfac)
# for (j in rank(-apply(prox_mat, 2, max), )) {
for (j in 1:control$nfac) {
    map_to[j] <- which.max(ifelse(1:control$nfac %in% map_to, -1, abs(t(loadings_main) %*% loadings_compare[, j])))
    # map_to[j] <- which.max(abs(t(loadings_main) %*% loadings_compare[, j]))
    # map_to[j] <- which.max(abs(t(loadings_main) %*% loadings_compare[, j]))
}

j=1
abs(t(loadings_main) %*% loadings_compare[, j])
loadings_compare[, map_to] <- loadings_compare
loadings_compare <- sweep(loadings_compare, 2, sign(colMeans(loadings_main * loadings_compare)), "*")
# loadings_compare <- loadings_compare[, apply(abs(t(loadings_main) %*% loadings_compare), 1, which.max)]
par(mfrow = c(1, 2))
image(t(loadings_main), zlim = c(-1, 1), col = rain)
image(t(loadings_compare), zlim = c(-1, 1))
# image(t(resl.comp[[control$mv_meth_nam_use]]$facs.varimax[rownames(loadings_main), ]))

image(abs(t(loadings_main) %*% loadings_compare))
plot(resl.comp[[control$mv_meth_nam_use]]$Sigcor.comb, Sig_corr_main)

plot(resl.comp[[control$mv_meth_nam_use]]$facs.varimax[rownames(loadings_main), 4],
      loadings_main[, 4])

dimnames(loadings_main)

par(mfrow = c(1, 2))
image(Sig_main)
image(Sig_compare)


str(fac.res.store)
head(fac.res.store[[seed.largest.kl]]$loadings)
head(resl.comp[[control$mv_meth_nam_use]]$facs.varimax)

head(fac.res.store[[1]]$loadings)
head(fac.res.store[[2]]$loadings)


#####################################################
# Checking concordance across splits when using LFSR as error rate control
compl <- lapply(compl, function(x){ x$signsigarr.lfsr <- as.numeric(x$lfsrarr < .025) * sign(x$mnarr); return(x) })
str(compl, m = 2)
table(c(compl$mash$signsigarr.lfsr[, , 1]), c(compl$mash$signsigarr.lfsr[, , 2]))
table(c(compl$impc_eb_1$signsigarr.lfsr[, , 1]), c(compl$impc_eb_1$signsigarr.lfsr[, , 2]))



resl.comp.largest.kl <- list(eb.glob = resl.out$eb, 
                             eb.furthest.kl = list(mn = compl$eb$mnarr[, , seed.largest.kl], sd = compl$eb$sdarr[, , seed.largest.kl]))
out.perm <- err.rate.control(resl = resl.comp.largest.kl, 
                             err.rate.meth = "perm", sep.imp.thresh = F,
                             test.stat = "z", linemap = linemap.use, phmap = phmap, Yhat = Yhat,
                             use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                             err.thresh = .05, centre.specific.thresh = F)

out.perm$restab
table(out.perm$resimp$eb.glob.perm.signsig, out.perm$resimp$eb.furthest.kl.perm.signsig)





# 
# #################################################
# # looking at concordance of methods across splits, based on permutation based FDR control
# resl.comp.largest.kl <- resl.comp
# file.base.largest.kl <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name, tab.all.seed[seed.largest.kl, var.in.name], sep = "_"), collapse = "_"))
# resl.file.largest.kl <- paste0(file.base.largest.kl, "_resl.RData")
# load(file = resl.file.largest.kl)
# resl.comp.largest.kl$eb.furthest <- resl.store[[1]][c("mn", "sd", "lfsr")]
# mash.file.base.largest.kl <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, tab.all.seed[seed.largest.kl, var.in.name.mash], sep = "_"), collapse = "_"))
# mash.resl.file.largest.kl <- paste0(mash.file.base.largest.kl, "_mash_resl.RData")
# load(file = mash.resl.file.largest.kl)
# resl.comp.largest.kl$mash.furthest <- resl.store[[1]][c("mn", "sd", "lfsr")]
# lines.largest.kl <- rownames(resl.comp.largest.kl$eb.furthest$mn)
# resl.comp.largest.kl <- rapply(resl.comp.largest.kl, function(x) return(x[lines.largest.kl, ph.use]), how = "replace")
# str(resl.comp.largest.kl,m=2)
# source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
# out <- err.rate.control(resl = resl.comp.largest.kl, err.rate.meth = "perm",
#                         calibrate = calibrate, sep.imp.thresh = sep.imp.thresh, test.stat = "z",
#                         linemap = linemap[match(lines.largest.kl, linemap$geno), ], phmap = phmap, Yhat = Yhat,
#                         use.upper.fp.est = use.upper.fp.est, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                         err.thresh = .05)
# out$restab
# table(out$resimp$eb.perm.signsig, out$resimp$eb.furthest.perm.signsig)
# table(out$resimp$mash.perm.signsig, out$resimp$mash.furthest.perm.signsig)





#