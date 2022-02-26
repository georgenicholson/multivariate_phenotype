xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))



load(file = uv.results.Y.S)
load(file = file.ebi.impc.results)
load(file = file.glob.res)







linemap.use <- linemap
linemap.use$line.type <- ifelse(linemap$line.type == "trueMut", "trueMutTes", "negConTes")
linemap.use <- linemap.use[linemap.use$geno %in% rownames(resl.comp$eb$mn), ]


Sig.comb <- resl.comp$impc_eb_1$Sig.comb
str(Sigll, m = 1)
#####################################################################
# Calculate KL divergence between split models and combined model
kl1 <- kl2 <- c()
for(seed in 1:n.seed.run.subsams){
  if(!is.null(Sigll$impc_eb_1[[seed]][[1]])){
    kl1[seed] <- .5 * (sum(diag(solve(Sig.comb) %*% Sigll$impc_eb_1[[seed]][[1]])) - P + 
                         determinant(Sig.comb, logarithm = T)$modulus - determinant(Sigll$impc_eb_1[[seed]][[1]], logarithm = T)$modulus)
    kl2[seed] <- .5 * (sum(diag(solve(Sigll$impc_eb_1[[seed]][[1]]) %*% Sig.comb)) - P + 
                         determinant(Sigll$impc_eb_1[[seed]][[1]], logarithm = T)$modulus - determinant(Sig.comb, logarithm = T)$modulus)
  }
}


plot(kl1, kl2)
seed.largest.kl <- which.max(kl1 + kl2)

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