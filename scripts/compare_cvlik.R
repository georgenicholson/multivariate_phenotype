rm(list = ls())
Data <- "eqtl"
xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
source(paste0(R.file.dir, "/impc_mv_paper_code/EM_fns_lean.R"))

load(file = file.resl.comp)
load(file = file.compl)
load(file = file.objl)
load(file = uv.results.Y.S)
truemuts <- linemap$geno[linemap$line.type == "trueMut"]
table(linemap$line.type)
negcons <- linemap$geno[linemap$line.type == "negCon"]

names(compl)

ebnam <- "impc_em.fit_nSig_1_fm_pca_bic_0_mash_FALSE"#"impc_em.fit_nSig_1_fm_pca_bic_0_mash_TRUE"#impc_em.fit_nSig_1_fm_fa"
mashnam <- "impc_mash_nSig_1"
ednam <- "impc_ed_nSig_1"
ebnam.eqtl <- "eqtl_em.fit_nSig_1_fm_fa"
mashnam.eqtl <- "eqtl_mash_nSig_1_fm_NA"
ednam.eqtl <- "eqtl_ed_nSig_1_fm_NA"

impc.llmean.splitmat.raw <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[truemuts, ], na.rm = T))
eqtl.llmean.splitmat.raw <- sapply(compl[grepl("eqtl", names(compl))], function(x) colMeans(x$llmat.raw[truemuts, ], na.rm = T))
# llmean.splitmat.zero <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.zero[truemuts, ], na.rm = T))
colMeans(impc.llmean.splitmat.raw, na.rm = T)
colMeans(eqtl.llmean.splitmat.raw, na.rm = T)



# colMeans(llmean.splitmat.zero, na.rm = T)
llmean.splitmat.raw <- sapply(compl[grepl("impc", names(compl))], function(x) apply(x$llmat.raw[truemuts, ], 2, function(v) median(v, na.rm = T)))
# llmean.splitmat.zero <- sapply(compl[grepl("impc", names(compl))], function(x) apply(x$llmat.zero[truemuts, ], na.rm = T))
colMeans(llmean.splitmat.raw, na.rm = T)
# colMeans(llmean.splitmat.zero, na.rm = T)

sapply(objl, function(lc) lc[[1]]$Ksig)

library(MCMCpack)
seed <- 10
ebnam <- "impc_em.fit_nSig_1_fm_fa_bic_1_mash_FALSE"#impc_em.fit_nSig_1_fm_pca_bic_0_mash_FALSE"#"impc_em.fit_nSig_1_fm_pca_bic_0_mash_TRUE"#impc_em.fit_nSig_1_fm_fa"
mashnam <- "impc_mash_nSig_1"
ednam <- "impc_ed_nSig_2"

# ldiwish(W = solve(objl[[ebnam]][[seed]]$Sighat), v = P + 1, S = diag(rep(1, P))) / N
# ldiwish(W = solve(objl[[ednam]][[seed]]$Sighat), v = P + 1, S = diag(rep(1, P))) / N
# ldiwish(W = solve(objl[[mashnam]][[seed]]$Sighat), v = P + 1, S = diag(rep(1, P))) / N
# ldiwish(W = solve(objl[[ebnam]][[seed]]$Sigl[[1]]), v = P + 1, S = diag(rep(1, P))) / N
# ldiwish(W = solve(objl[[ednam]][[seed]]$Sigl[[1]]), v = P + 1, S = diag(rep(1, P))) / N
# ldiwish(W = solve(objl[[mashnam]][[seed]]$Sigl$U.1), v = P + 1, S = diag(rep(1, P))) / N


plot(compl$impc_ed_nSig_2_fm_NA$llmat.raw, compl$impc_mash_nSig_1_fm_NA$llmat.raw)
abline(0, 1)
plot(compl$impc_em.fit_nSig_1_fm_fa$llmat.raw, compl$impc_mash_nSig_1_fm_NA$llmat.raw)
abline(0, 1)
plot(llmean.splitmat.raw[, mashnam], llmean.splitmat.raw[, ebnam])
abline(0, 1)
str(compl,m=2)
str(compl[[1]]$llmat.raw, m = 2)


plot(compl[[ebnam]]$llmat.raw[truemuts, seed], compl[[mashnam]]$llmat.raw[truemuts, seed]); abline(0, 1)
mean(compl[[ebnam]]$llmat.raw[truemuts, seed] > compl[[mashnam]]$llmat.raw[truemuts, seed], na.rm = T)
si <- 1
bo <- 1
N <- 2000
P <- 148
EDtol <- 1e-5
EDmeth <- c("justED", "mash")[1]
splitc <- seed <- 1
nSig <- 1
Data <- "impc"
file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
file.base.ed <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.ed, sapply(var.in.name.ed, get), sep = "_"), collapse = "_"))
mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
mash.raw.results.file.namc <- paste0(file.base.mash, "_mash_raw_results.RData")
load(mash.raw.results.file.namc)

sort(colMeans(res.mash.all.testing$posterior_weights))


str(objl[[ebnam]][[seed]], m =1)
str(objl[[ednam]][[seed]], m =1)
str(objl[[ednam]][[seed]], m = 2)

# load(file = paste0(meth.comp.output.dir, "/N_2000_P_148_nSig_1_seed_", seed, "_data_impc_sexspecific_FALSE_emout.RData"))
# load(file = paste0(meth.comp.output.dir, "/N_2000_P_148_nSig_1_seed_", seed, "_data_impc_sexspecific_FALSE_res.RData"))

plot(sqrt(diag(objl[[ebnam]][[seed]]$Sighat)), sqrt(diag(objl[[mashnam]][[seed]]$Sighat))); abline(0, 1)
plot(c(objl[[ebnam]][[seed]]$Sighat), c(objl[[mashnam]][[seed]]$Sighat)); abline(0, 1)
mean(abs(c(objl[[mashnam]][[seed]]$Sighat)) > abs(c(objl[[ebnam]][[seed]]$Sighat)))
mean(abs(diag(objl[[mashnam]][[seed]]$Sighat)) > abs(diag(objl[[ebnam]][[seed]]$Sighat)))
plot(diag(objl[[ebnam]][[seed]]$Sighat), diag(objl[[ednam]][[seed]]$Sighat)); abline(0, 1)
plot(diag(objl[[ebnam]][[seed]]$Sigl[[1]]), diag(objl[[ednam]][[seed]]$Sigl[[1]])); abline(0, 1)
plot(eb.diag, mash.diag); abline(0, 1)
str(objl[[ebnam]][[seed]])
mash.diag <- sqrt(diag(objl[[mashnam]][[seed]]$Sigl$U.1) * 
                    sum(objl[[mashnam]][[seed]]$omegaseq * objl[[mashnam]][[seed]]$pimat[, "U.1"] / sum(objl[[mashnam]][[seed]]$pimat[, "U.1"])))
eb.diag <- sqrt(diag(objl[[ebnam]][[seed]]$Sigl[[1]]) * sum(objl[[ebnam]][[seed]]$omegaseq * objl[[ebnam]][[seed]]$pimat[, 1]))


rownames(objl[[ebnam]][[seed]]$Sigl[[1]])[diag(objl[[ebnam]][[seed]]$Sigl[[1]]) > .015]
rownames(objl[[ebnam]][[seed]]$Sigl[[1]])[diag(objl[[ebnam]][[seed]]$Sigl[[1]]) < .015]

phsd <- apply(Yhat[truemuts, rownames(objl[[ebnam]][[seed]]$Sigl[[1]])], 2, function(v) sd(v, na.rm = T))
plot(diag(objl[[ebnam]][[seed]]$Sigl[[1]]), phsd); abline(0, 1)
plot(diag(objl[[ednam]][[seed]]$Sigl[[1]]), phsd); abline(0, 1)



# plot(diag(objl[[ebnam.eqtl]][[seed]]$Sighat), diag(objl[[mashnam.eqtl]][[seed]]$Sighat))
abline(0, 1)
mean(diag(objl[[ebnam]][[seed]]$Sighat) / diag(objl[[mashnam]][[seed]]$Sighat))
plot(diag(objl[[ebnam]][[seed]]$Sighat), diag(objl[[ednam]][[seed]]$Sighat))
abline(0, 1)
plot(diag(objl[[ebnam.eqtl]][[seed]]$Sighat), diag(objl[[ednam.eqtl]][[seed]]$Sighat))
abline(0, 1)
file.base.ed <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.ed, sapply(var.in.name.ed, get), sep = "_"), collapse = "_"))
bovy.output.file.namc <- gsub("justED", "mash", paste0(file.base.ed, "_bovy_output.RData"))
load(bovy.output.file.namc)
sams.for.testing <- setdiff(intersect(objl[[mashnam]][[seed]]$saml$sams.for.testing, objl[[ebnam]][[seed]]$saml$sams.for.testing),
                              rownames(mashdata.bovy$Bhat))
sams.for.testing <- objl[[mashnam]][[seed]]$saml$sams.for.model.fitting
# sams.for.testing <- objl[[ednam]][[seed]]$saml$sams.for.testing

# eb.scen <- 10
# file.base.seed <- gsub("XXX", seed, runtab[eb.scen, "file.base"])
# emout.file.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_emout.RData")
# emloaded <- load(file = emout.file.namc)
# dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
# npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
# file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
# suppressWarnings(load(file.in))
# emout.mix$Sigl <- lapply(emout.mix$Sigl, function(M){ dimnames(M) <- list(ph.use, ph.use); M})
# dimnames(emout.mix$Sig.mn) <- list(ph.use, ph.use)
# eb.objlc <- emout.mix
# eb.objlc$pimat <- eb.objlc$pi
eb.objlc <- objl[[ebnam]][[seed]]
mash.objlc <- objl[[mashnam]][[seed]]
ed.objlc <- objl[[ednam]][[seed]]
phord <- colnames(eb.objlc$Sigl[[1]])



# eb.Sigl <- lapply(eb.objlc$Sigl, function(M) M[phord, phord])
# ed.Sigl <- lapply(ed.objlc$Sigl, function(M) M[phord, phord])
# samlook <- sample(intersect(intersect(sams.for.testing, truemuts), rownames(res.store$mn)), 100)
samlook <- sample(intersect(sams.for.testing, truemuts), 100)
# samlook <- intersect(sams.for.testing, truemuts)
# samlook <- sample(intersect(sams.for.testing, negcons), 200)

# eb.Sigl <- lapply(emout.mix$Sigl, function(M) M[phord, phord])
# out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
#                                    Sigl = eb.Sigl, R = emout.mix$R[phord, phord], omegaseq = emout.mix$omegaseq,
#                                    pimat = emout.mix$pi, meth = "post.mn")
# eb.Sigl.adj <- lapply(eb.objlc$Sigl, function(M) M[phord, phord] * ((mash.diag / eb.diag) %o% (mash.diag / eb.diag)))
out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
                                   # Sigl = eb.Sigl.adj,
                                   Sigl = lapply(eb.objlc$Sigl, function(M) M[phord, phord]),
                                   R = eb.objlc$R[phord, phord], omegaseq = eb.objlc$omegaseq,
                                   pimat = eb.objlc$pimat, meth = "post.mn")
Y.em = Yhat[samlook, phord]; S.em = smat[samlook, phord]
# Sigl = eb.Sigl.adj,
Sigl = lapply(eb.objlc$Sigl, function(M) M[phord, phord]);
R = eb.objlc$R[phord, phord]; omegaseq = eb.objlc$omegaseq;
pimat = eb.objlc$pimat; meth = "just.obj"; prior.in.obj = T; bic.pen.mult = 0;
update.Sig = T; wish.pri = T; llmat = NULL; Ksig = eb.objlc$Ksig;
fac.model = c("fa", "pca")[2]
recalc.llmat = T
for(pri.in in c(T, F)){
  out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
                                     # Sigl = eb.Sigl.adj,
                                     Sigl = lapply(eb.objlc$Sigl, function(M) M[phord, phord]),
                                     R = eb.objlc$R[phord, phord], omegaseq = eb.objlc$omegaseq,
                                     pimat = eb.objlc$pimat, meth = "just.obj", prior.in.obj = pri.in, bic.pen.mult = 0,
                                     update.Sig = T, wish.pri = T, Ksig = eb.objlc$Ksig)
  mash.Sigl <- lapply(mash.objlc$Sigl[mashind], function(M) M[phord, phord] / 1)
  mash.Sigl$U.1 <- objl[[ednam]][[seed]]$Sigl[[1]]
  mash.out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
                                     Sigl = mash.Sigl, R = mash.objlc$R[phord, phord], omegaseq = mash.objlc$omegaseq,
                                     pimat = mash.objlc$pimat[, mashind, drop = F] / sum(mash.objlc$pimat[, mashind]), 
                                     meth = "just.obj", prior.in.obj = pri.in, bic.pen.mult = 0,
                                     update.Sig = T, wish.pri = T, Ksig = eb.objlc$Ksig)
  
  print(out.post.mix$obj)
  print(mash.out.post.mix$obj)
}
str(out.post.mix)
# colnames(pimat.mash)
pimat.mash <- objl[[mashnam]][[seed]]$pimat
pimat.eb <- objl[[ebnam]][[seed]]$pimat
mash.indsin <- 1
mashind <- colnames(pimat.mash)[order(-colSums(pimat.mash))][mash.indsin]
mashind <- setdiff(colnames(pimat.mash), phord)
# mashind <- paste0("U.", 1:8)
# mashind <- colnames(pimat.mash)[!grepl("U.", colnames(pimat.mash))]
# mashind <- setdiff(colnames(pimat.mash), phord)
# mashind <- setdiff(colnames(pimat.mash), paste0("U.", 1:3))
# mashind <- colnames(pimat.mash)
# mashind <- c("U.1", phord)
mashind <- c("U.1")
mashind <- intersect(mashind, colnames(pimat.mash)[colSums(pimat.mash) > 0])

mashind
mash.Sigl <- lapply(mash.objlc$Sigl[mashind], function(M) M[phord, phord] / 1)
mash.Sigl$U.1 <- objl[[ednam]][[seed]]$Sigl[[1]]
# mash.Sigl$U.1 <- objl[[ebnam]][[seed]]$Sigl[[1]] / objl[[ebnam]][[seed]]$Sigl[[1]][1, 1] * objl[[ednam]][[seed]]$Sigl[[1]][1, 1]
mash.out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
                                   Sigl = mash.Sigl, R = mash.objlc$R[phord, phord], omegaseq = mash.objlc$omegaseq,
                                   pimat = mash.objlc$pimat[, mashind, drop = F] / sum(mash.objlc$pimat[, mashind]), meth = "post.mn")
ed.out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
                                        Sigl = lapply(ed.objlc$Sigl, function(M) M[phord, phord]), 
                                        R = ed.objlc$R[phord, phord], omegaseq = ed.objlc$omegaseq,
                                        pimat = ed.objlc$pimat, meth = "post.mn")
eblik <- out.post.mix$loglikv[samlook]
edlik <- ed.out.post.mix$loglikv[samlook]
mashlik <- mash.out.post.mix$loglikv[samlook]
par(mfrow = c(2, 2))
plot(eblik, mashlik)
abline(0, 1, col = 2)
plot(eblik, edlik)
abline(0, 1, col = 2)
mean(eblik)
mean(edlik)
mean(mashlik)


ebnam <- "impc_em.fit_nSig_1_fm_pca_bic_0_mash_FALSE"#"impc_em.fit_nSig_1_fm_pca_bic_0_mash_TRUE"#impc_em.fit_nSig_1_fm_fa"
eb.objlc <- objl[[ebnam]][[seed]]

Sigc <- eb.objlc$Sighat
Rc <- eb.objlc$R

Pc <- 150
Rc <- rWishart(1, Pc, diag(rep(1, Pc)))[, , 1]
Sigc <- rWishart(1, Pc, diag(rep(1, Pc)))[, , 1]

str(Rc)

nits <- 200
system.time({
for(i in 1:nits)
  test <- solve(Sigc + Rc)
})
nits <- 200
system.time({
  for(i in 1:nits)
    test <- spdinv(Sigc + Rc)
})
K <- 20
Sig.low <- eigen(Sigc)$vectors[, 1:K]
R.low <- eigen(Rc)$vectors[, 1:K]
(R.low %*% t(R.low))[1:3, 1:3]
(R.low[1:3, ] %*% t(R.low[1:3, ]))[1:3, 1:3]


system.time({
  for(i in 1:nits){
    test1 <- Sig.low %*% solve(diag(rep(1, K)) + t(Sig.low) %*% Sig.low, t(Sig.low))
    test2 <- test1 + test1 %*% R.low %*% solve(diag(rep(1, K)) + t(R.low) %*% test1 %*% R.low, t(R.low))
  }
})





spl <- smooth.spline((mashlik + eblik) / 2, mashlik - eblik)
plot((mashlik + eblik) / 2, mashlik - eblik)
lines(spl$x, spl$y, col = 2)
abline(h = 0)
apply(mash.out.post.mix$rmat, 2:3, mean)
apply(out.post.mix$rmat, 2:3, mean)


apply(ed.out.post.mix$rmat, 2:3, mean)




hist(diag(eb.Sigl[[1]]))

par(mfrow = c(2, 2))
plot(diag(eb.Sigl[[1]]) / max(diag(eb.Sigl[[1]])), diag(mash.Sigl$U.1), xlim = c(0, 1), ylim = c(0, 1), xlab = "EB", ylab = "mash U.1")
abline(0, 1)
xpl <- diag(objl[[ebnam]][[seed]]$Sighat)
ypl <- diag(objl[[mashnam]][[seed]]$Sighat)
plot(xpl, ypl, xlim = c(0, max(xpl)), ylim = c(0, max(ypl)), xlab = "EB sighat", ylab = "mash sighat") 
abline(0, 1)
xpl <- diag(objl[[ebnam]][[seed]]$Sighat)
ypl <- diag(objl[[ednam]][[seed]]$Sigl[[1]])
plot(xpl, ypl, xlim = c(0, max(xpl)), ylim = c(0, max(ypl)), xlab = "EB sighat", ylab = "XD Sig")
abline(0, 1)


image(eb.Sigl[[1]], col = rain)
image(mash.Sigl$U.1, col = rain)

ordind <- 12
samdiff1 <- names(mashlik)[order(mashlik - eblik, decreasing = T)[ordind]]
samdiff2 <- names(mashlik)[order(eblik - mashlik, decreasing = T)[ordind]]
plot(Yhat[samdiff1, ] - compl[[ebnam]]$mnarr[samdiff1, , seed], Yhat[samdiff1, ] - compl[[mashnam]]$mnarr[samdiff1, , seed])
abline(0, 1)
plot(compl[[ebnam]]$mnarr[samdiff1, , seed], compl[[mashnam]]$mnarr[samdiff1, , seed])
abline(0, 1)



plot(Yhat[samdiff2, ] - compl[[ebnam]]$mnarr[samdiff2, , seed], Yhat[samdiff2, ] - compl[[mashnam]]$mnarr[samdiff2, , seed])
abline(0, 1)
plot(compl[[ebnam]]$mnarr[samdiff2, , seed], compl[[mashnam]]$mnarr[samdiff2, , seed])
abline(0, 1)


hist(apply(Yhat, 2, function(v) mad(v, na.rm = T)))
boxplot(Yhat)

, 2, function(v) sd(v, na.rm = T)))

matplot(objlc$omegaseq %*% t(diag(objlc$Sigl[[1]])[1:20]), objlc$pimat, ty = "l", log = "x")
plot(Yhat[samdiff2, ])

plot(objlc$omegaseq %*% t(mean(diag(objlc$Sigl[[1]]))), objlc$pimat, ty = "l", log = "x")
lines(mash.objlc$omegaseq %*% t(mean(diag(mash.objlc$Sigl[["U.1"]]))), mash.objlc$pimat[, "U.1"], col = 2)

mashlik[samdiff1]
eblik[samdiff1]

mashlik[samdiff2]
eblik[samdiff2]


sing.om.dist <- rowSums(mash.objlc$pimat[, phord])
sing.om.dist <- sing.om.dist / sum(sing.om.dist)
U.om.dist <- rowSums(mash.objlc$pimat[, grep("U", colnames(pimat.mash))])
U.om.dist <- U.om.dist / sum(U.om.dist)
plot(mash.objlc$omegaseq, sing.om.dist, log = "x", ty = "l")
lines(mash.objlc$omegaseq, U.om.dist)


colMeans(is.na(Yhat))
mean(is.na(Yhat))



mean(compl$impc_em.fit_nSig_1$llmat.raw[samlook, splitc])
mean(out.post.mix$loglikv[samlook])

par(mfrow = c(3, 3))
for(j in 1:9){
  # try({
    pimat.mash <- objl[[mashnam]][[seed]]$pimat
    pimat.eb <- objl[[ebnam]][[seed]]$pimat
    mashind <- colnames(pimat.mash)[order(-colSums(pimat.mash))][j]
    maxsig.mash <- max(diag(objl[[mashnam]][[seed]]$Sigl[[mashind]]))
    maxsig.eb <- max(diag(objl[[ebnam]][[seed]]$Sigl[[1]]))
    norm.omega.mash <- objl[[mashnam]][[seed]]$omegaseq * ifelse(maxsig.mash == 0, 1, maxsig.mash)
    norm.omega.eb <- objl[[ebnam]][[seed]]$omegaseq * ifelse(maxsig.eb == 0, 1, maxsig.eb)
    norm.pi.mash <- objl[[mashnam]][[seed]]$pimat[, mashind]
    norm.pi.mash <- norm.pi.mash / sum(norm.pi.mash, na.rm = T)
    plot(norm.omega.mash, norm.pi.mash, log = "x", ty = "l", main = paste(mashind, signif(colSums(pimat.mash)[mashind], 2)))
    lines(norm.omega.eb, pimat.eb, col = 2)
  # })
}

str(objl$impc_mash_nSig_1[[seed]]$Sigl)

names(objl)
str(resl.comp, m = 1)

str(res.mash.fitted.model, m = 2)
str(res.mash.all.testing, m = 2)
str(res.mash.all.testing$vloglik, m = 2)
str(Sigll)
str(compl, m = 2)

llmean.splitmat.raw <- sapply(compl, function(x) colMeans(x$llmat, na.rm = T))

par(mfrow = c(3, 3))
for(splitc in 1:9){
  plot(compl$impc_em.fit_nSig_1$llmat.raw[truemuts, splitc], compl$impc_mash_nSig_1$llmat.raw[truemuts, splitc])
  abline(0, 1, col = 2)
  mean(compl$impc_em.fit_nSig_1$llmat.raw[truemuts, splitc], na.rm = T)
  mean(compl$impc_mash_nSig_1$llmat.raw[truemuts, splitc], na.rm = T)
}
table(is.na(compl$impc_mash_nSig_1$llmat.raw[truemuts, splitc]), is.na(compl$impc_em.fit_nSig_1$llmat.raw[truemuts, splitc]))

table(compl$impc_mash_nSig_1$llmat.raw[truemuts, splitc] > compl$impc_em.fit_nSig_1$llmat.raw[truemuts, splitc])



plot(compl$impc_em.fit_nSig_1$llmat.zero[, 1], compl$impc_mash_nSig_1$llmat.zero[, 1])
abline(0, 1)
mean(compl$impc_em.fit_nSig_1$llmat.zero[, 1], na.rm = T)
mean(compl$impc_mash_nSig_1$llmat.zero[, 1], na.rm = T)

plot(res.mash.all.testing$vloglik[match(truemuts, rownames(res.mash.all.testing$vloglik)), 1], compl$impc_mash_nSig_1$llmat.zero[truemuts, splitc])

names(compl)

llmean.splitmat.look <- sapply(compl, function(x) x$llmat.raw[1:10, ])

colMeans(llmean.splitmat.look, na.rm = T)


llmean.splitmat.zero <- sapply(compl, function(x) colMeans(x$llmat.zero, na.rm = T))
llmean.splitmat.raw <- sapply(compl, function(x) colMeans(x$llmat.raw, na.rm = T))

mn.zero <- colMeans(llmean.splitmat.zero)
se.zero <- apply(llmean.splitmat.zero, 2, function(v) sd(v) / sqrt(length(v)))
mn.raw <- colMeans(llmean.splitmat.raw)
se.raw <- apply(llmean.splitmat.raw, 2, function(v) sd(v) / sqrt(length(v)))

formc <- "f"
ndig <- 1
cvlikdat <- data.frame(meth = names(mn.raw), 
         mn.se.zero = paste0(formatC(mn.zero, digits = ndig, format = formc), " (", formatC(se.zero, digits = ndig, format = formc), ")"),
         mn.se.raw = paste0(formatC(mn.raw, digits = ndig, format = formc), " (", formatC(se.raw, digits = ndig, format = formc), ")"))
cvlikdat
#


for(i in 1:length(default.parameters[[Data]]))
  assign(names(default.parameters[[Data]][i]), default.parameters[[Data]][[i]])

seed <- 4
nSig <- 1
K <- 74
nfit <- 250#length(cv.sams.curr)
si <- 1
bo <- 1
EDtol <- 1e-5
EDmeth <- "justED"
  
  
dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
suppressWarnings(load(file.in))

source(paste0(R.file.dir, "/impc_mv_paper_code/EM_fns_lean.R"))
load(file = uv.results.Y.S)
cv.sams.curr <- intersect(linemap$geno[linemap$line.type == "trueMut"], names(compl$impc_eb_1$llmat[, seed])[!is.na(compl$impc_eb_1$llmat[, seed])])


meth <- "eb"
sexspecific <- F
var.in.namec <- c(var.in.name.old, "sexspecific")
file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.namec, sapply(var.in.namec, get), sep = "_"), collapse = "_"))
# file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name, sapply(var.in.name, get), sep = "_"), collapse = "_"))
# resl.file.namc <- paste0(file.base, "_resl.RData")
emout.file.namc <- paste0(file.base, "_emout.RData")
res.store.namc <- paste0(file.base, "_res.RData")
loocv.res.store.namc <- paste0(file.base, "_loocv_res.RData")
fac.res.store.namc <- paste0(file.base, "_facres.RData")
load(file = emout.file.namc)


Y.em <- Yhat
S.em <- smat
# Y.em <- Y
# S.em <- S
# S.em[S.em == 10] <- 5
cv.sams.use <- sample(cv.sams.curr, nfit)
em.pi <- emout.mix$pi
# em.pi[] <- 1
em.pi <- em.pi / sum(em.pi)
out.post.mix.eb.obj <- em.update.function(Y.em = Y.em[cv.sams.use, ph.use], S.em = S.em[cv.sams.use, ph.use],
                                          Sigl = lapply(emout.mix$Sigl, function(M) M[ph.use, ph.use]), R = emout.mix$R[ph.use, ph.use],
                                          omegaseq = emout.mix$omegaseq, prior.in.obj = F,
                                          pimat = em.pi, meth = "just.obj")




var.in.name.ed <- setdiff(var.in.name.ed, "sexspecific")
# var.in.name.ed <- c(var.in.name.ed, "sexspecific")
EDtol <- 1e-4

file.base.ed <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.ed, sapply(var.in.name.ed, get), sep = "_"), collapse = "_"))
bovy.resl.file.namc <- paste0(file.base.ed, "_bovy_resl.RData")
load(file = bovy.resl.file.namc)
res.store <- resl.store$raw
phnam.ed <- colnames(res.store$mn)
ed.Sigl <- res.store$Sigl
out.post.mix.ed.obj <- em.update.function(Y.em = Y.em[cv.sams.use, ph.use], S.em = S.em[cv.sams.use, ph.use],
                                          Sigl = lapply(ed.Sigl, function(M) M[ph.use, ph.use]), R = emout.mix$R[ph.use, ph.use],
                                          omegaseq = res.store$omegaseq, prior.in.obj = F,
                                          pimat = res.store$pimat, meth = "just.obj")



EDtol <- 1e-5
file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
mash.raw.results.file.namc <- paste0(file.base.mash, "_mash_raw_results.RData")
load(mash.raw.results.file.namc)
phnam.mash <- colnames(res.mash.fitted.model$result$PosteriorMean)
mash.omegaseq <- res.mash.fitted.model$fitted_g$grid^2
mash.n.om <- length(mash.omegaseq)
null.Sigl <- list(diag(rep(0, P)))
names(null.Sigl) <- "null"samdiff1 <- which.max(mashlik - eblik)

mash.Sigl <- lapply(c(null.Sigl, res.mash.fitted.model$fitted_g$Ulist), function(M){ dimnames(M) <- list(phnam.mash, phnam.mash); M})
mash.nSig <- length(mash.Sigl)
mash.pimat.t <- matrix(0, mash.nSig, mash.n.om, dimnames = list(c("null", names(res.mash.fitted.model$fitted_g$Ulist)), as.character(mash.omegaseq)))
mash.piv <- res.mash.fitted.model$fitted_g$pi
mash.pimat.t[2:mash.nSig, ] <- mash.piv[2:length(mash.piv)]
mash.pimat.t[1, 1] <- mash.piv[1]
signam.keep <- names(mash.Sigl)#[!names(mash.Sigl) %in% ph.use] 
mash.Sigl.use <- mash.Sigl[signam.keep]
mash.pimat.t.use <- mash.pimat.t[signam.keep, , drop = F]
mash.pimat.t.use <- mash.pimat.t.use / sum(mash.pimat.t.use)
out.post.mix.mash.obj <- em.update.function(Y.em = Y.em[cv.sams.use, ph.use], S.em = S.em[cv.sams.use, ph.use],
                                            Sigl = lapply(mash.Sigl.use, function(M) M[ph.use, ph.use]), R = emout.mix$R[ph.use, ph.use],
                                            omegaseq = mash.omegaseq, prior.in.obj = F,
                                            pimat = t(mash.pimat.t.use), meth = "just.obj")

names(out.post.mix.eb.obj$llikv) <- names(out.post.mix.ed.obj$llikv) <- names(out.post.mix.mash.obj$llikv) <- cv.sams.use
# names(out.post.mix.eb.obj$llikv) <- names(out.post.mix.mash.obj$llikv) <- cv.sams.use
mean(out.post.mix.eb.obj$llikv[cv.sams.use], na.rm = T)
mean(out.post.mix.ed.obj$llikv[cv.sams.use], na.rm = T)
mean(out.post.mix.mash.obj$llikv[cv.sams.use], na.rm = T)


plot(out.post.mix.eb.obj$llikv[cv.sams.use], out.post.mix.mash.obj$llikv[cv.sams.use])
abline(0, 1)
plot(out.post.mix.eb.obj$llikv[cv.sams.use], out.post.mix.ed.obj$llikv[cv.sams.use])
abline(0, 1)
var.in.name <- c(var.in.name.old, "sexspecific")

load(file = res.store.namc)

cv.sams.curr <- intersect(linemap.in$geno[linemap.in$line.type == "trueMut"], names(compl$eb$llmat[, seed])[!is.na(compl$eb$llmat[, seed])])


nfit <- 500#length(cv.sams.curr)
source(paste0(R.file.dir, "/impc_mv_paper_code/EM_fns_lean.R"))
load(file = uv.results.Y.S)
Y.em <- Yhat
S.em <- smat
# Y.em <- Y
# S.em <- S
# S.em[S.em == 10] <- 5

cv.sams.use <- sample(cv.sams.curr, nfit)
em.pi <- emout.mix$pi
# em.pi[] <- 1
em.pi <- em.pi / sum(em.pi)
out.post.mix.eb.obj <- em.update.function(Y.em = Y.em[cv.sams.use, ph.use], S.em = S.em[cv.sams.use, ph.use],
                                          Sigl = lapply(emout.mix$Sigl, function(M) M[ph.use, ph.use]), R = emout.mix$R[ph.use, ph.use],
                                          omegaseq = emout.mix$omegaseq, prior.in.obj = F,
                                          pimat = em.pi, meth = "just.obj")





graphics.off()
par(mfrow = c(2, 2))
siglook.eb <- lapply(emout.mix$Sigl, function(M) M[ph.use, ph.use])[[1]]
siglook.mash <- lapply(mash.Sigl.use, function(M) M[ph.use, ph.use])[[1]]
image(siglook.eb)
image(siglook.mash)
plot(siglook.eb, siglook.mash)


graphics.off()
par(mfrow = c(3, 3))
for(indlook in 1:9){
  ylimc <- range(c(emout.mix$pi), c(mash.pimat.t.use))
  xlimc <- range(c(emout.mix$omegaseq * siglook.eb[indlook, indlook], mash.omegaseq * siglook.mash[indlook, indlook]))
  plot(emout.mix$omegaseq * siglook.eb[indlook, indlook], emout.mix$pi,  log = "x", xlim = xlimc, ylim = ylimc, ty = "l")
  lines(mash.omegaseq * siglook.mash[indlook, indlook], c(mash.pimat.t.use), col = 2)
}

sum(emout.mix$pi)
sum(c(mash.pimat.t.use))
mash.pidat[order(-mash.pidat$pitot)[1:30], ]


plot(out.post.mix.eb.obj$llikv, out.post.mix.mash.obj$llikv)
abline(0, 1)
plot(out.post.mix.eb.obj$llikv, out.post.mix.ed.obj$llikv)
abline(0, 1)


plot(res.mash.all.testing$vloglik[cv.sams.use, ], out.post.mix.mash.obj$llikv)
abline(0, 1)





sum(compl$mash$llmat[cv.sams.use, seed], na.rm = T)

plot(out.post.mix.eb.obj$llikv[cv.sams.use], compl$mash$llmat[cv.sams.use, seed])
abline(0, 1)

methdiff <- out.post.mix.eb.obj$llikv[cv.sams.use] - compl$mash$llmat[cv.sams.use, seed]
plot(methdiff, rowSums(Y[cv.sams.use, ] == 0))
samlook <- cv.sams.use[order(-1 * (out.post.mix.eb.obj$llikv[cv.sams.use] - compl$mash$llmat[cv.sams.use, seed]))[3]]
plot(Y[samlook, ])

plot(out.post.mix.eb.obj$llikv, compl$ed$llmat[cv.sams.use, seed])
abline(0, 1)


obj.out.curr <- em.update.function(Y.em = Y[sams.for.lik.cross.val, ph.use], S.em = S[sams.for.lik.cross.val, ph.use],
                                   Sigl = Sigll[[seed]], R = emout.mix$R, omegaseq = emout.mix$omegaseq,
                                   pimat = emout.mix$pi, meth = "just.obj", prior.in.obj = F)



str(compl)

#

plot(compl$ed$llmat[true.curr, 1], compl$eb$llmat[true.curr, 1])

plot(compl$mash$llmat[true.curr, 1], compl$eb$llmat[true.curr, 1])
abline(0, 1)
#


true.use <- unique(resimp[resimp$line.type == "trueMut", "geno"])
table(sign(resl.comp$eb.fac$mn)[rownames(resl.comp$eb.fac$mn) %in% true.use, "f.1"])
table(sign(resl.out$eb.fac$mn)[rownames(resl.out$eb.fac$mn) %in% true.use, "f.1"])





graphics.off()
matplot(mash.pimat.t, ty = "l")
mash.pidat <- data.frame(signam = rownames(mash.pimat.t), pitot = rowSums(mash.pimat.t))
mash.pidat[order(-mash.pidat$pitot)[1:30], ]



# sapply(compl, function(x) x$crossval.llv)
# true.curr <- linemap.in$geno[linemap.in$line.type == "trueMut"]
# plot(compl$mash$llmat[true.curr, 1], compl$eb$llmat[true.curr, 1])
# abline(0, 1)
# 
# colSums(compl$eb$llmat[true.curr, ], na.rm = T)
# colSums(compl$mash$llmat[true.curr, ], na.rm = T)
# colSums(compl$ed$llmat[true.curr, ], na.rm = T)
# 
# plot(colSums(compl$eb$llmat[true.curr, ], na.rm = T) - colSums(compl$ed$llmat[true.curr, ], na.rm = T), 
#      compl$eb$crossval.llv - compl$ed$crossval.llv)
# mean(colSums(compl$eb$llmat[true.curr, ], na.rm = T) - colSums(compl$ed$llmat[true.curr, ], na.rm = T)) - 
#   mean(compl$eb$crossval.llv - compl$ed$crossval.llv)
# compl$eb$crossval.llv - compl$mash$crossval.llv
# compl$eb$crossval.llv
# compl$mash$crossval.llv
# compl$ed$crossval.llv
# 

