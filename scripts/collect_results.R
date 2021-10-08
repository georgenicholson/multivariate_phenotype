rm(list = ls())
# Data <- c("eqtl", "impc")[2]
xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
source(paste0(R.file.dir, "/impc_mv_paper_code/EM_fns_lean.R"))
load(file = file.runtab)
runtab
dimnaml <- list()
for(Data in c("impc", "eqtl")){
  dimnaml[[Data]] <- list()
  if(Data == "eqtl"){
    mash.in <- readRDS(file = paste0(R.file.dir, "/urbut_code/stephenslab-gtexresults-961b969/data/MatrixEQTLSumStats.Portable.ld2.Z.rds"))
    Y <- Yhat <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
    S <- smat <- Yhat / rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
  }
  if(Data == "impc"){
    load(file = uv.results.Y.S)
    linemap.in <- linemap
  }
  dimnaml[[Data]]$sam.names <- rownames(Yhat)
  dimnaml[[Data]]$meas.names <- colnames(Yhat)
  dimnaml[[Data]]$N.all <- length(dimnaml[[Data]]$sam.names)
  dimnaml[[Data]]$P.all <- length(dimnaml[[Data]]$meas.names)
}


resl.comp <- compl <- Sigll <- Ksigl <- pimatl <- Sighatl <- Rl <- omegaseql <- samll <- objl <- list()
load(file = uv.results.Y.S)
resl.comp$uv <- list(mn = Yhat[dimnaml$impc$sam.names, dimnaml$impc$meas.names], sd = smat[dimnaml$impc$sam.names, dimnaml$impc$meas.names])
fac.methv <- "varimax"
runtab
ntimes
for(scen in 1:nrow(runtab)){
  print(runtab[scen, c("N", "P", "nSig", "me", "meth", "da")])
  for(j in 1:ncol(runtab))
    assign(colnames(runtab)[j], runtab[scen, j])
  for(j in 1:length(dimnaml[[Data]]))
    assign(names(dimnaml[[Data]])[j], dimnaml[[Data]][[j]]) #sam.names, meas.names, n.all, P.all
  array1 <- array(NA, dim = c(N.all, P.all, ntimes), dimnames = list(sam.names, meas.names, 1:ntimes))
  array2 <- array(NA, dim = c(N.all, ntimes), dimnames = list(sam.names, 1:ntimes))
  array.fac <- array(NA, dim = c(N.all, nfac, ntimes), dimnames = list(sam.names, facnam, 1:ntimes))
  crossval.llv <- rep(NA, ntimes)
  shared.formatl <- list(mnarr = array1, sdarr = array1, loocv.mnarr = array1, loocv.sdarr = array1, 
                         # varimax = list(fac.mnarr = array.fac, fac.sdarr = array.fac), 
                         # promax = list(fac.mnarr = array.fac, fac.sdarr = array.fac), 
                         lfsrarr = array1, llmat = array2, llmat.raw = array2, llmat.zero = array2)
  if(meth == "eb"){
    namc <- paste(Data, me, "nSig", nSig, "fm", EMfm, "bic", EMbic, "EMK", EMK, sep = "_")
  } else {
    namc <- paste(Data, me, "nSig", nSig, sep = "_")
  }
  compl[[namc]] <- shared.formatl
  objl[[namc]] <- list()
  pimatl[[namc]] <- Ksigl[[namc]] <- Sigll[[namc]] <- Sighatl[[namc]] <- Rl[[namc]] <- omegaseql[[namc]] <- samll[[namc]] <- list()
  for(seed in 1:ntimes){
    objl[[namc]][[seed]] <- list()
    file.base.seed <- gsub("XXX", seed, file.base)
    dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
    npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
    data.subset.file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
    suppressWarnings(objvecc <- load(data.subset.file.in))
    samll[[namc]][[seed]] <- objl[[namc]][[seed]]$saml <- list(sams.for.cor.est = sams.for.cor.est, sams.for.model.fitting = sams.for.model.fitting, 
                                  sams.for.testing = sams.for.testing, sams.for.lik.cross.val = sams.for.lik.cross.val)
    if(meth == "eb"){
      emout.file.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_emout.RData")
      res.store.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_res.RData")
      loocv.res.store.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_loocv_res.RData")
      fac.res.store.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_facres.RData")
      emloaded <- resloaded <- loocv.loaded <- fac.loaded <- NULL
      if(file.exists(emout.file.namc))
        emloaded <- load(file = emout.file.namc)
      if(file.exists(res.store.namc))
        resloaded <- load(file = res.store.namc)
      if(file.exists(loocv.res.store.namc))
        loocv.loaded <- load(file = loocv.res.store.namc)
      if(file.exists(fac.res.store.namc))
        fac.loaded <- load(file = fac.res.store.namc)
      if(is.null(emloaded) | is.null(resloaded))
        warning(paste0(namc, " seed ", seed, " not loaded"))
      res.store <- resl.store$raw
      Ksigl[[namc]][[seed]] <- objl[[namc]][[seed]]$Ksig <- emout.mix$Ksig
      if(me == "em.fit" & nSig == 1 & Data == "impc"){
        if(!is.null(loocv.loaded)){
          compl[[namc]]$loocv.mnarr[match(rownames(loocv.store$mn), sam.names), 
                                    match(colnames(loocv.store$mn), meas.names), seed] <- loocv.store$mn
          compl[[namc]]$loocv.sdarr[match(rownames(loocv.store$sd), sam.names), 
                                    match(colnames(loocv.store$sd), meas.names), seed] <- loocv.store$sd
        }
      }
    }
    if(meth == "mash"){
      mash.resl.file.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_mash_resl.RData")
      mashloaded <- NULL
      mashloaded <- load(file = mash.resl.file.namc)
      datatypec <- switch(Data, eqtl = "raw", impc = "zero")
      if(datatypec  %in% names(resl.store)){
        res.store <- resl.store[[datatypec]]
      } else {
        res.store <- resl.store[[1]]
        warning("Need to re-run mash on this data set")
      }
      if(is.null(mashloaded))
        warning(paste0(namc, " seed ", seed, " not loaded"))
    }
    if(meth == "XD"){
      bovy.resl.file.namc <- paste0(meth.comp.output.dir, "/", file.base.seed, "_bovy_resl.RData")
      xdloaded <- NULL
      xdloaded <- load(file = bovy.resl.file.namc)
      datatypec <- switch(Data, eqtl = "raw", impc = "zero")
      res.store <- resl.store[[datatypec]]
      if(is.null(xdloaded))
        warning(paste0(namc, " seed ", seed, " not loaded"))
    }
    Sigll[[namc]][[seed]] <- objl[[namc]][[seed]]$Sigl <- lapply(res.store$Sigl, function(M) M[match(meas.names, rownames(M)), match(meas.names, colnames(M))])
    Rl[[namc]][[seed]] <- objl[[namc]][[seed]]$R <- res.store$R[match(meas.names, rownames(res.store$R)), match(meas.names, colnames(res.store$R))]
    pimatl[[namc]][[seed]] <- objl[[namc]][[seed]]$pimat <- res.store$pimat
    omegaseql[[namc]][[seed]] <- objl[[namc]][[seed]]$omegaseq <- res.store$omegaseq
    Sighatl[[namc]][[seed]] <- objl[[namc]][[seed]]$Sighat <- 0
    for(sc in 1:length(Sigll[[namc]][[seed]])){
      Sighatl[[namc]][[seed]] <- Sighatl[[namc]][[seed]] + sum(pimatl[[namc]][[seed]][, sc] * omegaseql[[namc]][[seed]]) *  Sigll[[namc]][[seed]][[sc]]
      objl[[namc]][[seed]]$Sighat <- objl[[namc]][[seed]]$Sighat + sum(pimatl[[namc]][[seed]][, sc] * omegaseql[[namc]][[seed]]) *  Sigll[[namc]][[seed]][[sc]]
    }
    compl[[namc]]$mnarr[match(rownames(res.store$mn), sam.names), match(colnames(res.store$mn), meas.names), seed] <- res.store$mn
    compl[[namc]]$sdarr[match(rownames(res.store$sd), sam.names), match(colnames(res.store$sd), meas.names), seed] <- res.store$sd
    compl[[namc]]$lfsrarr[match(rownames(res.store$lfsr), sam.names), match(colnames(res.store$lfsr), meas.names), seed] <- res.store$lfsr
    compl[[namc]]$llmat[match(names(res.store$loglikv), sam.names), seed] <- res.store$loglikv
    compl[[namc]]$llmat.raw[match(names(resl.store$raw$loglikv), sam.names), seed] <- resl.store$raw$loglikv
    compl[[namc]]$llmat.zero[match(names(resl.store$zero$loglikv), sam.names), seed] <- resl.store$zero$loglikv
  }
  suppressWarnings(lmat.norm <- exp(compl[[namc]]$llmat - apply(compl[[namc]]$llmat, 1, function(v) max(v, na.rm = T))))
  pmix <- compl[[namc]]$pmix <- lmat.norm / rowSums(lmat.norm, na.rm = T)
  Sigmix.p <- colSums(pmix, na.rm = T) / sum(pmix, na.rm = T)
  compl[[namc]]$Sig.comb <- compl[[namc]]$R.comb <- 0
  for(seed in 1:ntimes){
    compl[[namc]]$Sig.comb <- compl[[namc]]$Sig.comb + Sighatl[[namc]][[seed]][meas.names, meas.names] * Sigmix.p[seed]
    compl[[namc]]$R.comb <- compl[[namc]]$R.comb + Rl[[namc]][[seed]][meas.names, meas.names] * Sigmix.p[seed]
  }
  resl.comp[[namc]]$Sig.comb <- compl[[namc]]$Sig.comb
  resl.comp[[namc]]$Sigcor.comb <- t(resl.comp[[namc]]$Sig.comb / sqrt(diag(resl.comp[[namc]]$Sig.comb))) / sqrt(diag(resl.comp[[namc]]$Sig.comb))
  eigc <- eigen(resl.comp[[namc]]$Sigcor.comb)
  resl.comp[[namc]]$facs.varimax <- varimax(eigc$vectors[, 1:nfac])$loadings
  resl.comp[[namc]]$facs.promax <- promax(eigc$vectors[, 1:nfac])$loadings
  dimnames(resl.comp[[namc]]$facs.varimax) <- dimnames(resl.comp[[namc]]$facs.promax) <- list(meas.names, facnam)
  wt.mnarr <- apply(sweep(compl[[namc]]$mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  wt.varr <- apply(sweep(compl[[namc]]$sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  wt.mnsqarr <- apply(sweep(compl[[namc]]$mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  wt.lfsrarr <- apply(sweep(compl[[namc]]$lfsrarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
  compl[[namc]]$mn.comb <- wt.mnarr
  compl[[namc]]$sd.comb <- sqrt(wt.varr + wt.mnsqarr - wt.mnarr^2)
  compl[[namc]]$lfsr.comb <- wt.lfsrarr
  if(me == "em.fit" & nSig == 1 & Data == "impc"){
    loocv.wt.mnarr <- apply(sweep(compl[[namc]]$loocv.mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
    loocv.wt.varr <- apply(sweep(compl[[namc]]$loocv.sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
    loocv.wt.mnsqarr <- apply(sweep(compl[[namc]]$loocv.mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
    compl[[namc]]$loocv.mn.comb <- loocv.wt.mnarr
    compl[[namc]]$loocv.sd.comb <- sqrt(loocv.wt.varr + loocv.wt.mnsqarr - loocv.wt.mnarr^2)
  }
  resl.comp[[namc]] <- c(resl.comp[[namc]], list(mn = compl[[namc]]$mn.comb[sam.names, meas.names], sd = compl[[namc]]$sd.comb[sam.names, meas.names], 
                            loocv.mn = compl[[namc]]$loocv.mn.comb[sam.names, meas.names], loocv.sd = compl[[namc]]$loocv.sd.comb[sam.names, meas.names],
                            lfsr = compl[[namc]]$lfsr.comb[sam.names, meas.names]))
}

save(resl.comp, file = file.resl.comp, version = 2)
save(compl, file = file.compl, version = 2)
save(objl, file = file.objl, version = 2)


resl.comp.fac <- list()
Data <- "impc"
namc <- mv.meth.nam.use
fac.meth <- "varimax"#fac.methv[1]
ntimes <- dim(compl[[namc]]$mnarr)[3]
ncore <- min(20, ntimes)
require(doParallel)
if(!"clust" %in% ls())
  clust <- makeCluster(rep("localhost", ncore), type = "SOCK")
registerDoParallel(clust)
fac.res.store <- foreach(seed = 1:ntimes, .verbose = T) %dopar% {
  sams.for.testing <- objl[[namc]][[seed]]$saml$sams.for.testing
  meas.names <- rownames(objl[[namc]][[seed]]$Sigl[[1]])
  facs <- switch(fac.meth, varimax = resl.comp[[namc]]$facs.varimax, promax = resl.comp[[namc]]$facs.promax)
  fac.out.post.mix <- em.update.function(Y.em = Yhat[sams.for.testing, meas.names], 
                                         S.em = smat[sams.for.testing, meas.names],
                                         Sigl = objl[[namc]][[seed]]$Sigl, R = objl[[namc]][[seed]]$R, 
                                         omegaseq = objl[[namc]][[seed]]$omegaseq,
                                         pimat = objl[[namc]][[seed]]$pimat, meth = "post.mn.fac", loadings = facs[meas.names, facnam])
  out <- list(mn = fac.out.post.mix$mnmat[sams.for.testing, facnam], sd = fac.out.post.mix$sdmat[sams.for.testing, facnam],
              loglikv = fac.out.post.mix$loglikv[sams.for.testing],
              lfsr = fac.out.post.mix$lfsrmat[sams.for.testing, facnam], 
              loadings = facs)
  return(out)
}

suppressWarnings(lmat.norm <- exp(compl[[namc]]$llmat - apply(compl[[namc]]$llmat, 1, function(v) max(v, na.rm = T))))
pmix <- lmat.norm / rowSums(lmat.norm, na.rm = T)
fac.mnarr <- fac.sdarr <- fac.lfsrarr <- fac.wt.mnarr <- fac.wt.varr <- fac.wt.mnsqarr <- fac.wt.lfsrarr <- 
  array(NA, dim = c(dimnaml[[Data]]$N.all, nfac, ntimes), dimnames = list(dimnaml[[Data]]$sam.names, facnam, 1:ntimes))
for(seed in 1:ntimes){
  fac.mnarr[samll[[namc]][[seed]]$sams.for.testing, , seed] <- fac.res.store[[seed]]$mn
  fac.sdarr[samll[[namc]][[seed]]$sams.for.testing, , seed] <- fac.res.store[[seed]]$sd
  fac.lfsrarr[samll[[namc]][[seed]]$sams.for.testing, , seed] <- fac.res.store[[seed]]$lfsr
}
fac.wt.mnarr <- apply(sweep(fac.mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
fac.wt.varr <- apply(sweep(fac.sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
fac.wt.mnsqarr <- apply(sweep(fac.mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
fac.wt.lfsrarr <- apply(sweep(fac.lfsrarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
resl.comp.fac[[namc]] <- list(mn = fac.wt.mnarr, sd = sqrt(fac.wt.varr + fac.wt.mnsqarr - fac.wt.mnarr^2), lfsr = fac.wt.lfsrarr)
save(resl.comp.fac, file = file.resl.comp.fac)


str(resl.comp.fac, m = 2)
#


















# #Comparing likelihoods 
# load(file = emout.file.namc)
# sub.sams.for.testing <- sams.for.lik.cross.val#[1:100]
# out.post.mix <- em.update.function(Y.em = Y[sub.sams.for.testing, ph.use], S.em = S[sub.sams.for.testing, ph.use],
#                                    Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
#                                    pimat = emout.mix$pi, meth = "post.mn")
# str(out.post.mix)
# eblik <- out.post.mix$loglikv[sub.sams.for.testing]
# mashlik <- compl$mash$llmat[sub.sams.for.testing, seed]
# plot(eblik, mashlik)
# abline(0, 1)
# 
# par(mfrow = c(4, 4))
# for(j in 1:10){
#   geno.look <- names(eblik)[order((eblik - mashlik))[j]]
#   plot(Yhat[geno.look, ])
#   abline(h = 0)
# }
# 
# 
# mean(eblik > mashlik)
# eblik[geno.look]
# mashlik[geno.look]
# exp(mean(mashlik) - mean(eblik))
# 
# sum(out.post.mix$loglikv[sub.sams.for.testing])
# sum(compl$mash$llmat[sub.sams.for.testing, seed])
# obj.out <- em.update.function(Y.em = Y[sams.for.lik.cross.val, ph.use], S.em = S[sams.for.lik.cross.val, ph.use],
#                               Sigl = emout.mix$Sigl, R = Rtest, omegaseq = emout.mix$omegaseq,
#                               pimat = emout.mix$pi, meth = "just.obj", prior.in.obj = F)
# 
# }













# 
# resl.out <- list()
# resimp <- data.frame()
# if(Data == "impc"){
#   ###############################################################################################
#   # add in estimate of correlation matrix from synthetic null lines
#   # sams.for.cor.est <- linemap.in[linemap.in$line.type == "negCon", "geno"]
#   # resl.comp$impc_eb_1$R <- cor((Y / S)[sams.for.cor.est, ph.use], use = "p", meth = "p")
#   # if("impc_eb.ss_1" %in% names(resl.comp))
#   #   resl.comp$impc_eb.ss_1$R <- cor((Y.genosex / S.genosex)[sams.for.cor.est, ph.use], use = "p", meth = "p")
#   # 
#   # ##################################################################
#   # # Calculate loadings from eigendecomposition of Sig to feed into EM_run
#   # Sig <- resl.comp$impc_eb_1$Sig.comb[ph.use, ph.use]
#   # R <- resl.comp$impc_eb_1$R[ph.use, ph.use]
#   # Sigdiag <- diag(Sig)
#   # Sigcor <- t(Sig / sqrt(Sigdiag)) / sqrt(Sigdiag)
#   # eigc <- eigen(Sigcor)
#   # facs.varimax <- varimax(eigen(Sigcor)$vectors[, 1:nfac])$loadings
#   # facs.promax <- promax(eigen(Sigcor)$vectors[, 1:nfac])$loadings
#   # dimnames(facs.promax) <- dimnames(facs.varimax) <- list(rownames(Sigcor), facnam)
#   # facs <- facs.varimax
#   # if(n.seed.run.subsams == 50)
#   #   save(facs.varimax, facs.promax, facs, file = file.glob.loadings, version = 2)
#   linemap <- linemap.in
#   linemap$line.type <- ifelse(linemap.in$line.type == "trueMut", "trueMutTes", "negConTes")
#   linemap <- linemap[linemap$geno %in% rownames(resl.comp$impc_eb_1$mn), ]
#   resl.comp <- resl.comp[!names(resl.comp) %in% c("varimax", "promax")]
#   out.perm <- err.rate.control(resl = resl.comp[!names(resl.comp) %in% c("uv.ss")], 
#                                err.rate.meth = "perm", sep.imp.thresh = F,
#                           test.stat = "z", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                           use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                           err.thresh = .05, centre.specific.thresh = F)
#   out.perm$restab
#   out.perm.lfsr <- err.rate.control(resl = resl.comp[!names(resl.comp) %in% c("uv", "uv.ss")],
#                                err.rate.meth = "perm", sep.imp.thresh = F,
#                                test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                                use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                                err.thresh = .05, centre.specific.thresh = F)
#   out.perm.lfsr$restab
#   out.lfsr <- err.rate.control(resl = resl.comp[!names(resl.comp) %in% c("uv", "uv.ss")],#[names(resl.comp) %in% c("eb", "ed", "mash")],
#                                err.rate.meth = "lfsr", sep.imp.thresh = F,
#                                test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                                use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                                err.thresh = .05, centre.specific.thresh = F)
#   out.lfsr$restab
#   colkeep <- c("meth", "err.rate.meth", "test.stat", "mvhitimp", "mvhitnonimp", "line.fdr.est", 
#                   "fdr.est", "fdr.est.imp", "fdr.est.nonimp")
#   restaball <- rbind(out.perm$restab[, colkeep], out.perm.lfsr$restab[, colkeep], out.lfsr$restab[, colkeep])
#   restaball.out <- restaball
#   for(j in 1:ncol(restaball)){
#     if(is.numeric(restaball[, j]))
#       restaball.out[, j] <- round(restaball[, j] * 100, 1)
#   }
#   restaball.out
#   resl.out$uv <- out.perm$resl$uv
#   resl.out$eb <- out.perm$resl$impc_eb_1
#   if(collect.factors){
#     for(fac.meth in fac.methv)
#       resl.out[[fac.meth]] <- out.perm$resl[[fac.meth]]
#   }
#   resimp <- out.perm$resimp
#   resimp$line.type <- ifelse(out.perm$resimp$line.type == "trueMutTes", "trueMut", "negCon")
# }
# 
# 
# 
# save(resl.out, resimp, Ksigl, Sigll, resl.comp, compl, file = file.glob.res, version = 2)
# 
# 
# 
# sapply(Ksigl, mean)
# 

# load(uv.resmat.file.post.qc)
# resmat.qc$uv.ss.t <-  resmat.qc$uv.mn.inter.genosex.sc / resmat.qc$uv.sd.inter.genosex.sc
# resmat.qc$uv.t <-  resmat.qc$uv.mn.sc / resmat.qc$uv.sd.sc
# resmat.qc0 <- resmat.qc[resmat.qc$line.type == "negCon", ]
# boxplot(uv.ss.t ~ mut.n.female, data = resmat.qc0)
# boxplot(uv.t ~ mut.n.female, data = resmat.qc0)
# length(which(resmat.qc0$uv.ss.t < -5))
# resmat.qc0[which(resmat.qc0$uv.ss.t < -5)[1:10], ]
# sort(table(resmat.qc0[which(resmat.qc0$uv.ss.t < -5), "ph"]))
# sort(unique(resmat.qc[which(resmat.qc$uv.ss.t < -10), "geno"]))
# table(resmat.qc0$uv.sd.inter.genosex.sc < resmat.qc0$uv.sd.inter.geno.sc)
# 
# seed <- 1
# resl.seed <- list()
# resl.seed$uv <- list(mn = Yhat[all.lines, ph.use], sd = smat[all.lines, ph.use])
# for(namc in names(compl)){
#   resl.seed[[namc]] <- list(mn = compl[[namc]]$mnarr[all.lines, ph.use, seed], sd = compl[[namc]]$sdarr[all.lines, ph.use, seed], 
#                                lfsr = compl[[namc]]$lfsrarr[all.lines, ph.use, seed])
# }
# 
# out.perm.seed <- err.rate.control(resl = resl.seed[!names(resl.seed) %in% c("uv.ss")], 
#                              err.rate.meth = "perm", sep.imp.thresh = F,
#                         test.stat = "z", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                         use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                         err.thresh = .05, centre.specific.thresh = F)
# out.perm.seed$restab
# out.perm.lfsr.seed <- err.rate.control(resl = resl.seed[!names(resl.seed) %in% c("uv", "uv.ss")],
#                              err.rate.meth = "perm", sep.imp.thresh = F,
#                              test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                              use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                              err.thresh = .05, centre.specific.thresh = F)
# out.perm.lfsr.seed$restab
# out.lfsr.seed <- err.rate.control(resl = resl.seed[!names(resl.seed) %in% c("uv", "uv.ss")],#[names(resl.comp) %in% c("eb", "ed", "mash")],
#                              err.rate.meth = "lfsr", sep.imp.thresh = F,
#                              test.stat = "lfsr", linemap = linemap, phmap = phmap, Yhat = Yhat,
#                              use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                              err.thresh = .05, centre.specific.thresh = F)
# out.lfsr.seed$restab
# 



# 
# 
# resimpc <- out.perm.lfsr$resimp
# resimpc0 <- resimpc[resimpc$line.type == "negConTes", ]
# resimpc1 <- resimpc[resimpc$line.type == "trueMutTes", ]
# str(resimpc0)
# table(resimpc$line.type)
# 
# str(resimpc0)
# 
# length(unique(resimpc0$geno))
# table(resimpc0$geno[resimpc0$impc_mash_1.perm.signsig != 0])
# hist(table(resimpc1$geno[resimpc1$impc_mash_1.perm.signsig != 0]))
# mean(is.na(resimpc0$impc_mash_1.perm.signsig))
# 
# str(resl.comp,m=1)
# 
# 
# compl$impc_mash_1$mn.comb["10061_0_0_negcon_1", ]
# compl$impc_mash_1$sd.comb["10061_0_0_negcon_1", ]
# compl$impc_mash_1$lfsr.comb["10061_0_0_negcon_1", ]
# str(resl.comp$impc_mash_1$lfsr["10061_0_0_negcon_1", ])
# str(compl$impc_mash_1$lfsrarr["10061_0_0_negcon_1", , 1])
# str(compl$impc_mash_1$mnarr["10061_0_0_negcon_1", , 1])
# str(compl$impc_mash_1$sdarr["10061_0_0_negcon_1", , 1])
# str(compl$impc_mash_1, m = 1)
# hist(resl.comp$mash$mn / resl.comp$mash$sd)
# hist(resl.comp$mash.no.imp$mn / resl.comp$mash.no.imp$sd)
# par(mfrow = c(2, 1))
# hist(-log10(resl.comp$mash$lfsr))
# hist(-log10(resl.comp$mash.no.imp$lfsr))

# mean(is.na(resl.comp$mash.no.imp$lfsr))
# hist(resl.comp$mash$lfsr[is.na(resl.comp$mash.no.imp$lfsr)])
# hist(resl.comp$mash$lfsr[!is.na(resl.comp$mash.no.imp$lfsr)])
# mn / resl.comp$mash$sd)
# hist(resl.comp$mash.no.imp$mn / resl.comp$mash.no.imp$sd)
      # file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name, sapply(var.in.name, get), sep = "_"), collapse = "_"))
      # emout.file.namc <- paste0(file.base, "_emout.RData")
      # emout.file.namc.with.ss <- paste0(file.base, "_sexspecific_", sexspecific, "_emout.RData")
      # res.store.namc.with.ss <- paste0(file.base, "_sexspecific_", sexspecific, "_res.RData")
      # loocv.res.store.namc <- paste0(file.base, "_loocv_res.RData")
      # fac.res.store.namc <- paste0(file.base, "_facres.RData")
      # load(file = emout.file.namc)
      # load(res.store.namc.with.ss)
      # load(file = res.store.namc)
      # file.info(res.store.namc)
# unique(phmap$procnam)
# proclook <- "Acoustic Startle and Pre-pulse Inhibition (PPI)"
# phvlook <- phmap[phmap$procnam == proclook, "ph"]
# plot(facs[, which.max(colMeans(abs(facs[phvlook, ])))])


# resl = resl.comp; err.rate.meth = "perm"; sep.imp.thresh = F;
# test.stat = "z"; linemap = linemap; phmap = phmap; Yhat = Yhat;
# use.upper.fp.est = F; control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1]
# err.thresh = .05; centre.specific.thresh = F

# str(resl.comp)
# # if(collect.factors)
# #   resl.comp$eb.fac <- list(mn = resl.comp$eb$fac.mn, sd = resl.comp$eb$fac.sd)
# resl.comp <- resl.comp
# resl.err.rate <- resl.comp[names(resl.comp) %in% c("uv", "uv.ss", "eb", "mash", "varimax", "promax")]

# load(file = "C:/Temp/genobad.RData")
# for(j in 1:length(resl.err.rate)){
#   resl.err.rate[[j]]$mn <- resl.err.rate[[j]]$mn[!rownames(resl.err.rate[[j]]$mn) %in% c(geno.bad.ko, geno.bad.wt), ]
#   resl.err.rate[[j]]$sd <- resl.err.rate[[j]]$sd[!rownames(resl.err.rate[[j]]$sd) %in% c(geno.bad.ko, geno.bad.wt), ]
# }


# linemap <- linemap[!linemap$geno %in% c(geno.bad.ko, geno.bad.wt), ]
# str(resl.err.rate)

# resl = resl.err.rate; err.rate.meth = "perm"; sep.imp.thresh = F;
# test.stat = "z"; linemap = linemap; phmap = phmap; Yhat = Yhat;
# use.upper.fp.est = F; control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1];
# err.thresh = .05; centre.specific.thresh = F

# str(resl.comp)
# resl.comp$mash.no.imp <- resl.comp$mash
# resl.comp$mash.no.imp$mn[is.na(resl.comp$uv$mn)] <- NA
# resl.comp$mash.no.imp$mn[is.na(resl.comp$uv$mn)] <- NA
# resl.comp$mash.no.imp$lfsr[is.na(resl.comp$uv$mn)] <- NA
# 
# 
# mean(is.na(resl.comp$mash.no.imp$sd))
# mean(is.na(resl.comp$mash$sd))
# tab.all.seed <- data.frame()
# for(i in 1:n.seed.run.subsams)
#   tab.all.seed <- rbind(tab.all.seed, tab[scen, ])
# tab.all.seed$seed <- 1:n.seed.run.subsams
# compl <- Sigll <- Ksigl <- list()
# for(mv.meth in mv.methv){
#   compl[[mv.meth]] <- shared.formatl
#   Sigll[[mv.meth]] <- list()
#   Ksigl[[mv.meth]] <- c()
# }

# vc.type <- c("vari", "pro")[1]
# if(vc.type == "vari")
#   facs <- varimax(eigen(Sigcor)$vectors[, 1:nfac])$loadings
# if(vc.type == "pro")
#   facs <- promax(eigen(Sigcor)$vectors[, facnam])$loadings
# 
# fcheckfields <- c("emout.ok", "emout.st", "emout.en", "res.ok", "res.st", "res.en", "loocv.ok", "loocv.st", "loocv.en", "fac.ok", "fac.st", "fac.en")
# runtab[, fcheckfields] <- NA
# #############################################
# # Check files are up to date
# for(scen in 1:nrow(runtab)){
#   for(j in 1:ncol(runtab))
#     assign(colnames(runtab)[j], runtab[scen, j])
#   Data <- da
#   file.base <- c()
#   if(grepl("em.fit", me)){
#     var.in.name.use <- var.in.name
#     if(rand)
#       var.in.name.use <- c(var.in.name, "rand")
#     meth <- "eb"
#   }
#   if(me == "mash"){
#     var.in.name.use <- var.in.name.mash
#     meth <- "mash"
#   }
#   if(me == "ed"){
#     var.in.name.use <- var.in.name.ed
#     meth <- "XD"
#   }
#   
#   for(seed in 1:ntimes)
#     file.base <- c(file.base, paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.use, sapply(var.in.name.use, get), sep = "_"), collapse = "_")))
#   # emout.file.namc <- paste0(file.base, "_emout.RData")
#   # res.store.namc <- paste0(file.base, "_res.RData")
#   loocv.res.store.namc <- paste0(file.base, "_loocv_res.RData")
#   fac.res.store.namc <- paste0(file.base, "_facres.RData")
#   # mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
#   # bovy.resl.file.namc <- paste0(file.base.ed, "_bovy_resl.RData")
#   res.store.namc <- paste0(file.base, switch(meth, eb = "_res.RData", mash = "_mash_resl.RData", XD = "_bovy_resl.RData"))
#   emout.file.namc <- paste0(file.base, switch(meth, eb = "_emout.RData", mash = NA, XD = NA))
#   runtab[scen, fcheckfields] <- list(all(file.exists(emout.file.namc)),
#                                                                   format(min(file.info(emout.file.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   format(max(file.info(emout.file.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   all(file.exists(res.store.namc)),
#                                                                   format(min(file.info(res.store.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   format(max(file.info(res.store.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   all(file.exists(loocv.res.store.namc)),
#                                                                   format(min(file.info(loocv.res.store.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   format(max(file.info(loocv.res.store.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   all(file.exists(fac.res.store.namc)),
#                                                                   format(min(file.info(fac.res.store.namc)[, "ctime"]), format = "%d-%b"),
#                                                                   format(max(file.info(fac.res.store.namc)[, "ctime"]), format = "%d-%b"))
#   
# }
# runtab
# 
# 

# # si <- 1
# # bo <- 1
# # EDmeth <- "justED"
# # EDtol <- 1e-4
# # mv.methv <- c("mash", "eb", "ed", "eb.ss", "eb.fac")[2]
# # fac.methv <- c("varimax", "promax")
# # collect.factors <- T
# # n.seed.run.subsams <- 10
# # if(full.analysis){
# #   for(i in 1:length(default.parameters[[Data]]))
# #     assign(names(default.parameters[[Data]][i]), default.parameters[[Data]][[i]])
# # }
# P <- 148#default.parameters[[Data]]$P
# N <- 500#default.parameters[[Data]]$N
# # K <- floor(P / 2)
# compdat.default <- data.frame(Data = Data, meth = c("ed", "ed", "ed", "mash", "eb", "eb", "eb"), nSig = c(1:3, 1, 1:3), 
#                       EMtol = c(1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4), 
#                       EDtol = c(1e-4, 1e-4, 1e-4, 1e-5, 1e-4, 1e-4, 1e-4))
# compdat.default$N <- default.parameters[[Data]]$N
# compdat.default$P <- default.parameters[[Data]]$P
# other.par.combs <- 
# compdat.extra <- compdat.default
# compdat[compdat$meth == "eb", "EDtol"] <- NA
# compdat[compdat$meth != "eb", "EMtol"] <- NA
# compdat <- compdat[order(compdat$meth != "eb"), ]
# if(N == 500)
#   compdat <- compdat[compdat$meth == "eb" & compdat$nSig == 1, ]
# compdat <- compdat[!(compdat$meth == "eb" & compdat$nSig > 1), ]
# compdat <- compdat[compdat$meth != "ed", ]
# compdat <- compdat[!(compdat$meth == "eb" & compdat$nSig > 1), ]























# for(scen in 1:nrow(runtab)){
#   print(runtab[scen, c("N", "P", "me", "meth", "da")])
#   for(j in 1:ncol(runtab))
#     assign(colnames(runtab)[j], runtab[scen, j])
#   for(j in 1:length(dimnaml[[Data]]))
#     assign(names(dimnaml[[Data]])[j], dimnaml[[Data]][[j]]) #sam.names, meas.names, n.all, P.all
#   array1 <- array(NA, dim = c(N.all, P.all, n.seed.run.subsams), dimnames = list(sam.names, meas.names, 1:n.seed.run.subsams))
#   array2 <- array(NA, dim = c(N.all, n.seed.run.subsams), dimnames = list(sam.names, 1:n.seed.run.subsams))
#   array.fac <- array(NA, dim = c(N.all, nfac, n.seed.run.subsams), dimnames = list(sam.names, facnam, 1:n.seed.run.subsams))
#   crossval.llv <- rep(NA, n.seed.run.subsams)
#   shared.formatl <- list(mnarr = array1, sdarr = array1, loocv.mnarr = array1, loocv.sdarr = array1, 
#                          varimax = list(fac.mnarr = array.fac, fac.sdarr = array.fac), 
#                          promax = list(fac.mnarr = array.fac, fac.sdarr = array.fac), 
#                          lfsrarr = array1, 
#                          llmat = array2, llmat.raw = array2, llmat.zero = array2)
#   compl <- Sigll <- Ksigl <- pimatl <- list()
#   namc <- paste(Data, meth, nSig, sep = "_")
#   compl[[namc]] <- shared.formatl
#   pimatl[[namc]] <- Ksigl[[namc]] <- Sigll[[namc]] <- list()
#   # for(namc in mv.methv){
#   for(seed in 1:n.seed.run.subsams){
#     if(meth %in% c("eb", "eb.ss")){
#       sexspecific <- switch(meth, eb = F, eb.ss = T)
#       var.in.name.tryl <- list(c(var.in.name.old, "sexspecific"), var.in.name.old, var.in.name)[3]
#       Ktry <- c(P, P / 2)
#       emloaded <- resloaded <- loocv.loaded <- fac.loaded <- F
#       for(K in Ktry){
#         for(var.in.namec in var.in.name.tryl){
#           file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.namec, sapply(var.in.namec, get), sep = "_"), collapse = "_"))
#           emout.file.namc <- paste0(file.base, "_emout.RData")
#           res.store.namc <- paste0(file.base, "_res.RData")
#           loocv.res.store.namc <- paste0(file.base, "_loocv_res.RData")
#           fac.res.store.namc <- paste0(file.base, "_facres.RData")
#           if(file.exists(emout.file.namc)){
#             emloaded <- T
#             load(file = emout.file.namc)
#           }
#           if(file.exists(res.store.namc)){
#             resloaded <- T
#             load(file = res.store.namc)
#           }
#           if(file.exists(fac.res.store.namc)){
#             load(file = fac.res.store.namc)
#             fac.loaded <- T
#           }
#         }
#       }
#       if(emloaded == F | resloaded == F )
#         warning(paste0(namc, " seed ", seed, " not loaded"))
#       res.store <- resl.store$raw
#       Ksigl[[namc]][[seed]] <- emout.mix$Ksig
#     }
#     if(meth == "mash"){
#       var.in.name.mash <- setdiff(var.in.name.mash, "sexspecific")
#       file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
#       # file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, tab.all.seed[seed, var.in.name.mash], sep = "_"), collapse = "_"))
#       mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
#       load(file = mash.resl.file.namc)
#       if("zero" %in% names(resl.store)){
#         res.store <- resl.store$zero
#       } else {
#         res.store <- resl.store[[1]]
#         warning("Need to re-run mash on this data set")
#       }
#     }
#     if(meth == "ed"){
#       var.in.name.ed <- setdiff(var.in.name.ed, "sexspecific")
#       file.base.ed <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.ed, sapply(var.in.name.ed, get), sep = "_"), collapse = "_"))
#       bovy.resl.file.namc <- paste0(file.base.ed, "_bovy_resl.RData")
#       load(file = bovy.resl.file.namc)
#       res.store <- resl.store$zero
#     }
#     Sigll[[namc]][[seed]] <- lapply(res.store$Sigl, function(M) M[match(meas.names, rownames(M)), match(meas.names, colnames(M))])
#     pimatl[[namc]][[seed]] <- res.store$pimat
#     compl[[namc]]$mnarr[match(rownames(res.store$mn), sam.names), 
#                         match(colnames(res.store$mn), meas.names), seed] <- res.store$mn
#     compl[[namc]]$sdarr[match(rownames(res.store$sd), sam.names), 
#                         match(colnames(res.store$sd), meas.names), seed] <- res.store$sd
#     compl[[namc]]$lfsrarr[match(rownames(res.store$lfsr), sam.names), 
#                           match(colnames(res.store$lfsr), meas.names), seed] <- res.store$lfsr
#     compl[[namc]]$llmat[match(names(res.store$loglikv), sam.names), seed] <- res.store$loglikv
#     compl[[namc]]$llmat.raw[match(names(resl.store$raw$loglikv), sam.names), seed] <- resl.store$raw$loglikv
#     compl[[namc]]$llmat.zero[match(names(resl.store$zero$loglikv), sam.names), seed] <- resl.store$zero$loglikv
#     if(meth == "eb" & nSig == 1){
#       if(loocv.loaded){
#         # load(file = loocv.res.store.namc)
#         compl[[namc]]$loocv.mnarr[match(rownames(loocv.store$mn), sam.names), 
#                                   match(colnames(loocv.store$mn), meas.names), seed] <- loocv.store$mn
#         compl[[namc]]$loocv.sdarr[match(rownames(loocv.store$sd), sam.names), 
#                                   match(colnames(loocv.store$sd), meas.names), seed] <- loocv.store$sd
#       }
#       if(fac.loaded){
#         # load(fac.res.store.namc)
#         # table(sign(fac.res.store$mn[rownames(fac.res.store$mn) %in% true.use, "f.1"]))
#         # fac.res.store$loadings[lean.ph, "f.1"]
#         for(fac.meth in fac.methv){
#           # compl[[namc]][[fac.meth]] <- list()
#           compl[[namc]][[fac.meth]]$fac.mnarr[match(rownames(fac.res.store[[fac.meth]]$mn), sam.names), 
#                                               match(colnames(fac.res.store[[fac.meth]]$mn), facnam), seed] <- fac.res.store[[fac.meth]]$mn
#           compl[[namc]][[fac.meth]]$fac.sdarr[match(rownames(fac.res.store[[fac.meth]]$sd), sam.names), 
#                                               match(colnames(fac.res.store[[fac.meth]]$sd), facnam), seed] <- fac.res.store[[fac.meth]]$sd
#         }
#       }
#     }
#   }
# }

# str(compl, m = 2)
# ######################################################
# # Bayesian model combination
# resl.comp <- list()
# # resl.comp$uv.ss <- list(mn = Yhat.genosex[all.lines, ph.use], sd = smat.genosex[all.lines, ph.use])
# resl.comp$uv <- list(mn = Yhat[all.lines, ph.use], sd = smat[all.lines, ph.use])
# # for(mv.meth in mv.methv){
# for(scen in 1:nrow(compdat)){
#   for(j in 1:ncol(compdat))
#     assign(colnames(compdat)[j], compdat[scen, j])
#   namc <- paste(Data, meth, nSig, sep = "_")
#   suppressWarnings(lmat.norm <- exp(compl[[namc]]$llmat - apply(compl[[namc]]$llmat, 1, function(v) max(v, na.rm = T))))
#   pmix <- lmat.norm / rowSums(lmat.norm, na.rm = T)
#   if(meth %in% c("eb", "eb.ss")){
#     Sigmix.p <- colSums(pmix, na.rm = T) / sum(pmix, na.rm = T)
#     Sig.comb <- 0
#     for(seed in 1:n.seed.run.subsams){
#       if(!is.null(Sigll[[namc]][[seed]][[1]])){
#         Sig.comb <- Sig.comb + Sigll[[namc]][[seed]][[1]][ph.use, ph.use] * c(Sigmix.p[seed])
#       } else {
#         warning(paste0("Analysis not complete for all seed ", seed))
#       }
#     }
#   }
#   wt.mnarr <- apply(sweep(compl[[namc]]$mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#   wt.varr <- apply(sweep(compl[[namc]]$sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#   wt.mnsqarr <- apply(sweep(compl[[namc]]$mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#   wt.lfsrarr <- apply(sweep(compl[[namc]]$lfsrarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#   compl[[namc]]$mn.comb <- wt.mnarr
#   compl[[namc]]$sd.comb <- sqrt(wt.varr + wt.mnsqarr - wt.mnarr^2)
#   compl[[namc]]$lfsr.comb <- wt.lfsrarr
#   if(meth == "eb" & nSig == 1){
#     loocv.wt.mnarr <- apply(sweep(compl[[namc]]$loocv.mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#     loocv.wt.varr <- apply(sweep(compl[[namc]]$loocv.sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#     loocv.wt.mnsqarr <- apply(sweep(compl[[namc]]$loocv.mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#     compl[[namc]]$loocv.mn.comb <- loocv.wt.mnarr
#     compl[[namc]]$loocv.sd.comb <- sqrt(loocv.wt.varr + loocv.wt.mnsqarr - loocv.wt.mnarr^2)
#     if(collect.factors){
#       for(fac.meth in fac.methv){
#         fac.wt.mnarr <- apply(sweep(compl[[namc]][[fac.meth]]$fac.mnarr, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#         fac.wt.varr <- apply(sweep(compl[[namc]][[fac.meth]]$fac.sdarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#         fac.wt.mnsqarr <- apply(sweep(compl[[namc]][[fac.meth]]$fac.mnarr^2, c(1, 3), pmix, '*'), 1:2, function(v) sum(v, na.rm = T))
#         compl[[namc]][[fac.meth]] <- list(fac.mn.comb = fac.wt.mnarr, 
#                                           fac.sd.comb = sqrt(fac.wt.varr + fac.wt.mnsqarr - fac.wt.mnarr^2))
#       }
#       # compl[[namc]][[fac.meth]]$fac.mn.comb <- fac.wt.mnarr
#       # compl[[namc]][[fac.meth]]$fac.sd.comb <- sqrt(fac.wt.varr + fac.wt.mnsqarr - fac.wt.mnarr^2)
#     }
#   }
#   resl.comp[[namc]] <- list(mn = compl[[namc]]$mn.comb[all.lines, ph.use], sd = compl[[namc]]$sd.comb[all.lines, ph.use], 
#                             loocv.mn = compl[[namc]]$loocv.mn.comb[all.lines, ph.use], loocv.sd = compl[[namc]]$loocv.sd.comb[all.lines, ph.use],
#                             lfsr = compl[[namc]]$lfsr.comb[all.lines, ph.use])
#   if(collect.factors){
#     for(fac.meth in fac.methv){
#       resl.comp[[fac.meth]] <- list(mn = compl$eb[[fac.meth]]$fac.mn.comb[all.lines, facnam], 
#                                     sd = compl$eb[[fac.meth]]$fac.sd.comb[all.lines, facnam])
#     }
#   }
#   if(meth %in% c("eb", "eb.ss"))
#     resl.comp[[namc]]$Sig.comb <- Sig.comb
# }

# 
# save(resl.comp, file = file.resl.comp, version = 2)
# save(resl.comp, file = file.resl.comp, version = 2)
# 
# file.compl <- paste0(global.res.dir, "/global_compl.RData")
# file.resl.comp <- paste0(global.res.dir, "/global_reslcomp.RData")
# file.resl.comp.fac <- paste0(global.res.dir, "/global_reslcompfac.RData")
# file.objl <- paste0(global.res.dir, "/global_objl.RData")
# if(!is.null(fac.loaded)){
#   for(fac.meth in fac.methv){
#     compl[[namc]][[fac.meth]]$fac.mnarr[match(rownames(fac.res.store[[fac.meth]]$mn), sam.names), 
#                                         match(colnames(fac.res.store[[fac.meth]]$mn), facnam), seed] <- fac.res.store[[fac.meth]]$mn
#     compl[[namc]][[fac.meth]]$fac.sdarr[match(rownames(fac.res.store[[fac.meth]]$sd), sam.names), 
#                                         match(colnames(fac.res.store[[fac.meth]]$sd), facnam), seed] <- fac.res.store[[fac.meth]]$sd
#   }
# }
