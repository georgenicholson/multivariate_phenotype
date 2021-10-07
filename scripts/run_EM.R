rm(list = ls())
xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))


#######################################
#Gather arguments from command line
arguments <- commandArgs()
if("--args" %in% arguments){
  argnam.in <- data.frame(nam = c("seed", "N", "P", "run.bovy", "include.singletons", 
                                  "Data", "methods", "max.num.sams.for.testing", "full.analysis", "sexspecific", "nSig", 
                                  "EDmeth", "EDtol", "EMtol", "EMmash",
                                  "EMfm", "EMbic", "EMwish", 
                                  "EMK", "EMKup", "rand", "loocv"),
                          coersion.fn = c("as.numeric", "as.numeric", "as.numeric", "as.numeric", "as.numeric", 
                                          "as.character", "as.character", "as.numeric", "as.logical", "as.logical", "as.numeric", 
                                          "as.character", "as.numeric", "as.numeric", "as.logical", 
                                          "as.character", "as.numeric", "as.logical", 
                                          "as.numeric", "as.logical", "as.logical", "as.logical"), 
                          stringsAsFactors = F)
  for(i in 1:nrow(argnam.in))
    assign(argnam.in$nam[i], eval(call(argnam.in$coersion.fn[i], arguments[grep("--args", arguments) + i])))
  print(arguments)
  print(paste0("nSig = ", nSig))
} else {
  seed <- 10   # seed <- as.numeric(format(Sys.time(), "%s")) %% 49 + 1
  N <- 200
  P <- 20
  run.bovy <- 1
  include.singletons <- 0
  Data <- c("eqtl", "impc")[2]
  methods <- c("em.fit", "em.test", "em.loocv", "em.fac", "mash", "em.all.test", "ed")[c(1, 2, 5, 7)[1]]
  max.num.sams.for.testing <- Inf
  full.analysis <- F
  sexspecific <- F
  EDmeth <- c("mash", "justED")[1]
  EDtol <- 1e-5 # Default
  EMtol <- 1e-4 # Default
  EMfm <- c("fa", "pca")[1]
  EMbic <- 0
  EMmash <- F
  EMwish <- F
  EMKup <- F
  EMK <- 10
  nSig <- 1
  rand <- F
  loocv <- F
}
if(full.analysis)
  P <- switch(Data, impc = 148, eqtl = 44)
EDnSig <- switch(EDmeth, mash = 3, justED = nSig)



create_table_of_analyses()


source(paste0(R.file.dir, "/impc_mv_paper_code/EM_algo_lean.R"))
source(paste0(R.file.dir, "/impc_mv_paper_code/EM_fns_lean.R"))
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
dat.typev <- switch(Data, impc = c("raw", "zero"), eqtl = "raw")

###################################
# Load data
dir.dataset <- paste0(sub.data.sets.dir, "/", Data)
npdir <- paste0(dir.dataset, "/N_", N, "_P_", P)
file.in <- paste0(npdir, "/", Data, "_N_", N, "_P_", P, "_seed_", seed, ".RData")
suppressWarnings(load(file.in))
if(Data == "eqtl"){
  mash.in <- readRDS(file = paste0(R.file.dir, "/urbut_code/stephenslab-gtexresults-961b969/data/MatrixEQTLSumStats.Portable.ld2.Z.rds"))
  Y <- Yhat <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
  S <- smat <- Yhat / rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
}
if(Data == "impc"){
  load(file = uv.results.Y.S)
  if(sexspecific){
    Y <- Y.genosex
    S <- S.genosex
    Yhat <- Yhat.genosex
    smat <- smat.genosex
  }
  ref <- read.csv(paste0(base.dir, "/data_in/reference_line_genotypeIds.csv"))
}

Y.eml <- S.eml <- list()
Y.eml$raw <- Yhat
Y.eml$zero <- Y
S.eml$raw <- smat
S.eml$zero <- S

##############################################
#Initialise R
ident.eps <- .05
cor.type <- c("weighted", "unweighted")[1]
if(Data == "impc"){
  if(cor.type == "weighted")
    R.init <- cor((Y / S)[sams.for.cor.est, ph.use], use = "p", meth = "p")
  if(cor.type == "unweighted")
    R.init <- cor(Y[sams.for.cor.est, ph.use], use = "p", meth = "p")
  R.init[is.na(R.init)] <- 0
  diag(R.init) <- 1
  if(qr(R.init)$rank < P)
    R.init <- (1 - ident.eps) * R.init + ident.eps * diag(rep(1, P))
}
if(Data == "eqtl"){
  snps.est.cor <- snpmap.sub$snp[snpmap.sub$task == "random.train.est.cor"]
  data.in <- list(Bhat = Y[sams.for.cor.est, ph.use], Shat = S[sams.for.cor.est, ph.use])
  R.init <- estimate_null_correlation_simple(data = data.in, z_thresh = 2)
  # Sig.init <- cov(Y[sams.for.model.fitting, ph.use], use = "p", meth = "p")
}
dimnames(R.init) <- list(ph.use, ph.use)

##############################################
#Initialise list of cov. mat's Sigl.em.init
if(!rand){
  if(nSig > 1){
    require(mclust)
    mcl.out <- Mclust(data = Y[sams.for.model.fitting, ph.use], G = nSig)
    cluster.membership <- apply(mcl.out$z, 1, which.max)
  } else {
    cluster.membership <- rep(1, length(sams.for.model.fitting))
  }
  names(cluster.membership) <- sams.for.model.fitting
  Sigl.em.init <- list()
  for(j in 1:nSig){
    sams.in <- names(cluster.membership)[which(cluster.membership == j)]
    Sig.init <- cov(Y[sams.in, ph.use], use = "p", meth = "p")
    Sig.init[is.na(Sig.init)] <- 0
    diag(Sig.init)[is.na(diag(Sig.init))] <- 1
    if(qr(Sig.init)$rank < P)
      Sig.init <- (1 - ident.eps) * Sig.init + ident.eps * diag(rep(1, P))
    dimnames(Sig.init) <- list(ph.use, ph.use)
    Sigl.em.init[[j]] <- Sig.init
  }
} else {
  Sigl.em.init <- list()
  for(j in 1:nSig){
    Sig.temp <- solve(rWishart(1, P, diag(rep(1, P)))[, , 1])
    Sigl.em.init[[j]] <- t(Sig.temp / sqrt(diag(Sig.temp))) / sqrt(diag(Sig.temp))
  }
}
update.Sig <- rep(T, nSig)
names(update.Sig) <- paste0("Sig", 1:nSig)

if(EMmash){
  add.singletons <- T
  add.ident <- T
  singl <- list()
  if(add.singletons){
    for(j in 1:P){
      namc <- paste0(ph.use[j], " singleton")
      diagc <- rep(0, P)
      diagc[j] <- 1
      singl.add <- list(diag(diagc))
      update.Sig.add <- F
      names(update.Sig.add) <- names(singl.add) <- namc
      singl <- c(singl, singl.add)
      update.Sig <- c(update.Sig, update.Sig.add)
    }
    Sigl.em.init <- c(Sigl.em.init, singl)
  }
  if(add.ident){
    Sigl.em.init <- c(Sigl.em.init, list(identity = diag(rep(1, P))))
    update.Sig <- c(update.Sig, c(identity = F))
  }
  Sigl.em.init <- lapply(Sigl.em.init, function(M){ dimnames(M) <- list(ph.use, ph.use); M})
}


K.init <- ifelse(update.Sig, EMK, NA)
update.K <- ifelse(update.Sig, EMKup, F)


##############################
# Run EM mixture
omegaseq <- exp((-10):10 * log(2))
# update.Sig <- T
# update.K <- T
# var.in.name <- c(var.in.name.old, "sexspecific")
var.in.name.use <- var.in.name
if(rand)
  var.in.name.use <- c(var.in.name.use, "rand")
if(EMmash)
  var.in.name.use <- c(var.in.name.use, "EMmash")
file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.use, sapply(var.in.name.use, get), sep = "_"), collapse = "_"))
emout.file.namc <- paste0(file.base, "_emout.RData")
res.store.namc <- paste0(file.base, "_res.RData")
resl.store.namc <- paste0(file.base, "_resl.RData")
fac.res.store.namc <- paste0(file.base, "_facres.RData")
loocv.res.store.namc <- paste0(file.base, "_loocv_res.RData")
print(emout.file.namc)
print(res.store.namc)
# if(grepl("em.fit",  any(c("em.fit", "em.all") %in% methods)){
if(any(grepl("em.fit", methods))){
  rerun.em <- T
  if(rerun.em | !file.exists(emout.file.namc)){
    print("Running EM algorithm")
    # Y.em = Yhat[sams.for.model.fitting, ph.use]; S.em = smat[sams.for.model.fitting, ph.use];
    # Sigl.em.init = Sigl.em.init; R.em.init = R.init; omegaseq = omegaseq; K.init = K.init; seed = seed;
    # update.Sig = update.Sig; update.K = update.K; pi.init = NULL; auto.tol.const = EMtol;
    # fac.model = EMfm; bic.pen.mult = EMbic; wish.pri = EMwish
    emout.mix <- EM_algo_mixture_multi_Sig(Y.em = Yhat[sams.for.model.fitting, ph.use], S.em = smat[sams.for.model.fitting, ph.use], 
                            Sigl.em.init = Sigl.em.init, R.em.init = R.init, omegaseq = omegaseq, K.init = K.init, seed = seed, 
                            update.Sig = update.Sig, update.K = update.K, pi.init = NULL, auto.tol.const = EMtol, 
                            fac.model = EMfm, bic.pen.mult = EMbic, wish.pri = EMwish)
    save(emout.mix, file = emout.file.namc)
  } else {
    print("Loading previously run EM output")
    load(file = emout.file.namc)
    emout.mix$Sigl <- lapply(emout.mix$Sigl, function(M){ dimnames(M) <- list(ph.use, ph.use); M})
    dimnames(emout.mix$Sig.mn) <- list(ph.use, ph.use)
  }
  print("Calculating posterior means")
  resl.store <- list()
  dat.type <- "raw"#c("zero", "raw")[2]#for(dat.type in dat.typev){
  # for(dat.type in dat.typev){
  # Y.em = Y.eml[[dat.type]][sams.for.testing, ph.use]; S.em = S.eml[[dat.type]][sams.for.testing, ph.use]
  # Sigl = emout.mix$Sigl; R = emout.mix$R; omegaseq = emout.mix$omegaseq;
  # pimat = emout.mix$pi; meth = "post.mn"; fac.model = EMfm; bic.pen.mult = EMbic
  # llmat = NULL; recalc.llmat = T
    out.post.mix <- em.update.function(Y.em = Y.eml[[dat.type]][sams.for.testing, ph.use], S.em = S.eml[[dat.type]][sams.for.testing, ph.use],
                                       Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
                                       pimat = emout.mix$pi, meth = "post.mn", fac.model = EMfm, bic.pen.mult = EMbic, wish.pri = EMwish)
    resl.store[[dat.type]] <- list(mn = out.post.mix$mnmat[sams.for.testing, ph.use], sd = out.post.mix$sdmat[sams.for.testing, ph.use],
                                   loglikv = out.post.mix$loglikv[sams.for.testing],
                                   lfsr = out.post.mix$lfsrmat[sams.for.testing, ph.use], 
                                   Sigl = lapply(emout.mix$Sigl, function(M) M[ph.use, ph.use]), 
                                   Sig.mn = emout.mix$Sigmn[ph.use, ph.use],
                                   Ksig = emout.mix$Ksig, R = emout.mix$R[ph.use, ph.use], pimat = emout.mix$pi, omegaseq = omegaseq)
  #}
  save(resl.store, file = res.store.namc)
  calc.loocv <- F
  if(nSig == 1 & "em.fit" %in% methods & calc.loocv){
    print("Calculating leave-one-procedure-out predictions")
    load(file = emout.file.namc)
    procun <- unique(phmap$procnam)
    matout <- matrix(NA, length(sams.for.testing), length(ph.use), dimnames = list(sams.for.testing, ph.use))
    loocv.store <- list(mn = matout, sd = matout)
    for(procc in procun){
      ph.leave.out <- phmap$ph[phmap$procnam == procc]
      Yhat.left.out <- Yhat
      Yhat.left.out[, ph.leave.out] <- NA
      smat.left.out <- smat
      smat.left.out[, ph.leave.out] <- NA
      post.mn.loocv <- em.update.function(Y.em = Yhat.left.out[sams.for.testing, ph.use], S.em = smat.left.out[sams.for.testing, ph.use],
                                          Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
                                          pimat = emout.mix$pi, meth = "post.mn", fac.model = EMfm, bic.pen.mult = EMbic, wish.pri = EMwish)
      loocv.store$mn[sams.for.testing, ph.leave.out] <- post.mn.loocv$mnmat[sams.for.testing, ph.leave.out]
      loocv.store$sd[sams.for.testing, ph.leave.out] <- post.mn.loocv$sdmat[sams.for.testing, ph.leave.out]
    }
    save(loocv.store, file = loocv.res.store.namc)
  }
}



##################################################
# Fit models using Extreme Deconvolution and/or MASH
if(any(c("mash", "ed") %in% methods) & Sys.info()["sysname"] != "Windows"){
  si <- include.singletons
  bo <- run.bovy
  file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
  file.base.ed <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.ed, sapply(var.in.name.ed, get), sep = "_"), collapse = "_"))
  mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
  mash.raw.results.file.namc <- paste0(file.base.mash, "_mash_raw_results.RData")
  bovy.output.file.namc <- paste0(file.base.ed, "_bovy_output.RData")
  bovy.resl.file.namc <- paste0(file.base.ed, "_bovy_resl.RData")
  library(mashr)
  normU <- T
  use.pt.mass <- T
  ############################################################
  # Get big effects covariance matrices for MASH
  if(Data == "impc"){
    #Choose strongest effects from UV analysis
    uv.t.mat <- (Y / S)[sams.for.model.fitting, ]
    uv.t.mat <- uv.t.mat[order(-apply(abs(uv.t.mat), 1, function(v) max(v, na.rm = T))), ]
    max.abs.t <- apply(abs(uv.t.mat), 1, function(v) max(v, na.rm = T))
    t.th <- 4
    min.sams.for.Ztil <- 100
    lines.strong.use <- names(max.abs.t)[max.abs.t > t.th]
    if(length(lines.strong.use) < min.sams.for.Ztil)
      lines.strong.use <- names(max.abs.t)[1:min.sams.for.Ztil]
    big.eff.use <- lines.strong.use
    Ztil <- scale((Y / S)[big.eff.use, ph.use], scale = F)
    if(any(colMeans((Ztil == 0)) == 1)) #If phen's completely missing among big effects then use all samples
      Ztil <- scale((Y / S)[sams.for.model.fitting, ph.use], scale = F)
  }
  if(Data == "eqtl"){
    table(snpmap.sub$task)
    snps.strong.use <- snpmap.sub$snp[snpmap.sub$task == "strong.est.cov"]
    big.eff.use <- snps.strong.use
    Ztil <- scale((Y / S)[big.eff.use, ph.use], scale = F)
  }
  Sigl.init.bigeff.mash <- list()
  Sigl.init.bigeff.mash[[1]] <- cov(Ztil)
  svd.Ztil <- svd(Ztil)
  Sigl.init.bigeff.mash[[2]] <- 0
  K1 <- min(3, K)
  K2 <- min(5, K)
  for(i in 1:K1)
    Sigl.init.bigeff.mash[[2]] <- Sigl.init.bigeff.mash[[2]] + svd.Ztil$d[i] * svd.Ztil$v[, i] %*% t(svd.Ztil$v[, i])
  Sigl.init.bigeff.mash[[3]] <- 0
  for(i in 1:K2)
    Sigl.init.bigeff.mash[[3]] <- Sigl.init.bigeff.mash[[3]] + svd.Ztil$d[i] * svd.Ztil$v[, i] %*% t(svd.Ztil$v[, i])
  names(Sigl.init.bigeff.mash) <- paste0("U.", 1:3)
  
  mashdata.for.model.fitting <- mash_set_data(Bhat = Y[sams.for.model.fitting, ph.use], 
                                              Shat = S[sams.for.model.fitting, ph.use], V = R.init)
  mashdata.big.effects.for.ED <- mash_set_data(Bhat = Y[big.eff.use, ph.use], 
                                               Shat = S[big.eff.use, ph.use], V = R.init)
  
  ##############################################
  # Run Extreme Deconvolution (ED)
  if("ed" %in% methods){
    if(EDmeth == "mash"){
      mashdata.bovy <- mashdata.big.effects.for.ED
      Ul.bovy.init <- Sigl.init.bigeff.mash
    }
    if(EDmeth == "justED"){
      mashdata.bovy <- mashdata.for.model.fitting
      Ul.bovy.init <- Sigl.em.init
    }
    print("Running Extreme Deconvolution")
    replace.bovy <- T
    if(!file.exists(bovy.output.file.namc) | replace.bovy){
      bovy.out <- ed_wrapper(data = mashdata.bovy, Ulist_init = Ul.bovy.init, tol = EDtol)
      Sigl.ED <- bovy.out$Ulist <- lapply(bovy.out$Ulist, function(M){ dimnames(M) <- list(ph.use, ph.use); M})
      pimat.ED <- t(bovy.out$pi)
      save(bovy.out, mashdata.bovy, Sigl.ED, pimat.ED, file = bovy.output.file.namc)
      print(bovy.output.file.namc)
    } else {
      load(file = bovy.output.file.namc)
    }
    Sigl.ED <- lapply(Sigl.ED, function(M) M[ph.use, ph.use])
    omegaseq.ED <- 1
    resl.store <- list()
    for(dat.type in dat.typev){
      out.post.mix.ed <- em.update.function(Y.em = Y.eml[[dat.type]][sams.for.testing, ph.use], S.em = S.eml[[dat.type]][sams.for.testing, ph.use],
                                            Sigl = Sigl.ED, R = R.init, omegaseq = omegaseq.ED,
                                            pimat = pimat.ED, meth = "post.mn")
      resl.store[[dat.type]] <- list(mn = out.post.mix.ed$mnmat[sams.for.testing, ph.use], sd = out.post.mix.ed$sdmat[sams.for.testing, ph.use],
                                     loglikv = out.post.mix.ed$loglikv[sams.for.testing],
                                     lfsr = out.post.mix.ed$lfsrmat[sams.for.testing, ph.use], Sigl = Sigl.ED, Sig.mn = NULL,
                                     Ksig = NULL, R = R.init, pimat = pimat.ED, omegaseq = omegaseq.ED)
    }
    save(resl.store, file = bovy.resl.file.namc)
  }
  
  ##############################################
  # Run MASH
  if("mash" %in% methods){
    load(file = file.objl)
    Sigl.for.mash <- objl$impc_em.fit_nSig_1[[seed]]$Sigl[1]
    names(Sigl.for.mash) <- "eb"
    
    replace.bovy <- F
    if(!file.exists(bovy.output.file.namc) | replace.bovy){
      print("Running Extreme Deconvolution")
      mashdata.bovy <- mashdata.big.effects.for.ED
      Ul.bovy.init <- Sigl.init.bigeff.mash
      bovy.out <- ed_wrapper(data = mashdata.bovy, Ulist_init = Ul.bovy.init, tol = EDtol)
      Sigl.ED <- bovy.out$Ulist <- lapply(bovy.out$Ulist, function(M){ dimnames(M) <- list(ph.use, ph.use); M})
      pimat.ED <- t(bovy.out$pi)
      save(bovy.out, mashdata.bovy, Sigl.ED, pimat.ED, file = bovy.output.file.namc)
      print(bovy.output.file.namc)
    } else {
      print("Loading Extreme Deconvolution results")
      load(file = bovy.output.file.namc)
    }
    Ul.bovy.use <- bovy.out$Ulist
    if(K2 > 1){
      sfa <- varimax(svd.Ztil$v[, 1:K2])
      sfa.F <- t(svd.Ztil$v[, 1:K2] %*% sfa$rotmat)
      sfa.L <- (svd.Ztil$u[, 1:K2, drop = F] %*% diag(svd.Ztil$d[1:K2], nrow = K2, ncol = K2)) %*% sfa$rotmat
    } else {
      sfa.F <- t(svd.Ztil$v[, 1:K2])
      sfa.L <- (svd.Ztil$u[, 1:K2, drop = F] %*% diag(svd.Ztil$d[1:K2], nrow = K2, ncol = K2))
    }
    Ul.rank1 <- list()
    for(i in 1:K2)
      Ul.rank1[[i]] <- t(sfa.L[, i] %*% sfa.F[i, , drop = F]) %*% sfa.L[, i] %*% sfa.F[i, , drop = F] / nrow(Ztil)
    Ul.data <- c(Ul.bovy.use, Ul.rank1)
    names(Ul.data) <- paste0("U.", 1:(3 + K2))
    cov.meth <- c("identity", "equal_effects", "simple_het")
    if(include.singletons == 1)
      cov.meth <- c(cov.meth, "singletons")
    Ul.canon <- cov_canonical(mashdata.for.model.fitting, cov_methods = cov.meth)
    # Ul.eb <- list(Sig)
    mashmethv <- list(c("data", "canon"), c("data"), c("eb"), c("canon"), c("data", "canon", "eb"))[1]#[c(1, 2, 3, 4)]
    # #############
    # mashmethv <- list(c("data", "canon"), c("data"), c("eb"), c("canon"), c("data", "canon", "eb"))[3]#[c(1, 2, 3, 4)]
    # mashmeth <- mashmethv[[1]]
    # priorc <- "nullbiased"
    # #####
    if(!"resl" %in% ls())
      resl <- list()
    for(priorc in c("nullbiased", "uniform")[1]){#null.bias in c(T, F)){#priorc <- "nullbiased"#
    for(mashmeth in mashmethv){#mashmeth <- mashmethv[[1]]#
        if("eb" %in% mashmeth && (!"Sigl.for.mash" %in% ls() | is.null(Sigl.for.mash)))
          next
        Ul.all <- list()
        if("data" %in% mashmeth)
          Ul.all <- c(Ul.all, Ul.data)
        if("canon" %in% mashmeth)
          Ul.all <- c(Ul.all, Ul.canon)
        if("eb" %in% mashmeth)
          Ul.all <- c(Ul.all, Sigl.for.mash[1])#emout.mix$Sigl)
        res.mash.fitted.model <- mash(data = mashdata.for.model.fitting, Ulist = Ul.all, outputlevel = 2, usepointmass = use.pt.mass,
                                      prior = priorc, normalizeU = normU, add.mem.profile = F)#, grid = sqrt(omegaseq))
        if("eb" %in% mashmeth)
          eb.model <- res.mash.fitted.model
        mashdata.all.testing <- mash_set_data(Bhat = Y[sams.for.testing, ph.use], Shat = S[sams.for.testing, ph.use], alpha = 0, V = R.init)
        res.mash.all.testing <- mash(mashdata.all.testing, g = res.mash.fitted.model$fitted_g, fixg = T, outputlevel = 2, add.mem.profile = F)
        dimnames(res.mash.all.testing$vloglik) <- list(sams.for.testing, "llik")
        mashnam <- paste0("mash_", paste(mashmeth, collapse = "+"), "_prior_", priorc)
        resl[[mashnam]] <- list(mn = res.mash.all.testing$result$PosteriorMean[sams.for.testing, ph.use],
                                sd = res.mash.all.testing$result$PosteriorSD[sams.for.testing, ph.use],
                                loglikv = res.mash.all.testing$vloglik,
                                lfdr = res.mash.all.testing$result$lfdr[sams.for.testing, ph.use],
                                lfsr = res.mash.all.testing$result$lfsr[sams.for.testing, ph.use],
                                loglik = sum(res.mash.all.testing$vloglik[sams.for.lik.cross.val, ]))
        names(resl[[mashnam]]$loglikv) <- sams.for.testing
      }
    }
    phnam.mash <- colnames(res.mash.fitted.model$result$PosteriorMean)
    mash.omegaseq <- res.mash.fitted.model$fitted_g$grid^2
    mash.n.om <- length(mash.omegaseq)
    null.Sigl <- list(diag(rep(0, P)))
    names(null.Sigl) <- "null"
    mash.Sigl <- lapply(c(null.Sigl, res.mash.fitted.model$fitted_g$Ulist), function(M){ dimnames(M) <- list(ph.use, ph.use); M})
    mash.nSig <- length(mash.Sigl)
    mash.pimat.t <- matrix(0, mash.nSig, mash.n.om, 
                         dimnames = list(c("null", names(res.mash.fitted.model$fitted_g$Ulist)), 
                                                                    as.character(mash.omegaseq)))
    mash.piv <- res.mash.fitted.model$fitted_g$pi
    mash.pimat.t[2:mash.nSig, ] <- mash.piv[2:length(mash.piv)]
    mash.pimat.t[1, 1] <- mash.piv[1]
    mash.pimat.t.use <- mash.pimat.t[, , drop = F]
    mash.pimat.t.use <- mash.pimat.t.use / sum(mash.pimat.t.use)
    resl.store <- resl[names(resl) != "uv"]
    for(dat.type in dat.typev){#dat.type <- "raw"#
      out.post.mix.mash <- em.update.function(Y.em = Y.eml[[dat.type]][sams.for.testing, ph.use], 
                                                S.em = S.eml[[dat.type]][sams.for.testing, ph.use],
                                                Sigl = mash.Sigl, R = R.init,
                                                omegaseq = mash.omegaseq, prior.in.obj = F,
                                                pimat = t(mash.pimat.t.use), meth = "post.mn")
      
      resl.store[[dat.type]] <- list(mn = out.post.mix.mash$mnmat[sams.for.testing, ph.use], sd = out.post.mix.mash$sdmat[sams.for.testing, ph.use],
                                     loglikv = out.post.mix.mash$loglikv[sams.for.testing],
                                     lfsr = out.post.mix.mash$lfsrmat[sams.for.testing, ph.use], Sigl = mash.Sigl, Sig.mn = NULL,
                                     Ksig = NULL, R = R.init, pimat = t(mash.pimat.t.use), omegaseq = mash.omegaseq)
    }
    save(resl.store, file = mash.resl.file.namc)
    save(res.mash.all.testing, res.mash.fitted.model, file = mash.raw.results.file.namc)
  }
}



# 
# truemuts <- linemap$geno[linemap$line.type == "trueMut"]
# samlook <- sample(intersect(sams.for.testing, truemuts), 100)
# samlook <- intersect(sams.for.testing, truemuts)[1:100]
# load(file = file.objl)
# res.mash.fitted.model$fitted_g$grid
# mash.omegaseq
# ebnam <- "impc_em.fit_nSig_1"
# mashnam <- "impc_mash_nSig_1"
# ednam <- "impc_ed_nSig_1"
# 
# objlc <- objl[[ebnam]][[seed]]
# mash.objlc <- objl[[mashnam]][[seed]]
# phord <- colnames(objlc$Sigl[[1]])
# 
# pimat.mash <- t(mash.pimat.t.use)
# # mash.Sigl <- c(list(diag(rep(0))), 
#                
# 
# # eb.Sigl <- lapply(emout.mix$Sigl, function(M) M[phord, phord])
# # out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
# #                                    Sigl = eb.Sigl, R = emout.mix$R[phord, phord], omegaseq = emout.mix$omegaseq,
# #                                    pimat = emout.mix$pi, meth = "post.mn")
# eb.Sigl <- lapply(objlc$Sigl, function(M) M[phord, phord])
# out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
#                                    Sigl = eb.Sigl, R = objlc$R[phord, phord], omegaseq = objlc$omegaseq,
#                                    pimat = objlc$pimat, meth = "post.mn")
# mash.indsin <- 1
# # colnames(pimat.mash)
# pimat.mash <- objl[[mashnam]][[seed]]$pimat
# pimat.eb <- objl[[ebnam]][[seed]]$pimat
# mashind <- colnames(pimat.mash)[order(-colSums(pimat.mash))][mash.indsin]
# # mashind <- setdiff(colnames(pimat.mash), phord)
# # mashind <- paste0("U.", 1:8)
# # mashind <- colnames(pimat.mash)[colSums(pimat.mash) > 0]
# # mashind <- colnames(pimat.mash)[!grepl("U.", colnames(pimat.mash))]
# mashind
# 
# # mash.Sigl <- lapply(mash.objlc$Sigl, function(M) M[phord, phord] / 1)
# # mash.Sigl$U.1 <- eb.Sigl[[1]] / max(diag(eb.Sigl[[1]]))
# mash.out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
#                                         Sigl = lapply(mash.Sigl, function(M) M[phord, phord] / 1), 
#                                         R = R.init[phord, phord], omegaseq = mash.omegaseq,
#                                         pimat = pimat.mash, meth = "post.mn")
# dim(pimat.mash)
# str(mash.Sigl)
# # mash.out.post.mix <- em.update.function(Y.em = Yhat[samlook, phord], S.em = smat[samlook, phord],
# #                                         Sigl = mash.objlc$Sigl, R = mash.objlc$R, omegaseq = mash.objlc$omegaseq,
# #                                         pimat = mash.objlc$pimat, meth = "post.mn")
# # plot(compl$impc_em.fit_nSig_1$llmat.raw[samlook, splitc], out.post.mix$loglikv[samlook])
# eblik <- out.post.mix$loglikv[samlook]
# mashlik <- mash.out.post.mix$loglikv[samlook]
# plot(eblik, mashlik)
# abline(0, 1, col = 2)
# mean(eblik)
# mean(mashlik)
# 
# plot(mash.omegaseq, pimat.mash[, 2], ty = "l", log = "x")
# lines(objlc$omegaseq * max(diag(eb.Sigl[[1]])), objlc$pimat[, 1], col = 2)
# 
# 
# 
# # rm(list = ls())
# Data <- "eqtl"
# xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
# source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
# source(paste0(R.file.dir, "/impc_mv_paper_code/EM_fns_lean.R"))
# 
# load(file = file.resl.comp)
# load(file = file.compl)
# load(file = uv.results.Y.S)
# 
# 
# 
# 
# 
# llmean.splitmat.raw <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[truemuts, ], na.rm = T))
# llmean.splitmat.zero <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.zero[truemuts, ], na.rm = T))
# colMeans(llmean.splitmat.raw, na.rm = T)
# colMeans(llmean.splitmat.zero, na.rm = T)
# plot(llmean.splitmat.raw[, mashnam], llmean.splitmat.raw[, ebnam])
# abline(0, 1)
# str(compl,m=2)
# str(compl[[1]]$llmat.raw, m = 2)
# 
# 
# si <- 1
# bo <- 1
# N <- 2000
# P <- 148
# EDtol <- 1e-5
# EDmeth <- c("justED", "mash")[2]
# splitc <- seed <- 5
# nSig <- 1
# Data <- "impc"
# file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
# file.base.ed <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.ed, sapply(var.in.name.ed, get), sep = "_"), collapse = "_"))
# mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
# mash.raw.results.file.namc <- paste0(file.base.mash, "_mash_raw_results.RData")
# load(mash.raw.results.file.namc)
# 
# # str(objl$impc_mash_nSig_1[[seed]], m = 1)
# # sort(colSums(objl$impc_mash_nSig_1[[seed]]$pimat))
# # plot(sort(colSums(objl$impc_mash_nSig_1[[seed]]$pimat)))
# # plot(sort(colSums(objl$impc_em.fit_nSig_1[[seed]]$pimat)))
# 
# ebnam <- "impc_em.fit_nSig_1"
# mashnam <- "impc_mash_nSig_1"
# ednam <- "impc_ed_nSig_1"
# 
# 
# load(file = paste0(meth.comp.output.dir, "/N_2000_P_148_nSig_1_seed_", seed, "_data_impc_sexspecific_FALSE_emout.RData"))
# load(file = paste0(meth.comp.output.dir, "/N_2000_P_148_nSig_1_seed_", seed, "_data_impc_sexspecific_FALSE_res.RData"))
# 
# str(res.store)
# str(emout.mix)
# str(emout.mix)
# all(samlook %in% rownames(res.store$mn))
# # plot(objl[[ebnam]][[seed]]$Sighat, objl[[mashnam]][[seed]]$Sighat)
# # abline(0, 1)
# plot(diag(objl[[ebnam]][[seed]]$Sighat), diag(objl[[mashnam]][[seed]]$Sighat))
# abline(0, 1)
# mean(diag(objl[[ebnam]][[seed]]$Sighat) / diag(objl[[mashnam]][[seed]]$Sighat))
# plot(diag(objl[[ebnam]][[seed]]$Sighat), diag(objl[[ednam]][[seed]]$Sighat))
# abline(0, 1)
# sams.for.testing <- objl$impc_em.fit_nSig_1[[seed]]$saml$sams.for.testing
# # sams.for.testing <- objl[[ednam]][[seed]]$saml$sams.for.testing
# 
# 
# objlc <- objl[[ebnam]][[seed]]
# mash.objlc <- objl[[mashnam]][[seed]]
# phord <- colnames(objlc$Sigl[[1]])
# # samlook <- sample(intersect(intersect(sams.for.testing, truemuts), rownames(res.store$mn)), 100)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# plot(diag(objlc$Sigl[[1]][phord, phord]), diag(Sig.init[phord, phord]))
# 
# 











# 
# ###########################################################################################
# # Note: must have run "collect_results.R" prior to running factor score estimation
# if(any(c("em.fac", "em.all.test") %in% methods)){
#   print("Running em.fac")
#   load(file = emout.file.namc)
#   # load(file = file.glob.res) # From "collect_results.R"
#   load(file = file.glob.loadings) # 'facs' From "collect_results.R"
#   fac.res.store <- list()
#   for(fac.meth in c("varimax", "promax")){
#     facs <- switch(fac.meth, varimax = facs.varimax, promax = facs.promax)
#     fac.out.post.mix <- em.update.function(Y.em = Yhat[sams.for.testing, ph.use], S.em = smat[sams.for.testing, ph.use],
#                                            Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
#                                            pimat = emout.mix$pi, meth = "post.mn.fac", loadings = facs[ph.use, facnam])
#     fac.res.store[[fac.meth]] <- list(mn = fac.out.post.mix$mnmat[sams.for.testing, facnam], sd = fac.out.post.mix$sdmat[sams.for.testing, facnam],
#                                       loglikv = fac.out.post.mix$loglikv[sams.for.testing],
#                                       lfsr = fac.out.post.mix$lfsrmat[sams.for.testing, facnam], 
#                                       loadings = facs[ph.use, facnam])
#   }
#   save(fac.res.store, file = fac.res.store.namc)
# }
# 



# 
# str(res.mash.fitted.model)
# 
# 
# 
# str(mashd)
# 
# 
# si <- include.singletons
# bo <- run.bovy
# file.base.mash <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
# mash.resl.file.namc <- paste0(file.base.mash, "_mash_resl.RData")
# mash.raw.results.file.namc <- paste0(file.base.mash, "_mash_raw_results.RData")
# load(mash.raw.results.file.namc)
# 
# str(res.mash.fitted.model)
# cvsams.look <- sams.for.lik.cross.val[1:3]
# mashdata.all.testing.sub <- mash_set_data(Bhat = Y[cvsams.look, ph.use], Shat = S[cvsams.look, ph.use], alpha = 0, V = R.init)
# xUlist <- mashr:::expand_cov(res.mash.fitted.model$fitted_g$Ulist, grid = res.mash.fitted.model$fitted_g$grid, usepointmass = T)
# lm <- mashr:::calc_relative_lik_matrix(data = mashdata.all.testing.sub, Ulist = xUlist, algorithm.version = c("R", "Rcpp")[2])
# vloglik <- mashr:::compute_vloglik_from_matrix_and_pi(pi_s = res.mash.fitted.model$fitted_g$pi, lm = lm, mashdata.all.testing.sub$Shat_alpha)
# 
# mash.omegaseq <- res.mash.fitted.model$fitted_g$grid^2
# mash.n.om <- length(mash.omegaseq)
# mash.Sigl <- c(list(diag(rep(0, P))), lapply(res.mash.fitted.model$fitted_g$Ulist, function(M) M[ph.use, ph.use]))
# mash.nSig <- length(mash.Sigl)
# mash.pimat.t <- matrix(0, mash.nSig, mash.n.om, dimnames = list(c("null", names(res.mash.fitted.model$fitted_g$Ulist)), as.character(mash.omegaseq)))
# mash.piv <- res.mash.fitted.model$fitted_g$pi
# mash.pimat.t[2:mash.nSig, ] <- mash.piv[2:length(mash.piv)]
# mash.pimat.t[1, 1] <- mash.piv[1]
# 
# 
# str(res.mash.fitted.model$fitted_g$Ulist)
# out.post.mix.mash <- em.update.function(Y.em = Y[cvsams.look, ph.use], S.em = S[cvsams.look, ph.use],
#                                         Sigl = mash.Sigl, R = R.init[ph.use, ph.use],
#                                         omegaseq = mash.omegaseq, prior.in.obj = F,
#                                         pimat = t(mash.pimat.t), meth = "just.obj")
# 
# out.post.mix.mash$llikv
# 
# 
# plot(vloglik, out.post.mix.mash$llikv)
# plot(vloglik, res.mash.all.testing$vloglik[cvsams.look, ])
# 
# str(xUlist)
# 
# resl[[mashnam]]$loglik
# plot(obj.out$llikv, res.mash.all.testing$vloglik[sams.for.lik.cross.val, ])
# abline(0, 1)
# plot(obj.out$llikv, obj.out.ed$llikv)
# abline(0, 1)
# plot(obj.out.ed$llikv, res.mash.all.testing$vloglik[sams.for.lik.cross.val, ],
#      col = ifelse(res.mash.all.testing$vloglik[sams.for.lik.cross.val, ] > obj.out.ed$llikv, 2, 1))
# abline(0, 1)
# plot(obj.out$llikv, res.mash.all.testing$vloglik[sams.for.lik.cross.val, ],
#      col = ifelse(res.mash.all.testing$vloglik[sams.for.lik.cross.val, ] > obj.out$llikv, 2, 1))
# abline(0, 1)
# 
# 
# str(res.mash.fitted.model)
# 
# 
# sum(obj.out$llikv)
# sum(obj.out.ed$llikv)
# sum(res.mash.all.testing$vloglik[sams.for.lik.cross.val, ])
# 
# 
# 
# Y.em = Yhat[sams.for.lik.cross.val, ph.use]; S.em = smat[sams.for.lik.cross.val, ph.use]
# Sigl = mash.Sigl; R = R.init;
# omegaseq = mash.omegaseq;
# pimat = t(mash.pimat.t); meth = "just.obj"
# 
# out.post.mix.mash <- em.update.function(Y.em = Yhat[sams.for.lik.cross.val, ph.use], S.em = smat[sams.for.lik.cross.val, ph.use],
#                                    Sigl = mash.Sigl, R = R.init,
#                                    omegaseq = mash.omegaseq, prior.in.obj = F,
#                                    pimat = t(mash.pimat.t), meth = "just.obj")
# out.post.mix.eb <- em.update.function(Y.em = Yhat[sams.for.lik.cross.val, ph.use], S.em = smat[sams.for.lik.cross.val, ph.use],
#                                         Sigl = emout.mix$Sigl, R = R.init,
#                                         omegaseq =emout.mix$omegaseq, prior.in.obj = F,
#                                         pimat = emout.mix$pi, meth = "just.obj")
# 
# load(file = uv.results.Y.S)
# S[S == 10] <- 5
# out.post.mix.mash <- em.update.function(Y.em = Y[sams.for.lik.cross.val, ph.use], S.em = S[sams.for.lik.cross.val, ph.use],
#                                         Sigl = mash.Sigl, R = R.init,
#                                         omegaseq = mash.omegaseq, prior.in.obj = F,
#                                         pimat = t(mash.pimat.t), meth = "just.obj")
# out.post.mix.eb <- em.update.function(Y.em = Y[sams.for.lik.cross.val, ph.use], S.em = S[sams.for.lik.cross.val, ph.use],
#                                       Sigl = emout.mix$Sigl, R = R.init,
#                                       omegaseq =emout.mix$omegaseq, prior.in.obj = F,
#                                       pimat = emout.mix$pi, meth = "just.obj")
# str(out.post.mix.mash)
# plot(out.post.mix.eb$llikv, out.post.mix.mash$llikv,
#      col = ifelse(out.post.mix.mash$llikv > out.post.mix.eb$llikv, 2, 1))
# abline(0, 1)
# mean(out.post.mix.eb$llikv, na.rm = T)
# mean(out.post.mix.mash$llikv, na.rm = T)
# str(res.mash.fitted.model$fitted_g$Ulist)
# 
# mash.pimat.t
# 
# 
# plot(obj.out$llikv, res.mash.all.testing$vloglik[sams.for.lik.cross.val, ],
#   col = ifelse(res.mash.all.testing$vloglik[sams.for.lik.cross.val, ] > obj.out$llikv, 2, 1))
# abline(0, 1)
# 
# mean(out.post.mix.eb$llikv - out.post.mix.mash$llikv, na.rm = T)
# mean(obj.out$llikv - res.mash.all.testing$vloglik[sams.for.lik.cross.val, ], na.rm = T)

# 
# fvc <- list.files("/mnt/c/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/RData_files/mv_results/methods_comparison")
# fvcc <- fvc[grep("bovy", fvc)]
# load(file = "/mnt/c/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/RData_files/mv_results/methods_comparison/N_500_P_20_si_1_bo_1_seed_42_data_impc_bovy_output.RData")
# str(bovy.out)


# ##################################################
# # Fit models using extreme deconvolution
# if("ed" %in% methods & Sys.info()["sysname"] != "Windows"){
#   # Sigl.for.mash <- unlist(lapply(resl[names(resl) != "uv"], function(x) x$Sigl), recursive = F)
#   # Sigl.for.mash <- unlist(lapply(resl[1][names(resl[1]) != "uv"], function(x) x$Sigl), recursive = F)
#   # file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
#   
#   ed.ulist.init <- list(Sigl.em.init)
#   library(mashr)
#   mashdata.for.model.fitting <- mash_set_data(Bhat = Y[sams.for.model.fitting, ph.use], Shat = S[sams.for.model.fitting, ph.use],
#                                               alpha = 0, V = R.init)#diag(rep(1, P)))
#   ed.out <- ed_wrapper(data = mashdata.for.model.fitting, Ulist_init = ed.ulist.init) 
#   ed.out <- cov_ed(data = mashdata.for.model.fitting, Ulist_init = ed.ulist.init)
#   
#   if (is.null(subset)) {
#     subset = 1:n_effects(data)
#   }
#   K = length(Ulist_init)
#   R = n_conditions(data)
#   pi_init = rep(1/K, K)
#   D = ncol(data$V)
#   if (all(data$V == diag(D))) {
#     ed.res = extreme_deconvolution(data$Bhat[subset, ], data$Shat[subset, 
#                                                                   ]^2, xamp = pi_init, xmean = matrix(0, nrow = K, 
#                                                                                                       ncol = R), xcovar = Ulist_init, fixmean = TRUE, ...)
#   }
#   else {
#     ycovar = lapply(subset, function(i) data$Shat[i, ] * 
#                       t(data$V * data$Shat[i, ]))
#     ed.res = extreme_deconvolution(data$Bhat[subset, ], ycovar, 
#                                    xamp = pi_init, xmean = matrix(0, nrow = K, ncol = R), 
#                                    xcovar = Ulist_init, fixmean = TRUE, ...)
#   }
#   
#   
  # 
  # 
  # 
  # ycovar.in <- array(NA, dim = c(N, P, P), dimnames = list(sams.for.model.fitting, ph.use, ph.use))
  # ycovar.in <- list()
  # for(samc in sams.for.model.fitting){
  #   
  #   ycovar.in[[samc]] <-  t(R.init * S[samc, ph.use]) * S[samc, ph.use]
  # 
  # extreme_deconvolution(ydata = Y[sams.for.testing, ph.use], ycovar = ycovar.in, 
  #                       xamp = 1, xmean = rep(0, P), xcovar = Sig.init,
  #                       projection = NULL, weight = NULL,
  #                       fixamp = NULL,fixmean=NULL,fixcovar=NULL,
  #                       tol=1.e-6,maxiter=1e9,w=0,logfile=NULL,
  #                       splitnmerge=0,maxsnm=FALSE,likeonly=FALSE,
  #                       logweight=FALSE)
  
  




























# 
# 
# 
# corout <- t(t(Sig.mn) / sqrt(diag(Sig.mn))) / sqrt(diag(Sig.mn))
# nfac.pl <- min(floor(P / 4), nfac)
# loadmat <- varimax(svd(corout)$v[, 1:nfac.pl])$loadings
# loadmat <- sweep(loadmat, 2, apply(loadmat, 2, function(v) v[which.max(abs(v))]), '/')
# loadmat <- loadmat[, order(colMeans(abs(loadmat) * 1:nrow(loadmat)))]



# 
# if(fac){
#   ##################################################################
#   # Calculate loadings from eigendecomposition of Sig
#   Sig <- emout.mix$Sig.mn
#   R <- emout.mix$R
#   dimnames(R) <- dimnames(Sig) <- list(ph.use, ph.use)
#   Sigdiag <- diag(Sig)
#   Sigcor <- t(Sig / sqrt(Sigdiag)) / sqrt(Sigdiag)
#   eigc <- eigen(Sigcor)
#   cumpropvar <- cumsum(eigc$values) / sum(eigc$values)
#   nfac <- match(T, cumpropvar > .9)
#   vc.type <- c("vari", "pro")[1]
#   if(vc.type == "vari")
#     facs <- varimax(eigen(Sigcor)$vectors[, 1:nfac])$loadings
#   if(vc.type == "pro")
#     facs <- promax(eigen(Sigcor)$vectors[, 1:nfac])$loadings
#   facnam <- paste("f", 1:nfac, sep = ".")
#   dimnames(facs) <- list(ph.use, facnam)
#   out.post.mix.fac <- em.update.function(Yin[sams.for.testing, ph.use], S.em = Sin[sams.for.testing, ph.use],
#                            Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq, pimat = emout.mix$pi,
#                            meth = "post.mn.fac", loadings = facs)
# }
# 























# Sigl.for.mash <- unlist(lapply(resl[names(resl) != "uv"], function(x) x$Sigl), recursive = F)
# Sigl.for.mash <- unlist(lapply(resl[1][names(resl[1]) != "uv"], function(x) x$Sigl), recursive = F)
# file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, sapply(var.in.name.mash, get), sep = "_"), collapse = "_"))
# file.base <- paste0(meth.comp.output.dir, "/", paste(paste(var.in.name.mash, tab[1, var.in.name.mash], sep = "_"), collapse = "_"))
# file.namc <- paste0(meth.comp.output.dir, paste(paste(var.in.name, tab[scen, var.in.name], sep = "_"), collapse = "_"), ".RData")
# library(RhpcBLASctl)
# blas_set_num_threads(1)
# Sig.init <- cov(Y[sams.for.model.fitting, ph.use], use = "p", meth = "p")
# Sig.init[is.na(Sig.init)] <- 0
# diag(Sig.init)[is.na(diag(Sig.init))] <- 1
# if(qr(Sig.init)$rank < P)
#   Sig.init <- (1 - ident.eps) * Sig.init + ident.eps * diag(rep(1, P))
# table(snpmap.sub$task)
# snpmap.sub <- snpmap.sub[order(snpmap.sub$random.ordering), ]
# sams.for.cor.est <- snpmap.sub$snp[snpmap.sub$task == "random.train.est.cor"]
# sams.for.strong.cov.est <- snpmap.sub$snp[snpmap.sub$task == "strong.est.cov "]
# sams.for.model.fitting <- snpmap.sub$snp[snpmap.sub$task == "random.train.fit.model"]
# sams.for.testing.random <- snpmap.sub$snp[which(snpmap.sub$task == "random.test.model")]
# sams.for.testing.random <- sams.for.testing.random[1:min(length(sams.for.testing.random), floor(max.num.sams.for.testing / 2))]
# sams.for.testing.strong <- snpmap.sub$snp[which(snpmap.sub$task == "strong.test.model")]
# sams.for.testing.strong <- sams.for.testing.strong[1:min(length(sams.for.testing.strong), floor(max.num.sams.for.testing / 2))]
# sams.for.testing <- c(sams.for.testing.random, sams.for.testing.strong)
# sams.for.lik.cross.val <- sams.for.testing[snpmap.sub[match(sams.for.testing, snpmap.sub$snp), "task"] == "random.test.model"]
# sams.for.cor.est <- linemap.sub$geno[linemap.sub$line.type == "negConTra"]
# sams.for.model.fitting <- linemap.sub$geno[linemap.sub$line.type == "trueMutTra"]
# sams.for.testing <- linemap.sub$geno[linemap.sub$line.type %in% c("negConVal", "trueMutVal", "negConTes", "trueMutTes")]
# sams.for.testing <- sams.for.testing[1:min(c(length(sams.for.testing), max.num.sams.for.testing))]
# sams.for.lik.cross.val <- sams.for.testing[linemap.sub[match(sams.for.testing, linemap.sub$geno), "line.type"] %in%
#                                              c("trueMutVal", "trueMutTes")]

# # ###########
# poss.existing.files <- c(bovy.output.file.namc, gsub("si_0", "si_1", bovy.output.file.namc), gsub("si_1", "si_0", bovy.output.file.namc))
# f.existing <- file.exists(poss.existing.files)
# if(any(f.existing)){
#   print("Loaded Bovy output")
#   fc <- poss.existing.files[match(T, f.existing)]
#   load(file = fc)
# } else {
#   stop("Need to run Extreme Deconvolution for this parameter combination, before running mash!!")
# }
# bol[[as.character(EDtol)]] <- bovy.out
# }
# str(bol)
# max(bol[[1]]$pi - bol[[2]]$pi)
# max(bol[[2]]$pi - bol[[3]]$pi)
# max(bol[[1]]$Ulist$U.1 - bol[[2]]$Ulist$U.1)
# max(bol[[2]]$Ulist$U.1 - bol[[3]]$Ulist$U.1)
# max(bol[[1]]$Ulist$U.2 - bol[[2]]$Ulist$U.2)
# max(bol[[2]]$Ulist$U.3 - bol[[3]]$Ulist$U.3)
# max(bol[[1]]$Ulist$U.3 - bol[[2]]$Ulist$U.3)
# max(bol[[2]]$Ulist$U.3 - bol[[3]]$Ulist$U.3)
# # #########

#   raw.out.post.mix.ed <- em.update.function(Y.em = Yhat[sams.for.testing, ph.use], S.em = smat[sams.for.testing, ph.use],
#                                    Sigl = bovy.out$Ulist, R = R.init, omegaseq = omegaseq.ED,
#                                    pimat = pimat.ED, meth = "post.mn")
# zerofilled.out.post.mix.ed <- em.update.function(Y.em = Y[sams.for.testing, ph.use], S.em = S[sams.for.testing, ph.use],
#                                           Sigl = bovy.out$Ulist, R = R.init, omegaseq = omegaseq.ED,
#                                           pimat = pimat.ED, meth = "post.mn")
# # raw.obj.out.ed <- em.update.function(Y.em = Yhat[sams.for.lik.cross.val, ph.use], S.em = smat[sams.for.lik.cross.val, ph.use],
# #                               Sigl = bovy.out$Ulist, R = R.init, omegaseq = omegaseq.ED,
# #                               pimat = pimat.ED, meth = "just.obj", prior.in.obj = F)
# # zero.filled.obj.out.ed <- em.update.function(Y.em = Y[sams.for.lik.cross.val, ph.use], S.em = S[sams.for.lik.cross.val, ph.use],
# #                                  Sigl = bovy.out$Ulist, R = R.init, omegaseq = omegaseq.ED,
# #                                  pimat = pimat.ED, meth = "just.obj", prior.in.obj = F)
# res.store.raw <- list(mn = raw.outpost.mix.ed$mnmat[sams.for.testing, ph.use], sd = raw.out.post.mix.ed$sdmat[sams.for.testing, ph.use],
#                       raw.loglikv = raw.out.post.mix.ed$loglikv[sams.for.testing],
#                       lfsr = raw.out.post.mix.ed$lfsrmat[sams.for.testing, ph.use], Sigl = bovy.out$Ulist, Sig.mn = NULL,
#                       Ksig = P, R = R.init, piv = pimat.ED, omegaseq = omegaseq.ED)
# res.store.zerofilled <- list(mn = zerofilled.outpost.mix.ed$mnmat[sams.for.testing, ph.use], sd = zerofilled.out.post.mix.ed$sdmat[sams.for.testing, ph.use],
#                              zerofilled.loglikv = zerofilled.out.post.mix.ed$loglikv[sams.for.testing],
#                              lfsr = zerofilled.out.post.mix.ed$lfsrmat[sams.for.testing, ph.use], Sigl = bovy.out$Ulist, Sig.mn = NULL,
#                              Ksig = P, R = R.init, piv = pimat.ED, omegaseq = omegaseq.ED)
# res.store.zerofilled <- list(mn = out.post.mix.ed$mnmat[sams.for.testing, ph.use], sd = out.post.mix.ed$sdmat[sams.for.testing, ph.use],
#                       raw.loglikv = out.post.mix.ed$loglikv[sams.for.testing],
#                       lfsr = out.post.mix.ed$lfsrmat[sams.for.testing, ph.use], Sigl = bovy.out$Ulist, Sig.mn = NULL,
#                       Ksig = P, R = R.init, piv = pimat.ED, omegaseq = omegaseq.ED, 
#                       zerofilled.loglikv = obj.out.ed$llikv, loglik = -obj.out.ed$obj)
# mashdata.lik.cross.val.testing <- mash_set_data(Bhat = Y[sams.for.lik.cross.val, ph.use],
#                                                 Shat = S[sams.for.lik.cross.val, ph.use], alpha = 0, V = R.init)
# res.mash.fitted.model$fitted_g$pi[1] <- 0
# res.mash.fitted.model$fitted_g$pi <- res.mash.fitted.model$fitted_g$pi / sum(res.mash.fitted.model$fitted_g$pi)
# res.mash.lik.cross.val.testing <- mash(mashdata.lik.cross.val.testing, g = res.mash.fitted.model$fitted_g, fixg = T, outputlevel = 2)
# if("eb" %in% mashmeth)
#   eb.model.cross.val <- res.mash.lik.cross.val.testing
# if(mashnam == "mash_eb")
#   my.pi.sub.res.mash <- res.mash.true.use.train
# ed_wrapper
# extreme_deconvolution

# if(any(c("em.test", "em.all.test", "em.all") %in% methods)){
#   print(nSig)
#   print("Running em.test")
#   load(file = emout.file.namc)
#   resl.store <- list()
#   for(dat.type in dat.typev){
#     out.post.mix <- em.update.function(Y.em = Y.eml[[dat.type]][sams.for.testing, ph.use], S.em = S.eml[[dat.type]][sams.for.testing, ph.use],
#                                           Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
#                                           pimat = emout.mix$pi, meth = "post.mn")
#     resl.store[[dat.type]] <- list(mn = out.post.mix$mnmat[sams.for.testing, ph.use], sd = out.post.mix$sdmat[sams.for.testing, ph.use],
#                       loglikv = out.post.mix$loglikv[sams.for.testing],
#                       lfsr = out.post.mix$lfsrmat[sams.for.testing, ph.use], 
#                       Sigl = lapply(emout.mix$Sigl, function(M) M[ph.use, ph.use]), 
#                       Sig.mn = emout.mix$Sigmn[ph.use, ph.use],
#                       Ksig = emout.mix$Ksig, R = emout.mix$R[ph.use, ph.use], pimat = emout.mix$pi, omegaseq = omegaseq)
#   }
#   save(resl.store, file = res.store.namc)
#   # out.post.mix <- em.update.function(Y.em = Yhat[sams.for.testing, ph.use], S.em = smat[sams.for.testing, ph.use],
#   #                                    Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
#   #                                    pimat = emout.mix$pi, meth = "post.mn")
#   # obj.out <- em.update.function(Y.em = Y[sams.for.lik.cross.val, ph.use], S.em = S[sams.for.lik.cross.val, ph.use],
#   #                               Sigl = emout.mix$Sigl, R = emout.mix$R, omegaseq = emout.mix$omegaseq,
#   #                               pimat = emout.mix$pi, meth = "just.obj", prior.in.obj = F)
#   # res.store <- list(mn = out.post.mix$mnmat[sams.for.testing, ph.use], sd = out.post.mix$sdmat[sams.for.testing, ph.use],
#   #                   loglikv = out.post.mix$loglikv[sams.for.testing],
#   #                   lfsr = out.post.mix$lfsrmat[sams.for.testing, ph.use], Sigl = emout.mix$Sigl, Sig.mn = emout.mix$Sigmn,
#   #                   Ksig = emout.mix$Ksig, R = emout.mix$R, pimat = emout.mix$pi, omegaseq = omegaseq, 
#   #                   zerofilled.loglikv = obj.out$llikv)
#   # save(res.store, file = res.store.namc)
# }

#Load ED output
# poss.existing.files <- c(bovy.output.file.namc, gsub("si_0", "si_1", bovy.output.file.namc), gsub("si_1", "si_0", bovy.output.file.namc))
# f.existing <- file.exists(poss.existing.files)
# if(any(f.existing)){
#   print("Loaded Bovy output")
#   fc <- poss.existing.files[match(T, f.existing)]
#   load(file = fc)
# } else {
#   stop("Need to run Extreme Deconvolution for this parameter combination, before running mash!!")
# }
