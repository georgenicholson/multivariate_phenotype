for(dataset in c("impc", "eqtl")){
  if(dataset == "impc"){
    Yhat <- readRDS("data/impc/Yhatmat.RDS")
    smat <- readRDS("data/impc/Shatmat.RDS")
    linemap <- readRDS("data/impc/linemap.RDS")
    reflinemap <- readRDS("data/impc/reflinemap.RDS")
    phmap <- readRDS("data/impc/phmap.RDS")
    P_full <- ncol(Yhat)
    N_full <- nrow(Yhat)
    N_max_train <- floor(sum(linemap$line.type == "trueMut") / 2)
  }
  if(dataset == "eqtl"){
    mash.in <- readRDS("data/eqtl/MatrixEQTLSumStats.Portable.ld2.Z.rds")
    Yhat <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
    P_full <- ncol(Yhat)
    N_full <- nrow(Yhat)
    N_max_train <- nrow(mash.in$random.z) / 2
  }
  Pseqsub <- c(control$Pseq[control$Pseq <= P_full], control$default_parameters[[dataset]]$P)
  Nseqsub <- c(control$Nseq[control$Nseq <= N_max_train], control$default_parameters[[dataset]]$N)
  train_test_list <- list()
  train_test_list$sams_for_cor_est <- train_test_list$sams_for_model_training <- 
    train_test_list$sams_for_model_testing <- train_test_list$sams_for_lik_cross_val <- 
    matrix(FALSE, nrow = N_full, ncol = control$n_subsamples, dimnames = list(rownames(Yhat), 1:control$n_subsamples))
  train_test_list$phens_to_use <- matrix(FALSE, nrow = P_full, ncol = control$n_subsamples, dimnames = list(colnames(Yhat), 1:control$n_subsamples))
  for(N in Nseqsub){#N <- max(Nseqsub)#Nseqsub[1]#
    for(P in c(P_full, Pseqsub)){#P <- P_full#  
      for(seed_curr in 1:control$n_subsamples){
          if(dataset == "eqtl"){
            Yhat <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
            smat <- Yhat / rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
            Y <- Yhat
            S <- smat
            Y[is.na(Yhat)] = 0
            S[is.na(Yhat)] = control$prior.sd.on.unobserved.thetahat
            data.type <- rep(c("strong", "random.train", "random.test"), 
                             times = c(nrow(mash.in$strong.z), nrow(mash.in$random.z), nrow(mash.in$random.test.z)))
            snpmap <- data.frame(snp = rownames(Yhat), data.type = data.type)
            # tasks <- c("random.est.cor", "strong.est.cov", "random.fit.model", "random.test.model", "strong.test.model")
            set.seed(seed_curr)
            snpmap <- snpmap[sample(1:nrow(snpmap)), ]
            snpmap$random.ordering <- 1:nrow(snpmap)
            # snpmap[snpmap$data.type == "strong", "task"] <- "random.est.cor"
            strong.all <- snpmap$snp[snpmap$data.type == "strong"]
            strong.est.cov <- sample(strong.all)[1:floor(length(strong.all) / 2)]
            strong.test.model <- setdiff(strong.all, strong.est.cov)
            random.train <- snpmap$snp[snpmap$data.type == "random.train"]
            random.train.est.cor <- sample(random.train)[1:floor(length(random.train) / 2)]
            random.train.fit.model <- setdiff(random.train, random.train.est.cor)
            snpmap$task[snpmap$snp %in% strong.est.cov] <- "strong.est.cov"
            snpmap$task[snpmap$snp %in% strong.test.model] <- "strong.test.model"
            snpmap$task[snpmap$snp %in% random.train.est.cor] <- "random.train.est.cor"
            snpmap$task[snpmap$snp %in% random.train.fit.model] <- "random.train.fit.model"
            snpmap$task[snpmap$data.type == "random.test"] <- "random.test.model"
            snpmap <- snpmap[rowSums(is.na(snpmap)) == 0, ]
            # Subsample model fit data to sample size N
            strong.est.cov.sub <- sample(snpmap$snp[snpmap$task == "strong.est.cov"], min(N, sum(snpmap$task == "strong.est.cov")))
            random.train.fit.model.sub <- sample(snpmap$snp[snpmap$task == "random.train.fit.model"], min(N, sum(snpmap$task == "random.train.fit.model")))
            random.train.est.cor.sub <- sample(snpmap$snp[snpmap$task == "random.train.est.cor"], min(N, sum(snpmap$task == "random.train.est.cor")))
            # snpmap$snp[snpmap$task == "random.train.est.cor"]#
            # Not currently subsampling test sets
            random.test.model.sub <- snpmap$snp[snpmap$task == "random.test.model"]#sample(snpmap$snp[snpmap$task == "random.test.model"], min(N, sum(snpmap$task == "random.test.model")))
            strong.test.model.sub <- snpmap$snp[snpmap$task == "strong.test.model"]#sample(snpmap$snp[snpmap$task == "random.test.model"], min(N, sum(snpmap$task == "random.test.model")))
            subsnps <- sample(c(strong.est.cov.sub, random.train.est.cor.sub, random.train.fit.model.sub, random.test.model.sub, strong.test.model.sub))
            snpmap.sub <- snpmap[match(subsnps, snpmap$snp), ]
            set.seed(seed_curr)
            ph.use <- sample(colnames(Y), P)
            snpmap.sub <- snpmap.sub[order(snpmap.sub$random.ordering), ]
            sams.for.cor.est <- snpmap.sub$snp[snpmap.sub$task == "random.train.est.cor"]
            sams.for.strong.cov.est <- snpmap.sub$snp[snpmap.sub$task == "strong.est.cov "]
            sams.for.model.fitting <- snpmap.sub$snp[snpmap.sub$task == "random.train.fit.model"]
            sams.for.testing.random <- snpmap.sub$snp[which(snpmap.sub$task == "random.test.model")]
            # sams.for.testing.random <- sams.for.testing.random[1:min(length(sams.for.testing.random), floor(max.num.sams.for.testing / 2))]
            sams.for.testing.strong <- snpmap.sub$snp[which(snpmap.sub$task == "strong.test.model")]
            # sams.for.testing.strong <- sams.for.testing.strong[1:min(length(sams.for.testing.strong), floor(max.num.sams.for.testing / 2))]
            sams.for.testing <- c(sams.for.testing.random, sams.for.testing.strong)
            sams.for.lik.cross.val <- sams.for.testing[snpmap.sub[match(sams.for.testing, snpmap.sub$snp), "task"] == "random.test.model"]
            # suppressWarnings(save(snpmap.sub, ph.use, 
            #                       sams.for.cor.est, sams.for.model.fitting, sams.for.testing, sams.for.lik.cross.val, 
            #                       file = file.out, version = 1))
          }
          if(dataset == "impc"){
            # load(file = uv.results.Y.S)
            ##############################################
            #Arrange data, and appropriately fill in values to be imputed
            Y <- Yhat
            S <- smat
            Y[is.na(Yhat)] = 0
            S[is.na(Yhat)] = control$prior.sd.on.unobserved.thetahat
            linemap$old.line.type <- linemap$line.type
            linemap <- linemap[match(rownames(Yhat), linemap$geno), ]
            linemap$geno.nozyg <- sapply(strsplit(linemap$geno, spl = "_"), function(x) x[1])

            ###############################################
            # Randomly arrange data into four categories, trueMutTra, trueMutTes, negConTra, negConTes
            set.seed(seed_curr)
            negcon.genzyg <- sample(linemap[linemap$old.line.type != "trueMut", "geno"])
            truemut.genzyg <- sample(linemap[linemap$old.line.type == "trueMut", "geno"])
            
            #Keep ref lines for testing true lines
            refline.genos <- linemap[linemap$geno.nozyg %in% reflinemap$genotype_id & linemap$geno %in% truemut.genzyg, "geno"]
            truemut.genzyg <- truemut.genzyg[order(truemut.genzyg %in% refline.genos)]#  (move ref lines to end)
          
            n.null.lines <- length(negcon.genzyg)
            linemap[match(negcon.genzyg[1:N], linemap$geno), "line_type_subsample"] <- "negConTra"
            
            
            sum(is.na(match(negcon.genzyg[1:N], linemap$geno)))
            
            linemap[match(negcon.genzyg[(N + 1):n.null.lines], linemap$geno), "line_type_subsample"] <- "negConTes"
            n.true.lines <- length(truemut.genzyg)
            linemap[match(truemut.genzyg[1:N], linemap$geno), "line_type_subsample"] <- "trueMutTra"
            linemap[match(truemut.genzyg[(N + 1):n.true.lines], linemap$geno), "line_type_subsample"] <- "trueMutTes"
            linemap <- linemap[rowSums(is.na(linemap)) == 0, ]

            neg.tra.use <- linemap$geno[linemap$line_type_subsample == "negConTra"]
            true.tra.use <- linemap$geno[linemap$line_type_subsample == "trueMutTra"]
            neg.tes.use <- linemap$geno[linemap$line_type_subsample == "negConTes"]
            true.tes.use <- linemap$geno[linemap$line_type_subsample == "trueMutTes"]
            
            # Randomly order lines
            lines.use <- sample(c(neg.tra.use, true.tra.use, neg.tes.use, true.tes.use))
            linemap.sub <- linemap[match(lines.use, linemap$geno), ]
            if(P < P_full){
              # Make sure no phen's are completely missing at sampled genes
              ok.phs <- phmap$ph[colMeans(Y[lines.use, phmap$ph] == 0) < .9]
              ph.use <- sample(ok.phs, P)
            } else {
              ph.use <- phmap$ph
            }
            ph.use <- ph.use[order(phmap[match(ph.use, phmap$ph), "procnam"])]
            sams.for.cor.est <- linemap.sub$geno[linemap.sub$line_type_subsample == "negConTra"]
            sams.for.model.fitting <- linemap.sub$geno[linemap.sub$line_type_subsample == "trueMutTra"]
            sams.for.testing <- linemap.sub$geno[linemap.sub$line_type_subsample %in% c("negConTes", "trueMutTes")]
            sams.for.lik.cross.val <- sams.for.testing[linemap.sub[match(sams.for.testing, linemap.sub$geno), "line_type_subsample"] %in%
                                                         c("trueMutTes")]
            # suppressWarnings(save(linemap.sub, ph.use, 
            #                       sams.for.cor.est, sams.for.model.fitting, sams.for.testing, sams.for.lik.cross.val, 
            #                       file = file.out, version = 1))
          }   
        train_test_list$sams_for_cor_est[sams.for.cor.est, seed_curr] <- TRUE
        train_test_list$sams_for_model_training[sams.for.model.fitting, seed_curr] <- TRUE
        train_test_list$sams_for_model_testing[sams.for.testing, seed_curr] <- TRUE
        train_test_list$sams_for_lik_cross_val[sams.for.lik.cross.val, seed_curr] <- TRUE
        train_test_list$phens_to_use[ph.use, seed_curr] <- TRUE
        
      }
      saveRDS(object = train_test_list, file = file.path(control$train_test_samples_dir, paste0(dataset, "_Ntrain_", N, "_P_", P, ".RDS")))
  
    }
  }
}

    # 
    # true.tra.use <- snpmap.sub$snp[snpmap.sub$task == "random.train.fit.model"]
    # lines.val.test.use <- snpmap.sub$snp[snpmap.sub$task == "random.test.model"]

# linemap[linemap$geno %in% refline.genos, "line.type"] <- "trueMutTes"

# #Add matching het/homs to testing true lines
# tabtemp <- table(linemap$geno.nozyg[linemap$old.line.type == "trueMut"])
# geno.nozyg.hethom <- names(tabtemp)[tabtemp == 2]
# geno.hethom <- linemap[linemap$geno.nozyg %in% geno.nozyg.hethom & linemap$old.line.type == "trueMut", "geno"]
# geno.hethom.use <- geno.hethom[!geno.hethom %in% true.tra.use]
# linemap[linemap$geno %in% geno.hethom.use, "line.type"] <- "trueMutTes"

#Add ref lines and hets/homs to testing true lines
# true.tes.use <- unique(c(true.tes.use, refline.genos, geno.hethom.use))
# set.seed(seed_curr)
# neg.tra.use <- sample(neg.tra.all, min(N, length(neg.tra.all)))
# true.tra.use <- sample(true.tra.all, min(N, length(true.tra.all)))
# Not currently subsampling validation/test sets here
# neg.tes.use <- neg.tes.all#sample(neg.tes.all, min(ceiling(max.lines.val.test / 4), length(neg.tes.all)))
# true.tes.use <- true.tes.all#sample(true.tes.all, min(ceiling(max.lines.val.test / 4), length(true.tes.all)))
# lines.val.test.use <- c(neg.val.use, neg.tes.use, true.val.use, true.tes.use)

# #######################################
# #Gather arguments from command line
# arguments <- commandArgs()
# if("--args" %in% arguments){
#   seed_curr <- as.numeric(arguments[grep("--args", arguments) + 1])
#   N <- as.numeric(arguments[grep("--args", arguments) + 2])
#   P <- as.numeric(arguments[grep("--args", arguments) + 3])
#   run.bovy <- as.numeric(arguments[grep("--args", arguments) + 4])
#   include.singletons <- as.numeric(arguments[grep("--args", arguments) + 5])
#   file.out <- as.character(arguments[grep("--args", arguments) + 6])
#   run.specific.scenario <- F
#   run.mash <- T
#   max.lines.val.test <- 2000
#   full.analysis <- F
#   do.plots <- F
# } else {
#   run.specific.scenario <- T
# }
# 
# if(run.specific.scenario){
#   # seed_curr <- 7575
#   seed_curr <- as.numeric(format(Sys.time(), "%s")) %% 30000
#   do.plots <- T
#   full.analysis <- F
#   if(full.analysis){
#     N <- 1500
#     P <- 148
#     max.lines.val.test <- 1000
#     nSig <- 1
#     run.mash <- F
#     run.bovy <- 0
#     include.singletons <- 0
#   }
#   if(!full.analysis){
#     # N <- 1500
#     # P <- 148
#     N <- 1000
#     P <- 40
#     max.lines.val.test <- 1000
#     run.mash <- T
#     run.bovy <- 1
#     include.singletons <- 1
#     nSig <- 1
#   }
#   K <- floor(P / 2)
# }

# estimate_null_correlation_simple <- function (data, z_thresh = 2, est_cor = TRUE) 
# {
#   z = data$Bhat/data$Shat
#   max_absz = apply(abs(z), 1, max)
#   nullish = which(max_absz < z_thresh)
#   if (length(nullish) < ncol(z)) {
#     stop("not enough null data to estimate null correlation")
#   }
#   nullish_z = z[nullish, ]
#   Vhat = cov(nullish_z)
#   if (est_cor) {
#     Vhat = cov2cor(Vhat)
#   }
#   return(Vhat)
# }

# if(Sys.info()["nodename"] != "TW")
#   source(paste0(ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:"),
#                 "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/EM_obj_grad_functions.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/EM_algo_mixture.R"))
# 
# ##############################################
# #Initial values for R, Sig
# initialize.at <- c("Identity", "Summary.stats")[2]
# if(initialize.at == "Identity"){
#   R.init <- Sig.init <- diag(rep(1, P))
# }
# if(initialize.at == "Summary.stats"){
#   estimate_null_correlation_simple
#   if(dataset == "impc"){
#     if(cor.type == "weighted"){
#       R.init <- cor(Y[neg.tra.use, ph.use] / S[neg.tra.use, ph.use], use = "p", meth = "p")
#     }
#     if(cor.type == "unweighted"){
#       R.init <- cor(Y[neg.tra.use, ph.use], use = "p", meth = "p")
#     }
#     if(cor.type == "identity")
#       R.init <- diag(1, P, P)
#     R.init[is.na(R.init)] <- 0
#     diag(R.init) <- 1
#     ident.eps <- .05
#     if(qr(R.init)$rank < P)
#       R.init <- (1 - ident.eps) * R.init + ident.eps * diag(rep(1, P))
#     Sig.init[is.na(Sig.init)] <- 0
#     diag(Sig.init)[is.na(diag(Sig.init))] <- 1
#     if(qr(Sig.init)$rank < P)
#       Sig.init <- (1 - ident.eps) * Sig.init + ident.eps * diag(rep(1, P))
#     Sig.init <- cov(Y[true.tra.use, ph.use], use = "p", meth = "p")
#   }
#   if(dataset == "eqtl"){
#     snps.est.cor <- snpmap.sub$snp[snpmap.sub$task == "random.train.est.cor"]
#     data.in <- list(Bhat = Y[snps.est.cor, ph.use], Shat = S[snps.est.cor, ph.use])
#     R.init <- estimate_null_correlation_simple(data = data.in, z_thresh = 2)
#     Sig.init <- cov(Y[snps.est.cor, ph.use], use = "p", meth = "p")
#   }
# }
# dimnames(R.init) <- dimnames(Sig.init) <- list(ph.use, ph.use)
# Sig.em.init <- Sig.init[ph.use, ph.use]
# R.em.init <- R.init[ph.use, ph.use]
# 
# ##############################
# # Run EM mixture
# if(!"resl" %in% ls())
#   resl <- list()
# resl$uv <- list(mn = Yhat[lines.val.test.use, ph.use], sd = smat[lines.val.test.use, ph.use])
# resl$uv$lfsr <- pnorm(-abs(resl$uv$mn) / resl$uv$sd)
# if(!full.analysis){
#   # Kseq <- c(1, 2, 5, 10, 20, 50)[1]
#   nSigseq <- c(1, 2, 5, 25)[1]#5#c(1, 2, 5, 10, 15)[c(5)]
#   Kseq <- floor(P / 4)#c(Kseq[Kseq < P], P)
#   omegamixseq <- c(T, F)[1]
#   calcfacseq <- c(T, F)[2]
#   miss.seq <- c("em", "prior")[1]
#   alphaminseq <- c(.1, 1, 2, 10)[2]
#   Rident.trainseq <- c(T, F)[2]
#   Rident.testseq <- c(T, F)[2]
#   methseq <- c("rank1", "full", "robust")[2]
#   random.init <- F
#   pi.null <- .2
# }
# if(full.analysis){
#   Kseq <- P
#   nSigseq <- nSig
#   omegamixseq <- c(T, F)[1]
#   calcfacseq <- c(T, F)[2]
#   miss.seq <- c("em", "prior")[1]
#   alphaminseq <- c(.1, 1, 5, 10)[4]
#   Rident.trainseq <- c(T, F)[2]
#   Rident.testseq <- c(T, F)[2]
#   methseq <- c("rank1", "full")[2]
#   random.init <- F
#   pi.null <- .2
# }
# tab <- expand.grid(N = N, P = P, K = Kseq, omix = omegamixseq, nSig = nSigseq, 
#                    miss = miss.seq, fac = calcfacseq, amin = alphaminseq, Rid.tr = Rident.trainseq, 
#                    Rid.te = Rident.testseq, meth = methseq, stringsAsFactors = F)
# 
# tab
# source(paste0(R.file.dir, "/impc_mv_analysis/EM_algo_mixture_multi_Sig.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/EM_algo_mixture_multi_Sig_functions.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/EM_algo_mixture_multi_Sig_robust.R"))
# source(paste0(R.file.dir, "/impc_mv_analysis/EM_algo_mixture_multi_Sig_rank1.R"))
# omegaseq.full <- 1 * exp(-10:10 * log(sqrt(2)))
# # omegaseq.full <- 1 * exp(-3:3 * log(2))
# scen=1#for(scen in 1:nrow(tab)){
# for(j in 1:ncol(tab))
#   assign(names(tab)[j], tab[scen, j])
# namc <- paste(c(t(cbind(names(tab), unlist(tab[scen, ])))), collapse = "_")
# if(omix){
#   omegaseq <- omegaseq.full
# } else {
#   omegaseq <- 1
# }
# if(miss == "em"){
#   Yin <- Yhat
#   Sin <- smat
# }
# if(miss == "prior"){
#   Yin <- Y
#   Sin <- S
# }
# if(dataset == "eqtl"){
# }
# keep.parameter.paths <- F
# diralpha <- 1
# Ksig <- rep(Kseq, length.out = nSig)
# Y.em <- Yin[true.tra.use, ph.use]
# S.em <- Sin[true.tra.use, ph.use]
# # Y.em <- Yin[c(true.tra.use, neg.tra.use), ph.use]
# # S.em <- Sin[c(true.tra.use, neg.tra.use), ph.use]
# rel.tol <- 1e-3 
# it.start.rel.tol <- 5
# alpha.start <- amin# * 10
# sd.large <- 10
# monitor.fdr <- T
# Sig.em.init = Sig.em.init;
# R.em.init = R.em.init; diralpha = diralpha; 
# # omegaseq = omegaseq.use;
# rel.tol = rel.tol; it.start.rel.tol = it.start.rel.tol; do.plots = do.plots; update.Sig = T;
# missing.data.treat = miss; nSig = nSig; Ksig = Ksig; seed = NULL;
# random.init = random.init; alpha.start = alpha.start; alpha.min = amin
# if(Rid.tr){
#   Rtrain <- diag(rep(1, P))
# } else {
#   Rtrain <- R.em.init
# }
# if(Rid.tr){
#   Rtest <- diag(rep(1, P))
# } else {
#   Rtest <- R.em.init
# }