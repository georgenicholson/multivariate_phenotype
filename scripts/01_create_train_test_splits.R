
get_test_train_split <- function(control, Data, N, P, phens_to_use) { 

# for(Data in c("impc", "eqtl")){
  if(Data == "impc"){
    Yhat <- readRDS("data/impc/Yhatmat.RDS")
    smat <- readRDS("data/impc/Shatmat.RDS")
    linemap <- readRDS("data/impc/linemap.RDS")
    reflinemap <- readRDS("data/impc/reflinemap.RDS")
    phmap <- readRDS("data/impc/phmap.RDS")
    P_full <- ncol(Yhat)
    N_full <- nrow(Yhat)
    N_max_train <- floor(sum(linemap$line.type == "trueMut") / 2)
  }
  if(Data == "eqtl"){
    mash.in <- readRDS("data/eqtl/MatrixEQTLSumStats.Portable.ld2.Z.rds")
    Yhat <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
    smat <- rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
    P_full <- ncol(Yhat)
    N_full <- nrow(Yhat)
    N_max_train <- nrow(mash.in$random.z) / 2
  }
  # Pseqsub <- c(control$Pseq[control$Pseq <= P_full], control$default_parameters[[Data]]$P)
  # Nseqsub <- c(control$Nseq[control$Nseq <= N_max_train], control$default_parameters[[Data]]$N)
  train_test_list <- list()
  train_test_list$sams_for_cor_est <- train_test_list$sams_for_model_training <- 
    train_test_list$sams_for_model_testing <- train_test_list$sams_for_lik_cross_val <- 
    matrix(FALSE, nrow = N_full, ncol = control$n_subsamples, dimnames = list(rownames(Yhat), 1:control$n_subsamples))
  # train_test_list$phens_to_use <- matrix(FALSE, nrow = P_full, ncol = control$n_subsamples, dimnames = list(colnames(Yhat), 1:control$n_subsamples))
  # for(N in Nseqsub){#N <- max(Nseqsub)#Nseqsub[1]#
    # for(P in c(P_full, Pseqsub)){#P <- P_full#  
      for(subsample_seed in 1:control$n_subsamples){
          if(Data == "eqtl"){
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
            set.seed(subsample_seed)
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
            set.seed(subsample_seed)
            # phens_to_use <- sample(colnames(Y), P)
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
            # suppressWarnings(save(snpmap.sub, phens_to_use, 
            #                       sams.for.cor.est, sams.for.model.fitting, sams.for.testing, sams.for.lik.cross.val, 
            #                       file = file.out, version = 1))
          }
          if(Data == "impc"){
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
            set.seed(subsample_seed)
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
            # if(P < P_full){
            #   # Make sure no phen's are completely missing at sampled genes
            #   ok.phs <- phmap$ph[colMeans(Y[lines.use, phmap$ph] == 0) < .9]
            #   phens_to_use <- sample(ok.phs, P)
            # } else {
            #   phens_to_use <- phmap$ph
            # }
            # phens_to_use <- phens_to_use[order(phmap[match(phens_to_use, phmap$ph), "procnam"])]
            sams.for.cor.est <- linemap.sub$geno[linemap.sub$line_type_subsample == "negConTra"]
            sams.for.model.fitting <- linemap.sub$geno[linemap.sub$line_type_subsample == "trueMutTra"]
            sams.for.testing <- linemap.sub$geno[linemap.sub$line_type_subsample %in% c("negConTes", "trueMutTes")]
            sams.for.lik.cross.val <- sams.for.testing[linemap.sub[match(sams.for.testing, linemap.sub$geno), "line_type_subsample"] %in%
                                                         c("trueMutTes")]
            # suppressWarnings(save(linemap.sub, phens_to_use, 
            #                       sams.for.cor.est, sams.for.model.fitting, sams.for.testing, sams.for.lik.cross.val, 
            #                       file = file.out, version = 1))
          }   
        train_test_list$sams_for_cor_est[sams.for.cor.est, subsample_seed] <- TRUE
        train_test_list$sams_for_model_training[sams.for.model.fitting, subsample_seed] <- TRUE
        train_test_list$sams_for_model_testing[sams.for.testing, subsample_seed] <- TRUE
        train_test_list$sams_for_lik_cross_val[sams.for.lik.cross.val, subsample_seed] <- TRUE
        # train_test_list$phens_to_use[phens_to_use, subsample_seed] <- TRUE
        # train_test_list$sams_for_cor_est[sams.for.cor.est, subsample_seed] <- TRUE
        # train_test_list$sams_for_model_training[sams.for.model.fitting, subsample_seed] <- TRUE
        # train_test_list$sams_for_model_testing[sams.for.testing, subsample_seed] <- TRUE
        # train_test_list$sams_for_lik_cross_val[sams.for.lik.cross.val, subsample_seed] <- TRUE
        # train_test_list$phens_to_use[phens_to_use, subsample_seed] <- TRUE
        
      }
      # saveRDS(object = train_test_list, file = file.path(control$train_test_samples_dir, paste0(Data, "_Ntrain_", N, "_P_", P, ".RDS")))
    return(train_test_list)
    }
#   }
# }
