
get_train_test_split <- function(control, Data_all, Data, N, P, n_subsamples) { 
  ###################################
  # Randomly choose phen subset for demo
  train_test_list <- list()
  if (P < Data_all[[Data]]$P_all) {
    phens_ok_to_use <- colnames(Data_all[[Data]]$Y_raw)[colMeans(is.na(Data_all[[Data]]$Y_raw)) < .8]
    if (length(phens_ok_to_use) < P) {
      phens_ok_to_use <- Data_all[[Data]]$meas_names
    }
    set.seed(1)
    phens_to_use <- sample(x = phens_ok_to_use, size = P)
    train_test_list$phens_to_use <- phens_to_use
  } else {
    train_test_list$phens_to_use <- Data_all[[Data]]$meas_names
  }
  empty_matrix <- matrix(FALSE, 
                         nrow = Data_all[[Data]]$N_all, 
                         ncol = n_subsamples, 
                         dimnames = list(rownames(Data_all[[Data]]$Y_raw), 
                                         1:n_subsamples))
  names_to_init <- c("sams_for_cor_est", 
                     "sams_for_model_training", 
                     "sams_for_model_testing", 
                     "sams_for_lik_cross_val", 
                     "sams_for_strong_cov_est")
  for (nam_init in names_to_init) {
    train_test_list[[nam_init]] <- empty_matrix
  }

  for(subsample_seed in 1:n_subsamples){
      if(Data == "eqtl"){
        set.seed(subsample_seed)
        snpmap <- Data_all$eqtl$snpmap
        snpmap <- snpmap[sample(1:nrow(snpmap)), ]
        snpmap$random.ordering <- 1:nrow(snpmap)
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
        # Not currently subsampling test sets
        random.test.model.sub <- snpmap$snp[snpmap$task == "random.test.model"]#sample(snpmap$snp[snpmap$task == "random.test.model"], min(N, sum(snpmap$task == "random.test.model")))
        strong.test.model.sub <- snpmap$snp[snpmap$task == "strong.test.model"]#sample(snpmap$snp[snpmap$task == "random.test.model"], min(N, sum(snpmap$task == "random.test.model")))
        subsnps <- sample(c(strong.est.cov.sub, random.train.est.cor.sub, random.train.fit.model.sub, random.test.model.sub, strong.test.model.sub))
        snpmap.sub <- snpmap[match(subsnps, snpmap$snp), ]
        set.seed(subsample_seed)
        snpmap.sub <- snpmap.sub[order(snpmap.sub$random.ordering), ]
        sams.for.cor.est <- snpmap.sub$snp[snpmap.sub$task == "random.train.est.cor"]
        sams.for.strong.cov.est <- snpmap.sub$snp[snpmap.sub$task == "strong.est.cov"]
        sams.for.model.fitting <- snpmap.sub$snp[snpmap.sub$task == "random.train.fit.model"]
        sams.for.testing.random <- snpmap.sub$snp[which(snpmap.sub$task == "random.test.model")]
        sams.for.testing.strong <- snpmap.sub$snp[which(snpmap.sub$task == "strong.test.model")]
        sams.for.testing <- c(sams.for.testing.random, sams.for.testing.strong)
        sams.for.lik.cross.val <- sams.for.testing[snpmap.sub[match(sams.for.testing, snpmap.sub$snp), "task"] == "random.test.model"]
        train_test_list$sams_for_strong_cov_est[sams.for.strong.cov.est, subsample_seed] <- TRUE
      }
      if(Data == "impc"){
        linemap <- Data_all$impc$linemap
        linemap$old.line.type <- linemap$line.type
        linemap <- linemap[match(rownames(Data_all[[Data]]$Y_raw), linemap$geno), ]
        linemap$geno.nozyg <- sapply(strsplit(linemap$geno, spl = "_"), function(x) x[1])

        ###############################################
        # Randomly arrange data into four categories, trueMutTra, trueMutTes, negConTra, negConTes
        set.seed(subsample_seed)
        negcon.genzyg <- sample(linemap[linemap$old.line.type != "trueMut", "geno"])
        truemut.genzyg <- sample(linemap[linemap$old.line.type == "trueMut", "geno"])
        #Keep ref lines for testing true lines
        refline.genos <- linemap[linemap$geno.nozyg %in% Data_all$impc$reflinemap$genotype_id & linemap$geno %in% truemut.genzyg, "geno"]
        truemut.genzyg <- truemut.genzyg[order(truemut.genzyg %in% refline.genos)]#  (move ref lines to end)
        n.null.lines <- length(negcon.genzyg)
        
        
        linemap[match(negcon.genzyg[1:N], linemap$geno), "line_type_subsample"] <- "negConTra"
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
        sams.for.cor.est <- linemap.sub$geno[linemap.sub$line_type_subsample == "negConTra"]
        sams.for.model.fitting <- linemap.sub$geno[linemap.sub$line_type_subsample == "trueMutTra"]
        sams.for.testing <- linemap.sub$geno[linemap.sub$line_type_subsample %in% c("negConTes", "trueMutTes")]
        sams.for.lik.cross.val <- sams.for.testing[linemap.sub[match(sams.for.testing, linemap.sub$geno), "line_type_subsample"] %in%
                                                     c("trueMutTes")]
        train_test_list$linemap.sub <- linemap.sub
      }   
    train_test_list$sams_for_cor_est[sams.for.cor.est, subsample_seed] <- TRUE
    train_test_list$sams_for_model_training[sams.for.model.fitting, subsample_seed] <- TRUE
    train_test_list$sams_for_model_testing[sams.for.testing, subsample_seed] <- TRUE
    train_test_list$sams_for_lik_cross_val[sams.for.lik.cross.val, subsample_seed] <- TRUE
  }
  return(train_test_list)
}
