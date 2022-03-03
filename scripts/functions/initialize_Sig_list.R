initialize_Sig_list <- function (control, Y, nSig, random_init = F) {
  sams_for_model_training <- rownames(Y)
  phens_to_use <- colnames(Y)
  if (random_init) {
    Sigl.em.init <- list()
    for(j in 1:nSig){
      Sig.temp <- solve(rWishart(1, P, diag(rep(1, P)))[, , 1])
      Sigl.em.init[[j]] <- t(Sig.temp / sqrt(diag(Sig.temp))) / sqrt(diag(Sig.temp))
    }
  } else {
    require(mclust)
    if (nSig > 1) {
      mcl.out <- mclust::Mclust(data = Y, G = nSig)
      cluster.membership <- apply(mcl.out$z, 1, which.max)
    } else {
      cluster.membership <- rep(1, length(sams_for_model_training))
    }
    names(cluster.membership) <- sams_for_model_training
    Sigl.em.init <- list()
    for (j in 1:nSig) {
      sams.in <- names(cluster.membership)[which(cluster.membership == j)]
      Sig.init <- cov(Y[sams.in, ], use = "p", meth = "p")
      Sig.init[is.na(Sig.init)] <- 0
      diag(Sig.init)[is.na(diag(Sig.init))] <- 1
      if(qr(Sig.init)$rank < P)
        Sig.init <- (1 - control$rank_deficient_R_eps) * Sig.init + control$rank_deficient_R_eps * diag(rep(1, P))
      dim(Sig.init)
      dimnames(Sig.init) <- list(phens_to_use, phens_to_use)
      Sigl.em.init[[j]] <- Sig.init
    }
  }
  Sigl.em.init
}
