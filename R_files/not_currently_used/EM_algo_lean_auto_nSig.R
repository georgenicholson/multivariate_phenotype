EM_algo_mixture_multi_Sig <- function(Y.em, S.em, Sigl.em.init, R.em.init, omegaseq = 1,#exp((-10):10 * log(2)), 
                            rel.tol = .001, it.start.rel.tol = 5, K.init = NULL,
                            seed = 1, update.Sig = T, update.K = T, update.nSig = F, pi.init = NULL){
  M <- length(omegaseq)
  R <- R.em.init
  N <- nrow(Y.em)
  ph.use <- colnames(Y.em)
  P <- length(ph.use)
  auto.tol.const.while.optimizing.nSig <- 1e-3
  auto.tol.const.final <- 5e-6
  new.mix.component.pi.th <- 0.05
  
  
  #################################
  # Initialise
  Sigl <- Sigl.em.init
  nSig <- length(Sigl)
  if(is.null(pi.init)){
    pimat <- matrix(1, M, nSig)
    pimat <- pimat / sum(pimat)
  } else {
    pimat <- pi.init
  }
  if(is.null(K.init)){
    Ksig <- rep(P, nSig)
  } else {
    Ksig <- K.init
  }
  
  ###################################
  # Start EM algorithm
  objv <- c()
  itnum <- 1
  # new.nSig <- nSig
  final.converged <- F
  if(update.nSig){
    stagec <- "opt.nSig"
  }  else {
    stagec <- "final"
  }
  while(!final.converged){
  # for(stagec in c("opt.nSig", "final")){
    auto.tol.const <- switch(stagec, opt.nSig = auto.tol.const.while.optimizing.nSig, final = auto.tol.const.final)
    converged <- F
    while(!converged){
      #############################################################
      # Compute updated Sigl, rmat
      em.up.out <- em.update.function(Y.em = Y.em, S.em = S.em, Sigl = Sigl, R = R, omegaseq = omegaseq, 
                                 pimat = pimat, meth = "update.Sigl", Ksig = Ksig, update.K = update.K)
      if(update.Sig)
        Sigl <- em.up.out$Sigl
      if(update.K)
        Ksig <- em.up.out$Ksig
      rmat <- em.up.out$rmat
      
      ############################################
      # Optimize wrt pi
      diralpha.new <- apply(rmat, 2:3, sum)
      # pimat[] <- diralpha.new / sum(diralpha.new)
      diralpha.shrink <- 0#.125
      pimat[] <- pmax(diralpha.new - diralpha.shrink, 0) / sum(pmax(diralpha.new - diralpha.shrink, 0))
      
      #############################
      # Evaluate objective function
      calc.obj.out <- em.update.function(Y.em = Y.em, S.em = S.em, Sigl = Sigl, R = R, omegaseq = omegaseq,
                         pimat = pimat, Ksig = Ksig, meth = "just.obj", update.Sig = update.Sig)
      objv[itnum] <- calc.obj.out$obj
      llikv <- calc.obj.out$llikv
      
      mix.keep <- which(colSums(pimat) > new.mix.component.pi.th)
      if(length(mix.keep) < nSig){
        stagec <- "final"
        auto.tol.const <- auto.tol.const.final
      }
      pimat <- pimat[, mix.keep, drop = F]
      Sigl <- Sigl[mix.keep]
      Ksig <- Ksig[mix.keep]
      rmat <- rmat[, , mix.keep, drop = F]
      nSig <- length(mix.keep)
      
      ##################################
      # Check convergence and print trace
      if(length(objv) > 1){
        most.recent.change <- objv[itnum] - objv[itnum - 1]
        llik.scale <- mad(llikv, na.rm = T)
        tolerance.eps <- llik.scale * N * auto.tol.const
        if(most.recent.change < 0 & most.recent.change > -tolerance.eps){
          converged <- TRUE
          if(stagec == "final")
            final.converged <- T
        }
        dob <- objv[length(objv)] - objv[length(objv) - 1]
        print(paste("Iteration = ", itnum))
        print(paste0("Ksig = ", Ksig))
        print(paste0("Obj = ", objv[length(objv)]))
        print(paste0("Change in obj = ", dob))
        print("Pi = ")
        print(round(t(pimat[, order(-colSums(pimat))]) * 100))
        print(paste0("llik.scale = ", llik.scale))
        print(paste0("tolerance.eps = ", tolerance.eps))
        print(paste0("most.recent.change = ", most.recent.change))
      }
      itnum <- itnum + 1
      
    }
    if(stagec == "opt.nSig"){
    # mix.keep <- which(pisum > new.mix.component.pi.th)
    # if(length(mix.keep) < nSig){
    #   pimat <- pimat[, mix.keep]
    #   Sigl <- Sigl[mix.keep]
    #   Ksig <- Ksig[mix.keep]
    #   rmat <- rmat[, , mix.keep]
    #   nSig <- length(mix.keep)
    #   stagec <- "final"
    #   # next
    # } else {
      nSig <- nSig + 1
      # Seed new mixture component
      lq.export <- .25
      uq.export <- .5
      prob.export <- .5
      llik.qv <- quantile(llikv, c(lq.export, uq.export), na.rm = T)
      sams.to.export <- which(llikv >= llik.qv[1] & llikv <= llik.qv[2])
      Sig.init <- cov(Y.em[sams.to.export, ], use = "p", meth = "p")
      Sig.init[is.na(Sig.init)] <- 0
      diag(Sig.init)[is.na(diag(Sig.init))] <- 1
      if(qr(Sig.init)$rank < P)
        Sig.init <- (1 - ident.eps) * Sig.init + ident.eps * diag(rep(1, P))
      Sigl[[nSig]] <- Sig.init
      Ksig[nSig] <- max(Ksig)
      rmat.new.component <- array(0, dim = c(N, M, nSig), dimnames = list(rownames(Y.em), 1:M, 1:nSig))
      for(i in 1:N){
        if(i %in% sams.to.export){
          rmat.new.component[i, , 1:(nSig - 1)] <- rmat[i, , ] * (1 - prob.export)
          rmat.new.component[i, , nSig] <- prob.export / M
        } else {
          rmat.new.component[i, , 1:(nSig - 1)] <- rmat[i, , ]
        }
      }
      pimat.new.component <- cbind(pimat, rep(0, M))
      diralpha.new.component <- apply(rmat.new.component, 2:3, sum)
      pimat <- diralpha.new.component / sum(diralpha.new.component)
      # Y.em = Y.em; S.em = S.em; Sigl = Sigl; R = R; omegaseq = omegaseq;
      # pimat = NULL; meth = "new.component"; Ksig = Ksig; update.K = update.K; rmat = rmat.new.component
      em.up.out <- em.update.function(Y.em = Y.em, S.em = S.em, Sigl = Sigl, R = R, omegaseq = omegaseq, 
                                      pimat = pimat, meth = "new.component", Ksig = Ksig, 
                                      update.K = F, rmat = rmat.new.component)
      Sigl <- em.up.out$Sigl
      str(Sigl)
    }
  }
  Sig.mn <- 0
  for(sc in 1:nSig){
    for(m in 1:M)
      Sig.mn <- Sig.mn + pimat[m, sc] * Sigl[[sc]] * omegaseq[m]
  }
  return(list(Sigl = Sigl, Sig.mn = Sig.mn, Ksig = Ksig, R = R, pi = pimat, omegaseq = omegaseq, objv = objv))
}






































































