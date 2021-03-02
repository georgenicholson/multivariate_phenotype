calc.log.mvnorm.dens <- function(yg, sg, Sig, R){
  obsinds <- which(!is.na(yg[]))
  if(length(obsinds) > 0){
    vc <- Sig[obsinds, obsinds, drop = F] + t(t(R[obsinds, obsinds, drop = F]) * sg[obsinds]) * sg[obsinds]
    logmvndensc <- -.5 * (determinant(vc, logarithm = T)$modulus + length(obsinds) * log(2 * pi)
                          + t(yg[obsinds]) %*% solve(vc) %*% yg[obsinds])
    return(logmvndensc)
  } else {
    return(NA)
  }
}

em.update.function <- function(Y.em, S.em, Sigl, R, omegaseq, pimat, Ksig, 
                                   meth = c("just.obj", "update.Sigl", "post.mn", "post.mn.fac")[2],
                                   loadings = NULL, prior.in.obj = T, update.Sig, update.K = F, max.K.change = .1,
                                   rmat = NULL, fac.model = c("fa", "pca")[2], bic.pen.mult = .5, llmat = NULL,
                                   recalc.llmat = T, wish.pri = T){
  nSig <- length(Sigl)
  M <- length(omegaseq)
  P <- ncol(Y.em)
  N <- nrow(Y.em)
  include.bic <- T
  wish.pri.nu <- P
  if(is.null(llmat)){
    recalc.llmat <- T
    Sig.inds.calc <- 1:nSig
    llmat <- array(-Inf, dim = c(N, M, nSig), dimnames = list(rownames(Y.em), 1:M, 1:nSig))
  } else {
    Sig.inds.calc <- which(update.Sig)
  }
  if(meth %in% c("update.Sigl", "post.mn", "post.mn.fac")){
    if(wish.pri){
      Ctl <- rep(list(diag(rep(1, P))), nSig)
      btl <- rep(list(2 * wish.pri.nu + 1), nSig)
    } else {
      Ctl <- rep(list(diag(rep(0, P))), nSig)
      btl <- rep(list(0), nSig)
    }
    rmat <- array(NA, dim = c(N, M, nSig), dimnames = list(rownames(Y.em), 1:M, 1:nSig))
    loglikv <- array(NA, dim = c(N), dimnames = list(rownames(Y.em)))
    if(meth == "post.mn.fac"){
      if(is.null(loadings))
        stop("Argument 'loadings' must be supplied if meth == 'post.mn.fac'")
      K <- ncol(loadings)
      mnmat <- sdmat <- lfsrmat <- varmat <- array(NA, dim = c(N, K), dimnames = list(rownames(Y.em), colnames(loadings)))
    } else {
      mnmat <- sdmat <- lfsrmat <- varmat <- array(NA, dim = dim(Y.em), dimnames = dimnames(Y.em))
    }
    t0 <- Sys.time()
    pi.thresh <- 1e-10
    for(i in 1:nrow(Y.em)){
      propdone <- i / nrow(Y.em)
      cat(paste0("\rCalculating likelihood matrix; estimated time remaining = ", 
                 ceiling(difftime(Sys.time(), t0, units = "mins") * (1 - propdone) / propdone), " mins "))
      mn.add <- var.add <- moment2.add <- lfsrmat.add <- 0
      if(any(!is.na(Y.em[i, ]))){
        # logmvndensc <- matrix(-Inf, M, nSig)
        for(sc in Sig.inds.calc){        # for(sc in 1:nSig){
          for(m in 1:M){
            if(pimat[m, sc] > pi.thresh){
              Sigc <- Sigl[[sc]] * omegaseq[m]
              llmat[i, m, sc] <- calc.log.mvnorm.dens(yg = Y.em[i, ], sg = S.em[i, ], Sig = Sigc, R = R)
            }
          }
        }
        logmvndensc <- llmat[i, , ]
        maxlogdens <- max(logmvndensc, na.rm = T)
        pi.normlik.mat <- pimat * exp(logmvndensc - maxlogdens)
        loglikv[i] <- log(sum(pi.normlik.mat, na.rm = T)) + maxlogdens
        rmat[i, , ] <- pi.normlik.mat / sum(pi.normlik.mat)
        for(sc in Sig.inds.calc){
          for(m in 1:M){
            if(rmat[i, m, sc] > 0){
              Sigc <- Sigl[[sc]] * omegaseq[m]
              obsinds <- which(!is.na(Y.em[i, ]))
              if(length(obsinds) == 0){
                mnc <- rep(0, P)
                V <- Sigc
              } else {
                Vinvobs <- solve(Sigc[obsinds, obsinds, drop = F] + t(R[obsinds, obsinds, drop = F] * S.em[i, obsinds]) * S.em[i, obsinds])
                tmp <- Sigc[, obsinds, drop = F] %*% Vinvobs
                mnc <- tmp %*% Y.em[i, obsinds]
                V <- Sigc - tmp %*% Sigc[obsinds, , drop = F]
              }
              if(meth == "post.mn"){
                mn.add <- mn.add + rmat[i, m, sc] * mnc
                moment2.add <- moment2.add + rmat[i, m, sc] * mnc^2
                var.add <- var.add + rmat[i, m, sc] * diag(V)
                tailv.below <- pnorm(0, mean = mnc, sd = sqrt(diag(V)))
                tailv.below[sqrt(diag(V)) == 0] <- .5
                lfsrmat.add <- lfsrmat.add + rmat[i, m, sc] * cbind(tailv.below, 1 - tailv.below)
              }
              if(meth == "post.mn.fac"){
                facmn <- t(loadings) %*% mnc
                facvar <- t(loadings) %*% V %*% loadings
                mn.add <- mn.add + rmat[i, m, sc] * facmn
                moment2.add <- moment2.add + rmat[i, m, sc] * facmn^2
                var.add <- var.add + rmat[i, m, sc] * diag(facvar)
                tailv.below <- pnorm(0, mean = facmn, sd = sqrt(diag(facvar)))
                tailv.below[sqrt(diag(facvar)) == 0] <- .5
                lfsrmat.add <- lfsrmat.add + rmat[i, m, sc] * cbind(tailv.below, 1 - tailv.below)
              }
              if(meth  == "update.Sigl"){
                Ctl[[sc]] <- Ctl[[sc]] + rmat[i, m, sc] * (V + mnc %*% t(mnc)) / omegaseq[m]
                btl[[sc]] <- btl[[sc]] + rmat[i, m, sc]
              }
            }
          }
        }
        if(meth %in% c("post.mn", "post.mn.fac")){
          mnmat[i, ] <- mn.add
          sdmat[i, ] <- sqrt(var.add + moment2.add - mn.add^2)
          lfsrmat[i, ] <- pmin(lfsrmat.add[, 1], lfsrmat.add[, 2])
        }
      }
    }
    if(meth == "update.Sigl"){
      Sigoutl <- Sigl
      Ksig.new <- Ksig
      fa.failed <- F
      # if(update.K){
      for(sc in Sig.inds.calc){
        obj.K.v <- c()
        eff.n <- sum(pimat[, sc]) * N
        fa.startc <- NULL
        if(update.K[sc]){
            dim.check.fa <- switch(fac.model, pca = 1:P,
                               fa = unique(pmin(c(1, 2, 3, 5, 10, 15, 20, 25, 30, 35, 40, 50, 100, 200, 300, 500), P)))
        } else {
          dim.check.fa <- Ksig[sc]
        }
        ndim.check <- length(dim.check.fa)
        for(j in 1:(ndim.check + 1)){
          if(j <= ndim.check)
            Kc <- dim.check.fa[j]
          if(fac.model == "fa")
            max.K.change <- 1
          if(j == ndim.check + 1)
            K.new <- Kc <- max(c(1, floor(Ksig[sc] * (1 - max.K.change)), dim.check.fa[which.min(obj.K.v)]))
          # K.new <- Kc <- max(c(1, floor(Ksig[sc] * (1 - max.K.change)), dim.check.fa[which.min(obj.K.v)]))
          # K.new <- Kc <- max(c(1, floor(Ksig[sc] * (1 - max.K.change)), (1:P)[which.min(obj.K.v)] + 1))
          if(Kc <= P){
            if(fac.model == "pca"){
              eigc <- eigen(Ctl[[sc]] / btl[[sc]])
              sig2hat <- ifelse(Kc == P, 0, mean(eigc$values[(Kc + 1):P]))
              Sigoutl[[sc]] <- eigc$vectors[, 1:Kc, drop = F] %*% 
                as.matrix(diag(eigc$values[1:Kc], nrow = Kc, ncol = Kc) - sig2hat * diag(rep(1, Kc), nrow = Kc, ncol = Kc)) %*% 
                t(eigc$vectors[, 1:Kc, drop = F]) + sig2hat * diag(rep(1, P))
            }
            if(fac.model == "fa"){
              # print(Kc)
              fa.covmat <- Ctl[[sc]] / btl[[sc]]
              lower.psi <- .005
              control <- list(nstart = 1, trace = F, lower = lower.psi, opt = list(maxit = 10000, factr = 1e7), rotate = NULL)
              fa.tryer <- try({
                fa.out <- factanal(factors = Kc, covmat = fa.covmat, n.obs = N, start = fa.startc, control = control)
              })
              fa.failed <- inherits(fa.tryer, "try-error")
              if(!fa.failed){
                fa.startc <- as.matrix(pmax(lower.psi, fa.out$uniquenesses))
                fa.corout <- fa.out$loadings %*% t(fa.out$loadings) + diag(fa.out$uniquenesses)
                # image(fa.corout)
                Sigoutl[[sc]] <- fa.corout * (sqrt(diag(fa.covmat)) %o% sqrt(diag(fa.covmat)))
              } else { #Do pca optimization
                eigc <- eigen(Ctl[[sc]] / btl[[sc]])
                # sig2hat <- mean(eigc$values[(Kc + 1):P])
                sig2hat <- ifelse(Kc == P, 0, mean(eigc$values[(Kc + 1):P]))
                Sigoutl[[sc]] <- eigc$vectors[, 1:Kc, drop = F] %*% 
                  as.matrix(diag(eigc$values[1:Kc], nrow = Kc, ncol = Kc) - sig2hat * diag(rep(1, Kc), nrow = Kc, ncol = Kc)) %*% 
                  t(eigc$vectors[, 1:Kc, drop = F]) + sig2hat * diag(rep(1, P))
              }
            }
          } else {
            Sigoutl[[sc]] <- Ctl[[sc]] / btl[[sc]]
          }
          if(j <= ndim.check){
            if(fac.model == "pca" | fa.failed)
              npar <- Kc * P + 1 - .5 * Kc * (Kc - 1)
            if(fac.model == "fa")
              npar <- Kc * P + P - .5 * Kc * (Kc - 1)
            # npar <- 2 * (Kc * P + 1 - .5 * Kc * (Kc - 1))
            # obj.K.v[j] <- 2 * npar / 2 * log(eff.n) + 
            obj.K.v[j] <- btl[[sc]] * determinant(Sigoutl[[sc]], logarithm = T)$modulus + sum(diag(solve(Sigoutl[[sc]], Ctl[[sc]]))) 
            if(include.bic)
              obj.K.v[j] <- obj.K.v[j] + npar * log(eff.n) * bic.pen.mult
          }
        }
        # dimnames(Sigoutl[[sc]]) <- list(ph.use, ph.use)
        Ksig.new[sc] <- K.new
        if(update.K[sc])
          plot(dim.check.fa, obj.K.v, ty = "l")
      }
      Ksig <- Ksig.new
      # }
    # } else {
    #   for(sc in Sig.inds.calc){
    #     Kc <- Ksig[sc]
    #     if(Kc < P){
    #       eigc <- eigen(Ctl[[sc]] / btl[[sc]])
    #       sig2hat <- mean(eigc$values[(Kc + 1):P])
    #       Sigoutl[[sc]] <- eigc$vectors[, 1:Kc, drop = F] %*% as.matrix(diag(eigc$values[1:Kc], nrow = Kc, ncol = Kc) - 
    #               sig2hat * diag(rep(1, Kc), nrow = Kc, ncol = Kc)) %*% t(eigc$vectors[, 1:Kc, drop = F]) + sig2hat * diag(rep(1, P))
    #     } else {
    #       Sigoutl[[sc]] <- Ctl[[sc]] / btl[[sc]]
    #     }
    #     # dimnames(Sigoutl[[sc]]) <- list(ph.use, ph.use)
    #   }
    # }
      return(list(Sigl = Sigoutl, Ksig = Ksig, Ctl = Ctl, btl = btl, rmat = rmat, loglikv = loglikv, llmat = llmat))
    }
    if(meth %in% c("post.mn", "post.mn.fac"))
      return(list(mnmat = mnmat, sdmat = sdmat, lfsrmat = lfsrmat, rmat = rmat, loglikv = loglikv, llmat = llmat))
  }
  if(meth == "just.obj"){
    negobj <- 0
    if(prior.in.obj){
      for(sc in which(update.Sig)){
        eff.n <- N
        Kc <- Ksig[sc]
        if(fac.model == "pca")
          npar <- Kc * P + 1 - .5 * Kc * (Kc - 1)
        if(fac.model == "fa")
          npar <- Kc * P + P - .5 * Kc * (Kc - 1)
        # npar <- Kc * P + 1 - .5 * Kc * (Kc - 1)
        # npar <- 2 * (Kc * P + 1 - .5 * Kc * (Kc - 1))
        # negobj <- negobj - 2 * npar / 2 * log(eff.n)
        if(include.bic)
          negobj <- negobj - npar * log(eff.n) * bic.pen.mult / 2
        if(wish.pri)
          negobj <- negobj - .5 * sum(diag(solve(Sigl[[sc]]))) - .5 * (2 * wish.pri.nu + 1) * determinant(Sigl[[sc]], logarithm = TRUE)$modulus
      }
    }
    llikv <- rep(NA, nrow(Y.em))
    names(llikv) <- rownames(Y.em)
    # cat(paste0("Initializing likelihood matrix\n"))
    t0 <- Sys.time()
    for(i in 1:nrow(Y.em)){
      propdone <- i / nrow(Y.em)
      cat(paste0("\rCalculating likelihood matrix; estimated time remaining = ", 
                 ceiling(difftime(Sys.time(), t0, units = "mins") * (1 - propdone) / propdone), " mins "))
      if(any(!is.na(Y.em[i, ]))){
        if(recalc.llmat){
          # logmvndensc <- matrix(-Inf, M, nSig)
          for(sc in Sig.inds.calc){
            for(m in 1:M){
              if(pimat[m, sc] > 0){
                llmat[i, m, sc] <- calc.log.mvnorm.dens(yg = Y.em[i, ], sg = S.em[i, ], Sig = Sigl[[sc]] * omegaseq[m], R = R)
              }
            }
          }
        }
        logmvndensc <- llmat[i, , ]
        maxlogdens <- max(logmvndensc, na.rm = T)
        llikv[i] <- log(sum(pimat * exp(logmvndensc - maxlogdens), na.rm = T)) + maxlogdens
        negobj <- negobj + llikv[i]
      }
    }
    cat("\n")
    return(list(obj = c(-negobj), llikv = llikv, llmat = llmat))
  }
}















# round(Ctl[[1]][1:10, 1:10], 2)







# var.explained.init <- cumsum(eigen(Sig.init)$values) / sum(eigen(Sig.init)$values)
# tail.area <- .05
# percent.start.tail <- .8
# tail.start <- match(T, percent.start.tail < var.explained.init)
# geom.p <- .5#1 - exp(log(tail.area) / tail.start)
# 
# for(Kc in seq(10, 140, by = 10)){#                Kc = 6
#   print(Kc)
#   fa.covmat <- Ctl[[sc]] / btl[[sc]]
#   lower.psi <- .005
#   control <- list(nstart = 1, trace = TRUE, lower = lower.psi, opt = list(maxit = 10000, factr = 1e8), rotate = NULL)
#   fa.out <- factanal(factors = Kc, covmat = fa.covmat, n.obs = N, start = startc, control = control)
#   startc <- as.matrix(pmax(lower.psi, fa.out$uniquenesses))
# }
# # , start = diag(Ctl[[sc]] / btl[[sc]]))
# Kc <- floor(P / 32)
# fa.out <- factanal(factors = Kc, covmat = Ctl[[sc]] / btl[[sc]])
# str(test)
# Sigoutl[[sc]] <- fa.out$loadings %*% t(fa.out$loadings) + diag(fa.out$uniquenesses)
# eigc <- eigen(Ctl[[sc]] / btl[[sc]])
# sig2hat <- mean(eigc$values[(Kc + 1):P])
# Sigoutl[[sc]] <- eigc$vectors[, 1:Kc, drop = F] %*% 
#   as.matrix(diag(eigc$values[1:Kc], nrow = Kc, ncol = Kc) - sig2hat * diag(rep(1, Kc), nrow = Kc, ncol = Kc)) %*% 
#   t(eigc$vectors[, 1:Kc, drop = F]) + sig2hat * diag(rep(1, P))
