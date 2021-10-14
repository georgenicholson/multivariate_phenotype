#Function for estimating Fsr from a concordance table
fdr.est.tab <- function(tab){
  require(Hmisc)
  namv <- seq(-1, 1, by = 1)
  tab <- tab[match(namv, rownames(tab)), match(namv, colnames(tab))]
  dimnames(tab) <- list(namv, namv)
  tab[is.na(tab)] <- 0
  x <- tab["1", "-1"] + tab["-1", "1"] 
  n <- tab["1", "-1"] + tab["-1", "1"] + tab["-1", "-1"] + tab["1", "1"]
  ci <- binconf(x = x, n = n)
  ci[ci > 2 / 3] <- 2 / 3 # 2/3 is maximum q compatible with model
  fsr.est <- .5 * (1 - sqrt(1 - 2 * ci))
  return(list(est = fsr.est[1], ci = fsr.est[2:3], n = n, x = x))
}

