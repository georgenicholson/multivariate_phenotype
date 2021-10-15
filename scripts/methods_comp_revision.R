rm(list = ls())
# data <- "impc"
Data <- "eqtl"
xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
source(paste0(R.file.dir, "/impc_mv_analysis/err_rate_control.R"))
load(file = file.glob.res)
load(file = uv.results.Y.S)
linemap.use <- linemap
linemap.use$line.type <- ifelse(linemap$line.type == "trueMut", "trueMutTes", "negConTes")
linemap.use <- linemap.use[linemap.use$geno %in% rownames(resl.comp$eb$mn), ]
# resl = resl.comp[names(resl.comp) %in% c("uv", "eb", "mash", "mash.no.imp", "varimax")]
# err.rate.meth = "perm"; sep.imp.thresh = F
# test.stat = "z"; 
# use.upper.fp.est = F; control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1];
# err.thresh = .05; centre.specific.thresh = F
str(resl.comp, m = 1)

out.perm <- err.rate.control(resl = resl.comp[names(resl.comp) %in% c("uv", "eb", "mash", "mash.no.imp", "varimax")], 
                             err.rate.meth = "perm", sep.imp.thresh = F,
                             test.stat = "z", linemap = linemap.use, phmap = phmap, Yhat = Yhat,
                             use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                             err.thresh = .05, centre.specific.thresh = F)
out.perm$restab
out.perm.lfsr <- err.rate.control(resl = resl.comp[names(resl.comp) %in% c("eb", "mash", "mash.no.imp")], 
                                  err.rate.meth = "perm", sep.imp.thresh = F,
                                  test.stat = "lfsr", linemap = linemap.use, phmap = phmap, Yhat = Yhat,
                                  use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                                  err.thresh = .05, centre.specific.thresh = F)
out.perm.lfsr$restab
out.lfsr <- err.rate.control(resl = resl.comp[names(resl.comp) %in% c("eb", "mash")], 
                             err.rate.meth = "lfsr", sep.imp.thresh = F,
                             test.stat = "lfsr", linemap = linemap.use, phmap = phmap, Yhat = Yhat,
                             use.upper.fp.est = F, control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                             err.thresh = .05, centre.specific.thresh = F)
out.lfsr$restab


#######################################################
# Extract cross-validated likelihoods and add to table


str(compl, m = 2)


dput(colnames(out.perm$restab))

colkeep <- c("meth", "err.rate.meth",  
             "test.stat", "mvhitimp", "mvhitnonimp", "fdr.est.imp", 
             "fdr.est.nonimp")
tabout <- rbind(out.perm$restab[out.perm$restab$meth %in% c("uv", "eb", "mash"), colkeep],
      out.perm.lfsr$restab[out.perm.lfsr$restab$meth %in% c("uv", "eb", "mash"), colkeep],
      out.lfsr$restab[out.lfsr$restab$meth %in% c("mash", "eb"), colkeep])
for(j in 1:ncol(tabout)){
  if(is.numeric(tabout[, j]))
  tabout[, j] <- formatC(100 * tabout[, j], digits = 1, format = "f")
}
tabout[order(tabout$meth), ]


str(compl)


plot(compl$eb$crossval.llv, compl$mash$crossval.llv)
abline(0, 1)



hist(resl.comp$mash$mn / resl.comp$mash$sd)
hist(resl.comp$mash.no.imp$mn / resl.comp$mash.no.imp$sd)


par(mfrow = c(2, 1))
hist(-log10(resl.comp$mash$lfsr))
hist(-log10(resl.comp$mash.no.imp$lfsr))


rownames(resl.comp$mash$lfsr)[which(rowSums(resl.comp$mash.no.imp$lfsr < 1e-5, na.rm = T) > 0)]
str(out.perm$resimp, m = 1)
sum()



mean(is.na(resl.comp$mash.no.imp$lfsr))
hist(resl.comp$mash$lfsr[is.na(resl.comp$mash.no.imp$lfsr)])
hist(resl.comp$mash$lfsr[!is.na(resl.comp$mash.no.imp$lfsr)])

mn / resl.comp$mash$sd)
hist(resl.comp$mash.no.imp$mn / resl.comp$mash.no.imp$sd)


resl.comp <- list()
resl.comp$uv <- out$resl$uv
resl.comp$uv.ss <- out$resl$uv.ss
resl.comp$eb <- out$resl$eb
resl.comp$eb.ss <- out$resl$eb.ss
if(collect.factors){
  for(fac.meth in fac.methv)
    resl.comp[[fac.meth]] <- out$resl[[fac.meth]]
}
resimp <- out$resimp
resimp$line.type <- ifelse(out$resimp$line.type == "trueMutTes", "trueMut", "negCon")
# resimp <- resimp[, !grepl("mash", colnames(resimp))]
sapply(Ksigl, mean)



c("meth", "err.rate.meth", "centre.specific.thresh", "sep.imp.thresh", 
  "test.stat", "use.upper", "mvhitimp", "mvhitnonimp", "fdr.est.imp", 
  "fdr.est.nonimp", "ref.lines.post", "ref.lines.postcomp", "ref.lines.postimp", 
  "ref.lines.postimpcomp", "ref.lines.post.l", "ref.lines.postcomp.l", 
  "ref.lines.postimp.l", "ref.lines.postimpcomp.l", "ref.lines.post.u", 
  "ref.lines.postcomp.u", "ref.lines.postimp.u", "ref.lines.postimpcomp.u", 
  "ref.lines.post.x", "ref.lines.postcomp.x", "ref.lines.postimp.x", 
  "ref.lines.postimpcomp.x", "ref.lines.post.n", "ref.lines.postcomp.n", 
  "ref.lines.postimp.n", "ref.lines.postimpcomp.n", "het.hom.post", 
  "het.hom.postcomp", "het.hom.postimp", "het.hom.postimpcomp", 
  "het.hom.post.l", "het.hom.postcomp.l", "het.hom.postimp.l", 
  "het.hom.postimpcomp.l", "het.hom.post.u", "het.hom.postcomp.u", 
  "het.hom.postimp.u", "het.hom.postimpcomp.u", "het.hom.post.x", 
  "het.hom.postcomp.x", "het.hom.postimp.x", "het.hom.postimpcomp.x", 
  "het.hom.post.n", "het.hom.postcomp.n", "het.hom.postimp.n", 
  "het.hom.postimpcomp.n")