# rm(list = ls())
data <- "impc"
source(paste0(ifelse(Sys.info()["sysname"] == "Linux", "/mnt/c", "C:"),
              "/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
###########################################################
#Load UV results
load(uv.resmat.file.post.qc)
str(resmat.qc)

genuse <- unique(resmat.qc$geno)
phuse <- unique(resmat.qc$ph)
true.use <- unique(resmat.qc[resmat.qc$line.type == "trueMut", "geno"])
neg.tra.use <- unique(resmat.qc[resmat.qc$line.type == "negConTra", "geno"])
neg.tes.use <- unique(resmat.qc[resmat.qc$line.type == "negConTes", "geno"])
em.geno.all <- c(true.use, neg.tra.use)

################################################
#Create wide-format matrices of data from long format
ng <- length(genuse)
nph <- length(phuse)
Yhat <- smat <- Yhat.genosex <- smat.genosex <- matrix(NA, ng, nph, dimnames = list(genuse, phuse))

# Use main genotype effect from two-way sex-genotype interaction model
resmat.qc$uv.main.geno.eff.mn.sc <- resmat.qc$uv.mn.inter.geno.sc
resmat.qc$uv.main.geno.eff.sd.sc <- resmat.qc$uv.sd.inter.geno.sc
# When interaction model not feasible due to only one sex being available, then just use main effect from additive sex plus genotype model
resmat.qc$uv.main.geno.eff.mn.sc[is.na(resmat.qc$uv.main.geno.eff.mn.sc)] <- 
  resmat.qc$uv.mn.just.geno.sc[is.na(resmat.qc$uv.main.geno.eff.mn.sc)]
resmat.qc$uv.main.geno.eff.sd.sc[is.na(resmat.qc$uv.main.geno.eff.sd.sc)] <- 
  resmat.qc$uv.sd.just.geno.sc[is.na(resmat.qc$uv.main.geno.eff.sd.sc)]

Yhat[cbind(match(resmat.qc$geno, rownames(Yhat)), match(resmat.qc$ph, colnames(Yhat)))] <- resmat.qc$uv.main.geno.eff.mn.sc
smat[cbind(match(resmat.qc$geno, rownames(smat)), match(resmat.qc$ph, colnames(smat)))] <- resmat.qc$uv.main.geno.eff.sd.sc
Yhat.genosex[cbind(match(resmat.qc$geno, rownames(Yhat.genosex)), match(resmat.qc$ph, colnames(Yhat.genosex)))] <- resmat.qc$uv.mn.inter.genosex.sc
smat.genosex[cbind(match(resmat.qc$geno, rownames(smat.genosex)), match(resmat.qc$ph, colnames(smat.genosex)))] <- resmat.qc$uv.sd.inter.genosex.sc


# graphics.off()  
# # plot(resmat.qc$uv.main.geno.eff.mn.sc, resmat.qc$uv.mn.sc)
# plot(resmat.qc$uv.main.geno.eff.mn.sc / resmat.qc$uv.main.geno.eff.sd.sc, resmat.qc$uv.mn.sc / resmat.qc$uv.sd.sc)
# abline(0, 1)
# Yhat[cbind(match(resmat.qc$geno, rownames(Yhat)), match(resmat.qc$ph, colnames(Yhat)))] <- resmat.qc$uv.mn.sc
# smat[cbind(match(resmat.qc$geno, rownames(smat)), match(resmat.qc$ph, colnames(smat)))] <- resmat.qc$uv.sd.sc

# ##########################################
# # Checking new analysis
# t.unadjusted <- resmat.qc$uv.mn.sc / resmat.qc$uv.sd.sc
# t.new.adjusted <- resmat.qc$uv.mn.inter.geno.sc / resmat.qc$uv.sd.inter.geno.sc
# t.new.genosex.adjusted <- resmat.qc$uv.mn.inter.genosex.sc / resmat.qc$uv.sd.inter.genosex.sc
# t.new <- resmat.qc$uv.mn.just.geno.sc / resmat.qc$uv.sd.just.geno.sc
# 
# t.new.adjusted[is.na(t.new.adjusted)] <- t.new[is.na(t.new.adjusted)]
# 
# mean(is.na(t.new.adjusted))
# mean(is.na(t.unadjusted))
# 
# graphics.off()  
# plot(t.new, t.unadjusted)  
# abline(0, 1)
# plot(t.new, t.new.adjusted)  
# abline(0, 1)
# plot(t.new.adjusted, t.new.genosex.adjusted)  
# abline(0, 1)
# th <- 1.5
# table(abs(t.new) > th, abs(t.new.adjusted) > th)
# str(resmat.qc)
# str(d)
# table(resmat.qc$mut.n)
# mean(resmat.qc$mut.n >= 14)
# 
# ##################################

# #####################################################################################
# #Recalibrate UV results, if specified by scale.uv.pre, and calculate z-statistics
# uv.t.temp <- Yhat[genuse, phuse] / smat[genuse, phuse]
# smat.rescaled <- smat[genuse, phuse]
# if(scale.uv.pre){
#   smat.rescaled[] = NA
#   for(cenc in cennamun){
#     cen.all.geno <- unique(resmat.qc$geno[which(resmat.qc$cenlong == cenc)])
#     uv.cen.scale.t <- apply(uv.t.temp[neg.tra.use[neg.tra.use %in% cen.all.geno], ], 2, function(v) sd(v, na.rm = T))
#     smat.rescaled[cen.all.geno, ] <- t(t(smat[cen.all.geno, phuse]) * uv.cen.scale.t)
#   }
# }

# plot(resmat.qc$uv.mn.sc, resmat.qc$uv.mn.just.geno.sc)
# hist(resmat.qc$uv.mn.inter.genosex.sc / resmat.qc$uv.sd.inter.genosex.sc)
# mean(abs(resmat.qc$uv.mn.inter.genosex.sc / resmat.qc$uv.sd.inter.genosex.sc)[resmat.qc$line.type == "trueMut"] > 2, na.rm = T)
# mean(abs(resmat.qc$uv.mn.sc / resmat.qc$uv.sd.sc)[resmat.qc$line.type != "trueMut"] > 2, na.rm = T)
# str(resmat.qc)
# npl <- 100000
# graphics.off()
# plot(resmat.qc$uv.mn.inter.genosex.sc[1:npl], resmat.qc$uv.mn.sc[1:npl])


linemap <- unique(resmat.qc[, c("geno", "line.type", "cen")])
phmap <- pout[match(colnames(Yhat), pout$ph), ]
Y <- Yhat
S <- smat
Y[is.na(Yhat)] = 0
S[is.na(Yhat)] = prior.sd.on.unobserved.thetahat
Y.genosex <- Yhat.genosex
S.genosex <- smat.genosex
Y.genosex[is.na(Yhat.genosex)] = 0
S.genosex[is.na(Yhat.genosex)] = prior.sd.on.unobserved.thetahat
plot((Y.genosex / S.genosex)[1:10000], (Y / S)[1:10000])


##################################
#Shorten a few phen names
phmap$nam <- gsub("Forelimb grip strength measurement mean", "Forelimb grip strength", phmap$nam)
phmap$nam <- gsub("Forelimb and hindlimb grip strength measurement mean", "Forelimb and hindlimb grip strength", phmap$nam)
phmap$nam <- gsub("Forelimb grip strength normalised against body weight", "Forelimb grip strength (BWnorm)", phmap$nam)
phmap$nam <- gsub("Forelimb and hindlimb grip strength normalised against body weight", "Forelimb and hindlimb grip strength (BWnorm)", phmap$nam)
phmap$nam <- gsub("Bone Mineral Density \\(excluding skull\\)", "Bone Mineral Density", phmap$nam)
phmap$nam <- gsub("Bone Mineral Content \\(excluding skull\\)", "Bone Mineral Content", phmap$nam)

save(Yhat, smat, Y, S, Yhat.genosex, smat.genosex, Y.genosex, S.genosex, linemap, phmap, file = uv.results.Y.S, version = 2)
























