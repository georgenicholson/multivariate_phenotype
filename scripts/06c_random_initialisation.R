
#####################################################################
# Calculate KL divergence between split models and combined model for Factor Sensitivity Analysis
kl1 <- kl2 <- c()
Sig.comb <- resl.comp[[control$mv_meth_nam_use]]$Sig.comb
P <- control$default_parameters$impc$P
names(objl)
subsam_meth_name_rand <- "impc_MVphen_rand_nSig_1_K_20"
subsam_meth_name_fixed <- "impc_MVphen_nSig_1_K_20"
symm_KL <- function(Sig1, Sig2) {
  kl1 <- .5 * (sum(diag(solve(Sig1) %*% Sig2)) - P + 
                       determinant(Sig1, logarithm = T)$modulus - determinant(Sig2, logarithm = T)$modulus)
  kl2 <- .5 * (sum(diag(solve(Sig2) %*% Sig1)) - P + 
                       determinant(Sig2, logarithm = T)$modulus - determinant(Sig1, logarithm = T)$modulus)
  return(kl1 + kl2)
}
KL_split_rand_fixed <- KL_split_rand_comb <- KL_split_fixed_comb <- c()
for(seed in 1:control$n_subsamples_benchmark){
  Sig_for_curr_data_split_rand <- objl[[subsam_meth_name_rand]][[seed]]$Sigl[[1]]
  Sig_for_curr_data_split_fixed <- objl[[subsam_meth_name_fixed]][[seed]]$Sigl[[1]]
  KL_split_rand_fixed[seed] <- symm_KL(Sig_for_curr_data_split_rand, Sig_for_curr_data_split_fixed)
  KL_split_fixed_comb[seed] <- symm_KL(Sig_for_curr_data_split_fixed, Sig.comb)
  KL_split_rand_comb[seed] <- symm_KL(Sig_for_curr_data_split_rand, Sig.comb)
  
  # if(!is.null(Sig_for_curr_data_split)){
  #   kl1[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split_rand) %*% Sig_for_curr_data_split)) - P + 
  #                        determinant(Sig_for_curr_data_split_rand, logarithm = T)$modulus - determinant(Sig_for_curr_data_split, logarithm = T)$modulus)
  #   kl2[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split) %*% Sig_for_curr_data_split_rand)) - P + 
  #                        determinant(Sig_for_curr_data_split, logarithm = T)$modulus - determinant(Sig_for_curr_data_split_rand, logarithm = T)$modulus)
  # }
}

str(all_sig)
all_sig <- lapply(c(objl[[subsam_meth_name_rand]][1:10], objl[[subsam_meth_name_fixed]][1:10]), function(x) x$Sigl[[1]])
dmat <- matrix(NA, 20, 20)
for (i in 1:20) {
  for(j in 1:20) {
    dmat[i, j] <- symm_KL(all_sig[[i]], all_sig[[j]])
  }
}
image(dmat)

# Cross lik validation
#OR images of cov matrices

test <- outer(all_sig, all_sig, FUN = function(x, y) symm_KL(x, y))


plot(KL_split_fixed_comb, KL_split_rand_comb)


seed.largest.kl <- which.max(kl1 + kl2)
resl_rand <- list(mn = compl[[subsam_meth_name_rand]]$mnarr[, , seed.largest.kl], 
                         sd = compl[[subsam_meth_name_rand]]$sdarr[, , seed.largest.kl], 
                         lfsr = compl[[subsam_meth_name_rand]]$lfsrarr[, , seed.largest.kl])
resl_fixed <- list(mn = compl[[subsam_meth_name_fixed]]$mnarr[, , seed.largest.kl], 
                  sd = compl[[subsam_meth_name_fixed]]$sdarr[, , seed.largest.kl], 
                  lfsr = compl[[subsam_meth_name_fixed]]$lfsrarr[, , seed.largest.kl])

resl_into_err_rate_fn <- list(uv = resl.comp$uv, mv_subsam_rand = resl_rand, mv_subsam_fixed = resl_fixed)
out.perm <- err.rate.control(control = control, 
                             resl = resl_into_err_rate_fn,
                             err.rate.meth = "perm", 
                             test.stat = "z", 
                             linemap = Data_all$impc$linemap, 
                             reflinemap = Data_all$impc$reflinemap, 
                             phmap = Data_all$impc$phmap, 
                             cenmap = Data_all$impc$cenmap,
                             Yhat = Data_all$impc$Y_raw,
                             control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
                             err.thresh = .05, 
                             p.complete.null.true = 1, 
                             p.test.null.true = 1)
resimp_subsam <- out.perm$resimp[!grepl("fac", out.perm$resimp$ph), ]

all_lines_in_subsam_test <- rownames(resl_rand$mn)[rowSums(is.na(resl_rand$mn)) == 0]
ko_lines_in_subsam_test <- all_lines_in_subsam_test[
  Data_all$impc$linemap[match(all_lines_in_subsam_test, Data_all$impc$linemap$geno), "line.type"] == "trueMut"]
tests_compare <- resimp_subsam[resimp_subsam$geno %in% ko_lines_in_subsam_test, "ph_geno"]
tab_rand_init <- table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_rand.perm.signsig"],
                       resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_fixed.perm.signsig"])

table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_rand.perm.signsig"])
table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_fixed.perm.signsig"])


# resimp_full <- readRDS(file = control$file.resimp)
# tab_rand_init <- table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam.perm.signsig"],
#                   resimp_full[match(tests_compare, resimp_full$ph_geno), paste0(control$mv_meth_nam_use, ".perm.signsig")])


tab_rand_init[] <- prettyNum(tab_rand_init, big.mark = ",")
tab_out_rand_init <- print(xtable::xtable(tab_rand_init, label = "tab:rand_init", 
                          caption = ""), floating = FALSE, caption.placement = 'top')
cat(tab_out_rand_init, file = paste(control$dropbox_table_dir, "/tab_rand_init.txt", sep = ""))




# 
# 
# par(mfrow = c(1, 2))
# Cor_for_curr_data_split_rand <- t(Sig_for_curr_data_split_rand / sqrt(diag(Sig_for_curr_data_split_rand)))  / sqrt(diag(Sig_for_curr_data_split_rand))
# Cor_for_curr_data_split_fixed <- t(Sig_for_curr_data_split_fixed / sqrt(diag(Sig_for_curr_data_split_fixed)))  / sqrt(diag(Sig_for_curr_data_split_fixed))
# image(Cor_for_curr_data_split_rand)
# image(Cor_for_curr_data_split_fixed)
# 
# # fsr est not valid here (no conditional indep between outputs as same data)
# # fdr.est.tab(tab_curr) 
# 
# 
# 
# #####################################################################
# # Calculate KL divergence between split models and combined model for Factor Sensitivity Analysis
# kl1 <- kl2 <- c()
# Sig.comb <- resl.comp[[control$mv_meth_nam_use]]$Sig.comb
# P <- control$default_parameters$impc$P
# names(objl)
# subsam_meth_name_rand <- "impc_MVphen_rand_nSig_1_K_20"
# subsam_meth_name_fixed <- "impc_MVphen_N_500_nSig_1_K_20"
# for(seed in 1:control$n_subsamples_benchmark){
#   Sig_for_curr_data_split_rand <- objl[[subsam_meth_name_rand]][[seed]]$Sigl[[1]]
#   Sig_for_curr_data_split_fixed <- objl[[subsam_meth_name_fixed]][[seed]]$Sigl[[1]]
#   if(!is.null(Sig_for_curr_data_split)){
#     kl1[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split_rand) %*% Sig_for_curr_data_split)) - P + 
#                          determinant(Sig_for_curr_data_split_rand, logarithm = T)$modulus - determinant(Sig_for_curr_data_split, logarithm = T)$modulus)
#     kl2[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split) %*% Sig_for_curr_data_split_rand)) - P + 
#                          determinant(Sig_for_curr_data_split, logarithm = T)$modulus - determinant(Sig_for_curr_data_split_rand, logarithm = T)$modulus)
#   }
# }
# seed.largest.kl <- which.max(kl1 + kl2)
# resl_rand <- list(mn = compl[[subsam_meth_name_rand]]$mnarr[, , seed.largest.kl], 
#                   sd = compl[[subsam_meth_name_rand]]$sdarr[, , seed.largest.kl], 
#                   lfsr = compl[[subsam_meth_name_rand]]$lfsrarr[, , seed.largest.kl])
# resl_fixed <- list(mn = compl[[subsam_meth_name_fixed]]$mnarr[, , seed.largest.kl], 
#                    sd = compl[[subsam_meth_name_fixed]]$sdarr[, , seed.largest.kl], 
#                    lfsr = compl[[subsam_meth_name_fixed]]$lfsrarr[, , seed.largest.kl])
# 
# resl_into_err_rate_fn <- list(uv = resl.comp$uv, mv_subsam_rand = resl_rand, mv_subsam_fixed = resl_fixed)
# out.perm <- err.rate.control(control = control, 
#                              resl = resl_into_err_rate_fn,
#                              err.rate.meth = "perm", 
#                              test.stat = "z", 
#                              linemap = Data_all$impc$linemap, 
#                              reflinemap = Data_all$impc$reflinemap, 
#                              phmap = Data_all$impc$phmap, 
#                              cenmap = Data_all$impc$cenmap,
#                              Yhat = Data_all$impc$Y_raw,
#                              control.level = c("line.fdr", "line.fwer", "phcen.fdr")[1],
#                              err.thresh = .05, 
#                              p.complete.null.true = 1, 
#                              p.test.null.true = 1)
# resimp_subsam <- out.perm$resimp[!grepl("fac", out.perm$resimp$ph), ]
# 
# all_lines_in_subsam_test <- rownames(resl_rand$mn)[rowSums(is.na(resl_rand$mn)) == 0]
# ko_lines_in_subsam_test <- all_lines_in_subsam_test[
#   Data_all$impc$linemap[match(all_lines_in_subsam_test, Data_all$impc$linemap$geno), "line.type"] == "trueMut"]
# tests_compare <- resimp_subsam[resimp_subsam$geno %in% ko_lines_in_subsam_test, "ph_geno"]
# tab_rand_init <- table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_rand.perm.signsig"],
#                        resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_fixed.perm.signsig"])
# 
# table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_rand.perm.signsig"])
# table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam_fixed.perm.signsig"])
