
#####################################################################
# Calculate KL divergence between split models and combined model for Factor Sensitivity Analysis
kl1 <- kl2 <- c()
Sig.comb <- resl.comp[[control$mv_meth_nam_use]]$Sig.comb
P <- control$default_parameters$impc$P
names(objl)
subsam_meth_name <- c("impc_MVphen_N_500_nSig_1_K_20", "impc_MVphen_rand_nSig_1_K_20")[1]
for(seed in 1:control$n_subsamples_benchmark){
  Sig_for_curr_data_split <- objl[[subsam_meth_name]][[seed]]$Sigl[[1]]
  if(!is.null(Sig_for_curr_data_split)){
    kl1[seed] <- .5 * (sum(diag(solve(Sig.comb) %*% Sig_for_curr_data_split)) - P + 
                         determinant(Sig.comb, logarithm = T)$modulus - determinant(Sig_for_curr_data_split, logarithm = T)$modulus)
    kl2[seed] <- .5 * (sum(diag(solve(Sig_for_curr_data_split) %*% Sig.comb)) - P + 
                         determinant(Sig_for_curr_data_split, logarithm = T)$modulus - determinant(Sig.comb, logarithm = T)$modulus)
  }
}
seed.largest.kl <- which.max(kl1 + kl2)
resl.comp.subset <- list(mn = compl[[subsam_meth_name]]$mnarr[, , seed.largest.kl], 
                         sd = compl[[subsam_meth_name]]$sdarr[, , seed.largest.kl], 
                         lfsr = compl[[subsam_meth_name]]$lfsrarr[, , seed.largest.kl])

resl_into_err_rate_fn <- list(uv = resl.comp$uv, mv_subsam = resl.comp.subset)
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
all_lines_in_subsam_test <- rownames(resl.comp.subset$mn)[rowSums(is.na(resl.comp.subset$mn)) == 0]
ko_lines_in_subsam_test <- all_lines_in_subsam_test[
  Data_all$impc$linemap[match(all_lines_in_subsam_test, Data_all$impc$linemap$geno), "line.type"] == "trueMut"]
tests_compare <- resimp_subsam[resimp_subsam$geno %in% ko_lines_in_subsam_test, "ph_geno"]
resimp_full <- readRDS(file = control$file.resimp)
tab_curr <- table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam.perm.signsig"],
                  resimp_full[match(tests_compare, resimp_full$ph_geno), paste0(control$mv_meth_nam_use, ".perm.signsig")])


tab.uv.ebi.comp[] <- prettyNum(tab.uv.ebi.comp, big.mark = ",")
tabout.uv <- print(xtable(tab.uv.ebi.comp, label = "tab:uv_ebi_comparison", 
                          caption = "Comparison of signed phenotype hits between our UV model (left) and the existing phenotype calls in the IMPC database (top)"),
                   caption.placement = "top")
cat(tabout.uv, file = paste(dir.save, "/uv_ebi_comp_tab.txt", sep = ""))

tab_curr
# fdr.est.tab(tab_curr)



rin <- readRDS(file = file_preproc)
resimp <- readRDS(file = paste0(control$global_res_dir, "/resimp_comb.RDS"))  
resimp <- resimp[!grepl("fac", resimp$ph), ]
mv_signsig_name <- paste0(control$mv_meth_nam_use, ".perm.signsig")
resimp[which(is.na(resimp[, mv_signsig_name]))[1:10], ]
resimp[, paste0("ebi.", c("mn", "se", "p"))] <- rin[match(resimp$testid, rin$testid), 
                                                    c("genotype_effect_parameter_estimate", "genotype_effect_stderr_estimate", "p_value")]
resimp$ebi.t <- resimp$ebi.mn / resimp$ebi.se
resimp$ebi.signsig <- sign(resimp$ebi.mn) * (resimp$ebi.p < 1e-4)
tab.eb.ebi.comp <- table(resimp[, mv_signsig_name], resimp$ebi.signsig)
tab.uv.ebi.comp <- table(resimp$uv.perm.signsig, resimp$ebi.signsig)
control$mv_meth_nam_use
n.eb.ebi.disagree <- tab.eb.ebi.comp["-1", "1"] + tab.eb.ebi.comp["1", "-1"]
prop.eb.ebi.disagree <- n.eb.ebi.disagree / (n.eb.ebi.disagree + tab.eb.ebi.comp["-1", "-1"] + tab.eb.ebi.comp["1", "1"])
ebi.hit.rate <- mean(resimp$ebi.signsig != 0, na.rm = T)
eb.hit.rate <- mean(resimp[, mv_signsig_name] != 0, na.rm = T)
dir.save <- control$dropbox_text_numbers_dir
save.prop <- c("ebi.hit.rate", "prop.eb.ebi.disagree")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
save.num <- c("n.eb.ebi.disagree")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)

library(xtable)
tab.uv.ebi.comp[] <- prettyNum(tab.uv.ebi.comp, big.mark = ",")
tabout.uv <- print(xtable(tab.uv.ebi.comp, label = "tab:uv_ebi_comparison", 
                          caption = "Comparison of signed phenotype hits between our UV model (left) and the existing phenotype calls in the IMPC database (top)"),
                   caption.placement = "top")
cat(tabout.uv, file = paste(dir.save, "/uv_ebi_comp_tab.txt", sep = ""))
tab.eb.ebi.comp[] <- prettyNum(tab.eb.ebi.comp, big.mark = ",")
tabout.mv <- print(xtable(tab.eb.ebi.comp, label = "tab:eb_ebi_comparison", 
                          caption = "Comparison of signed phenotype hits between our MV model (left) and the existing phenotype calls in the IMPC database (top)"),
                   caption.placement = "top")
cat(tabout.mv, file = paste(dir.save, "/eb_ebi_comp_tab.txt", sep = ""))
