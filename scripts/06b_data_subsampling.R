
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
tab_curr<-tab_curr_num <- table(resimp_subsam[match(tests_compare, resimp_subsam$ph_geno), "mv_subsam.perm.signsig"],
                  resimp_full[match(tests_compare, resimp_full$ph_geno), paste0(control$mv_meth_nam_use, ".perm.signsig")])
subsam_n_disagree <- (tab_curr_num["-1", "1"] + tab_curr_num["1", "-1"])
subsam_n_both_call <- tab_curr_num["-1", "1"] + tab_curr_num["1", "-1"] + tab_curr_num["1", "1"] + tab_curr_num["-1", "-1"]

save.num <- c("subsam_n_disagree", "subsam_n_both_call")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)

tab_curr[] <- prettyNum(tab_curr, big.mark = ",")
tabout_subsam <- print(xtable(tab_curr, label = "tab:subsampling_comp", 
                          caption = "Comparison of signed phenotype hits between our UV model (left) and the existing phenotype calls in the IMPC database (top)"),
                   caption.placement = "top",
                   floating = FALSE)
cat(tabout_subsam, file = paste(control$dropbox_table_dir, "/subsampling_comp.txt", sep = ""))


