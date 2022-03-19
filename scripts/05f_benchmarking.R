

linemap <- Data_all$impc$linemap
methnam_code <- c("MVphen", "ComposeMV")[1]
methnam_out <- c("MVphen", "ComposeMV")[2]
truemuts <- linemap$geno[linemap$line.type == control$nam.truemut]
negcons <- linemap$geno[linemap$line.type == control$nam.negcon]
llmean.splitmat.true <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[truemuts, 1:10], na.rm = T))
llmean.splitmat.null <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[negcons, 1:10], na.rm = T))
true.mns <- colMeans(llmean.splitmat.true, na.rm = T)
true.ses <- apply(llmean.splitmat.true, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))
null.mns <- colMeans(llmean.splitmat.null, na.rm = T)
null.ses <- apply(llmean.splitmat.null, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))

err.data.type <- "comb"
formatfn <- function(v)  formatC(v * 100, digits = 1, format = "f")
format.cvlik.fn <- function(v) formatC(v, digits = 1, format = "f")
caption <- switch(err.data.type, comb = "IMPC data: hit rate and error rate comparison across methods", 
                  single = "IMPC data: hit rate and error rate comparison across methods (single CV split)")
colkeep1 <- c("meth", "err.rate.meth", "test.stat", 
              "mvhitimp", "hit.rate.imp.ci.l", "hit.rate.imp.ci.u", 
              "mvhitnonimp", "hit.rate.nonimp.ci.l", "hit.rate.nonimp.ci.u",
              "line.fdr.est", "line.fdr.ci.l", "line.fdr.ci.u", 
              "fdr.est", "fdr.ci.l", "fdr.ci.u", "fdr.est.imp", "fdr.est.nonimp", 
              "ref.lines.post", "ref.lines.post.l", "ref.lines.post.u")
restaball <- rbind(restabl$perm[, colkeep1], restabl$perm.lfsr[, colkeep1], restabl$lfsr[, colkeep1])
restaball <- restaball[!(grepl("N", restaball$meth) | grepl("rand", restaball$meth)), ]
restaball <- restaball[restaball$meth != "varimax", ]
names(restaball)
restaball$hit.imp.ci <- paste0(formatfn(restaball$mvhitimp), " (", formatfn(restaball$hit.rate.imp.ci.l), ", ", formatfn(restaball$hit.rate.imp.ci.u), ")")
restaball$hit.nonimp.ci <- paste0(formatfn(restaball$mvhitnonimp), " (", formatfn(restaball$hit.rate.nonimp.ci.l), ", ", formatfn(restaball$hit.rate.nonimp.ci.u), ")")
restaball$line.fdr.ci <- paste0(formatfn(restaball$line.fdr.est), " (", formatfn(restaball$line.fdr.ci.l), ", ", formatfn(restaball$line.fdr.ci.u), ")")
restaball$fdr.ci <- paste0(formatfn(restaball$fdr.est), " (", formatfn(restaball$fdr.ci.l), ", ", formatfn(restaball$fdr.ci.u), ")")
restaball$fsr.ci <- paste0(formatfn(restaball$ref.lines.post), " (", formatfn(restaball$ref.lines.post.l), ", ", formatfn(restaball$ref.lines.post.u), ")")
restaball$Method <- NA
restaball$facmod <- "FA"
restaball$Method[grepl(methnam_code, restaball$meth)] <- methnam_out
restaball$Method[grepl("mash", restaball$meth)] <- "mash"
restaball$Method[grepl("XD", restaball$meth)] <- "XD"
restaball$Method[grepl("uv", restaball$meth)] <- "UV"
restaball$mvhitimp[grepl("uv", restaball$meth)] <- NA
restaball.out <- restaball

rename_tests <- c('1A' = '(a)', '1B' = '(b)', '2B' = '(c)')
restaball$Test <- rename_tests[paste0(ifelse(restaball$err.rate.meth == "perm", 1, 2), ifelse(restaball$test.stat == "z", "A", "B"))]

restaball$S <- sapply(strsplit(restaball$meth, spl = "_nSig_"), function(x) substr(x[2], 1, 1))
restaball$S[restaball$Method == "mash"] <- control$default_parameters$impc$P + 10
restaball$S[restaball$Method == "UV"] <- ""
restaball$K <- sapply(strsplit(restaball$meth, spl = "_K_"), function(x) substr(x[2], 1, 2))
restaball$fm <- "fa"
restaball <- restaball[which(restaball$fm != "pca" | is.na(restaball$fm)), ]
restaball$Control <- ifelse(restaball$err.rate.meth == "perm", "$\\mathrm{Fdr}_{\\mathrm{complete}}\\leq 5\\%$", "$\\text{lfsr} \\leq 5\\%$")
restaball$Statistic <- ifelse(restaball$test.stat == "z", "$z$", "$\\text{lfsr}$")
restaball$cvlik <- true.mns[restaball$meth]
restaball$cvlik.l <- true.mns[restaball$meth] - 2 * true.ses[restaball$meth]
restaball$cvlik.u <- true.mns[restaball$meth] + 2 * true.ses[restaball$meth]
restaball$cvlik.ci <- paste0(format.cvlik.fn(restaball$cvlik), " (", format.cvlik.fn(restaball$cvlik.l), ", ", format.cvlik.fn(restaball$cvlik.u), ")")
restaball <- restaball[order(match(restaball$Method, c("UV", "XD", "mash", methnam_code))), ]
rates.format <- c("mvhitimp", "mvhitnonimp", "line.fdr.est", "fdr.est", "fdr.est.imp", "fdr.est.nonimp")
for(j in rates.format){
  if(is.numeric(restaball[, j]))
    restaball[, j] <- formatC(restaball[, j] * 100, digits = 1, format = "f")
}
tail(restaball)

pdf(file = paste0(control$dropbox_figure_dir, "/rand_init_cv_lik.pdf"), 6, 6)
par(oma = c(1, 1, 1, 1), mar = c(4, 4, 4, 4))
plot(llmean.splitmat.true[, "impc_MVphen_nSig_1_K_20"], llmean.splitmat.true[, "impc_MVphen_rand_nSig_1_K_20"],
     main = "Cross-validated log likelihood comparison",
     xlab = "Sample-covariance initialised",
     ylab = "Randomly initialised")
abline(0, 1)
dev.off()

#############################################################################
#Output power and type I error table
library(xtable)
colkeepmap <- data.frame(keep = c("meth", "err.rate.meth", "test.stat", "hit.imp.ci", "hit.nonimp.ci", "line.fdr.est", "line.fdr.ci",
                                  "fdr.est", "fdr.ci", "ref.lines.post", "fsr.ci", "fdr.est.imp", "fdr.est.nonimp", "Method", "Test", 
                                  "S", "K", "Control", "Statistic", "cvlik.ci"),
                         keep.as = c("meth", "err.rate.meth", "test.stat", "missing", "measured", 
                                     "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{complete}}$", "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{complete}}$",
                                     "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{single}}$", "$\\widehat{\\mathrm{Fdr}}_{\\mathrm{single}}$",
                                     "$\\widehat{\\mathrm{Fsr}}_{\\mathrm{replicate}}$", "$\\widehat{\\mathrm{Fsr}}_{\\mathrm{replicate}}$",
                                     "fdr.est.imp", "fdr.est.nonimp", "Method", 
                                     "Test", "S", "K", "Control", "Statistic", "CV Log Likelihood"))
colkeep <- c("Test", "Control", "Statistic", "Method", "S", "K", "hit.nonimp.ci", "hit.imp.ci", "line.fdr.ci",
             "fdr.ci", "fsr.ci")
restaball.out <- restaball[, colkeep]
colnames(restaball.out) <- colkeepmap[match(colkeep, colkeepmap$keep), "keep.as"]
colnames(restaball)


tabout <- print(xtable(restaball.out, label = "tab:hitrates", align = rep("r", ncol(restaball.out) + 1),
                       caption = "Hit rates and error rate comparison across methods"),
                caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
cat(tabout, file = paste(control$dropbox_text_numbers_dir, "/hitrates_combined.txt", sep = ""))

hitList <- split(restaball.out[, 4:ncol(restaball.out)], f = restaball.out$Test)
hitList[[1]][hitList[[1]]$Method == "UV", "missing"] <- ""


which_row_main_meth <- which(hitList[[1]]$Method == methnam_out & hitList[[1]]$S == 1 & hitList[[1]]$K == 20)
which_row_best_measured <- which.max(sapply(strsplit(hitList[[1]]$measured, split = " "), function(x) x[1]))
which_row_best_missing <- which.max(sapply(strsplit(hitList[[1]]$missing, split = " "), function(x) x[1]))
hitList[[1]][which_row_main_meth, ] <- paste0("\\textbf{", hitList[[1]][which_row_main_meth, ], "}")
hitList[[1]][which_row_best_measured, "measured"] <- paste0("\\underline{", hitList[[1]][which_row_best_measured, "measured"], "}")
hitList[[1]][which_row_best_missing, "missing"] <- paste0("\\underline{", hitList[[1]][which_row_best_missing, "missing"], "}")



attr(hitList, "subheadings") <- paste0(names(hitList), " Controlling ", 
                                       restaball.out$Control[match(names(hitList), restaball.out$Test)],
                                       " using ", restaball.out$Statistic[match(names(hitList), restaball.out$Test)],  " statistic ")
xList <- xtableList(hitList, label = paste0("tab:hitrates.", err.data.type), 
                    align = gsub(" llllll", "llllll|", paste0(" ", paste0(rep("l", ncol(hitList[[1]]) + 1), collapse = ""))),
                    caption = caption)

tabout.test <- print(xList, 
                     table.placement = "hp", 
                     colnames.format = "multiple", 
                     add.to.row = addtorow,
                     caption.placement = "top", 
                     sanitize.text.function = function(x){x}, 
                     include.rownames = F, 
                     floating = FALSE)
tabout.test.with.header <- gsub("\n\\\\hline\nMe", 
                                paste0("\n\\\\hline\n& & &\\\\multicolumn\\{2\\}\\{c\\}\\{Hit rate  in \\\\% when data are\\}&",
                                       "\\\\multicolumn\\{3\\}\\{|c\\}\\{Estimated error rate in \\\\% (95\\\\% CI) \\}\\\\\\\\\n",
                                       # "\\\\cline\\{3\\-4\\}\\\\cline\\{5\\-7\\}\\\nMe"), tabout.test)
                                       "\\\nMe"), tabout.test)
cat(tabout.test.with.header, file = paste(control$dropbox_text_numbers_dir, "/hitrates_", err.data.type, ".txt", sep = ""))


colkeep.lik <- c("Method", "S", "K", "cvlik.ci")#, "hit.nonimp.ci", "hit.imp.ci")
restaball.lik <- unique(restaball[, colkeep.lik])
restaball.lik <- restaball.lik[restaball.lik$Method != "UV", ]
colnames(restaball.lik) <- colkeepmap[match(colkeep.lik, colkeepmap$keep), "keep.as"]
table_for_order <- restaball.lik
restaball.lik
which_row_main_meth <- which(restaball.lik$Method == methnam_out & restaball.lik$S == 1 & restaball.lik$K == 20)
which_row_best_cvlik <- which.max(sapply(strsplit(restaball.lik$`CV Log Likelihood`, split = " "), function(x) x[1]))
restaball.lik[which_row_best_cvlik, "CV Log Likelihood"] <- paste0("\\underline{", restaball.lik[which_row_best_cvlik, "CV Log Likelihood"], "}")
restaball.lik[which_row_main_meth, ] <- paste0("\\textbf{", restaball.lik[which_row_main_meth, ], "}")

tabout.lik <- print(xtable(restaball.lik, label = "tab:cvlik_impc", align = rep("l", ncol(restaball.lik) + 1),
                            caption = "Comparison of cross-validated log likelihood across MV methods"),
                            caption.placement = "top", 
                            sanitize.text.function = function(x){x},
                            include.rownames = F, 
                            floating = FALSE)
cat(tabout.lik, file = paste(control$dropbox_text_numbers_dir, "/cvlik_table.txt", sep = ""))



# table_for_order$Method[table_for_order$Method == methnam_out] <- methnam_code
# nam_order <- paste0("eqtl_", table_for_order$Method, "_nSig_", table_for_order$S)
# restaball.lik
# dput(names(true.mns))


res_to_fill <- restaball.lik
Data <- "eqtl"
model_ordered <- c("eqtl_XD_nSig_1", "eqtl_XD_nSig_2", "eqtl_mash_nSig_1", "eqtl_MVphen_nSig_1_K_15", "eqtl_MVphen_nSig_1_K_20", 
                   "eqtl_MVphen_nSig_1_K_30", "eqtl_MVphen_nSig_1_K_40", "eqtl_MVphen_nSig_2_K_15", 
                   "eqtl_MVphen_nSig_2_K_20", "eqtl_MVphen_nSig_2_K_30", "eqtl_MVphen_nSig_2_K_40") 
llmean.splitmat.true <- sapply(compl[grepl(Data, names(compl))], function(x) colMeans(x$llmat.raw[, 1:10], na.rm = T))
llmean.splitmat.true <- llmean.splitmat.true[, model_ordered]
# llmean.splitmat.true <- llmean.splitmat.true[, order(match(colnames(llmean.splitmat.true), c("UV", "XD", "mash", methnam_code)))]
true.mns <- colMeans(llmean.splitmat.true, na.rm = T)
true.ses <- apply(llmean.splitmat.true, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))
true_mn_with_cis <- true.mns
format.cvlik.fn_eqtl <- function(v) formatC(v, digits = 2, format = "f")
which_row_best_cvlik_eqtl <- which.max(true.mns)

res_to_fill$`CV Log Likelihood` <- paste0(format.cvlik.fn_eqtl(true.mns),
                                          " (", format.cvlik.fn_eqtl(true.mns - 2 * true.ses), ", ", format.cvlik.fn_eqtl(true.mns + 2 * true.ses), ")")
res_to_fill$`CV Log Likelihood`[which_row_best_cvlik_eqtl] <- paste0("\\underline{", res_to_fill$`CV Log Likelihood`[which_row_best_cvlik_eqtl], "}")

tabout_eqtl <- print(xtable(res_to_fill, label = "tab:cvlik_impc", align = rep("l", ncol(restaball.lik) + 1),
                           caption = "Comparison of cross-validated log likelihood across MV methods"),
                    caption.placement = "top", 
                    sanitize.text.function = function(x){x},
                    include.rownames = F, 
                    floating = FALSE)
cat(tabout_eqtl, file = paste(control$dropbox_text_numbers_dir, "/eqtl_cvlik_table.txt", sep = ""))

# length(true.mns)

# }

