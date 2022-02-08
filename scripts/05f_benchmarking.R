# load(file = file.runtab)
# load(file = file.resl.comp)
# load(file = file.compl)
# load(file = uv.results.Y.S)
# load(file = paste0(global.res.dir, "/resimpl_comb.RData"))  

# for(i in 1:length(default.parameters[[Data]]))
#   assign(names(default.parameters[[Data]][i]), default.parameters[[Data]][[i]])

linemap <- Data_all$impc$linemap
methnam <- "MVphen"
truemuts <- linemap$geno[linemap$line.type == control$nam.truemut]
negcons <- linemap$geno[linemap$line.type == control$nam.negcon]
llmean.splitmat.true <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[truemuts, 1:10], na.rm = T))
llmean.splitmat.null <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.raw[negcons, 1:10], na.rm = T))
# llmean.splitmat.zero <- sapply(compl[grepl("impc", names(compl))], function(x) colMeans(x$llmat.zero[truemuts, ], na.rm = T))
true.mns <- colMeans(llmean.splitmat.true, na.rm = T)
true.ses <- apply(llmean.splitmat.true, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))
null.mns <- colMeans(llmean.splitmat.null, na.rm = T)
null.ses <- apply(llmean.splitmat.null, 2, function(v) sd(v, na.rm = T) / sqrt(length(v)))

err.data.type <- "comb"#for(err.data.type in c("comb", "single")){
# load(file = paste0(global.res.dir, "/restabl_", err.data.type, ".RData"))  
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
restaball$hit.imp.ci <- paste0(formatfn(restaball$mvhitimp), " (", formatfn(restaball$hit.rate.imp.ci.l), "-", formatfn(restaball$hit.rate.imp.ci.u), ")")
restaball$hit.nonimp.ci <- paste0(formatfn(restaball$mvhitnonimp), " (", formatfn(restaball$hit.rate.nonimp.ci.l), "-", formatfn(restaball$hit.rate.nonimp.ci.u), ")")
restaball$line.fdr.ci <- paste0(formatfn(restaball$line.fdr.est), " (", formatfn(restaball$line.fdr.ci.l), "-", formatfn(restaball$line.fdr.ci.u), ")")
restaball$fdr.ci <- paste0(formatfn(restaball$fdr.est), " (", formatfn(restaball$fdr.ci.l), "-", formatfn(restaball$fdr.ci.u), ")")
restaball$fsr.ci <- paste0(formatfn(restaball$ref.lines.post), " (", formatfn(restaball$ref.lines.post.l), "-", formatfn(restaball$ref.lines.post.u), ")")
restaball$Method <- NA
restaball$facmod <- "FA"#ifelse(grepl("_fa_", restaball$meth), "FA", "PCA")
restaball$Method[grepl(methnam, restaball$meth)] <- methnam
# restaball$Method[grepl("em.fit", restaball$meth)] <- paste0(methnam, " (", restaball$facmod[grepl("em.fit", restaball$meth)], ")")
restaball$Method[grepl("mash", restaball$meth)] <- "mash"
restaball$Method[grepl("XD", restaball$meth)] <- "XD"
restaball$Method[grepl("uv", restaball$meth)] <- "UV"
restaball$mvhitimp[grepl("uv", restaball$meth)] <- NA
restaball.out <- restaball
restaball$Test <- paste0(ifelse(restaball$err.rate.meth == "perm", 1, 2), ifelse(restaball$test.stat == "z", "A", "B"))
restaball$S <- sapply(strsplit(restaball$meth, spl = "_"), function(x) x[4])
restaball$S[restaball$Method == "mash"] <- control$default_parameters$impc$P + 10
restaball$S[restaball$Method == "UV"] <- ""
restaball$K <- sapply(strsplit(restaball$meth, spl = "_"), function(x) x[6])
restaball$fm <- "fa"#sapply(strsplit(restaball$meth, spl = "_"), function(x) x[6])
restaball <- restaball[which(restaball$fm != "pca" | is.na(restaball$fm)), ]
restaball$Control <- ifelse(restaball$err.rate.meth == "perm", "$\\mathrm{Fdr}_{\\mathrm{complete}}\\leq 5\\%$", "$lfsr \\leq 5\\%$")
restaball$Statistic <- ifelse(restaball$test.stat == "z", "$z$", "$lfsr$")
restaball$cvlik <- true.mns[restaball$meth]
restaball$cvlik.l <- true.mns[restaball$meth] - 2 * true.ses[restaball$meth]
restaball$cvlik.u <- true.mns[restaball$meth] + 2 * true.ses[restaball$meth]
restaball$cvlik.ci <- paste0(format.cvlik.fn(restaball$cvlik), " (", format.cvlik.fn(restaball$cvlik.l), ",", format.cvlik.fn(restaball$cvlik.u), ")")
restaball <- restaball[order(match(restaball$Method, c("UV", "XD", "mash", methnam))), ]
rates.format <- c("mvhitimp", "mvhitnonimp", "line.fdr.est", "fdr.est", "fdr.est.imp", "fdr.est.nonimp")
for(j in rates.format){
  if(is.numeric(restaball[, j]))
    restaball[, j] <- formatC(restaball[, j] * 100, digits = 1, format = "f")
}
tail(restaball)
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
restaball.out

tabout <- print(xtable(restaball.out, label = "tab:hitrates", align = rep("r", ncol(restaball.out) + 1),
                       caption = "Hit rates and error rate comparison across methods"),
                caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
cat(tabout, file = paste(control$dropbox_text_numbers_dir, "/hitrates_combined.txt", sep = ""))


hitList <- split(restaball.out[, 4:ncol(restaball.out)], f = restaball.out$Test)
attr(hitList, "subheadings") <- paste0(names(hitList), ". Controlling ", 
                                       restaball.out$Control[match(names(hitList), restaball.out$Test)],
                                       " using ", restaball.out$Statistic[match(names(hitList), restaball.out$Test)],  " statistic ")
xList <- xtableList(hitList, label = paste0("tab:hitrates.", err.data.type), 
                    align = gsub(" llllll", "llllll|", paste0(" ", paste0(rep("l", ncol(hitList[[1]]) + 1), collapse = ""))),
                    caption = caption)
tabout.test <- print(xList, table.placement = "hp", colnames.format = "multiple", add.to.row = addtorow,
                     caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
tabout.test.with.header <- gsub("\n\\\\hline\nMe", 
                                paste0("\n\\\\hline\n& & &\\\\multicolumn\\{2\\}\\{c\\}\\{Hit rate  in \\\\% when data are\\}&",
                                       "\\\\multicolumn\\{3\\}\\{|c\\}\\{Estimated error rate in \\\\% (95\\\\% CI) \\}\\\\\\\\\n",
                                       # "\\\\cline\\{3\\-4\\}\\\\cline\\{5\\-7\\}\\\nMe"), tabout.test)
                                       "\\\nMe"), tabout.test)
tabout.test.with.header
cat(tabout.test.with.header, file = paste(control$dropbox_text_numbers_dir, "/hitrates_", err.data.type, ".txt", sep = ""))


colkeep.lik <- c("Method", "S", "K", "cvlik.ci")#, "hit.nonimp.ci", "hit.imp.ci")
restaball.lik <- unique(restaball[, colkeep.lik])
restaball.lik <- restaball.lik[restaball.lik$Method != "UV", ]
colnames(restaball.lik) <- colkeepmap[match(colkeep.lik, colkeepmap$keep), "keep.as"]
restaball.lik

tabout.lik <- print(xtable(restaball.lik, label = "tab:cvlik_impc", align = rep("l", ncol(restaball.lik) + 1),
                           caption = "Comparison of cross-validated log likelihood across MV methods"),
                    caption.placement = "top", sanitize.text.function = function(x){x}, include.rownames = F)
cat(tabout.lik, file = paste(control$dropbox_text_numbers_dir, "/cvlik_table.txt", sep = ""))

# }

