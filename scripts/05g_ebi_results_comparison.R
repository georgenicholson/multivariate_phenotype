# data <- "impc"
# xdir <- ifelse(Sys.info()["sysname"] == "Linux", "/mnt/x", "X:")
# source(paste0(xdir, "/projects/impc_mv_analysis/R_files/impc_mv_parameters.R"))
load.raw.data <- F
if(load.raw.data){
  rin <- read.csv(file = paste0(xdir, "/projects/impc_mv_analysis/data_in/impc_results_ebi/IMPC_ALL_statistical_results.csv"))
#   save(rin, file = paste0(processed.data.dir, "/ebi_impc_results_v11.RData"))
# }
  
  load(file = paste0(processed.data.dir, "/ebi_impc_results_v11.RData"))
  load(file = uv.results.Y.S)
  load(file = file.glob.res)
  genemap <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/gene_info_impc.tsv"), sep = "\t", header = T)
  # rin$impc_genotype_id_2 <- genemap[match(rin$marker_symbol, genemap$gene_symbol), "genotype_id"]
  # table(is.na(rin$impc_genotype_id_2), is.na(rin$impc_genotype_id))
  # mean(is.na(rin$impc_genotype_id_2) & !is.na(rin$impc_genotype_id))
  # mean(is.na(rin$impc_genotype_id) & !is.na(rin$impc_genotype_id_2))
  
  head(rin$marker_symbol)
  str(genemap)
  resimp <- resimp[!resimp$ph %in% facnam & resimp$line.type == "trueMut", ]
  ebi.cenmap <- data.frame(ebinam = c("WTSI", "UC Davis", "TCP", "CCP-IMG", "KMPC", "MARC", "JAX", 
    "MRC Harwell", "ICS", "HMGU", "RBRC", "BCM"),
                            mynam = c("Wtsi", "Ucd", "Tcp", "CCP-IMG", "KMPC", "MARC", "J", "H", "Ics", "Gmc", "Rbrc", "Bcm"))
  ebi.cenmap$mynum <- cenmap[match(ebi.cenmap$mynam, cenmap$nam), "cen"]
  rin$mycen <- ebi.cenmap[match(rin$phenotyping_center, ebi.cenmap$ebinam), "mynum"]
  genemap$cen <- resimp[match(genemap$genotype_id, resimp$geno.sh), "cen"]
  genemap$mgi_cen <- paste(genemap$gene_id, genemap$cen, sep = "_")
  rin$mgi_cen <- paste(rin$marker_accession_id, rin$mycen, sep = "_")
  rin$impc_genotype_id <- genemap[match(rin$colony_id, genemap$genotype), "genotype_id"]
  rin$zyg <- ifelse(rin$zygosity == "hemizygote", 2, ifelse(rin$zygosity == "heterozygote", 0, 1))
  rin$geno.sh <- rin$impc_genotype_id
  rin$geno <- paste(rin$geno.sh, rin$zyg, sep = "_")
  rin$ph <- rin$parameter_stable_id
  rin$testid <- paste(rin$mycen, rin$ph, rin$geno.sh, rin$zyg, sep = "_")
  rin$testid_nocen <- paste(rin$ph, rin$geno.sh, rin$zyg, sep = "_")
  resimp$testid_nocen <- paste(resimp$ph, resimp$geno.sh, resimp$zyg, sep = "_")
  mean(resimp$testid[!resimp$imputed] %in% rin$testid)
  mean(resimp$testid_nocen %in% rin$testid_nocen)
  # prob.res <- resimp[sample(which(!resimp$testid %in% rin$testid &!resimp$imputed), 100), ]
  prob.res <- resimp[which(!resimp$testid %in% rin$testid &!resimp$imputed), ]
  fields.to.hugh <- c("ph", "geno", "cen", "cenlong", "geno.sh", "zyg")
  # just.genofields.to.hugh <- c("ph", "geno", "cen", "cenlong", "geno.sh", "zyg")
  geno.not.there <- prob.res[!prob.res$geno.sh %in% rin$geno.sh, fields.to.hugh]
  geno.there.but.missing.params <- prob.res[prob.res$geno.sh %in% rin$geno.sh, fields.to.hugh]
  # genemap
  str(geno.not.there)
  str(geno.there.but.missing.params)
  str(rin)
  save(rin, file = file.ebi.impc.results)
  # save(genemap, geno.there.but.missing.params, geno.not.there, 
  #      file = "C:/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/data_out/records_george_still_cannot_find_in_ebi_v11.RData")
}


load(file = file.ebi.impc.results)
load(file = file.glob.res)
resimp[, paste0("ebi.", c("mn", "se", "p"))] <- rin[match(resimp$testid, rin$testid), 
                                                    c("genotype_effect_parameter_estimate", "genotype_effect_stderr_estimate", "p_value")]
resimp$ebi.t <- resimp$ebi.mn / resimp$ebi.se
resimp$ebi.signsig <- sign(resimp$ebi.mn) * (resimp$ebi.p < 1e-4)


tab.eb.ebi.comp <- table(resimp$eb.perm.signsig, resimp$ebi.signsig)
tab.uv.ebi.comp <- table(resimp$uv.perm.signsig, resimp$ebi.signsig)
# str(tab.uv.ebi.comp)
# 
# table(resimp$line.type)
# str(resimp)
# resimp.comp <- resimp[!is.na(resimp$uv.perm.signsig) & !is.na(resimp$eb.perm.signsig) & !is.na(resimp$ebi.signsig), ]
# ebi.hit.rate.p <- mean(resimp$ebi.signsig != 0, na.rm = T)
# uv.hit.rate.p <- mean(resimp$uv.perm.signsig != 0, na.rm = T)
# eb.hit.rate.p <- mean(resimp$eb.perm.signsig[!is.na(resimp$uv.perm.signsig)] != 0, na.rm = T)

prettyNum(tab.eb.ebi.comp, big.mark = ",")

tab.eb.ebi.comp

n.eb.ebi.disagree <- tab.eb.ebi.comp["-1", "1"] + tab.eb.ebi.comp["1", "-1"]
prop.eb.ebi.disagree <- n.eb.ebi.disagree / (n.eb.ebi.disagree + tab.eb.ebi.comp["-1", "-1"] + tab.eb.ebi.comp["1", "1"])
ebi.hit.rate <- mean(resimp$ebi.signsig != 0, na.rm = T)
dir.save <- revision.text.numbers
save.prop <- c("ebi.hit.rate", "prop.eb.ebi.disagree")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 1, format = "f"), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
save.num <- c("n.eb.ebi.disagree")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)



save.num <- c("n.sig.eb.non", "n.sig.eb.imp", "n.sig.uv", "n.sig.eb.tot", "n.ko.line.measured", "n.ph.measured")
for(numc in save.num)
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)
# save.num2 <- c("mean.s.measured.phens", "sd.y.measured.phens")
# for(numc in save.num2)
#   write.table(formatC(eval(as.name(numc)), format = "f", digits = 2), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
#               col.names = F, row.names = F, quote = F)
# 
# 
# save.fold <- c("eb.fold.increase.", "imp.fold.increase", "fold.increase")
# for(numc in save.fold)
#   write.table(formatC(eval(as.name(numc)), digits = 1, format = "f"), file = paste(dir.save, "/", numc, ".txt", sep = ""), 
#               col.names = F, row.names = F, quote = F)
library(xtable)
# addtorow <- addtocol.mv <- list()
# addtorow$pos <- addtocol.mv$pos <- list(-1)
# addtorow$command <- '&\\multicolumn{3}{c}{Existing IMPC call}\\\\'

tabout.uv <- print(xtable(prettyNum(tab.uv.ebi.comp, big.mark = ","), label = "tab:uv_ebi_comparison", 
                       caption = "Comparison of signed phenotype hits between our UV model (left) and the existing phenotype calls in the IMPC database (top)"),
                   caption.placement = "top")
cat(tabout.uv, file = paste(revision.text.numbers, "/uv_ebi_comp_tab.txt", sep = ""))
tabout.mv <- print(xtable(tab.eb.ebi.comp, label = "tab:eb_ebi_comparison", 
                          caption = "Comparison of signed phenotype hits between our MV model (left) and the existing phenotype calls in the IMPC database (top)"),
                   caption.placement = "top")
cat(tabout.mv, file = paste(revision.text.numbers, "/eb_ebi_comp_tab.txt", sep = ""))


tab.uv.eb.comp <- table(resimp$uv.perm.signsig, resimp$eb.perm.signsig)

resimp.disagree <- resimp[which(resimp$eb.perm.signsig != 0 & resimp$ebi.signsig != 0 & resimp$eb.perm.signsig != resimp$ebi.signsig), ]
resimp.disagree.uvmv <- resimp[which(resimp$eb.perm.signsig != 0 & resimp$uv.perm.signsig != 0 & resimp$eb.perm.signsig != resimp$uv.perm.signsig), ]
# resimp.all.agree <- resimp[which(resimp$uv.perm.signsig == resimp$ebi.signsig), ]
table(resimp.all.agree$eb.perm.signsig)


resimp.all.agree[resimp.all.agree$ph == "IMPC_ECH_011_001", ]

#resimp.disagree$prior.p.pos <- sapply(resimp.disagree$ph, function(ph) sum(resimp$ebi.signsig[resimp$ph == ph] == 1, na.rm = T) / sum(resimp$ebi.signsig[resimp$ph == ph] != 0, na.rm = T))
# resimp.disagree$prior.p.pos <- sapply(resimp.disagree$ph, function(ph) sum(resimp$eb.perm.signsig[resimp$ph == ph] == 1, na.rm = T) / sum(resimp$eb.perm.signsig[resimp$ph == ph] != 0, na.rm = T))
# resimp.disagree$prior.p.pos <- sapply(resimp.disagree$ph, function(ph) sum(resimp$uv.perm.signsig[resimp$ph == ph] == 1, na.rm = T) / sum(resimp$uv.perm.signsig[resimp$ph == ph] != 0, na.rm = T))
# 


resimp.prior <- resimp#.all.agree
resimp.disagree$prior.p.pos <- sapply(resimp.disagree$ph, function(ph) 
  sum(resimp.prior$ebi.signsig[resimp.prior$ph == ph] == 1 | resimp.prior$uv.perm.signsig[resimp.prior$ph == ph] == 1, na.rm = T) / 
    sum(resimp.prior$ebi.signsig[resimp.prior$ph == ph] != 0 | resimp.prior$uv.perm.signsig[resimp.prior$ph == ph] != 0, na.rm = T))


# resimp.disagree.uvmv$prior.p.pos <- sapply(resimp.disagree.uvmv$ph, function(ph) 
#   sum(resimp.prior$ebi.signsig[resimp.prior$ph == ph] == 1 | resimp.prior$uv.perm.signsig[resimp.prior$ph == ph] == 1, na.rm = T) / 
#     sum(resimp.prior$ebi.signsig[resimp.prior$ph == ph] != 0 | resimp.prior$uv.perm.signsig[resimp.prior$ph == ph] != 0, na.rm = T))

# resimp.disagree$prior.p.pos <- sapply(resimp.disagree$ph, function(ph) sum(resimp.all.agree$uv.perm.signsig[resimp.all.agree$ph == ph] == 1, na.rm = T) /
#                                         sum(resimp.all.agree$uv.perm.signsig[resimp.all.agree$ph == ph] != 0, na.rm = T))

p.eb <- resimp.disagree$prior.p.pos^as.numeric(resimp.disagree$eb.perm.signsig == 1) * 
  (1 - resimp.disagree$prior.p.pos)^as.numeric(resimp.disagree$eb.perm.signsig == -1)
p.ebi <- resimp.disagree$prior.p.pos^as.numeric(resimp.disagree$ebi.signsig == 1) * 
  (1 - resimp.disagree$prior.p.pos)^as.numeric(resimp.disagree$ebi.signsig == -1)


# p.uv.uvmv <- resimp.disagree.uvmv$prior.p.pos^as.numeric(resimp.disagree.uvmv$uv.perm.signsig == 1) * 
#   (1 - resimp.disagree.uvmv$prior.p.pos)^as.numeric(resimp.disagree.uvmv$uv.perm.signsig == -1)
# p.eb.uvmv <- resimp.disagree.uvmv$prior.p.pos^as.numeric(resimp.disagree.uvmv$eb.perm.signsig == 1) * 
#   (1 - resimp.disagree.uvmv$prior.p.pos)^as.numeric(resimp.disagree.uvmv$eb.perm.signsig == -1)


table(p.eb > p.ebi)

plot(p.eb, p.ebi, xlim = c(0, 1))
abline(v = .5)
bf.eb.ebi <- prod(p.eb) / prod(p.ebi)
save.num <- c("bf.eb.ebi")
for(numc in save.num)
  write.table(prettyNum(round(eval(as.name(numc)), 0), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


post.eb.mod.prob <- prod(p.eb) / (prod(p.eb) + prod(p.ebi))
save.num4 <- c("post.eb.mod.prob")
for(numc in save.num4)
  write.table(prettyNum(round(eval(as.name(numc)), 4), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)




mean(is.na(resimp$ebi.signsig))


fdr.est.tab(tab.eb.ebi.comp)



table(resimp$eb.perm.signsig, resimp$ebi.signsig)
tab.imp <- table(resimp$eb.perm.signsig[is.na(resimp$uv.perm.signsig)], resimp$ebi.signsig[is.na(resimp$uv.perm.signsig)])
fdr.est.tab(tab.imp)

table(resimp$eb.perm.signsig[is.na(resimp$uv.perm.signsig)], resimp$ebi.signsig[is.na(resimp$uv.perm.signsig)])
table(resimp$uv.perm.signsig[!is.na(resimp$uv.perm.signsig) & !is.na(resimp$ebi.signsig)])
table(resimp$ebi.signsig[!is.na(resimp$uv.perm.signsig) & !is.na(resimp$ebi.signsig)])


# resimp[resimp$]








annos <- rin
unique(annos[annos$colony_id == genemap[genemap$genotype_id == geno.not.there[1,"geno.sh"],"genotype"],"marker_accession_id"])
mean(resimp$geno %in% rin$geno)
length(unique(resimp$geno[!resimp$geno %in% rin$geno]))
str(rin$geno)
str(geno.not.there)
str(geno.there.but.missing.params)
table(geno.not.there$ph)
table(geno.there.but.missing.params$ph)
npl <- 10000
plot(resimp$ebi.t[1:npl], resimp$uv.t[1:npl])




tab.eb.ebi.imp <- table(resimp$eb.perm.signsig[resimp$imputed], resimp$ebi.signsig[resimp$imputed])
tab.eb.ebi.non <- table(resimp$eb.perm.signsig[!resimp$imputed], resimp$ebi.signsig[!resimp$imputed])
tab.eb.ebi.imp
fdr.est.tab(tab.eb.ebi.imp)
tab.eb.ebi.non
fdr.est.tab(tab.eb.ebi.non)



table(resimp$eb.perm.signsig, resimp$uv.perm.signsig)

rsub <- rin[rin$parameter_stable_id %in% ph.use & rin$impc_genotype_id %in% resimp$geno.sh, ]

mean(resimp$geno.sh %in% rin$impc_genotype_id)
unique(resimp$ph[!resimp$ph %in% rin$parameter_stable_id])
mean(resimp$cen %in% rin$mycen)
table(rin$zyg)
mean(resimp$testid %in% rin$testid)
table(resimp[!resimp$testid %in% rin$testid, "cen"])
head(rin$testid)
head(resimp$testid)

mean(resimp$geno.sh %in% rin$impc_genotype_id)
mean(resimp$cen %in% rin$mycen)


dput(unique(rsub$phenotyping_center))


    
table(rsub$phenotyping_center)
cenmap

resimp$ebi.mn[1:100]


plot()
table(rsub$phenotyping_center_id)
cenmap

sort(table(rsub$parameter_stable_id[rsub$interaction_significant == "true"]))

colnames(rsub)[grep("sig", colnames(rsub))]
str(rsub)
rsub$
str(resimp)
mean(resimp$geno.sh %in% rin$mutant_biological_model_id)
str(rin)

mean(!is.na(rin$impc_genotype_id))
ph.use
mean(genemap$gene_id %in% rin$marker_accession_id)
table(table(rsub$allele_accession_id))
namsort <- sort(colnames(rin))
str(rin[, namsort])

rin$z
range(rin$mutant_biological_model_id)
range(as.numeric(resimp$geno.sh))
str(genemap)

for(j in 1:ncol(rin))
  print(ph.use[1] %in% rin[, j])
#
