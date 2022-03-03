file_preproc <- file.path(control$data_dir, "impc", "IMPC_results_preprocessed.RDS")
import.raw.data <- TRUE
if(import.raw.data){
  temp <- tempfile()
  download.file(url = "https://www.dropbox.com/s/j7a0b3ey855704o/IMPC_ALL_statistical_results.csv.gz?dl=1", temp, mode = "wb")
  unzipped_file <- file.path(control$data_dir, "impc", "IMPC_ALL_statistical_results.csv")
  R.utils::decompressFile(temp, destname = unzipped_file, ext = "gz", FUN=gzfile)
  unlink(temp)
  rin <- read.csv(file = unzipped_file)
  Data_all <- readRDS(file = control$Data_all_file)
  resimp <- readRDS(file = paste0(control$global_res_dir, "/resimp_comb.RDS"))  
  load(file = uv.results.Y.S)
  load(file = file.glob.res)
  genemap <- Data_all$impc$genemap
  cenmap <- Data_all$impc$cenmap
  resimp <- resimp[!grepl("fac", resimp$ph) & resimp$line.type == "trueMut", ]
  ebi.cenmap <- data.frame(ebinam = c("WTSI", "UC Davis", "TCP", "CCP-IMG", "KMPC", "MARC", "JAX", 
                                      "MRC Harwell", "ICS", "HMGU", "RBRC", "BCM"),
                            mynam = c("Wtsi", "Ucd", "Tcp", "CCP-IMG", "KMPC", "MARC", "J", 
                                      "H", "Ics", "Gmc", "Rbrc", "Bcm"))
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
  prob.res <- resimp[which(!resimp$testid %in% rin$testid &!resimp$imputed), ]
  fields.to.hugh <- c("ph", "geno", "cen", "cenlong", "geno.sh", "zyg")
  # These are the genes missing or incomplete in EBI data
  geno.not.there <- prob.res[!prob.res$geno.sh %in% rin$geno.sh, fields.to.hugh]
  geno.there.but.missing.params <- prob.res[prob.res$geno.sh %in% rin$geno.sh, fields.to.hugh]
  saveRDS(rin, file = file_preproc)
  # save(rin, file = file.ebi.impc.results)
  # save(genemap, geno.there.but.missing.params, geno.not.there, 
  #      file = "C:/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/data_out/records_george_still_cannot_find_in_ebi_v11.RData")
}

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


tab.uv.eb.comp <- table(resimp$uv.perm.signsig, resimp[, mv_signsig_name])

resimp.disagree <- resimp[which(resimp[, mv_signsig_name] != 0 & 
                                  resimp$ebi.signsig != 0 & 
                                  resimp[, mv_signsig_name] != resimp$ebi.signsig), ]
resimp.disagree.uvmv <- resimp[which(resimp[, mv_signsig_name] != 0 & 
                                       resimp$uv.perm.signsig != 0 & 
                                       resimp[, mv_signsig_name] != resimp$uv.perm.signsig), ]
table(resimp.all.agree$eb.perm.signsig)

resimp.prior <- resimp#.all.agree
resimp.disagree$prior.p.pos <- sapply(resimp.disagree$ph, function(ph) 
  sum(resimp.prior$ebi.signsig[resimp.prior$ph == ph] == 1 | resimp.prior$uv.perm.signsig[resimp.prior$ph == ph] == 1, na.rm = T) / 
    sum(resimp.prior$ebi.signsig[resimp.prior$ph == ph] != 0 | resimp.prior$uv.perm.signsig[resimp.prior$ph == ph] != 0, na.rm = T))

p.eb <- resimp.disagree$prior.p.pos^as.numeric(resimp.disagree[, mv_signsig_name] == 1) * 
  (1 - resimp.disagree$prior.p.pos)^as.numeric(resimp.disagree[, mv_signsig_name] == -1)
p.ebi <- resimp.disagree$prior.p.pos^as.numeric(resimp.disagree$ebi.signsig == 1) * 
  (1 - resimp.disagree$prior.p.pos)^as.numeric(resimp.disagree$ebi.signsig == -1)
table(p.eb > p.ebi)

plot(p.eb, p.ebi, xlim = c(0, 1))
abline(v = .5)
bf.eb.ebi <- prod(p.eb) / prod(p.ebi)
save.num <- c("bf.eb.ebi")
for(numc in save.num)
  write.table(prettyNum(round(eval(as.name(numc)), 2), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)


post.eb.mod.prob <- prod(p.eb) / (prod(p.eb) + prod(p.ebi))
save.num4 <- c("post.eb.mod.prob")
for(numc in save.num4)
  write.table(prettyNum(round(eval(as.name(numc)), 4), big.mark = ","), file = paste(dir.save, "/", numc, ".txt", sep = ""),
              col.names = F, row.names = F, quote = F)

install.packages("plotrix")

res_ebi <- resimp[!is.na(resimp$uv.perm.signsig), 
                  c("ebi.signsig", "uv.perm.signsig", mv_signsig_name)]
res_missing <- resimp[which(resimp$imputed &
                              !is.na(resimp[, mv_signsig_name])), mv_signsig_name,
                      drop = FALSE]
meths <- c("ebi", "uv", "mv")
names(res_ebi) <- meths
n_ebi_tests <- nrow(res_ebi)
n_missing_mvphen_tests <- nrow(res_missing)
p_missing <- n_missing_mvphen_tests / (n_ebi_tests + n_missing_mvphen_tests)
p_measured <- 1 - p_missing


methnamv <- c(ebi = "ebi.signsig", uv = "uv.perm.signsig", mv = mv_signsig_name)
p_pos <- rad_pos <- list()
for (meth in meths) {
  # p_pos[[meth]] <- mean(res_ebi[, meth] != 0)
  p_pos[[meth]] <- mean(!is.na(resimp[, methnamv[meth]]) & resimp[, methnamv[meth]] != 0)
  rad_pos[[meth]] <- sqrt(p_pos[[meth]] / pi)
}
p_miss_and_pos <- mean(res_missing != 0)


n_obs <- sum(!resimp$imputed)
n_miss <- n_tot - n_obs
n_tot <- nrow(resimp)
n_ebi_pos <- sum(!is.na(resimp[, methnamv["ebi"]]) & resimp[, methnamv["ebi"]] != 0)
n_uv_pos <- sum(!is.na(resimp[, methnamv["uv"]]) & resimp[, methnamv["uv"]] != 0)
n_mv_pos_obs <- sum(!resimp$imputed & resimp[, methnamv["mv"]] != 0)
n_ebi_pos_mv_neg <- sum(!resimp$imputed & resimp[, methnamv["ebi"]] != 0 & resimp[, methnamv["mv"]] == 0 )


n_ebi_pos_mv_pos <- sum(na.omit(!resimp$imputed & resimp[, methnamv["ebi"]] != 0 & resimp[, methnamv["mv"]] != 0 ))
n_uv_pos_mv_pos <- sum(na.omit(!resimp$imputed & resimp[, methnamv["uv"]] != 0 & resimp[, methnamv["mv"]] != 0 ))
n_uv_pos_ebi_pos <- sum(na.omit(!resimp$imputed & resimp[, methnamv["uv"]] != 0 & resimp[, methnamv["ebi"]] != 0 ))

p23 <- n_ebi_pos_mv_pos / n_tot
p12 <- n_uv_pos_ebi_pos / n_tot
p13 <- n_uv_pos_mv_pos / n_tot
p1 <- n_uv_pos / n_tot
p2 <- n_ebi_pos / n_tot
p3 <- n_mv_pos_obs / n_tot
r1 <- sqrt(n_uv_pos / n_tot / pi)
r2 <- sqrt(n_ebi_pos / n_tot / pi)
r3 <- sqrt(n_mv_pos_obs / n_tot / pi)
d = .02
area <- function (r1, r2, d) { 
  l = sqrt(((r1 + r2)^2 - d^2) * (d^2 - (r1 - r2)^2)) / d
  alp <- asin(l / 2 / r1)
  bet <- asin(l / 2 / r2)
  area <- (alp - sin(alp) * cos(alp)) * r1^2 + (bet - sin(bet) * cos(bet)) * r2^2
  return(area)
}

area <- function (r, R, d) { 
  area <- r^2 * acos((d^2 + r^2 - R^2) / (2 * d * r)) + 
          R^2 * acos((d^2 + R^2 - r^2) / (2 * d * R)) -
          .5 * sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R))  
  return(area)
}

r1
r2
d12seq <- seq(abs(r1 - r2), r1 + r2, by = .0001)
d12 <- d12seq[which.min(abs(area(r1, r2, d12seq) - p12))]
d13seq <- seq(abs(r1 - r3), r1 + r3, by = .0001)
d13 <- d13seq[which.min(abs(area(r1, r3, d13seq) - p13))]
d23seq <- seq(abs(r2 - r3), r2 + r3, by = .0001)
d23 <- d23seq[which.min(abs(area(r2, r3, d23seq) - p23))]

plot(abs(area(r2, r3, d23seq) - p23))
rr1 <- 1  
rr2 <- 2  
dseq <- seq(abs(rr1 - rr2), rr1 + rr2, by = .0001)
plot(dseq, area(rr1, rr2, dseq))

x3 <- p_obs / 2
y3 <- .5
x1 <- x3 + d13
y1 <- y3
cd <- d13
ad <- d12
bd <- d23                                                       
x2 <- x3 + (ad^2 - bd^2 - cd^2) / (-2 * cd)
y2 <- y3 + sqrt(bd^2 - (x2 - x3)^2)
n_mv_pos_miss <- sum(resimp$imputed & resimp[, methnamv["mv"]] != 0)
p_miss <- n_miss / n_tot
p_obs <- n_obs / n_tot
table(resimp[, methnamv["ebi"]], resimp[, methnamv["mv"]])

pdf(file = paste0(control$dropbox_figure_dir, "/venn_diagrams.pdf"), 12, 8.5)
par(oma = c(1, 1, 1, 1), mar = c(0, 0, 6, 20))
col_use <- rgb(red = c(1, 0.2, 0), green = c(0, 0.2, 0), blue = c(0, 0.2, 1), alpha = .45)
# viridis::rocket(n = 3, alpha = .6)#[sample(1:3, 3)]
#col_use <- viridis::turbo(n = 3, alpha = .5)#[sample(1:3, 3)]
names(col_use) <- meths
plot(0, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", ty = "l", xaxs = "i", yaxs = "i")
polygon(x = c(0, rep(p_obs, 2), 0, 0), y = c(0, 0, 1, 1, 0), col = "white", border = 1)
plotrix::draw.circle(x = x3, y = y3, radius = r3, nv = 10000, border = NA, 
                     col = col_use["mv"], lty = 1, density = NULL, angle = 45, lwd = 1)
plotrix::draw.circle(x = x1, y = y1, radius = r1, nv = 10000, border = NA, 
                     col = col_use["uv"], lty = 1, density = NULL, angle = 45, lwd = 1)
plotrix::draw.circle(x = x2, y = y2, radius = r2, nv = 10000, border = NA, 
                     col = col_use["ebi"], lty = 1, density = NULL, angle = 45, lwd = 1)
plotrix::draw.circle(x = p_obs + p_miss / 2, y = .5, radius = sqrt(n_mv_pos_miss / n_tot / pi), nv = 10000, border = NA, 
                     col = col_use["mv"], lty = 1, density = NULL, angle = 45, lwd = 1)
cexax <- 1.4
mtext(side = 3, text = c("Observed data", "Missing data"), at = c(p_obs / 2, p_obs + p_miss / 2), line = 1, cex = cexax)
par(xpd = NA)
# legend(x = .5, y = 1.09, legend = c("IMPC database", "UV method", "MV method"), col = col_use, pch = 19, pt.cex = 2, xjust = .5)
legend(x = 1.03, y = .75, legend = c("IMPC database", "UV method", "MV method"), col = col_use,
        pch = 19, pt.cex = 2, xjust = 0, cex = cexax)
par(xpd = F)
dev.off()

# plotrix::draw.circle(x = p_obs / 2, y = .5, radius = sqrt(n_uv_pos / n_tot / pi), nv = 10000, border = NULL, 
#                      col = col_use["uv"], lty = 1, density = NULL, angle = 45, lwd = 1)
for (meth in c("mv", "ebi", "uv")) {
  plotrix::draw.circle(x = p_obs / 2, y = .5, radius = rad_pos[[meth]], nv = 10000, border = col_use[meth], 
                     col = col_use[meth], lty = 1, density = NULL, angle = 45, lwd = 1)
  plotrix::draw.circle(x = p_obs + p_miss / 2, y = .5, radius = rad_pos[[meth]], nv = 10000, border = col_use[meth], 
                       col = col_use[meth], lty = 1, density = NULL, angle = 45, lwd = 1)
}
































































