# .libPaths("X:/R_libraries")
# source("scripts/05_generate_results.R")
####################################
# Control parameters

sd_mult_big_eff_thresh <- 2
dirc <- "bo"
nperm <- 1000
fwer.th <- .05
methc <- c("ph", "proc", "fac")[1]
rerun.GO.analysis <- F
use.just.homs <- T   # To be want to do the analysis using just homozygotes?
fac.meth <- c("varimax", "promax")[1]



resimp <- readRDS(file = control$file.resimp)
resl.comp.fac <- readRDS(file = control$file_raw_factor_results)
phmap <- Data_all$impc$phmap
library("foreach")

ph.use <- Data_all$impc$phord
if(use.just.homs){
  true.use <- unique(resimp[resimp$line.type == "trueMut" & resimp$zyg == 1, "geno"])
} else {
  true.use <- unique(resimp[resimp$line.type == "trueMut", "geno"])
}
procord <- Data_all$impc$procord
nph <- length(ph.use)
ng <- length(true.use)
dirv <- c("up", "do", "bo")


sum(resimp$uv.sig, na.rm = T)
sum(resimp$eb.sig, na.rm = T)

resl.out <- list()
resl.out$eb <- resll$perm[[control$mv_meth_nam_use]]
resl.out$uv <- resll$perm$uv
resl.out[[fac.meth]] <- resl.comp.fac[[control$mv_meth_nam_use]]
facnam <- colnames(resl.out[[fac.meth]]$mn)
resl.out[[fac.meth]]$t <- resl.out[[fac.meth]]$mn / resl.out[[fac.meth]]$sd
resl.out[[fac.meth]]$th <- resimp[, paste0(fac.meth, ".th.final")][1]

ph_high <- names(sort(abs(resl.comp[[control$mv_meth_nam_use]]$facs.varimax[, "fac_7"]), decreasing = T))

Data_all$impc$phmap[match(ph_high[1:5], Data_all$impc$phmap$ph), ]

mvphen_effect_size_thresh <- outer(rep(1, nrow(resl.out$eb$mn)), apply(resl.out$eb$mn, 2, sd, na.rm = T)) * sd_mult_big_eff_thresh
uv_effect_size_thresh <- outer(rep(1, nrow(resl.out$uv$mn)), apply(resl.out$uv$mn, 2, sd, na.rm = T)) * sd_mult_big_eff_thresh
fac_effect_size_thresh <- outer(rep(1, nrow(resl.out[[fac.meth]]$mn)), apply(resl.out[[fac.meth]]$mn, 2, sd, na.rm = T)) * sd_mult_big_eff_thresh
resl.out$eb$signsig.bigeff <- (abs(resl.out$eb$t) > resl.out$eb$th & abs(resl.out$eb$mn) >= mvphen_effect_size_thresh) * sign(resl.out$eb$t)
resl.out$uv$signsig.bigeff <- (abs(resl.out$uv$t) > resl.out$uv$th & abs(resl.out$uv$mn) >= uv_effect_size_thresh) * sign(resl.out$uv$t)
resl.out[[fac.meth]]$signsig.bigeff <- (abs(resl.out[[fac.meth]]$t) > resl.out[[fac.meth]]$th & 
                                          abs(resl.out[[fac.meth]]$mn) >= fac_effect_size_thresh) * 
                                            sign(resl.out[[fac.meth]]$t)

sigl <- list()
sigl$uv <- sigl$uv.bigeff <- sigl$mv <- sigl$mv.bigeff <- sigl$f <- list()
sigl$uv$up <- resl.out$uv$signsig[true.use, ph.use] == 1
sigl$uv$do <- resl.out$uv$signsig[true.use, ph.use] == -1
sigl$uv$bo <- resl.out$uv$signsig[true.use, ph.use] != 0
sigl$mv$up <- resl.out$eb$signsig[true.use, ph.use] == 1
sigl$mv$do <- resl.out$eb$signsig[true.use, ph.use] == -1
sigl$mv$bo <- resl.out$eb$signsig[true.use, ph.use] != 0
sigl$mv.bigeff$up <- resl.out$eb$signsig.bigeff[true.use, ph.use] == 1
sigl$mv.bigeff$do <- resl.out$eb$signsig.bigeff[true.use, ph.use] == -1
sigl$mv.bigeff$bo <- resl.out$eb$signsig.bigeff[true.use, ph.use] != 0
sigl$uv.bigeff$up <- resl.out$uv$signsig.bigeff[true.use, ph.use] == 1
sigl$uv.bigeff$do <- resl.out$uv$signsig.bigeff[true.use, ph.use] == -1
sigl$uv.bigeff$bo <- resl.out$uv$signsig.bigeff[true.use, ph.use] != 0
sigl$f$up <- resl.out[[fac.meth]]$signsig.bigeff[true.use, facnam] == 1
sigl$f$do <- resl.out[[fac.meth]]$signsig.bigeff[true.use, facnam] == -1
sigl$f$bo <- resl.out[[fac.meth]]$signsig.bigeff[true.use, facnam] != 0

colSums(sigl$f$bo)
colSums(sigl$mv.bigeff$bo)


gene.all.num <- unique(sapply(strsplit(true.use, spl = "_"), function(v) v[1]))
genemap.in <- Data_all$impc$genemap[Data_all$impc$genemap$genotype_id %in% gene.all.num, ]
sym2eg <- BiocGenerics::toTable(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
sym2eg <- sym2eg[order(as.integer(sym2eg$gene_id)), ]

#Note taking the first Entrez ID in case of one-to-many sym2eg mapping 
genemap.in$entrez <- sym2eg[match(genemap.in$gene_symbol, sym2eg$symbol), "gene_id"]

############################################
# Note that there are currently some Symbols not being mapped to Entrez IDs
#Could not match symbol in map to Entrez for these
genemap.in[is.na(genemap.in$entrez), ]
genemap2 <- genemap.in[!is.na(genemap.in$entrez), ]
egGenomicLoc <- BiocGenerics::toTable(org.Mm.eg.db::org.Mm.egCHRLOC)
sym2eg[, c("loc", "chr")] <- egGenomicLoc[match(sym2eg$gene_id, egGenomicLoc$gene_id), c("start_location", "Chromosome")]
sym2eg[which(sym2eg$chr == "X"), "chr"] <- 20
sym2eg[which(sym2eg$chr == "Y"), "chr"] <- 21
sym2eg$chr <- as.integer(sym2eg$chr)
sym2eg <- sym2eg[rowSums(is.na(sym2eg)) == 0, ]
genfacl <- list()
for(restype in c("uv", "mv", "mv.bigeff", "uv.bigeff", "f")[1:5]){#restype <- "uv"#
  genfacl[[restype]] <- list()
  for(dirc in dirv){#dirc <- dirv[1]#
    genfacl[[dirc]] <- list()
    if(restype %in% c("uv", "mv", "mv.bigeff", "uv.bigeff")){
      phvc <- ph.use
      procv <- Data_all$impc$procord
    }
    if(restype == "f") {
      phvc <- facnam
    }
    for(phc in phvc){#phc <- ph.use[1]#
      impc.idc <- unique(sapply(strsplit(rownames(sigl[[restype]][[dirc]])[which(sigl[[restype]][[dirc]][, phc])], spl = "_"), 
                                function(v) v[1]))
      impc.idc.meas <- unique(sapply(strsplit(rownames(sigl[[restype]][[dirc]])[which(!is.na(sigl[[restype]][[dirc]][, phc]))], spl = "_"), 
                                function(v) v[1]))
      hitsc <- unique(c(na.omit(genemap2[match(impc.idc, genemap2$genotype_id), "entrez"])))
      measc <- unique(c(na.omit(genemap2[match(impc.idc.meas, genemap2$genotype_id), "entrez"])))
      phnamc <- phmap[match(phc, phmap$ph), "nam"]
      genfacl[[restype]][[dirc]][[phc]]$hits <- hitsc
      genfacl[[restype]][[dirc]][[phc]]$measured <- measc
      if(restype %in% c("uv", "mv", "mv.bigeff", "uv.bigeff")){
        procc <- phmap[match(phc, phmap$ph), "procnam"]
        genfacl[[restype]][[dirc]][[procc]]$hits <- unique(c(genfacl[[restype]][[dirc]][[procc]]$hits, hitsc))
        genfacl[[restype]][[dirc]][[procc]]$measured <- unique(c(genfacl[[restype]][[dirc]][[procc]]$measured, measc))
      }
    }
  }
}


#############################################
# Gene enrichment analysis
go_file_use <- gsub("go", paste0("go_sdmult_th_", sd_mult_big_eff_thresh), control$file.go.results)
if (rerun.GO.analysis | !file.exists(go_file_use)) {
  cl <- parallel::makeCluster(20)
  doParallel::registerDoParallel(cl = cl)
  featvc <- c(ph.use, facnam)
  parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  outl <- foreach(phc = featvc, .verbose = T, .errorhandling = "pass", .packages = c("GOfuncR", "Mus.musculus")) %dopar% {
    tryer <- try({
      allresl <- list()
      if (phc %in% c(ph.use, Data_all$impc$procord)) {
        restypev <- c("uv", "mv", "uv.bigeff", "mv.bigeff")[c(3, 4)]
      }
      if(phc %in% facnam) {
        restypev <- "f"
      }
      for(restype in restypev){
          myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
          allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
          myInterestingGenes.symbol <- c(na.omit(sym2eg[match(myInterestingGenes, sym2eg$gene_id), "symbol"]))
          allPossibleGenes.symbol <- c(na.omit(sym2eg[match(allPossibleGenes, sym2eg$gene_id), "symbol"]))
          geneList.bin <- as.integer(allPossibleGenes.symbol %in% myInterestingGenes.symbol)
          names(geneList.bin) <- geneList.bin
          input_hyper <- data.frame(allPossibleGenes.symbol, is_candidate = geneList.bin)
          set.seed(1)
          res_hyper <- GOfuncR::go_enrich(input_hyper, n_randset = nperm, organismDb = 'Mus.musculus', 
                                 domains = c('cellular_component', 'biological_process', 'molecular_function')[2])
          stats <- res_hyper[[1]]
          allRes <- stats[stats$FWER_overrep < fwer.th, ]
          allresl[[restype]] <- allRes
      }
    })
    if(!inherits(tryer, "try-error")) {
      return(allresl)
    } else {
      return(list(ph = phc, error = tryer))
    }
  }
  parallel::stopCluster(cl)
  featvc_full_ph_nam <- Data_all$impc$phmap[match(featvc, Data_all$impc$phmap$ph), "nam"]
  featvc_full_ph_nam[is.na(featvc_full_ph_nam)] <- featvc[is.na(featvc_full_ph_nam)]
  names(outl) <- featvc_full_ph_nam
  for (j in 1:length(outl)) {
    if (NROW(outl[[j]]$uv.bigeff) > 0) {
      outl[[j]]$uv.bigeff$ph <- names(outl)[j]
    }
    if (NROW(outl[[j]]$mv.bigeff) > 0) {
      outl[[j]]$mv.bigeff$ph <- names(outl)[j]
    }
  }
  
  str(outl[[1]])
  gonamv <- unlist(lapply(outl, function(x) c(x$uv$node_name, x$mv$node_name, x$mv.bigeff$node_name, x$uv.bigeff$node_name)))
  goidv <- unlist(lapply(outl, function(x) c(x$uv$node_id, x$mv$node_id, x$mv.bigeff$node_id, x$uv.bigeff$node_id)))
  gonamidmap <- unique(data.frame(nam = gonamv, id = goidv))
  go.gene.dat <- GOfuncR::get_anno_genes(go_ids = goidv, database = 'Mus.musculus', genes = NULL, annotations = NULL,
                                term_df = NULL, graph_path_df = NULL, godir = NULL)
  go.gene.dat$entrez <- sym2eg[match(go.gene.dat$gene, sym2eg$symbol), "gene_id"]
  go.gene.dat$go_name <- gonamidmap[match(go.gene.dat$go_id, gonamidmap$id), "nam"]
  outl.sub <- outl
  
  save(go.gene.dat, outl.sub, outl, file = go_file_use)
} else {
  load(file = go_file_use)
}


# library(help = "Mus.musculus")
genfacl[[restype]][[dirc]][[phc]]$measured

names(outl)
length(unique(go.gene.dat$go_id))

str(go.gene.dat)
names(outl.sub)
sort(phmap$nam)

co_enrich <- list(list(impc = "Locomotor activity", go = "locomotory behavior"),
                  list(impc = "Bone Area", go = "locomotory behavior"),
                  list(impc = "Click-evoked ABR threshold", go = "sensory perception of sound"),
                  list(impc = "% Pre-pulse inhibition - PPI4", go = "sensory perception of sound"))






control$dropbox_small_table_dir <- file.path(control$dropbox_table_dir, "little_go_tables")
dir.create(control$dropbox_small_table_dir, showWarnings = FALSE)

# go_inlab <- "GO in"
# go_outlab <- "GO out"
# impc_inlab <- "IMPC in"
# impc_outlab <- "IMPC out"
# go_inlab <- "In gene set"
# go_outlab <- "Not in"
# impc_inlab <- "In"
# impc_outlab <- "Not in"
# inoutvec <- c("In", "Out")
inoutvec <- c("Yes", "No")
go_inlab <- inoutvec[1]
go_outlab <- inoutvec[2]
impc_inlab <- inoutvec[1]
impc_outlab <- inoutvec[2]
for (pair_num in 1:length(co_enrich)) {#pair_num <- 2
  impc_ph_name_curr <- co_enrich[[pair_num]]$impc
  impc_ph_code_curr <- phmap[match(impc_ph_name_curr, phmap$nam), "ph"]
  go_term_name_curr <- co_enrich[[pair_num]]$go
  go_term_code_curr <- go.gene.dat[match(go_term_name_curr, go.gene.dat$go_name), "go_id"]
  
  impc_ph_mv_hits <- genfacl$mv.bigeff$bo[[impc_ph_code_curr]]$hits
  impc_ph_uv_hits <- genfacl$uv.bigeff$bo[[impc_ph_code_curr]]$hits
  impc_ph_measured_mv <- genfacl$mv.bigeff$bo[[impc_ph_code_curr]]$measured
  impc_ph_measured_uv <- genfacl$uv.bigeff$bo[[impc_ph_code_curr]]$measured
  go_term_tab <- go.gene.dat[go.gene.dat$go_name == go_term_name_curr, ]
  # go_term_measured <- intersect(impc_ph_measured, go_term_tab$entrez)
  go_enrich_df_mv <- data.frame(all_measured = impc_ph_measured_mv, 
                                impc_hit = ifelse(impc_ph_measured_mv %in% impc_ph_mv_hits, 2, 1),
                                go_hit = ifelse(impc_ph_measured_mv %in% go_term_tab$entrez, 2, 1))
  go_enrich_df_uv <- data.frame(all_measured = impc_ph_measured_uv, 
                                impc_hit = ifelse(impc_ph_measured_uv %in% impc_ph_uv_hits, 2, 1),
                                go_hit = ifelse(impc_ph_measured_uv %in% go_term_tab$entrez, 2, 1))
  
  tab_mv <- table(go_enrich_df_mv$go_hit, go_enrich_df_mv$impc_hit)
  tab_uv <- table(go_enrich_df_uv$go_hit, go_enrich_df_uv$impc_hit)
  dimnames(tab_mv) <- dimnames(tab_uv) <- list(c(go_outlab, go_inlab), c(impc_outlab, impc_inlab))
  print(co_enrich[[pair_num]])
  pval_uv <- fisher.test(tab_uv)$p.value
  pval_mv <- fisher.test(tab_mv)$p.value
  write.table(formatC(pval_uv, digits = 2, format = "g"), file = paste0(control$dropbox_small_table_dir, "/uv_pval_", pair_num, ".txt"), 
              col.names = F, row.names = F, quote = F)
  write.table(formatC(pval_mv, digits = 2, format = "g"), file = paste0(control$dropbox_small_table_dir, "/mv_pval_", pair_num, ".txt"), 
              col.names = F, row.names = F, quote = F)
  
  print(tab_mv)
  print(tab_uv)
  print(pval_mv)
  print(pval_uv)
  other_tab <- outl.sub[[impc_ph_name_curr]]$mv.bigeff
  print(other_tab[other_tab$node_name == go_term_name_curr, ])
  tabout_mv <- xtable::print.xtable(xtable::xtable(tab_mv, caption = "MV analysis"), floating = FALSE, caption.placement = 'top') 
  cat(tabout_mv, file = paste0(control$dropbox_small_table_dir, "/go_tab_mv_", pair_num, ".txt"))
  print(other_tab[other_tab$node_name == go_term_name_curr, ])
  tabout_uv <- xtable::print.xtable(xtable::xtable(tab_uv, caption = "UV analysis"), floating = FALSE, caption.placement = 'top') 
  cat(tabout_uv, file = paste0(control$dropbox_small_table_dir, "/go_tab_uv_", pair_num, ".txt"))
  # cat(tabout, file = paste0(control$dropbox_small_table_dir, "/go_tab_", impc_ph_code_curr, "_", gsub(":", "_", go_term_code_curr), ".txt"))
  impc_ph_name_curr_out <- gsub("%", "\\\\%", impc_ph_name_curr)
  write.table(impc_ph_name_curr_out, file = paste0(control$dropbox_small_table_dir, "/impc_ph_name_", pair_num, ".txt"), 
              col.names = F, row.names = F, quote = F)
  write.table(go_term_name_curr, file = paste0(control$dropbox_small_table_dir, "/go_term_name_", pair_num, ".txt"), 
              col.names = F, row.names = F, quote = F)
  
}#



gonamv <- unlist(lapply(outl, function(x) c(x$uv$node_name, x$mv$node_name, x$mv.bigeff$node_name, x$uv.bigeff$node_name)))
goidv <- unlist(lapply(outl, function(x) c(x$uv$node_id, x$mv$node_id, x$mv.bigeff$node_id, x$uv.bigeff$node_id)))
gonamidmap <- unique(data.frame(nam = gonamv, id = goidv))
go.gene.dat <- GOfuncR::get_anno_genes(go_ids = goidv, database = 'Mus.musculus', genes = NULL, annotations = NULL,
                                       term_df = NULL, graph_path_df = NULL, godir = NULL)
go.gene.dat$entrez <- sym2eg[match(go.gene.dat$gene, sym2eg$symbol), "gene_id"]
go.gene.dat$go_name <- gonamidmap[match(go.gene.dat$go_id, gonamidmap$id), "nam"]


##########################################
# Numbers for text
uniquenams <- unique(go.gene.dat$go_name)
ngo <- length(uniquenams)
ph.use.nam <- phmap[match(ph.use, phmap$ph), "nam"]
uvmat <- uvmat.nup <- uvmat.ndo <- mvmat <- mvmat.nup <- mvmat.ndo <- matrix(NA, ngo, nph, dimnames = list(uniquenams, ph.use.nam))
for(phc in ph.use.nam){# phc <- ph.use.nam[1]#
  for(restype in c("mv.bigeff", "uv.bigeff")){# restype <- "mv.bigeff"#
    gotabc <- outl[[phc]][[restype]]
    if(NROW(gotabc) > 0) {
      for(j in 1:nrow(gotabc)){
        idc <- gotabc[j, "node_id"]
        namc <- gotabc[j, "node_name"]
        genesc <- go.gene.dat[which(go.gene.dat$go_id == idc), "entrez"]
        ph.id <- phmap[match(phc, phmap$nam), "ph"]
        ndo <- sum(genfacl[[restype]]$do[[ph.id]]$hits %in% genesc)
        nup <- sum(genfacl[[restype]]$up[[ph.id]]$hits %in% genesc)
        propup <- (nup / (nup + ndo) - .5) * 2
        if(restype == "mv.bigeff"){
          mvmat[namc, phc] <- propup
          mvmat.nup[namc, phc] <- nup
          mvmat.ndo[namc, phc] <- ndo
        }
        if(restype == "uv.bigeff"){
          uvmat[namc, phc] <- propup
          uvmat.nup[namc, phc] <- nup
          uvmat.ndo[namc, phc] <- ndo
        }
      }
    }
  }
}

n.go.annot.mv <- sum(!is.na(mvmat))
n.go.annot.uv <- sum(!is.na(uvmat))
n.mv.greater <- sum(colSums(!is.na(mvmat)) > colSums(!is.na(uvmat)))
n.uv.greater <- sum(colSums(!is.na(uvmat)) > colSums(!is.na(mvmat)))
prop.mv.greater <- mean(colSums(!is.na(mvmat)) > colSums(!is.na(uvmat)))
prop.uv.greater <- mean(colSums(!is.na(uvmat)) > colSums(!is.na(mvmat)))
save.num <- c("n.go.annot.mv", "n.go.annot.uv", "n.mv.greater", "n.uv.greater")
for(numc in save.num){
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
}
save.prop <- c("prop.mv.greater", "prop.uv.greater")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 0, format = "f"), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
                                                col.names = F, row.names = F, quote = F)

##########################################
# Heatmap of GO terms for each phenotype perturbing gene set
ngo.plth <- 3#5
nph.plth <- 3#3
clustlinkage <- c("ward.D", "ward.D2", "single", "complete", "average")[1]
dist_method <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")[3]
mvmat.pl <- mvmat
uvmat.pl <- uvmat
# mvmat.pl <- mvmat.pl[rowSums(!is.na(mvmat.pl)) >= ngo.plth, ] 
for(j in 1:10){
  mvmat.pl <- mvmat.pl[rowSums(!is.na(mvmat.pl)) >= ngo.plth, ] 
  mvmat.pl <- mvmat.pl[, colSums(!is.na(mvmat.pl)) >= nph.plth] 
}
mvmat.pl.nona <- mvmat.pl
mvmat.pl.nona[!is.na(mvmat.pl.nona)] <- 1
mvmat.pl.nona[is.na(mvmat.pl.nona)] <- 0
uvmat.pl.nona <- uvmat.pl
uvmat.pl.nona[!is.na(uvmat.pl.nona)] <- 1
uvmat.pl.nona[is.na(uvmat.pl.nona)] <- 0
# for (j in 1:10) {
  mvmat.pl <- mvmat.pl[hclust(dist(mvmat.pl.nona, method = dist_method), method = clustlinkage)$order, ]
  mvmat.pl <- mvmat.pl[, hclust(dist(t(mvmat.pl.nona), method = dist_method), method = clustlinkage)$order]
# }
uvmat.pl <- uvmat[rownames(mvmat.pl), colnames(mvmat.pl)]
mode(mvmat.pl) <- mode(uvmat.pl) <- "numeric"
mvmat.cl <- mvmat.pl
uvmat.pl.nona <- uvmat
uvmat.pl.nona[!is.na(uvmat.pl.nona)] <- 1
uvmat.pl.nona[is.na(uvmat.pl.nona)] <- 0
# for (j in 1:10) {
  uvmat.cl <- uvmat[hclust(dist(uvmat.pl.nona, method = dist_method), method = clustlinkage)$order, ]
  uvmat.cl <- uvmat.cl[, hclust(dist(t(uvmat.pl.nona), method = dist_method), method = clustlinkage)$order]
# }
mv.go.ph <- lapply(rownames(mvmat.cl), function(x) colnames(mvmat.cl)[which(!is.na(mvmat.cl[x, ]))])
uv.go.ph <- lapply(rownames(uvmat.cl), function(x) colnames(uvmat.cl)[which(!is.na(uvmat.cl[x, ]))])
names(mv.go.ph) <- rownames(mvmat.cl)
names(uv.go.ph) <- rownames(uvmat.cl)
mv.go.proc <- lapply(mv.go.ph, function(x) unique(phmap[match(x, phmap$nam), "procnam"]))
uv.go.proc <- lapply(uv.go.ph, function(x) unique(phmap[match(x, phmap$nam), "procnam"]))

subplot.go.terms <- c("equilibrioception",  "nervous system process", 
                      "sensory perception",  "inner ear receptor cell development", 
                      "sensory perception of sound", 
                      "negative regulation of metabolic process", 
                      "developmental growth", "regulation of cell development", 
                      "regulation of neurogenesis", 
                      "regulation of steroid metabolic process", 
                      "multicellular organism growth", "regulation of cell differentiation", 
                      "growth hormone secretion", 
                      "response to nutrient levels", 
                      "leptin-mediated signaling pathway", 
                      "regulation of gluconeogenesis", "physiological cardiac muscle hypertrophy", 
                      "regulation of lipid metabolic process",
                      "muscle system process", 
                      "regulation of cell projection organization", 
                      "positive regulation of circadian sleep/wake cycle, non-REM sleep", 
                      "positive regulation of developmental growth", 
                      "positive regulation of multicellular organism growth", 
                      "regulation of signaling", 
                      "cell population proliferation", "hematopoietic or lymphoid organ development", 
                      "immune system development",  
                      "B cell homeostasis", "lymphocyte homeostasis", "homeostasis of number of cells",  
                      "circadian rhythm", "growth", "negative regulation of transport", 
                      "signal release", "regulation of biological quality", "regulation of hormone levels", 
                      "lipid biosynthetic process", "response to leptin", 
                      "regulation of circadian rhythm", 
                      "system development", "cell differentiation", 
                      "synaptic signaling", 
                      "nervous system development", 
                      "locomotory behavior", "behavior",  "neuron development", 
                      "neuron differentiation")

subplot.go.terms <- subplot.go.terms[subplot.go.terms %in% rownames(mvmat.pl)]
allout <- data.frame()
for(j in 1:length(outl.sub))
  allout <- rbind(allout, outl.sub[[j]]$mv.bigeff)
allout <- allout[order(allout$raw_p_overrep), ]
allsub <- allout[match(subplot.go.terms, allout$node_name), ]
str(allout)
gotab.paper <- allsub[, c("node_name", "ph", "raw_p_overrep")]
gotab.paper$proc <- phmap[match(gotab.paper$ph, phmap$nam), "procnam"]
gotab.paper$nup <- mvmat.nup[cbind(gotab.paper$node_name, gotab.paper$ph)]
gotab.paper$ndo <- mvmat.ndo[cbind(gotab.paper$node_name, gotab.paper$ph)]
gotab.paper$raw_p_overrep <- signif(gotab.paper$raw_p_overrep, 2)
for(j in 1:ncol(gotab.paper)){
  gotab.paper[, j] <- gsub("positive regulation of", "up", gotab.paper[, j])
  gotab.paper[, j] <- gsub("negative regulation of", "down", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Acoustic Startle and Pre-pulse Inhibition \\(PPI\\)", "Acoustic Startle", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Intraperitoneal glucose tolerance test \\(IPGTT\\)", "IPGTT", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Body Composition \\(DEXA lean/fat\\)", "DEXA", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Intraperitoneal glucose tolerance test \\(IPGTT\\)", "IPGTT", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Intraperitoneal glucose tolerance test \\(IPGTT\\)", "IPGTT", gotab.paper[, j])
  nchmax <- 40
  if(!is.numeric(mode(gotab.paper[, j])))
    gotab.paper[, j] <- paste0(substr(gotab.paper[, j], 1, nchmax), ifelse(nchar(gotab.paper[, j]) > nchmax, "...", ""))
}
gotab.paper


ph_unique <- sort(unique(allout$ph))
go_unique <- sort(unique(allout$node_name))
dput(names(allout))
fields_out <- data.frame(old = c("node_name", "ph", "raw_p_overrep"), 
                         new = c("GO gene set", "IMPC gene set", "Co-enrich p"))

n_go_pval_tables <- 8
go_terms_tables <- data.frame(set = c("sensory perception of sound", "locomotory behavior", 
                                      "regulation of lipid biosynthetic process", "brain development", 
                                      "chemical synaptic transmission", "growth", "circulatory system development",
                                      "anatomical structure development",
                                      "Fat/Body weight", "Insulin"),
                              impcgo = c("go", "go", "go", "go", "go",  "go", "go", "go", "impc", "impc"))[1:n_go_pval_tables, ]
go_terms_tables <- go_terms_tables[order(match(go_terms_tables$set, rev(rownames(mvmat.pl)))), ]
for (j in 1:nrow(go_terms_tables)) {
  set_curr <- go_terms_tables$set[j]
  if (go_terms_tables$impcgo[j] == "impc") {
    fields_out_use <- fields_out[fields_out$new != "IMPC gene set", ]
    tab_curr <- allout[allout$ph == set_curr, fields_out_use$old]
  }
  if (go_terms_tables$impcgo[j] == "go") {
    fields_out_use <- fields_out[fields_out$new != "GO gene set", ]
    tab_curr <- allout[allout$node_name == set_curr, fields_out_use$old]
  }
  # max_nchar <- 35
  # tab_curr[, 1] <- paste0(substr(tab_curr[, 1], 1, max_nchar), ifelse(nchar(tab_curr[, 1]) > max_nchar, "...", ""))
  # tab_curr
  name_out <- paste0("(", letters[j], ") ", toupper(go_terms_tables$impcgo[j]), " gene set \\emph{", set_curr, "}")
  name_out_short <- set_curr
  tab_curr$raw_p_overrep <- formatC(tab_curr$raw_p_overrep, digits = 2, format = "g")
  names(tab_curr) <- fields_out_use$new
  # tab_curr[, 1] <- paste0("\\textit{", tab_curr[, 1], "}")
  tabout <- xtable::print.xtable(xtable::xtable(tab_curr, caption = "", align = c("l", "p{5cm}", "r")), floating = FALSE, 
                                 caption.placement = 'top', include.rownames = FALSE) 
  cat(tabout, file = paste0(control$dropbox_small_table_dir, "/go_impc_enrich_", j, ".txt"))
  
  # name_out <- set_curr
  write.table(name_out, file = paste0(control$dropbox_small_table_dir, "/go_impc_enrich_name_", j, ".txt"), 
              col.names = F, row.names = F, quote = F)
  write.table(name_out_short, file = paste0(control$dropbox_small_table_dir, "/go_impc_enrich_name_short_", j, ".txt"), 
              col.names = F, row.names = F, quote = F)
}

phmap <- Data_all$impc$phmap
genemap <- Data_all$impc$genemap
y1 <- resl.out$eb$signsig[, phmap[phmap$nam == "Insulin", "ph"]]
y2 <- resl.out$eb$signsig[, phmap[phmap$nam == "Fat/Body weight", "ph"]]
table(y1, y2)

procord
impc_gene_ids <- sapply(strsplit(names(which(y1 == -1 & y2 == 1)), split = "_"), function(x) x[[1]])
genemap[match(impc_gene_ids, genemap$genotype_id), ]
str(genemap)
impc_gene_ids %in% genemap$genotype_id
phmap[phmap$nam == "Fat mass", "ph"]
str(resl.out$eb$signsig)


for(go.plot.type in c("all", "sub")){
  for(mod.type in c("mv", "uv")){
    go.use <- switch(go.plot.type, all = rownames(mvmat.pl), sub = subplot.go.terms)
    fnamc <- paste0(mod.type, "_", go.plot.type, "_go_heatmap_th_", sd_mult_big_eff_thresh, ".jpg")
    zpl <- switch(mod.type, mv = mvmat.pl[go.use, ], uv = uvmat.pl[go.use, ])
    jpeg(filename = paste(control$figure_dir, "/", fnamc, sep = ""), 
         width = 8, height = ifelse(go.plot.type == "all", 12, 6), 
         units = "in", res = 500)
    par(oma = c(12, 13, 0, 3), mar = c(0, 0, 0, 0))
    cexlab <- .6
    image(x = 1:ncol(zpl), y = 1:nrow(zpl), z = t(zpl), xlab = "", ylab = "", xaxt = "n", yaxt = "n", col = control$heat_col_palette, zlim = c(-1, 1))
    ncut <- 48
    axis(side = 2, labels = NA, at = 1:nrow(zpl), las = 2, cex.axis = cexlab)
    # axis(side = 3, labels = c("(c)", "(d)"), at = match(go_terms_tables$set[3:4], colnames(zpl)), las = 1, cex.axis = cexlab)
    axis(side = 4, labels = paste0("(", letters[1:n_go_pval_tables], ")"), 
         at = match(go_terms_tables$set, rownames(zpl)), las = 1, cex.axis = cexlab)
    labcut <- paste0(substr(rownames(zpl), 1, ncut), ifelse(nchar(rownames(zpl)) > ncut, "...", ""))
    mtext(side = 2, line = 1, text = labcut, at = 1:nrow(zpl), las = 2, cex = cexlab)
    axis(side = 1, labels = NA, at = 1:ncol(zpl), las = 2, cex.axis = cexlab)
    labcut <- paste0(substr(colnames(zpl), 1, ncut), ifelse(nchar(colnames(zpl)) > ncut, "...", ""))
    lincol <- "grey"
    abline(v = 1:ncol(zpl) - .5, col = lincol)
    abline(v = 1:ncol(zpl) - 1.5, col = lincol)
    abline(h = 1:nrow(zpl) - .5, col = lincol)
    mtext(side = 4, line = 2, text = "Rows labeled (a-h) are shown in Figure 9", cex = .8)
    phmap.sub <- phmap[phmap$nam %in% colnames(zpl), ]
    proctab <- sort(table(phmap.sub$procnam), decreasing = T)
    minperprocname <- 0
    names.leg <- names(proctab)[proctab > minperprocname]
    legeps <- .2
    textcoltab <- data.frame(phnam = colnames(zpl))
    textcoltab$procnam <- phmap[match(textcoltab$phnam, phmap$nam), "procnam"]
    colv <- control$heat_col_palette[(1:(length(names.leg) + 1)) / (length(names.leg) + 1) * 1000]
    textcoltab$col <- colv[match(textcoltab$procnam, names.leg)]
    textcoltab$col[is.na(textcoltab$col)] <- colv[length(names.leg) + 1]
    mtext(side = 1, line = 1, text = labcut, at = 1:ncol(zpl), las = 2, cex = cexlab, col = textcoltab[match(colnames(zpl), textcoltab$phnam), "col"])
    par(xpd = NA)
    legeps = .02
    if(minperprocname > 0){
      legend(x = -nrow(zpl) * legeps, y = -ncol(zpl) * legeps, xjust = 1, text.col = colv, legend = c(names.leg, "Other"), cex = .6)
    } else {
      legend(x = -nrow(zpl) * legeps, y = -ncol(zpl) * legeps, xjust = 1, text.col = colv[1:(length(colv) - 1)], title.col = 1,
             legend = names.leg, cex = .6, title = "Procedure color-code for phenotypes")
    }
    par(xpd = F)
    dev.off()
    file.copy(from = paste(control$figure_dir, "/", fnamc, sep = ""),
              to = paste(control$dropbox_figure_dir, "/", fnamc, sep = ""), overwrite = TRUE)
  }
}


###################################################################
# Table of GO term enrichment comparing UV and MV models
uvrows <- sapply(outl, function(x) NROW(x$uv.bigeff))
mvrows <- sapply(outl, function(x) NROW(x$mv.bigeff))
tabuvmvcomp <- table(uvrows, mvrows)
# binl <- list(0, 1, 2, 3:5, 6:10, 11:20, 20:50)
binl <- list(0, 1:5, 6:10, 11:20, 20:50)
binnam <- sapply(binl, function(x) if(length(x) == 1){ x} else {paste(min(x), max(x), sep = "-")})
nb <- length(binl)
tabcomp <- matrix(NA, nb, nb, dimnames = list(binnam, binnam))
for(j in 1:nb){
  for(k in 1:nb){
    tabcomp[j, k] <- sum(tabuvmvcomp[rownames(tabuvmvcomp) %in% as.character(binl[[j]]), colnames(tabuvmvcomp) %in% as.character(binl[[k]])])
  }
}
tabcomp
library(xtable)
tabout <- print(xtable(tabcomp, label = "tab:uvmvgocounts", 
                       caption = "Number of co-enriched  GO term gene sets for each IMPC gene set for the UV (left) and MV (top) models"),
                caption.placement = "top")
cat(tabout, file = paste(control$dropbox_table_dir, "/mvugotab_th_", sd_mult_big_eff_thresh, ".txt", sep = ""))


###################################################################
# Table of GO term enrichment comparing UV and MV models


str(uvmat)
uvcols <- rowSums(!is.na(uvmat))#sapply(outl, function(x) NCOL(x$uv.bigeff))
mvcols <- rowSums(!is.na(mvmat))#sapply(outl, function(x) NROW(x$mv.bigeff))
tabuvmvcomp <- table(uvcols, mvcols)
# binl <- list(0, 1, 2, 3:5, 6:10, 11:20, 20:50)
binl <- list(0, 1:5, 6:10, 11:20, 20:50)
binnam <- sapply(binl, function(x) if(length(x) == 1){ x} else {paste(min(x), max(x), sep = "-")})
nb <- length(binl)
tabcomp <- matrix(NA, nb, nb, dimnames = list(binnam, binnam))
for(j in 1:nb){
  for(k in 1:nb){
    tabcomp[j, k] <- sum(tabuvmvcomp[rownames(tabuvmvcomp) %in% as.character(binl[[j]]), colnames(tabuvmvcomp) %in% as.character(binl[[k]])])
  }
}
tabcomp
library(xtable)
tabout <- print(xtable(tabcomp, label = "tab:uvmv_enriched_counts_num_impc_per_go", 
                       caption = "Number of co-enriched IMPC gene sets for each GO term gene set for the UV (left) and MV (top) models"),
                caption.placement = "top")
cat(tabout, file = paste(control$dropbox_table_dir, "/uvmv_enriched_counts_num_impc_per_go_th_", sd_mult_big_eff_thresh, ".txt", sep = ""))







# 
# ############################################################
# # Just checking KEGG for enrichment
# # keggtab <- toTable(org.Mm.egPATH)
# keggtab <- toTable(org.Mm.egGO)
# keggtab2 <- keggtab[keggtab$gene_id %in% genemap2$entrez, ]
# n.min.at.term <- 1
# str(keggtab2)
# idc <- c("go_id", "path_id")[1]
# # keggtabc <- table(keggtab2$path_id)
# keggtabc <- table(keggtab2[, idc])
# keggkeep <- names(keggtabc)[keggtabc >= n.min.at.term]
# keggtab3 <- keggtab2[keggtab2[, idc] %in% keggkeep, ]
# keggentun <- unique(keggtab3$gene_id)
# keggpathun <- unique(keggtab3[, idc])
# entrez2kegg <- lapply(keggentun, function(x) keggtab3[keggtab3$gene_id == x, "path_id"])
# names(entrez2kegg) <- keggentun
# ensun <- unique(names(entrez2kegg))
# ng <- length(ensun)
# nk <- length(keggpathun)
# keggensmat <- matrix(0, ng, nk, dimnames = list(ensun, keggpathun))
# for(ensc in names(entrez2kegg)){
#   keggensmat[ensc, which(colnames(keggensmat) %in% entrez2kegg[[ensc]])] <- 1
# }
# perml <- list()
# methc <- "mv.bigeff"
# pvalmatl <- list()
# nperm <- 1
# for(i in 1:(nperm + 1)){
#   if(i == 1)
#     keggensmat.use <- keggensmat
#   if(i > 1)
#     keggensmat.use[] <- c(keggensmat[sample(1:nrow(keggensmat)), ])
#     # keggensmat.use[] <- c(keggensmat[, sample(1:ncol(keggensmat))])
#   pvalmat <- matrix(NA, nph, nk, dimnames = list(ph.use, keggpathun))
#   for(phc in ph.use){
#     print(phc)
#     ensmeasv <- genfacl[[methc]][[dirc]][[phc]]$measured
#     enshitv.use <- enshitv <- genfacl[[methc]][[dirc]][[phc]]$hits
#     # if(i == 1)
#     #   enshitv.use <- enshitv
#     # if(i > 1)
#     #   enshitv.use <- sample(ensmeasv, length(enshitv))
#     n.hit.in.kegg <- colSums(keggensmat.use[rownames(keggensmat.use) %in% enshitv.use, , drop = F])
#     n.hit.not.in.kegg <- colSums(1 - keggensmat.use[rownames(keggensmat.use) %in% enshitv.use, , drop = F])
#     n.measnothit.in.kegg <- colSums(keggensmat.use[rownames(keggensmat.use) %in% ensmeasv, , drop = F]) - n.hit.in.kegg
#     n.measnothit.not.in.kegg <- colSums(1 - keggensmat.use[rownames(keggensmat.use) %in% ensmeasv, , drop = F]) - n.hit.not.in.kegg
#     for(keggc in keggpathun){
#       fisher.tab <- matrix(c(n.hit.in.kegg[keggc], n.hit.not.in.kegg[keggc], n.meas.in.kegg[keggc], n.meas.not.in.kegg[keggc]), 2, 2)
#       fishout <- fisher.test(fisher.tab)
#       pvalmat[phc, keggc] <- fishout$p.value
#     }
#     pvalmatl[[i]] <- pvalmat
#   }
# }  
# 

# #############################################
# # Compare to mousenet network
# download.file(url = "https://www.inetbio.org/mousenet/dl_conv.php?f=MouseNetV2&s=0&t=0",
#     destfile = paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2.txt"))
# download.file(url = "https://www.inetbio.org/mousenet/dl_conv.php?f=MouseNetV2_GS&s=0&t=0",
#               destfile = paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_gs.txt"))
# download.file(url = "https://www.inetbio.org/mousenet/dl_conv.php?f=MM-LC&s=0&t=0",
#               destfile = paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_lc.txt"))
# download.file(url = "https://www.inetbio.org/mousenet/dl_conv.php?f=MM-CX&s=0&t=0",
#               destfile = paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_cx.txt"))
# download.file(url = "https://www.inetbio.org/mousenet/dl_conv.php?f=MM-GN&s=0&t=0",
#               destfile = paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_gn.txt"))
# download.file(url = "https://www.inetbio.org/mousenet/dl_conv.php?f=MM-PG&s=0&t=0",
#               destfile = paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_pg.txt"))
# mousenet <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2.txt"), sep = "\t")
# mousenet.gs <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_gs.txt"), sep = "\t")
# mousenet.lc <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_lc.txt"), sep = "\t")
# mousenet.cx <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_cx.txt"), sep = "\t")
# mousenet.gn <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_gn.txt"), sep = "\t")
# mousenet.pg <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_pg.txt"), sep = "\t")
# colnames(mousenet) <- colnames(mousenet.gs) <- colnames(mousenet.lc) <-
#   colnames(mousenet.cx) <- colnames(mousenet.gn) <- colnames(mousenet.pg) <- c("eg1", "eg2", "lr")
# 
# if("cl" %in% ls())
#   stopCluster(cl)
# cl <- makeCluster(12)
# registerDoParallel(cl = cl)
# dirc <- "bo"
# methc <- c("ph", "proc", "fac")[1]
# featvc <- switch(methc, ph = ph.use, proc = procord, fac = facnam)
# 
# outl <- foreach(phc = featvc, .verbose = T, .errorhandling = "pass") %dopar% {
#   # phc = ph.use[1]
#   pvall <- list()
#   for(restype in c("uv.bigeff", "mv.bigeff")){
#     myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
#     allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
#     all(myInterestingGenes %in% allPossibleGenes)
#     # myInterestingGenes.symbol <- c(na.omit(sym2eg[match(myInterestingGenes, sym2eg$gene_id), "symbol"]))
#     # allPossibleGenes.symbol <- c(na.omit(sym2eg[match(allPossibleGenes, sym2eg$gene_id), "symbol"]))
#     # mousenet.sub <- mousenet[(mousenet$eg1 %in% allPossibleGenes | mousenet$eg2 %in% allPossibleGenes), ]
#     mousenet.use <- list(mousenet.gs, mousenet.lc, mousenet.cx, mousenet.gn, mousenet.pg, mousenet)[[1]]
#     # if(any(!is.na(mousenet.use$lr))){
#     #   mousenet.sub <- mousenet.use[which((mousenet.use$eg1 %in% allPossibleGenes & mousenet.use$eg2 %in% allPossibleGenes) & mousenet.use$lr > 0), ]
#     # } else {
#     mousenet.sub <- mousenet.use[(mousenet.use$eg1 %in% allPossibleGenes & mousenet.use$eg2 %in% allPossibleGenes), ]
#     # }
#     str(mousenet.sub)
#     npair <- sum(mousenet.sub$eg1 %in% myInterestingGenes & mousenet.sub$eg2 %in% myInterestingGenes)
#     perm.npair <- c()
#     for(i in 1:1000){
#       myInterestingGenes.null <- sample(allPossibleGenes, length(myInterestingGenes))
#       perm.npair[i] <- sum(mousenet.sub$eg1 %in% myInterestingGenes.null & mousenet.sub$eg2 %in% myInterestingGenes.null)
#     }
#     pval.est <- mean(perm.npair >= npair)
#     pvall[[restype]] <- pval.est
#   }
#   return(list(pval = pvall, ph = phc))
# }
# stopCluster(cl)
# if(methc == "ph")
#   names(outl) <- phmap[match(ph.use, phmap$ph), "nam"]
# if(methc %in% c("proc", "fac"))
#   names(outl) <- featvc#phmap[match(ph.use, phmap$ph), "nam"]
# pvalv.uv <- sapply(outl, function(v) v$pval$uv.bigeff)
# pvalv.mv <- sapply(outl, function(v) v$pval$mv.bigeff)
# qvalv.uv <- p.adjust(pvalv.uv, meth = "BH")
# qvalv.mv <- p.adjust(pvalv.mv, meth = "BH")
# thc <- .05
# table(pvalv.uv < thc, pvalv.mv < thc)
# table(qvalv.uv < thc, qvalv.mv < thc)
# names(qvalv.mv)[qvalv.mv < thc]
# 
# mean(p.adjust(pvalv.mv, meth = "BH") < thc)
# 
# mean(pvalv.mv < thc)
# mean(pvalv.uv < thc)
# 
# 
# 
# ###################
# restype <- "mv.bigeff"
# methc <- c("ph", "proc", "fac")[2]
# featvc <- switch(methc, ph = ph.use, proc = procord, fac = facnam)
# pvalv <- c()
# for(phc in featvc){
#   myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
#   allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
#   geneList.bin <- as.integer(allPossibleGenes %in% myInterestingGenes)
#   maptab <- data.frame(eg = allPossibleGenes, sig = geneList.bin)
#   maptab[, c("loc", "chr")] <- sym2eg[match(allPossibleGenes, sym2eg$gene_id), c("loc", "chr")]
#   maptab <- maptab[order(maptab$chr, maptab$loc), ]
#   maptab$sym <- sym2eg[match(maptab$eg, sym2eg$gene_id), "symbol"]
#   chrun <- unique(maptab$chr)
#   ntot <- nrow(maptab)
#   nsig <- sum(maptab$sig)
#   # str(maptab)
#   table(maptab$chr)
#   vecsig <- table(maptab$chr[maptab$sig == 1])
#   vecsig <- vecsig[match(chrun, names(vecsig))]
#   vecsig[is.na(vecsig)] <- 0
#   vecnon <- table(maptab$chr[maptab$sig == 0])
#   vecnon <- vecnon[match(chrun, names(vecnon))]
#   vecnon[is.na(vecnon)] <- 0
#   
#   tab.chi <- cbind(vecsig, vecnon)
#   tab.chi <- tab.chi[rowSums(tab.chi) > 0, ]
#   chi.out <- chisq.test(tab.chi)
#   pvalv[phc] <- chi.out$p.value
# }
# 
# sort(pvalv)
# sort(p.adjust(pvalv, meth = "BH"))
# #################
# 
# 
# ############################################
# # Positional analysis
# methc <- c("ph", "proc", "fac")[2]
# featvc <- switch(methc, ph = ph.use, proc = Data_all$impc$procord, fac = facnam)
# dirc <- "bo"
# if("cl" %in% ls())
#   stopCluster(cl)
# cl <- makeCluster(25)
# registerDoParallel(cl = cl)
# par(mfrow = c(5, 5))
# phsub <- sample(ph.use, 25)
# nperm <- 500
# outl <- foreach(phc = featvc, .verbose = T, .errorhandling = "pass") %dopar% {#for(phc in phsub){#phc <- ph.use[20]#
#   library(topGO)
#   library(GOfuncR)
#   library(org.Mm.eg.db)
#   allresl <- list()
#   for(restype in c("uv.bigeff", "mv.bigeff")){
#     myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
#     allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
#     geneList.bin <- as.integer(allPossibleGenes %in% myInterestingGenes)
#     maptab <- data.frame(eg = allPossibleGenes, sig = geneList.bin)
#     maptab[, c("loc", "chr")] <- sym2eg[match(allPossibleGenes, sym2eg$gene_id), c("loc", "chr")]
#     maptab <- maptab[order(maptab$chr, maptab$loc), ]
#     maptab$sym <- sym2eg[match(maptab$eg, sym2eg$gene_id), "symbol"]
#     chrun <- unique(maptab$chr)
#     ntot <- nrow(maptab)
#     nsig <- sum(maptab$sig)
#     str(maptab)
#     table(maptab$chr)
#     chisq.test(cbind(table(maptab$chr[maptab$sig == 1]), table(maptab$chr[maptab$sig == 0])))
#     filthalfwidv <- c(3, 5, 7, 10, 15, 20, 30)[2]##[1:3]
#     filtshapel <- list()
#     for(filtwidc in filthalfwidv){
#       filtshapel[[as.character(filtwidc)]] <- dnorm(seq(-3, 3, len = 2 * filtwidc + 1))
#       filtshapel[[as.character(filtwidc)]] <- filtshapel[[as.character(filtwidc)]] / sum(filtshapel[[as.character(filtwidc)]])
#     }
#     
#     
#     
#     overlap.inds <- c()
#     act <- rep(0, ntot)
#     for(filtwidc in filthalfwidv){
#       filtc <- filter(maptab$sig, filter = filtshapel[[as.character(filtwidc)]])
#       # overlap.inds <- outer(match(chrun, maptab$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
#       # overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
#       # filtc[overlap.inds] <- NA
#       act <- pmax(act, filtc, na.rm = T)
#     }
#     permmat <- matrix(NA, ntot, nperm)
#     for(it in 1:nperm){
#       maptab.perm <- maptab
#       maptab.perm$sig <- sample(maptab.perm$sig)
#       permc <- rep(0, ntot)
#       for(filtwidc in filthalfwidv){
#         filtc <- filter(maptab.perm$sig, filter = filtshapel[[as.character(filtwidc)]])
#         # overlap.inds <- outer(match(chrun, maptab.perm$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
#         # overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
#         # filtc[overlap.inds] <- NA
#         permc <- pmax(permc, filtc, na.rm = T)
#       }
#       permmat[, it] <- permc#filter(maptab.perm$sig, filter = rep(1, filtwid))
#     }
#     allresl[[restype]] <-  list(permmat = permmat, act = act)
#   }
#   return(allresl)
# }
# names(outl) <- featvc
# stopCluster(cl)
# rm(cl)
# 
# str(outl)
# 
# 
# relmatl <- list()
# fwerc <- .05
# for(restype in c("uv.bigeff", "mv.bigeff")){
#   relmatl[[restype]] <- sapply(outl, function(x){
#                   thr.up <- quantile(apply(x[[restype]]$permmat, 2, function(v) max(v, na.rm = T)), 1 - fwerc, na.rm = T)
#                 return(x[[restype]]$act / thr.up)
#     })
# }
# 
# max.rel <- apply(relmatl$mv.bigeff, 2, function(v) max(v, na.rm = T))
# nplph <- 10
# plph <- colnames(relmatl$mv.bigeff)[order(-max.rel)[1:nplph]]
# graphics.off()
# matplot(relmatl$mv.bigeff[, plph], ty = "l")
# abline(h = 1)
# 
# phord <- ph.use[order(phmap[match(ph.use, phmap$ph), "procnam"])]
# graphics.off()
# par(mfrow = c(2, 1))
# image(relmatl$mv.bigeff[rowSums(is.na(relmatl$mv.bigeff)) == 0, phord])
# image(relmatl$mv.bigeff[rowSums(is.na(relmatl$mv.bigeff)) == 0, phord] > 1)
# 
# ############################################################################
# # Combined positional analysis (aggregating the number of hits across phenotypes within gene)
# # ssmat <- (abs(ebl$t) > ebl$th) * sign(ebl$t)
# ssmat <- resl.out$eb$signsig.bigeff[true.use, ]
# maptab <- data.frame(genzyg = rownames(ssmat))
# maptab$impc.id <- sapply(strsplit(maptab$genzyg, spl = "_"), function(v) v[1])
# maptab$zyg <- sapply(strsplit(maptab$genzyg, spl = "_"), function(v) v[2])
# maptab$eg <- genemap2[match(maptab$impc.id, genemap2$genotype_id), "entrez"]
# maptab[, c("loc", "chr")] <- sym2eg[match(maptab$eg, sym2eg$gene_id), c("loc", "chr")]
# maptab <- maptab[order(match(maptab$zyg, c(1, 0, 2))), ]
# genshun <- unique(maptab$eg)
# maptab <- maptab[match(genshun, maptab$eg), ]
# maptab <- maptab[order(maptab$chr, maptab$loc), ]
# chrun <- unique(maptab$chr)
# phord <- ph.use[order(phmap[match(ph.use, phmap$ph), "procnam"])]
# ssmat <- ssmat[maptab$genzyg, phord]
# ntot <- nrow(ssmat)
# agghits <- rowSums(abs(ssmat[, , drop = F]))
# filthalfwidv <- c(5)#10, 15, 20, 30, 50)[1]#[1:3]
# filtshapel <- list()
# for(filtwidc in filthalfwidv){
#   filtshapel[[as.character(filtwidc)]] <- dnorm(seq(-3, 3, len = 2 * filtwidc + 1))
#   filtshapel[[as.character(filtwidc)]] <- filtshapel[[as.character(filtwidc)]] / sum(filtshapel[[as.character(filtwidc)]])
# }
# overlap.inds <- c()
# act <- rep(0, ntot)
# for(filtwidc in filthalfwidv){
#   filtc <- filter(agghits, filter = filtshapel[[as.character(filtwidc)]])
#   overlap.inds <- outer(match(chrun, maptab$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
#   overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
#   filtc[overlap.inds] <- NA
#   act <- pmax(act, filtc, na.rm = T)
# }
# plot(act)
# permmat <- matrix(NA, ntot, nperm)
# nperm <- 500
# 
# for(it in 1:nperm){
#   agghits.perm <- sample(agghits)
#   permc <- rep(0, ntot)
#   for(filtwidc in filthalfwidv){
#     filtc <- filter(agghits.perm, filter =  filtshapel[[as.character(filtwidc)]])
#     overlap.inds <- outer(match(chrun, maptab$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
#     overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
#     filtc[overlap.inds] <- NA
#     permc <- pmax(permc, filtc, na.rm = T)
#   }
#   permmat[, it] <- permc#filter(maptab.perm$sig, filter = rep(1, filtwid))
# }
# thr.up <- quantile(apply(permmat, 2, function(v) max(v, na.rm = T)), 1 - fwerc, na.rm = T)
# plot(act / thr.up, ty = "l")
# # plot(act)
# 
# 
# 
# 
# 
