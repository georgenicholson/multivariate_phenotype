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

resl.out <- list()
resl.out$eb <- resll$perm[[control$mv_meth_nam_use]]
resl.out$uv <- resll$perm$uv
resl.out[[fac.meth]] <- resl.comp.fac[[control$mv_meth_nam_use]]
facnam <- colnames(resl.out[[fac.meth]]$mn)
resl.out[[fac.meth]]$t <- resl.out[[fac.meth]]$mn / resl.out[[fac.meth]]$sd
resl.out[[fac.meth]]$th <- resimp[, paste0(fac.meth, ".th.final")][1]

ph_high <- names(sort(abs(resl.comp[[control$mv_meth_nam_use]]$facs.varimax[, "fac_7"]), decreasing = T))

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
          length(allPossibleGenes)
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


outl <- outl[!grepl("fac", names(outl))]

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
  tabout_mv <- xtable::print.xtable(xtable::xtable(tab_mv, caption = "(b)     MV analysis"), floating = FALSE, caption.placement = 'top') 
  cat(tabout_mv, file = paste0(control$dropbox_small_table_dir, "/go_tab_mv_", pair_num, ".txt"))
  print(other_tab[other_tab$node_name == go_term_name_curr, ])
  tabout_uv <- xtable::print.xtable(xtable::xtable(tab_uv, caption = "(a)     UV analysis"), floating = FALSE, caption.placement = 'top') 
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

###########################################
# Check the # GO terms in BP subgraph for gene background set
gotab <- as.data.frame(org.Mm.eg.db::org.Mm.egGO)
allPossibleGenes <- genfacl$mv.bigeff$bo[[1]]$measured
n_hom_genes <- length(allPossibleGenes)
n_bp_go_terms <- length(unique(gotab$go_id[gotab$Ontology == "BP" & gotab$gene_id %in% allPossibleGenes]))

n.go.annot.mv <- sum(!is.na(mvmat))
n.go.annot.uv <- sum(!is.na(uvmat))
n.mv.greater <- sum(colSums(!is.na(mvmat)) > colSums(!is.na(uvmat)))
n.uv.greater <- sum(colSums(!is.na(uvmat)) > colSums(!is.na(mvmat)))
prop.mv.greater <- mean(colSums(!is.na(mvmat)) > colSums(!is.na(uvmat)))
prop.uv.greater <- mean(colSums(!is.na(uvmat)) > colSums(!is.na(mvmat)))
save.num <- c("n.go.annot.mv", "n.go.annot.uv", "n.mv.greater", "n.uv.greater", "n_bp_go_terms", "n_hom_genes")
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
save.num <- c("ngo.plth", "nph.plth")
for(numc in save.num){
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(control$dropbox_text_numbers_dir, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
}

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

red_blue_col_palette <- c(rgb(0, 0, 1, alpha = seq(1, 0, len = 500)), rgb(1, 0, 0, alpha = seq(0, 1, len = 500)))
for(go.plot.type in c("all", "sub")[1]){
  for(mod.type in c("mv", "uv")){
    go.use <- switch(go.plot.type, all = rownames(mvmat.pl), sub = subplot.go.terms)
    fnamc <- paste0(mod.type, "_", go.plot.type, "_go_heatmap_th_", sd_mult_big_eff_thresh, ".jpg")
    zpl <- switch(mod.type, mv = mvmat.pl[go.use, ], uv = uvmat.pl[go.use, ])
    jpeg(filename = paste(control$figure_dir, "/", fnamc, sep = ""), 
         width = 8, height = ifelse(go.plot.type == "all", 12, 6), 
         units = "in", res = 500)
    
    wid_eps <- .1
    hei_eps <- .035
    # layout(mat = matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2), widths = c(1 - wid_eps, wid_eps), c(hei_eps, 1 - hei_eps))
    layout(mat = matrix(c(2, 1), nrow = 2, ncol = 1), heights = c(hei_eps, 1 - hei_eps))
    par(oma = c(12, 13, 3, 3), mar = c(0, 0, 0, 0))
    cexlab <- .6
    image(x = 1:ncol(zpl), y = 1:nrow(zpl), z = t(zpl), xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
          col = red_blue_col_palette, zlim = c(-1, 1))
    ncut <- 48
    axis(side = 2, labels = NA, at = 1:nrow(zpl), las = 2, cex.axis = cexlab)
    # axis(side = 3, labels = c("(c)", "(d)"), at = match(go_terms_tables$set[3:4], colnames(zpl)), las = 1, cex.axis = cexlab)
    axis(side = 4, labels = paste0("(", letters[1:n_go_pval_tables], ")"), 
         at = match(go_terms_tables$set, rownames(zpl)), las = 1, cex.axis = .8)
    labcut <- paste0(substr(rownames(zpl), 1, ncut), ifelse(nchar(rownames(zpl)) > ncut, "...", ""))
    mtext(side = 2, line = 1, text = labcut, at = 1:nrow(zpl), las = 2, cex = cexlab)
    axis(side = 1, labels = NA, at = 1:ncol(zpl), las = 2, cex.axis = cexlab)
    labcut <- paste0(substr(colnames(zpl), 1, ncut), ifelse(nchar(colnames(zpl)) > ncut, "...", ""))
    lincol <- "light grey"
    abline(v = 1:ncol(zpl) - .5, col = lincol)
    abline(v = 1:ncol(zpl) - 1.5, col = lincol)
    abline(h = 1:nrow(zpl) - .5, col = lincol)
    for (i in 1:ncol(zpl)) {
      for (j in 1:nrow(zpl)) {
        if (!is.na(t(zpl)[i, j])) {
          polygon(x = i + c(-1, 1, 1, -1, -1) * .4, y = j + c(-1, -1, 1, 1, -1) * .4, lwd = 1.5, lend = 2)
        }
      }
    }
    mtext(side = 4, line = 2, text = "Rows labeled (a-h) are shown in Table 3", cex = .8)
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
             legend = names.leg, cex = .6, title = "Procedure colour-code for phenotypes")
    }
    par(xpd = F)
    par(mar = c(.5, 13, .5, 0), xpd = F)
    cexax <- .8
    image(z = (as.matrix(1:1000)), y = 1, x = seq(0, 100, len = 1000), col = red_blue_col_palette, xaxt = "n", yaxt = "n", 
          ylab = "")
    axis(side = 3, las = 2, cex.axis = cexax, las = 1)#, labels = c("< -1.0", -0.5, 0.0, 0.5, "> 1.0"), at = seq(0, 1, by = .25))
    mtext(side = 3, text = "% genes with phenotype perturbed upwards", line = 2, cex = cexax)
    
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
binl <- list(0, 1:5, 6:10, 11:20, 21:1000)
binnam <- sapply(binl, function(x) if(length(x) == 1){ x} else {paste(min(x), max(x), sep = "-")})
nb <- length(binl)
tabcomp <- matrix(NA, nb, nb, dimnames = list(binnam, binnam))
for(j in 1:nb){
  for(k in 1:nb){
    tabcomp[j, k] <- sum(tabuvmvcomp[rownames(tabuvmvcomp) %in% as.character(binl[[j]]), colnames(tabuvmvcomp) %in% as.character(binl[[k]])])
  }
}
rownames(tabcomp)[rownames(tabcomp) == "21-1000"] <- ">20"
colnames(tabcomp)[colnames(tabcomp) == "21-1000"] <- ">20"
tabcomp
tabcomp[] <- prettyNum(tabcomp, big.mark = ",")

sum(tabcomp)
library(xtable)
tabout <- print(xtable(tabcomp, label = "tab:uvmvgocounts", 
                       caption = "Number of co-enriched  GO term gene sets for each IMPC gene set for the UV (left) and MV (top) models"),
                caption.placement = "top", 
                floating = FALSE)
cat(tabout, file = paste(control$dropbox_table_dir, "/mvugotab_th_", sd_mult_big_eff_thresh, ".txt", sep = ""))


###################################################################
# Table of GO term enrichment comparing UV and MV models


str(uvmat)
uvcols <- rowSums(!is.na(uvmat))#sapply(outl, function(x) NCOL(x$uv.bigeff))
mvcols <- rowSums(!is.na(mvmat))#sapply(outl, function(x) NROW(x$mv.bigeff))
tabuvmvcomp <- table(uvcols, mvcols)
# binl <- list(0, 1, 2, 3:5, 6:10, 11:20, 20:50)
binl <- list(0, 1:5, 6:10, 11:20, 21:1000)
binnam <- sapply(binl, function(x) if(length(x) == 1){ x} else {paste(min(x), max(x), sep = "-")})
nb <- length(binl)
tabcomp <- matrix(NA, nb, nb, dimnames = list(binnam, binnam))
for(j in 1:nb){
  for(k in 1:nb){
    tabcomp[j, k] <- sum(tabuvmvcomp[rownames(tabuvmvcomp) %in% as.character(binl[[j]]), colnames(tabuvmvcomp) %in% as.character(binl[[k]])])
  }
}
tabcomp_with_zeroes <- tabcomp
tabcomp_with_zeroes["0", "0"] <- n_bp_go_terms - sum(tabcomp)
rownames(tabcomp_with_zeroes)[rownames(tabcomp_with_zeroes) == "21-1000"] <- ">20"
colnames(tabcomp_with_zeroes)[colnames(tabcomp_with_zeroes) == "21-1000"] <- ">20"
tabcomp_with_zeroes[] <- prettyNum(tabcomp_with_zeroes, big.mark = ",")
tabcomp_with_zeroes
library(xtable)
tabout <- print(xtable(tabcomp_with_zeroes, label = "tab:uvmv_enriched_counts_num_impc_per_go", 
                       caption = "Number of co-enriched IMPC gene sets for each GO term gene set for the UV (left) and MV (top) models"),
                caption.placement = "top", 
                floating = FALSE)
cat(tabout, file = paste(control$dropbox_table_dir, "/uvmv_enriched_counts_num_impc_per_go_th_", sd_mult_big_eff_thresh, ".txt", sep = ""))





