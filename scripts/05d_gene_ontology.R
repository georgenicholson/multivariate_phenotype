.libPaths("X:/R_libraries")
resimp <- readRDS(file = control$file.resimp)
resl.comp.fac <- readRDS(file = control$file_raw_factor_results)
phmap <- Data_all$impc$phmap
library("foreach")

ph.use <- Data_all$impc$phord
use.just.homs <- F   # To be want to do the analysis using just homozygotes?
fac.meth <- c("varimax", "promax")[1]
if(use.just.homs){
  true.use <- unique(resimp[resimp$line.type == "trueMut" & resimp$zyg == 1, "geno"])
} else {
  true.use <- unique(resimp[resimp$line.type == "trueMut", "geno"])
}
nph <- length(ph.use)
ng <- length(true.use)
dirv <- c("up", "do", "bo")


resl.out <- list()
resl.out$eb <- resll$perm[[control$mv_meth_nam_use]]
resl.out$uv <- resll$perm$uv
resl.out[[fac.meth]] <- resl.comp.fac[[control$mv_meth_nam_use]]
facnam <- colnames(resl.out[[fac.meth]]$mn)
resl.out[[fac.meth]]$t <- resl.out[[fac.meth]]$mn / resl.out[[fac.meth]]$sd

mn.th.big <- 1
resl.out$eb$signsig.bigeff <- (abs(resl.out$eb$t) > resl.out$eb$th & abs(resl.out$eb$mn) >= mn.th.big) * sign(resl.out$eb$t)
resl.out$uv$signsig.bigeff <- (abs(resl.out$uv$t) > resl.out$uv$th & abs(resl.out$uv$mn) >= mn.th.big) * sign(resl.out$uv$t)
resl.out[[fac.meth]]$signsig.bigeff <- (abs(resl.out[[fac.meth]]$t) > resl.out[[fac.meth]]$th & 
                                          abs(resl.out[[fac.meth]]$mn) >= mn.th.big) * sign(resl.out[[fac.meth]]$t)

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
sigl$f$up <- resl.out[[fac.meth]]$signsig[true.use, facnam] == 1
sigl$f$do <- resl.out[[fac.meth]]$signsig[true.use, facnam] == -1
sigl$f$bo <- resl.out[[fac.meth]]$signsig[true.use, facnam] != 0


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
for(restype in c("uv", "mv", "mv.bigeff", "uv.bigeff", "f")[1:4]){#restype <- "uv"#
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

dirc <- "bo"
nperm <- 1000
fwer.th <- .05
methc <- c("ph", "proc", "fac")[1]
run.GO.analysis <- T

#############################################
# Gene enrichment analysis
if (run.GO.analysis) {
  if("cl" %in% ls())
    parallel::stopCluster(cl)
  cl <- parallel::makeCluster(20)
  doParallel::registerDoParallel(cl = cl)
  featvc <- switch(methc, ph = ph.use, proc = procord, fac = facnam)
  # library(GOfuncR)
  parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  outl <- foreach(phc = featvc, .verbose = T, .errorhandling = "pass", .packages = "GOfuncR") %dopar% {
    tryer <- try({
      allresl <- list()
      if (phc %in% c(ph.use, Data_all$impc$procord)) {
        restypev <- c("uv", "mv", "uv.bigeff", "mv.bigeff")[c(3, 4)]
      }
      # if(phc %in% facnam)
      #   restypev <- "f"
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
      return(tryer)
    }
  }
  parallel::stopCluster(cl)
  if (methc == "ph") {
    names(outl) <- phmap[match(ph.use, phmap$ph), "nam"]
    for (j in 1:length(outl)) {
      if (NROW(outl[[j]]$uv.bigeff) > 0) {
        outl[[j]]$uv.bigeff$ph <- names(outl)[j]
      }
      if (NROW(outl[[j]]$mv.bigeff) > 0) {
        outl[[j]]$mv.bigeff$ph <- names(outl)[j]
      }
    }
  }
  if (methc %in% c("proc", "fac")) {
    names(outl) <- featvc
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
  save(go.gene.dat, outl.sub, outl, file = file.go.results)
} else {
  load(file = file.go.results)
}

# str(outl[[1]])
# str(gotermtab)
# 
# gotab <- toTable(org.Mm.egGO)
# gotab2 <- gotab[gotab$gene_id %in% genemap2$entrez, ]
# tabc <- table(gotab2$go_id)
# gokeep <- names(tabc)#[tabc > 0]
# gotab3 <- gotab2[gotab2$go_id %in% gokeep, ]
# entun <- unique(gotab3$gene_id)
# goun <- unique(gotab3$go_id)
# entrez2go <- lapply(entun, function(x) gotab3[gotab3$gene_id == x, "go_id"])
# go2entrez <- lapply(goun, function(x) gotab3[gotab3$go_id == x, "gene_id"])
# names(entrez2go) <- entun
# names(go2entrez) <- goun
# 
# str(gotab)
# 
# mean(uniqueids %in% gotab$go_id)

##########################################
# Numbers for text
uniquenams <- unique(go.gene.dat$go_name)
ngo <- length(uniquenams)
ph.use.nam <- phmap[match(ph.use, phmap$ph), "nam"]
uvmat <- uvmat.nup <- uvmat.ndo <- mvmat <- mvmat.nup <- mvmat.ndo <- matrix(NA, ngo, nph, dimnames = list(uniquenams, ph.use.nam))
for(phc in ph.use.nam){
  for(restype in c("mv.bigeff", "uv.bigeff")){
    gotabc <- outl[[phc]][[restype]]
    if(NROW(gotabc) == 0)
      next
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


n.go.annot.mv <- sum(!is.na(mvmat))
n.go.annot.uv <- sum(!is.na(uvmat))
n.mv.greater <- sum(colSums(!is.na(mvmat)) > colSums(!is.na(uvmat)))
n.uv.greater <- sum(colSums(!is.na(uvmat)) > colSums(!is.na(mvmat)))
prop.mv.greater <- mean(colSums(!is.na(mvmat)) > colSums(!is.na(uvmat)))
prop.uv.greater <- mean(colSums(!is.na(uvmat)) > colSums(!is.na(mvmat)))
save.num <- c("n.go.annot.mv", "n.go.annot.uv", "n.mv.greater", "n.uv.greater")
for(numc in save.num){
  write.table(prettyNum(eval(as.name(numc)), big.mark = ","), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
              col.names = F, row.names = F, quote = F)
}
save.prop <- c("prop.mv.greater", "prop.uv.greater")
for(numc in save.prop)
  write.table(formatC(100 * eval(as.name(numc)), digits = 0, format = "f"), file = paste(revision.text.numbers, "/", numc, ".txt", sep = ""), 
                                                col.names = F, row.names = F, quote = F)



# str(mvmat)
# phsigncheck <- c("Sodium", "Creatine kinase", "Fat mass", "Insulin", "Glycerol")[1]
# table(resl.out$eb$signsig[true.use, phmap[match(phsigncheck, phmap$nam), "ph"]])
# phmap[grep("Bone", phmap$nam), ]

# str(phmap)
##########################################
# Heatmap of GO terms for each phenotype perturbing gene set
ngo.plth <- 2
nph.plth <- 3
clustlinkage <- c("ward.D", "ward.D2", "single", "complete", "average")[1]
mvmat.pl <- mvmat
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

mvmat.pl <- mvmat.pl[hclust(dist(mvmat.pl.nona), method = clustlinkage)$order, ]
mvmat.pl <- mvmat.pl[, hclust(dist(t(mvmat.pl.nona)), method = clustlinkage)$order]


uvmat.pl <- uvmat[rownames(mvmat.pl), colnames(mvmat.pl)]
# uvmat.pl <- mvmat.pl
# for(phc in colnames(uvmat.pl))
#   uvmat.pl[, phc] <- rownames(uvmat.pl) %in% outl[[phc]]$uv.bigeff$node_name
mode(mvmat.pl) <- mode(uvmat.pl) <- "numeric"


mvmat.cl <- mvmat.pl

uvmat.pl.nona <- uvmat
uvmat.pl.nona[!is.na(uvmat.pl.nona)] <- 1
uvmat.pl.nona[is.na(uvmat.pl.nona)] <- 0
uvmat.cl <- uvmat[hclust(dist(uvmat.pl.nona), method = clustlinkage)$order, ]
uvmat.cl <- uvmat.cl[, hclust(dist(t(uvmat.pl.nona)), method = clustlinkage)$order]


mv.go.ph <- lapply(rownames(mvmat.cl), function(x) colnames(mvmat.cl)[which(!is.na(mvmat.cl[x, ]))])
uv.go.ph <- lapply(rownames(uvmat.cl), function(x) colnames(uvmat.cl)[which(!is.na(uvmat.cl[x, ]))])
names(mv.go.ph) <- rownames(mvmat.cl)
names(uv.go.ph) <- rownames(uvmat.cl)
mv.go.proc <- lapply(mv.go.ph, function(x) unique(phmap[match(x, phmap$nam), "procnam"]))
uv.go.proc <- lapply(uv.go.ph, function(x) unique(phmap[match(x, phmap$nam), "procnam"]))
mv.go.ph
mv.go.proc

str(mvmat.nup)



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

gotab.paper <- allsub[, c("node_name", "ph", "raw_p_overrep")]
gotab.paper$proc <- phmap[match(gotab.paper$ph, phmap$nam), "procnam"]
gotab.paper$nup <- mvmat.nup[cbind(gotab.paper$node_name, gotab.paper$ph)]
gotab.paper$ndo <- mvmat.ndo[cbind(gotab.paper$node_name, gotab.paper$ph)]
gotab.paper$raw_p_overrep <- signif(gotab.paper$raw_p_overrep, 2)
# gotab.paper <- gsub("positive regultion", "pos. reg.", gotab.paper)
for(j in 1:ncol(gotab.paper)){
  # gotab.paper[, j] <- gsub("positive regulation", "pos. reg.", gotab.paper[, j])
  gotab.paper[, j] <- gsub("positive regulation of", "up", gotab.paper[, j])
  gotab.paper[, j] <- gsub("negative regulation of", "down", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Acoustic Startle and Pre-pulse Inhibition \\(PPI\\)", "Acoustic Startle", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Intraperitoneal glucose tolerance test \\(IPGTT\\)", "IPGTT", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Body Composition \\(DEXA lean/fat\\)", "DEXA", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Intraperitoneal glucose tolerance test \\(IPGTT\\)", "IPGTT", gotab.paper[, j])
  gotab.paper[, j] <- gsub("Intraperitoneal glucose tolerance test \\(IPGTT\\)", "IPGTT", gotab.paper[, j])
  
  nchmax <- 40
  if(!is.numeric(mode(gotab.paper[, j])))
    # gotab.paper[, j] <- substr(gotab.paper[, j], 1, 30)  
    gotab.paper[, j] <- paste0(substr(gotab.paper[, j], 1, nchmax), ifelse(nchar(gotab.paper[, j]) > nchmax, "...", ""))
}
gotab.paper


unique(phmap$procnam)
labcut <- paste0(substr(rownames(zpl), 1, ncut), ifelse(nchar(rownames(zpl)) > ncut, "...", ""))

dput(names(allsub))

for(go.plot.type in c("all", "sub")){
  for(mod.type in c("mv", "uv")){
    go.use <- switch(go.plot.type, all = rownames(mvmat.pl), sub = subplot.go.terms)
    fnamc <- paste0(mod.type, "_", go.plot.type, "_go_heatmap.jpg")
  #   
  # for(fnamc in c("mv_go_heatmap.jpg", "uv_go_heatmap.jpg")){
    zpl <- switch(mod.type, mv = mvmat.pl[go.use, ], uv = uvmat.pl[go.use, ])
    # zpl <- zpl[, colMeans(is.na;(zpl)) < 1]
    # jpeg(paste(go.plot.dir, "/", fnamc, sep = ""), 8, 7, units = "in", res = 500)
    jpeg(paste(go.plot.dir, "/", fnamc, sep = ""), 8, ifelse(go.plot.type == "all", 12, 6), units = "in", res = 500)
    par(oma = c(11, 13, .5, 0.5), mar = c(0, 0, 0, 0))
    cexlab <- .6
    image(x = 1:ncol(zpl), y = 1:nrow(zpl), z = t(zpl), xlab = "", ylab = "", xaxt = "n", yaxt = "n", col = rain, zlim = c(-1, 1))
    ncut <- 48
    axis(side = 2, labels = NA, at = 1:nrow(zpl), las = 2, cex.axis = cexlab)
    labcut <- paste0(substr(rownames(zpl), 1, ncut), ifelse(nchar(rownames(zpl)) > ncut, "...", ""))
    # phlabcut <- paste0(substr(colnames(zpl), 1, ncut), ifelse(nchar(colnames(zpl)) > ncut, "...", ""))
    mtext(side = 2, line = 1, text = labcut, at = 1:nrow(zpl), las = 2, cex = cexlab)
    axis(side = 1, labels = NA, at = 1:ncol(zpl), las = 2, cex.axis = cexlab)
    labcut <- paste0(substr(colnames(zpl), 1, ncut), ifelse(nchar(colnames(zpl)) > ncut, "...", ""))
    # mtext(side = 3, line = 1, text = labcut, at = 1:ncol(zpl), las = 2, cex = cexlab)
    lincol <- "grey"
    abline(v = 1:ncol(zpl) - .5, col = lincol)
    abline(h = 1:nrow(zpl) - .5, col = lincol)
    phmap.sub <- phmap[phmap$nam %in% colnames(zpl), ]
    proctab <- sort(table(phmap.sub$procnam), decreasing = T)
    minperprocname <- 0
    names.leg <- names(proctab)[proctab > minperprocname]
    legeps <- .2
    textcoltab <- data.frame(phnam = colnames(zpl))
    textcoltab$procnam <- phmap[match(textcoltab$phnam, phmap$nam), "procnam"]
    colv <- rain[(1:(length(names.leg) + 1)) / (length(names.leg) + 1) * 1000]
    textcoltab$col <- colv[match(textcoltab$procnam, names.leg)]
    textcoltab$col[is.na(textcoltab$col)] <- colv[length(names.leg) + 1]
    mtext(side = 1, line = 1, text = labcut, at = 1:ncol(zpl), las = 2, cex = cexlab, col = textcoltab[match(colnames(zpl), textcoltab$phnam), "col"])
    par(xpd = NA)
    legeps = .02
    if(minperprocname > 0){
      legend(x = -nrow(zpl) * legeps, y = -ncol(zpl) * legeps, xjust = 1, text.col = colv, legend = c(names.leg, "Other"), cex = .6)
    } else {
      legend(x = -nrow(zpl) * legeps, y = -ncol(zpl) * legeps, xjust = 1, text.col = colv[1:(length(colv) - 1)], legend = names.leg, cex = .6)
    }
    par(xpd = F)
    dev.off()
    file.copy(from = paste(go.plot.dir, "/", fnamc, sep = ""),
              to = paste(revision.paper.figures, "/", fnamc, sep = ""), overwrite = TRUE)
  }
}

###################################################################
# Table of GO term enrichment comparing UV and MV models
uvrows <- sapply(outl, function(x) NROW(x$uv.bigeff))
mvrows <- sapply(outl, function(x) NROW(x$mv.bigeff))
tabuvmvcomp <- table(uvrows, mvrows)
binl <- list(0, 1, 2, 3:5, 6:10, 11:20, 20:50)
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
tabout <- print(xtable(tabcomp, label = "tab:uvmvgocounts", caption = "Number of enriched BP GO terms per phenotype for the UV (left) and MV (top) models"),
                caption.placement = "top")
cat(tabout, file = paste(revision.text.numbers, "/mvugotab.txt", sep = ""))





longphnamv <- colnames(mvmat)[col(mvmat)]
longprocnamv <- phmap[match(longphnamv, phmap$nam), "procnam"]
fulltab <- data.frame(gonam = rownames(mvmat)[row(mvmat)], 'Procedure' = longprocnamv, 
                      'Phenotype' = longphnamv, Up = c(mvmat.nup), Down = c(mvmat.ndo), z = c(mvmat))
                      # 'Perturbation Direction' = paste0(c(mvmat.nup), '+', c(mvmat.ndo), '-'), z = c(mvmat))

subtab.go.terms <- c("lipid biosynthetic process", "growth", "B cell homeostasis", 
                     "immune system development", "regulation of lipid metabolic process",
                     "positive regulation of circadian sleep/wake cycle", "negative regulation of metabolic process", "equilibrioception")

fulltab <- fulltab[!is.na(fulltab$z), ]
rownames(fulltab) <- NULL
head(fulltab)
subtab <- fulltab[fulltab$gonam %in% subtab.go.terms, ]
mtcarsList <- split(subtab[, !names(subtab) %in% c("gonam", "z")], f = subtab$gonam)
attr(mtcarsList, "subheadings") <- paste0("\\emph{",
                                          names(mtcarsList), "}")
# attr(mtcarsList, "message") <- c("Line 1 of Message",
#                                  "Line 2 of Message")
xList <- xtableList(mtcarsList)
print.xtableList(xList)
tabout.test <- print.xtableList(xList, colnames.format = "multiple")
cat(tabout.test, file = paste(revision.text.numbers, "/mvugotab_test.txt", sep = ""))

goList <- split(fulltab[, !names(fulltab) %in% c("gonam", "z")], f = fulltab$gonam)
attr(goList, "subheadings") <- paste0(names(goList))
tabout.test <- ""
go.tab.dir <- paste0(revision.text.numbers, "/go_tables")
dir.create(go.tab.dir, showWarnings = F)
for(namc in names(goList)){
  try({
    xList <- xtable(goList[[namc]], caption = paste0("\\emph{", namc, "} genes perturb these phenotypes"), 
                    label = paste0("tab:", namc))
    tabout.test <- print(xList, include.rownames = F, caption.placement = "top")
    cat(tabout.test, file = paste0(go.tab.dir, "/gotab_", namc, ".txt"))
    
})  # print.xtableList(xList)
  # tabout.test <- paste(tabout.test, print.xtableList(xList, colnames.format = "multiple", include.rownames = F))
}


dput(rownames(mvmat.pl))


cat(tabout.test, file = paste(revision.text.numbers, "/mvugotab_test.txt", sep = ""))

str(xList)

data(mtcars)
mtcars <- mtcars[, 1:6]
mtcarsList <- split(mtcars, f = mtcars$cyl)
attr(mtcarsList, "subheadings") <- paste0("Number of cylinders = ",
                                          names(mtcarsList))
attr(mtcarsList, "message") <- c("Line 1 of Message",
                                 "Line 2 of Message")
xList <- xtableList(mtcarsList)
print.xtableList(xList)
tabout.test <- print.xtableList(xList, colnames.format = "multiple")
cat(tabout.test, file = paste(revision.text.numbers, "/mvugotab_test.txt", sep = ""))




############################################################
# Just checking KEGG for enrichment
# keggtab <- toTable(org.Mm.egPATH)
keggtab <- toTable(org.Mm.egGO)
keggtab2 <- keggtab[keggtab$gene_id %in% genemap2$entrez, ]
n.min.at.term <- 1
str(keggtab2)
idc <- c("go_id", "path_id")[1]
# keggtabc <- table(keggtab2$path_id)
keggtabc <- table(keggtab2[, idc])
keggkeep <- names(keggtabc)[keggtabc >= n.min.at.term]
keggtab3 <- keggtab2[keggtab2[, idc] %in% keggkeep, ]
keggentun <- unique(keggtab3$gene_id)
keggpathun <- unique(keggtab3[, idc])
entrez2kegg <- lapply(keggentun, function(x) keggtab3[keggtab3$gene_id == x, "path_id"])
names(entrez2kegg) <- keggentun
ensun <- unique(names(entrez2kegg))
ng <- length(ensun)
nk <- length(keggpathun)
keggensmat <- matrix(0, ng, nk, dimnames = list(ensun, keggpathun))
for(ensc in names(entrez2kegg)){
  keggensmat[ensc, which(colnames(keggensmat) %in% entrez2kegg[[ensc]])] <- 1
}
perml <- list()
methc <- "mv.bigeff"
pvalmatl <- list()
nperm <- 1
for(i in 1:(nperm + 1)){
  if(i == 1)
    keggensmat.use <- keggensmat
  if(i > 1)
    keggensmat.use[] <- c(keggensmat[sample(1:nrow(keggensmat)), ])
    # keggensmat.use[] <- c(keggensmat[, sample(1:ncol(keggensmat))])
  pvalmat <- matrix(NA, nph, nk, dimnames = list(ph.use, keggpathun))
  for(phc in ph.use){
    print(phc)
    ensmeasv <- genfacl[[methc]][[dirc]][[phc]]$measured
    enshitv.use <- enshitv <- genfacl[[methc]][[dirc]][[phc]]$hits
    # if(i == 1)
    #   enshitv.use <- enshitv
    # if(i > 1)
    #   enshitv.use <- sample(ensmeasv, length(enshitv))
    n.hit.in.kegg <- colSums(keggensmat.use[rownames(keggensmat.use) %in% enshitv.use, , drop = F])
    n.hit.not.in.kegg <- colSums(1 - keggensmat.use[rownames(keggensmat.use) %in% enshitv.use, , drop = F])
    n.measnothit.in.kegg <- colSums(keggensmat.use[rownames(keggensmat.use) %in% ensmeasv, , drop = F]) - n.hit.in.kegg
    n.measnothit.not.in.kegg <- colSums(1 - keggensmat.use[rownames(keggensmat.use) %in% ensmeasv, , drop = F]) - n.hit.not.in.kegg
    for(keggc in keggpathun){
      fisher.tab <- matrix(c(n.hit.in.kegg[keggc], n.hit.not.in.kegg[keggc], n.meas.in.kegg[keggc], n.meas.not.in.kegg[keggc]), 2, 2)
      fishout <- fisher.test(fisher.tab)
      pvalmat[phc, keggc] <- fishout$p.value
    }
    pvalmatl[[i]] <- pvalmat
  }
}  






#############################################
# Compare to mousenet network
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
mousenet <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2.txt"), sep = "\t")
mousenet.gs <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_gs.txt"), sep = "\t")
mousenet.lc <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_lc.txt"), sep = "\t")
mousenet.cx <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_cx.txt"), sep = "\t")
mousenet.gn <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_gn.txt"), sep = "\t")
mousenet.pg <- read.table(paste0(platdir, "/projects/impc_mv_analysis/data_in/mousenet_v2_pg.txt"), sep = "\t")
colnames(mousenet) <- colnames(mousenet.gs) <- colnames(mousenet.lc) <- 
  colnames(mousenet.cx) <- colnames(mousenet.gn) <- colnames(mousenet.pg) <- c("eg1", "eg2", "lr")

if("cl" %in% ls())
  stopCluster(cl)
cl <- makeCluster(12)
registerDoParallel(cl = cl)
dirc <- "bo"
methc <- c("ph", "proc", "fac")[1]
featvc <- switch(methc, ph = ph.use, proc = procord, fac = facnam)

outl <- foreach(phc = featvc, .verbose = T, .errorhandling = "pass") %dopar% {
  # phc = ph.use[1]
  pvall <- list()
  for(restype in c("uv.bigeff", "mv.bigeff")){
    myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
    allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
    all(myInterestingGenes %in% allPossibleGenes)
    # myInterestingGenes.symbol <- c(na.omit(sym2eg[match(myInterestingGenes, sym2eg$gene_id), "symbol"]))
    # allPossibleGenes.symbol <- c(na.omit(sym2eg[match(allPossibleGenes, sym2eg$gene_id), "symbol"]))
    # mousenet.sub <- mousenet[(mousenet$eg1 %in% allPossibleGenes | mousenet$eg2 %in% allPossibleGenes), ]
    mousenet.use <- list(mousenet.gs, mousenet.lc, mousenet.cx, mousenet.gn, mousenet.pg, mousenet)[[1]]
    # if(any(!is.na(mousenet.use$lr))){
    #   mousenet.sub <- mousenet.use[which((mousenet.use$eg1 %in% allPossibleGenes & mousenet.use$eg2 %in% allPossibleGenes) & mousenet.use$lr > 0), ]
    # } else {
    mousenet.sub <- mousenet.use[(mousenet.use$eg1 %in% allPossibleGenes & mousenet.use$eg2 %in% allPossibleGenes), ]
    # }
    str(mousenet.sub)
    npair <- sum(mousenet.sub$eg1 %in% myInterestingGenes & mousenet.sub$eg2 %in% myInterestingGenes)
    perm.npair <- c()
    for(i in 1:1000){
      myInterestingGenes.null <- sample(allPossibleGenes, length(myInterestingGenes))
      perm.npair[i] <- sum(mousenet.sub$eg1 %in% myInterestingGenes.null & mousenet.sub$eg2 %in% myInterestingGenes.null)
    }
    pval.est <- mean(perm.npair >= npair)
    pvall[[restype]] <- pval.est
  }
  return(list(pval = pvall, ph = phc))
}
stopCluster(cl)
if(methc == "ph")
  names(outl) <- phmap[match(ph.use, phmap$ph), "nam"]
if(methc %in% c("proc", "fac"))
  names(outl) <- featvc#phmap[match(ph.use, phmap$ph), "nam"]
pvalv.uv <- sapply(outl, function(v) v$pval$uv.bigeff)
pvalv.mv <- sapply(outl, function(v) v$pval$mv.bigeff)
qvalv.uv <- p.adjust(pvalv.uv, meth = "BH")
qvalv.mv <- p.adjust(pvalv.mv, meth = "BH")
thc <- .05
table(pvalv.uv < thc, pvalv.mv < thc)
table(qvalv.uv < thc, qvalv.mv < thc)
names(qvalv.mv)[qvalv.mv < thc]

mean(p.adjust(pvalv.mv, meth = "BH") < thc)

mean(pvalv.mv < thc)
mean(pvalv.uv < thc)














###################



restype <- "mv.bigeff"
methc <- c("ph", "proc", "fac")[2]
featvc <- switch(methc, ph = ph.use, proc = procord, fac = facnam)
pvalv <- c()
for(phc in featvc){
  myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
  allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
  geneList.bin <- as.integer(allPossibleGenes %in% myInterestingGenes)
  maptab <- data.frame(eg = allPossibleGenes, sig = geneList.bin)
  maptab[, c("loc", "chr")] <- sym2eg[match(allPossibleGenes, sym2eg$gene_id), c("loc", "chr")]
  maptab <- maptab[order(maptab$chr, maptab$loc), ]
  maptab$sym <- sym2eg[match(maptab$eg, sym2eg$gene_id), "symbol"]
  chrun <- unique(maptab$chr)
  ntot <- nrow(maptab)
  nsig <- sum(maptab$sig)
  # str(maptab)
  table(maptab$chr)
  vecsig <- table(maptab$chr[maptab$sig == 1])
  vecsig <- vecsig[match(chrun, names(vecsig))]
  vecsig[is.na(vecsig)] <- 0
  vecnon <- table(maptab$chr[maptab$sig == 0])
  vecnon <- vecnon[match(chrun, names(vecnon))]
  vecnon[is.na(vecnon)] <- 0
  
  tab.chi <- cbind(vecsig, vecnon)
  tab.chi <- tab.chi[rowSums(tab.chi) > 0, ]
  chi.out <- chisq.test(tab.chi)
  pvalv[phc] <- chi.out$p.value
}

sort(pvalv)
sort(p.adjust(pvalv, meth = "BH"))







#################


############################################
# Positional analysis
methc <- c("ph", "proc", "fac")[2]
featvc <- switch(methc, ph = ph.use, proc = procord, fac = facnam)
dirc <- "bo"
if("cl" %in% ls())
  stopCluster(cl)
cl <- makeCluster(25)
registerDoParallel(cl = cl)
par(mfrow = c(5, 5))
phsub <- sample(ph.use, 25)
nperm <- 500
outl <- foreach(phc = featvc, .verbose = T, .errorhandling = "pass") %dopar% {#for(phc in phsub){#phc <- ph.use[20]#
  library(topGO)
  library(GOfuncR)
  library(org.Mm.eg.db)
  allresl <- list()
  for(restype in c("uv.bigeff", "mv.bigeff")){
    myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
    allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
    geneList.bin <- as.integer(allPossibleGenes %in% myInterestingGenes)
    maptab <- data.frame(eg = allPossibleGenes, sig = geneList.bin)
    maptab[, c("loc", "chr")] <- sym2eg[match(allPossibleGenes, sym2eg$gene_id), c("loc", "chr")]
    maptab <- maptab[order(maptab$chr, maptab$loc), ]
    maptab$sym <- sym2eg[match(maptab$eg, sym2eg$gene_id), "symbol"]
    chrun <- unique(maptab$chr)
    ntot <- nrow(maptab)
    nsig <- sum(maptab$sig)
    str(maptab)
    table(maptab$chr)
    chisq.test(cbind(table(maptab$chr[maptab$sig == 1]), table(maptab$chr[maptab$sig == 0])))
    filthalfwidv <- c(3, 5, 7, 10, 15, 20, 30)[2]##[1:3]
    filtshapel <- list()
    for(filtwidc in filthalfwidv){
      filtshapel[[as.character(filtwidc)]] <- dnorm(seq(-3, 3, len = 2 * filtwidc + 1))
      filtshapel[[as.character(filtwidc)]] <- filtshapel[[as.character(filtwidc)]] / sum(filtshapel[[as.character(filtwidc)]])
    }
    
    
    
    overlap.inds <- c()
    act <- rep(0, ntot)
    for(filtwidc in filthalfwidv){
      filtc <- filter(maptab$sig, filter = filtshapel[[as.character(filtwidc)]])
      # overlap.inds <- outer(match(chrun, maptab$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
      # overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
      # filtc[overlap.inds] <- NA
      act <- pmax(act, filtc, na.rm = T)
    }
    permmat <- matrix(NA, ntot, nperm)
    for(it in 1:nperm){
      maptab.perm <- maptab
      maptab.perm$sig <- sample(maptab.perm$sig)
      permc <- rep(0, ntot)
      for(filtwidc in filthalfwidv){
        filtc <- filter(maptab.perm$sig, filter = filtshapel[[as.character(filtwidc)]])
        # overlap.inds <- outer(match(chrun, maptab.perm$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
        # overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
        # filtc[overlap.inds] <- NA
        permc <- pmax(permc, filtc, na.rm = T)
      }
      permmat[, it] <- permc#filter(maptab.perm$sig, filter = rep(1, filtwid))
    }
    allresl[[restype]] <-  list(permmat = permmat, act = act)
  }
  return(allresl)
}
names(outl) <- featvc
stopCluster(cl)
rm(cl)

str(outl)

# plot(outl[[1]]$mv.bigeff$act)
# plot(outl[[1]]$mv.bigeff$permmat[, 1])

relmatl <- list()
fwerc <- .05
for(restype in c("uv.bigeff", "mv.bigeff")){
  relmatl[[restype]] <- sapply(outl, function(x){
                  thr.up <- quantile(apply(x[[restype]]$permmat, 2, function(v) max(v, na.rm = T)), 1 - fwerc, na.rm = T)
                return(x[[restype]]$act / thr.up)
    })
}

# colnames(relmatl$mv.bigeff) <- ph.use
# plot(relmatl[[1]][[1]])


# apply(relmatl$uv.bigeff, 2, function(v) max(v, na.rm = T))

str(relmatl)
max.rel <- apply(relmatl$mv.bigeff, 2, function(v) max(v, na.rm = T))
nplph <- 10
plph <- colnames(relmatl$mv.bigeff)[order(-max.rel)[1:nplph]]
graphics.off()
matplot(relmatl$mv.bigeff[, plph], ty = "l")
abline(h = 1)

phord <- ph.use[order(phmap[match(ph.use, phmap$ph), "procnam"])]
graphics.off()
par(mfrow = c(2, 1))
image(relmatl$mv.bigeff[rowSums(is.na(relmatl$mv.bigeff)) == 0, phord])
image(relmatl$mv.bigeff[rowSums(is.na(relmatl$mv.bigeff)) == 0, phord] > 1)


############################################################################
# Combined positional analysis (aggregating the number of hits across phenotypes within gene)
# ssmat <- (abs(ebl$t) > ebl$th) * sign(ebl$t)
ssmat <- resl.out$eb$signsig.bigeff[true.use, ]
maptab <- data.frame(genzyg = rownames(ssmat))
maptab$impc.id <- sapply(strsplit(maptab$genzyg, spl = "_"), function(v) v[1])
maptab$zyg <- sapply(strsplit(maptab$genzyg, spl = "_"), function(v) v[2])
maptab$eg <- genemap2[match(maptab$impc.id, genemap2$genotype_id), "entrez"]
maptab[, c("loc", "chr")] <- sym2eg[match(maptab$eg, sym2eg$gene_id), c("loc", "chr")]
maptab <- maptab[order(match(maptab$zyg, c(1, 0, 2))), ]
genshun <- unique(maptab$eg)
maptab <- maptab[match(genshun, maptab$eg), ]
maptab <- maptab[order(maptab$chr, maptab$loc), ]
chrun <- unique(maptab$chr)
phord <- ph.use[order(phmap[match(ph.use, phmap$ph), "procnam"])]
ssmat <- ssmat[maptab$genzyg, phord]
ntot <- nrow(ssmat)
agghits <- rowSums(abs(ssmat[, , drop = F]))
filthalfwidv <- c(5)#10, 15, 20, 30, 50)[1]#[1:3]
filtshapel <- list()
for(filtwidc in filthalfwidv){
  filtshapel[[as.character(filtwidc)]] <- dnorm(seq(-3, 3, len = 2 * filtwidc + 1))
  filtshapel[[as.character(filtwidc)]] <- filtshapel[[as.character(filtwidc)]] / sum(filtshapel[[as.character(filtwidc)]])
}
overlap.inds <- c()
act <- rep(0, ntot)
for(filtwidc in filthalfwidv){
  filtc <- filter(agghits, filter = filtshapel[[as.character(filtwidc)]])
  overlap.inds <- outer(match(chrun, maptab$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
  overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
  filtc[overlap.inds] <- NA
  act <- pmax(act, filtc, na.rm = T)
}
plot(act)
permmat <- matrix(NA, ntot, nperm)
nperm <- 500

for(it in 1:nperm){
  agghits.perm <- sample(agghits)
  permc <- rep(0, ntot)
  for(filtwidc in filthalfwidv){
    filtc <- filter(agghits.perm, filter =  filtshapel[[as.character(filtwidc)]])
    overlap.inds <- outer(match(chrun, maptab$chr), seq(-2 * filtwidc, 2 * filtwidc, by = 1), '+')
    overlap.inds <- overlap.inds[overlap.inds >= 1 & overlap.inds < ntot]
    filtc[overlap.inds] <- NA
    permc <- pmax(permc, filtc, na.rm = T)
  }
  permmat[, it] <- permc#filter(maptab.perm$sig, filter = rep(1, filtwid))
}
thr.up <- quantile(apply(permmat, 2, function(v) max(v, na.rm = T)), 1 - fwerc, na.rm = T)
plot(act / thr.up, ty = "l")
# plot(act)






str(ssmat)

#



myInterestingGenes <- genfacl[[restype]][[dirc]][[phc]]$hits
allPossibleGenes <- genfacl[[restype]][[dirc]][[phc]]$measured
geneList.bin <- as.integer(allPossibleGenes %in% myInterestingGenes)
maptab <- data.frame(eg = allPossibleGenes, sig = geneList.bin)
maptab[, c("loc", "chr")] <- sym2eg[match(allPossibleGenes, sym2eg$gene_id), c("loc", "chr")]
maptab <- maptab[order(maptab$chr, maptab$loc), ]
maptab$sym <- sym2eg[match(maptab$eg, sym2eg$gene_id), "symbol"]
ssmat <- resl.out$eb$signsig.bigeff[true.use, ]
# ssmat <- resl.out$uv$signsig.bigeff[true.use, ]
agghits <- rowSums(abs(ssmat[, , drop = F]))
maptab$chrloc <- 1e10 * maptab$chr + maptab$loc
maptab$sym <- sym2eg[match(maptab$eg, sym2eg$gene_id), "symbol"]
maptab$impc.gen <- genemap[match(maptab$sym, genemap$gene_symbol), "genotype_id"]
maptab$agg <- agghits[match(paste0(maptab$impc.gen, "_1"), names(agghits))]
maptab <- maptab[!is.na(maptab$agg), ]
maptab <- maptab[maptab$agg >= 1, ]
maptab <- maptab[order(maptab$chrloc), ]
maptab$dist <- c(NA, diff(maptab$chrloc))
head(maptab)
sort(maptab$dist)[1:10]
str(maptab)
maptab[maptab$dist %in% sort(maptab$dist)[1:60], ]




agghits <- rowSums(abs(ssmat[, , drop = F]))
str(genemap)
ssmat.eg <- sapply()








# phc <- "IMPC_OFD_012_001"
# phc <- "IMPC_IPG_012_001"
# phc <- "IMPC_CBC_017_001"
# phc <- "IMPC_CBC_025_001"
# phc <- "IMPC_CBC_004_001"
# for(restype in c("uv", "mv"))
#   print(genfacl[[restype]][[dirc]][[phc]])
# phc <- "IMPC_DXA_002_001"#ph.use[1]#for(facc in facord){#facc <- 3

# BiocInstaller::biocLite(c("GOfuncR"))
# BiocInstaller::biocLite(c("org.Mm.eg.db"))
# ## top GO-categories per domain
# by(stats, stats$ontology, head, n=3)
# 
# 
# GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = entrez2go)
# test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
# resultFisher <- getSigGroups(GOdata, test.stat)
# allRes <- GenTable(GOdata, Fis = resultFisher, topNodes = 50)
# print(phmap[match(phc, phmap$ph), ])
# print(allRes)
# myInterestingGenes.symbol <- c(na.omit(sym2eg[match(myInterestingGenes, sym2eg$gene_id), "symbol"]))
# allPossibleGenes.symbol <- c(na.omit(sym2eg[match(allPossibleGenes, sym2eg$gene_id), "symbol"]))
# geneList <- factor(as.integer(allPossibleGenes %in% myInterestingGenes))
# maptab <- maptab[order(match(maptab$zyg, c(1, 0, 2))), ]
# genshun <- unique(maptab$eg)
# maptab <- maptab[match(genshun, maptab$eg), ]
# 
# maptab <- data.frame(genzyg = rownames(ssmat))
# maptab$impc.id <- sapply(strsplit(maptab$genzyg, spl = "_"), function(v) v[1])
# maptab$zyg <- sapply(strsplit(maptab$genzyg, spl = "_"), function(v) v[2])
# maptab$eg <- genemap2[match(maptab$impc.id, genemap2$genotype_id), "entrez"]
# maptab[, c("loc", "chr")] <- sym2eg[match(maptab$eg, sym2eg$gene_id), c("loc", "chr")]
# maptab <- maptab[order(match(maptab$zyg, c(1, 0, 2))), ]
# genshun <- unique(maptab$eg)
# maptab <- maptab[match(genshun, maptab$eg), ]
# maptab <- maptab[order(maptab$chr, maptab$loc), ]
# nrow(maptab)
# egun <- unique(maptab$eg)
# length(egun)
# maptab <- maptab[match(egun, maptab$eg), ]
############################################
# Note that there are currently some Symbols not being mapped to Entrez IDs
# str(sym2eg)
# "52808" %in% sym2eg$gene_id
# table(table(sym2eg$symbol))
# table(table(sym2eg$gene_id))
# all(genemap.in$gene_symbol %in% sym2eg$symbol)
# genemap.in$gene_symbol[! genemap.in$gene_symbol %in% sym2eg$symbol]
#####################

# resimp.truelines <- resimp[resimp$geno %in% true.use, ]
# uvl <- lapply(uvl, function(M) M[true.use, ph.use])
# ebl <- lapply(ebl, function(M) M[true.use, ph.use])
# uvl$th <- ebl$th <- matrix(NA, ng, nph, dimnames = list(true.use, ph.use))
# uvl$th[cbind(resimp.truelines$geno, resimp.truelines$ph)] <- resimp.truelines$uv.th.final
# ebl$th[cbind(resimp.truelines$geno, resimp.truelines$ph)] <- resimp.truelines$eb.th.final
# sigl <- list()
# sigl$uv <- sigl$mv <- sigl$f <- list()
# sigl$uv$up <- (uvl$t > uvl$th)[true.use, ]
# sigl$uv$do <- (uvl$t < -uvl$th)[true.use, ]
# sigl$uv$bo <- (abs(uvl$t) > uvl$th)[true.use, ]
# sigl$mv$up <- (ebl$t > ebl$th)[true.use, ]
# sigl$mv$do <- (ebl$t < -ebl$th)[true.use, ]
# sigl$mv$bo <- (abs(ebl$t) > ebl$th)[true.use, ]
# sigl$f$up <- (fl$t > fl$th)[true.use, ]
# sigl$f$do <- (fl$t < -fl$th)[true.use, ]
# sigl$f$bo <- (abs(fl$t) > fl$th)[true.use, ]

# 
# fl <- list()
# fl$t <- resl.out$eb.fac$mn / resl.out$eb.fac$sd
# str(fl)
# str(resimp)
# 
# fl$th <- fl$t
# if(length(unique(resimp$eb.fac.th.final)) > 1)
#   stop("Multiple factor sig thresh, but code written for single thresh")
# fl$th[] <- resimp$eb.fac.th.final[1]
# facdat <- resl$f$facdat

# nfac <- 24
# vc.type <- c("vari", "pro")[2]
# load(file = paste0(meth.comp.output.dir, "/ebmix_results_nf_", nfac, "_vt_", vc.type, ".RData"))
# facdat <- resl$f$facdat
# facs <- resl$f$loadings.ord
# facs[, 21][order(abs(facs[, 21]))]
# resl$f$ref.lines
# str(resl$f, m = 2)

# myInterestingGenes.uv <- myInterestingGenes
# allPossibleGenes.uv <- allPossibleGenes
# myInterestingGenes.symbol.uv <- myInterestingGenes.symbol
# allPossibleGenes.symbol.uv <- allPossibleGenes.symbol
# 
# myInterestingGenes.mv <- myInterestingGenes
# allPossibleGenes.mv <- allPossibleGenes
# myInterestingGenes.symbol.mv <- myInterestingGenes.symbol
# allPossibleGenes.symbol.mv <- allPossibleGenes.symbol
# 
# 
# length(intersect(myInterestingGenes.uv, myInterestingGenes.mv))
# length(setdiff(myInterestingGenes.uv, myInterestingGenes.mv))
# length(setdiff(myInterestingGenes.mv, myInterestingGenes.uv))
# 
# length(intersect(allPossibleGenes.uv, allPossibleGenes.mv))
# length(setdiff(allPossibleGenes.uv, allPossibleGenes.mv))
# length(setdiff(allPossibleGenes.mv, allPossibleGenes.uv))
# 
# res_hyper.uv <- res_hyper
# res_hyper.mv <- res_hyper
# 
# res.mv <- res_hyper.mv$results
# res.uv <- res_hyper.uv$results
# res.mv[res.mv$node_name == "sensory perception of sound", ]
# res.uv[res.uv$node_name == "sensory perception of sound", ]
# sound.ents <- get_anno_genes(go_ids = "GO:0007605", database = 'Mus.musculus', genes = NULL, annotations = NULL,
#                              term_df = NULL, graph_path_df = NULL, godir = NULL)
# 
# sound.ents$eg <- sym2eg[match(sound.ents$gene, sym2eg$symbol), "gene_id"]
# sound.ents$mv.sig <- sound.ents$eg %in% myInterestingGenes.mv
# sound.ents$uv.sig <- sound.ents$eg %in% myInterestingGenes.uv
# sound.ents$in.mv <- sound.ents$eg %in% allPossibleGenes.mv
# sound.ents$in.uv <- sound.ents$eg %in% allPossibleGenes.uv
# resimp$sym <- genemap[match(resimp$geno.sh, genemap$genotype_id), "gene_symbol"]
# resimp$eg <- sym2eg[match(resimp$sym, sym2eg$symbol), "gene_id"]
# resimp$sound.gene <- resimp$eg %in% sound.ents$eg
# resimp$proc <- phmap[match(resimp$ph, phmap$ph), "procnam"]
# 
# resimpsub <- resimp[which(resimp$proc == "Auditory Brain Stem Response" & resimp$eb.perm.signsig != 0), ]
# plot(resimpsub$eb.t, col = ifelse(resimpsub$sound.gene, 2, 1))
# 
# boxplot(abs(eb.t) ~ sound.gene, data = resimpsub)
# 
# summary(abs(resimpsub$eb.t[resimpsub$sound.gene]))
# summary(abs(resimpsub$eb.t[!resimpsub$sound.gene]))
# summary(abs(resimpsub$uv.t[resimpsub$sound.gene]))
# summary(abs(resimpsub$uv.t[!resimpsub$sound.gene]))
# 
# sound.ents$cen <- resimp[match(sound.ents$eg, resimp$eg), "cenlong"]
# sound.ents.sub <- sound.ents[sound.ents$in.mv, ]
# table(sound.ents$mv.sig, sound.ents$uv.sig)
# table(sound.ents.sub$mv.sig, sound.ents.sub$uv.sig)
# 
# mvmat <- matrix(c(16, 22, length(myInterestingGenes.mv), length(allPossibleGenes.mv) - length(myInterestingGenes.mv)), 2, 2)
# uvmat <- matrix(c(11, 22, length(myInterestingGenes.uv), length(allPossibleGenes.uv) - length(myInterestingGenes.uv)), 2, 2)
# mvmat
# fisher.test(mvmat)
# uvmat
# fisher.test(uvmat)

