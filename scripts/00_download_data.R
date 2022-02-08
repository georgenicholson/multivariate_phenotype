rm(list = ls())
##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control contains parameter settings
control <- get_control_parameters_mv()

##########################################
# Download data
if (!"data/impc" %in% list.dirs(control$data_dir)) {
  temp <- tempfile()
  download.file(url = "https://www.dropbox.com/s/hpsyomy7y9paxll/multivariate_phenotype_paper_supporting_data.zip?dl=1", temp, mode = "wb")
  unzip(temp, exdir = control$data_dir)
  unlink(temp)
}

###################################
# Organise data
Data_all <- list()
for (Data in c("eqtl", "impc")) {
  Data_all[[Data]] <- list()
  if (Data == "eqtl") {
    mash.in <- readRDS("data/eqtl/MatrixEQTLSumStats.Portable.ld2.Z.rds")
    Data_all[[Data]]$Y_raw <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
    Data_all[[Data]]$S_raw <- Data_all[[Data]]$Y_raw / rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
    data.type <- rep(c("strong", "random.train", "random.test"), 
                     times = c(nrow(mash.in$strong.z), nrow(mash.in$random.z), nrow(mash.in$random.test.z)))
    Data_all[[Data]]$snpmap <- data.frame(snp = rownames(Data_all[[Data]]$Y_raw), data.type = data.type)
    
  }
  if (Data == "impc") {
    Data_all[[Data]]$Y_raw <- readRDS("data/impc/Yhatmat.RDS")
    Data_all[[Data]]$S_raw <- readRDS("data/impc/Shatmat.RDS")
    Data_all[[Data]]$linemap <- readRDS("data/impc/linemap.RDS")
    Data_all[[Data]]$reflinemap <- readRDS("data/impc/reflinemap.RDS")
    Data_all[[Data]]$phmap <- readRDS("data/impc/phmap.RDS")
    Data_all[[Data]]$cenmap <- readRDS("data/impc/cenmap.RDS")
    Data_all[[Data]]$genemap <- readRDS("data/impc/genemap.RDS")
  }
  
  Data_all[[Data]]$Y_zeroed <- Data_all[[Data]]$Y_raw
  Data_all[[Data]]$S_zeroed <- Data_all[[Data]]$S_raw
  Data_all[[Data]]$Y_zeroed[is.na(Data_all[[Data]]$Y_raw)] = 0
  Data_all[[Data]]$S_zeroed[is.na(Data_all[[Data]]$Y_raw)] = control$prior.sd.on.unobserved.thetahat
  # Data_all[[Data]]$Y_raw <- Y_raw
  # Data_all[[Data]]$S_raw <- S_raw
  # Data_all[[Data]]$Y_zeroed <- Y_zeroed
  # Data_all[[Data]]$S_zeroed <- S_zeroed
  Data_all[[Data]]$sam_names <- rownames(Data_all[[Data]]$Y_raw)
  Data_all[[Data]]$meas_names <- colnames(Data_all[[Data]]$Y_raw)
  Data_all[[Data]]$N_all <- length(Data_all[[Data]]$sam_names)
  Data_all[[Data]]$P_all <- length(Data_all[[Data]]$meas_names)
}




#############################################################
#create and save ordering of procedures for plotting, so that labels don't overlap
###################################
phmap <- Data_all$impc$phmap
procall <- unique(phmap$procnam)
phenall <- unique(phmap$ph)
proctab <- table(phmap[match(phenall, phmap$ph), "procnam"])
proctabord <- sort(proctab)
procord <- c()
smallbig <- 0
while(length(proctabord) > 0){
  if(smallbig == 0){
    procord <- c(procord, names(proctabord)[1])
    proctabord <- proctabord[-1]
  }
  if(smallbig == 1){
    procord <- c(procord, names(proctabord)[length(proctabord)])
    proctabord <- proctabord[-length(proctabord)]
  }
  smallbig <- 1 - smallbig
}
inds.swap <- which(procord %in% c("Intraperitoneal glucose tolerance test (IPGTT)", "Heart Weight"))
procord[inds.swap] <- procord[rev(inds.swap)]
procord[]
phord <- unique(phmap[order(match(phmap$procnam, procord)), "ph"])
Data_all$impc$phord <- phord
Data_all$impc$procord <- procord


saveRDS(Data_all, file = control$Data_all_file)



