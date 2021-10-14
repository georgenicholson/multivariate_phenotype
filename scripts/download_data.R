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
  download.file(url = "https://www.dropbox.com/s/wz6mg35n502au6a/multivariate_phenotype_paper_supporting_data.zip?dl=1", temp, mode = "wb")
  unzip(temp, exdir = control$data_dir)
  unlink(temp)
}


Data_all <- list()
for (Data in c("eqtl", "impc")) {
  Data_all[[Data]] <- list()
  ###################################
  # Load data
  if (Data == "eqtl") {
    mash.in <- readRDS("data/eqtl/MatrixEQTLSumStats.Portable.ld2.Z.rds")
    Y_raw <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
    S_raw <- Y_raw / rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
  }
  if (Data == "impc") {
    Y_raw <- readRDS("data/impc/Yhatmat.RDS")
    S_raw <- readRDS("data/impc/Shatmat.RDS")
    Data_all[[Data]]$linemap <- readRDS("data/impc/linemap.RDS")
    Data_all[[Data]]$reflinemap <- readRDS("data/impc/reflinemap.RDS")
    Data_all[[Data]]$phmap <- readRDS("data/impc/phmap.RDS")
    Data_all[[Data]]$cenmap <- readRDS("data/impc/cenmap.RDS")
  }
  
  
  Y_zeroed <- Y_raw
  S_zeroed <- S_raw
  
  
  Y_zeroed[is.na(Y_raw)] = 0
  S_zeroed[is.na(Y_raw)] = control$prior.sd.on.unobserved.thetahat
  
  Data_all[[Data]]$Y_raw <- Y_raw
  Data_all[[Data]]$S_raw <- S_raw
  Data_all[[Data]]$Y_zeroed <- Y_zeroed
  Data_all[[Data]]$S_zeroed <- S_zeroed
  Data_all[[Data]]$sam_names <- rownames(Y_raw)
  Data_all[[Data]]$meas_names <- colnames(Y_raw)
  Data_all[[Data]]$N_all <- length(Data_all[[Data]]$sam_names)
  Data_all[[Data]]$P_all <- length(Data_all[[Data]]$meas_names)
}

saveRDS(Data_all, file = control$Data_all_file)

#   Y.eml <- S.eml <- list()
#   Y.eml$raw <- Yhat
#   Y.eml$zero <- Y
#   S.eml$raw <- smat
#   S.eml$zero <- S
#   dimnaml[[Data]]$sam.names <- rownames(Yhat)
#   dimnaml[[Data]]$meas.names <- colnames(Yhat)
#   dimnaml[[Data]]$N.all <- length(dimnaml[[Data]]$sam.names)
#   dimnaml[[Data]]$P.all <- length(dimnaml[[Data]]$meas.names)
#   
# 
# ##########################################
# # Get dimnames from data sets
# dimnaml <- list()
# for(Data in c("impc", "eqtl")){
#   dimnaml[[Data]] <- list()
#   if (Data == "eqtl") {
#     mash.in <- readRDS("data/eqtl/MatrixEQTLSumStats.Portable.ld2.Z.rds")
#     Yhat <- rbind(mash.in$strong.b, mash.in$random.b, mash.in$random.test.b)
#     smat <- Yhat / rbind(mash.in$strong.z, mash.in$random.z, mash.in$random.test.z)
#   }
#   if (Data == "impc") {
#     Yhat <- readRDS("data/impc/Yhatmat.RDS")
#     smat <- readRDS("data/impc/Shatmat.RDS")
#     linemap <- readRDS("data/impc/linemap.RDS")
#     reflinemap <- readRDS("data/impc/reflinemap.RDS")
#     phmap <- readRDS("data/impc/phmap.RDS")
#   }
#   dimnaml[[Data]]$sam.names <- rownames(Yhat)
#   dimnaml[[Data]]$meas.names <- colnames(Yhat)
#   dimnaml[[Data]]$N.all <- length(dimnaml[[Data]]$sam.names)
#   dimnaml[[Data]]$P.all <- length(dimnaml[[Data]]$meas.names)
# }
# 
