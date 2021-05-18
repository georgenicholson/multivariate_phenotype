rm(list = ls())
# source("R_files/impc_mv_parameters.R")

fn.source.out <- sapply(list.files("R_files/functions", full.names = T), source)

control <- get_control_parameters_mv()
# undebug({
# debugonce(create_table_of_analyses)
runtab <- create_table_of_analyses(check_status = T)
# })


save(runtab, file = file.runtab)
scen.do <- which(!runtab$res.ok)
runtab[scen.do, ]
ntimes <- runtab$ntimes[scen]
partc <- c("debug", "debug11", "debug6", "debug8")[2]
file.copy(from = paste0(R.file.dir, "/impc_mv_paper_code/run_EM.R"),
          to = paste0(R.file.dir, "/impc_mv_analysis/run_EM_algorithm_temp_sims.R"), overwrite = T)
for(scen in scen.do){#scen <- 1#
  ntimes <- runtab$ntimes[scen]
  for(seedc in 1:ntimes){#seedc <- 1#
    # for(seedc in 1:4){#seedc <- 1#
    var.in.name.bash <- unique(c(var.in.name, var.in.name.ed, var.in.name.mash))
    var.in.name.bash <- var.in.name.bash[order(match(var.in.name.bash, colnames(runtab)))]
    # var.in.name.bash <- setdiff(colnames(runtab), c("ma", "fu", "ntimes", "mem"))[1:8]
    var.in.name.use <- var.in.name.bash[var.in.name.bash %in% names(runtab)]
    var.in.name.use <- var.in.name.use[!is.na(runtab[scen, match(var.in.name.use, names(runtab))])]
    var.in.name.use <- c(var.in.name.use, "seed")
    namc <- gsub("XXX", seedc, runtab[scen, "file.base"])
    #paste(c(paste(var.in.name.use, runtab[scen, var.in.name.use], sep = "_"), paste0("seed_", seedc)), collapse = "_")
    bashscript <- paste0(bash.file.dir, "/test_script_", namc, ".sh")
    # file.out <- paste0(prefix, namc, ".RData")#paste0("restab_", namc, ".RData")
    # file.out.restab <- paste0(prefix, "_restab_", namc, ".RData")
    # file.out.resl <- paste0(prefix, "_res_", namc, ".RData")
    # if((!file.exists(paste0(meth.comp.output.dir, "/", file.out.restab))) | replace.files){
    output.file <- file(bashscript, "wb")
    cat(file = output.file, append = F, sep = "\n",
        c("#!/bin/bash",
          "#SBATCH --ntasks=1",
          "#SBATCH --output=/mnt/s/job_output_files/job_%j.out",
          "#SBATCH --error=/mnt/s/job_output_files/job_%j.err",
          paste0("#SBATCH --partition ", partc),
          paste0("#SBATCH --job-name=", namc),
          paste0("#SBATCH --mem-per-cpu=", runtab$mem[scen], "MB"),
          paste("srun Rscript --no-restore --no-save",
              "/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/run_EM_algorithm_temp_sims.R",
              seedc, runtab$N[scen], runtab$P[scen], runtab$bo[scen], runtab$si[scen], 
              runtab$da[scen], runtab$me[scen], runtab$ma[scen], runtab$fu[scen], 
              runtab$ss[scen], runtab$nSig[scen], runtab$EDmeth[scen], runtab$EDtol[scen], 
              runtab$EMtol[scen], runtab$EMmash[scen], runtab$EMfm[scen], runtab$EMbic[scen], 
              runtab$EMwish[scen], runtab$EMK[scen], runtab$EMKup[scen], 
              runtab$rand[scen], runtab$loocv[scen],
              "--profile=task", sep = " ")))
    if(Sys.info()["sysname"] == "Windows"){
      system(paste("wsl sbatch", gsub("C:/", "/mnt/c/", bashscript)))
      # system(paste("wsl sbatch", gsub("C:/Users/nicho/Documents/bauer_sync", "/mnt/x/", bashscript)))
    } else {
      system(paste("sbatch", bashscript))
    }
    close(output.file)
  }
  # Sys.sleep(2)
}






# "#SBATCH --output=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%j.out",
# "#SBATCH --error=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%j.err",












readLines(bashscript)

scen <- 1
for(seedc in 1:2){#seedc <- 1#
  namc <- paste(c(paste(var.in.name, runtab[scen, var.in.name], sep = "_"), paste0("seed_", seedc)), collapse = "_")
  bashscript <- paste0(bash.file.dir, "/test_script_", namc, ".sh")
  file.out <- paste0(prefix, namc, ".RData")#paste0("restab_", namc, ".RData")
  file.out.restab <- paste0(prefix, "_restab_", namc, ".RData")
  file.out.resl <- paste0(prefix, "_res_", namc, ".RData")
load(file = paste0(meth.comp.output.dir, "/", file.out.resl))
load(file = paste0(meth.comp.output.dir, "/", file.out.restab))
restab
object.size(resl.store)

}

# dataset <- as.character(arguments[grep("--args", arguments) + 7])
# methods <- as.character(arguments[grep("--args", arguments) + 8])
# max.num.sams.for.testing <- as.character(arguments[grep("--args", arguments) + 9])
# full.analysis <- as.character(arguments[grep("--args", arguments) + 10])


comp.typev <- c("ref.lines.post", "het.hom.post")
pr.grpv <- c("", "comp", "imp", "impcomp")
tab.temp <- expand.grid(comp.typev, pr.grpv)
grpnam <- paste0(tab.temp[, 1], tab.temp[, 2])
nth <- 0
######################################
#collect results
restabl <- list()
restab.checkl <- list()
scen = 7#for(scen in 1:nrow(runtab)){#scen <- 1#
  ntimes <- runtab$ntimes[scen]
  restab.num.sum <- restab.num.sumsq <- restab.denom.sum <- 0
  nsuccess <- 0
  for(seedc in 1:ntimes){#seedc <- 1#
    namc <- paste(c(paste(var.in.name, runtab[scen, var.in.name], sep = "_"), paste0("seed_", seedc)), collapse = "_")
    bashscript <- paste0(bash.file.dir, "/test_script_", namc, ".sh")
    file.out <- paste0(prefix, namc, ".RData")#paste0("restab_", namc, ".RData")
    full.file.out <- paste0(meth.comp.output.dir, "/", file.out)
    if(!file.exists(full.file.out))
      next
# for(scen in 1:nrow(runtab)){#scen <- 1#
#   N <- runtab$N[scen]
#   P <- runtab$P[scen]
#   ntimes <- runtab$ntimes[scen]
  # for(seedc in 1:ntimes){#seedc <- 1#
    
    try({
      # namc <- paste(c(paste(var.in.name, runtab[scen, var.in.name], sep = "_"), paste0("seed_", seedc)), collapse = "_")
      # file.out <- paste0(meth.comp.output.dir, "/", namc, ".RData")
      # namc <- paste(c(paste(var.in.name, runtab[scen, var.in.name], sep = "_"), paste0("seed_", seedc)), collapse = "_")
      # fc <- paste0(meth.comp.output.dir, "/restab_N_", N, "_P_", P, "_seed_", seedc, ".RData")
      # if(difftime(Sys.time(), file.info(full.file.out)$mtime, units = "hours") < Inf){
        load(file = full.file.out)
        restab.checkl[[seedc]] <- restab
        for(grpnamc in grpnam){
          nnam <- paste0(grpnamc, ".n")
          nanam <- paste0(grpnamc, c("", ".u", ".l"))
          restab[which(restab[, nnam] < nth), nanam] <- NA
        }
        restab[]
        # restabl[[namc]] <- restab
        restab.num <- as.matrix(restab)
        
        suppressWarnings(mode(restab.num) <- "double")
        restab.denom <- restab.num
        restab.num[] <- ifelse(is.finite(restab.num), restab.num, 0)
        restab.denom[] <- ifelse(is.finite(restab.num), 1, 0)
        restab.num.sum <- restab.num.sum + restab.num
        restab.num.sumsq <- restab.num.sumsq + (restab.num)^2
        restab.denom.sum <- restab.denom.sum + restab.denom
        nsuccess <- nsuccess + 1
      # }      
    })
  }
  nc <- restab.denom.sum
  mnc <- restab.num.sum / nc
  sec <- sqrt((restab.num.sumsq / nc - mnc^2) * nc / (nc - 1)) / sqrt(nc)
  lc <- ifelse(as.numeric(mnc) < 0, mnc - 2 * sec, pmax(0, mnc - 2 * sec))
  # uc <- pmin(1, mnc + 2 * sec)
  uc <- mnc + 2 * sec
  nrnd <- 3
  mnc.pr <- round(mnc, nrnd)# * 100
  lc.pr <- round(lc, nrnd)# * 100
  uc.pr <- round(uc, nrnd)# * 100
  restab.out <- mnc
  restab.out[] <- paste0(mnc.pr, " (", lc.pr, ", ", uc.pr, ")")
  restab.out <- as.data.frame(restab.out)
  # restab.out.mn <- 
  # restab.out.se <- sqrt((restab.num.sumsq / restab.denom.sum - restab.out.mn^2) * 
  #   (restab.denom.sum / (restab.denom.sum - 1)))
  # as.data.frame(round(restab.num.sum / restab.denom.sum, 3) * 100)
  cols.do.not.agg <- c("meth", "err.rate.meth", "calibrate", "sep.imp.thresh", "test.stat", "use.upper")
  if("restab" %in% ls())  
    restab.out[, cols.do.not.agg] <- restab[, cols.do.not.agg]
  restab.out <- restab.out[, order(!colnames(restab.out) %in% cols.do.not.agg)]
  restabl[[namc]] <- restab.out
  print(runtab[scen, ])
  print(paste0("nsuccess = ", nsuccess))
  restab.show <- restabl[[namc]]#[!restabl[[scen]]$sep.imp.thresh & restabl[[scen]]$test.stat == "z", ]
  # restab.show <- restabl[[scen]][restabl[[scen]]$test.stat == "z", ]
  # print(restab.show[, !colnames(restab.show) %in% c(outer(grpnam, c(".l", ".u"), paste0))])
  print(restab.show[, !colnames(restab.show) %in% c(outer(grpnam, c(".n", ".x", ".l", ".u"), paste0))])
  # print(restabl[[scen]])
  # print(restabl[[scen]][restabl[[scen]]$err.rate.meth == "lfsr", ])
}


methv


colnames(restabl[[1]])
str(restab.checkl)
pl <- sapply(restab.checkl, function(x) x$het.hom.post)

pl <- sapply(restab.checkl, function(x) x$ref.lines.postimpcomp)
table(pl[47, ])
apply(t(pl), 2, function(v) median(v, na.rm = T))

boxplot(t(pl))



options(max.print = 5000)








# var.in.name <- c("N", "P")
# for(scen in 1:nrow(runtab)){
#   namc <- paste(paste(var.in.name, runtab[scen, var.in.name], sep = "_"), collapse = "_")
#   bashscript <- paste0(bash.file.dir, "/test_script_", namc, ".sh")
#   cat(file = bashscript, append = F, sep = "\n",
#       c("#!/bin/bash",
#         "#SBATCH --ntasks=1",
#         "#SBATCH --output=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%x_%a.out",
#         "#SBATCH --error=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%x_%a.err",
#         "#SBATCH --partition debug",
#         paste0("#SBATCH --array=1-", runtab$ntimes[scen]),
#         paste0("#SBATCH --job-name=", namc),
#         paste0("#SBATCH --mem-per-cpu=", runtab$mem[scen], "MB"),
#         paste("Rscript --no-restore --no-save /mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/run_EM_algorithm_temp.R", 
#               runtab$N[scen], runtab$P[scen], sep = " ")))
#   print(readLines(bashscript))
#   system(paste("sbatch", bashscript))
# }

#paste0("#SBATCH --ncpus-per-task=", runtab$ntimes[scen]),
#     ))
# namc <- paste(runtab[scen, ], collapse = "_")
# cat(file = bashscript, append = T, sep = " ",
#     c("srun",
#       paste0("--job-name=", namc),
#       paste0("--mem-per-cpu=", runtab$mem[scen], "MB"),
#       paste("Rscript --no-restore --no-save /mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/run_EM_algorithm_temp.R", 
#       runtab$N, runtab$P, sep = " ")))
# cat(file = bashscript, append = T, "\n")


# "#SBATCH --output=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%A_%a.out",
# "#SBATCH --error=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%A_%a.err",
# "#SBATCH --output=/mnt/c/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%j.out",
# "#SBATCH --error=/mnt/c/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%j.err",

























































# "#SBATCH --ntasks-per-core=1 --overcommit",
# "/mnt/c/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/R_files/impc_mv_analysis/run_EM_algorithm_temp.R",
# 
# 
# for(scen in 1:nrow(runtab)){#scen <- 1#
#   ntimes <- runtab$ntimes[scen]
#   namc <- paste(c(paste(var.in.name, runtab[scen, var.in.name], sep = "_")), collapse = "_")
#   bashscript <- paste0(bash.file.dir, "/test_script_", namc, ".sh")
#   if((!file.exists(paste0(meth.comp.output.dir, "/", file.out))) | replace.files){
#     cat(file = bashscript, append = F, sep = "\n",
#         c("#!/bin/bash",
#           "#SBATCH --ntasks=1",
#           "#SBATCH --output=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%j.out",
#           "#SBATCH --error=/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/testing/test_%j.err",
#           "#SBATCH --partition debug",
#           paste0("#SBATCH --job-name=", namc),
#           paste0("#SBATCH --mem-per-cpu=", runtab$mem[scen], "MB")))
#     for(seedc in 1:ntimes){#seedc <- 1#
#       file.out <- paste0(prefix, namc, "_seed_", seedc, ".RData")
#       cat(file = bashscript, append = T, sep = "\n",
#           paste("srun Rscript --no-restore --no-save",
#               "/mnt/x/projects/impc_mv_analysis/R_files/impc_mv_analysis/run_EM_algorithm_temp_sims.R",
#               seedc, runtab$N[scen], runtab$P[scen], runtab$bo[scen], runtab$si[scen], 
#               file.out, "--profile=task", sep = " "))
#     }
#     system(paste("sbatch", bashscript))
#   }
# }
# add.all.data.job <- F
# if(add.all.data.job){
#   runtab.add <- data.frame(N = 4548, P = 148, bo = 0, si = 0)
#   runtab <- rbind(runtab.add, runtab)
# }
# runtab <- runtab[!(runtab$bo == 1 & (runtab$P > 20 | runtab$N > 500)) & 
#                    !(runtab$si == 1 & (runtab$P > 20 | runtab$N > 500)), ]
# runtab <- runtab[runtab$bo == runtab$si, ]
# prefix <- c("test_server", "K_auto_restab", "phcen_fdr_em_repairing_restab", "cen_line_fdr_repairing_restab", "repairing_restab", "revamped_retab", "redu_rank_omega_pinull_cenlinefdr_restab", "omega_pinull_cenlinefdr_restab", "omega_pinull_linefdr_restab", "omega_is_back_w_pinull_restab", "omega_is_back__restab", 
#             "random2_start_upto_5_full_rank_restab", "robust2_upto_5_full_rank_restab", "rewritten_upto_5_full_rank_restab", 
#             "single_full_rank_restab", "back_to_basics_restab", "R_trates_ident_restab", 
#             "null_atom_no_om_w_fullrank_restab", "null_atom_w_om_restab", "null_atom_restab", "rank1_mix_restab", "newest_mix_restab","newer_mix_restab", "new_mix_restab", 
#             "w_calib_restab_", "w_more_hethom_restab_", "deeper_restab_", "deepest_restab_")[1]
# # replace.files <- T
# t.ago.start <- 0
# t.ago.end <- Inf
# remove.earlier.files <- F
# if(remove.earlier.files){
#   fvc <- list.files(meth.comp.output.dir, full.names = T)
#   finfoc <- file.info(fvc)
#   tdiff <- as.numeric(difftime(Sys.time(), finfoc$ctime, units = "hours"))
#   hist(tdiff)
#   abline(v = c(t.ago.start, t.ago.end), col = "red")
#   fv.delete <- fvc[tdiff > t.ago.start & tdiff < t.ago.end]
#   file.remove(fv.delete)
# }
# 
# ####################################################
# # General parameters to pass to run_EM.R
# daseq <- c("impc", "eqtl")[1:2]
# methodseq <- c("em.fit", "ed", "mash", "em.fit.rand", "em.fit.N.500", 
#                "em.test", "em.loocv", "em.fac", "em.all.test", "em.all")[1:5]
# # Nseq <- c(100, 200, 500, 1000, 2000, 5000)[3]
# # Pseq <- c(10, 20, 40, 60, 100, 148)[6]
# Nseq <- Pseq <- NA
# full.analysis.seq <- c(T)
# sexspecific.seq <- F#c(F, T)
# nSig.seq <- 1:2
# n.seed.run.current <- 10
# max.test.seq <- 20000
# seedc <- 3
# EMfm.seq <- c("fa", "pca")[1]
# EMmash.seq <- F#c(T, F)
# EMbic.seq <- c(0, .5, 1, 2, 3)[1]
# EMwish.seq <- F#c(0, .25, .5)
# EMK.seq <- c(15, 20, 30, 40)
# K.use <- 20
# EMKup.seq <- F
# 
# ###############################################
# # Our EM specific parameters
# EMtol.seq <- c(1e-4, 1e-5, 5e-6)[1]
# 
# ###############################
# # MASH specific parameters
# bovyseq <- 1
# singleseq <- 1
# mash.EDtol <- 1e-5
# 
# ########################################
# # Extreme Deconvolution specific  parameters
# EDmeth.seq <- c("mash", "justED")[2]
# EDtol.seq <- c(1e-3, 1e-4, 1e-5, 1e-6)[2]
# 
# #############################
# # Create runtab
# runtab <- expand.grid(N = Nseq, P = Pseq, bo = bovyseq, si = singleseq, da = daseq, 
#                       me = methodseq, ma = max.test.seq, fu = full.analysis.seq, 
#                       ss = sexspecific.seq, nSig = nSig.seq, EDmeth = EDmeth.seq, 
#                       EDtol = EDtol.seq, EMtol = EMtol.seq, EMfm = EMfm.seq, 
#                       EMbic = EMbic.seq, EMwish = EMwish.seq, EMmash = EMmash.seq, 
#                       EMK = EMK.seq, EMKup = EMKup.seq, stringsAsFactors = F)
# runtab[!grepl("em", runtab$me), c("EMtol", "EMfm", "EMbic", "EMmash", "EMK", "EMKup")] <- NA
# # runtab[which(runtab$EMbic == 0), c("EMfm")] <- "pca"
# # runtab[!grepl("em", runtab$me), c("EMfm")] <- NA
# runtab[!runtab$me %in% c("ed", "mash"), c("EDtol", "EDmeth")] <- NA
# runtab[runtab$me != "mash", c("bo", "si")] <- NA
# runtab[runtab$me == "mash", "EDmeth"] <- "mash"
# runtab[runtab$me == "mash", "EDtol"] <- mash.EDtol
# runtab[runtab$me == "ed" & runtab$EDmeth == "mash", "EDtol"] <- mash.EDtol
# runtab[runtab$me == "ed" & runtab$EDmeth == "mash", c("nSig")] <- 1
# runtab[runtab$me == "mash", c("nSig")] <- 1
# runtab[!(runtab$me == "em.fit" & runtab$nSig == 1), "fu"] <- T
# runtab[runtab$me == "em.fit" & runtab$nSig == 1, "EMtol"] <- 1e-4#5e-6
# runtab[runtab$me == "em.fit" & !runtab$EMKup, "EMbic"] <- 0
# 
# # runtab <- runtab[runtab$N > 3 * runtab$P, ]
# 
# runtab[runtab$fu & runtab$da == "impc", "P"] <- 148
# runtab[runtab$fu & runtab$da == "impc", "N"] <- 2000
# runtab[runtab$fu & runtab$da == "eqtl", "P"] <- 44
# runtab[runtab$fu & runtab$da == "eqtl", "N"] <- 5000
# runtab <- unique(runtab)
# runtab$ntimes <- n.seed.run.current
# runtab[runtab$me == "em.fit.N.500", "N"] <- 500
# runtab <- runtab[!(runtab$me == "em.fit.N.500" & (runtab$da == "eqtl" | runtab$nSig > 1)), ]
# runtab <- runtab[!(runtab$me == "em.fit.rand" & (runtab$da == "eqtl")), ]
# runtab <- runtab[!(runtab$me %in% c("em.fit.rand", "em.fit.N.500") & runtab$EMK != K.use), ]
# runtab <- runtab[!(runtab$me == "em.fit" & runtab$nSig > 1 & runtab$EMK != K.use), ]
# runtab <- runtab[!(runtab$me == "em.fit.rand" & runtab$nSig > 1), ]
# 
# 
# ###############################################################
# #Specify memory requirements for each type of job
# runtab$mem <- 2300
# runtab[runtab$da %in% c("impc", "eqtl") & runtab$me == "em.fit", "mem"] <- 2000
# runtab[runtab$da  == "impc" & runtab$me == "mash" & runtab$fu, "mem"] <- 6000
# runtab[runtab$da == "eqtl" & runtab$me == "mash" & runtab$fu, "mem"] <- 4600
# 
# runtab <- unique(runtab)
# # runtab <- runtab[order(runtab$da, runtab$me, runtab$nSig), ]
# runtab <- runtab[order(runtab$da, match(runtab$me, methodseq), runtab$nSig), ]
# runtab
# runtab$rand <- ifelse(grepl("rand", runtab$me), T, F)
# runtab$loocv <- ifelse(runtab$me == "em.fit" & runtab$nSig == 1 & runtab$N == 2000, T, F)
# rownames(runtab) <- 1:nrow(runtab)
# # runtab$meth <- 
# fcheckfields <- c("emout.ok", "emout.st", "emout.en", "res.ok", "res.st", "res.en", "loocv.ok", 
#                   "loocv.st", "loocv.en", "fac.ok", "fac.st", "fac.en")
# runtab[, fcheckfields] <- NA
# runtab[, c("file.base", "meth")] <- NA
