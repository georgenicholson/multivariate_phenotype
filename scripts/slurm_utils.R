rm(list = ls())
##########################################
# Super parameters
run_type <- c("demo", "main", "benchmark", "test_benchmark")[1]

##########################################
# Source function files
fns_to_source <- list.files("scripts/functions", full.names = TRUE)
for (file_curr in fns_to_source) {
  source(file_curr)
}

##########################################
# control to contain parameters
control <- get_control_parameters_mv()

##########################################
# Create directory structure
for(dirc in c(control$output_dir, control$methods_comp_dir, control$global_res_dir, control$data_dir, control$train_test_samples_dir))
  dir.create(dirc, showWarnings = F)

# ##########################################
# # Download data
# if (!"data/impc" %in% list.dirs(control$data_dir)) {
#   temp <- tempfile()
#   download.file(url = "https://www.dropbox.com/s/wz6mg35n502au6a/multivariate_phenotype_paper_supporting_data.zip?dl=1", temp, mode = "wb")
#   unzip(temp, exdir = control$data_dir)
#   unlink(temp)
# }

##########################################
# Get table of analyses
analysis_table <- create_table_of_analyses(control = control, check_status = T, run_type = run_type)
# analysis_table <- analysis_table[order(analysis_table$Data, decreasing = TRUE), ]
scen_to_run <- 1:nrow(analysis_table)
# network_path_to_MVphen_git_repo <- "/mnt/x/projects/impc_mv_analysis/github_multivariate_phenotype/multivariate_phenotype"
network_path_to_MVphen_git_repo <- "/mnt/c/Users/nicho/Dropbox/GitHub_Dropbox/multivariate_phenotype"
# network_path_to_MVphen_git_repo <- "//TW/Users/nicho/Documents/bauer_sync/projects/impc_mv_analysis/github_multivariate_phenotype/multivariate_phenotype"
dir.create(".job", showWarnings = FALSE)
partc <- c("debug", "debug11", "debug6", "debug8", "debug4")[5]
file.copy(from = "scripts/01_model_fitting_wrapper.R",
          to = "scripts/01_model_fitting_wrapper_temp.R", overwrite = T)
for (scen in scen_to_run) {#scen <- 1#
  for (subsamseed in 1:analysis_table$n_subsamples[scen]) {#subsamseed <- 1#
    file_name_curr <- gsub("XXX", subsamseed, analysis_table$file_core_name[scen])
    bashscript <- paste0(".job/", file_name_curr, ".sh")
    output.file <- file(bashscript, "wb")
    # partc <- ifelse(analysis_table$mem[scen] > 3000, "debug4", "debug11")
    cat(file = output.file, append = F, sep = "\n",
        c("#!/bin/bash\n",
          "#SBATCH --ntasks=1", 
          paste0("#SBATCH --output=/mnt/s/job_output_files/", file_name_curr, ".out"),
          paste0("#SBATCH --error=/mnt/s/job_output_files/", file_name_curr, ".err"),
          paste0("#SBATCH --partition ", partc),
          paste0("#SBATCH --job-name=", file_name_curr),
          paste0("#SBATCH --mem-per-cpu=", analysis_table$mem[scen], "MB"),
          paste0("sudo mount -t drvfs //TW/Users/nicho/Documents/bauer_sync /mnt/x"),
          paste0("cd ", network_path_to_MVphen_git_repo),
          paste("srun Rscript --no-restore --no-save scripts/01_model_fitting_wrapper_temp.R",
                run_type, scen, subsamseed, sep = " ")))
    if(Sys.info()["sysname"] == "Windows"){
      # system(paste("wsl sudo dos2unix -n", bashscript, bashscript))
      system(paste("wsl sbatch", file.path(network_path_to_MVphen_git_repo, bashscript)))
    } else {
      system(paste("sbatch", file.path(network_path_to_MVphen_git_repo, bashscript)))
    }
    close(output.file)
  }
  # Sys.sleep(2)
}




# system(paste0("wsl cat ", "/mnt/s/job_output_files/", file_name_curr, ".err"))
# system(paste0("wsl cat ", "/mnt/s/job_output_files/", file_name_curr, ".out"))
# Data_impc_Meth_mash_N_100_P_10_subsamseed_1_nSig_1_MVphen_K_3_XDmeth_mash.err
# sinfo
# scancel -u root
 # squeue
# scontrol update NodeName=HPZ420-4 State=IDLE
# sacct
# sudo mount -t drvfs //HPZ420-0/Users/nicho/Documents/slurm /mnt/s
# sudo mount -t drvfs //TW/bauer_sync /mnt/x

# C:\Users\nicho\OneDrive\Documents\slurm_setup_files\shell_files


