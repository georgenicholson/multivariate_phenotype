
source("X:/projects/impc_mv_analysis/R_files/impc_mv_parameters.R")

source(paste(R.file.dir, "/impc_import_data/impc_import_data.R", sep = ""))
source(paste(R.file.dir, "/impc_import_data/create_negative_controls.R", sep = ""))
#note this file is quite heavy on RAM in parallel processing stage...
#If you have time you could recode to only pass subsets of d to slave cores
source(paste(R.file.dir, "/impc_uv_analysis/univariate_analysis.R", sep = ""))
source(paste(R.file.dir, "/impc_uv_analysis/add_sig_thresh_to_uv_results.R", sep = ""))
source(paste(R.file.dir, "/impc_uv_analysis/identify_outlying_uv.R", sep = ""))
source(paste(R.file.dir, "/impc_uv_analysis/qc+scale_uv_data.R", sep = ""))

# use.resmat.qc <- T
# source("X:/projects/impc_mv_analysis/R_files/impc_uv_analysis/add_sig_thresh_to_uv_results.R")
# source(paste(R.file.dir, "/impc_mv_analysis/run_EM_algorithm_temp.R", sep = ""))
source(paste(R.file.dir, "/impc_mv_paper_code/run_EM.R", sep = ""))
source(paste(R.file.dir, "/impc_mv_paper_code/collect_results.R", sep = ""))
source(paste(R.file.dir, "/impc_mv_paper_code/calculate_hit_rates.R", sep = ""))

source(paste(R.file.dir, "/impc_mv_paper_plots_code/table_hit_rate_err_rate_impc.R", sep = ""))
# hitrates_comb.txt: table comparing hit rates and error rates
# cvlik_table.txt: CV lik comparison for IMPC

source(paste(R.file.dir, "/impc_mv_paper_plots_code/plot_power_comparisons_revision.R", sep = ""))
# combined_power_plot.jpg: multipanel power plot

# source(paste(R.file.dir, "/impc_mv_analysis/process_EM_output.R", sep = ""))
# em_convergence.jpg: EM convergence plot
# source(paste(R.file.dir, "/impc_mv_analysis/add_sig_thresholds_to_uv+mv_results.R", sep = ""))
source(paste(R.file.dir, "/impc_mv_paper_code/factor_model.R", sep = ""))

# source(paste(R.file.dir, "/impc_plots/factor_interpretation_plots.R", sep = ""))
# factor_interpretation_plot.jpg: multipanel factor interpretation plot
# cumulative_correlation_explained.jpg: Cumulative proportion of correlation structure explained by eigenvectors

source(paste(R.file.dir, "/impc_plots/plot_power_comparisons.R", sep = ""))
# uv_mv_scatter.jpg: supplementary figure examining concordance between UV and MV

source(paste(R.file.dir, "/impc_plots/plot_ref_lines.R", sep = ""))
# ref_lines_heatmap.jpg: Replicability heat map
# ref_lines_z_scatter.jpg: replicability scatterplots
# hethom_scatter.jpg: Heterozygote/homozygote concordance scatterplot

source(paste(R.file.dir, "/impc_plots/plot_annotation_heatmaps.R", sep = ""))
# paper_correlation_heatmaps_figure.jpg: Paper heat maps
# paper_correlation_heatmaps_figure.jpg: Paper covariance matrices

source(paste(R.file.dir, "/impc_plots/post_eb_assess_imputation.R", sep = ""))
# looEB_uv_scatter.jpg: Scatterplots examining concordance of LOO-MV with UV/MV
# power_MV_looMV_UV_by_phenotype.jpg: The phenotype-specific proportion of lines annotated

source(paste(R.file.dir, "/impc_plots/plot_impc_design_cartoon.R", sep = ""))
# impc_design_cartoon.jpc: Cartoon for paper








# line_by_line_uv_mv_comp.jpg
# fnamc <- paste("estmeth_", estimation.meth, "_nIts_", mv.nits, "_annotation_heatmap_justSig_", 
# justSig, "_compuvmv_", incl, ".jpg", sep = "")
# power_MV_looMV_UV.jpg: The procedure-specific proportion of lines annotated

library("jpeg")
library("tiff")
library("png")
install.packages("png")
img <- readTIFF("origin.tiff", native=TRUE)
writeJPEG(img, target = "Converted.jpeg", quality = 1)



dropfigdir <- "C:/Users/nicho/Dropbox/impc_mv_paper/tex_files/figures"
plfv <- list.files(dropfigdir)
plfv <- plfv[grepl(".jpg", plfv) | grepl(".png", plfv)]

figmap <- data.frame(file = plfv, fig = NA, supp = NA)
figmap[figmap$file == "impc_design_cartoon.jpg", c("fig", "supp")] <- c(1, 0)
figmap[figmap$file == "paper_annotation_heatmaps_figure.jpg", c("fig", "supp")] <- c(2, 0)
figmap[figmap$file == "paper_correlation_heatmaps_figure.jpg", c("fig", "supp")] <- c(3, 0)
figmap[figmap$file == "combined_power_plot.jpg", c("fig", "supp")] <- c(4, 0)
figmap[figmap$file == "ref_lines_heatmap.jpg", c("fig", "supp")] <- c(5, 0)
figmap[figmap$file == "ref_lines_z_scatter.jpg", c("fig", "supp")] <- c(6, 0)
figmap[figmap$file == "factor_interpretation_plot.jpg", c("fig", "supp")] <- c(7, 0)

figmap[figmap$file == "impc_pipeline_figure.png", c("fig", "supp")] <- c(1, 1)
figmap[figmap$file == "uv_mv_scatter.jpg", c("fig", "supp")] <- c(2, 1)
figmap[figmap$file == "looEB_uv_scatter.jpg", c("fig", "supp")] <- c(3, 1)
figmap[figmap$file == "power_MV_looMV_UV_by_phenotype.jpg", c("fig", "supp")] <- c(4, 1)
figmap[figmap$file == "hethom_scatter.jpg", c("fig", "supp")] <- c(5, 1)
figmap[figmap$file == "cumulative_correlation_explained.jpg", c("fig", "supp")] <- c(6, 1)
figmap[figmap$file == "em_convergence.jpg", c("fig", "supp")] <- c(7, 1)
figmap[figmap$file == "uv_qc_plot_for_paper.jpg", c("fig", "supp")] <- c(8, 1)

figmap <- figmap[!is.na(figmap$fig), ]
figmap$full.name <- paste0(ifelse(figmap$supp, "Supplementary ", ""), "Figure ", figmap$fig)
figmap$short.name <- paste0(ifelse(figmap$supp, "Supp", ""), "Fig", figmap$fig)
figmap <- figmap[order(figmap$supp, figmap$fig), ]



tiffoutdir <- paste0(dropfigdir, "/tiffForPaper")
jpgoutdir <- paste0(dropfigdir, "/jpgForPaper")
dir.create(tiffoutdir, showWarnings = F)
dir.create(jpgoutdir, showWarnings = F)
for(j in 1:nrow(figmap)){
  if(grepl("jpg", figmap$file[j]))
    img <- readJPEG(paste0(dropfigdir, "/", figmap$file[j]), native = TRUE)
  if(grepl("png", figmap$file[j]))
    img <- readPNG(paste0(dropfigdir, "/", figmap$file[j]), native = TRUE)
  # writeTIFF(img, where = paste0(tiffoutdir, "/", figmap$short.name[j], ".tif"))
  writeJPEG(img, target = paste0(jpgoutdir, "/", figmap$short.name[j], ".jpg"))
}


