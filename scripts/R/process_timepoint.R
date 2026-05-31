#!/usr/bin/env Rscript
# Usage: Rscript process_timepoint.R week3
suppressPackageStartupMessages({
library(yaml)
})

if (interactive()) {
  # Manually set these for running in R Studio
  selected_timepoint <- "3_weeks" # e.g., 3_weeks, 7_weeks, 18.4_week_dams, 7_months, 18_months
  base_dir           <- "C:/Users/David/Documents/Research/PhD/scRNAseq/" # update to your base directory
  config_path        <- "configs/scrnaseq_qc_config.yaml" # update to your config path
} else {
  # Command line logic
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 3) {
    stop("Usage: Rscript process_timepoint.R <timepoint> <base_dir> <config_path>", call. = FALSE)
  }
  selected_timepoint <- args[1]
  base_dir           <- args[2] # e.g. "C:/Users/Name/Documents/Research/PhD/scRNAseq/"
  config_path        <- args[3] # e.g. "configs/config.yaml"
}

# format base_dir to ennsure it doesn't end with a slash, to avoid issues with file.path later
if (endsWith(base_dir, "/")) base_dir <- sub("/$", "", base_dir)
if (!file.exists(base_dir)) stop(sprintf("Cannot find base directory at %s", base_dir), call. = FALSE)

configs <- read_yaml(config_path)

# Use valid_timepoints from configs
if(!is.null(configs$valid_timepoints)) {
  valid_timepoints <- configs$valid_timepoints
} else {
  stop("Configuration must define valid_timepoints", call. = FALSE)
}

# ensure selected_timepoint is valid
if (!selected_timepoint %in% valid_timepoints) stop(sprintf("Error: '%s' is not a valid timepoint.\nAllowed values are: %s", 
               selected_timepoint, 
               paste(configs$valid_timepoints, collapse = ", ")), 
               call. = FALSE)

# load libraries
suppressPackageStartupMessages({
library(Seurat)
library(SeuratObject)
library(scCustomize)
library(tidyverse)
library(scDblFinder)
library(DropletUtils)
library(rhdf5)
library(hdf5r)
library(SoupX)
library(qs2) # for fast data saving/loading
library(plyr)
library(zeallot) # for unpacking multiple values from functions
library(scales)
library(BiocParallel)
})

# source utility functions
source(here::here("scripts/R/utils.R"))

# Apply the global settings
set.seed(configs$random_seed)
options(future.globals.maxSize <- configs$future_globals_maxSize)
nthreads <- parallel::detectCores(logical=FALSE) - 2 # leave two physical cores free for system processes
c(nthreads, bp_param) %<-% configure_parallelism(rng_seed = configs$random_seed)

# Directories (Update these to match your exact environment)
# base_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/"
# base_dir_external <- "D:/PhD/scRNAseq/"
# base_dir_turbo <- "//sph-colacino-win.turbo.storage.umich.edu/sph-colacino/aguilada/scRNAseq"
# cellbender_fulloutput_dir <- paste0(base_dir,"cellbender/cellbender_h5_outputs/FPR_0.0/")
# cellbender_dir <- paste0(base_dir, "cellbender/cellbender_h5_outputs_for_seurat/FPR_0.0/")
# cellranger_dir <- ifelse(selected_timepoint == 'week3',
#                            paste0(base_dir,"8651-JC-mammary/10x_analysis_JC-8651_cellranger_v.9.0.0/"),
#                            paste0(base_dir,"10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/")
#                           )

# analysis type, cellbender or cellranger ----
analysis = configs$analysis$cellbender # e.g., "cellbender_analysis" or "cellranger.v.9.0.0_analysis" for creating appropriate folders
# analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analysis' # we set emptydrops-minimum-umis = 100, and default values for all else

# Directories where we'll save our generated data and QC figures
data_dir <- file.path(base_dir, "analysis", selected_timepoint, analysis, "FPR_0.0", "data")
fig_dir  <- file.path(base_dir, "analysis", selected_timepoint, analysis, "FPR_0.0", "qc")

dirs_to_make <- c(
  data_dir, 
  fig_dir,
  file.path(data_dir, "emptyDrops"),
  file.path(fig_dir, "emptyDrops"),
  file.path(data_dir, "soupx"),
  file.path(fig_dir, "soupx"),
  file.path(data_dir, "cellbender_ambient_profiles")
)

walk(dirs_to_make, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))

# Initialize empty lists to store Seurat objects
cb.fltrd.list <- list() # for cellbender data filtered by emptyDrops calls
cb.list <- list() # for cellbender data. NO emptyDrops filtering.
raw.list <- list() # for cellranger raw data
fltrd.list <- list() # for cellranger filtered data
joined.cb.list <- list() # for joined cellbender and raw cellranger data (to validate cellbender on same barcodes)
soupx.list <- list() # for soupx filtered cell ranger matrices
joined.sx.list <- list() # for joined cellranger w/ soupx and raw cellranger data (to validate soupx on same barcodes)

# load the sample sheet, assuming it's in the configs directory
if (!file.exists('configs/sample_sheet.csv')) stop(sprintf("Cannot find 'configs/sample_sheet.csv'), make sure it is in the working directory."), call. = FALSE)
samples <- read.csv(file.path('configs/sample_sheet.csv'),header=TRUE)

# subset the sample sheet to only include rows for the selected time point
samples_to_process <- samples[samples$timepoint==selected_timepoint,]

cat("=================================================\n")
cat("Processing data for:", selected_timepoint, "\n")
cat("=================================================\n")

# --- Core Processing Loop ---
processed_seurat <- lapply(seq_along(nrow(samples_to_process)), function(row_idx) {
  row <- samples_to_process[row_idx, ]

  cat("\nProcessing ", selected_timepoint, " Sample:", row$sample_id, "\n")
  
  # path  if pool directory is specified, use that, otherwise use project directory
  raw_pool_file_path <- file.path(base_dir, row$cellranger_dir, row$pool_dir, row$core_sample_id, row$cellranger_raw_h5_filename)
  fltrd_pool_file_path <- file.path(base_dir, row$cellranger_dir, row$pool_dir, row$core_sample_id, row$cellranger_filtered_h5_filename)

  # path without pool
  raw_file_path <- file.path(base_dir, row$cellranger_dir, row$core_sample_id, row$cellranger_raw_h5_filename)
  fltrd_file_path <- file.path(base_dir, row$cellranger_dir, row$core_sample_id, row$cellranger_filtered_h5_filename)

  # cellbender path
  cb_file_path <- file.path(base_dir, row$cellbender_seurat_h5_dir, row$fpr_dir, row$cellbender_seurat_h5_filename)
  # cellbender full output path (for metadata)
  cb_full_h5_path <- file.path(base_dir, row$cellbender_h5_dir, row$fpr_dir, row$cellbender_h5_filename)

  # Resolve which paths actually exists on this file system
  if (file.exists(raw_pool_file_path) && file.exists(fltrd_pool_file_path)) {
    raw_file_path <- raw_pool_file_path
    fltrd_file_path <- fltrd_pool_file_path
  } else if (file.exists(raw_file_path)) {
    # raw_file_path is already set to the non-pool path, so we can just check if it exists
    if (!file.exists(fltrd_file_path)) {
      stop(sprintf("\nCRITICAL ERROR: Cannot find filtered matrix for %s.\nChecked Location 1: %s\nChecked Location 2: %s", row$sample_id, fltrd_pool_file_path, fltrd_file_path))
    }
  } else {
    stop(sprintf("\nCRITICAL ERROR: Cannot find raw matrix for %s.\nChecked Location 1: %s\nChecked Location 2: %s", row$sample_id, raw_pool_file_path, raw_file_path))
  }

  if (!file.exists(cb_file_path)) stop(paste("Cannot find CellBender matrix at:", cb_file_path))
  if (!file.exists(cb_full_h5_path)) stop(paste("Cannot find CellBender full output matrix at:", cb_full_h5_path))


  ## Read the counts ----
  cat('\nloading from',raw_file_path, 'for ',selected_timepoint,'sample:',row$core_sample_id,'\n')
  raw_matrix <- Read10X_h5(filename=raw_file_path)

  cat('\nloading from',fltrd_file_path, 'for',selected_timepoint,'sample:',row$core_sample_id,'\n')
  fltrd_matrix <- Read10X_h5(filename=fltrd_file_path)

  cat('\nloading from',cb_file_path, 'for',selected_timepoint,'sample:',row$core_sample_id,'\n')
  cb_matrix <- Read10X_h5(filename=cb_file_path)
  
  ## run EmptyDrops on raw sample ----
  # set lower UMI bound based on yaml config, with timepoint-specific override if it exists.
  lower = get_param(configs, "emptydrops", "lower_umi_bound", selected_timepoint)

  cat("---------------------------------------------\n")
  cat('\nRunning EmptyDrops for',selected_timepoint,'sample:',row$core_sample_id,'with lower = ', lower, '\n')
  cat("---------------------------------------------\n")

  e.out <- emptyDrops(raw_matrix,
                        lower= lower, # lower bound of total UMIs. default = 100
                        test.ambient=TRUE,
                        niters=10000, # n iterations for monte carlo p-value calculation
                        retain=Inf,
                        BPPARAM = bp_param
  )

  gc() # clean up memory after emptyDrops, which can be intensive on large matrices

  fdr_cutoff = configs$emptydrops$default$fdr_cutoff # start at 0.01, can adjust later based on results
  is.cell = e.out$FDR <= fdr_cutoff
  e.out$is.cell = is.cell # droplets with FDR <= fdr_cutoff considered non-empty. Can adjust later with df.
  cat("Number of cells detected by EmptyDrops at FDR <=0.01:", sum(is.cell, na.rm = TRUE), "\n",
      "Number of limited droplets:", sum(e.out$Limited, na.rm = TRUE), "\n")
  print(table(Limited=e.out$Limited, Significant=is.cell))

  e.calls <- data.frame(
    row.names = rownames(e.out),
    emptydrops_FDR = e.out$FDR, # NAs included for droplets that didn't pass with < lower UMI bound.
    emptydrops_logprob = e.out$LogProb,
    emptydrops_pvalue = e.out$PValue,
    emptydrops_Limited = e.out$Limited,
    emptydrops_is_cell = is.cell
  )

  write.csv(e.calls, file.path(data_dir, "emptyDrops", paste0(selected_timepoint, "_sample_", row$core_sample_id, "_emptyDrops_calls.csv")), row.names = TRUE)

  ## histogram of emptyDrops null p-values ----
  cat('\nPlotting EmptyDrops graphs for sample:',selected_timepoint,'sample:',row$core_sample_id,'\n')
  png(filename = paste0(fig_dir,'emptyDrops/',selected_timepoint,'_sample_',row$core_sample_id,'_histogram_emptyDrops_null_pvalues.png'), units='in', width=10, height=10, res=180)
  hist(e.out[which(e.out$Total <= configs$emptydrops$default$lower_umi_bound & e.out$Total > 0),]$PValue,breaks = 20,
        main = 'EmptyDrops null hypothesis p-values',
        sub = paste0('(cells < ', configs$emptydrops$default$lower_umi_bound, ' UMI counts)'),
        xlab='p-value')
  dev.off()
  
  ## Create emptyDrops classification plot ----
  ggplot(e.out, aes(x = Total, y = -LogProb)) +
    # Add points with color based on cell classification
    geom_jitter(aes(color = is.cell), alpha = 0.6, size = 1.2) +
    
    # Use log10 scale for x-axis with comma formatting for readability
    scale_x_log10(labels = comma_format(), n.breaks = 6) +
    scale_y_log10(labels = comma_format(), n.breaks = 6) +
    annotation_logticks(color='black') +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("TRUE" = "Cell", "FALSE" = "Background"),
      name = "Classification"
    ) +
    labs(
      title = "Cell Classification by UMI Count and Probability",
      subtitle = paste0("Cells with < ", configs$emptydrops$default$lower_umi_bound, " UMI counts considered empty"),
      x = "Total UMI count (log scale)",
      y = "-Log Probability"
    ) +
    
    # Add a vertical reference line at a threshold value
    geom_vline(xintercept = configs$emptydrops$default$lower_umi_bound,
                linetype = "dashed", color = "blue", alpha = 0.7) +
    theme_classic2() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  ggsave(paste0(fig_dir,'emptyDrops/',selected_timepoint,'_sample_',row$core_sample_id,'_emptyDrops_cell_calls.png'), units='in', width=10, height=10, dpi=320)
    
  
  ## Create Seurat objects ----
  cat('\ncreating seurat objects for ',selected_timepoint,' sample:',row$core_sample_id,'\n')
  
  raw_seurat_obj <- CreateSeuratObject(
    counts  = raw_matrix,
    assay   = "RNA",
    project = paste0(row$core_sample_id, "_raw")
  )
  
  # joint seurat object with cellbender and raw cell ranger counts for cellbender QC.
  joined_cb_and_raw <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_matrix,
                                            raw_counts_matrix = raw_matrix)
  joined_cb_and_raw <- AddMetaData(joined_cb_and_raw, metadata=e.calls) # add emptyDrops calls to metadata for later comparison
  
  fltrd_seurat_obj <- CreateSeuratObject(
    counts  = fltrd_matrix,
    assay   = "RNA",
    project = paste0("sample", row$core_sample_id, "_filtered")
  )

  cb_seurat_obj <- CreateSeuratObject(
    counts  = cb_matrix,
    assay   = "RNA",
    project = paste0("sample", row$core_sample_id, "_cellbender")
  )
  cb_seurat_obj <- AddMetaData(cb_seurat_obj, metadata=e.calls) # add emptyDrops calls to metadata

  # read in cellbender metadata, i.e. cell probabilites, cell size, 
  # droplet efficiency (reverse transcription efficiency within the droplet), 
  # background fraction
  
  cat('\nMerging cellbedender metadata with seurat objects for sample:',row$core_sample_id,'\n')
  cellbender_h5_output_path <- file.path(base_dir, row$cellbender_h5_dir, row$fpr_dir, row$cellbender_h5_filename)

  ## grab cellbender metadata from hdf5 file ----
  metadata <- h5read(cellbender_h5_output_path, 'metadata')
  droplet_latents <- h5read(cellbender_h5_output_path, 'droplet_latents')
  droplet_latents <- as.data.frame(droplet_latents[-c(2,6)]) # don't want gene_expression encoding or barcode indices
  rownames(droplet_latents) <- metadata$barcodes_analyzed
  
  ## grab cellbender ambient profile
  matrices <- h5read(cellbender_h5_output_path, 'matrix')
  global_latents <- h5read(cellbender_h5_output_path, 'global_latents')
  global_latents <- as.data.frame(global_latents)[1]
  global_latents$gene <- matrices$features$name
  global_latents$ensemble_ids <- matrices$features$id
  global_latents <- global_latents[,c(2,1,3)]
  
  # save ambient profile
  write.csv(global_latents, paste0(data_dir,'cellbender_ambient_profiles/',selected_timepoint,'_sample_',row$core_sample_id,'_ambient_expression.csv'), row.names = FALSE)

  # add cellbender metadata to all seurat objects
  raw_seurat_obj <- AddMetaData(raw_seurat_obj, metadata=droplet_latents)
  joined_cb_and_raw <- AddMetaData(joined_cb_and_raw, metadata=droplet_latents)
  fltrd_seurat_obj <- AddMetaData(fltrd_seurat_obj, metadata=droplet_latents)
  cb_seurat_obj <- AddMetaData(cb_seurat_obj, metadata=droplet_latents)

})


# --- Merging ---
cat("\nMerging Seurat Objects...\n")
merged_cb <- merge(processed_data[[1]], processed_data[-1], add.cell.ids = paste0("tp_", samples_to_process)) %>% 
  JoinLayers()

# Add Custom Metadata
merged_cb$cells <- rownames(merged_cb@meta.data)
merged_cb$sample <- str_extract(merged_cb$cells, '^[^_]+_[^_]+')

# --- Doublet Detection ---
cat("\nRunning Doublet Detection on Unfiltered Merged Data...\n")
sce_cb <- as.SingleCellExperiment(merged_cb)
sce_cb <- scDblFinder(sce_cb, samples = "sample", multiSampleMode = 'split', clusters = TRUE, 
                      dbr.sd = 1, dbr.per1k = 0.008, BPPARAM = SnowParam(workers = nthreads))

merged_cb <- AddMetaData(merged_cb, metadata = as.data.frame(colData(sce_cb)) %>% dplyr::select(scDblFinder.class, scDblFinder.score))

# --- Save Checkpoints ---
cat("\nSaving Unfiltered Checkpoints...\n")
qs_save(merged_cb, paste0(data_dir, selected_timepoint, "_unfiltered_checkpoint.qs"), nthreads = nthreads)
qs_save(merged_cb@meta.data, paste0(data_dir, selected_timepoint, "_unfiltered_metadata.qs"), nthreads = nthreads)

cat("\nCheckpoint Generation Complete!\n")