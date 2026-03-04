#!/usr/bin/env Rscript
# Usage: Rscript process_timepoint.R week3
suppressPackageStartupMessages({
library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
valid_timepoints <- c("week3", "week10", "week18.4dam", "month7", "month18")
if (length(args) < 3) {
  stop("Usage: Rscript 01_process_checkpoint.R <timepoint> <base_dir> <config_path>", call. = FALSE)
}

selected_timepoint <- args[1]
base_dir <- args[2] # e.g. "C:/Users/Name/Documents/Research/PhD/scRNAseq/"
config_path <- args[3] # e.g. "config.yaml"

if (!selected_timepoint %in% valid_timepoints) stop(sprintf("Error: '%s' is not a valid timepoint.\nAllowed values are: %s", 
               selected_timepoint, 
               paste(valid_timepoints, collapse = ", ")), 
               call. = FALSE)
if (!endsWith(base_dir, "/")) base_dir <- paste0(base_dir, "/")

configs <- read_yaml(config_path)

# Apply the global settings
set.seed(configs$random_seed)
options(future.globals.maxSize = configs$future_globals_maxSize)
nthreads <- configs$nthreads

suppressPackageStartupMessages({
library(SeuratObject)
library(scCustomize)
library(tidyverse)
library(scDblFinder)
library(DropletUtils)
library(rhdf5)
library(SoupX)
library(qs2) # for fast data saving/loading
library(plyr)
library(BiocParallel)
})

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

## analysis type, cellbender or cellranger ----
# analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analysis' # we set emptydrops-minimum-umis = 100, and default values for all else

# Directories where we'll save our generated data and QC figures
data_dir <- file.path(base_dir, "analysis", selected_timepoint, analysis, "FPR_0.0_MAD", "data")
fig_dir  <- file.path(base_dir, "analysis", selected_timepoint, analysis, "FPR_0.0_MAD", "qc")

dirs_to_make <- c(
  data_dir, fig_dir,
  file.path(data_dir, "ambient_expression_profiles"),
  file.path(data_dir, "emptyDrops"),
  file.path(fig_dir, "emptyDrops"),
  file.path(data_dir, "soupx"),
  file.path(data_dir, "cellbender_ambient_profiles")
)

walk(dirs_to_make, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))


samples_to_process <- timepoints[[selected_timepoint]]

## Utility function to calculate metrics ----
calc_metrics = function(seurat.list) {
    
    lapply(seurat.list, function(x) {
    log10genes = log10(x$nFeature_RNA)
    log10counts = log10(x$nCount_RNA)
    
    # mean and standard deviation of log10 ngenes and ncount
    mean_ngene = mean(log10genes)
    mean_ncount = mean(log10counts)
    std_ngene = sd(log10genes)
    std_ncount = sd(log10counts)
    
    # z-scores of log10 ngenes and ncounts
    x$log10_genes_zscored = (log10genes - mean_ngene) / std_ngene 
    x$log10_counts_zscored = (log10counts - mean_ncount) / std_ncount
    
    # calculate log10 genes per log10 UMI for each cell 
    # acts as a novelty score to measure complexity
    x$complexity = log10genes/log10counts
    
    # Might as well do mitochondrial expression percentage here, along with ribosomal and hemoglobin
    x = PercentageFeatureSet(x, pattern = "^mt-", col.name = "percent.mt") # percent mitochondrial genes
    x = PercentageFeatureSet(x, pattern = "^Rp[sl]", col.name = "percent.rb") # percent ribosomal genes
    x = PercentageFeatureSet(x, pattern = "^Hb[^(P|E|S)]", col.name = "percent.hb") # percent hemoglobin genes
    
  })
}


cat("=================================================\n")
cat("Processing data for:", selected_timepoint, "\n")
cat("=================================================\n")


# Initialize empty lists to store Seurat objects
cb.fltrd.list <- list() # for cellbender data filtered by emptyDrops calls
cb.list <- list() # for cellbender data. NO emptyDrops filtering.
raw.list <- list() # for cellranger raw data
fltrd.list <- list() # for cellranger filtered data
joined.list <- list() # for joined cellbender and raw cellranger data (to validate cellbender on same barcodes)
soupx.list <- list() # for soupx filtered cell ranger matrices

# --- Core Processing Loop ---
processed_data <- lapply(samples_to_process, function(i) {
  cat("\nProcessing ", selected_timepoint, " Sample:", i, "\n")

  # create directories
  dirs = c(cellbender_dir, cellranger_dir, data_dir, fig_dir)
  for (dir in dirs) {
    if (!dir.exists(dir)) { dir.create(dir,
                                       recursive=TRUE) }
  }

  if (!dir.exists(paste0(data_dir,'ambient_expression_profiles/'))) {
    dir.create(paste0(data_dir,'ambient_expression_profiles/'),recursive=TRUE) }
  
  if (!dir.exists(paste0(data_dir,'emptyDrops/'))) {
    dir.create(paste0(data_dir,'emptyDrops/'),recursive=TRUE) }
  
  if (!dir.exists(paste0(fig_dir,'emptyDrops/'))) {
    dir.create(paste0(fig_dir,'emptyDrops/'),recursive=TRUE) }
  
  if (!dir.exists(paste0(data_dir,'soupx/'))) {
    dir.create(paste0(data_dir,'soupx/'),recursive=TRUE) }
  
  if (!dir.exists(paste0(data_dir,'cellbender_ambient_profiles/'))) {
    dir.create(paste0(data_dir,'cellbender_ambient_profiles/'),recursive=TRUE) }
  
  if (selected_timepoint == 'week3') {
      # Construct the file path for each sample directory
      raw_file_path <- paste0(cellranger_dir,
                              "8651-JC-", i, "/count/sample_raw_feature_bc_matrix.h5")
      
      fltrd_file_path <- paste0(cellranger_dir,
                                "8651-JC-", i, "/count/sample_filtered_feature_bc_matrix.h5")
      
      cb_file_path <- paste0(cellbender_dir,
                             "output_cellbender_8651-JC-", i, "_FPR_0.0_filtered_seurat.h5") # filtered version contains cells w/ >50% posterior probability of being a cell
      # "output_cellbender_8651-JC-", i, "_FPR_0.01_filtered_seurat.h5")
      
    } else { # load other time points from turbo
      if (cellranger_dir == "//sph-colacino-win.turbo.storage.umich.edu/sph-colacino/aguilada/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/") {
        raw_file_path <- paste0(cellranger_dir, timepoint_pools[[selected_timepoint]],
                                "11389-DA-", i, "/count/sample_raw_feature_bc_matrix.h5")
  
        fltrd_file_path <- paste0(cellranger_dir, timepoint_pools[[selected_timepoint]],
                                  "11389-DA-", i, "/count/sample_filtered_feature_bc_matrix.h5")
  
        cb_file_path <- paste0(cellbender_dir,
                               "output_cellbender_11389-DA-", i, "_FPR_0.0_filtered_seurat.h5")
        # "output_cellbender_8651-JC-", i, "_FPR_0.01_filtered_seurat.h5")
      } else { # otherwise load other time points locally
        raw_file_path <- paste0(cellranger_dir,
                                "11389-DA-", i, "/count/sample_raw_feature_bc_matrix.h5")
  
        fltrd_file_path <- paste0(cellranger_dir,
                                  "11389-DA-", i, "/count/sample_filtered_feature_bc_matrix.h5")
  
        cb_file_path <- paste0(cellbender_dir,
                               "output_cellbender_11389-DA-", i, "_FPR_0.0_filtered_seurat.h5")
    }
  }
  
  
  ## Read the counts ----
  cat('\nloading from',raw_file_path, 'for ',selected_timepoint,' sample:',i,'\n')
  raw_matrix <- Read10X_h5(filename=raw_file_path)

  cat('\nloading from',fltrd_file_path, 'for ',selected_timepoint,' sample:',i,'\n')
  fltrd_matrix <- Read10X_h5(filename=fltrd_file_path)

  cat('\nloading from',cb_file_path, 'for ',selected_timepoint,' sample:',i,'\n')
  cb_matrix <- Read10X_h5(filename=cb_file_path)
  
  ## run EmptyDrops on raw sample ----
  cat("---------------------------------------------\n")
  cat('\nRunning EmptyDrops for ',selected_timepoint,' sample:',i,'\n')
  cat("---------------------------------------------\n")

  e.out <- emptyDrops(raw_matrix,
                        lower=100, # lower bound of total UMIs. default.
                        test.ambient=TRUE,
                        niters=10000, # n iterations for monte carlo p-value calculation
                        retain=Inf,
                        BPPARAM = (workers = nthreads))
  fdr_cutoff = configs$emptydrops$default$fdr_cutoff
  is.cell = e.out$FDR <= configs$emptydrops$default$fdr_cutoff
  e.out$is.cell = is.cell # droplets with FDR <= fdr_cutoff considered non-empty. Can adjust later with df.
  cat("Number of cells detected by EmptyDrops at FDR <=0.01:", sum(is.cell, na.rm = TRUE), "\n")
  print(table(Limited=e.out$Limited, Significant=is.cell))
  cat("Number of cells detected by EmptyDrops at FDR <=0.001:", sum(e.out$FDR <= 0.001, na.rm = TRUE), "\n")

  e.calls <- data.frame(
    row.names = rownames(e.out),
    emptydrops_FDR = e.out$FDR,
    emptydrops_logprob = e.out$LogProb,
    emptydrops_pvalue = e.out$PValue,
    emptydrops_Limited = e.out$Limited,
    emptydrops_is_cell = is.cell
  )

  ## histogram of emptyDrops null p-values ----
  cat('\nPlotting EmptyDrops graphs for sample:',i,'\n')
  png(filename = paste0(fig_dir,'emptyDrops/',selected_timepoint,'_sample',i,'_histogram_emptyDrops_null_pvalues.png'), units='in', width=10, height=10, res=180)
  hist(e.out[which(e.out$Total < 101),]$PValue,breaks = 100,
        main = 'EmptyDrops null hypothesis p-values',
        sub = '(cells < 100 UMI counts)',
        xlab='p-value')
  dev.off()
  
  ## Create emptyDrops plot ----
  ggplot(e.out, aes(x = Total, y = -LogProb)) +
    # Add points with color based on cell classification
    geom_jitter(aes(color = is.cell), alpha = 0.6, size = 1.2) +
    
    # Use log10 scale for x-axis with comma formatting for readability
    scale_x_log10(labels = comma, n.breaks = 6) +
    scale_y_log10(labels = comma, n.breaks = 6) +
    annotation_logticks(color='black') +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("TRUE" = "Cell", "FALSE" = "Background"),
      name = "Classification"
    ) +
    labs(
      title = "Cell Classification by UMI Count and Probability",
      subtitle = "Cells with < 100 UMI counts considerd empty",
      x = "Total UMI count (log scale)",
      y = "-Log Probability"
    ) +
    
    # Add a vertical reference line at a threshold value
    geom_vline(xintercept = 100, 
                linetype = "dashed", color = "blue", alpha = 0.7) +
    theme_classic2() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  ggsave(paste0(fig_dir,'emptyDrops/',selected_timepoint,'_sample',i,'_emptyDrops_cell_calls.png'), units='in', width=10, height=10, dpi=320)
    
    
  ## Create Seurat objects ----
  cat('\ncreating seurat objects for ',selected_timepoint,' sample:',i,'\n')
  
  raw_seurat_obj <- CreateSeuratObject(
    counts  = raw_matrix,
    assay   = "RNA",
    project = paste0("sample", i)
  )
  
  # joint seurat object with cellbender and raw cell range counts for cellbender QC.
  joined_cb_and_raw <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_matrix,
                                            raw_counts_matrix = raw_matrix)
  joined_cb_and_raw <- AddMetaData(joined_cb_and_raw, metadata=e.calls) # add emptyDrops calls to metadata for later comparison

  fltrd_seurat_obj <- CreateSeuratObject(
    counts  = fltrd_matrix,
    assay   = "RNA",
    project = paste0("sample", i)
  )

  cb_seurat_obj <- CreateSeuratObject(
    counts  = cb_matrix,
    assay   = "RNA",
    project = paste0("sample", i)
  )
  cb_seurat_obj <- AddMetaData(cb_seurat_obj, metadata=e.calls) # add emptyDrops calls to metadata

  # read in cellbender metadata, i.e. cell probabilites, cell size, 
  # droplet efficiency (reverse transcription efficiency within the droplet), 
  # background fraction
  
  cat('\nMerging cellbedender metadata with seurat objects for sample:',i,'\n')
  cellbender_h5_output_path = ifelse(selected_timepoint == 'week3',
                                      paste0(cellbender_fulloutput_dir, "output_cellbender_8651-JC-",i,"_FPR_0.0.h5"),
                                      paste0(cellbender_fulloutput_dir, "output_cellbender_11389-DA-",i,"_FPR_0.0.h5")
  )
  ## grab cellbender metadata from hdf5 file ----
  metadata = h5read(cellbender_h5_output_path, 'metadata')
  droplet_latents = h5read(cellbender_h5_output_path, 'droplet_latents')
  droplet_latents = as.data.frame(droplet_latents[-c(2,6)]) # don't want gene_expression encoding or barcode indices
  rownames(droplet_latents) = metadata$barcodes_analyzed
  
  ## grab cellbender ambient profile
  matrices = h5read(cellbender_h5_output_path, 'matrix')
  global_latents = h5read(cellbender_h5_output_path, 'global_latents')
  global_latents = as.data.frame(global_latents)[1]
  global_latents$gene = matrices$features$name
  global_latents$ensemble_ids = matrices$features$id
  global_latents = global_latents[,c(2,1,3)]
  
  # save ambient profile
  write.csv(global_latents, paste0(data_dir,'cellbender_ambient_profiles/',selected_timepoint,'_sample',i,'_ambient_expression.csv'), row.names = FALSE)

  # add cellbender metadata to all seurat objects
  raw_seurat_obj = AddMetaData(raw_seurat_obj, metadata=droplet_latents)
  joined_cb_and_raw = AddMetaData(joined_cb_and_raw, metadata=droplet_latents)
  fltrd_seurat_obj = AddMetaData(fltrd_seurat_obj, metadata=droplet_latents)
  cb_seurat_obj = AddMetaData(cb_seurat_obj, metadata=droplet_latents)
})


# --- Merging ---
cat("\nMerging Seurat Objects...\n")
merged_cb <- merge(processed_data[[1]], processed_data[-1], add.cell.ids = paste0("tp_", samples_to_process)) %>% JoinLayers()

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