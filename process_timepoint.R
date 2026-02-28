#!/usr/bin/env Rscript
# Usage: Rscript 01_process_checkpoint.R week3

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Provide a timepoint (e.g., week3)", call. = FALSE)
selected_timepoint <- args[1]

cat("=================================================\n")
cat("Generating QC Checkpoint for:", selected_timepoint, "\n")
cat("=================================================\n")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(scDblFinder)
  library(DropletUtils)
  library(rhdf5)
  library(SoupX)
 library(qs2) # for fast data saving/loading
  library(plyr)
  library(BiocParallel)
})

set.seed(281330800)
options(future.globals.maxSize = 1024^3)
nthreads <- 4

# Directories (Update these to match your exact environment)
base_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/"
cellbender_dir <- paste0(base_dir, "cellbender/cellbender_h5_outputs_for_seurat/FPR_0.0/")
cellranger_dir <- paste0(base_dir, "10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/")
data_dir <- paste0(base_dir, "analysis/", selected_timepoint, "/data/checkpoints/")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Sample Mappings (Truncated for brevity)
timepoints <- list("week3" = 21:32, "week10" = 33:40) 
samples_to_process <- timepoints[[selected_timepoint]]

# Utility: Calculate Base Metrics
calc_base_metrics <- function(seurat_obj) {
  seurat_obj$complexity <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt-", col.name = "percent.mt")
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]", col.name = "percent.rb")
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[^(P|E|S)]", col.name = "percent.hb")
  return(seurat_obj)
}

# --- Core Processing Loop ---
processed_data <- lapply(samples_to_process, function(i) {
  cat("\nProcessing Sample:", i, "\n")
  
  # Load CellBender and Raw Matrices (Add your specific pathing logic here)
  raw_mat <- Read10X_h5(paste0(cellranger_dir, "11389-DA-", i, "/count/sample_raw_feature_bc_matrix.h5"))
  cb_mat <- Read10X_h5(paste0(cellbender_dir, "output_cellbender_11389-DA-", i, "_FPR_0.0_filtered_seurat.h5"))
  
  # 1. Run EmptyDrops (Tag, do not subset)
  e.out <- emptyDrops(raw_mat, lower=100, test.ambient=TRUE, niters=10000, BPPARAM = SnowParam(workers = nthreads))
  e_calls <- data.frame(
    row.names = rownames(e.out),
    emptyDrops_FDR = e.out$FDR,
    emptyDrops_is_cell = e.out$FDR <= 0.01
  )
  
  # 2. Create Object and Add Metadata
  cb_obj <- CreateSeuratObject(counts = cb_mat, assay = "RNA", project = paste0("sample", i))
  cb_obj <- AddMetaData(cb_obj, metadata = e_calls)
  cb_obj <- calc_base_metrics(cb_obj)
  
  return(cb_obj)
})
names(processed_data) <- paste0("sample", samples_to_process)

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