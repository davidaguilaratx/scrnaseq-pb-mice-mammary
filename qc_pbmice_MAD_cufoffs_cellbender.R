library(Seurat)
library(SeuratObject)
library(scCustomize)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(scDblFinder)
library(DropletUtils)
library(rhdf5)
library(SoupX)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(qs2) # for fast data saving/loading
library(BiocParallel)
# library(metap)
# library(RCurl)
# library(Polychrome)

# set number of threads to use for saving qs formatted data
nthreads <- 2 

# set seed
# set.seed(2024)
set.seed(281330800)

# create polychrome color palette ----

# p12 = createPalette(12, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
# swatch(p12)
# names(p12) = NULL

# read in data ----

# data directories

## default directory
default_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/"
default_dir_ext = "D:/PhD/scRNAseq/"
default_dir_turbo = "//sph-colacino-win.turbo.storage.umich.edu/sph-colacino/aguilada/scRNAseq/"

## cellbender data directories----
# for original h5 outputs
# cellbender_fulloutput_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs/FPR_0.0/"
cellbender_fulloutput_dir <- paste0(default_dir,"/cellbender/cellbender_h5_outputs/FPR_0.0/")

# for h5 outputs processed to be seurat friendly
cellbender_dir <- paste0(default_dir,"/cellbender/cellbender_h5_outputs_for_seurat/FPR_0.0/")
# cellbender_dir <- paste0(default_dir,"/cellbender/cellbender_h5_outputs/FPR_0.0/")
# cellbender_dir <- paste0(default_dir,"/cellbender/cellbender_h5_outputs_for_seurat/FPR_0.01/")

# cell ranger data directory
# cellranger_dir <- "//sph-colacino-win.turbo.storage.umich.edu/sph-colacino/aguilada/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/"


## analysis type, cellbender or cellranger ----
analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analysis' # we set emptydrops-minimum-umis = 100, and default values for all else

## directories for saving based on timepoint being analyzed ----
selected_timepoint <- 'week3'
# selected_timepoint <- 'week10'
# selected_timepoint <- 'month7'
# selected_timepoint <- 'month18'

# cellranger_dir <- ifelse(selected_timepoint == 'week3', 
#                          "D:/PhD/scRNAseq/8651-JC-mammary/10x_analysis_JC-8651_cellranger_v.9.0.0/",
#                          "D:/PhD/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/")

# Directories where we'll save our generated data and QC figures
timepoint_dirs <- list(
  "week3" = list(
    fig_dir = paste0(default_dir,"/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/qc/"),
    data_dir = paste0(default_dir,"/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/data/")
    ),
  "week10" = list(
    fig_dir = paste0(default_dir,"/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/qc/"),
    data_dir = paste0(default_dir,"/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/data/")  
    ),
  "month7" = list(
    fig_dir = paste0(default_dir,"/analysis/7_months/",analysis,"/FPR_0.0_MAD/qc/"),
    data_dir = paste0(default_dir,"/analysis/7_months/",analysis,"/FPR_0.0_MAD/data/")
    ),
  "month18" = list(
    fig_dir = paste0(default_dir,"/analysis/18_months/",analysis,"/FPR_0.0_MAD/qc/"),
    data_dir = paste0(default_dir,"/analysis/18_months/",analysis,"/FPR_0.0_MAD/data/")
  )
)

# timepoint_dirs <- list(
#   "week3" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/qc/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/"
#   ),
#   "week10" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/qc/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/data/"
#   ),
#   "month7" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/qc/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/data/"
#   ),
#   "month18" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/qc/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/data/"
#   )
# )

list2env(timepoint_dirs[[selected_timepoint]], envir = environment()) # save directory variables to environment

## Create directories to save results if they don't already exist: ----
# dirs = c(cellbender_dir, cellranger_dir, data_dir, fig_dir)
# 
# for (dir in dirs) {
#   if (!dir.exists(dir)) { dir.create(dir,
#                                      recursive=TRUE) }
# }

## initialize empty list to store seurat objects ----
# mammary.list <- list() # for cellbender data filtered by emptyDrops calls
# cb.list <- list() # for cellbender data. NO emptyDrops filtering.
# raw.list <- list() # for cellranger filtered data
# fltrd.list <- list() # for cellranger filtered data
# joined.list <- list() # for joined cellbender and raw cellranger data (to validate cellbender on same barcodes)
# soupx.list <- list() # for soupx filtered cell ranger matrices

## timepoint sample number mappings ----
timepoints <- list(
  "week3" = 21:32,
  "week10" = 33:40,
  "month7" = 1:8,
  "month18" = 25:32
)
# Create treatment and control group mapping
treatment_groups <- list(
  "week3" = list(
    "ctrl" = 21:26,
    "pb" = 27:32
  ),
  "week10" = list(
    "ctrl" = 33:36,
    "pb" = 37:40
  ),
  "month7" = list(
    "ctrl" = 1:4,
    "pb" = 5:8
  ),
  "month18" = list(
    "ctrl" = 25:28,
    "pb" = 29:32
  )
)

# directories sample data outputs
timepoint_pools = list(
  "week3" = NA,
  "week10" = '11389-DA_Pool03/outs/per_sample_outs/',
  "month7" = '11389-DA_Pool01/outs/per_sample_outs/',
  "month18" = '11389-DA_Pool02/outs/per_sample_outs/'
)

## create extra directories for emptyDrops and soupx data ----
# if (!dir.exists(paste0(data_dir,'ambient_expression_profiles/'))) { 
#   dir.create(paste0(data_dir,'ambient_expression_profiles/'),recursive=TRUE) }
# 
# if (!dir.exists(paste0(data_dir,'emptyDrops/'))) { 
#   dir.create(paste0(data_dir,'emptyDrops/'),recursive=TRUE) }
# 
# if (!dir.exists(paste0(fig_dir,'emptyDrops/'))) { 
#   dir.create(paste0(fig_dir,'emptyDrops/'),recursive=TRUE) }
# 
# if (!dir.exists(paste0(data_dir,'soupx/'))) { 
#   dir.create(paste0(data_dir,'soupx/'),recursive=TRUE) }
# 
# if (!dir.exists(paste0(data_dir,'cellbender_ambient_profiles/'))) { 
#   dir.create(paste0(data_dir,'cellbender_ambient_profiles/'),recursive=TRUE) }

options(future.globals.maxSize = 1024^3)  # Set global memory limits to 1 GB
for (selected_timepoint in names(timepoints)) {
  # create directories
  list2env(timepoint_dirs[[selected_timepoint]], envir = environment())  # save directory variables to environment

  # cellranger_dir <- ifelse(selected_timepoint == 'week3',
  #                          "D:/PhD/scRNAseq/8651-JC-mammary/10x_analysis_JC-8651_cellranger_v.9.0.0/",
  #                          "D:/PhD/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/")
  
  cellranger_dir <- ifelse(selected_timepoint == 'week3',
                           paste0(default_dir_ext,"/8651-JC-mammary/10x_analysis_JC-8651_cellranger_v.9.0.0/"),
                           paste0(default_dir_turbo,"/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/"))

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


  mammary.list <- list() # for cellbender data filtered by emptyDrops calls
  cb.list <- list() # for cellbender data. NO emptyDrops filtering.
  raw.list <- list() # for cellranger raw data
  fltrd.list <- list() # for cellranger filtered data
  joined.list <- list() # for joined cellbender and raw cellranger data (to validate cellbender on same barcodes)
  soupx.list <- list() # for soupx filtered cell ranger matrices

  # Loop over numbers in file directory names and read in data into seurat objects ----
  for (i in timepoints[[selected_timepoint]]) {
  
    # if (selected_timepoint == 'week3') {
    # # Construct the file path for each sample directory
    # raw_file_path <- paste0(cellranger_dir,
    #   "8651-JC-", i, "/count/sample_raw_feature_bc_matrix.h5")
    # 
    # fltrd_file_path <- paste0(cellranger_dir,
    #   "8651-JC-", i, "/count/sample_filtered_feature_bc_matrix.h5")
    # 
    # cb_file_path <- paste0(cellbender_dir,
      # "output_cellbender_8651-JC-", i, "_FPR_0.0_filtered_seurat.h5")

    if (selected_timepoint == 'week3') {
      # Construct the file path for each sample directory
      raw_file_path <- paste0(cellranger_dir,
                              "8651-JC-", i, "/count/sample_raw_feature_bc_matrix.h5")
      
      fltrd_file_path <- paste0(cellranger_dir,
                                "8651-JC-", i, "/count/sample_filtered_feature_bc_matrix.h5")
      
      cb_file_path <- paste0(cellbender_dir,
                             "output_cellbender_8651-JC-", i, "_FPR_0.0_filtered_seurat.h5")
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
        # "output_cellbender_8651-JC-", i, "_FPR_0.01_filtered_seurat.h5")
  
        # Read the counts
        # cat('\nloading from',raw_file_path, 'for sample:',i,'\n')
        # raw_matrix <- Read10X_h5(filename=raw_file_path)
        # cat('\nloading from',fltrd_file_path, 'for sample:',i,'\n')
        # fltrd_matrix <- Read10X_h5(filename=fltrd_file_path)
        # cat('\nloading from',cb_file_path, 'for sample:',i,'\n')
        # cb_matrix <- Read10X_h5(filename=cb_file_path)
        #
        # Create a Seurat objects
        # joined_cb_and_raw <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_matrix,
        # raw_counts_matrix = raw_matrix)
        #
        #
        # fltrd_seurat_obj <- CreateSeuratObject(
        #   counts  = fltrd_matrix,
        #   assay   = "RNA",
        #   project = paste0("sample", i)
        # )
        #
        #
        # cb_seurat_obj <- CreateSeuratObject(
        #   counts  = cb_matrix,
        #   assay   = "RNA",
        #   project = paste0("sample", i)
        # )
        #
        # # Create a joined seurat object, which will
        #
        # # Store in a named list
        # # Example: list element named "count21", "count22", etc.
        # mammary.list[[paste0("sample", i)]] <- cb_seurat_obj
        # fltrd.list[[paste0("sample", i)]] <- fltrd_seurat_obj
        # joined.list[[paste0('sample',i)]] <- joined_cb_and_raw
        #
        # cat('Finished loading sample',i)
      }
    }
  
    ## Read the counts ----
    cat('\nloading from',raw_file_path, 'for sample:',i,'\n')
    raw_matrix <- Read10X_h5(filename=raw_file_path)
    
    cat('\nloading from',fltrd_file_path, 'for sample:',i,'\n')
    fltrd_matrix <- Read10X_h5(filename=fltrd_file_path)
    
    cat('\nloading from',cb_file_path, 'for sample:',i,'\n')
    cb_matrix <- Read10X_h5(filename=cb_file_path)
    
    
    ## run EmptyDrops on raw sample ----
    cat('\nRunning EmptyDrops for sample:',i,'\n')
  
    e.out <- emptyDrops(raw_matrix,
                        lower=100, # lower bound of total UMIs. default.
                        test.ambient=TRUE,
                        niters=10000, # n iterations for monte carlo p-value calculation
                        retain=Inf,
                        BPPARAM = SnowParam(workers = 4))
    register(SerialParam()) # close extra nodes/workers used
    bpstop()
    gc() # clean up RAM
    
    is.cell = e.out$FDR <= 0.01
    e.out$is.cell = is.cell # droplets with FDR <= 0.01 considered non-empty.
    sum(e.out$is.cell, na.rm = TRUE)
    table(Limited=e.out$Limited, Significant=is.cell)
  
    qs_save(e.out, paste0(data_dir,'emptyDrops/',selected_timepoint,'_sample',i,'.qs'), nthreads=2) # save file
    
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
    
    
    ## Create a Seurat objects ----
    cat('\ncreating seurat objects for sample:',i,'\n')
    
    raw_seurat_obj <- CreateSeuratObject(
      counts  = raw_matrix,
      assay   = "RNA",
      project = paste0("sample", i)
    )
    
    # joint seurat object with cellbender and raw cell range counts for cellbender QC.
    joined_cb_and_raw <- Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_matrix,
                                             raw_counts_matrix = raw_matrix)
  
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
    
    write.csv(global_latents, paste0(data_dir,'cellbender_ambient_profiles/',selected_timepoint,'_sample',i,'_ambient_expression.csv'), row.names = FALSE)
  
    raw_seurat_obj = AddMetaData(raw_seurat_obj, metadata=droplet_latents)
    joined_cb_and_raw = AddMetaData(joined_cb_and_raw, metadata=droplet_latents)
    fltrd_seurat_obj = AddMetaData(fltrd_seurat_obj, metadata=droplet_latents)
    cb_seurat_obj = AddMetaData(cb_seurat_obj, metadata=droplet_latents)
    
    
    ## subset cellbender seurat object by EmptyDrop cell calls ----
    cb_ed_seurat_obj = subset(cb_seurat_obj, cells=WhichCells(cb_seurat_obj, cells=rownames(e.out[which(e.out$is.cell==T),])))
  
    ## run soupx on filtered cellranger data ----
    # need same number of genes in raw and filtered matrices
    raw_matrix_subset <- raw_matrix[rownames(raw_matrix) %in% rownames(fltrd_matrix), ]
    soup.channel = SoupChannel(raw_matrix_subset, fltrd_matrix)
    
    ### normalize and cluster filtered cell ranger data to feed clusters into soupx ----
    cat('\nClustering filtered cell ranger data for sample:',i,'\n')
    fltrd_seurat_obj <- SCTransform(fltrd_seurat_obj, verbose = T)
    fltrd_seurat_obj <- RunPCA(fltrd_seurat_obj, verbose = T)
    fltrd_seurat_obj <- FindNeighbors(fltrd_seurat_obj, dims = 1:50, verbose = T)
    fltrd_seurat_obj <- FindClusters(fltrd_seurat_obj, verbose = T)
    fltrd_seurat_obj <- RunUMAP(fltrd_seurat_obj, dims = 1:50, verbose = T)
    
    ### estimate soup ----
    cat('\nRunning SoupX for sample:',i,'\n')
    soup.channel = setClusters(soup.channel,
                               setNames(fltrd_seurat_obj@meta.data$seurat_clusters,
                                        rownames(fltrd_seurat_obj@meta.data)))
    soup.channel = setDR(soup.channel, fltrd_seurat_obj@reductions$umap@cell.embeddings)
    soup.channel = autoEstCont(soup.channel)
    
    qs_save(soup.channel, paste0(data_dir,'soupx/',selected_timepoint,'_sample',i,'_soup.channel.qs'), nthreads=2) # save file
    
    soup = soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ]
    head(soup, n = 20)
    
    soup$gene = rownames(soup)
    soup = soup[,c(3,1,2)]
    write.csv(soup, paste0(data_dir,'soupx/ambient_soupx_profile_',selected_timepoint,'_sample',i,'.csv'), row.names = F)
    
    adj.matrix = adjustCounts(soup.channel, roundToInt = T)
    
    soupx_seurat_obj = CreateSeuratObject(
      counts  = adj.matrix,
      assay   = "RNA",
      project = paste0("sample", i)
    )
    
    # add cellbender droplet metadata for comparison
    soupx_seurat_obj = AddMetaData(soupx_seurat_obj, metadata=droplet_latents)
    
    # remove clustering data from fltrd_seurat_obj.
    DefaultAssay(fltrd_seurat_obj) = 'RNA'
    fltrd_seurat_obj = DietSeurat(fltrd_seurat_obj,assays='RNA')
    
    # Store in a named list
    # Example: list element named "count21", "count22", etc.
    mammary.list[[paste0("sample", i)]] <- cb_ed_seurat_obj
    cb.list[[paste0("sample", i)]] <- cb_seurat_obj
    raw.list[[paste0("sample", i)]] <- raw_seurat_obj
    fltrd.list[[paste0("sample", i)]] <- fltrd_seurat_obj
    joined.list[[paste0('sample',i)]] <- joined_cb_and_raw
    soupx.list[[paste0('sample',i)]] <- soupx_seurat_obj
  
    cat('Finished loading sample',i)
  }
  
  rm(cb_matrix,raw_matrix,fltrd_matrix,adj.matrix,joined_cb_and_raw,cb_seurat_obj, # clean up RAM
     cb_ed_seurat_obj, raw_seurat_obj, fltrd_seurat_obj, soupx_seurat_obj,
     soup.channel,soup, raw_matrix_subset, droplet_latents, matrices, global_latents,
     metadata, e.out) 
  gc() # free RAM
  
  
  # calculate log10 z-scores; percent mitochondrial, ribosomal, and hemoglobin gene expression; and complexity ----
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
  
  mammary.list = calc_metrics(mammary.list)
  cb.list = calc_metrics(cb.list)
  raw.list = calc_metrics(raw.list)
  fltrd.list = calc_metrics(fltrd.list)
  joined.list = calc_metrics(joined.list)
  soupx.list = calc_metrics(soupx.list)
  
  
  ## create list of sample names, separating count data into ctrl and pb and add timepoint prefix ----
  # cell_ids = unlist(lapply(names(mammary.list), function(name) {
  #   n = as.numeric(gsub('[^0-9]', '', name))
  #   n = ifelse(n < 27, paste0('ctrl',n), paste0('pb',n))
  # }))
  
  timepoint_sample_prefix = list(
    "week3" = 'wk3',
    "week10" = 'wk10',
    "month7" = 'mn7',
    "month18" = 'mn18'
  )
  
  prefix = timepoint_sample_prefix[[selected_timepoint]]
  cell_ids = unlist(lapply(names(mammary.list), function(name) {
    n = as.numeric(gsub('[^0-9]', '', name)) # sample number
    paste0(ifelse(n %in% treatment_groups[[selected_timepoint]]$ctrl, paste0(prefix,"_ctrl",n), paste0(prefix,"_pb",n)))
  }))
  cell_ids
  
  # merge seurat object lists into one object ----
  cat('\nMerging seurat objects for',selected_timepoint,'\n')
  merged_seurat_cb = merge(mammary.list[[1]], mammary.list[-1], add.cell.ids=cell_ids)
  merged_seurat_cb_orig = merge(cb.list[[1]], cb.list[-1], add.cell.ids=cell_ids)
  merged_seurat_raw = merge(raw.list[[1]], raw.list[-1], add.cell.ids=cell_ids)
  merged_seurat_fltrd = merge(fltrd.list[[1]], fltrd.list[-1], add.cell.ids=cell_ids)
  merged_seurat_joined = merge(joined.list[[1]], joined.list[-1], add.cell.ids=cell_ids)
  merged_seurat_soupx = merge(soupx.list[[1]], soupx.list[-1], add.cell.ids=cell_ids)
  
  rm(mammary.list, raw.list, fltrd.list, joined.list, cb.list, soupx.list) # clean up RAM
  gc()
  
  # add more metadata, such as sample name, cells, and condition ----
  add_metadata = function(seurat_object) {
    # create metadata dataframe
    metadata = seurat_object@meta.data
    
    # add a cell IDs column to metadata
    metadata$cells = rownames(metadata)
    
    # add sample type column to metadata
    metadata$sample = str_extract(metadata$cells, '^[^_]+_[^_]+') # ^start of text. [^_]+ match any charcters that are not underscore.
    
    
    # add timepoint column to metadata
    metadata$timepoint = str_extract(metadata$cells, '^[^_]+')
    metadata$timepoint = factor(metadata$timepoint,
                        levels=c('wk3','wk10','mn7','mn18'),
                        ordered = TRUE)
    
    # week 3 samples run on different date than rest of timepoints. 
    # That being said, each timepint was run on a different flow cell. 
    # Will add a batch column for the run date
    metadata$batch = ifelse(metadata$timepoint =='wk3','Experiment_1','Experiment_2')
    
    # add condition and exposed columns to metadata
    metadata$condition = str_extract(metadata$cells, '(?<=_)[^0-9_]+') #?<= checks for underscore. capture non-number or underscore text after
    metadata$exposed = ifelse(metadata$condition == 'pb', 1, 0)
    
    
    # save metadata dataframe back to merged seurat data
    seurat_object@meta.data = metadata
    
    return(seurat_object)
  }
  
  # add metadata 
  cat('\nadding metadata to seurat objects for',selected_timepoint,'\n')
  merged_seurat_cb = add_metadata(merged_seurat_cb)
  merged_seurat_cb_orig = add_metadata(merged_seurat_cb_orig)
  merged_seurat_raw = add_metadata(merged_seurat_raw)
  merged_seurat_fltrd = add_metadata(merged_seurat_fltrd)
  merged_seurat_joined = add_metadata(merged_seurat_joined)
  merged_seurat_soupx = add_metadata(merged_seurat_soupx)
  
  # add cell call column to label whether cell was called by cellbender and cellranger
  # or cellbender only
  cellbender_cells = colnames(merged_seurat_cb)
  cellranger_cells = colnames(merged_seurat_fltrd)
  all_cells = union(cellbender_cells, cellranger_cells)
  cell_calls = ifelse(all_cells %in% cellbender_cells & all_cells %in% cellranger_cells,
                      "cellbender and cellranger",
                      ifelse(all_cells %in% cellranger_cells,
                             "cellranger only",
                             "cellbender only"))
  
  table(cell_calls)
  t = intersect(cellbender_cells, cellranger_cells)
  # cellbender and cellranger           cellbender only           cellranger only 
  # 46080                               24125                        12 
  
  # cellbender and cellranger           cellbender only           cellranger only #v.9.0.0
  # 49620                     20585                      1227 
  
  merged_seurat_cb$cell_calls = cell_calls
  table(merged_seurat_cb$cell_calls)
  # cellbender and cellranger           cellbender only 
  # 46080                     24125 
  
  # cellbender and cellranger           cellbender only # v.9.0.0
  # 49620                     20585 
  
  # join seurat layers since we can group by sample in metadata ----
  cat('\njoining layers in seurat objects for',selected_timepoint,'\n')
  merged_seurat_cb = JoinLayers(merged_seurat_cb)
  merged_seurat_cb_orig = JoinLayers(merged_seurat_cb_orig)
  # merged_seurat_raw = JoinLayers(merged_seurat_raw)
  merged_seurat_fltrd = JoinLayers(merged_seurat_fltrd)
  merged_seurat_joined[['RNA']] = JoinLayers(merged_seurat_joined[['RNA']]) # RAW assay already has only one layer
  merged_seurat_soupx = JoinLayers(merged_seurat_soupx)
  gc() # free RAM
  
  # # save merged seurat objects and metadata
  qs_save(merged_seurat_cb, paste0(data_dir,'/merged_seurat_cb.qs'),nthreads=nthreads)
  qs_save(merged_seurat_cb_orig, paste0(data_dir,'/merged_seurat_cb_orig.qs'),nthreads=nthreads)
  qs_save(merged_seurat_raw, paste0(data_dir,'/merged_seurat_raw.qs'),nthreads=nthreads)
  qs_save(merged_seurat_fltrd, paste0(data_dir,'/merged_seurat_fltrd.qs'),nthreads=nthreads)
  qs_save(merged_seurat_joined, paste0(data_dir,'/merged_seurat_joined.qs'),nthreads=nthreads)
  qs_save(merged_seurat_soupx, paste0(data_dir,'/merged_seurat_soupx.qs'),nthreads=nthreads)
  
  
  qs_save(merged_seurat_cb@meta.data, paste0(data_dir,'/metadata_cb.qs'),nthreads=nthreads)
  qs_save(merged_seurat_cb_orig@meta.data, paste0(data_dir,'/metadata_cb_orig.qs'),nthreads=nthreads)
  qs_save(merged_seurat_raw@meta.data, paste0(data_dir,'/metadata_raw.qs'),nthreads=nthreads)
  qs_save(merged_seurat_fltrd@meta.data, paste0(data_dir,'/metadata_fltrd.qs'),nthreads=nthreads)
  qs_save(merged_seurat_joined@meta.data, paste0(data_dir,'/metadata_joined.qs'),nthreads=nthreads)
  qs_save(merged_seurat_soupx@meta.data, paste0(data_dir,'/metadata_soupx.qs'),nthreads=nthreads)
  
  rm(merged_seurat_raw) # free RAM
  gc()
  
  # merged_seurat_cb = qs_read(paste0(data_dir, 'merged_seurat_cb.qs'),nthreads=nthreads)
  # merged_seurat_fltrd = qs_read(paste0(data_dir, 'merged_seurat_fltrd.qs'),nthreads=nthreads)
  # merged_seurat_joined = qs_read(paste0(data_dir, 'merged_seurat_joined.qs'),nthreads=nthreads)
  
  
  # number of cells per dataset ----
  table(Idents(merged_seurat_cb))
  # sample21 sample22 sample23 sample24 sample25 sample26 sample27 sample28 sample29 sample30 sample31 sample32 
  # 3965     5140     6121     3249     6475     7187     8335     5946     6505     8230     4551     4501 
  
  table(Idents(merged_seurat_fltrd))
  # sample21 sample22 sample23 sample24 sample25 sample26 sample27 sample28 sample29 sample30 sample31 sample32 
  # 2142     3035     4567     1720     6151     5671     7976     4396     1249     6560     2096      529 
  
  # sample21 sample22 sample23 sample24 sample25 sample26 sample27 sample28 sample29 sample30 sample31 sample32 # v.9.0.0
  # 2633     3345     4822     2078     6414     6157     8206     4597     2056     7059     2362     1118 
  
  Idents(merged_seurat_joined) = merged_seurat_joined$sample
  table(Idents(merged_seurat_joined))
  # ctrl21 ctrl22 ctrl23 ctrl24 ctrl25 ctrl26   pb27   pb28   pb29   pb30   pb31   pb32 
  # 3015   3630   4969   2551   4260   5636   6189   3694   4780   5834   3691   3008 
  
  
  
  
  # proposed initial QC cutoffs ----
  percent.mt.cutoff = 8
  gene_cutoff_low = 150
  gene_cutoff_high = 5000
  umi_cutoff_low = 150
  umi_cutoff_high = 10000
  complexity_cutoff = 0.8 # log10genes / log10UMIs
  min_cells_per_gene = 3 # filtered genes must be expressed in at least this # of cells
  
  # QC plots ----
  plot_qcs = function(seurat_object, percent.mt.cutoff=8, gene_cutoff_low=150, gene_cutoff_high=5000,
                      umi_cutoff_low=150, umi_cutoff_high=5000, complexity_cutoff=0.8, min_cells_per_gene=3,
                      fig_dir=fig_dir) {
    
    ## plot number of cells per dataset ----
    cat('\nPlotting number of cells per sample')
    ggplot(seurat_object@meta.data, aes(x=sample, fill=sample)) + 
      geom_bar() +
      theme_classic2() +
      NoLegend()+
      labs(title="Number of cells per experiment",
           subtitle = paste0('Total cells = ',dim(seurat_object@meta.data)[1],', # ctrl cells = ',
                             sum(seurat_object@meta.data$condition =='ctrl'),', # pb cells = ',
                             sum(seurat_object@meta.data$condition == 'pb')))+
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      scale_y_continuous(n.breaks=8)
    ggsave(paste0(fig_dir,'/ncells.png'), dpi=320, width=10, height= 10)
    
    ## plot proportion of mitochondrial genes per sample ----
    cat('\nplotting percent.mt per sample')
    VlnPlot(seurat_object, features = "percent.mt", group.by='sample', pt.size = 0.05)+
      NoLegend()+
      labs(title='Percent mitochondrial genes',
           subtitle=paste0('cutoff above ', percent.mt.cutoff,'%'),
           x='sample',
           y='percent.mt')+
      scale_y_continuous(n.breaks=20, limits=c(0,NA))+
      theme(axis.text.x = element_text(angle = 0, hjust=0.5),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      geom_hline(yintercept=percent.mt.cutoff, linetype='dashed',color='red')+
      annotate('text', x=7.5, y=10.5, label='cutoff', size=5, color='red')
    ggsave(paste0(fig_dir,'/percent_mt.png'), dpi=320, width=10, height= 10)
    
    ## plot novelty score (complexity = log10genes / log10UMI) ----
    cat('\nPlotting complexity per sample')
    ggplot(seurat_object@meta.data, aes(color=sample, x=complexity, fill=sample)) + 
      geom_density(alpha = 0.1)+
      theme_classic2()+
      labs(title='Overal complexity via genes detected per UMI (novelty score)',
           x='log10(nGenes) / log10(nUMI)',
           y='Cell density')+
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      geom_vline(xintercept = complexity_cutoff, linetype='dashed', color='red')+
      annotate('text', x=complexity_cutoff-0.02, y=40, label='cutoff', size=5, color='red')
    ggsave(paste0(fig_dir,'/complexity_novelty_scores.png'), dpi=320, width=10, height= 10)
    
    ## plot number of reads per cell as vln plots ----
    cat('\nPlotting UMIs per cell per sample as vln plots')
    VlnPlot(seurat_object, features = "nCount_RNA", group.by='sample', pt.size = 0.05)+
      theme_classic2()+
      NoLegend()+
      labs(title='Number of transcript UMIs per cell',
           subtitle=paste0('cutoffs below ',umi_cutoff_low,' and above ',umi_cutoff_high, ' UMIs'),
           x='sample',
           y='nUMI')+
      scale_y_log10(n.breaks=10,labels = comma_format())+
      annotation_logticks(sides = "l")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      geom_hline(yintercept=umi_cutoff_low, linetype='dashed',color='red')+
      geom_hline(yintercept=umi_cutoff_high, linetype='dashed',color='red')
    ggsave(paste0(fig_dir,'/umi_counts.png'), dpi=320, width=10, height= 10)
    
    # plot number of reads per cell as density plots ----
    cat('\nPlotting UMIs per cell per sample as density plots')
    ggplot(seurat_object@meta.data, aes(color=sample, x=nCount_RNA, fill=sample))+ 
      geom_density(alpha = 0.1)+
      labs(title='Number of transcript UMIs per cell',
           subtitle=paste0('cutoffs below ',umi_cutoff_low,' and above ',umi_cutoff_high, ' UMIs'),
           x='nUMI',
           y='Cell density')+
      theme_classic2()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      scale_x_log10(n.breaks=10, labels=comma_format())+
      annotation_logticks(sides = "b")+
      geom_vline(xintercept = umi_cutoff_low, linetype='dashed', color='red')+
      geom_vline(xintercept = umi_cutoff_high, linetype='dashed', color='red')
    ggsave(paste0(fig_dir,'/umi_counts_density.png'), dpi=320, width=10, height= 10)
    
    ## plot number of genes per cell as vln plots ----
    cat('\nPlotting genes per cell per sample as vln plots')
    VlnPlot(seurat_object, features = "nFeature_RNA", group.by='sample', pt.size = 0.05)+
      theme_classic2()+
      NoLegend()+
      labs(title='Number of genes detected per cell',
           subtitle=paste0('cutoffs below ',gene_cutoff_low,' and above ',gene_cutoff_high, ' genes'),
           x='sample',
           y='nGenes')+
      scale_y_log10(n.breaks=10, labels = comma_format())+
      annotation_logticks(sides = "l")+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')
    ggsave(paste0(fig_dir,'/gene_counts.png'), dpi=320, width=10, height= 10)
    
    ## plot number of genes per cell as density plots ----
    cat('\nPlotting genes per cell per sample as density plots')
    ggplot(seurat_object@meta.data, aes(color=sample, x=nFeature_RNA, fill=sample))+
      geom_density(alpha = 0.1)+
      labs(title='Number of genes detected per cell',
           subtitle=paste0('cutoffs below ',gene_cutoff_low,' and above ',gene_cutoff_high, ' genes'),
           x='nGene',
           y='Cell density')+
      theme_classic2()+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      scale_x_log10(n.breaks=10, labels=comma_format())+
      annotation_logticks(sides = "b")+
      geom_vline(xintercept = gene_cutoff_low, linetype='dashed', color='red')+
      geom_vline(xintercept = gene_cutoff_high, linetype='dashed', color='red')
    ggsave(paste0(fig_dir,'/gene_counts_density.png'), dpi=320, width=10, height= 10)
    
    ## plot count depth (number of UMIs) vs number of genes, colored by sample name (left) and percent mitochondrial gene counts (right) ----
    p1 = FeatureScatter(seurat_object, feature1='nCount_RNA', feature2='nFeature_RNA', group.by='sample', pt.size=0.5, shuffle=TRUE)+
      DarkTheme()+
      labs(title='# of UMIs vs # of genes',
           subtitle='colored by sample',
           x='# of UMI counts (log scale)',
           y='# of genes (log scale)',
           color='sample')+
      scale_x_log10(n.breaks=8, labels=comma_format())+
      scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
      annotation_logticks(sides = "bl", colour='white', size=1)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18),
            plot.background=element_rect(color='black'))+
      NoGrid()+
      stat_smooth(method=lm, linewidth=0.5, colour='#61e4e8')+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='red')
    
    p2 = ggplot(seurat_object@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
      geom_point(size=0.5)+
      DarkTheme()+
      labs(title='# of UMIs vs # of genes',
           subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
           x='# of UMI counts (log scale)',
           y='# of genes (log scale)',
           color='% mt counts')+
      scale_x_log10(n.breaks=8, labels=comma_format())+
      scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
      annotation_logticks(sides = "bl", colour='white', size=1)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18),
            plot.background=element_rect(color='black'))+
      NoGrid()+
      stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='red')+
      scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#292a2d')
    
    cat('\nPlotting UMIs vs genes colored by sample (left) and mt.percent (right)')
    plot_grid(p1, p2)
    ggsave(paste0(fig_dir,'/ncounts_vs_ngenes_coloredby_percentmt.png'), dpi=320, width=20, height= 10)
    
    ## plot count depth (number of reads) vs number of genes, colored by exposure condition and mitochondrial gene percentage ----
    p1 = FeatureScatter(seurat_object, feature1='nCount_RNA', feature2='nFeature_RNA', group.by='condition', pt.size=0.5, shuffle=TRUE)+
      DarkTheme()+
      labs(title='# of UMIs vs # of genes',
           subtitle='colored by condition',
           x='# of UMI counts (log scale)',
           y='# of genes (log scale)',
           color='condition')+
      scale_x_log10(n.breaks=8, labels=comma_format())+
      scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
      annotation_logticks(sides = "bl", colour='white', size=1)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18),
            plot.background=element_rect(color='black'))+
      NoGrid()+
      stat_smooth(method=lm, linewidth=0.5, colour='#61e4e8')+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='red')
    
    p2 = ggplot(seurat_object@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
      geom_point(size=0.5)+
      DarkTheme()+
      labs(title='# of UMIs vs # of genes',
           subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
           x='# of UMI counts (log scale)',
           y='# of genes (log scale)',
           color='% mt counts')+
      scale_x_log10(n.breaks=8, labels=comma_format())+
      scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
      annotation_logticks(sides = "bl", colour='white', size=1)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18),
            plot.background=element_rect(color='black'))+
      NoGrid()+
      stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='red')+
      scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#292a2d')
    
    cat('\nPlotting UMIs vs genes colored by condition (left) and mt.percent (right)')
    plot_grid(p1, p2)
    ggsave(paste0(fig_dir,'/ncounts_vs_ngenes_coloredby_percentmt_and_exposure_condition.png'), dpi=320, width=20, height= 10)
    
    ## more vln plots for UMIs, genes, and percent.mt by exposure condition ----
    cat('\nPlotting genes per condition')
    QC_Plots_Genes(seurat_object, group.by = 'condition', low_cutoff = gene_cutoff_low, high_cutoff = gene_cutoff_high)
    ggsave(paste0(fig_dir,'/genes_per_cell_by_exposure_condition.png'), dpi=320, width=10, height= 10)
    
    cat('\nPlotting UMIs per condition')
    QC_Plots_UMIs(seurat_object, group.by = 'condition', low_cutoff = umi_cutoff_low, high_cutoff = umi_cutoff_high, y_axis_log = TRUE)+
      scale_y_log10(n.breaks=10, labels=comma_format())+
      annotation_logticks(sides = "l", size=.5)
    ggsave(paste0(fig_dir,'/UMIs_per_cell_by_exposure_condition.png'), dpi=320, width=10, height= 10)
    
    cat('\nPlotting percent.mt per condition')
    QC_Plots_Mito(seurat_object, mito_name='percent.mt',group.by = 'condition', low_cutoff = percent.mt.cutoff)
    ggsave(paste0(fig_dir,'/percent.mt_per_cell_by_exposure_condition.png'), dpi=320, width=10, height= 10)
    
  }
  
  cat('\nplotting qcs for',selected_timepoint,'\n')
  plot_qcs(merged_seurat_cb)
  
  
  # What would median absolute deviation cutoffs look like? ----
  #' Calculate MAD thresholds for a metric in a Seurat object
  #' @param seurat_obj Seurat object
  #' @param metric_col Name of metric column in meta.data
  #' @param sample_col Name of sample column in meta.data (default: "sample")
  #' @return Data frame with MAD thresholds for each sample
  #' 
  calculate_mad_thresholds <- function(seurat_obj, 
                                       metric_col, 
                                       sample_col = "sample") {
    # Input validation
    if (!metric_col %in% names(seurat_obj@meta.data)) {
      stop(sprintf("Metric column '%s' not found in Seurat object", metric_col))
    }
    if (!sample_col %in% names(seurat_obj@meta.data)) {
      stop(sprintf("Sample column '%s' not found in Seurat object", sample_col))
    }
    
    # Get unique samples
    samples <- unique(seurat_obj$sample)
    
    # Calculate thresholds for each sample
    threshold_data <- lapply(samples, function(samp) {
      sample_metrics <- seurat_obj@meta.data[seurat_obj@meta.data[[sample_col]] == samp, metric_col]
      log1p_sample_metrics = log1p(sample_metrics)
      
      
      med_val <- median(sample_metrics,na.rm=T)
      mad_val <- mad(sample_metrics,na.rm=T)
      # med_val <- median(log1p_sample_metrics)
      # mad_val <- mad(log1p_sample_metrics)
      
      data.frame(
        sample = samp,
        median = med_val,
        mad = mad_val,
        mad3_lower = med_val - 3 * mad_val,
        mad3_upper = med_val + 3 * mad_val,
        mad5_lower = med_val - 5 * mad_val,
        mad5_upper = med_val + 5 * mad_val
      )
    })
    
    do.call(rbind, threshold_data)
  }
  
  #' Create violin plot with MAD threshold lines and hard lower cutoffs
  #' @param seurat_obj Seurat object
  #' @param metric_col Name of metric to plot
  #' @param mad_thresholds Data frame with MAD thresholds from calculate_mad_thresholds()
  #' @param title Plot title (optional)
  #' @param ylabel Y-axis label (optional)
  #' @param use_log_scale Use log10 scale for y-axis (default: TRUE)
  #' @param mad3_color Color for 3 MAD lines (default: "cyan")
  #' @param mad5_color Color for 5 MAD lines (default: "red")
  #' @param cutoff_color Color for hard cutoff lines (default: "magenta")
  #' @param hard_cutoffs Named list of hard cutoffs to display (optional)
  #' @return ggplot object
  plot_metric_with_mad <- function(seurat_obj,
                                   metric_col,
                                   mad_thresholds,
                                   gene_cutoff_low=150,
                                   umi_cutoff_low=150,
                                   title = NULL,
                                   ylabel = NULL,
                                   use_log_scale = TRUE,
                                   mad3_color = "cyan",
                                   mad5_color = "magenta",
                                   cutoff_color = "red",
                                   hard_cutoffs = list(
                                     nCount_RNA = 150,
                                     nFeature_RNA = 150,
                                     percent.mt = 8
                                   )) {
    
    # Create base violin plot
    p <- VlnPlot(seurat_obj,
                 features = metric_col,
                 group.by = 'sample',
                 pt.size = 0.05) +
      NoLegend() +
      labs(title = title %||% sprintf("Distribution of %s", metric_col),
           subtitle = "Per-sample MAD cutoffs shown",
           x = "Sample",
           y = ylabel %||% metric_col)
    
    # Add log scale if requested
    if (use_log_scale) {
      p <- p + scale_y_log10(n.breaks = 10, labels = comma_format())+
        annotation_logticks(sides = "l")
    } else {
      p <- p + scale_y_continuous(n.breaks=20, limits=c(0,NA))
    }
    
    # Prepare threshold data
    samples <- unique(seurat_obj$sample)
    threshold_data <- data.frame(
      sample = rep(samples, 4),
      y = c(mad_thresholds$mad3_lower,
            mad_thresholds$mad3_upper,
            mad_thresholds$mad5_lower,
            mad_thresholds$mad5_upper),
      type = rep(c("mad3_lower", "mad3_upper", "mad5_lower", "mad5_upper"),
                 each = length(samples))
    )
    
    # Add threshold lines
    p <- p +
      # 3 MAD thresholds
      geom_segment(data = threshold_data[threshold_data$type %in% 
                                           c("mad3_lower", "mad3_upper"), ],
                   mapping = aes(x = as.numeric(factor(sample)) - 0.4,
                                 xend = as.numeric(factor(sample)) + 0.4,
                                 y = y,
                                 yend = y),
                   inherit.aes = FALSE,
                   color = mad3_color,
                   linewidth = 1.5,
                   linetype = "dashed") +
      # 5 MAD thresholds
      geom_segment(data = threshold_data[threshold_data$type %in% 
                                           c("mad5_lower", "mad5_upper"), ],
                   mapping = aes(x = as.numeric(factor(sample)) - 0.4,
                                 xend = as.numeric(factor(sample)) + 0.4,
                                 y = y,
                                 yend = y),
                   inherit.aes = FALSE,
                   color = mad5_color,
                   linewidth = 1.5,
                   linetype = "dashed") 
    
    
    # Add hard cutoff line ONLY for non-log1p metrics
    if (!grepl("log1p_", metric_col) && metric_col %in% names(hard_cutoffs)) {
      cutoff_value <- hard_cutoffs[[metric_col]]
      
      # For percent.mt, the cutoff is an upper bound
      if (metric_col == "percent.mt") {
        p <- p + geom_hline(yintercept = cutoff_value, 
                            color = cutoff_color, 
                            linetype = "dashed",
                            linewidth = 1.0)
      } 
      # For other metrics, the cutoff is a lower bound
      else {
        p <- p + geom_hline(yintercept = cutoff_value, 
                            color = cutoff_color, 
                            linetype = "dashed",
                            linewidth = 1.0)
      }
    }
    
    # Add legend text explaining the lines
    legend_text <- sprintf("%s lines: ±3 MAD\n%s lines: ±5 MAD", 
                           mad3_color, mad5_color)
    
    # Add cutoff information to legend if applicable (only for non-log1p metrics)
    if (!grepl("log1p_", metric_col) && metric_col %in% names(hard_cutoffs)) {
      cutoff_value <- hard_cutoffs[[metric_col]]
      
      if (metric_col == "percent.mt") {
        legend_text <- paste0(legend_text, 
                              sprintf("\n%s line: upper hard cutoff at %g%%", 
                                      cutoff_color, cutoff_value))
      } else {
        legend_text <- paste0(legend_text, 
                              sprintf("\n%s line: lower hard cutoff at %g", 
                                      cutoff_color, cutoff_value))
      }
    }
    
    p <- p + annotate("text", x = -Inf, y = Inf,
                      label = legend_text,
                      hjust = -0.1, vjust = 1.2,
                      size = 3)
    
    return(p)
  }
  
  #' Analyze metric distributions and MAD thresholds
  #' @param seurat_obj Seurat object
  #' @param metrics Vector of metric names to analyze
  #' @param output_dir Directory for saving plots (optional)
  #' @param prefix Prefix for output files (optional)
  #' @param hard_cutoffs Named list of hard cutoffs to display (optional)
  #' @return List threshold data frames for each metric
  analyze_mad_metrics <- function(seurat_obj, 
                                  metrics, 
                                  output_dir = NULL,
                                  prefix = "MAD_",
                                  use_log_scale = TRUE,
                                  hard_cutoffs = list(
                                    nCount_RNA = 150,
                                    nFeature_RNA = 150,
                                    percent.mt = 8,
                                    complexity = 0.8
                                  )) {
    
    # Initialize results list
    results <- list()
    
    for (metric in metrics) {
      # Calculate MAD thresholds
      thresholds <- calculate_mad_thresholds(seurat_obj, metric)
      
      # Create plot
      p <- plot_metric_with_mad(
        seurat_obj = seurat_obj,
        metric_col = metric,
        mad_thresholds = thresholds,
        use_log_scale = !grepl("percent", metric),  # Don't log scale percentages
        hard_cutoffs = hard_cutoffs
      )
      
      # Save plot if output directory is specified
      if (!is.null(output_dir)) {
        plot_filename <- paste0(output_dir, 
                                sprintf("/%s%s", prefix, metric),'_cutoffs_per_sample.png')
        ggsave(plot_filename, plot = p, dpi = 320, width = 10, height = 10)
      }
      
      # Store results
      results[[metric]] <- list(
        thresholds = thresholds
      )
    }
    
    return(results)
  }
  
  # Example usage:
  # metrics_to_analyze <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  
  norm_log1p = function(seurat_obj) {
    seurat_obj$log1p_nCount_RNA = log1p(seurat_obj$nCount_RNA)
    seurat_obj$log1p_nFeature_RNA = log1p(seurat_obj$nFeature_RNA)
    seurat_obj$log1p_percent.mt = log1p(seurat_obj$percent.mt)
    
    return(seurat_obj)
  }
  
  merged_seurat_cb = norm_log1p(merged_seurat_cb)
  
  metrics_to_analyze <- c("nCount_RNA", "nFeature_RNA", "percent.mt", 
                          "log1p_nCount_RNA", "log1p_nFeature_RNA", "log1p_percent.mt")
  
  # Run analysis
  cat('\nanalyzing MAD metrics for',selected_timepoint,'\n')
  MAD_results <- analyze_mad_metrics(
    seurat_obj = merged_seurat_cb,
    metrics = metrics_to_analyze,
    output_dir = fig_dir,
    hard_cutoffs = list(
      nCount_RNA = 150,
      nFeature_RNA = 150,
      percent.mt = 8
    )
  )
  # 
  # MAD_results <- analyze_mad_metrics(
  #   seurat_obj = merged_seurat_fltrd,
  #   metrics = metrics_to_analyze,
  #   output_dir = fig_dir,
  #   hard_cutoffs = list(
  #     nCount_RNA = 150,
  #     nFeature_RNA = 150,
  #     percent.mt = 8
  #   )
  # )
  # 
  # MAD_results <- analyze_mad_metrics(
  #   seurat_obj = merged_seurat_soupx,
  #   metrics = metrics_to_analyze,
  #   output_dir = fig_dir,
  #   hard_cutoffs = list(
  #     nCount_RNA = 150,
  #     nFeature_RNA = 150,
  #     percent.mt = 8
  #   )
  # )
  
  
  
  # MAD_results <- analyze_mad_metrics(
  #   seurat_obj = merged_seurat_cb,
  #   metrics = metrics_to_analyze,
  #   output_dir = fig_dir,
  #   use_log_scale = FALSE
  # )
  qs_save(MAD_results, paste0(data_dir,'/MAD_thresholds.qs'),nthreads=nthreads)
  
  # Access results for a specific metric
  nCount_thresholds <- MAD_results$nCount_RNA$thresholds
  
  # Print thresholds if desired
  print(nCount_thresholds)
  
  
  if (!dir.exists(paste0(fig_dir,'cellbender_validation/'))) { dir.create(paste0(fig_dir,'cellbender_validation/'),
                                                        recursive=TRUE) }
  
  
  # Cellbender validation ----
  dim(merged_seurat_joined)
  # [1] 13815 51257
  
  cat('\nbuilding cellbender validation plots for',selected_timepoint,'\n')
  
  merged_seurat_joined <- Add_CellBender_Diff(seurat_object = merged_seurat_joined, raw_assay_name = "RAW", cell_bender_assay_name = "RNA")
  
  head(merged_seurat_joined@meta.data, 5)
  
  median_stats <- Median_Stats(seurat_object = merged_seurat_joined, 
                               group_by_var = "sample", 
                               median_var = c("nCount_Diff", "nFeature_Diff", 
                                              'complexity',"percent.mt", 
                                              'percent.hb','percent.rb'))
  
  feature_diff <- CellBender_Feature_Diff(seurat_object = merged_seurat_joined, raw_assay = "RAW", cell_bender_assay = "RNA")
  feature_diff$gene <- rownames(feature_diff)
  summary(feature_diff)
  # View(feature_diff)
  # genes = c("Slc27a2", "Cfd", "Depp1", "Slc4a4", "Apoe", "Adap1", "Peg3", 
  #           "Fabp3", "Marc1", "Pgm3", "Cidea", "Ppp2r5b", "Gpx3", "Tmem82",
  #           "Fads3", "Tmem120b", "Zim1", "Cdo1", "Rgs2", "Orm1", "Fam13a",
  #           "Slc25a34", "Lctl", "Hspa4l", "Aldh1a7", "Dhrs7", "Dock9", "Adcy3",
  #           "Sgpl1", "Tlcd1", "Prkce", "Zfyve21", "Bsg", "Hsph1", "Uck1",
  #           "Acss2", "Fam110c", "Bnip3", "Podn", "Ccdc69", "Acly", "Unc119",
  #           "Sorbs1", "Zc3h12a", "Mmd", "Sc5d", "Gsn", "Gk", "Raf1", "Plin4",
  #           "Ppm1b", "Stip1", "Frat2", "Sorbs3", "Apold1", "Klf2", "Tcf15",
  #           "Slc9a3r2", "Aqp1", "Dusp1", "Id1", "Bcl2l1", "Klf4", "Flt1",
  #           "Timp4", "Tcim", "Jund", "Bmp6", "Cdkn1a", "Rapgef1", "Cldn5",
  #           "Mcl1", "Kdm6b", "Esam", "Sgk1", "Cavin1", "Irak2", "Ddit4",
  #           "Fabp4", "Mcf2l", "Amotl2", "Snrk", "Tmcc3", "Map3k6", "Cldn15",
  #           "Itpkb", "Trp53i11", "Tsc22d1", "Rasgrp2", "Ahnak", "Tob1", "Csf1",
  #           "Pimreg", "Tpx2", "Prnd", "Ptms", "Ahsa2", "mt-Cytb", "Clu",
  #           "mt-Atp6", "Snx2", "Trim25", "Ucp2", "Lyve1", "Rassf2", "Dab2",
  #           "Cyth1", "Per1", "Ms4a4a", "Stk17b", "Sult1a1", "Hnrnpu", "Spry1",
  #           "Oaz1", "Ctso", "Haao", "Nab2", "Timeless", "Hsp90ab1", "Cd180",
  #           "Mt1", "Mef2c", "Rubcnl", "Tatdn1", "St13")
  p1 <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold=0.001)
  p2 <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold = 5)
  p3 <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold = 10)
  p4 <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold = 15)
  p5 <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold = 20)
  p6 <- CellBender_Diff_Plot(feature_diff_df = feature_diff, pct_diff_threshold = 25)
  
  patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 2);
  ggsave(paste0(fig_dir,'cellbender_validation/feature_diffs.png'), dpi = 320, width = 20, height = 20)
  
  
  
  # Function to compare CellBender and raw counts
  # compare_cellbender_raw <- function(seurat_obj) {
  #   # Get counts from both assays
  #   raw_counts <- LayerData(seurat_obj, assay = "RAW", layer = "counts")
  #   rna_counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
  #   
  #   # Calculate differences
  #   diff_counts <- raw_counts - rna_counts
  #   
  #   # Create a new assay for the differences
  #   seurat_obj[["DIFF"]] <- CreateAssayObject(counts = diff_counts)
  #   
  #   # Calculate summary statistics
  #   total_removed <- Matrix::colSums(diff_counts)
  #   genes_affected <- Matrix::rowSums(diff_counts > 0)
  #   
  #   # Create summary plots
  #   p1 <- ggplot(data.frame(removed = total_removed), aes(x = removed)) +
  #     geom_histogram(bins = 100) +
  #     labs(title = "Distribution of Removed Counts per Cell",
  #          x = "Number of Counts Removed",
  #          y = "Number of Cells")+
  #     theme_classic2()
  #   # ggsave(paste0(fig_dir,'cellbender_validation/removed_counts_per_cell_cellbender.png'), plot = p1, dpi = 320, width = 10, height = 10)
  #   
  #   p2 <- ggplot(data.frame(affected = genes_affected), aes(x = affected)) +
  #     geom_histogram(bins = 200) +
  #     xlim(0,2500)+
  #     ylim(0,1100)+
  #     labs(title = "Distribution of Affected Genes",
  #          x = "Number of Cells where Gene was Affected",
  #          y = "Number of Genes")+
  #     theme_classic()
  #   # ggsave(paste0(fig_dir,'cellbender_validation/n_affeced_cells_per_gene_cellbender.png'), plot = p2, dpi = 320, width = 10, height = 10)
  #   plot_grid(p1,p2)
  #   ggsave(paste0(fig_dir,'cellbender_validation/distribution_of_affected_genes_and_cells_cellbender.png'), dpi = 320, width = 10, height = 10)
  #   
  # }
  
  # compare_cellbender_raw(merged_seurat_joined)
  
  ## visualize cellbender filtering ----
  cat('\nclustering joined cellbender for',selected_timepoint,'\n')
  
  merged_seurat_joined = NormalizeData(merged_seurat_joined, assay='RNA')
  merged_seurat_joined = NormalizeData(merged_seurat_joined, assay='RAW')
  
  merged_seurat_joined = FindVariableFeatures(merged_seurat_joined)
  merged_seurat_joined = ScaleData(merged_seurat_joined)
  merged_seurat_joined = RunPCA(merged_seurat_joined)
  merged_seurat_joined = FindNeighbors(merged_seurat_joined)
  merged_seurat_joined = FindClusters(merged_seurat_joined)
  merged_seurat_joined = RunUMAP(merged_seurat_joined,
                 dims = 1:30,
                 reduction='pca')
  
  qs_save(merged_seurat_joined, paste0(data_dir, 'merged_seurat_joined_UMAP.qs'),nthreads=nthreads)
  
  # merged_seurat_joined = qs_read(paste0(data_dir, 'merged_seurat_joined_UMAP.qs'),nthreads=nthreads)
  
  grep("^Hb[^(P|E|S)]", rownames(merged_seurat_joined),value = T)
  grep("^mt-", rownames(merged_seurat_joined),value = T)
  
  # ambient RNA
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('mt-Co1','mt-Nd4l')) # mitochondrial genes
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_ambient_RNA_cleaning_mt_genes.png'), dpi = 320, width = 20, height = 10)
  
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Hbb-bt','Hbb-bs')) # hemoglobin genes
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_ambient_RNA_cleaning_hb_genes.png'), dpi = 320, width = 20, height = 10)
  
  # adipocyte markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Fabp4','Plin1'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_adipocyte_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  # b cell markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Cd19','Pax5'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_bcell_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  # endothelial markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Esam','Angptl4'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_endothelial_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  # t cell markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Lck','Lat'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_tcell_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  # fibroblast markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Col1a1','Col3a1'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_fibroblast_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  # epitheliali markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Krt7','Krt18'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_epithelial_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  # macrophage markers
  FeaturePlot_DualAssay(merged_seurat_joined, features=c('Adgre1','Csf1r'))
  ggsave(paste0(fig_dir,'cellbender_validation/cellbender_macrophage_marker_cleaning.png'), dpi = 320, width = 20, height = 10)
  
  
  # cell-level filtering with mean absolute deviations (MADs) ----
  
  # Function to identify outliers using MAD with optional manual lower bounds
  is_outlier <- function(seurat_obj, metric, nmads, sample_col = "sample", 
                         log1p_transform = FALSE, lower_bound = NULL) {
    # Input validation
    if (!metric %in% names(seurat_obj@meta.data)) {
      stop(sprintf("Metric '%s' not found in Seurat object metadata", metric))
    }
    
    # Initialize results vector
    outlier_vector <- logical(length = ncol(seurat_obj))
    names(outlier_vector) <- colnames(seurat_obj)
    
    # Process each sample separately
    for (samp in unique(seurat_obj@meta.data[[sample_col]])) {
      # Get sample indices - explicitly check the metadata
      sample_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[sample_col]] == samp, ])
      sample_idx <- colnames(seurat_obj) %in% sample_cells
      
      # Get metric values for this sample
      M <- seurat_obj@meta.data[sample_cells, metric]
      
      cat(sprintf("\nDEBUG %s - Sample %s:", metric, samp))
      cat(sprintf("\n  Number of cells: %d", length(M)))
      cat(sprintf("\n  Raw median: %.2f", median(M)))
      cat(sprintf("\n  Raw MAD: %.2f", mad(M)))
      
      
      raw_M <- M  # Keep raw values for lower bound comparison
      
      # Apply log1p transform if requested (only for MAD calculations)
      if (log1p_transform) {
        M <- log1p(M)
      }
      
      # Calculate median and MAD
      med <- median(M)
      mad_val <- mad(M)
      
      # Identify upper outliers using MAD
      upper_boundary = med + nmads * mad_val
      upper_outliers <- M > upper_boundary
      
      # Print boundary for this sample
      display_boundary <- if (log1p_transform) {
        expm1(upper_boundary)
      } else {
        upper_boundary
      }
      cat(sprintf("\n  Final upper boundary: %.2f", display_boundary))
      cat(sprintf("\n  Number of upper outliers: %d\n", sum(upper_outliers)))
      
      # Combine with manual lower bound if provided
      if (!is.null(lower_bound)) {
        lower_outliers <- raw_M < lower_bound
        outliers <- upper_outliers | lower_outliers
      } else {
        # If no lower bound, use MAD for both upper and lower
        lower_outliers <- M < med - nmads * mad_val
        outliers <- upper_outliers | lower_outliers
      }
      
      # Store results
      outlier_vector[sample_idx] <- outliers
    }
    cat("\n")  # Add newline after printing all samples
    
    return(outlier_vector)
  }
  
  # Function to apply MAD filtering to multiple metrics with manual lower bounds
  apply_mad_filters <- function(seurat_obj, 
                                metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt", "complexity"),
                                nmads = 5,
                                mt_nmads = 3,
                                mt_hard_cutoff = 8,
                                complexity_cutoff = 0.8,
                                sample_col = "sample",
                                log1p_transform = FALSE,
                                lower_bounds = list(
                                  nCount_RNA = 150,
                                  nFeature_RNA = 150,
                                  percent.mt = NULL,
                                  complexity = NULL
                                )) {
    
    # Initialize combined outlier vector
    combined_outliers <- logical(length = ncol(seurat_obj))
    
    # Process each metric
    for (metric in metrics) {
      # Get lower bound for this metric if specified
      lower_bound <- lower_bounds[[metric]]
      
      # Special handling for mitochondrial percentage - only hard cutoff
      if (metric == "percent.mt") {
        mt_outliers <- seurat_obj[[metric]] > mt_hard_cutoff
        combined_outliers <- combined_outliers | mt_outliers
        seurat_obj$mt_outlier <- mt_outliers
        
        # Print summary for mt outliers
        n_mt_outliers <- sum(mt_outliers, na.rm=T)
        cat(sprintf("\nMitochondrial outliers (> %g%%): %d cells (%.1f%%)\n", 
                    mt_hard_cutoff, n_mt_outliers, 
                    100 * n_mt_outliers/ncol(seurat_obj)))
        
        # Special handling for complexity
      } else if (metric == "complexity") {
        complexity_outliers <- seurat_obj[[metric]] < complexity_cutoff
        combined_outliers <- combined_outliers | complexity_outliers
        seurat_obj$complexity_outlier <- complexity_outliers
        
        # Print summary for complexity outliers
        n_complexity_outliers <- sum(complexity_outliers, na.rm=T)
        cat(sprintf("\nComplexity outliers (< %g): %d cells (%.1f%%)\n", 
                    complexity_cutoff, n_complexity_outliers, 
                    100 * n_complexity_outliers/ncol(seurat_obj)))
        
      } else {
        # Regular metrics
        metric_outliers <- is_outlier(seurat_obj, metric, nmads, 
                                      sample_col, log1p_transform, lower_bound)
        combined_outliers <- combined_outliers | metric_outliers
        
        # Add individual outlier status to metadata
        seurat_obj[[paste0(metric, "_outlier")]] <- metric_outliers
        
        # Print summary for this metric
        n_metric_outliers <- sum(metric_outliers)
        cat(sprintf("\n%s outliers: %d cells (%.1f%%)", 
                    metric, n_metric_outliers, 
                    100 * n_metric_outliers/ncol(seurat_obj)))
        
        if (!is.null(lower_bound)) {
          n_lower <- sum(seurat_obj[[metric]] < lower_bound)
          cat(sprintf("\n  Below lower bound (%g): %d cells\n", 
                      lower_bound, n_lower))
        }
      }
    }
    
    # Add combined outlier status to metadata
    seurat_obj$outlier <- combined_outliers
    
    # Print total outliers
    n_total_outliers <- sum(combined_outliers)
    cat(sprintf("\n\nTotal outliers: %d cells (%.1f%%)\n", 
                n_total_outliers, 100 * n_total_outliers/ncol(seurat_obj)))
    
    # Return updated Seurat object
    return(seurat_obj)
  }
  
  # Function to filter Seurat object based on outlier status
  filter_mad_outliers <- function(seurat_obj, 
                                  nmads = 5,
                                  mt_nmads = 3,
                                  mt_hard_cutoff = 8,
                                  complexity_cutoff = 0.8,
                                  sample_col = "sample",
                                  log1p_transform = FALSE,
                                  lower_bounds = list(
                                    nCount_RNA = 150,
                                    nFeature_RNA = 150
                                  )) {
    # calculate outliers
    seurat_obj = apply_mad_filters(seurat_obj,
                                   nmads=nmads,
                                   mt_nmads = mt_nmads,
                                   mt_hard_cutoff = mt_hard_cutoff,
                                   complexity_cutoff = complexity_cutoff,
                                   sample_col = sample_col,
                                   log1p_transform = log1p_transform,
                                   lower_bounds=lower_bounds)
    
    # Print summary of filtering
    n_total <- ncol(seurat_obj)
    
    n_outliers <- sum(seurat_obj$outlier)
    cat(sprintf("\nRemoving %d outlier cells (%.1f%% of total)", 
                n_outliers, 100 * n_outliers/n_total), 'from', n_total, 'cells.',
        '\n',n_total - n_outliers, 'cells remain.')
    
    # Filter object
    filtered_obj <- subset(seurat_obj, subset = outlier == FALSE)
    return(filtered_obj)
  }
  
  # Example usage:
  # seurat_obj <- apply_mad_filters(
  #   seurat_obj,
  #   metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt", "complexity"),
  #   lower_bounds = list(
  #     nCount_RNA = 100,
  #     nFeature_RNA = 100,
  #     percent.mt = NULL,
  #     complexity = NULL  # Handled separately with complexity_cutoff
  #   ),
  #   complexity_cutoff = 0.8
  # )
  # filtered_seurat <- filter_mad_outliers(seurat_obj)
  
  
  # test = apply_mad_filters(merged_seurat_cb,
  #                          nmads = 5,
  #                          lower_bounds = list(
  #                            nCount_RNA = 100,
  #                            nFeature_RNA = 100
  #                          ))
  
  # test = filter_mad_outliers(merged_seurat_cb,
  #                            nmads = 5,
  #                            lower_bounds = list(
  #                              nCount_RNA = 100,
  #                              nFeature_RNA = 100
  #                            ))
  
  
  
  # # cell-level filtering of all samples ----
  # filter_cells = function(seurat_obj, percent.mt.cutoff=10, gene_cutoff_low=150, gene_cutoff_high=5000,
  #                         umi_cutoff_low=150, umi_cutoff_high=20000, complexity_cutoff=0.8) {
  #   
  #   filtered_seurat_obj = subset(x = seurat_obj, 
  #                                subset=(nCount_RNA >= umi_cutoff_low) & 
  #                                  (nCount_RNA <= umi_cutoff_high) & 
  #                                  (nFeature_RNA >= gene_cutoff_low) &
  #                                  (nFeature_RNA <= gene_cutoff_high) &
  #                                  (complexity >= complexity_cutoff) &
  #                                  (percent.mt < percent.mt.cutoff))
  #   return(filtered_seurat_obj)
  # }
  
  
  # Find cells with NA in complexity
  na_cells <- rownames(merged_seurat_cb@meta.data)[is.na(merged_seurat_cb@meta.data$complexity)]
  
  # Look at the metadata for these specific cells
  cell_metadata <- merged_seurat_cb@meta.data[na_cells, ]
  
  # Print the metadata for these cells
  print(cell_metadata)
  
  
  ncells = length(merged_seurat_cb$orig.ident) # for tracking # cells before filtering
  
  # filtered_seurat_cb = filter_cells(merged_seurat_cb)
  # filtered_seurat_fltrd = filter_cells(merged_seurat_fltrd)
  # filtered_seurat_joined = filter_cells(merged_seurat_joined)
  
  
  cat('\nfiltering cells for',selected_timepoint,'\n')
  
  filtered_seurat_cb = filter_mad_outliers(merged_seurat_cb,
                                           nmads = 5,
                                           mt_hard_cutoff = 8,
                                           complexity_cutoff = 0.8,
                                           lower_bounds = list(
                                             nCount_RNA = 150,
                                             nFeature_RNA = 150
                                           ))
  
  filtered_seurat_cb_orig = filter_mad_outliers(merged_seurat_cb_orig,
                                           nmads = 5,
                                           mt_hard_cutoff = 8,
                                           complexity_cutoff = 0.8,
                                           lower_bounds = list(
                                             nCount_RNA = 150,
                                             nFeature_RNA = 150
                                           ))
  
  filtered_seurat_fltrd = filter_mad_outliers(merged_seurat_fltrd,
                                              nmads = 5,
                                              mt_hard_cutoff = 8,
                                              complexity_cutoff = 0.8,
                                              lower_bounds = list(
                                                nCount_RNA = 150,
                                                nFeature_RNA = 150
                                              ))
  filtered_seurat_joined = filter_mad_outliers(merged_seurat_joined,
                                              nmads = 5,
                                              mt_hard_cutoff = 8,
                                              complexity_cutoff = 0.8,
                                              lower_bounds = list(
                                                nCount_RNA = 150,
                                                nFeature_RNA = 150
                                              ))
  
  filtered_seurat_soupx = filter_mad_outliers(merged_seurat_soupx,
                                               nmads = 5,
                                               mt_hard_cutoff = 8,
                                               complexity_cutoff = 0.8,
                                               lower_bounds = list(
                                                 nCount_RNA = 150,
                                                 nFeature_RNA = 150
                                               ))
  
  
  
  ## number of cell before and after cell-level filtering ----
  ncells_post_cell_filter = length(filtered_seurat_cb$orig.ident)
  cat(paste0(ncells,' cells before cell-level filtering\n',
             ncells_post_cell_filter,' cells left after cell-level filtering\n',
             ncells - ncells_post_cell_filter,' cells removed after cell-level filtering'))
  
  # cellranger processed cells
  # 44182 cells before cell-level filtering
  # 44102 cells left after cell-level filtering
  # 80 cells removed after cell-level filtering
  
  # cellbender processed cells
  # 70205 cells before cell-level filtering
  # 54787 cells left after cell-level filtering
  # 15418 cells removed after cell-level filtering
  
  
  
  
  
  # Plot number of cells with non-zero expression per gene ----
  
  
  
  
  
  # gene-level filtering ----
  # according to best practices, want to filter out genes that
  # are not detected in at least 20 cells.
  # https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
  filter_genes = function(seurat_obj, min_cells_per_gene=3) {
    
    # extract counts from joint layers
    jointcounts = JoinLayers(seurat_obj)
    counts = LayerData(jointcounts, assay='RNA', layer='counts')
    
    # Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
    nonzero = counts > 0
    
    # Sums all TRUE values and returns TRUE if more than 20 TRUE values per gene
    keep_genes = Matrix::rowSums(nonzero) >= min_cells_per_gene
    
    # Only keeping those genes expressed in more than 20 cells
    filtered_counts = counts[keep_genes, ]
    
    # Reassign to filtered Seurat object and filtered metadata
    filtered_seurat = CreateSeuratObject(filtered_counts, meta.data = seurat_obj@meta.data)
    
    ## number of genes before and after filtering
    ngenes = sum(Matrix::rowSums(nonzero) > 0)
    ngenes_post_filtering = sum(Matrix::rowSums(nonzero) >min_cells_per_gene)
    print(paste0(ngenes,' genes before gene-level filtering.'))
    print(paste0(ngenes_post_filtering,' genes left after gene-level filtering.'))
    print(paste0(ngenes - ngenes_post_filtering,' genes removed after gene-level filtering.'))
    
    return(filtered_seurat)
    
  }
  
  ## filter genes and show number of genes before and after filtering ----
  cat('\nfiltering genes for',selected_timepoint,'\n')
  
  filtered_seurat_cb = filter_genes(filtered_seurat_cb)
  # [1] "17443 genes before gene-level filtering." ### week3
  # [1] "15044 genes left after gene-level filtering."
  # [1] "2399 genes removed after cell-level filtering."
  
  
  # [1] "17984 genes before gene-level filtering."  # week10
  # [1] "16351 genes left after gene-level filtering."
  # [1] "1633 genes removed after cell-level filtering."
  
  # [1] "17906 genes before gene-level filtering."  # 7month
  # [1] "16449 genes left after gene-level filtering."
  # [1] "1457 genes removed after cell-level filtering."
  
  # [1] "17642 genes before gene-level filtering."    # 18month
  # [1] "16577 genes left after gene-level filtering."
  # [1] "1065 genes removed after cell-level filtering."
  
  filtered_seurat_cb_orig = filter_genes(filtered_seurat_cb_orig)
  
  
  filtered_seurat_fltrd = filter_genes(filtered_seurat_fltrd)
  # [1] "17123 genes before gene-level filtering."
  # [1] "14863 genes left after gene-level filtering."
  # [1] "2260 genes removed after cell-level filtering."
  filtered_seurat_joined = filter_genes(filtered_seurat_joined)
  # [1] "13706 genes before gene-level filtering."
  # [1] "13570 genes left after gene-level filtering."
  # [1] "136 genes removed after cell-level filtering."
  
  filtered_seurat_soupx = filter_genes(filtered_seurat_soupx)
  
  
  
  ## remove any deprecated genes in dataset  ----
  remove_deprecated_genes <- function(seurat_obj) {
    # Find deprecated genes
    deprecated_genes <- grep("^DEPRECATED-", rownames(seurat_obj), value = TRUE)
    
    # If deprecated genes exist, remove them
    if (length(deprecated_genes) > 0) {
      message(paste("Removing", length(deprecated_genes), "deprecated genes from the Seurat object"))
      # Keep all genes EXCEPT the deprecated ones
      seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), deprecated_genes))
      
      message(paste('Keeping', length(rownames(seurat_obj)),'genes'))
    } else {
      message("No deprecated genes found in this Seurat object")
    }
    
    return(seurat_obj)
  }
  
  # remove the deprecated genes
  filtered_seurat_cb = remove_deprecated_genes(filtered_seurat_cb)
  filtered_seurat_cb_orig = remove_deprecated_genes(filtered_seurat_cb_orig)
  filtered_seurat_fltrd = remove_deprecated_genes(filtered_seurat_fltrd)
  filtered_seurat_joined = remove_deprecated_genes(filtered_seurat_joined)
  filtered_seurat_soupx = remove_deprecated_genes(filtered_seurat_soupx)
  
  
  
  qs_save(filtered_seurat_cb, paste0(data_dir,'/filtered_seurat_cb.qs'),nthreads=nthreads)
  qs_save(filtered_seurat_cb_orig, paste0(data_dir,'/filtered_seurat_cb_orig.qs'),nthreads=nthreads)
  qs_save(filtered_seurat_fltrd, paste0(data_dir,'/filtered_seurat_fltrd.qs'),nthreads=nthreads)
  qs_save(filtered_seurat_joined, paste0(data_dir,'/filtered_seurat_joined.qs'),nthreads=nthreads)
  qs_save(filtered_seurat_soupx, paste0(data_dir,'/filtered_seurat_soupx.qs'),nthreads=nthreads)
  
  gc()  # clean up RAM
  
  
  
  # filtered_seurat_cb = qs_read(paste0(data_dir,'filtered_seurat_cb.qs'),nthreads=nthreads)
  # filtered_seurat_fltrd = qs_read(paste0(data_dir,'filtered_seurat_fltrd.qs'),nthreads=nthreads)
  # filtered_seurat_joined = qs_read(paste0(data_dir,'filtered_seurat_joined.qs'),nthreads=nthreads)
  # filtered_seurat_soupx = qs_read(paste0(data_dir,'filtered_seurat_soupx.qs'),nthreads=nthreads)
  
  
  # doublet detection ----
  # to get scDblFinder to work:
  # installed Matrix version 1.6-5 followed by the lates irlba pacakge version 2.3.5.1 as per
  # https://github.com/plger/scDblFinder/issues/91
  
  cat('\nrunning doublet detection for',selected_timepoint,'\n')
  
  # set up paralle worker settings
  # bp <- SnowParam(
  #   workers = 4,
  #   stop.on.error = FALSE,  # This helps see partial results if some samples succeed
  #   log = TRUE  # Enable logging to see more details about what's happening
  # )
  
  sce_cb = as.SingleCellExperiment(filtered_seurat_cb)
  # sce_cb$clusters = fastcluster(sce_cb)
  sce_cb = scDblFinder(sce_cb, samples = "sample", multiSampleMode = 'split',
                       clusters = TRUE,
                       # clusters = 'clusters',
                       dbr.sd = 1, # lets the dbr be set emperically, makes setting dbr irrelevant
                       dbr.per1k = 0.008, # default is 0.008
                       BPPARAM=SnowParam(4, RNGseed=2024))
  # dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
  # which, translates to 0.05% for a 16-plex experiment according to 
  # Howitt etal, 2024, DOI: 10.1101/2024.10.03.616596
  # They calculate a fitted fraction of 0.12% for 16-plex
  # "when running multiple samples we recommend to first cluster all samples 
  # together (for example using sce$cluster <- fastcluster(sce)) and then 
  # provide the clusters to scDblFinder." It's also recommended to set dbr.sd = 1
  # for FLEX data, both according to Howitt paper and scDblFinder docs
  # https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
  
  
  register(SerialParam()) # close extra nodes/workers used
  bpstop()
  table(sce_cb$scDblFinder.class)
  gc()
  # default
  # singlet doublet 
  # 50590    5060
  
  # cluster = TRUE
  # singlet doublet 
  # 50731    4919 
  
  # clusters = TRUE
  # singlet doublet 
  # 50771    4879 
  
  
  # dbr.sd = 1
  # singlet doublet 
  # 49891    5759 
  
  # dbr.sd = 1, dbr.per1k = 0.0012
  # singlet doublet 
  # 49768    5882 
  
  
  # clusters = TRUE, dbr.sd = 1
  # singlet doublet 
  # 50542    5108 
  
  # clusters = TRUE, dbr.sd = 1, dbr.per1k = 0.0012
  # singlet doublet 
  # 50614    5036 
  
  # fastcluster() cluster = 'cluster', dbr.sd = 1, dbr.per1k = 0.0012
  # singlet doublet 
  # 50447    5203 
  
  # fastcluster() clusters = 'cluster', dbr.sd = 1, dbr.per1k = 0.0012
  # singlet doublet 
  # 50436    5214
  
  # dbr.sd = 1, dbr.per1k = 0.0012
  # singlet doublet 
  # 50436    5214
  
  
  
  png(filename = paste0(fig_dir,'histogram_scDblFinder_scores.png'), units='in', width=10, height=10, res=320)
  hist(sce_cb$scDblFinder.score, breaks = 100)
  dev.off()
  
  # hist(sce_cb$scDblFinder.weighted, breaks = 100)
  # hist(sce_cb$scDblFinder.difficulty, breaks = 100)
  # hist(sce_cb$scDblFinder.cxds_score, breaks = 100)
  # table(sce_cb$scDblFinder.mostLikelyOrigin)
  # table(sce_cb$scDblFinder.originAmbiguous)
  
  
  
  
  gc() # free RAM
  
  ### add doublets to cellbender seurat object ----
  meta_sce_cb <- sce_cb@colData@listData %>% as.data.frame() %>% 
    dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
  rownames(meta_sce_cb) = sce_cb@colData@rownames
  filtered_seurat_cb = AddMetaData(filtered_seurat_cb, metadata = meta_sce_cb %>% dplyr::select('scDblFinder.class',
                                                                                                'scDblFinder.score'))
  head(filtered_seurat_cb@meta.data)
  table(filtered_seurat_cb$scDblFinder.class)
  qs_save(sce_cb, paste0(data_dir,'sce_cb_scDblFinder.qs'),nthreads=nthreads)
  rm(list = c('meta_sce_cb', 'sce_cb')) # clean up RAM
  gc()
  qs_save(filtered_seurat_cb, paste0(data_dir,'filtered_seurat_cb.qs'),nthreads=nthreads)
  
  
  ## cellranger data doublets ----
  sce_fltrd = as.SingleCellExperiment(filtered_seurat_fltrd)
  # sce_fltrd$clusters = fastcluster(sce_fltrd)
  sce_fltrd = scDblFinder(sce_fltrd, samples = "sample", multiSampleMode = 'split',
                          clusters = TRUE,
                          # clusters = 'clusters',
                          dbr.sd = 1, # lets the dbr be set emperically
                          dbr.per1k = 0.0012, # default is 0.008
                          BPPARAM=SnowParam(4)) #dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
  register(SerialParam()) # close extra nodes/workers used
  bpstop()
  table(sce_fltrd$scDblFinder.class)
  gc()
  
  ### add doublets to filtered cellranger seurat object ----
  meta_sce_fltrd <- sce_fltrd@colData@listData %>% as.data.frame() %>% 
    dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
  rownames(meta_sce_fltrd) = sce_fltrd@colData@rownames
  filtered_seurat_fltrd = AddMetaData(filtered_seurat_fltrd, metadata = meta_sce_fltrd  %>% dplyr::select('scDblFinder.class',
                                                                                                          'scDblFinder.score'))
  head(filtered_seurat_fltrd@meta.data)
  table(filtered_seurat_fltrd$scDblFinder.class)
  qs_save(sce_fltrd, paste0(data_dir,'sce_fltrd_scDblFinder.qs'),nthreads=nthreads)
  rm(list = c('meta_sce_fltrd', 'sce_fltrd')) # clean up RAM
  gc()
  qs_save(filtered_seurat_fltrd, paste0(data_dir,'filtered_seurat_fltrd.qs'),nthreads=nthreads)
  
  
  
  ## soupx cleaned cellranger data doublets
  sce_soupx = as.SingleCellExperiment(filtered_seurat_soupx)
  # sce_soupx$clusters = fastcluster(sce_soupx)
  sce_soupx = scDblFinder(sce_soupx, samples = "sample", multiSampleMode = 'split',
                          clusters = TRUE,
                          # clusters = 'clusters',
                          dbr.sd = 1, # lets the dbr be set emperically. Makes dbr.per1k irrelevant
                          dbr.per1k = 0.0012, # default is 0.008
                          BPPARAM=SnowParam(4)) #dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
  register(SerialParam()) # close extra nodes/workers used
  bpstop()
  table(sce_soupx$scDblFinder.class)
  gc()
  
  # add doublets to filtered cellranger seurat object
  meta_sce_soupx <- sce_soupx@colData@listData %>% as.data.frame() %>% 
    dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
  rownames(meta_sce_soupx) = sce_soupx@colData@rownames
  filtered_seurat_soupx = AddMetaData(filtered_seurat_soupx, metadata = meta_sce_soupx  %>% dplyr::select('scDblFinder.class',
                                                                                                          'scDblFinder.score'))
  head(filtered_seurat_soupx@meta.data)
  table(filtered_seurat_soupx$scDblFinder.class)
  qs_save(sce_soupx, paste0(data_dir,'sce_soupx_scDblFinder.qs'),nthreads=nthreads)
  rm(list = c('meta_sce_soupx', 'sce_soupx')) # clean up RAM
  gc()
  qs_save(filtered_seurat_soupx, paste0(data_dir,'filtered_seurat_soupx.qs'),nthreads=nthreads)
  
  
  ## qc plots for doublets
  plot_doublets = function(seurat_object,fig_dir=fig_dir) {
    metrics = c("nCount_RNA", "nFeature_RNA", 'complexity',"percent.mt", "percent.hb")
    # Check how doublets and singlets differ in QC measures per sample
    cat('\nPlotting doublet metrics by sample\n')
    VlnPlot(seurat_object, group.by = 'sample', split.by = "scDblFinder.class",
            features = metrics, 
            ncol = 3, pt.size = 0) + theme(legend.position = 'right')
    ggsave(paste0(fig_dir,'/doublets_by_sample.png'), dpi=320, width=20, height= 20)
    
    # Check how doublets and singlets differ in QC measures by condition
    cat('\nPlotting doublet metrics by condition\n')
    VlnPlot(seurat_object, group.by = 'condition', split.by = "scDblFinder.class",
            features = metrics, 
            ncol = 3, pt.size = 0) + theme(legend.position = 'right')
    ggsave(paste0(fig_dir,'/doublets_by_condition.png'), dpi=320, width=20, height= 20)
  }
  
  cat('\nplotting doublet metrics for',selected_timepoint,'\n')
  plot_doublets(filtered_seurat_cb)
  
  # plot_doublets(filted_seurat_fltrd)
  # plot_doublets(filterd_seurat_soupx)
  percent.mt.cutoff = 8
  gene_cutoff_low = 150
  gene_cutoff_high = 5000
  umi_cutoff_low = 150
  umi_cutoff_high = 10000
  complexity_cutoff = 0.8 # log10genes / log10UMIs
  min_cells_per_gene = 3
  
  ## pre and post filtering plots ----
  # plot_before_and_after = function(seurat_obj, filtered_seurat_obj, percent.mt.cutoff=8, gene_cutoff_low=100, gene_cutoff_high=5000,
  #                                  umi_cutoff_low=150, umi_cutoff_high=5000, complexity_cutoff=0.8, min_cells_per_gene=3) {
  plot_before_and_after = function(seurat_obj, filtered_seurat_obj, percent.mt.cutoff=percent.mt.cutoff, gene_cutoff_low=gene_cutoff_low, gene_cutoff_high=gene_cutoff_high,
                                   umi_cutoff_low=umi_cutoff_low, umi_cutoff_high=umi_cutoff_high, complexity_cutoff=complexity_cutoff, min_cells_per_gene=min_cells_per_gene) {  
    metadata = seurat_obj@meta.data
    filtered_metadata = filtered_seurat_obj@meta.data
    
    max_y <- max(table(metadata$sample),table(filtered_metadata$sample)) # calculate max y limit
    max_y <- ceiling(max_y / 1000) * 1000 # round up to the nearest 1000
  
    p1 = ggplot(metadata, aes(x=sample, fill=sample)) + 
          geom_bar() +
          theme_classic2() +
          NoLegend()+
          labs(title="Number of cells per sample",
               subtitle = paste0('Total cells = ',dim(metadata)[1],', # ctrl cells = ',
                                 sum(metadata$condition =='ctrl'),', # pb cells = ',
                                 sum(metadata$condition == 'pb')))+
          theme(plot.title = element_text(hjust=0.5, face="bold"),
                plot.subtitle = element_text(hjust=0.5),
                axis.title = element_text(size=12), 
                axis.text = element_text(size=12),
                title = element_text(size=18))+
          scale_y_continuous(n.breaks=10, limits=c(0,max_y))
        
    p2 = ggplot(filtered_metadata, aes(x=sample, fill=sample)) + 
          geom_bar() +
          theme_classic2() +
          NoLegend()+
          labs(title="Number of cells per sample after filtering",
               subtitle = paste0('Total cells = ',dim(filtered_metadata)[1],', # ctrl cells = ',
                                 sum(filtered_metadata$condition =='ctrl'),', # pb cells = ',
                                 sum(filtered_metadata$condition == 'pb')))+
          theme(plot.title = element_text(hjust=0.5, face="bold"),
                plot.subtitle = element_text(hjust=0.5),
                axis.title = element_text(size=12), 
                axis.text = element_text(size=12),
                title = element_text(size=18))+
          scale_y_continuous(n.breaks=10, limits=c(0,max_y))
    
    plot_grid(p1, p2)
    ggsave(paste0(fig_dir,'/ncells_before_and_after_filtering.png'), dpi=320, width=20, height= 10)
    
    
    ## plot count depth (number of UMIs) vs number of genes, colored by percent mitochondrial gene counts, before and after filtering ----
    p1 = ggplot(metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
      geom_point(size=0.5)+
      DarkTheme()+
      labs(title='# of UMIs vs # of genes',
           subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
           x='# of UMI counts (log scale)',
           y='# of genes (log scale)',
           color='% mt counts')+
      scale_x_log10(n.breaks=8, labels=comma_format())+
      scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
      annotation_logticks(sides = "bl", colour='white', size=1)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18),
            plot.background=element_rect(color='black'))+
      NoGrid()+
      stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      # geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='red')+
      # geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='red')+
      scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#4b4d53')
    
    p2 = ggplot(filtered_metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
      geom_point(size=0.5)+
      DarkTheme()+
      labs(title='Count depth vs # of genes after filtering',
           subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
           x='# of UMI counts (log scale)',
           y='# of genes (log scale)',
           color='% mt counts')+
      scale_x_log10(n.breaks=8, labels=comma_format())+
      scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
      annotation_logticks(sides = "bl", colour='white', size=1)+
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18),
            plot.background=element_rect(color='black'))+
      NoGrid()+
      stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
      geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
      # geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
      geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='red')+
      # geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='red')+
      scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff,NA), na.value ='#4b4d53')
    
    plot_grid(p1, p2)
    ggsave(paste0(fig_dir,'/ncounts_vs_ngenes_coloredby_percentmt_before_and_after_filtering.png'), dpi=320, width=20, height= 10)
  }
  
  cat('\nplotting before and after filtering plots for',selected_timepoint,'\n')
  plot_before_and_after(merged_seurat_cb_orig, filtered_seurat_cb)
  
  plot_before_and_after(merged_seurat_cb_orig, filtered_seurat_cb)
  
  
  # plot_before_and_after(merged_seurat_fltrd, filtered_seurat_fltrd)
  # plot_before_and_after(merged_seurat_soupx, filtered_seurat_soupx)
  
  
  # free RAM for next timepoint loop
  rm(merged_seurat_cb,merged_seurat_cb_orig,merged_seurat_fltrd,merged_seurat_joined,
     merged_seurat_soupx, filtered_seurat_soupx, filtered_seurat_fltrd, filtered_seurat_cb,
     filtered_seurat_cb_orig, filtered_seurat_joined)
  gc()
  
  cat('\nfinished qc for',selected_timepoint,'\n')

}

#### sessionInfo() ----
sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] viridis_0.6.5               viridisLite_0.4.2           plyr_1.8.9                  ggpubr_0.6.0               
# [5] ggthemes_5.1.0              BiocParallel_1.38.0         scales_1.3.0                lubridate_1.9.4            
# [9] forcats_1.0.0               stringr_1.5.1               dplyr_1.1.4                 purrr_1.0.4                
# [13] readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1              
# [17] tidyverse_2.0.0             cowplot_1.1.3               patchwork_1.3.0             scCustomize_3.0.1          
# [21] Seurat_5.2.1                SeuratObject_5.0.2          sp_2.2-0                    scDblFinder_1.20.2         
# [25] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0 Biobase_2.64.0              GenomicRanges_1.56.1       
# [29] GenomeInfoDb_1.42.3         IRanges_2.38.1              S4Vectors_0.42.1            BiocGenerics_0.52.0        
# [33] MatrixGenerics_1.18.1       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0    bitops_1.0-9             httr_1.4.7               RColorBrewer_1.1-3       backports_1.5.0         
# [6] tools_4.4.1              sctransform_0.4.1        R6_2.6.1                 mgcv_1.9-1               lazyeval_0.2.2          
# [11] uwot_0.2.2               withr_3.0.2              gridExtra_2.3            progressr_0.15.1         cli_3.6.4               
# [16] textshaping_1.0.0        spatstat.explore_3.3-4   fastDummies_1.7.5        labeling_0.4.3           spatstat.data_3.1-4     
# [21] ggridges_0.5.6           pbapply_1.7-2            Rsamtools_2.22.0         systemfonts_1.2.1        scater_1.34.0           
# [26] parallelly_1.42.0        limma_3.60.6             rstudioapi_0.17.1        generics_0.1.3           shape_1.4.6.1           
# [31] BiocIO_1.16.0            ica_1.0-3                spatstat.random_3.3-2    car_3.1-3                Matrix_1.7-0            
# [36] ggbeeswarm_0.7.2         abind_1.4-8              lifecycle_1.0.4          yaml_2.3.10              edgeR_4.2.1             
# [41] carData_3.0-5            snakecase_0.11.1         SparseArray_1.6.1        Rtsne_0.17               paletteer_1.6.0         
# [46] grid_4.4.1               promises_1.3.2           dqrng_0.4.1              crayon_1.5.3             miniUI_0.1.1.1          
# [51] lattice_0.22-6           beachmat_2.22.0          pillar_1.10.1            metapod_1.14.0           rjson_0.2.23            
# [56] xgboost_1.7.8.1          future.apply_1.11.3      codetools_0.2-20         glue_1.8.0               spatstat.univar_3.1-1   
# [61] data.table_1.16.4        vctrs_0.6.5              png_0.1-8                spam_2.11-1              gtable_0.3.6            
# [66] rematch2_2.1.2           S4Arrays_1.4.1           mime_0.12                survival_3.8-3           statmod_1.5.0           
# [71] bluster_1.16.0           fitdistrplus_1.2-2       ROCR_1.0-11              nlme_3.1-167             RcppAnnoy_0.0.22        
# [76] irlba_2.3.5.1            vipor_0.4.7              KernSmooth_2.23-26       colorspace_2.1-1         ggrastr_1.0.2           
# [81] tidyselect_1.2.1         compiler_4.4.1           curl_6.2.1               BiocNeighbors_2.0.1      DelayedArray_0.32.0     
# [86] plotly_4.10.4            rtracklayer_1.66.0       lmtest_0.9-40            digest_0.6.37            goftest_1.2-3           
# [91] spatstat.utils_3.1-2     XVector_0.44.0           htmltools_0.5.8.1        pkgconfig_2.0.3          fastmap_1.2.0           
# [96] rlang_1.1.5              GlobalOptions_0.1.2      htmlwidgets_1.6.4        UCSC.utils_1.2.0         shiny_1.10.0            
# [101] farver_2.1.2             zoo_1.8-12               jsonlite_1.9.0           BiocSingular_1.22.0      RCurl_1.98-1.16         
# [106] magrittr_2.0.3           Formula_1.2-5            scuttle_1.16.0           GenomeInfoDbData_1.2.13  dotCall64_1.2           
# [111] munsell_0.5.1            Rcpp_1.0.14              reticulate_1.40.0        stringi_1.8.4            zlibbioc_1.50.0         
# [116] MASS_7.3-64              parallel_4.4.1           listenv_0.9.1            ggrepel_0.9.6            deldir_2.0-4            
# [121] Biostrings_2.74.1        splines_4.4.1            tensor_1.5               hms_1.1.3                circlize_0.4.16         
# [126] locfit_1.5-9.11          igraph_2.1.4             spatstat.geom_3.3-5      ggsignif_0.6.4           RcppHNSW_0.6.0          
# [131] reshape2_1.4.4           ScaledMatrix_1.14.0      XML_3.99-0.18            scran_1.34.0             ggprism_1.0.5           
# [136] tzdb_0.4.0               httpuv_1.6.15            RANN_2.6.2               polyclip_1.10-7          future_1.34.0           
# [141] scattermore_1.2          rsvd_1.0.5               janitor_2.2.1            broom_1.0.7              xtable_1.8-4            
# [146] restfulr_0.0.15          RSpectra_0.16-2          rstatix_0.7.2            later_1.4.1              ragg_1.3.3              
# [151] snow_0.4-4               beeswarm_0.4.0           GenomicAlignments_1.42.0 cluster_2.1.8            timechange_0.3.0        
# [156] globals_0.16.3   