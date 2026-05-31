library(Seurat)
library(SeuratObject)
library(scCustomize)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(scDblFinder)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(BiocParallel)

set.seed(2024)

timepoints = c('3_weeks','10_weeks','7_months','18_months')

for (time in timepoints) {
  # directories for saving
  
  output_data_dir = paste0('C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/',time,'/cellbender_analysis/FPR_0.0_MAD/data/')
  output_fig_dir = paste0('C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/',time,'/cellbender_analysis/FPR_0.0_MAD/qc/')

  # load data
  filtered_seurat_cb = readRDS(paste0(output_data_dir,'/filtered_seurat_cb.rds'))
  filtered_seurat_fltrd = readRDS(paste0(output_data_dir,'/filtered_seurat_fltrd.rds'))
  
  sce_cb = as.SingleCellExperiment(filtered_seurat_cb)
  sce_cb$cluster = fastcluster(sce_cb,rdname='PCA', ndims=50)
  sce_cb = scDblFinder(sce_cb, samples = "sample", 
                       # clusters = TRUE,
                       clusters = 'cluster',
                       dbr.sd = 1, # lets the dbr be set emperically
                       dbr.per1k = 0.0012, # default is 0.008
                       BPPARAM=SnowParam(4)) #dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
  #dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
  # which, translates to 0.05% for a 16-plex experiment according to DOI,
  # 10.1101/2024.10.03.616596
  # they calculate a fitted fraction of 0.12% for 16-plex
  # when running multiple samples we recommend to first cluster all samples 
  # together (for example using sce$cluster <- fastcluster(sce)) and then 
  # provide the clusters to scDblFinder.
  
  
  register(SerialParam()) # close extra nodes/workers used
  bpstop()
  table(sce_cb$scDblFinder.class)
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
  
  
  hist(sce_cb$scDblFinder.score, breaks = 100)
  png(filename = paste0(output_fig_dir,'histogram_scDblFinder_scores.png'), units='in', width=10, height=10, res=320)
  dev.off()
  
  hist(sce_cb$scDblFinder.cxds_score, breaks = 100)
  table(sce_cb$scDblFinder.mostLikelyOrigin)
  table(sce_cb$scDblFinder.originAmbiguous)
  
  
  
  
  gc() # free RAM
  
  # add doublets to cellbendener seurat object
  meta_sce_cb <- sce_cb@colData@listData %>% as.data.frame() %>% 
    dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
  rownames(meta_sce_cb) = sce_cb@colData@rownames
  filtered_seurat_cb = AddMetaData(filtered_seurat_cb, metadata = meta_sce_cb %>% dplyr::select('scDblFinder.class'))
  head(filtered_seurat_cb@meta.data)
  table(filtered_seurat_cb$scDblFinder.class)
  saveRDS(sce_cb, paste0(output_data_dir,'sce_cb.rds'))
  rm(list = c('meta_sce_cb', 'sce_cb')) # clean up RAM
  gc()
  saveRDS(filtered_seurat_cb, paste0(output_data_dir,'filtered_seurat_cb.rds'))
  
  # cellranger data doublets
  sce_fltrd = as.SingleCellExperiment(filtered_seurat_fltrd)
  sce_fltrd = scDblFinder(sce_fltrd, samples = "sample", dbr=0.0008, BPPARAM=SnowParam(4))
  register(SerialParam()) # close extra nodes/workers used
  bpstop()
  table(sce_fltrd$scDblFinder.class)
  gc()
  
  # add doublets to cellbedner seurat object
  meta_sce_fltrd <- sce_fltrd@colData@listData %>% as.data.frame() %>% 
    dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
  rownames(meta_sce_fltrd) = sce_fltrd@colData@rownames
  filtered_seurat_fltrd = AddMetaData(filtered_seurat_fltrd, metadata = meta_sce_fltrd  %>% dplyr::select('scDblFinder.class'))
  head(filtered_seurat_fltrd@meta.data)
  table(filtered_seurat_fltrd$scDblFinder.class)
  saveRDS(sce_fltrd, paste0(output_data_dir,'sce_fltrd.rds'))
  rm(list = c('meta_sce_fltrd', 'sce_fltrd')) # clean up RAM
  gc()
  saveRDS(filtered_seurat_fltrd, paste0(output_data_dir,'filtered_seurat_fltrd.rds'))
  
  
  ## qc plots for doublets
  plot_doublets = function(seurat_object) {
    metrics = c("nCount_RNA", "nFeature_RNA", 'complexity',"percent.mt", "percent.hb")
    # Check how doublets and singlets differ in QC measures per sample
    cat('\nPlotting doublet metrics by sample\n')
    VlnPlot(seurat_object, group.by = 'sample', split.by = "scDblFinder.class",
            features = metrics, 
            ncol = 3, pt.size = 0) + theme(legend.position = 'right')
    ggsave(paste0(output_fig_dir,'/doublets_by_sample.png'), dpi=320, width=20, height= 20)
    
    # Check how doublets and singlets differ in QC measures by condition
    cat('\nPlotting doublet metrics by condition\n')
    VlnPlot(seurat_object, group.by = 'condition', split.by = "scDblFinder.class",
            features = metrics, 
            ncol = 3, pt.size = 0) + theme(legend.position = 'right')
    ggsave(paste0(output_fig_dir,'/doublets_by_condition.png'), dpi=320, width=20, height= 20)
    
    plot_doublets(filtered_seurat_cb)
  }
}
  
  
