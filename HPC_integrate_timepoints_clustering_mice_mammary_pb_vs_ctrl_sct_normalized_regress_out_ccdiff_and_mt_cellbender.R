library(Seurat)
library(SeuratObject)
library(sctransform)
library(harmony)
library(patchwork)
library(cowplot)
library(tidyverse)
# library(Matrix)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(glmGamPoi)
# library(metap)
# library(RCurl)
# library(Polychrome)
library(pheatmap)
library(presto)
library(qs2) # for fast data saving/loading
library(BiocParallel)
library(scDblFinder)

# set number of threads to use for saving qs formatted data
nthreads = 4 

# set seed
set.seed(2024)

# create polychrome color palette
# p12 = createPalette(12, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
# swatch(p12)
# names(p12) = NULL

# sessionInfo()

# read in data ----

# analysis type, cellbender or cellranger
analysis <- 'cellbender_analysis'
# analysis <- 'cellranger.v.9.0.0_analysis'

# Set timepoint to analyze for this script
# selected_timepoint <- 'week3'
# selected_timepoint <- 'week10'
# selected_timepoint <- 'month7'
# selected_timepoint <- 'month18'


## directories for timepoint data ----
timepoint_dirs <- list(
  "week3" = list(
    fig_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/annotation/")
  ),
  "week10" = list(
    fig_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/annotation/")
  ),
  "month7" = list(
    fig_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/7_months/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/7_months/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/7_month/",analysis,"/FPR_0.0_MAD/annotation/")
  ),
  "month18" = list(
    fig_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/annotation/")
  )
)
# timepoint_dirs <- list(
#   "week3" = list(
#     fig_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   ),
#   "week10" = list(
#     fig_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   ),
#   "month7" = list(
#     fig_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/7_month/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   ),
#   "month18" = list(
#     fig_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   )
# )
# list2env(timepoint_dirs[[selected_timepoint]], envir = environment())

## directories for integrated data ----
fig_dir_integrated = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir_integrated = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir_integrated = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")

# Create directories to save results if they don't already exist:
dirs = c(fig_dir_integrated, data_dir_integrated, annot_dir_integrated)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# load in Qc'd data ----
seurat.list = list()
for (timepoint in names(timepoint_dirs)) {
  data_dir = timepoint_dirs[[timepoint]]$data_dir
  integrated_seurat = qs_read(paste0(data_dir,'filtered_seurat_cb.qs'),nthreads=nthreads)
  seurat.list[[timepoint]] = integrated_seurat
}
rm(integrated_seurat, data_dir) # free RAM
gc()

## merge timepoints ----
integrated_seurat = merge(seurat.list[[1]], seurat.list[-1])
integrated_seurat = JoinLayers(integrated_seurat)

## add ordinal timepoints metadata column ----
integrated_seurat$timepoint = str_extract(integrated_seurat$cells, '^[^_]+')
integrated_seurat$timepoint = factor(integrated_seurat$timepoint,
                                     levels=c('wk3','wk10','mn7','mn18'),
                                     ordered = TRUE)


rm(seurat.list) # free RAM
gc()

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat.qs'), nthreads=nthreads)

# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat.qs'), nthreads=nthreads)


# exploring sources of variation, cell cycle and mitochondrial gene expression----

# cell cycle genes, formatted to match genes in matrix with first letter upper-cased
s.genes = stringr::str_to_title(cc.genes$s.genes)
g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)


# Normalize the counts
integrated_seurat = NormalizeData(integrated_seurat)



# Score cells for cell cycle
integrated_seurat = CellCycleScoring(integrated_seurat, 
                                    g2m.features = g2m.genes, 
                                    s.features = s.genes,
                                    set.ident=T)

# From https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Seurat suggests regressing out the difference between S and G2M phase scores
# in tissues with differentiating cells, like our developing mammary gland.
integrated_seurat$cc.difference = integrated_seurat$S.Score - integrated_seurat$G2M.Score


# Check quartile values of mitochondrial gene expression
# summary(integrated_seurat@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1548  0.7305  1.2666  1.9231  8.0000 

# Turn percent.mt into categorical factor vector based on quartile values
integrated_seurat@meta.data$mito.level = cut(integrated_seurat@meta.data$percent.mt, 
                                        breaks=c(-Inf, 2, 4, 6, Inf), 
                                        labels=c("1st","2nd","3rd", "4th"),
                                        ordered_result=TRUE)


# Identify the most variable genes
seurat_pca = FindVariableFeatures(integrated_seurat, 
                                  selection.method = "vst", # default
                                  nfeatures = 3000, # default is 2000
                                  verbose = FALSE)

# Identify the 25 most highly variable genes
top25 = head(VariableFeatures(seurat_pca), 25)
# plot variable features with and without labels
p = VariableFeaturePlot(seurat_pca)
LabelPoints(plot=p, points = top25, repel = TRUE)+
  labs(title='Top 25 variably expressed genes')+
  theme_classic2()
ggsave(paste0(fig_dir_integrated,'top25_variable_genes_integrated_seurat.png'), dpi=320, width=20, height= 10)

# test = qs_read('/nfs/turbo/sph-colacino/aguilada/scRNAseq/8651-JC-mammary/3_weeks/data/annotated_seurat_res.0.8_df.qs',nthreads=nthreads)


# Scale the counts
seurat_pca = ScaleData(seurat_pca)

# Perform PCA on merged non-integrated data using cell cycle genes ----
seurat_pca = RunPCA(seurat_pca, features=c(s.genes,g2m.genes), npcs=50) # default params
gc() # free RAM

#### Plot the PCA colored by cell cycle phase using cell cycle genes ----
DimPlot(seurat_pca,
        reduction = "pca",
        group.by = "Phase",
        split.by = 'Phase')+
  labs(title='PCA by cell cycle phase',
       color='phase')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir_integrated,'/PCA_by_cellcycle_phase_genes_no_regressed_variables.png'), dpi=320, width=20, height= 10)
# there is a lot of variation between S and G2M phase. 
# Perhaps better to regress this out in our developing mammary tissue


# From https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Regress out the difference between S and G2M phase scores
# takes a long while to regress out. Pushing up on 20-30 mins
# seurat_pca = ScaleData(seurat_pca, vars.to.regress = "cc.difference", features = rownames(seurat_pca))
# gc() # free RAM
# qs_save(seurat_pca, paste0(data_dir,'seurat_pca_regressout_cc.diff.qs'),nthreads=nthreads)

# run PCA again w/ cc.diff now regressed out
# seurat_pca = RunPCA(seurat_pca, features=c(s.genes,g2m.genes), npcs=50) # default params

## What affect did regressing out the cc.diff have?
# DimPlot(seurat_pca,
#         reduction = "pca",
#         group.by= "Phase",
#         split.by='Phase')+
#   labs(title='PCA by cell cycle phase',
#        color='phase')+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.title = element_text(hjust=0.5, face="bold"),
#         plot.subtitle = element_text(hjust=0.5),
#         axis.title = element_text(size=12), 
#         axis.text = element_text(size=12),
#         title = element_text(size=18))
# ggsave(paste0(fig_dir,'/PCA_by_cellcycle_phase_genes_ccdiff_regressed_out.png'), dpi=320, width=20, height= 10)
# still some variability between phases

# run PCA on all genes with cc.diff regressed out
# seurat_pca = RunPCA(seurat_pca, npcs=50) # default params

# run pca on all features
seurat_pca = RunPCA(seurat_pca, npcs=50)
gc() # free RAM


DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "Phase",
        split.by='Phase',
        shuffle=T)+
  labs(title='PCA by cell cycle phase',
       color='phase')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir_integrated, '/PCA_by_cellcycle_phase_allgenes_.png'), dpi=320, width=20, height= 10)


### Plot the PCA colored by mitochondrial expression quartile ----
DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "mito.level",
        split.by = "mito.level")+
  labs(title='PCA by mitochondrial expression quartile',
       color='quartile')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir_integrated,'/PCA_by_mitochondrial_expression_quartile_all_genes.png'), dpi=320, width=20, height= 10)

qs_save(seurat_pca, paste0(data_dir_integrated,'integrated_seurat_rm29and32_regressout_cc.diff.qs'),nthreads=nthreads)# save integrated_seurat object
##looking at the plots there doesn't seem to be any significant source of variation from cell cycle phase or mitochondrial expression

### plot the PCA colored by condition ----
DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "condition",
        split.by='condition')+
  labs(title='PCA by condition',
       color='condition')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir_integrated,'/PCA_by_condition_all_genes.png'), dpi=320, width=20, height= 10)


rm(seurat_pca,p) # free RAM. Don't need any further
gc()


# no longer need data log-normalized layer, as we will normalize with SCT
integrated_seurat = DietSeurat(integrated_seurat,assays='RNA',layers='counts')
gc() # free RAM


qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat.qs'), nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat.qs'),nthreads=nthreads)
gc() # free RAM

# Split integrated_seurat by sample and apply sctransform to each individually ----
# Split seurat object by sample to perform cell cycle scoring and SCTransform on all samples
## will normalize integrated_seurat with SCTransform as don't need PCA exploration from seurat_pca object.
split_seurat = SplitObject(integrated_seurat, split.by = "sample")

gc() # save RAM

# SCTransform ----
# adjust object size limits of variables/objects in R
options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB

#### apply SCTransform to each of the 10 samples, then merge normalized samples and save object ----
# Do for each sample separately according to https://github.com/satijalab/seurat/issues/6003
# return all genes/features according to https://github.com/satijalab/seurat/issues/5205

for(i in seq_along(split_seurat)) { 
  cat('\n Beginning sample:', names(split_seurat[i]), '\n\n')
  # split_seurat[i] = NormalizeData(split_seurat[i])
  # split_seurat[i] = FindVariableFeatures(split_seurat[i])
  # split_seurat[i] = ScaleData(split_seurat[i])
  split_seurat[i] = SCTransform(split_seurat[[i]], assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                  vars.to.regress=c('cc.difference','percent.mt'),
                  return.only.var.genes=FALSE)
  # cat('/n running PCA for,', names(split_seurat[i]), '\n\n')
  # split_seurat[i] = RunPCA(split_seurat[[i]], assay='SCT', npcs=50)
  # split_seurat[i] = FindNeighbors(split_seurat[[i]], dims=1:50)
  # split_seurat[i] = FindClusters(object = split_seurat[[i]])
  # split_seurat[i] = RunUMAP(split_seurat[[i]], dims=1:50, reduction='pca')
  
}
gc() # free RAM

qs_save(split_seurat, paste0(data_dir_integrated,'split_seurat_ct_regressout_cc.diffandmt.qs'),nthreads=nthreads)
gc()

# split_seurat = qs_read(paste0(data_dir_integrated,'split_seurat_ct_regressout_cc.diffandmt.qs'),nthreads=nthreads)
# gc()

# merge sctransform normalized samples. Takes a little while to run
integrated_seurat = merge(split_seurat[[1]], split_seurat[-1], merge.data=TRUE,merge.dr=FALSE)
integrated_seurat = JoinLayers(integrated_seurat, assay='RNA')
gc() # free RAM


qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)

# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)
gc() #free up RAM

### continue analysis without integration using sctransform normalization ----
# select integration features according to https://github.com/satijalab/seurat/issues/5205
# and https://github.com/satijalab/seurat/issues/6185
# selecting 5000 features to account for the possible different rankings of highly variable genes
# between samples. Since we're normalizing with SCTtransform this # shouldn't change things much.
features = SelectIntegrationFeatures(split_seurat, nfeatures=5000, fvf.nfeatures=5000)
qs_save(features,paste0(data_dir_integrated,'variable_features.qs'),nthreads=nthreads)
VariableFeatures(integrated_seurat) = features

sct_results = SCTResults(integrated_seurat, slot = "feature.attributes")
qs_save(sct_results, paste0(data_dir_integrated,'sct_results_regressout_cc.diffandmt.qs'),nthreads=nthreads)
rm(sct_results)

# Variable features lost/unset after merging transformed objects. Need to set again.
# should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(integrated_seurat[["SCT"]]) = rownames(integrated_seurat[["SCT"]]@scale.data)

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)

# can plot variable genes by plotting gmean vs residual variance
# according to https://github.com/satijalab/seurat/issues/8321
# z = integrated_seurat@assays$SCT@SCTModel.list[[3]]@feature.attributes
# qplot(z$gmean, z$residual_variance) + scale_y_log10() + scale_x_log10()

# Identify the 25 most highly variable genes
top25_sct = head(VariableFeatures(integrated_seurat), 25)
top25_sct

sum(top25 %in% top25_sct) # = 0. completely different top 25 variable features between sctransform and log norm

# VariableFeatures() won't work for merged object after SCTtransform.
# Need to run for each sample individually
lapply(seq_along(split_seurat), function(x,n,i) {
  top25_sct = head(VariableFeatures(x[[i]]), 25)
  p = VariableFeaturePlot(x[[i]], selection.method='sct', assay='SCT') 
  LabelPoints(plot=p, points = top25_sct, repel = TRUE)+
    labs(title=paste0('Top 25 variably expressed genes for ', n[[i]]))+
    theme_classic2()
  ggsave(paste0(fig_dir_integrated,'top25_variable_genes_',n[[i]],'.png'), dpi=320, width=20, height= 10)
  
}, x=split_seurat, n=names(split_seurat))

rm(split_seurat) # free RAM
gc()


# split by timepoint so we can see batch effects and integrate them out later
# integrated_seurat[['SCT']] = split(integrated_seurat[["SCT"]], f = integrated_seurat$timepoint)

# Run PCA on non-integrated merged SCTransformed data ---- 
integrated_seurat = RunPCA(integrated_seurat, assay='SCT', npcs=50) # default params
gc() # free up RAM


qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)
  

### Visualize expression of most highly weighted genes for PCs 1-50 (loadings) ----
VizDimLoadings(integrated_seurat, dims = 1:12, ncol=4)
ggsave(paste0(fig_dir_integrated,'PCA_loadings_sct_1-12.png'), dpi=320, width=20, height= 15)

VizDimLoadings(integrated_seurat, dims = 13:24, ncol=4)
ggsave(paste0(fig_dir_integrated,'PCA_loadings_sct_13-24.png'), dpi=320, width=20, height= 15)

VizDimLoadings(integrated_seurat, dims = 25:36, ncol=4)
ggsave(paste0(fig_dir_integrated,'PCA_loadings_sct_25-36.png'), dpi=320, width=20, height= 15)

VizDimLoadings(integrated_seurat, dims = 37:43, ncol=4)
ggsave(paste0(fig_dir_integrated,'PCA_loadings_sct_37-43.png'), dpi=320, width=20, height= 15)

VizDimLoadings(integrated_seurat, dims = 44:50, ncol=4)
ggsave(paste0(fig_dir_integrated,'PCA_loadings_sct_44-50.png'), dpi=320, width=20, height= 15)

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC1-6_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 1:6,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC7-12_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 7:12,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC13-18_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 13:18,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC19-24_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 19:24,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC25-30_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 25:30,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC31-36_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 31:36,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC37-42_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 37:42,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir_integrated,'high_score_genes_per_PC43-50_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = integrated_seurat, 
           reduction = "pca", 
           dims = 43:50,
           ncol=4,
           cells=500,
           balanced = TRUE)
dev.off()

# ranked list (PC score) of genes for PC1
names(sort(Loadings(integrated_seurat, reduction="pca")[,1], decreasing=TRUE))


# Printing out the top 10 most variable genes driving PCs
print(x = integrated_seurat[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

### Plot the elbow plot to see variation captured by PCs ----
ElbowPlot(object = integrated_seurat, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir_integrated,'elbow_plot_sct.png'), dpi=320, width=20, height= 10)

# We can calculate where the principal components start to elbow by taking the smaller value of:
#
# 1)The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# 2)The point where the percent change in variation between the consecutive PCs is less than 0.1%.

### Determine percent of variation associated with each PC ----
pct = integrated_seurat[["pca"]]@stdev / sum(integrated_seurat[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 42
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 12
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 12

# Create a dataframe to visualize percent variation captured by PCs
plot_df = data.frame(pct = pct, 
                     cumu = cumu, 
                     rank = 1:length(pct))

### Elbow plot to visualize percent variation captured by PCs ----
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_classic2()+
  scale_x_continuous(n.breaks=10)+
  scale_y_continuous(n.breaks=10)+
  labs(title = 'Eblow Plot of percent variation explained by PCs',
       x='Cumulative percentage',
       y='Percent of variation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir_integrated,'quantitative_elbow_plot_sct.png'), dpi=320, width=20, height= 10)
rm(plot_df)


# cluster sctransformed non-integrated cells----
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)

# Run PCA on non-integrated merged SCTransformed data ---- 
# integrated_seurat = RunPCA(integrated_seurat, assay='SCT', npcs=50) # default params
gc() # free up RAM

resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)


integrated_seurat = FindNeighbors(integrated_seurat, dims=1:50, reduction='pca')
integrated_seurat = FindClusters(object = integrated_seurat, algorithm=4,
                                 resolution = resolutions)
#### run UMAP on non-integrated, merged, sct-normalized data ----
integrated_seurat = RunUMAP(integrated_seurat,
                            dims=1:50,
                            reduction='pca') 


qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_umap_non-integrated.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_umap_non-integrated.qs'),nthreads=nthreads)
gc() # free RAM


DimPlot(integrated_seurat,
        reduction = "umap",group.by = 'timepoint',
        label=F,
        label.size=6,
        shuffle=TRUE )

DimPlot(integrated_seurat,
        reduction = "umap",group.by = 'Phase',split.by = 'Phase',
        label=F,
        label.size=6,
        shuffle=TRUE )

DimPlot(integrated_seurat,
        reduction = "umap",group.by = 'condition', split.by='condition',
        label=F,
        label.size=6,
        shuffle=TRUE)

DimPlot(integrated_seurat,
        reduction = "umap",group.by = 'mito.level',split.by='mito.level',
        label=F,
        label.size=6,
        shuffle=TRUE)


Idents(integrated_seurat) = 'SCT_snn_res.0.1'
DimPlot(integrated_seurat,
        reduction = "umap",
        label=T,
        label.size=6)


p1 = DimPlot(integrated_seurat,
        reduction = "umap",group.by = 'timepoint',
        label=F,
        label.size=6,
        shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
        reduction = "umap",group.by='SCT_snn_res.0.1',
        label=T,
        label.size=6)
p1 + p2

features = c('Adipoq','Krt17','Sdc4','Krt19','Cd8a','Cd4','Pecam1',
             'Col1a1')

FeaturePlot(integrated_seurat,
            reduction='umap',
            features='G2M.Score',
            pt.size = 6,
            slot = 'data')

## Try integration w/ harmony across samples ----
integrated_seurat = IntegrateLayers(
  object = integrated_seurat,
  method = HarmonyIntegration,
  orig.reduction = 'pca',
  new.reduction = 'harmony',
  assay='SCT',
  normalization.method = "SCT",
  verbose = T
)


integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "harmony",
                                   graph.name=c('harmony_nn','harmony_snn'))
integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
                                  graph.name='harmony_snn')
integrated_seurat = RunUMAP(integrated_seurat,
                            dims=1:50,
                            reduction='harmony',
                            reduction.name = 'umap_harmony')

qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_sample.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_sample.qs'),nthreads=nthreads)
gc()

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'timepoint',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 # looks okay, could be better

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'scDblFinder.class',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'Phase',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 



# integrate across timepoint/flowcell and cluster sctransformed integrated cells ----

# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)

# split by timepoint so we can see batch effects and integrate them out later
# integrated_seurat[['RNA']] = split(integrated_seurat[["RNA"]], f = integrated_seurat$timepoint)


# adjust object size limits of variables/objects in R
# options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
# integrated_seurat = SCTransform(integrated_seurat)
# gc() # free RAM
# integrated_seurat = SCTransform(integrated_seurat,vst.flavor = 'vs2',
#                                 vars.to.regress=c('percenty.mt','cc.difference'),
#                                 return.only.var.genes=FALSE)

# integrated_seurat = RunPCA(integrated_seurat, assay='SCT',npcs=50)


# integrated_seurat = IntegrateLayers(
#   object = integrated_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = 'pca',
#   new.reduction = 'harmony',
#   assay='SCT',
#   normalization.method = "SCT",
#   verbose = T
# )

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "integrated.dr",
#                                  graph.name = c('idr_nn','idr_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                 graph.name = 'idr_snn')
# integrated_seurat = RunUMAP(integrated_seurat,
#                           dims=1:50,
#                           reduction='integrated.dr',
#                           reduction.name = 'umap_idr')

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "harmony",
#                                    graph.name=c('harmony_nn','harmony_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                   graph.name='harmony_snn')
# integrated_seurat = RunUMAP(integrated_seurat,
#                             dims=1:50,
#                             reduction='harmony',
#                             reduction.name = 'umap_harmony')
# 
# gc() # free RAM
# 
# qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_timepoint.qs'),nthreads=nthreads)
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'timepoint',split.by='timepoint',
#         label=F,
#         label.size=6,
#         shuffle=TRUE )
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'Phase', split.by='Phase',
#         label=F,
#         label.size=6,
#         shuffle=TRUE )     
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'sample',
#         label=F,
#         label.size=6,
#         shuffle=TRUE )  
# 
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'condition',
#         split.by='condition',
#         label=F,
#         label.size=6,
#         shuffle=TRUE)
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'mito.level',split.by='mito.level',
#         label=F,
#         label.size=6,
#         shuffle=TRUE)
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'cell_calls', split.by='cell_calls',
#         # split.by = 'cell_calls',
#         label=F,
#         label.size=6,
#         shuffle=TRUE)
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'scDblFinder.class',
#         # split.by = 'scDblFinder.class',
#         label=F,
#         label.size=6,
#         shuffle=TRUE)
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'harmony_snn_res.0.1',
#         label=T,
#         label.size=6,
#         shuffle=TRUE)
# 
# 
# DimPlot(integrated_seurat,
#         reduction = "umap_harmony",group.by = 'harmony_snn_res.0.2',
#         label=T,
#         label.size=6,
#         shuffle=TRUE)
# 
# 
# 
# 
# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'timepoint',
#              label=F,
#              label.size=6,
#              shuffle=TRUE )
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 
# 
# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'scDblFinder.class',
#              label=F,
#              label.size=6,
#              shuffle=TRUE )
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 
# 
# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'Phase',
#              label=F,
#              label.size=6,
#              shuffle=TRUE )
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 



# integrate across batch/experiment and cluster sctransformed integrated cells ----
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)
integrated_seurat$batch = ifelse(integrated_seurat$timepoint == 'wk3','Experiment_1','Experiment_2')
integrated_seurat$timepoint = factor(integrated_seurat$timepoint,
                                     levels=c('wk3','wk10','mn7','mn18'),
                                     ordered = TRUE)

# # split by batch so we can see batch effects and integrate them out later
# integrated_seurat[['RNA']] = split(integrated_seurat[["RNA"]], f = integrated_seurat$batch)


# adjust object size limits of variables/objects in R
# options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
# integrated_seurat = SCTransform(integrated_seurat)
# gc() # free RAM
# integrated_seurat = SCTransform(integrated_seurat,vst.flavor = 'vs2',
#                                 vars.to.regress=c('percenty.mt','cc.difference'),
#                                 return.only.var.genes=FALSE)
# gc()

# integrated_seurat = RunPCA(integrated_seurat, assay='SCT',npcs=50)


integrated_seurat = RunHarmony(integrated_seurat, 'batch')


# integrated_seurat = IntegrateLayers(
#   object = integrated_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = 'pca',
#   new.reduction = 'harmony',
#   assay='SCT',
#   normalization.method = "SCT",
#   verbose = T
# )

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "integrated.dr",
#                                  graph.name = c('idr_nn','idr_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                 graph.name = 'idr_snn')
# integrated_seurat = RunUMAP(integrated_seurat,
#                           dims=1:50,
#                           reduction='integrated.dr',
#                           reduction.name = 'umap_idr')

integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "harmony",
                                   graph.name=c('harmony_nn','harmony_snn'))
resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
                                  algorithm=4, graph.name='harmony_snn')
integrated_seurat = RunUMAP(integrated_seurat,
                            dims=1:50,
                            reduction='harmony',
                            reduction.name = 'umap_harmony')

gc() # free RAM

qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_batch.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_batch.qs'),nthreads=nthreads)
gc()

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'timepoint',
             label=F,
             label.size=6,
             shuffle=TRUE,
             cols=hue_pal()(length(unique(integrated_seurat$timepoint))))

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'batch',
             label=F,
             label.size=6,
             shuffle=TRUE)

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.25',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'scDblFinder.class',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'mito.level',
             label=F,
             label.size=6,
             shuffle=TRUE)

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'Phase',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'cell_calls',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'mito.level',
             label=F,
             label.size=6,
             shuffle=TRUE)


p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'cell_calls',
             label=F,
             label.size=6,
             shuffle=TRUE )
p1 + p2 

DimPlot(integrated_seurat,
        reduction = "umap_harmony",group.by = 'mito.level',
        split.by = 'timepoint',
        label=F,
        label.size=6,
        shuffle=TRUE )

DimPlot(integrated_seurat,
        reduction = "umap_harmony",group.by = 'condition',
        split.by = 'timepoint',
        label=F,
        label.size=6,
        shuffle=TRUE )



features = c('Adipoq','Krt17','Sdc4','Krt19','Cd8a','Cd4','Pecam1',
             'Col1a1')

FeaturePlot(integrated_seurat, 
            features = c('Col3a1'),
            reduction='umap_harmony',
            min.cutoff='q10',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

# integrate across batch/experiment and timepoint and cluster sctransformed integrated cells ----
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)
# integrated_seurat$batch = ifelse(integrated_seurat$timepoint == 'wk3','Experiment_1','Experiment_2')
# integrated_seurat$timepoint = factor(integrated_seurat$timepoint,
#                                      levels=c('wk3','wk10','mn7','mn18'),
#                                      ordered = TRUE)

# # split by batch so we can see batch effects and integrate them out later
# integrated_seurat[['RNA']] = split(integrated_seurat[["RNA"]], f = integrated_seurat$batch)


# adjust object size limits of variables/objects in R
# options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
# integrated_seurat = SCTransform(integrated_seurat)
# gc() # free RAM
# integrated_seurat = SCTransform(integrated_seurat,vst.flavor = 'vs2',
#                                 vars.to.regress=c('percenty.mt','cc.difference'),
#                                 return.only.var.genes=FALSE)
# gc()

# integrated_seurat = RunPCA(integrated_seurat, assay='SCT',npcs=50)


# integrated_seurat = RunHarmony(integrated_seurat, c('batch','timepoint'))


# integrated_seurat = IntegrateLayers(
#   object = integrated_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = 'pca',
#   new.reduction = 'harmony',
#   assay='SCT',
#   normalization.method = "SCT",
#   verbose = T
# )

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "integrated.dr",
#                                  graph.name = c('idr_nn','idr_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                 graph.name = 'idr_snn')
# integrated_seurat = RunUMAP(integrated_seurat,
#                           dims=1:50,
#                           reduction='integrated.dr',
#                           reduction.name = 'umap_idr')

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "harmony",
                                   # graph.name=c('harmony_nn','harmony_snn'))

# library(reticulate)
# py_config()  # This shows the Python configuration that R is using
# leidenalg_available <- py_module_available("leidenalg")
# print(paste("leidenalg available:", leidenalg_available))
# py_require("leidenalg")

# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolution,
#                                   graph.name='harmony_snn')

# resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
# integrated_seurat <- FindClusters(integrated_seurat, random.seed=1, resolution = resolutions,
#                                   algorithm=4, graph.name='harmony_snn') # use leiden clustering
# 
# integrated_seurat = RunUMAP(integrated_seurat,
#                             dims=1:50,
#                             reduction='harmony',
#                             reduction.name = 'umap_harmony')
# 
# gc() # free RAM

# qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_batch_and_timepoint_leiden.qs'),nthreads=nthreads)
# gc()

# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_batch_and_timepoint_leiden.qs'),nthreads=nthreads)

# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'timepoint',
#              label=F,
#              label.size=6,
#              shuffle=TRUE)
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 
# 
# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'batch',
#              label=F,
#              label.size=6,
#              shuffle=TRUE)
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 
# 
# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'scDblFinder.class',
#              label=F,
#              label.size=6,
#              shuffle=TRUE )
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 
# 
# p1 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by = 'Phase',
#              label=F,
#              label.size=6,
#              shuffle=TRUE )
# 
# p2 = DimPlot(integrated_seurat,
#              reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
#              label=T,
#              label.size=6)
# p1 + p2 
# 
# 
# features = c('Adipoq','Krt17','Sdc4','Krt19','Cd8a','Cd4','Pecam1',
#              'Col1a1')
# 
# FeaturePlot(integrated_seurat, 
#             features = c('Adipoq'),
#             reduction='umap_harmony',
#             min.cutoff='q10',
#             max.cutoff='q90',
#             cols=viridis(256)) &
#   DarkTheme()


## doublets QC ----
umap_doublets_qc = function(seurat_umap, reduction = 'umap',resolution, directory) {
  
  if (!dir.exists(paste0(directory, resolution,'/'))) {dir.create(paste0(directory, resolution,'/'),
                                                                  recursive=TRUE)}
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  reduction_name = toupper(reduction)
  
  # set umap resolution
  Idents(seurat_umap) = resolution
  seurat_umap@meta.data$resolution = Idents(seurat_umap)
  
  
  
  # plot side by side UMAP of clusters and doublet classification
  p1 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = resolution,
               label=T,
               repel=T,
               label.size=6,
               shuffle=TRUE)
  
  p2 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = 'scDblFinder.class',
               label=F,
               label.size=6,
               shuffle=TRUE)
  p1 | p2
  ggsave(filename=paste0(directory,resolution,'/',reduction_name,'_doublets.png'), dpi=320, width=20, height=10)
  
  
  # Extract identity and sample information from seurat object to determine 
  # the number of cells per cluster per sample
  n_cells = FetchData(seurat_umap, 
                      vars = c("ident", "scDblFinder.class")) %>%
    dplyr::count(ident, scDblFinder.class)
  
  # barplot of singlets vs doublets by cluster
  ggplot(n_cells, aes(x=ident, y=n, fill=scDblFinder.class))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_scDblFinder.png'), dpi=320, width=20, height= 10)
  
  # Barplot of proportion of doublets by cluster
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=scDblFinder.class), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5),
          legend.key.size = unit(1,'cm'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20))
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_scDblFinder.png'), dpi=320, width=20, height= 10)
  
}

### plot doublets on UMAP ----
# umap_doublets_qc(integrated_seurat, resolution='SCT_snn_res.0.8', directory=fig_dir_integrated)
# umap_doublets_qc(integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.2', directory=fig_dir_integrated)
umap_doublets_qc(integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.1', directory=fig_dir_integrated)
umap_doublets_qc(integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.25', directory=fig_dir_integrated)


### subset doublets out of dataset ----

# tabulate doublets by sample
table(integrated_seurat$sample, integrated_seurat$scDblFinder.class)


# tabulate total doublets
table(integrated_seurat$scDblFinder.class)

integrated_seurat = subset(integrated_seurat, subset=(scDblFinder.class == 'singlet'))
table(integrated_seurat$scDblFinder.class)

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_leiden_UMAP_doublets_removed_batch_integrated.qs'),nthreads=nthreads)
gc()

## re-do PCA and harmony integration with doublets removed, clustering on singlets ----

# # adjust object size limits of variables/objects in R
# options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
# integrated_seurat = SCTransform(integrated_seurat)
#                                 # vars.to.regress=c('cc.difference','percent.mt'))


integrated_seurat = RunPCA(integrated_seurat, assay='SCT',npcs=50,verbose=T)

integrated_seurat = RunHarmony(integrated_seurat, c('batch'))

# gc()
# 
# integrated_seurat = IntegrateLayers(
#   object = integrated_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = 'pca',
#   new.reduction = 'harmony',
#   assay='SCT',
#   normalization.method = "SCT",
#   verbose = T
# )

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "integrated.dr",
#                                  graph.name = c('idr_nn','idr_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                 graph.name = 'idr_snn')
# integrated_seurat = RunUMAP(integrated_seurat,
#                           dims=1:50,
#                           reduction='integrated.dr',
#                           reduction.name = 'umap_idr')

integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "harmony",
                                   graph.name=c('harmony_nn','harmony_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                   graph.name='harmony_snn')
resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

integrated_seurat <- FindClusters(integrated_seurat, random.seed=1, resolution = resolutions,
                                  algorithm=4, graph.name='harmony_snn') # use leiden clustering
integrated_seurat = RunUMAP(integrated_seurat,
                            dims=1:50,
                            reduction='harmony',
                            reduction.name = 'umap_harmony')

gc() # free RAM

integrated_seurat = JoinLayers(integrated_seurat,assay='RNA')

qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden.qs'),nthreads=nthreads)

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'timepoint',
             label=F,
             label.size=6,
             shuffle=TRUE,
             cols=hue_pal()(length(unique(integrated_seurat$timepoint)))
             )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.15',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'mito.level',
             label=F,
             label.size=6,
             shuffle=TRUE)

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.25',
             split.by='timepoint',
             label=T,
             label.size=6)

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'mito.level',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='cell_calls',
             label=F,
             label.size=6)

p1 + p2 



p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'Phase',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.1',
             label=T,
             label.size=6)
p1 + p2 

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'scDblFinder.class',
             label=F,
             label.size=6,
             shuffle=TRUE )

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.25',
             label=T,
             label.size=6)
p1 + p2 


features = c('Adipoq','Krt17','Sdc4','Krt19','Cd8a','Cd4','Pecam1',
             'Col1a1','H2-Ab1')


FeaturePlot(integrated_seurat,
            features = 'Krt8',
            reduction='umap_harmony',
            min.cutoff='q10',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(integrated_seurat,
            features = 'nCount_RNA',
            reduction='umap_harmony',
            min.cutoff='q10',
            max.cutoff='q90',
)



FeaturePlot(integrated_seurat,
            features = 'nFeature_RNA',
            reduction='umap_harmony',
            min.cutoff='q10',
            max.cutoff='q90',
)

FeaturePlot(integrated_seurat,
            features = 'percent.mt',
            reduction='umap_harmony',
            min.cutoff='q10',
            max.cutoff='q90',
)


FeaturePlot(integrated_seurat,
            features = 'Adipoq',
            reduction='umap_harmony',
            min.cutoff='q10',
            max.cutoff='q90',
            cols=viridis(256)) &
  DarkTheme()

DimPlot(integrated_seurat,
        reduction = "umap_harmony",group.by='harmony_snn_res.0.25',
        # split.by = 'timepoint',
        label=T,
        label.size=6)


table(integrated_seurat$scDblFinder.class, integrated_seurat$sample)
# mn18_ctrl25 mn18_ctrl26 mn18_ctrl27 mn18_ctrl28 mn18_pb29 mn18_pb30 mn18_pb31 mn18_pb32 mn7_ctrl1 mn7_ctrl2
# doublet         404         260         307         328       279       337       371       268      1099       799
# singlet        3140        1963        2766        2916      3306      3287      3058      2821      6656      3799
# 
# mn7_ctrl3 mn7_ctrl4 mn7_pb5 mn7_pb6 mn7_pb7 mn7_pb8 wk10_ctrl33 wk10_ctrl34 wk10_ctrl35 wk10_ctrl36 wk10_pb37
# doublet       858       880    1100     816     855     819         807         684         608         636       540
# singlet      4334      4883    5523    4939    5478    5013        4168        4171        3766        4493      3640
# 
# wk10_pb38 wk10_pb39 wk10_pb40 wk3_ctrl21 wk3_ctrl22 wk3_ctrl23 wk3_ctrl24 wk3_ctrl25 wk3_ctrl26 wk3_pb27
# doublet       621       782       971        296        458        430        195        301        544      647
# singlet      3856      5090      5154       2583       2376       4659       2429       4144       5312     5781
# 
# wk3_pb28 wk3_pb29 wk3_pb30 wk3_pb31 wk3_pb32
# doublet      344      532      586      346      363
# singlet     2731     4543     5381     3262     2800


table(integrated_seurat$scDblFinder.class, integrated_seurat$batch)
# Experiment_1 Experiment_2
# doublet         5042        15429
# singlet        46001        98220

table(integrated_seurat$scDblFinder.class, integrated_seurat$timepoint)
# wk3  wk10   mn7  mn18
# doublet  5042  5649  7226  2554
# singlet 46001 34338 40625 23257

table(integrated_seurat$scDblFinder.class)
# doublet singlet 
# 20471  144221 


### check cluster silhouette index ----
library(cluster)
library(parallelDist)

wk10 = subset(integrated_seurat, subset=(timepoint=='wk10'))

# Calculate distance matrix 
dist_matrix <- parDist(Embeddings(wk10[["harmony"]])[,1:10],
                       method = "euclidean", threads = nthreads)

Idents(wk10) = 'harmony_snn_res.0.25'
# Calculate silhouette scores
sil_scores <- silhouette(as.integer(Idents(wk10)), dist_matrix)

summary(sil_scores)
sil_df <- data.frame(
  cluster = factor(sil_scores[,1]),  # Convert to factor for better plotting
  neighbor = sil_scores[,2],
  sil_width = sil_scores[,3]) %>%
  arrange(cluster)


ggplot(sil_df, aes(x = cluster, y = sil_width, fill = cluster)) + 
  geom_violin(draw_quantiles = 0.5, 
              scale = "width", 
              size = 1, 
              show.legend = FALSE) + 
  labs(x = "Louvain Cluster", 
       y = "Silhouette Score") + 
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_line(color = "grey80")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 1)





### calculate doublets again ----
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden.qs'),nthreads=nthreads)
integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_integrated_layers_by_batch.qs'),nthreads=nthreads)
gc()
integrated_seurat$scDblFinder.class = NULL
integrated_seurat$scDblFinder.score = NULL

sce = as.SingleCellExperiment(integrated_seurat, assay='RNA')
gc()

bp = MulticoreParam(nthreads, RNGseed=2024)
bpstart(bp)

# label doublets. This takes a while to run.
sce = scDblFinder(sce, samples = "batch", multiSampleMode = 'split',
                     # clusters = T,
                     clusters = 'harmony_snn_res.0.25',
                     dbr.sd = 1, # lets the dbr be set empirically, makes setting dbr irrelevant
                     BPPARAM=bp)
bpstop(bp)
register(SerialParam()) # close extra nodes/workers used
table(sce$scDblFinder.class)
# singlet doublet 
# 143109   21583 
table(sce$scDblFinder.class, sce$sample)
table(sce$scDblFinder.class, sce$timepoint)

# Store first-round doublet scores
sce$doublet_score_round1 <- sce$scDblFinder.score
sce$doublet_class_round1 <- sce$scDblFinder.class

qs_save(sce, paste0(data_dir_integrated, 'sce_scDblFinder.qs'), nthreads=nthreads)
sce = qs_read(paste0(data_dir_integrated, 'sce_scDblFinder.qs'), nthreads=nthreads)
gc()

# create singlets subset
sce_singlets <- sce[, sce$scDblFinder.class == "singlet"]
qs_save(sce_singlets, paste0(data_dir_integrated, 'sce_singlets.qs'), nthreads=nthreads)
gc()

#### Full reprocessing of normalization and clustering of singlets for second round of doublet detection (like She et al. 2025) ----
seu_singlets <- as.Seurat(sce_singlets, counts='counts', data=NULL)
seu_singlets_split = SplitObject(seu_singlets, split.by = "sample")
gc() # save RAM
options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
for(i in seq_along(seu_singlets_split)) { 
  cat('\n Beginning sample:', names(seu_singlets_split[i]), '\n\n')
  seu_singlets_split[i] = SCTransform(seu_singlets_split[[i]], assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                                vars.to.regress=c('cc.difference','percent.mt'),
                                return.only.var.genes=FALSE)
  }
gc() # free RAM

qs_save(seu_singlets_split, paste0(data_dir_integrated,'seu_singlets_split_singlets_integrated_by_batch_regressout_cc.diffandmt.qs'),nthreads=nthreads)
seu_singlets_split = qs_read(paste0(data_dir_integrated,'seu_singlets_split_singlets_integrated_by_batch_regressout_cc.diffandmt.qs'),nthreads=nthreads)
gc()

# seu_singlets_split = qs_read(paste0(data_dir_integrated,'seu_singlets_split_ct_regressout_cc.diffandmt.qs'),nthreads=nthreads)
# gc()
seu_singlets = merge(seu_singlets_split[[1]], seu_singlets_split[-1], merge.data=TRUE,merge.dr=FALSE)
gc() # free RAM
qs_save(seu_singlets, paste0(data_dir_integrated,'seu_singlets_integrated_by_batch_regressout_cc.diffandmt.qs'),nthreads=nthreads)
# seu_singlets = qs_read(paste0(data_dir_integrated,'seu_singlets_regressout_cc.diffandmt.qs'),nthreads=nthreads)
gc() #free up RAM
features = SelectIntegrationFeatures(seu_singlets_split, nfeatures=5000, fvf.nfeatures=5000)
qs_save(features,paste0(data_dir_integrated,'variable_features_seu_singlet_SCT.qs'),nthreads=nthreads)
VariableFeatures(seu_singlets) = features

sct_results = SCTResults(seu_singlets, slot = "feature.attributes")
qs_save(sct_results, paste0(data_dir_integrated,'sct_results_seu_singlets_SCT.qs'),nthreads=nthreads)
rm(sct_results)

VariableFeatures(seu_singlets[["SCT"]]) = rownames(seu_singlets[["SCT"]]@scale.data)

qs_save(seu_singlets, paste0(data_dir_integrated,'seu_singlets_regressout_cc.diffandmt_SCT.qs'),nthreads=nthreads)
# seu_singlets = qs_read(paste0(data_dir_integrated,'seu_singlets_regressout_cc.diffandmt_SCT.qs'),nthreads=nthreads)
rm(seu_singlets_split)
gc()


seu_singlets <- RunPCA(seu_singlets, assay='SCT',npcs=50,verbose=T)
seu_singlets <- RunHarmony(seu_singlets, group.by = "batch")
seu_singlets <- FindNeighbors(seu_singlets, dims = 1:50, reduction = "harmony",
                              graph.name=c('harmony_nn','harmony_snn'))

resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

seu_singlets <- FindClusters(seu_singlets, random.seed=1, resolution = resolutions,
                                  algorithm=4, graph.name='harmony_snn') # use leiden clustering

seu_singlets = RunUMAP(seu_singlets,
                            dims=1:50,
                            reduction='harmony',
                            reduction.name = 'umap_harmony')

DimPlot(seu_singlets,
        group.by = 'harmony_snn_res.0.25',
        label=T,
        label.size=6)

DimPlot(seu_singlets,
        group.by = 'harmony_snn_res.0.25',
        split.by = 'timepoint',
        label=T,
        label.size=6)


p1 = DimPlot(seu_singlets,
             group.by = 'harmony_snn_res.0.25',
             label=T,
             label.size=6)


p2 = DimPlot(seu_singlets,
        group.by = 'scDblFinder.class',
        label=T,
        label.size=6)
p1 + p2

FeaturePlot(seu_singlets,
            features = 'Ly6d',
            min.cutoff = 'q20',
            max.cutoff = 'q80',
            label=F,
            cols=viridis(256))+
  DarkTheme()


qs_save(seu_singlets, paste0(data_dir_integrated,'seu_singlets_regressout_cc.diffandmt_SCT_UMAP.qs'),nthreads=nthreads)
gc()

# Convert back to sce data type for round 2 of scDblFinder
sce_singlets = as.SingleCellExperiment(seu_singlets, assay='RNA')

bp = MulticoreParam(nthreads, RNGseed=2025) # Use different seed for second round
bpstart(bp)

# Run second round on singlets only
sce_singlets = scDblFinder(sce_singlets, samples = "batch", multiSampleMode = 'split',
                            # clusters = T,
                            clusters = 'harmony_snn_res.0.25',
                            dbr.sd = 1, # lets the dbr be set empirically, makes setting dbr irrelevant
                            BPPARAM=bp)

bpstop(bp)
register(SerialParam()) # close extra nodes/workers used

table(sce_singlets$scDblFinder.class)
# singlet doublet 
# 134854    8000 
table(sce_singlets$scDblFinder.class, sce_singlets$sample)
table(sce_singlets$scDblFinder.class, sce_singlets$timepoint)



# Transfer second-round results back to original object
sce$doublet_score_round2 <- NA
sce$doublet_class_round2 <- NA
sce$doublet_score_round2[match(colnames(sce_singlets), colnames(sce))] <- 
  sce_singlets$scDblFinder.score
sce$doublet_class_round2[match(colnames(sce_singlets), colnames(sce))] <- 
  sce_singlets$scDblFinder.class # returns numeric labels
sce$doublet_class_round2 <- factor(sce$doublet_class_round2, # set back to singlet and doublet labels
                                   levels = c(1, 2), 
                                   labels = c("singlet", "doublet"))
table(sce_singlets$scDblFinder.class)
table(sce$doublet_class_round2)

# A cell is classified as a doublet if identified in either round
sce$final_doublet_class <- "doublet"
singlet_cells <- sce$doublet_class_round1 == "singlet" & 
  (is.na(sce$doublet_class_round2) | sce$doublet_class_round2 == "singlet")
sce$final_doublet_class[singlet_cells] <- "singlet"

table(sce$final_doublet_class)
# doublet singlet 
# 29838  134854 

qs_save(sce, paste0(data_dir_integrated, 'sce_doublets_x2.qs'), nthreads=nthreads)
gc()

# doublet singlet 
# 20471  144221 
# # 
# qs_save(sce, paste0(data_dir_integrated,'sce_scDblFinder_run_second_time_on_integrated_data_using_random_doublets.qs'),nthreads=nthreads)
# gc()
# 
# 
# 
#### add doublets to seurat object ----
meta_sce <- sce@colData@listData %>% as.data.frame() %>%
  dplyr::select(starts_with('scDblFinder'), # 'scDblFinder.class' and 'scDblFinder.score'
                'doublet_class_round1', 'doublet_score_round1',
                'doublet_class_round2', 'doublet_score_round2',
                'final_doublet_class'
  )
rownames(meta_sce) = sce@colData@rownames

# integrated_seurat = AddMetaData(integrated_seurat, metadata = meta_sce %>% dplyr::select('scDblFinder.class','scDblFinder.score_2'))
# test = integrated_seurat
# test$scDblFinder.class = NULL
# test$scDblFinder.score = NULL
# integrated_seurat = AddMetaData(integrated_seurat, metadata = meta_sce %>% dplyr::select('scDblFinder.class', 'scDblFinder.score'))
integrated_seurat = AddMetaData(integrated_seurat, metadata = meta_sce)
# head(integrated_seurat@meta.data)
# table(integrated_seurat$scDblFinder.class)
rm(meta_sce, sce, sce_singlets, seu_singlets) # clean up RAM
qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_by_batch_scDblFinderx2.qs'),nthreads=nthreads)
gc()

table(integrated_seurat$final_doublet_class)
# doublet singlet 
# 29838  134854 

table(integrated_seurat$final_doublet_class, integrated_seurat$sample)
# wk10_pb37 wk10_pb38 wk10_pb39 wk10_pb40 wk3_ctrl21 wk3_ctrl22 wk3_ctrl23
# doublet       809       760      1343      1439        418        598       1117
# singlet      3371      3717      4529      4686       2461       2236       3972
# 
# wk3_ctrl24 wk3_ctrl25 wk3_ctrl26 wk3_pb27 wk3_pb28 wk3_pb29 wk3_pb30
# doublet        451        374        960      837      349      969      792
# singlet       2173       4071       4896     5591     2726     4106     5175
# 
# wk3_pb31 wk3_pb32
# doublet      512      311
# singlet     3096     2852
# > table(integrated_seurat$final_doublet_class, integrated_seurat$sample)
# 
# mn18_ctrl25 mn18_ctrl26 mn18_ctrl27 mn18_ctrl28 mn18_pb29 mn18_pb30
# doublet         246          97         281         235       225       200
# singlet        3298        2126        2792        3009      3360      3424
# 
# mn18_pb31 mn18_pb32 mn7_ctrl1 mn7_ctrl2 mn7_ctrl3 mn7_ctrl4 mn7_pb5
# doublet       391       337      2038      1106      1195      1889    1165
# singlet      3038      2752      5717      3492      3997      3874    5458
# 
# mn7_pb6 mn7_pb7 mn7_pb8 wk10_ctrl33 wk10_ctrl34 wk10_ctrl35 wk10_ctrl36
# doublet    1222    1775    1319         894         997        1009        1178
# singlet    4533    4558    4513        4081        3858        3365        3951
# 
# wk10_pb37 wk10_pb38 wk10_pb39 wk10_pb40 wk3_ctrl21 wk3_ctrl22 wk3_ctrl23
# doublet       809       760      1343      1439        418        598       1117
# singlet      3371      3717      4529      4686       2461       2236       3972
# 
# wk3_ctrl24 wk3_ctrl25 wk3_ctrl26 wk3_pb27 wk3_pb28 wk3_pb29 wk3_pb30
# doublet        451        374        960      837      349      969      792
# singlet       2173       4071       4896     5591     2726     4106     5175
# 
# wk3_pb31 wk3_pb32
# doublet      512      311
# singlet     3096     2852

table(integrated_seurat$scDblFinder.class, integrated_seurat$sample)
# mn18_ctrl25 mn18_ctrl26 mn18_ctrl27 mn18_ctrl28 mn18_pb29
# singlet        3541        2223        3037        3239      3557
# doublet           3           0          36           5        28
# 
# mn18_pb30 mn18_pb31 mn18_pb32 mn7_ctrl1 mn7_ctrl2 mn7_ctrl3
# singlet      3617      3362      3062      5478      3758      4460
# doublet         7        67        27      2277       840       732
# 
# mn7_ctrl4 mn7_pb5 mn7_pb6 mn7_pb7 mn7_pb8 wk10_ctrl33 wk10_ctrl34
# singlet      3818    6034    5224    4934    5181        4702        4512
# doublet      1945     589     531    1399     651         273         343
# 
# wk10_ctrl35 wk10_ctrl36 wk10_pb37 wk10_pb38 wk10_pb39 wk10_pb40
# singlet        3898        4771      3686      4160      5302      5419
# doublet         476         358       494       317       570       706
# 
# wk3_ctrl21 wk3_ctrl22 wk3_ctrl23 wk3_ctrl24 wk3_ctrl25 wk3_ctrl26
# singlet       2844       2718       4568       2543       4418       5443
# doublet         35        116        521         81         27        413
# 
# wk3_pb27 wk3_pb28 wk3_pb29 wk3_pb30 wk3_pb31 wk3_pb32
# singlet     6156     3044     5002     5846     3549     3159
# doublet      272       31       73      121       59        4

table(integrated_seurat$scDblFinder.class, integrated_seurat$timepoint)
# wk3  wk10   mn7  mn18
# singlet 45267 34058 39041 24488
# doublet  5776  5929  8810  1323

p1 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by = 'final_doublet_class',
             label=F,
             label.size=6,
             shuffle=TRUE)

p2 = DimPlot(integrated_seurat,
             reduction = "umap_harmony",group.by='harmony_snn_res.0.25',
             label=T,
             label.size=6)

p1 + p2 

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_doublets2x.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_doublets2x.qs'),nthreads=nthreads)
gc()

# plot doublets ----
plot_doublets = function(seurat_object,fig_dir=fig_dir_integrated, doubletmetadata = 'scDblFinder.class') {
  metrics = c("nCount_RNA", "nFeature_RNA", 'complexity',"percent.mt", "percent.hb")
  # Check how doublets and singlets differ in QC measures per sample
  cat('\nPlotting doublet metrics by sample\n')
  VlnPlot(seurat_object, group.by = 'sample', split.by = doubletmetadata,
          features = metrics, 
          ncol = 3, pt.size = 0) + theme(legend.position = 'right')
  ggsave(paste0(fig_dir,'/doublets_by_sample_',doubletmetadata,'.png'), dpi=320, width=20, height= 20)
  
  # Check how doublets and singlets differ in QC measures by condition
  cat('\nPlotting doublet metrics by condition\n')
  VlnPlot(seurat_object, group.by = 'condition', split.by = doubletmetadata,
          features = metrics, 
          ncol = 3, pt.size = 0) + theme(legend.position = 'right')
  ggsave(paste0(fig_dir,'/doublets_by_condition_',doubletmetadata,'.png'), dpi=320, width=20, height= 20)
}

plot_doublets(integrated_seurat, doubletmetadata = 'doublet_class_round1')
plot_doublets(integrated_seurat, doubletmetadata = 'final_doublet_class')

## doublets QC ----
umap_doublets_qc = function(seurat_umap, reduction = 'umap_harmony',resolution, directory) {
  
  if (!dir.exists(paste0(directory, resolution,'/'))) {dir.create(paste0(directory, resolution,'/'),
                                                                  recursive=TRUE)}
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  reduction_name = toupper(reduction)
  
  # set umap resolution
  Idents(seurat_umap) = resolution
  seurat_umap@meta.data$resolution = Idents(seurat_umap)
  
  
  
  # plot side by side UMAP of clusters and doublet classification
  p1 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = resolution,
               label=T,
               repel=T,
               label.size=6,
               shuffle=TRUE)
  
  p2 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = 'final_doublet_class',
               label=F,
               label.size=6,
               shuffle=TRUE)
  p1 | p2
  ggsave(filename=paste0(directory,resolution,'/',reduction_name,'_doublets.png'), dpi=320, width=20, height=10)
  
  p1 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = resolution,
               label=T,
               repel=T,
               label.size=6,
               shuffle=TRUE)
  
  p2 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = 'doublet_class_round1',
               label=F,
               label.size=6,
               shuffle=TRUE,
               cols=c( "#00BFC4","#F8766D"))
  
  p3 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = 'doublet_class_round2',
               label=F,
               label.size=6,
               shuffle=TRUE,
               cols=c( "#00BFC4","#F8766D"))
  
  p4 = DimPlot(integrated_seurat,
               reduction = reduction,
               group.by = 'final_doublet_class',
               label=F,
               label.size=6,
               shuffle=TRUE)
  p1 + p2 + p3 + p4 + plot_layout(ncol=2)
  ggsave(filename=paste0(directory,resolution,'/',reduction_name,'_doublets_round1_and_round2.png'), dpi=320, width=20, height=20)
  
  
  # Extract identity and sample information from seurat object to determine 
  # the number of cells per cluster per sample
  n_cells = FetchData(seurat_umap, 
                      vars = c("ident", "final_doublet_class")) %>%
    dplyr::count(ident, final_doublet_class)
  
  # barplot of singlets vs doublets by cluster
  ggplot(n_cells, aes(x=ident, y=n, fill=final_doublet_class))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_scDblFinder.png'), dpi=320, width=20, height= 10)
  
  # Barplot of proportion of doublets by cluster
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=final_doublet_class), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5),
          legend.key.size = unit(1,'cm'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20))
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_scDblFinder.png'), dpi=320, width=20, height= 10)
  
}

### plot doublets on UMAP ----
# umap_doublets_qc(integrated_seurat, resolution='SCT_snn_res.0.8', directory=fig_dir_integrated)
# umap_doublets_qc(integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.2', directory=fig_dir_integrated)
umap_doublets_qc(integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.1', directory=fig_dir_integrated)
umap_doublets_qc(integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.25', directory=fig_dir_integrated)


### subset doublets out of dataset once more----

# tabulate doublets by sample
table(integrated_seurat$final_doublet_class, integrated_seurat$sample)


# tabulate total doublets
table(integrated_seurat$final_doublet_class)

integrated_seurat = subset(integrated_seurat, subset=(final_doublet_class == 'singlet'))

table(integrated_seurat$final_doublet_class)
# singlet 
# 134854

# qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_leiden_UMAP_doublets_removed_batch_integrated.qs'),nthreads=nthreads)
gc()

## re-do PCA and harmony integration with doublets removed x2, clustering on singlets, ----

split_seurat = SplitObject(integrated_seurat, split.by = "sample")

gc() # save RAM

# adjust object size limits of variables/objects in R
options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
for(i in seq_along(split_seurat)) { 
  cat('\n Beginning sample:', names(split_seurat[i]), '\n\n')
  split_seurat[i] = SCTransform(split_seurat[[i]], assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                                vars.to.regress=c('cc.difference','percent.mt'),
                                return.only.var.genes=FALSE)
  
}
gc() # free RAM

qs_save(split_seurat, paste0(data_dir_integrated,'split_seurat_doublets_2x.qs'),nthreads=nthreads)
gc()
integrated_seurat = merge(split_seurat[[1]], split_seurat[-1], merge.data=TRUE,merge.dr=FALSE)
integrated_seurat = JoinLayers(integrated_seurat, assay='RNA')
gc() # free RAM

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_doublets_2x.qs'),nthreads=nthreads)

# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_regressout_cc.diffandmt.qs'),nthreads=nthreads)
gc() #free up RAM

features = SelectIntegrationFeatures(split_seurat, nfeatures=5000, fvf.nfeatures=5000)
qs_save(features,paste0(data_dir_integrated,'variable_features_doublets_2x.qs'),nthreads=nthreads)
VariableFeatures(integrated_seurat) = features

sct_results = SCTResults(integrated_seurat, slot = "feature.attributes")
qs_save(sct_results, paste0(data_dir_integrated,'sct_results_doublets_2x.qs'),nthreads=nthreads)
rm(sct_results, split_seurat)
gc()

# Variable features lost/unset after merging transformed objects. Need to set again.
# should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(integrated_seurat[["SCT"]]) = rownames(integrated_seurat[["SCT"]]@scale.data)

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat_doublets_2x.qs'),nthreads=nthreads)
gc()

integrated_seurat = RunPCA(integrated_seurat, assay='SCT',npcs=50,verbose=T)

integrated_seurat = RunHarmony(integrated_seurat, c('batch'))

# gc()
# 
# integrated_seurat = IntegrateLayers(
#   object = integrated_seurat,
#   method = HarmonyIntegration,
#   orig.reduction = 'pca',
#   new.reduction = 'harmony',
#   assay='SCT',
#   normalization.method = "SCT",
#   verbose = T
# )

# integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "integrated.dr",
#                                  graph.name = c('idr_nn','idr_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                 graph.name = 'idr_snn')
# integrated_seurat = RunUMAP(integrated_seurat,
#                           dims=1:50,
#                           reduction='integrated.dr',
#                           reduction.name = 'umap_idr')

integrated_seurat <- FindNeighbors(integrated_seurat, dims = 1:50, reduction = "harmony",
                                   graph.name=c('harmony_nn','harmony_snn'))
# integrated_seurat <- FindClusters(integrated_seurat, resolution = resolutions,
#                                   graph.name='harmony_snn')
resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

integrated_seurat <- FindClusters(integrated_seurat, random.seed=1, resolution = resolutions,
                                  algorithm=4, graph.name='harmony_snn') # use leiden clustering
integrated_seurat = RunUMAP(integrated_seurat,
                            dims=1:50,
                            reduction='harmony',
                            reduction.name = 'umap_harmony')

gc() # free RAM

# reintroduce order to ordinal variables that were lost during doublet classification ----
integrated_seurat$timepoint = factor(integrated_seurat$timepoint,
                                     levels=c('wk3','wk10','mn7','mn18'),
                                     ordered = TRUE)
integrated_seurat$batch = factor(integrated_seurat$batch,
                                     levels=c('Experiment_1','Experiment_2'),
                                     ordered = FALSE)
integrated_seurat$mito.level = factor(integrated_seurat$mito.level,
                                      levels=c("1st","2nd","3rd", "4th"),
                                      ordered = TRUE)


qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden_doublets_x2.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden_doublets_x2.qs'),nthreads=nthreads)
gc()

DimPlot(integrated_seurat,
        reduction = "umap_harmony",group.by='harmony_snn_res.0.25',
        label=T,
        label.size=6)




# Visualizations of clustering----
## visualize clustering results with clustree ----
library(clustree)

# clustree(integrated_seurat, prefix='SCT_snn_res.')
# ggsave(paste0(fig_dir_integrated,'clustree.png'), dpi=320, width=20, height= 10)

clustree(integrated_seurat, prefix='harmony_snn_res.')
ggsave(paste0(fig_dir_integrated,'clustree_harmony_integrated.png'), dpi=320, width=20, height= 10)

# clustree(integrated_seurat, prefix = 'SCT_snn_res.',
#          node_colour = 'Pgr', node_colour_aggr = 'mean')



## function to plot initial umaps ----
plot_umaps = function (seurat_umap, reduction='umap',resolution, directory) {
  
  if (!dir.exists(paste0(directory, resolution))) { dir.create(paste0(directory, resolution,'/'),
                                                              recursive=TRUE) }
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # Assign umap resolution
  Idents(seurat_umap) = resolution
  
  reduction_name = toupper(reduction)
  
  # Plot the UMAP of non-integrated data
  DimPlot(seurat_umap,
          reduction = reduction,
          label=TRUE,
          repel=TRUE,
          label.size=6)
  ggsave(paste0(directory,resolution,'/',reduction_name,'.png'), dpi=320, width=10, height=10)
  
  
  # plot UMAP by sample and condition
  p1 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='sample',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='condition',
               shuffle=T)
  
  png(filename=paste0(directory,resolution,'/',reduction_name,'_by_sample_and_condition.png'), units='in', width=20, height=10, res=320)
  print(plot_grid(p1, p2))
  dev.off()
  
  # plot UMAP by sample and cell calls in filtered cellranger and cellbender
  p1 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='sample',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='cell_calls',
               shuffle=T)
  
  png(filename=paste0(directory,resolution,'/',reduction_name,'_by_sample_and_cell_calls_CELLRANGER_vs_CELLBENDER.png'), units='in', width=20, height=10, res=320)
  print(plot_grid(p1, p2))
  dev.off()
  
  
  # plot UMAP by Phase and mitochondrial expression quartile
  p1 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='Phase',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='mito.level',
               shuffle=T,
               cols = hue_pal()(length(unique(seurat_umap$mito.level))))
  
  p1 | p2
  ggsave(filename=paste0(directory,resolution,'/',reduction_name,'_by_cellcycle_and_mito.quartile.png'), dpi=320, width=20, height=10)

  # plot UMAP by cell calls and mitochondrial expression quartile
  p1 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='cell_calls',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = reduction,
               group.by='mito.level',
               shuffle=T,
               cols = hue_pal()(length(unique(seurat_umap$mito.level))))
  
  p1 | p2
  ggsave(filename=paste0(directory,resolution,'/',reduction_name,'_by_cell_calls_and_mito.quartile.png'), dpi=320, width=20, height=10)
  
}

### run plot_umaps ----
# plot_umaps(seurat_umap=integrated_seurat, resolution='SCT_snn_res.0.8', directory=fig_dir_integrated)
plot_umaps(seurat_umap=integrated_seurat, reduction='umap_harmony', resolution='harmony_snn_res.0.1', directory=fig_dir_integrated)
plot_umaps(seurat_umap=integrated_seurat, reduction='umap_harmony', resolution='harmony_snn_res.0.25', directory=fig_dir_integrated)


## UMAP clustering QC ----
metrics =  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt", "percent.hb") # metrics to overlay on umap

### create function to plot all QC plots and umap overlays ----
umap_clustering_qc = function(seurat_umap, reduction = 'umap', resolution, directory) {
  
  if (!dir.exists(paste0(directory, resolution))) { dir.create(paste0(directory, resolution,'/'),
                                                               recursive=TRUE) }
  DefaultAssay(seurat_umap) = 'SCT'
  
  reduction_name = toupper(reduction)
  
  # set umap resolution
  Idents(seurat_umap) = resolution
  seurat_umap@meta.data$resolution = Idents(seurat_umap)
  
  # Explore whether clusters segregate by cell cycle phase
  DimPlot(seurat_umap,
          reduction = reduction,
          label = TRUE, 
          split.by = "Phase")  + NoLegend()
  ggsave(paste0(directory,resolution,'/split_by_cell_cycle_',reduction_name,'.png'), dpi=320, width=30, height=20)
  
  # Metrics to plot on UMAP
  FeaturePlot(seurat_umap, 
              reduction = reduction, 
              features = metrics,
              pt.size = 0.4, 
              order = TRUE,
              min.cutoff = 'q10',
              label = TRUE)
  ggsave(paste0(directory,resolution,'/',reduction_name,'_counts_cellcycle_mt_metrics.png'), dpi=320, width=20, height=20)
  
  # Boxplot of nUMIs per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nCount_RNA, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nUMIs_per_cluster.png'), dpi=320, width=20, height= 10)
  
  # Boxplot of nGenes per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nFeature_RNA, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nGenes_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Boxplot of percent.mt per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=percent.mt, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/percent.mt_per_cluster.png'), dpi=320, width=20, height= 10)
  
  # Boxplot of percent.hb per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=percent.hb, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/percent.hb_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Boxplot of S.Score per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=S.Score, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/S.score_per_cluster.png'), dpi=320, width=20, height= 10)
  

  # Boxplot of G2M.Score per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=G2M.Score, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/G2M.score_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Extract identity and sample information from seurat object to determine 
  # the number of cells per cluster per sample
  n_cells = FetchData(seurat_umap, 
                       vars = c("ident", "sample")) %>%
    dplyr::count(ident, sample)
  
  # Barplot of number of cells per cluster by sample
  ggplot(n_cells, aes(x=ident, y=n, fill=sample))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
    
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_sample.png'), dpi=320, width=20, height= 10)
  
  
  # Extract identity and condition information from seurat object to determine 
  # the number of cells per cluster per condition
  n_cells = FetchData(seurat_umap, 
                      vars = c("ident", "condition")) %>%
    dplyr::count(ident, condition)
  
  # Barplot of number of cells per cluster by condition
  ggplot(n_cells, aes(x=ident, y=n, fill=condition))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_condition.png'), dpi=320, width=20, height= 10)
  
  # Barplot of proportion of cells in each cluster by sample
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=sample), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5),
          legend.key.size = unit(1,'cm'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20))
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_sample.png'), dpi=320, width=20, height= 10)
  
  
  # Barplot of proportion of cells in each cluster by condition
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=condition), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5),
          legend.key.size = unit(1,'cm'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20))
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_condition.png'), dpi=320, width=20, height= 10)

  
  # histogram by cluster for nCount_RNA
  ggplot(seurat_umap@meta.data, aes(x = nCount_RNA, fill = resolution)) +
    geom_histogram(position = "stack", bins = 100) +
    theme_classic2() +
    labs(x = "Total UMI counts", 
         y = "Frequency",
         title = paste0("Distribution of UMI counts by ", resolution)) +
    theme(text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(n.breaks = 20)
  ggsave(paste0(directory,resolution,"/histogram_nCount_RNA_by",resolution,'.png'), width = 20, height = 10)
  
  
  # histogram by cluster for nFeature_RNA
  ggplot(seurat_umap@meta.data, aes(x = nFeature_RNA, fill = resolution)) +
    geom_histogram(position = "stack", bins = 100) +
    theme_classic2() +
    labs(x = "Number of genes detected", 
         y = "Frequency",
         title = paste0("Distribution of genes detected by ", resolution))+
           theme(text = element_text(size = 12),
                 plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(n.breaks = 20)
  ggsave(paste0(directory,resolution,"/histogram_nFeature_RNA_by",resolution,'.png'), width = 20, height = 10)
  
  # histogram by cluster for percent.mt
  ggplot(seurat_umap@meta.data, aes(x = percent.mt, fill = resolution)) +
    geom_histogram(position = "stack", bins = 100) +
    theme_classic2() +
    labs(x = "Mitochondrial percentage", 
         y = "Frequency",
         title = paste0("Distribution of mitochondrial percentage by ", resolution))+
    theme(text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(n.breaks = 20)
  ggsave(paste0(directory,resolution,"/histogram_percent.mt_by",resolution,'.png'), width = 20, height = 10)
  
  ggplot(seurat_umap@meta.data, aes(x = complexity, fill = resolution)) +
    geom_histogram(position = "stack", bins = 100) +
    theme_classic2() +
    labs(x = "log10(genes) / log10(UMIs)", 
         y = "Frequency",
         title = paste0("Distribution of cell complexity by ", resolution))+
    theme(text = element_text(size = 12),
          plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(n.breaks = 20)
  ggsave(paste0(directory,resolution,"/histogram_complexity_by_",resolution,'.png'), width = 20, height = 10)

}


#### run umap_clustering_qc ----
# umap_clustering_qc(seurat_umap=integrated_seurat, resolution='SCT_snn_res.0.8', directory=fig_dir_integrated)
# umap_clustering_qc(seurat_umap=integrated_seurat, resolution='SCT_snn_res.0.8', directory=fig_dir_integrated)
umap_clustering_qc(seurat_umap=integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.1', directory=fig_dir_integrated)
umap_clustering_qc(seurat_umap=integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.25', directory=fig_dir_integrated)






## Which PCs drive different clusters? ----

# determine max number of PCs present in seurat object
max_pcs = length(grep("^PC_", colnames(integrated_seurat@reductions$pca@cell.embeddings)))

## overlay PCs on umap ----
# takes a while to run so will use parallel processing
# library(future.apply)
# library(parallel)
# 
# options(future.globals.maxSize = 8 * 1024^3)  # Set limit to 8 GB

## function to overlay PCs on umap
# takes a while to run so will use parallel processing with 3 workers
# pcs_umap_overlay = function(seurat_umap, reduction='umap', resolution, directory,max_pcs=50) {
#   
#   if (!dir.exists(paste0(directory, resolution,'/'))) {dir.create(paste0(directory, resolution,'/'),
#                                                                   recursive=TRUE)}
#   
#   DefaultAssay(seurat_umap) = 'SCT'
#   Idents(seurat_umap) = resolution
#   
#   reduction_name = toupper(reduction)
#   
#   # Extract all PC data at once
#   pc_data = FetchData(seurat_umap, vars = c("umap_1", "umap_2", paste0("PC_", 1:max_pcs), "ident"))
#   
#   # Adding cluster label data
#   umap_label = pc_data %>%
#     group_by(ident) %>%
#     dplyr::summarise(umap_1=mean(umap_1), umap_2=mean(umap_2))
#   
#   # Reshape to long format for faceting
#   pc_data_long = pc_data %>%
#     pivot_longer(cols = starts_with("PC_"), 
#                  names_to = "PC", 
#                  values_to = "PC_value")
#   
#   # Split into chunks of 12 PCs and plot
#   pc_chunks = split(unique(pc_data_long$PC), ceiling(seq_along(unique(pc_data_long$PC))/12))
#   
#   for(i in seq_along(pc_chunks)) {
#     pc_subset = pc_data_long %>% filter(PC %in% pc_chunks[[i]])
#     first_pc = pc_chunks[[i]][1]
#     last_pc = pc_chunks[[i]][length(pc_chunks[[i]])]
#     
#     p = ggplot(pc_subset, aes(umap_1, umap_2, color = PC_value)) +
#       geom_point(alpha = 0.7, size = 0.5) +
#       MetBrewer::scale_color_met_c('Hokusai2', guide='none') +
#       facet_wrap(~PC, ncol = 4) +
#       geom_text(data = umap_label, aes(label = ident), color = "black", size = 3) +
#       theme_classic() +
#       theme(strip.text = element_text(size = 14, face = "bold"))
#     
#     ggsave(paste0(directory, resolution, '/', reduction_name, '_by_PC_components_', 
#                   first_pc, '-', last_pc, '.png'), 
#            plot = p, dpi = 320, width = 20, height = 15)
#   }
# }

pcs_umap_overlay_optimized = function(seurat_umap, reduction='umap', resolution, directory, max_pcs=50) {
  if (!dir.exists(paste0(directory, resolution,'/'))) {dir.create(paste0(directory, resolution,'/'),
                                                                  recursive=TRUE)}
  
  DefaultAssay(seurat_umap) = 'SCT'
  Idents(seurat_umap) = resolution
  
  reduction_name = toupper(reduction)
  
  # Extract UMAP coordinates once
  umap_coords = as.data.frame(Embeddings(seurat_umap, reduction = reduction))
  colnames(umap_coords) = c("umap_1", "umap_2")
  
  # Add cluster IDs
  umap_coords$cluster = as.character(Idents(seurat_umap))
  
  # Pre-calculate cluster centers once - ensure we get a usable data frame
  cluster_centers_df = data.frame()
  for(clust in unique(umap_coords$cluster)) {
    subset_data = umap_coords[umap_coords$cluster == clust, ]
    x_mean = mean(subset_data$umap_1)
    y_mean = mean(subset_data$umap_2)
    cluster_centers_df = rbind(cluster_centers_df, 
                               data.frame(cluster = clust, 
                                          x = x_mean, 
                                          y = y_mean))
  }
  
  # Process PCs in chunks of 12 to limit memory usage
  pc_chunks = split(1:max_pcs, ceiling(seq_along(1:max_pcs)/12))
  
  # Process each chunk sequentially
  for(chunk_idx in seq_along(pc_chunks)) {
    pc_indices = pc_chunks[[chunk_idx]]
    first_pc = pc_indices[1]
    last_pc = pc_indices[length(pc_indices)]
    
    # Extract only the PCs needed for this chunk
    pc_data = as.data.frame(seurat_umap[["pca"]]@cell.embeddings[, paste0("PC_", pc_indices)])
    
    # Create a combined data frame for plotting
    plot_df = cbind(umap_coords, pc_data)
    
    # Create all plots for this chunk at once
    plots_list = list()
    
    for(i in seq_along(pc_indices)) {
      pc_num = pc_indices[i]
      pc_col = paste0("PC_", pc_num)
      
      # Calculate data-specific color limits for better visualization
      pc_min = min(plot_df[, pc_col])
      pc_max = max(plot_df[, pc_col])
      pc_abs_max = max(abs(c(pc_min, pc_max)))
      
      # Transform values to spread out extremes (keeping the sign)
      pc_transformed = paste0(pc_col, "_transformed")
      plot_df[[pc_transformed]] = sign(plot_df[[pc_col]]) * sqrt(abs(plot_df[[pc_col]]))
      
      # Calculate alpha values based on the absolute value of PC loading
      # Scale between 0.3 and 0.9 for better visibility
      alpha_col = paste0(pc_col, "_alpha")
      plot_df[[alpha_col]] = scales::rescale(abs(plot_df[[pc_col]]), to = c(0.3, 0.9))
      
      p = ggplot(plot_df, aes(umap_1, umap_2)) +
        # Replace geom_raster with geom_point
        geom_point(
          aes(
            color = !!sym(pc_transformed),
            alpha = !!sym(alpha_col)
          ), 
          size = 0.8  # Adjust point size as needed
        ) +
        MetBrewer::scale_color_met_c('Hokusai2') +
        guides(color = "none") +
        # scale_color_gradient2(
        #   low = "blue", 
        #   mid = "white", 
        #   high = "red", 
        #   midpoint = 0, 
        #   guide = "none"
        # ) +
        # Turn off alpha legend
        scale_alpha(guide = "none") +
        # Use the explicit data frame to avoid issues
        geom_text(
          data = cluster_centers_df,
          mapping = aes(x = x, y = y, label = cluster),
          size = 3.5, 
          color = "black",
          fontface = "bold"
        ) +
        ggtitle(pc_col) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid = element_blank()
        )
      
      plots_list[[i]] = p
    }
    
    # Calculate optimal number of columns
    n_plots = length(plots_list)
    if (n_plots <= 3) {
      n_cols = min(n_plots, 2)  # 1-3 plots: use 1-2 cols
    } else if (n_plots <= 8) {
      n_cols = 3  # 4-8 plots: use 3 cols
    } else {
      n_cols = 4  # 9+ plots: use 4 cols
    }
    
    # Combine all plots for this chunk with dynamic column count
    combined_plot = cowplot::plot_grid(plotlist = plots_list, ncol = n_cols)

    # Save the combined plot
    filename = paste0(directory, resolution, '/', reduction_name, '_by_PC_components_', first_pc, '-', last_pc, '.png')
    
    # Additionally, adjust the figure width based on the number of columns
    fig_width = 5 * n_cols  # Base width per column
    ggsave(filename, combined_plot, width = fig_width, height = min(10, 3 * ceiling(n_plots/n_cols)), dpi = 180)    
    # Clear memory
    rm(plot_df, plots_list, combined_plot)
    gc()
  }
}

## plot PCs on umap for resolution SCT_snn_res.0.8 ----
# takes a while to run
# pcs_umap_overlay(seurat_umap=integrated_seurat, resolution='SCT_snn_res.0.8', directory=fig_dir_integrated)
pcs_umap_overlay_optimized(seurat_umap=integrated_seurat, reduction='umap_harmony',resolution='harmony_snn_res.0.25', directory=fig_dir_integrated)


# PC_7 looks like it's driving cell cycle genes
# PC_13 looks like it's driving myoepithelial cells
# PC_14 marcophage mb
# PC_17 mystery cluster by endothelial cells
# PC_43 and PC_46 maybe cells transitioning between fibroblast and epithelial 


#### exploring cell marker genes ----
# library(readxl)
# 
# # read in cell markers from excel file
# cell_markers = read_excel('/nfs/turbo/sph-colacino/aguilada/scRNAseq/8651-JC-mammary/3_weeks/annotation/aggregate_sources_celltype_markers.xlsx')
# cell_markers = lapply(cell_markers, na.omit) # convert to list of lists, omit NAs
# 
# # set RNA counts slot to be default and normalize for visualizing genes
# # SCTransform normalization was performed only on the 3000 most variable genes, 
# # so many of our genes of interest may not be present in the SCT data. Hence,
# # we switch over to the original RNA count to explore cell type marker genes.
# DefaultAssay(integrated_seurat) = 'RNA'
# # integrated_seurat = NormalizeData(integrated_seurat, verbose = FALSE)
# # DefaultAssay(integrated_seurat) = 'SCT'
# # qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat.qs'),nthreads=nthreads)
# 
# # Want to only consider cell markers that are actually present in the dataset
# cell_markers = lapply(cell_markers, function(x){
#   x = unique(x[x %in% rownames(integrated_seurat)])
# })
# 
# # remove any cell types with zero markers in the dataset
# cell_markers = Filter(length, cell_markers)
# # Erythrocyte cell type removed since no gene markers present in this dataset
# 
# # process and save the gene markers actually present in the data into a new file
# max_length = max(sapply(cell_markers, length)) # max length lists of gene markers
# 
# # Pad shorter lists with NA to match max length so all gene lists have same length
# padded_markers = lapply(cell_markers, function(x) {
#   length(x) = max_length
#   return(x)
# })
# 
# # save gene markers actually present in the dataset to csv
# write.csv(as.data.frame(padded_markers), 
#           file='annotation/aggregate_sources_celltype_markers_actually_in_data.csv',
#           row.names=FALSE,
#           na='')


# # function to plot all gene markers on umap and dotplots for each cell type in batches of 16 genes ----
# # takes a while to run
# plot_gene_markers = function(seurat_umap, resolution, directory) {
#   
#   DefaultAssay(seurat_umap) = 'RNA'
#   
#   # set resolution
#   Idents(seurat_umap) = resolution
#   
#   lapply(seq_along(cell_markers), function(i) {
#     gene_chunks = split(cell_markers[[i]], ceiling(seq_along(cell_markers[[i]])/16)) # plot 16 gene markers at a time
#     lapply(seq_along(gene_chunks), function(j) {
#       FeaturePlot(seurat_umap, # overlay gene marker expression on umap
#                   reduction='umap',
#                   features=gene_chunks[[j]],
#                   order=T,
#                   min.cutoff='q10',
#                   label=T)+
#         plot_annotation(
#           title=paste0(names(cell_markers[i]),'-',j),
#           theme=theme(plot.title = element_text(hjust = 0.5,size=24,face='bold')))
#       ggsave(paste0(directory,resolution,'/annotation/',names(cell_markers[i]),'-',j,'.png'), 
#                     dpi=320, width=20, height=20)
#       
#       png(filename=paste0(directory,resolution,'/annotation/',names(cell_markers[i]),'-',j,'_dotplot.png'), 
#           units='in', width=20, height=10, res=320)
#       print(
#         DotPlot(seurat_umap, features=gene_chunks[[j]])+ # construct dotplot of gene markers for each cell type
#           RotatedAxis()+
#           ggtitle(paste0(names(cell_markers[i]),'-',j))+
#           theme(plot.title = element_text(hjust = 0.5, size=24, face='bold'))
#       )
#       graphics.off()
#     })
#   })
# }
#   
  
# ## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.1. ----
# # takes a while to run
# plot_gene_markers(seurat_umap=integrated_seurat, resolution='SCT_snn_res.0.1', directory=fig_dir_integrated)
# 
# # looks like we could benefit from increasing resolution
# # resolutions available. For reference
# resolutions = sprintf('SCT_snn_res.%.1f',c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 2.0))
# 
# # Assign identity of clusters to desired resolution
# Idents(integrated_seurat) = 'SCT_snn_res.1.4'
# 
# DimPlot(integrated_seurat,
#         reduction = "umap",
#         label=TRUE,
#         label.size=6)



DimPlot(integrated_seurat,
        group.by='harmony_snn_res.0.25',
        reduction='umap_harmony',
        label=T,
        repel=T,
        shuffle=T,
        )

DimPlot(integrated_seurat,
        group.by='Phase',
        reduction='umap_harmony',
        label=F,
        shuffle=T,
)

#### Find cluster markers ----


# join layers to be able to run FindAllMarkers()
# integrated_seurat = JoinLayers(integrated_seurat, assays='RNA')

# find_markers = function(seurat_obj, ident, assay='SCT') {}

Idents(integrated_seurat) = 'harmony_snn_res.0.25' # set cluster #s as ident

DefaultAssay(integrated_seurat) = 'SCT'

# Find marker genes ----

## Find all positive markers ----
integrated_seurat = PrepSCTFindMarkers(integrated_seurat)

qs_save(integrated_seurat, paste0(data_dir_integrated,'harmony_integrated_seurat_singlets_only_SCTprep.qs'))
integrated_seurat = qs_read(paste0(data_dir_integrated,'harmony_integrated_seurat_singlets_only_SCTprep.qs'), nthreads=nthreads)
gc()

all.pos.markers = FindAllMarkers(integrated_seurat,
                                 only.pos=TRUE,
                                 logfc.threshold = 0.25,
                                 min.pct=0.1,
                                 min.diff.pct = -Inf) 

# Combine markers with gene descriptions 
annotations = read.csv('/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/3_weeks/annotation/annotationhub_mice_genes.csv')
ann.pos.markers = left_join(x = all.pos.markers, 
                          y = annotations[, c("gene_name", "description")],
                          by = c("gene" = "gene_name")) %>%
  unique()

# add pct.diff col
ann.pos.markers$pct.diff = ann.pos.markers$pct.1 - ann.pos.markers$pct.2

# Rearrange the columns to be more intuitive
ann.pos.markers = ann.pos.markers[ , c(6, 7, 2:4, 9, 1, 5,8)]

# Order the rows by log2FC values
ann.pos.markers = ann.pos.markers %>%
  dplyr::arrange(cluster, desc(avg_log2FC))

# save markers
write.csv(ann.pos.markers,paste0(annot_dir_integrated,'annotated_pos_markers_res.0.25_df.csv'),row.names=FALSE)


# Extract top 10 markers per cluster
top20 = ann.pos.markers %>% 
  group_by(cluster) %>%
  filter(p_val_adj < 0.05, pct.diff > 0.5) %>%
  slice_max(order_by=tibble(avg_log2FC),n=20)
write.csv(top20,paste0(annot_dir_integrated,'top20_markers_res.0.25_df.csv'),row.names=FALSE)



## Find all positive conserved markers between conditions ----
# Create function to get conserved markers for any given cluster
get_conserved = function(cluster, assay = 'SCT', group.by='condition'){
  FindConservedMarkers(integrated_seurat, # need to specify exact seurat object here
                       ident.1 = cluster,
                       assay=assay,
                       grouping.var = group.by, # condition variable
                       only.pos=TRUE,
                       logfc.threshold = 0.25,
                       min.pct=0.1,
                       min.diff.pct = -Inf) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}




### get conserved markers across conditions
n = length(unique(integrated_seurat$harmony_snn_res.0.25))
# get conserved markers # adjust 0:## to (number of clusters - 1)
# conserved.pos.markers = map_dfr(c(0:n-1), get_conserved) # 28 clusters at resolution 0.8
conserved.pos.markers = map_dfr(c(1:n), get_conserved) # 21 clusters at harmony resolution 0.25. Leiden clustering starts from index 1


# add avg.pct.diff, max_adj_pval, and min.pct 1 cols
conserved.pos.markers = conserved.pos.markers %>%
  mutate(avg_log2FC = (ctrl_avg_log2FC + pb_avg_log2FC) / 2,
         max_adj_pval = ifelse(ctrl_p_val_adj > pb_p_val_adj,ctrl_p_val_adj,pb_p_val_adj),
         avg.pct.diff = ((pb_pct.1 - pb_pct.2) + (ctrl_pct.1 - ctrl_pct.2)) / 2,
         min.pct.1 = ifelse(ctrl_pct.1 < pb_pct.1, ctrl_pct.1,pb_pct.1),
         max.pct.2 = ifelse(ctrl_pct.2 > pb_pct.2, ctrl_pct.2,pb_pct.2))

# Rearrange the columns to be more intuitive
conserved.pos.markers = conserved.pos.markers[ , c(1:2,16:20,15,3:14)]

# Order the rows by log2FC values
conserved.pos.markers = conserved.pos.markers %>%
  dplyr::arrange(cluster_id, desc(avg_log2FC))

# save markers
write.csv(conserved.pos.markers, paste0(annot_dir_integrated,'conserved_pos_markers_across_condition_harmonyres.0.25_df.csv'),row.names=FALSE)


# Extract top 10 conserved markers per cluster
top20.conserved = conserved.pos.markers %>% 
  group_by(cluster_id) %>%
  filter(max_adj_pval < 0.05, avg.pct.diff > 0.5) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  # slice_max(order_by=tibble(avg.pct.diff),n=20)
slice_max(order_by=tibble(avg_log2FC,avg.pct.diff),n=20)
write.csv(top20.conserved, paste0(annot_dir_integrated,'top20.conserved_pos_markers_across_condition_harmonyres.0.25_df.csv'),row.names=FALSE)


# top10.conserved = top10.conserved[!top10.conserved$cluster_id == 26,] # remove because only one condition

gc() # free up RAM

### get conserved markers across timepoints ----
# get conserved markers # adjust 0:## to (number of clusters - 1)
# conserved.pos.markers_t = map_dfr(c(0:n-1), get_conserved(group.by='timepoint')) # 28 clusters at resolution 0.8
conserved.pos.markers_t = map_dfr(c(1:n), group.by='timepoint',get_conserved) # 21 clusters at harmony resolution 0.25. Leiden clustering starts from index 1

# add avg.pct.diff, max_adj_pval, and min.pct 1 cols
conserved.pos.markers_t = conserved.pos.markers_t %>%
  dplyr::mutate(
    avg_log2FC = (wk3_avg_log2FC + wk10_avg_log2FC + mn7_avg_log2FC + mn18_avg_log2FC) / 4,
    avg.pct.diff = ((wk3_pct.1 - wk3_pct.2) + (wk10_pct.1 - wk10_pct.2) + 
                      (mn7_pct.1 - mn7_pct.2) + (mn18_pct.1 - mn18_pct.2)) / 4
  ) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(
    max_adj_pval = max(c(wk3_p_val_adj, wk10_p_val_adj, mn7_p_val_adj, mn18_p_val_adj), na.rm = TRUE),
    min.pct.1 = min(c(wk3_pct.1, wk10_pct.1, mn7_pct.1, mn18_pct.1), na.rm = TRUE),
    max.pct.2 = max(c(wk3_pct.2, wk10_pct.2, mn7_pct.2, mn18_pct.2), na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Rearrange the columns to be more intuitive
conserved.pos.markers_t = conserved.pos.markers_t[ , c(1:2,26:30,25,3:24)]

# Order the rows by log2FC values
conserved.pos.markers_t = conserved.pos.markers_t %>%
  dplyr::arrange(cluster_id, desc(avg_log2FC))

# save markers
write.csv(conserved.pos.markers_t, paste0(annot_dir_integrated,'conserved_pos_markers_across_timepoints_harmonyres.0.25_df.csv'),row.names=FALSE)

top20.conserved_t = conserved.pos.markers_t %>% 
  group_by(cluster_id) %>%
  filter(max_adj_pval < 0.05, avg.pct.diff > 0.5) %>%
  # slice_max(order_by=tibble(avg.pct.diff),n=20)
slice_max(order_by=tibble(avg_log2FC),n=20)
write.csv(top20.conserved_t, paste0(annot_dir_integrated,'top20.conserved_pos_markers_across_timepoints_harmonyres.0.25_df.csv'),row.names=FALSE)

# # clusters 26 is made up of one condition.
# # Will treat individual values as avg and max to easily bind rows to top10.conserved
# topc26 = conserved.pos.markers[conserved.pos.markers$cluster_id == 26,] %>%
#   mutate(avg_fc = ctrl_avg_log2FC,
#          max_adj_pval = ctrl_p_val_adj,
#          avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
#   group_by(cluster_id) %>%
#   slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)
# 
# top10.conserved = rbind(top10.conserved, topc26)


# cell annotation ----
# integrated_seurat = qs_read(paste0(data_dir_integrated,'annotated_seurat_Li_etal_sctype.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden_doublets_x2.qs'),nthreads=nthreads)

# run T-SNE. Takes 20 minutes
# integrated_seurat = RunTSNE(integrated_seurat,
#                             reduction='harmony',
#                             dims=1:50,
#                             reduction.name='tsne_harmony')

# qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden_doublets_x2.qs'),nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden_doublets_x2.qs'),nthreads=nthreads)
# qs_save(integrated_seurat, paste0(data_dir_integrated,'annotated_seurat_Li_etal_sctype.qs'), nthreads=nthreads)
integrated_seurat = qs_read(paste0(data_dir_integrated,'annotated_seurat_Li_etal_sctype.qs'),nthreads=nthreads)
gc()

DimPlot(integrated_seurat,
        group.by='celltype',
        split.by='timepoint',
        reduction='umap_harmony',
        label=T,
        label.size=3,
        repel=T
)

DimPlot(integrated_seurat,
        group.by='gptcelltype_clusters',
        # split.by='timepoint',
        reduction='umap_harmony',
        label=T,
        label.size=5,
        repel=T
)


Idents(integrated_seurat) = integrated_seurat$harmony_snn_res.0.25
integrated_seurat = RenameIdents(integrated_seurat,
                                 '1' = 'T cells',
                                 '2' = 'B cells',
                                 '3' = 'Fibroblast',
                                 '4' = 'T cells',
                                 '5' = 'Endothelial',
                                 '6' = 'Adipocytes',
                                 '7' = 'Myeloid',
                                 '8' = 'Epithelial',
                                 '9' = 'T cells',
                                 '10' = 'Epithelial', 
                                 '11' = 'Myeloid',
                                 '12' = 'Epithelial',
                                 '13' = 'Myeloid',
                                 '14' = 'Pericytes',
                                 '15' = 'Proliferating cells',
                                 '16' = 'T cells',
                                 '17' = 'Fibroblast', 
                                 '18' = 'Endothelial',
                                 '19' = 'Myeloid', 
                                 '20' = 'Schwann cells', 
                                 '21' = 'Mast cells',
                                 '22' = 'Muscle'
)

integrated_seurat$major_celltype = Idents(integrated_seurat)
# qs_save(integrated_seurat, paste0(data_dir_integrated, 'integrated_seurat_by_batch_singlets_only_leiden_doublets_x2.qs'),nthreads=nthreads)


DimPlot(integrated_seurat,
        group.by='major_celltype',
        # split.by='timepoint',
        reduction='umap_harmony',
        label=T,
        label.size=5,
        repel=T
)
DimPlot(integrated_seurat,
        group.by='major_celltype',
        # split.by='timepoint',
        reduction='tsne_harmony',
        label=T,
        label.size=5,
        repel=T
)













#### sessionInfo ----
sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Red Hat Enterprise Linux 8.8 (Ootpa)
# 
# Matrix products: default
# BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRblas.so 
# LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: America/Detroit
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] clustree_0.5.1              ggraph_2.2.1               
# [3] scDblFinder_1.18.0          SingleCellExperiment_1.28.1
# [5] SummarizedExperiment_1.36.0 Biobase_2.66.0             
# [7] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
# [9] IRanges_2.40.1              S4Vectors_0.44.0           
# [11] BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
# [13] matrixStats_1.5.0           BiocParallel_1.40.0        
# [15] qs2_0.1.5                   presto_1.0.0               
# [17] data.table_1.17.0           pheatmap_1.0.12            
# [19] glmGamPoi_1.16.0            scales_1.3.0               
# [21] viridis_0.6.5               viridisLite_0.4.2          
# [23] plyr_1.8.9                  ggpubr_0.6.0               
# [25] ggthemes_5.1.0              lubridate_1.9.3            
# [27] forcats_1.0.0               stringr_1.5.1              
# [29] dplyr_1.1.4                 purrr_1.0.4                
# [31] readr_2.1.5                 tidyr_1.3.1                
# [33] tibble_3.2.1                ggplot2_3.5.1              
# [35] tidyverse_2.0.0             cowplot_1.1.3              
# [37] patchwork_1.3.0             harmony_1.2.3              
# [39] Rcpp_1.0.14                 sctransform_0.4.1          
# [41] Seurat_5.2.1                SeuratObject_5.0.2         
# [43] sp_2.2-0                   
# 
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0     bitops_1.0-9              httr_1.4.7               
# [4] RColorBrewer_1.1-3        numDeriv_2016.8-1.1       tools_4.4.0              
# [7] backports_1.5.0           DT_0.33                   R6_2.6.1                 
# [10] sn_2.1.1                  lazyeval_0.2.2            uwot_0.2.3               
# [13] withr_3.0.2               gridExtra_2.3             progressr_0.15.1         
# [16] textshaping_1.0.0         cli_3.6.4                 spatstat.explore_3.3-4   
# [19] fastDummies_1.7.5         sandwich_3.1-1            labeling_0.4.3           
# [22] mvtnorm_1.3-3             spatstat.data_3.1-4       ggridges_0.5.6           
# [25] pbapply_1.7-2             systemfonts_1.2.1         Rsamtools_2.20.0         
# [28] scater_1.32.1             parallelly_1.43.0         plotrix_3.8-4            
# [31] limma_3.62.2              rstudioapi_0.17.1         generics_0.1.3           
# [34] BiocIO_1.14.0             ica_1.0-3                 spatstat.random_3.3-2    
# [37] car_3.1-3                 Matrix_1.7-0              ggbeeswarm_0.7.2         
# [40] MetBrewer_0.2.0           abind_1.4-8               lifecycle_1.0.4          
# [43] multcomp_1.4-28           yaml_2.3.10               edgeR_4.2.2              
# [46] carData_3.0-5             mathjaxr_1.6-0            SparseArray_1.6.2        
# [49] Rtsne_0.17                grid_4.4.0                promises_1.3.2           
# [52] dqrng_0.4.1               crayon_1.5.3              miniUI_0.1.1.1           
# [55] lattice_0.22-6            msigdbr_10.0.0            beachmat_2.22.0          
# [58] pillar_1.10.1             metapod_1.12.0            rjson_0.2.23             
# [61] xgboost_1.7.8.1           future.apply_1.11.3       codetools_0.2-20         
# [64] mutoss_0.1-13             glue_1.8.0                leidenbase_0.1.32        
# [67] spatstat.univar_3.1-2     Rdpack_2.6.3              vctrs_0.6.5              
# [70] png_0.1-8                 spam_2.11-1               gtable_0.3.6             
# [73] cachem_1.1.0              rbibutils_2.3             S4Arrays_1.6.0           
# [76] mime_0.13                 tidygraph_1.3.1           survival_3.5-8           
# [79] statmod_1.5.0             bluster_1.14.0            TH.data_1.1-3            
# [82] fitdistrplus_1.2-2        ROCR_1.0-11               nlme_3.1-164             
# [85] RcppAnnoy_0.0.22          irlba_2.3.5.1             vipor_0.4.7              
# [88] KernSmooth_2.23-22        colorspace_2.1-1          mnormt_2.1.1             
# [91] ggrastr_1.0.2             tidyselect_1.2.1          compiler_4.4.0           
# [94] curl_6.2.2                BiocNeighbors_2.0.1       TFisher_0.2.0            
# [97] DelayedArray_0.32.0       plotly_4.10.4             stringfish_0.16.0        
# [100] rtracklayer_1.64.0        checkmate_2.3.2           lmtest_0.9-40            
# [103] digest_0.6.37             goftest_1.2-3             spatstat.utils_3.1-2     
# [106] RhpcBLASctl_0.23-42       XVector_0.46.0            htmltools_0.5.8.1        
# [109] pkgconfig_2.0.3           sparseMatrixStats_1.18.0  fastmap_1.2.0            
# [112] rlang_1.1.5               htmlwidgets_1.6.4         UCSC.utils_1.2.0         
# [115] shiny_1.10.0              DelayedMatrixStats_1.28.1 farver_2.1.2             
# [118] zoo_1.8-13                jsonlite_2.0.0            BiocSingular_1.22.0      
# [121] RCurl_1.98-1.17           magrittr_2.0.3            Formula_1.2-5            
# [124] scuttle_1.16.0            GenomeInfoDbData_1.2.13   dotCall64_1.2            
# [127] munsell_0.5.1             babelgene_22.9            reticulate_1.41.0.1      
# [130] stringi_1.8.7             zlibbioc_1.52.0           MASS_7.3-60.2            
# [133] parallel_4.4.0            listenv_0.9.1             ggrepel_0.9.6            
# [136] deldir_2.0-4              graphlayouts_1.2.2        Biostrings_2.72.1        
# [139] splines_4.4.0             multtest_2.60.0           tensor_1.5               
# [142] hms_1.1.3                 qqconf_1.3.2              locfit_1.5-9.12          
# [145] igraph_2.1.4              spatstat.geom_3.3-5       ggsignif_0.6.4           
# [148] RcppHNSW_0.6.0            reshape2_1.4.4            ScaledMatrix_1.14.0      
# [151] XML_3.99-0.18             metap_1.12                RcppParallel_5.1.10      
# [154] scran_1.32.0              tweenr_2.0.3              tzdb_0.4.0               
# [157] httpuv_1.6.15             RANN_2.6.2                polyclip_1.10-7          
# [160] future_1.34.0             scattermore_1.2           ggforce_0.4.2            
# [163] rsvd_1.0.5                broom_1.0.5               xtable_1.8-4             
# [166] restfulr_0.0.15           RSpectra_0.16-2           rstatix_0.7.2            
# [169] later_1.4.1               ragg_1.3.3                memoise_2.0.1            
# [172] beeswarm_0.4.0            GenomicAlignments_1.40.0  cluster_2.1.6            
# [175] timechange_0.3.0          globals_0.16.3