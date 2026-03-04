library(Seurat)
library(SeuratObject)
library(sctransform)
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

# set number of threads to use for saving qs formatted data
nthreads <- 2 

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
    fig_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/",analysis,"/FPR_0.0_MAD/annotation/")
  ),
  "week10" = list(
    fig_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/",analysis,"/FPR_0.0_MAD/annotation/")
  ),
  "month7" = list(
    fig_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_month/",analysis,"/FPR_0.0_MAD/annotation/")
  ),
  "month18" = list(
    fig_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/clustering/"),
    data_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/data/"),
    annot_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/annotation/")
  )
)
# timepoint_dirs <- list(
#   "week3" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   ),
#   "week10" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   ),
#   "month7" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_month/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   ),
#   "month18" = list(
#     fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/clustering/",
#     data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/data/",
#     annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/annotation/"
#   )
# )
# list2env(timepoint_dirs[[selected_timepoint]], envir = environment())

## directories for integrated data ----
fig_dir_integrated = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir_integrated = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir_intergrated = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")

# Create directories to save results if they don't already exist:
dirs = c(fig_dir_integrated, data_dir_integrated, annot_dir_intergrated)

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
integrated_seurat == JoinLayers(integrated_seurat)

## add timepoint metadata column ----
integrated_seurat$timepoint = str_extract(integrated_seurat$cells, '^[^_]+')

rm(seurat.list) # free RAM
gc()

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat.qs'), nthreads=nthreads)

# exploring sources of variation, cell cycle and mitochondrial gene expression----

# cell cycle genes, formatted to match genes in matrix with first letter upper-cased
s.genes = stringr::str_to_title(cc.genes$s.genes)
g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)


# Normalize the counts
integrated_seurat = NormalizeData(integrated_seurat)



# Score cells for cell cycle
integrated_seurat = JoinLayers(integrated_seurat, assay = 'RNA')
integrated_seurat = CellCycleScoring(integrated_seurat, 
                                    g2m.features = g2m.genes, 
                                    s.features = s.genes,
                                    set.ident=T)

# From https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Seurat suggests regressing out the difference between S and G2M phase scores
# in tissues with differentiating cells, like our developing mammary gland.
integrated_seurat$cc.difference = integrated_seurat$S.Score - integrated_seurat$G2M.Score


# Check quartile values of mitochondrial gene expression
summary(integrated_seurat@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.0000  0.4739  0.8186  1.1372  9.9042 

# Turn percent.mt into categorical factor vector based on quartile values
integrated_seurat@meta.data$mito.level = cut(integrated_seurat@meta.data$percent.mt, 
                                        breaks=c(-Inf, 0,0.4739, 1.1372, Inf), 
                                        labels=c("1st","2nd","3rd", "4th"))

# no longer need data log-normalized layer
integrated_seurat = DietSeurat(integrated_seurat,assays='RNA',layers='counts')
gc() # free RAM

qs_save(integrated_seurat, paste0(data_dir_integrated,'integrated_seurat.qs'), nthreads=nthreads)
# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat.qs'), nthreads=nthreads)
gc() # free RAM

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

# test = qs_read('C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/data/annotated_seurat_res.0.8_df.qs',nthreads=nthreads)


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
ggsave(paste0(fig_dir, '/PCA_by_cellcycle_phase_allgenes_.png'), dpi=320, width=20, height= 10)


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
ggsave(paste0(fig_dir,'/PCA_by_mitochondrial_expression_quartile_all_genes.png'), dpi=320, width=20, height= 10)

# qs_save(seurat_pca, paste0(data_dir,'integrated_seurat_rm29and32_regressout_cc.diff.qs'),nthreads=nthreads)# save integrated_seurat object
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
ggsave(paste0(fig_dir,'/PCA_by_condition_all_genes.png'), dpi=320, width=20, height= 10)


rm(seurat_pca,p) # free RAM. Don't need any further
gc()
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

# merge sctransform normalized samples. Takes a little while to run
integrated_seurat = merge(split_seurat[[1]], split_seurat[-1], merge.data=TRUE,merge.dr=FALSE)
integrated_seurat = JoinLayers(seurat_merged_sct, assay='RNA')
gc() # free RAM


qs_save(seurat_merged_sct, paste0(data_dir,'seurat_merged_sct_regressout_cc.diffandmt.qs'),nthreads=nthreads)

gc() #free up RAM

### continue analysis without integration using sctransform normalization ----
# select integration features according to https://github.com/satijalab/seurat/issues/5205
# and https://github.com/satijalab/seurat/issues/6185
# selecting 5000 features to account for the possible different rankings of highly variable genes
# between samples. Since we're normalizing with SCTtransform this # shouldn't change things much.
features = SelectIntegrationFeatures(split_seurat, nfeatures=5000, fvf.nfeatures=5000)
qs_save(features,paste0(data_dir,'variable_features.qs'),nthreads=nthreads)
VariableFeatures(seurat_merged_sct) = features

sct_results = SCTResults(seurat_merged_sct, slot = "feature.attributes")
qs_save(sct_results, paste0(data_dir,'sct_results_regressout_cc.diffandmt.qs'),nthreads=nthreads)
rm(sct_results)

# Variable features lost/unset after merging transformed objects. Need to set again.
# should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(seurat_merged_sct[["SCT"]]) = rownames(seurat_merged_sct[["SCT"]]@scale.data)

# can plot variable genes by plotting gmean vs residual variance
# according to https://github.com/satijalab/seurat/issues/8321
# z = seurat_merged_sct@assays$SCT@SCTModel.list[[3]]@feature.attributes
# qplot(z$gmean, z$residual_variance) + scale_y_log10() + scale_x_log10()

# Identify the 25 most highly variable genes
top25_sct = head(VariableFeatures(seurat_merged_sct), 25)
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
  ggsave(paste0(fig_dir,'top25_variable_genes_',n[[i]],'.png'), dpi=320, width=20, height= 10)
  
}, x=split_seurat, n=names(split_seurat))

rm(split_seurat) # free RAM
gc()

# Run PCA on non-integrated merged SCTransformed data ---- 
seurat_merged_sct = RunPCA(seurat_merged_sct, assay='SCT', npcs=50) # default params
gc() # free up RAM

qs_save(seurat_merged_sct, paste0(data_dir,'seurat_merged_sct_regressout_cc.diffandmt.qs'),nthreads=nthreads)


### Visualize expression of most highly weighted genes for PCs 1-50 (loadings) ----
VizDimLoadings(seurat_merged_sct, dims = 1:12, ncol=4)
ggsave(paste0(fig_dir,'PCA_loadings_sct_1-12.png'), dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 13:24, ncol=4)
ggsave(paste0(fig_dir,'PCA_loadings_sct_13-24.png'), dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 25:36, ncol=4)
ggsave(paste0(fig_dir,'PCA_loadings_sct_25-36.png'), dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 37:43, ncol=4)
ggsave(paste0(fig_dir,'PCA_loadings_sct_37-43.png'), dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 44:50, ncol=4)
ggsave(paste0(fig_dir,'PCA_loadings_sct_44-50.png'), dpi=320, width=20, height= 15)

png(filename = paste0(fig_dir,'high_score_genes_per_PC1-6_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 1:6,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC7-12_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 7:12,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC13-18_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 13:18,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC19-24_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 19:24,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC25-30_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 25:30,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC31-36_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 31:36,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC37-42_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 37:42,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = paste0(fig_dir,'high_score_genes_per_PC43-50_heatmap_sct.png'), units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 43:50,
           ncol=4,
           cells=500,
           balanced = TRUE)
dev.off()

# ranked list (PC score) of genes for PC1
names(sort(Loadings(seurat_merged_sct, reduction="pca")[,1], decreasing=TRUE))


# Printing out the top 10 most variable genes driving PCs
print(x = seurat_merged_sct[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

### Plot the elbow plot to see variation captured by PCs ----
ElbowPlot(object = seurat_merged_sct, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave(paste0(fig_dir,'elbow_plot_sct.png'), dpi=320, width=20, height= 10)

# We can calculate where the principal components start to elbow by taking the smaller value of:
#
# 1)The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# 2)The point where the percent change in variation between the consecutive PCs is less than 0.1%.

### Determine percent of variation associated with each PC ----
pct = seurat_merged_sct[["pca"]]@stdev / sum(seurat_merged_sct[["pca"]]@stdev) * 100
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
ggsave(paste0(fig_dir,'quantitative_elbow_plot_sct.png'), dpi=320, width=20, height= 10)
rm(plot_df)

# cluster sctransformed integrated cells----
resolutions = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0)

seurat_umap_sct = JoinLayers(seurat_umap_sct, assay='RNA')

seurat_umap_sct = IntegrateLayers(
  object = seurat_umap_sct,
  method = RPCAIntegration,
  orig.reduction = 'pca',
  new.reduction = 'integrated.dr',
  normalization.method = "SCT",
  verbose = T
)

seurat_umap_sct <- FindNeighbors(seurat_umap_sct, dims = 1:50, reduction = "integrated.dr",
                                 graph.name = c('idr_nn','idr_snn'))
seurat_umap_sct <- FindClusters(seurat_umap_sct, resolution = resolutions,
                                graph.name = 'idr_snn')
seurat_umap_sct = RunUMAP(seurat_umap_sct,
                          dims=1:50,
                          reduction='integrated.dr',
                          reduction.name = 'umap_idr') 
gc() # free RAM

qs_save(seurat_umap_sct, paste0(data_dir, 'seurat_umap_sct_integrated_layers.qs'),nthreads=nthreads)

# seurat_umap_sct = qs_read(paste0(data_dir, 'seurat_umap_sct_integrated_layers.qs'),nthreads=nthreads)


DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'Phase',
        label=F,
        label.size=6,
        shuffle=TRUE )     

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'sample',
        label=F,
        label.size=6,
        shuffle=TRUE )  


DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'condition',
        split.by='condition',
        label=F,
        label.size=6,
        shuffle=TRUE)

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'mito.level',
        label=F,
        label.size=6,
        shuffle=TRUE)

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'cell_calls',
        # split.by = 'cell_calls',
        label=F,
        label.size=6,
        shuffle=TRUE)

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'scDblFinder.class',
        # split.by = 'scDblFinder.class',
        label=F,
        label.size=6,
        shuffle=TRUE)

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'SCT_snn_res.0.1',
        label=F,
        label.size=6,
        shuffle=TRUE)


DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'idr_snn_res.0.8',
        label=F,
        label.size=6,
        shuffle=TRUE)



# cluster sctransformed non-integrated cells----
seurat_umap_sct = FindNeighbors(seurat_umap_sct, dims=1:50, reduction='pca')
seurat_umap_sct = FindClusters(object = seurat_umap_sct,
                           resolution = resolutions)

qs_save(seurat_umap_sct, paste0(data_dir,'seurat_umap_sct_regressout_cc.diffandmt.qs'),nthreads=nthreads)
rm(seurat_merged_sct) # free RAM
gc()


#### run UMAP on non-integrated, merged, sct-normalized data ----
seurat_umap_sct = RunUMAP(seurat_umap_sct,
                      dims=1:50,
                      reduction='pca') 



sce_cb = as.SingleCellExperiment(seurat_umap_sct)
sce_cb = scDblFinder(sce_cb, samples = "sample", 
                     clusters = 'SCT_snn_res.0.8',
                     # clusters = 'cluster',
                     dbr.sd = 1, # lets the dbr be set emperically
                     dbr.per1k = 0.008, # default is 0.008
                     BPPARAM=SnowParam(4)
) #dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
#dbr is 0.8% per 1000 cells according to 10x reps (0.08% = 0.0008)
# which, translates to 0.05% for a 16-plex experiment according to DOI: 0.1101/2024.10.03.616596

# they calculate a fitted fraction of 0.12% for 16-plex
# when running multiple samples we recommend to first cluster all samples 
# together (for example using sce$cluster <- fastcluster(sce)) and then 
# provide the clusters to scDblFinder.


# register(SerialParam()) # close extra nodes/workers used
# bpstop()
table(sce_cb$scDblFinder.class)
metadata(sce_cb)$scDblFinder.threshold


hist(sce_cb$scDblFinder.score, breaks = 100)

table(sce_cb$sample, sce_cb$scDblFinder.class)  # Rates should correlate with cell counts


hist(sce_cb$scDblFinder.weighted, breaks = 100)
hist(sce_cb$scDblFinder.difficulty, breaks = 100)
hist(sce_cb$scDblFinder.cxds_score, breaks = 100)
table(sce_cb$scDblFinder.mostLikelyOrigin)
table(sce_cb$scDblFinder.originAmbiguous)








qs_save(seurat_umap_sct, paste0(data_dir,'seurat_umap_sct_regressout_cc.diffandmt_include_doublets.qs'),nthreads=nthreads)
# seurat_umap_sct = qs_read(paste0(data_dir,'seurat_umap_sct_rm29and32_regressout_cc.diffandmt_.qs'),nthreads=nthreads)

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'Phase',
        label=F,
        label.size=6,
        shuffle=TRUE )

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'condition',
        label=F,
        label.size=6,
        shuffle=TRUE)

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'mito.level',
        label=F,
        label.size=6,
        shuffle=TRUE)


Idents(seurat_umap_sct) = 'SCT_snn_res.0.1'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.2'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.4'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.6'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.8'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1.2'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.2'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

# start at 0.2
# 0.4 maybe
# 0.8 looks like best to start with
# 1.2 is the highest res to use probably
# 1.4 has same clusters as 1.2, but slightly different shape
# between cluster 19 and 12, which are the same cluster comparatively speaking.
# 1.2 has a better defined cluster shape
# 1. just has one less B cell cluster than 1.2 and 1.4

# 
# going with 0.8

# Visualizations of clustering----
## visualize clustering results with clustree ----
library(clustree)

clustree(seurat_umap_sct, prefix='SCT_snn_res.')
ggsave(paste0(fig_dir,'clustree.png'), dpi=320, width=20, height= 10)

clustree(seurat_umap_sct, prefix='idr_snn_res.')
ggsave(paste0(fig_dir,'clustree_integrated.png'), dpi=320, width=20, height= 10)

# clustree(seurat_umap_sct, prefix = 'SCT_snn_res.',
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
               shuffle=T)
  
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
               shuffle=T)
  
  p1 | p2
  ggsave(filename=paste0(directory,resolution,'/',reduction_name,'_by_cellcycle_and_mito.quartile.png'), dpi=320, width=20, height=10)
  
}

### run plot_umaps ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory=fig_dir)
plot_umaps(seurat_umap=seurat_umap_sct, reduction='umap_idr', resolution='idr_snn_res.0.8', directory=fig_dir)


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
# umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory=fig_dir)
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory=fig_dir)
umap_clustering_qc(seurat_umap=seurat_umap_sct, reduction='umap_idr',resolution='idr_snn_res.0.8', directory=fig_dir)


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
  p1 = DimPlot(seurat_umap_sct,
               reduction = reduction,
               group.by = resolution,
               label=T,
               repel=T,
               label.size=6,
               shuffle=TRUE)
  
  p2 = DimPlot(seurat_umap_sct,
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
umap_doublets_qc(seurat_umap_sct, resolution='SCT_snn_res.0.8', directory=fig_dir)
umap_doublets_qc(seurat_umap_sct, reduction='umap_idr',resolution='idr_snn_res.0.8', directory=fig_dir)


### try scDblFidner on seraut clusters
sce_cb = as.SingleCellExperiment(integrated_seurat_cb)
sce_cb$cluster = fastcluster(sce_cb)
sce_cb = scDblFinder(sce_cb, samples = "sample", 
                     clusters = '',
                     # clusters = 'cluster',
                     # dbr.sd = 1, # lets the dbr be set emperically
                     # dbr.per1k = 0.0012, # default is 0.008
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
png(filename = paste0(fig_dir,'histogram_scDblFinder_scores.png'), units='in', width=10, height=10, res=320)
dev.off()

hist(sce_cb$scDblFinder.cxds_score, breaks = 100)
table(sce_cb$scDblFinder.mostLikelyOrigin)
table(sce_cb$scDblFinder.originAmbiguous)






### subset doublets out of dataset ----




# tabulate doublets by sample
table(seurat_umap_sct$sample, seurat_merged_sct$scDblFinder.class)


# tabulate total doublets
table(seurat_umap_sct$scDblFinder.class)

# seurat_merged_sct = subset(seurat_merged_sct, subset=dscDblFinder.class == 'singlet')





## Which PCs drive different clusters? ----

# determine max number of PCs present in seurat object
max_pcs = length(grep("^PC_", colnames(seurat_umap_sct@reductions$pca@cell.embeddings)))

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
# pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory=fig_dir)
pcs_umap_overlay_optimized(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory=fig_dir)


# PC_7 looks like it's driving cell cycle genes
# PC_13 looks like it's driving myoepithelial cells
# PC_14 marcophage mb
# PC_17 mystery cluster by endothelial cells
# PC_43 and PC_46 maybe cells transitioning between fibroblast and epithelial 


#### exploring cell marker genes ----
# library(readxl)
# 
# # read in cell markers from excel file
# cell_markers = read_excel('C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/annotation/aggregate_sources_celltype_markers.xlsx')
# cell_markers = lapply(cell_markers, na.omit) # convert to list of lists, omit NAs
# 
# # set RNA counts slot to be default and normalize for visualizing genes
# # SCTransform normalization was performed only on the 3000 most variable genes, 
# # so many of our genes of interest may not be present in the SCT data. Hence,
# # we switch over to the original RNA count to explore cell type marker genes.
# DefaultAssay(seurat_umap_sct) = 'RNA'
# # seurat_umap_sct = NormalizeData(seurat_umap_sct, verbose = FALSE)
# # DefaultAssay(seurat_umap_sct) = 'SCT'
# # qs_save(seurat_umap_sct, paste0(data_dir,'seurat_umap_sct.qs'),nthreads=nthreads)
# 
# # Want to only consider cell markers that are actually present in the dataset
# cell_markers = lapply(cell_markers, function(x){
#   x = unique(x[x %in% rownames(seurat_umap_sct)])
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
# plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory=fig_dir)
# 
# # looks like we could benefit from increasing resolution
# # resolutions available. For reference
# resolutions = sprintf('SCT_snn_res.%.1f',c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 2.0))
# 
# # Assign identity of clusters to desired resolution
# Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'
# 
# DimPlot(seurat_umap_sct,
#         reduction = "umap",
#         label=TRUE,
#         label.size=6)



DimPlot(seurat_umap_sct,
        group.by='SCT_snn_res.0.4',
        reduction='umap',
        label=T,
        repel=T,
        shuffle=T,
        )

#### annotate clusters ----
# will annotate at a resolution of 0.8 and subcluster as needed
Idents(seurat_umap_sct) = 'SCT_snn_res.0.8' # set cluster #s as ident

DefaultAssay(seurat_umap_sct) = "RNA"

# Find marker genes ----

## Find all positive markers ----
# join layers to be able to run FindAllMarkers()

all.pos.markers = FindAllMarkers(seurat_umap_sct,
                                 only.pos=TRUE,
                                 logfc.threshold = 0.25,
                                 min.pct=0.1,
                                 min.diff.pct = -Inf) 

# Combine markers with gene descriptions 
annotations = read.csv('C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/annotation/annotationhub_mice_genes.csv')
ann.pos.markers = left_join(x = all.pos.markers, 
                          y = annotations[, c("gene_name", "description")],
                          by = c("gene" = "gene_name")) %>%
  unique()

# add pct.diff col
ann.pos.markers$pct.diff = ann.pos.markers$pct.1 - ann.pos.markers$pct.2

# Rearrange the columns to be more intuitive
ann.pos.markers = ann.pos.markers[ , c(6, 7, 2:4, 9, 1, 5,8)]

# Order the rows by p-adjusted values
ann.pos.markers = ann.pos.markers %>%
  dplyr::arrange(cluster, p_val_adj)

# save markers
write.csv(ann.pos.markers,paste0(annot_dir,'annotated_pos_markers_res.0.8_df.csv'),row.names=FALSE)

# Extract top 10 markers per cluster
top10 = ann.pos.markers %>% 
  group_by(cluster) %>%
  slice_max(order_by=tibble(avg_log2FC,pct.diff,p_val_adj),n=10)


## Find all positive conserved markers between conditions ----
# Create function to get conserved markers for any given cluster
get_conserved = function(cluster){
  FindConservedMarkers(seurat_umap_sct, # need to specify exact seurat object here
                       ident.1 = cluster,
                       grouping.var = "condition", # condition variable
                       only.pos=TRUE,
                       logfc.threshold = 0.25,
                       min.pct=0.1,
                       min.diff.pct = -Inf) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}


# get conserved markers # adjust 0:## to (number of clusters - 1)
conserved.pos.markers = map_dfr(c(0:27), get_conserved) # 28 clusters at resolution 0.8

# add avg.pct.diff, max_adj_pval, and min.pct 1 cols
conserved.pos.markers = conserved.pos.markers %>%
  mutate(avg_log2FC = (ctrl_avg_log2FC + pb_avg_log2FC) / 2,
         max_adj_pval = ifelse(ctrl_p_val_adj > pb_p_val_adj,ctrl_p_val_adj,pb_p_val_adj),
         avg.pct.diff = ((pb_pct.1 - pb_pct.2) + (ctrl_pct.1 - ctrl_pct.2)) / 2,
         min.pct.1 = ifelse(ctrl_pct.1 < pb_pct.1, ctrl_pct.1,pb_pct.1))

# Rearrange the columns to be more intuitive
conserved.pos.markers = conserved.pos.markers[ , c(1:2,16:19,15,3:14)]

# Order the rows by p-adjusted values
conserved.pos.markers = conserved.pos.markers %>%
  dplyr::arrange(cluster_id, max_adj_pval)

# save markers
write.csv(conserved.pos.markers, paste0(annot_dir,'conserved_pos_markers_res.0.8_df.csv'),row.names=FALSE)


# Extract top 10 conserved markers per cluster
top10.conserved = conserved.pos.markers %>% 
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_log2FC,avg.pct.diff,max_adj_pval),n=10)


top10.conserved = top10.conserved[!top10.conserved$cluster_id == 26,] # remove because only one condition

gc() # free up RAM


#### sessionInfo ----
sessionInfo()
# 
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DoubletFinder_2.0.4 presto_1.0.0        data.table_1.16.0   Rcpp_1.0.13         pheatmap_1.0.12     Polychrome_1.5.1    RCurl_1.98-1.16     metap_1.11          glmGamPoi_1.14.3   
# [10] scales_1.3.0        viridis_0.6.5       viridisLite_0.4.2   plyr_1.8.9          ggpubr_0.6.0        ggthemes_5.1.0      Matrix_1.6-5        lubridate_1.9.3     forcats_1.0.0      
# [19] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2         readr_2.1.5         tidyr_1.3.1         tibble_3.2.1        ggplot2_3.5.1       tidyverse_2.0.0     cowplot_1.1.3      
# [28] patchwork_1.3.0     sctransform_0.4.1   Seurat_5.1.0        SeuratObject_5.0.2  sp_2.1-4           
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                polyclip_1.10-7             fastDummies_1.7.4          
# [7] lifecycle_1.0.4             Rdpack_2.6.1                rstatix_0.7.2               globals_0.16.3              lattice_0.22-6              MASS_7.3-60.0.1            
# [13] backports_1.5.0             magrittr_2.0.3              plotly_4.10.4               plotrix_3.8-4               qqconf_1.3.2                httpuv_1.6.15              
# [19] sn_2.1.1                    spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0           pbapply_1.7-2               RColorBrewer_1.1-3         
# [25] multcomp_1.4-26             abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1        BiocGenerics_0.48.1        
# [31] TH.data_1.1-2               sandwich_3.1-1              GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6              
# [37] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-0        TFisher_0.2.0               goftest_1.2-3               RSpectra_0.16-2            
# [43] spatstat.random_3.3-2       fitdistrplus_1.2-1          parallelly_1.38.0           leiden_0.4.3.1              codetools_0.2-20            DelayedArray_0.28.0        
# [49] tidyselect_1.2.1            farver_2.1.2                matrixStats_1.4.1           stats4_4.3.1                spatstat.explore_3.3-2      mathjaxr_1.6-0             
# [55] jsonlite_1.8.9              multtest_2.58.0             progressr_0.14.0            Formula_1.2-5               ggridges_0.5.6              survival_3.7-0             
# [61] tools_4.3.1                 ica_1.0-3                   glue_1.7.0                  mnormt_2.1.1                gridExtra_2.3               SparseArray_1.2.4          
# [67] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.8         numDeriv_2016.8-1.1         withr_3.0.1                 fastmap_1.2.0               fansi_1.0.6                
# [73] digest_0.6.37               timechange_0.3.0            R6_2.5.1                    mime_0.12                   colorspace_2.1-1            scattermore_1.2            
# [79] tensor_1.5                  spatstat.data_3.1-2         utf8_1.2.4                  generics_0.1.3              httr_1.4.7                  htmlwidgets_1.6.4          
# [85] S4Arrays_1.2.1              scatterplot3d_0.3-44        uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.5                lmtest_0.9-40              
# [91] XVector_0.42.0              htmltools_0.5.8.1           carData_3.0-5               dotCall64_1.1-1             Biobase_2.62.0              png_0.1-8                  
# [97] spatstat.univar_3.0-1       rstudioapi_0.16.0           tzdb_0.4.0                  reshape2_1.4.4              nlme_3.1-166                zoo_1.8-12                 
# [103] KernSmooth_2.23-24          parallel_4.3.1              miniUI_0.1.1.1              pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                
# [109] RANN_2.6.2                  promises_1.3.0              car_3.1-3                   xtable_1.8-4                cluster_2.1.6               mvtnorm_1.3-1              
# [115] cli_3.6.3                   compiler_4.3.1              rlang_1.1.4                 crayon_1.5.3                mutoss_0.1-13               future.apply_1.11.2        
# [121] ggsignif_0.6.4              stringi_1.8.4               deldir_2.0-4                munsell_0.5.1               lazyeval_0.2.2              spatstat.geom_3.3-3        
# [127] RcppHNSW_0.6.0              hms_1.1.3                   future_1.34.0               shiny_1.9.1                 SummarizedExperiment_1.32.0 rbibutils_2.2.16  