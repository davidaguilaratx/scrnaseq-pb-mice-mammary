setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(SeuratObject)
library(sctransform)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(glmGamPoi)
library(metap)
library(RCurl)
library(Polychrome)
library(pheatmap)
library(presto)
library(DESeq2)

# set seed
set.seed(2024)

# sessionInfo()

# load in filtered data
filtered_seurat = readRDS('data/filtered_seurat.rds')


### create polychrome color palette ----
p12 = createPalette(12, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
swatch(p12)
names(p12) = NULL

### exploring sources of variation, cell cycle and mitochondrial gene expression----

# cell cycle genes, formatted to match genes in matrix with first letter upper-cased
s.genes = stringr::str_to_title(cc.genes$s.genes)
g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)


# Normalize the counts
seurat_phase = NormalizeData(filtered_seurat)

# Score cells for cell cycle
seurat_phase = CellCycleScoring(seurat_phase, 
                                g2m.features = g2m.genes, 
                                s.features = s.genes)

# Identify the most variable genes
seurat_phase = FindVariableFeatures(seurat_phase, 
                                    selection.method = "vst", # default
                                    nfeatures = 2000, # default
                                    verbose = FALSE)

# Identify the 25 most highly variable genes
top25 = head(VariableFeatures(seurat_phase), 25)
# plot variable features with and without labels
p = VariableFeaturePlot(seurat_phase)
LabelPoints(plot=p, points = top25, repel = TRUE)+
  labs(title='Top 25 variably expressed genes')+
  theme_classic2()
ggsave('figures/clustering/top25_variable_genes_seurat_phase.png', dpi=320, width=20, height= 10)


# Scale the counts
seurat_phase = ScaleData(seurat_phase)

### Perform PCA on merged non-integrated data ----
seurat_phase = RunPCA(seurat_phase, npcs=50) # default params

#### Plot the PCA colored by cell cycle phase ----
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")+
  labs(title='PCA by cell cycle phase',
       color='phase')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering/PCA_by_cellcycle_phase.png', dpi=320, width=20, height= 10)


# Check quartile values of mitochondrial gene expression
summary(seurat_phase@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3610  0.7853  1.0592  1.4286  9.8712 

# Turn percent.mt into categorical factor vector based on quartile values
seurat_phase@meta.data$mito.level = cut(seurat_phase@meta.data$percent.mt, 
                                    breaks=c(-Inf, 0.3610, 0.7853, 1.0592, Inf), 
                                    labels=c("1st","2nd","3rd", "4th"))

### Plot the PCA colored by mitochondrial expression quartile ----
DimPlot(seurat_phase,
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
ggsave('figures/clustering/PCA_by_mitochondrial_expression_quartile.png', dpi=320, width=20, height= 10)

saveRDS(seurat_phase, 'data/seurat_phase.rds')# save seurat_phase object
##looking at the plots there doesn't seem to be any significant source of variation from cell cycle phase or mitochondrial expression

### plot the PCA colored by condition ----
DimPlot(seurat_phase,
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
ggsave('figures/clustering/PCA_by_condition.png', dpi=320, width=20, height= 10)

### Split seurat_phase by sample and apply sctransform to each individually ----
# Split seurat object by sample to perform cell cycle scoring and SCTransform on all samples
split_seurat = SplitObject(seurat_phase, split.by = "sample")

rm(seurat_phase) # save RAM
gc()

# adjust object size limits of variables/objects in R
options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB

#### apply SCTransform to each of the 12 samples, then merge normalized samples and save object ----
# takes a while to run. Do for each sample separately according to https://github.com/satijalab/seurat/issues/6003
split_seurat = lapply(split_seurat, function(x) {
  x = SCTransform(x, assay='RNA', vst.flavor ="v2", variable.features.n=3000) # default params
})
saveRDS(split_seurat, 'data/split_seurat.rds')

# merge sctransform normalized samples. Takes a while to run
seurat_merged_sct = merge(split_seurat[[1]], split_seurat[-1], merge.data=TRUE,merge.dr=FALSE) 
saveRDS(seurat_merged_sct, 'data/seurat_merged_sct.rds')

gc() #free up RAM

### continue analysis without integration using sctransform normalization ----
# select integration features according to https://github.com/satijalab/seurat/issues/5205
features = SelectIntegrationFeatures(split_seurat, nfeatures=5000, fvf.nfeatures=5000)
VariableFeatures(seurat_merged_sct) = features

SCTResults(seurat_merged_sct, slot = "feature.attributes")

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
  ggsave(paste0('figures/clustering/top25_variable_genes_seurat_phase_',n[[i]],'.png'), dpi=320, width=20, height= 10)
  
}, x=split_seurat, n=names(split_seurat))

rm(split_seurat) # free RAM
gc()

### Run PCA on non-integrated merged SCTransformed data ---- 
seurat_merged_sct = RunPCA(seurat_merged_sct, assay='SCT', npcs=50) # default params
gc() # free up RAM

saveRDS(seurat_merged_sct, 'data/seurat_merged_sct.rds')


### Visualize expression of most highly weighted genes for PCs 1-50 (loadings) ----
VizDimLoadings(seurat_merged_sct, dims = 1:12, ncol=4)
ggsave('figures/clustering/sct/PCA_loadings_sct.png', dpi=320, width=20, height= 15)


png(filename = 'figures/clustering/sct/high_score_genes_per_PC1-6_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 1:6,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC7-12_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 7:12,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC13-18_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 13:18,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC19-24_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 19:24,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC25-30_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 25:30,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC31-36_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 31:36,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC37-42_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 37:42,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/sct/high_score_genes_per_PC43-50_heatmap_sct.png', units='in', width=20, height=10, res=320)
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
ggsave('figures/clustering/sct/elbow_plot_sct.png', dpi=320, width=20, height= 10)

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
co2 # 19
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 19
pcs = 19

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
ggsave('figures/clustering/sct/quantitative_elbow_plot_sct.png', dpi=320, width=20, height= 10)


### cluster sctransformed non-integrated cells----
seurat_umap_sct = FindNeighbors(seurat_merged_sct, dims=1:50)
seurat_umap_sct = FindClusters(object = seurat_umap_sct,
                           resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0))

## clusters after running with dims=1:50 filtered nUMI<10000, removed low quality samples 29 and 31
# Resolution 0.1: 14 clusters
# Resolution 0.2: 18 clusters
# Resolution 0.4: 22 clusters
# Resolution 0.6: 26 clusters
# Resolution 0.8: 28 clusters
# Resolution 1.0: 30 clusters
# Resolution 1.2: 31 clusters
# Resolution 1.4: 32 clusters
# Resolution 2.0: 39 clusters


## clusters after running with dims=1:50. filtered nUMI<50000
# Resolution 0.1: 12 clusters
# Resolution 0.2: 17 clusters
# Resolution 0.4: 20 clusters
# Resolution 0.6: 25 clusters
# Resolution 0.8: 27 clusters
# Resolution 1.0: 30 clusters
# Resolution 1.2: 31 clusters
# Resolution 1.4: 33 clusters
# Resolution 2.0: 39 clusters



saveRDS(seurat_umap_sct, 'data/seurat_umap_sct.rds')
rm(seurat_merged_sct) # free RAM
gc()


### run UMAP on non-integrated, merged, sct-normalized data ----
seurat_umap_sct = RunUMAP(seurat_umap_sct,
                      dims=1:50,
                      reduction='pca') 

saveRDS(seurat_umap_sct, 'data/seurat_umap_sct.rds')
# seurat_umap_sct = readRDS('data/seurat_umap_sct.rds')

Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'
DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'Phase',
        label=TRUE,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'
DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'Phase',
        label=F,
        label.size=6)




## function to plot initial umaps ----
plot_umaps = function (seurat_umap, resolution, directory) {
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # Assign umap resolution
  Idents(seurat_umap) = resolution
  
  # Plot the UMAP of non-integrated data
  DimPlot(seurat_umap,
          reduction = "umap",
          label=TRUE,
          label.size=6)
  ggsave(paste0(directory,resolution,'/UMAP.png'), dpi=320, width=10, height=10)
  
  
  # plot UMAP by sample and condition
  p1 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='sample',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='condition',
               shuffle=T)
  
  png(filename=paste0(directory,resolution,'/UMAP_by_sample_and_condition.png'), units='in', width=20, height=10, res=320)
  print(plot_grid(p1, p2))
  dev.off()
  
  
  # plot UMAP by Phase and mitochondrial expression quartile
  p1 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='Phase',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='mito.level',
               shuffle=T)
  
  png(filename=paste0(directory,resolution,'/UMAP_by_cellcycle_and_mitoqrtl.png'), units='in', width=20, height=10, res=320)
  print(plot_grid(p1, p2))
  graphics.off()
}

## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.1 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering/sct/')


#### UMAP clustering QC ----
metrics =  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt") # metrics to overlay on umap

## create function to plot all QC plots and umap overlays ----
umap_clustering_qc = function(seurat_umap, resolution, directory) {
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # set umap resolution
  Idents(seurat_umap) = resolution
  seurat_umap@meta.data$resolution = Idents(seurat_umap)
  
  # Explore whether clusters segregate by cell cycle phase
  DimPlot(seurat_umap,
          label = TRUE, 
          split.by = "Phase")  + NoLegend()
  ggsave(paste0(directory,resolution,'/cell_cycle_split_umap.png'), dpi=320, width=30, height=20)
  
  # Metrics to plot on UMAP
  FeaturePlot(seurat_umap, 
              reduction = "umap", 
              features = metrics,
              pt.size = 0.4, 
              order = TRUE,
              min.cutoff = 'q10',
              label = TRUE)
  ggsave(paste0(directory,resolution,'/counts_cycle_mt_metrics.png'), dpi=320, width=20, height=20)
  
  # Boxplot of nUMIs per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nCount_RNA, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/nUMIs_per_cluster.png'), dpi=320, width=10, height= 10)
  
  # Boxplot of nGenes per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nFeature_RNA, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/nGenes_per_cluster.png'), dpi=320, width=10, height= 10)
  
  
  # Boxplot of percent.mt per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=percent.mt, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/percent.mt_per_cluster.png'), dpi=320, width=10, height= 10)
  
  
  # Boxplot of S.Score per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=S.Score, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/S.score_per_cluster.png'), dpi=320, width=10, height= 10)
  
  
  # Boxplot of G2M.Score per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=G2M.Score, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/G2M.score_per_cluster.png'), dpi=320, width=10, height= 10)
  
  
  # Extract identity and sample information from seurat object to determine 
  # the number of cells per cluster per sample
  n_cells = FetchData(seurat_umap, 
                       vars = c("ident", "sample")) %>%
    dplyr::count(ident, sample)
  
  # Barplot of number of cells per cluster by sample
  ggplot(n_cells, aes(x=ident, y=n, fill=sample))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_sample.png'), dpi=320, width=10, height= 10)
  
  
  # Extract identity and condition information from seurat object to determine 
  # the number of cells per cluster per condition
  n_cells = FetchData(seurat_umap, 
                      vars = c("ident", "condition")) %>%
    dplyr::count(ident, condition)
  
  # Barplot of number of cells per cluster by condition
  ggplot(n_cells, aes(x=ident, y=n, fill=condition))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_condition.png'), dpi=320, width=10, height= 10)
  
  # Barplot of proportion of cells in each cluster by sample
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=sample), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_sample.png'), dpi=320, width=10, height= 10)
  
  
  # Barplot of proportion of cells in each cluster by sample
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=condition), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_condition.png'), dpi=320, width=10, height= 10)
}


# clustering QC at SCT_snn_res.0.1 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering/sct/')
  
## Which PCs drive different clusters? ----

# determine max number of PCs present in seurat object
max_pcs = length(grep("^PC_", colnames(seurat_umap_sct@reductions$pca@cell.embeddings)))

# Create PC names based on the max PCs present
columns = c(paste0("PC_", 1:max_pcs),
            "ident",
            "umap_1", "umap_2")

## overlay PCs on umap ----
# takes a while to run so will use parallel processing
# library(future.apply)
# library(parallel)
# 
# options(future.globals.maxSize = 8 * 1024^3)  # Set limit to 8 GB

## function to overlay PCs on umap
# takes a while to run so will use parallel processing with 3 workers
pcs_umap_overlay = function(seurat_umap, resolution, directory) {
  
  # plan(list(sequential, tweak(multisession, workers=2)))
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # set resolution
  Idents(seurat_umap) = resolution
  
  # Extract columns from the seurat object. FetchData() can pull any data
  # from expression matrices, cell embeddings, or metadata
  pc_data = FetchData(seurat_umap, vars = columns)
  
  # Adding cluster label to center of cluster on UMAP
  umap_label = FetchData(seurat_umap, 
                          vars = c("ident", "umap_1", "umap_2"))  %>%
    group_by(ident) %>%
    dplyr::summarise(x=mean(umap_1), y=mean(umap_2))
  
  # Chunk PCs into blocks of 12
  pc_chunks = split(paste0("PC_", 1:max_pcs), ceiling(seq_along(1:max_pcs)/12))
  
  # Plotting a UMAP plot for each of the PCs
  lapply(seq_along(pc_chunks), function(chunk){
    
    pc_range = pc_chunks[[chunk]]
    first_pc = pc_range[1]
    last_pc = pc_range[length(pc_range)]
    
    lapply(pc_chunks[[chunk]], function(pc) {
      ggplot(pc_data, aes(umap_1, umap_2))+
        geom_point(aes(color=!!sym(pc)), alpha = 0.7)+
        scale_color_gradient(guide = 'none', low = "grey90", high = "blue")+
        geom_text(data=umap_label, aes(x=x, y=y, label=ident))+
        ggtitle(pc)+
        theme(plot.title = element_text(hjust = 0.5, size=24, face='bold'))+
        theme_classic2()
    }) %>% 
    plot_grid(plotlist = .)
    ggsave(paste0(directory,resolution,'/UMAP_by_PC_components_',first_pc,'-',last_pc,'.png'), dpi=320, width=20, height=15)
    
  })
  # on.exit(plan(sequential)) # restore to default non-parallel settings
}

## plot PCs on umap for resolution SCT_snn_res.0.1 ----
# takes a while to run
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering/sct/')


#### exploring cell marker genes ----
library(readxl)

# read in cell markers from excel file
cell_markers = read_excel('annotation/aggregate_sources_celltype_markers.xlsx')
cell_markers = lapply(cell_markers, na.omit) # convert to list of lists, omit NAs

# set RNA counts slot to be default and normalize for visualizing genes
# SCTransform normalization was performed only on the 3000 most variable genes, 
# so many of our genes of interest may not be present in the SCT data. Hence,
# we switch over to the original RNA count to explore cell type marker genes.
DefaultAssay(seurat_umap_sct) = 'RNA'
# seurat_umap_sct = NormalizeData(seurat_umap_sct, verbose = FALSE)
# DefaultAssay(seurat_umap_sct) = 'SCT'
# saveRDS(seurat_umap_sct, 'data/seurat_umap_sct.rds')

# Want to only consider cell markers that are actually present in the dataset
cell_markers = lapply(cell_markers, function(x){
  x = unique(x[x %in% rownames(seurat_umap_sct)])
})

# remove any cell types with zero markers in the dataset
cell_markers = Filter(length, cell_markers)
# Erythrocyte cell type removed since no gene markers present in this dataset

# process and save the gene markers actually present in the data into a new file
max_length = max(sapply(cell_markers, length)) # max length lists of gene markers

# Pad shorter lists with NA to match max length so all gene lists have same length
padded_markers = lapply(cell_markers, function(x) {
  length(x) = max_length
  return(x)
})

# save gene markers actually present in the dataset to csv
write.csv(as.data.frame(padded_markers), 
          file='annotation/aggregate_sources_celltype_markers_actually_in_data.csv',
          row.names=FALSE,
          na='')


# function to plot all gene markers on umap and dotplots for each cell type in batches of 16 genes ----
# takes a while to run
plot_gene_markers = function(seurat_umap, resolution, directory) {
  
  DefaultAssay(seurat_umap) = 'RNA'
  
  # set resolution
  Idents(seurat_umap) = resolution
  
  lapply(seq_along(cell_markers), function(i) {
    gene_chunks = split(cell_markers[[i]], ceiling(seq_along(cell_markers[[i]])/16)) # plot 16 gene markers at a time
    lapply(seq_along(gene_chunks), function(j) {
      FeaturePlot(seurat_umap, # overlay gene marker expression on umap
                  reduction='umap',
                  features=gene_chunks[[j]],
                  order=T,
                  min.cutoff='q10',
                  label=T)+
        plot_annotation(
          title=paste0(names(cell_markers[i]),'-',j),
          theme=theme(plot.title = element_text(hjust = 0.5,size=24,face='bold')))
      ggsave(paste0(directory,resolution,'/annotation/',names(cell_markers[i]),'-',j,'.png'), 
                    dpi=320, width=20, height=20)
      
      png(filename=paste0(directory,resolution,'/annotation/',names(cell_markers[i]),'-',j,'_dotplot.png'), 
          units='in', width=20, height=10, res=320)
      print(
        DotPlot(seurat_umap, features=gene_chunks[[j]])+ # construct dotplot of gene markers for each cell type
          RotatedAxis()+
          ggtitle(paste0(names(cell_markers[i]),'-',j))+
          theme(plot.title = element_text(hjust = 0.5, size=24, face='bold'))
      )
      graphics.off()
    })
  })
}
  
  
## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.1. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering/sct/')

# looks like we could benefit from increasing resolution
# resolutions available. For reference
resolutions = sprintf('SCT_snn_res.%.1f',c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 2.0))

# Assign identity of clusters to desired resolution
Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'

DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=TRUE,
        label.size=6)


#### Resolution SCT_snn_res.0.2 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.2 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering/sct/')

# clustering QC at SCT_snn_res.0.2 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering/sct/')

## plot PCs on umap for resolution SCT_snn_res.0.2 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.2. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering/sct/')


#### Resolution SCT_snn_res.0.4 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.4 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering/sct/')

# clustering QC at SCT_snn_res.0.4 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering/sct/')

## plot PCs on umap for resolution SCT_snn_res.0.4 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.4. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering/sct/')



#### Resolution SCT_snn_res.1.2 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.1.2 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering/sct/')

# clustering QC at SCT_snn_res.1.2 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering/sct/')

## plot PCs on umap for resolution SCT_snn_res.1.2 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.1.2 ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering/sct/')


#### Resolution SCT_snn_res.1.4 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.1.4 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering/sct/')

# clustering QC at SCT_snn_res.1.4 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering/sct/')

## plot PCs on umap for resolution SCT_snn_res.1.4 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.1.4 ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering/sct/')


#### annotate clusters ----
Idents(seurat_umap_sct) = 'SCT_snn_res.1.2' # set cluster #s as ident

DefaultAssay(seurat_umap_sct) = "RNA"

#### Find marker genes ----

## Find all positive markers
# join layers to be able to run FindAllMarkers()
combined_seurat = JoinLayers(seurat_umap_sct)

all.pos.markers = FindAllMarkers(combined_seurat,
                                 only.pos=TRUE,
                                 logfc.threshold = 0.25,
                                 min.pct=0.1,
                                 min.diff.pct = -Inf) 

# Combine markers with gene descriptions 
annotations = read.csv('annotation/annotationhub_mice_genes.csv')
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
write.csv(ann.pos.markers,'annotation/annotated_pos_markers.csv',row.names=FALSE)

# Extract top 10 markers per cluster
top10 = ann.pos.markers %>% 
  group_by(cluster) %>%
  slice_max(order_by=tibble(avg_log2FC,pct.diff,p_val_adj),n=10)


## Find all positive conserved markers between conditions ----
# Create function to get conserved markers for any given cluster
get_conserved = function(cluster){
  FindConservedMarkers(combined_seurat, # need to specify exact seurat object here
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

# get conserved markers
conserved.pos.markers = map_dfr(c(0:30), get_conserved) # 31 clusters at resolution 1.4

# add pct.diff col
conserved.pos.markers = conserved.pos.markers %>%
  mutate(avg_fc = (ctrl_avg_log2FC + pb_avg_log2FC) / 2,
         max_adj_pval = ifelse(ctrl_p_val_adj > pb_p_val_adj,ctrl_p_val_adj,pb_p_val_adj),
         avg.pct.diff = ((pb_pct.1 - pb_pct.2) + (ctrl_pct.1 - ctrl_pct.2)) / 2)

# Rearrange the columns to be more intuitive
conserved.pos.markers = conserved.pos.markers[ , c(1:2,16:18,15,3:14)]

# Order the rows by p-adjusted values
conserved.pos.markers = conserved.pos.markers %>%
  dplyr::arrange(cluster_id, max_adj_pval)

# save markers
write.csv(conserved.pos.markers,'annotation/conserved_pos_markers.csv',row.names=FALSE)


# Extract top 10 conserved markers per cluster
top10.conserved = conserved.pos.markers %>% 
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)
top10.conserved = top10.conserved[!(top10.conserved$cluster_id == 29 | 
                                           top10.conserved$cluster_id == 30),]

# clusters 29 and 30 are made up of one condition. 
# Will treat individual values as avg and max to easily bind rows to top10.conserved
topc29 = conserved.pos.markers[conserved.pos.markers$cluster_id == 29,] %>%
  mutate(avg_fc = pb_avg_log2FC,
         max_adj_pval = pb_p_val_adj,
         avg.pct.diff = pb_pct.1 - pb_pct.2) %>%
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)

topc30 = conserved.pos.markers[conserved.pos.markers$cluster_id == 30,] %>%
  mutate(avg_fc = ctrl_avg_log2FC,
         max_adj_pval = ctrl_p_val_adj,
         avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)

top10.conserved = rbind(top10.conserved, topc29, topc30)


# Rename cluster #'s to cell type annotations ----
seurat_umap_sct = RenameIdents(seurat_umap_sct,
                                 '0' = 'Naive CD4+ T cells',
                                 '1' = 'CD8+ NKT-like cells',
                                 '2' = 'Naive B cells',
                                 '3' = 'Naive B cells',
                                 '4' = 'Naive B cells',
                                 '5' = 'Naive B cells',
                                 '6' = 'Gamma Delta T cells',
                                 '7' = 'Vascular.endothelial',
                                 '8' = 'Effector CD8+ T cells',
                                 '9' = 'Adipocyte',
                                 '10' = 'Neutrophils?',
                                 '11' = 'Fibroblast',
                                 '12' = 'Immature B cells',
                                 '13' = 'Effector CD4+ T cells',
                                 '14' = 'Macrophage.Ma',
                                 '15' = 'T helper cells',
                                 '16' = 'Unkonwn1',
                                 '17' = 'Unknown2',
                                 '18' = 'T-B reticular cells',
                                 '19' = 'Naive B cells',
                                 '20' = 'Endothelial',
                                 '21' = 'Macrophage.Mb',
                                 '22' = 'Adipocyte',
                                 '23' = 'Myofibroblasts/Pericyte',
                                 '24' = 'Macrophage.Mb',
                                 '25' = 'Macrophage.Mb',
                                 '26' = 'Fibroblast',
                                 '27' = 'Luminal.HS',
                                 '28' = 'Unknown3',
                                 '29' = 'Naive CD8+ T cells',
                                 '30' = 'Myocyte')

# assign cluster names to meta data
seurat_umap_sct$seurat_clusters = Idents(seurat_umap_sct)
saveRDS(seurat_umap_sct, 'data/annotated_seurat.rds')