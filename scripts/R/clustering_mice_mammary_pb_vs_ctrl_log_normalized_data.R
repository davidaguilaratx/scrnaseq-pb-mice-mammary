setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(SeuratObject)
library(scCustomize)
library(sctransform)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(ggplot2)
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


### continue analysis without integration using log normalization ----

# project scores of each gene in the dataset (including genes not included in the PCA) based on their correlation with the calculated components.
seurat_phase = ProjectDim(object=seurat_phase, reduction='pca') 

# visualize PC loadings
VizDimLoadings(seurat_phase, dims = 1:12, ncol=4)
ggsave('figures/clustering/lognorm/PCA_loadings.png', dpi=320, width=20, height= 15)


### expression of most highly weighted genes per PC ----
png(filename = 'figures/clustering/lognorm/high_score_genes_per_PC1-6_heatmap.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_phase, 
           reduction = "pca", 
           dims = 1:6,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/lognorm/high_score_genes_per_PC7-12_heatmap.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_phase, 
           reduction = "pca", 
           dims = 7:12,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/lognorm/high_score_genes_per_PC13-18_heatmap.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_phase, 
           reduction = "pca", 
           dims = 13:18,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()


png(filename = 'figures/clustering/lognorm/high_score_genes_per_PC19-24_heatmap.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_phase, 
           reduction = "pca", 
           dims = 19:24,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering/lognorm/high_score_genes_per_PC25-30_heatmap.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_phase, 
           reduction = "pca", 
           dims = 25:30,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

# ranked list (PC score) of genes for PC1
names(sort(Loadings(seurat_phase, reduction="pca")[,1], decreasing=TRUE))


# Printing out the top 10 most variable genes driving PCs
print(x = seurat_phase[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

### Plot the elbow plot to see variation captures by PCs ----
ElbowPlot(object = seurat_phase, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering/lognorm/elbow_plot.png', dpi=320, width=20, height= 10)

# We can calculate where the principal components start to elbow by taking the smaller value of:
#
# 1)The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# 2)The point where the percent change in variation between the consecutive PCs is less than 0.1%.

# # Determine percent of variation associated with each PC
pct = seurat_phase[["pca"]]@stdev / sum(seurat_phase[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu > 90 & pct < 5)[1]
co1
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs
# [1] 21
pcs = 21

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
ggsave('figures/clustering/lognorm/quantitative_elbow_plot.png', dpi=320, width=20, height= 10)


#  Jackstraw analysis to estimate a p-value for each component and plot
# seurat_phase = JackStraw(object = seurat_phase, num.replicate = 100, dims=30); # takes around 30 minutes
# seurat_phase = ScoreJackStraw(object = seurat_phase, dims = 1:30)
# JackStrawPlot(seurat_phase, dims=1:30)

### cluster log-normalized non-integrated cells ----
seurat_umap = FindNeighbors(seurat_phase, dims=1:pcs)
seurat_umap = FindClusters(object = seurat_umap,
                           resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 2.0))

# Resolution 0.1: 12 clusters
# Resolution 0.2: 15 clusters
# Resolution 0.4: 18 clusters
# Resolution 0.6: 20 clusters
# Resolution 0.8: 24 clusters
# Resolution 1.0: 26 clusters
# Resolution 1.4: 33 clusters
# Resolution 2.0: 41 clusters

saveRDS(seurat_umap, 'data/seurat_umap.rds')



# resolutions
resolutions = sprintf('RNA_snn_res.%.1f',c(0.1, 0.4, 0.6, 0.8, 1.0, 1.4, 2.0))

# Assign identity of clusters to resolution 0.1
Idents(seurat_umap) = 'RNA_snn_res.0.1'

# run UMAP on non-intergated, merged, log-normalized data ----
seurat_umap = RunUMAP(seurat_umap,
                      dims=1:pcs, # pcs=21
                      reduction='pca') 

# saveRDS(seurat_umap, 'data/seurat_umap_umap_res0.1')

# Plot the UMAP of non-integrated data
DimPlot(seurat_umap,
        reduction = "umap")


# plot UMAP by sample and condition
p1 = DimPlot(seurat_umap,
             reduction = "umap",
             group.by='sample',
             shuffle=T)
p2 = DimPlot(seurat_umap,
             reduction = "umap",
             group.by='condition',
             shuffle=T)

plot_grid(p1, p2)

# plot UMAP by Phase and mitochondrial expression quartile
p1 = DimPlot(seurat_umap,
             reduction = "umap",
             group.by='Phase',
             shuffle=T)
p2 = DimPlot(seurat_umap,
             reduction = "umap",
             group.by='mitoFr',
             shuffle=T)

plot_grid(p1, p2)




marker_list = c("Trf", "Cited1", "Krt14", "Rgs5", "Csf1r", "Cd4", "Mmp12","Prox1","Epcam","Cdh5","Col1a1")

marker_list = c('Acta2','CD44','CK14','CK5','Etv5','Krt14') # basal

marker_list = c('Acta2','Oxtr','Krt4','Krt14') # myoepithelial

marker_list = c('Aldh1a3','c-Kit','Cd14','Cited1','Esr1','Pgr','Prlr','S100a6') # hormoning sensing proginitor

marker_list = c('CD29','CD24','CD49f','Sox2','Oct4','Sca-1','Nanog','CD44') # stem cells

marker_list = 'Pdgfra'


FeaturePlot(seurat_umap, features = marker_list) ###Hmmm... no epithelial cells


DotPlot(seurat_umap, features = marker_list) + RotatedAxis()



















# 
# lapply(resolutions, funciont(res) {
#   Idents(seurat_umap) = res # Assign identity of clusters to resolution res
#   seurat_umap = RunUMAP(seurat_umap,
#                         dims=1:pcs, # pcs=21
#                         reduction='pca') 
#   # Plot the UMAP of non-integrated data
#   DimPlot(seurat_umap,
#           reduction = "umap")
#   
#   # plot UMAP by sample and condition
#   p1 = DimPlot(seurat_umap,
#                reduction = "umap",
#                group.by='sample',
#                shuffle=T)
#   p2 = DimPlot(seurat_umap,
#                reduction = "umap",
#                group.by='condition',
#                shuffle=T)
#   
#   plot_grid(p1, p2)
#   
#   
#   # plot UMAP by Phase and mitochondrial expression quartile
#   p1 = DimPlot(seurat_umap,
#                reduction = "umap",
#                group.by='Phase',
#                shuffle=T)
#   p2 = DimPlot(seurat_umap,
#                reduction = "umap",
#                group.by='mitoFr',
#                shuffle=T)
#   
#   plot_grid(p1, p2)
#   
# })