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
library(Polychrome)
library(pheatmap)
library(presto)
library(qs2) # for fast data saving/loading


# set seed
set.seed(2024)

# sessionInfo()

# read in data ----
# Set timepoint to analyze for this script
selected_timepoint <- 'week3'
selected_timepoint <- 'week10'
selected_timepoint <- 'month7'
selected_timepoint <- 'month18'


# directories for saving
timepoint_dirs <- list(
  "week3" = list(
    fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/qc/",
    data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/",
    annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
  ),
  "week10" = list(
    fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/qc/",
    data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/data/",
    annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
  ),
  "month7" = list(
    fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/qc/",
    data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/data/",
    annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_month/cellbender_analysis/FPR_0.0_MAD/annotation/"
  ),
  "month18" = list(
    fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/qc/",
    data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/data/",
    annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/annotation/"
  )
)
list2env(timepoint_dirs[[selected_timepoint]], envir = environment())

# Create directories to save results if they don't already exist:
dirs = c(annot_dir, data_dir, fig_dir)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}


seurat_umap_sct = qs_read(paste0(data_dir,''))










# clusters 26 is made up of one condition.
# Will treat individual values as avg and max to easily bind rows to top10.conserved
topc26 = conserved.pos.markers[conserved.pos.markers$cluster_id == 26,] %>%
  mutate(avg_fc = ctrl_avg_log2FC,
         max_adj_pval = ctrl_p_val_adj,
         avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)

top10.conserved = rbind(top10.conserved, topc26)

## test subclustering t cells ----
test = subset(DietSeurat(seurat_umap_sct, assays='RNA'), subset=(seurat_clusters == 1) | 
                (seurat_clusters == 2) | (seurat_clusters == 4) | 
                (seurat_clusters == 10) | (seurat_clusters == 13) | 
                (seurat_clusters == 20))

test = SplitObject(test, split.by='sample')

test = lapply(test, function(x) {
  x = SCTransform(x, assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                  vars.to.regress=c('cc.difference','percent.mt'),
                  return.only.var.genes=FALSE)
})

features = SelectIntegrationFeatures(test, nfeatures=5000, fvf.nfeatures=5000)

# merge sctransform normalized samples. Takes a little while to run
test = merge(test[[1]], test[-1], merge.data=TRUE,merge.dr=FALSE)

VariableFeatures(test) = features

sct_results = SCTResults(test, slot = "feature.attributes")
qs_save(sct_results, 'figures/clustering_rm29and32_regress_ccdiffandmt/sct_results_test.qs',nthreads=2)
rm(sct_results)

# Variable features lost/unset after merging transformed objects. Need to set again.
# should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(test[["SCT"]]) = rownames(test[["SCT"]]@scale.data)

test = RunPCA(test, assay='SCT', npcs=50)
test = FindNeighbors(test, dims=1:50)
test = FindClusters(object=test, resolution=c(0.1,0.4,0.8,1.2), graph.name='SCT_snn')
test = RunUMAP(test, dims=1:50, reduction='pca')

Idents(test) = 'SCT_snn_res.0.4'
DimPlot(test,
        reduction='umap',
        label=TRUE,
        shuffle=TRUE)

table(test$condition, test$SCT_snn_res.0.4)




# subclustering cluster 19 and 25, epithelial cells ----
# normalization
epithelial = subset(DietSeurat(seurat_umap_sct, assays='RNA'), subset=(seurat_clusters == 19) | 
                      (seurat_clusters == 25))

epithelial = SplitObject(epithelial, split.by='sample')

epithelial = lapply(epithelial, function(x) {
  x = SCTransform(x, assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                  vars.to.regress=c('cc.difference','percent.mt'),
                  return.only.var.genes=FALSE)
})

features = SelectIntegrationFeatures(epithelial, nfeatures=5000, fvf.nfeatures=5000)

# merge sctransform normalized samples. Takes a little while to run
epithelial = merge(epithelial[[1]], epithelial[-1], merge.data=TRUE,merge.dr=FALSE)

VariableFeatures(epithelial) = features

epithelial_sct_results = SCTResults(epithelial, slot = "feature.attributes")
qs_save(epithelial_sct_results, 'figures/clustering_rm29and32_regress_ccdiffandmt/sct_results_epithelial_epithelial.qs',nthreads=2)
rm(sct_results)

# # Variable features lost/unset after merging transformed objects. Need to set again.
# # should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(epithelial[["SCT"]]) = rownames(epithelial[["SCT"]]@scale.data)

epithelial = RunPCA(epithelial, npcs=50)

ElbowPlot(object = epithelial, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))

### Determine percent of variation associated with each PC
pct = epithelial[["pca"]]@stdev / sum(epithelial[["pca"]]@stdev) * 100
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


epithelial = FindNeighbors(epithelial, dims=1:50)
epithelial = FindClusters(object=epithelial, resolution=c(0.1,0.4,0.8,1.0,1.2))
epithelial = RunUMAP(epithelial, dims=1:50, reduction='pca')

qs_save(epithelial,paste0(data_dir,'epithelial_19&25_umap_sct.qs'),nthreads=2)

Idents(epithelial) = 'RNA_snn_res.0.4'
DimPlot(epithelial,
        reduction='umap',
        label=TRUE,
        shuffle=TRUE)

table(epithelial$condition, epithelial$RNA_snn_res.0.4)
# 0   1   2   3   4
# ctrl 182 174 117  37   2
# pb   114  77  71  37   1


Idents(epithelial) = 'SCT_snn_res.0.8'
DimPlot(epithelial,
        reduction='umap',
        label=TRUE,
        shuffle=TRUE)
table(epithelial$condition, epithelial$RNA_snn_res.0.4)
# 0   1   2
# ctrl 174 119  52
# pb   149  94  65

DimPlot(epithelial,
        reduction='umap',
        group.by = 'condition',
        label=FALSE,
        shuffle=TRUE)

table(epithelial$condition, epithelial$SCT_snn_res.1.2)
# 0   1   2   3   4
# ctrl 182 174 117  37   2
# pb   114  77  71  37   1

DoHeatmap(epithelial, features[1:100],size=4,angle=90)+NoLegend()

## plot umap and umap qcs for epithelial cells ----
if (!dir.exists(paste0(fig_dir,'epithelial/SCT_snn_res.0.8/'))) { dir.create(paste0(fig_dir,'epithelial/SCT_snn_res.0.8/'),
                                                                             recursive=TRUE) }
plot_umaps(epithelial, resolution='SCT_snn_res.0.8', directory=paste0(fig_dir,'epithelial/'))
umap_clustering_qc(epithelial, resolution='SCT_snn_res.0.8', directory=paste0(fig_dir,'epithelial/'))


## myoepithelial  cells ----
# cluster 2
Idents(epithelial) = 'SCT_snn_res.0.8' # set cluster #s as ident

DefaultAssay(epithelial) = "RNA"

epithelial.markers = FindMarkers(epithelial,
                                 ident.1 = 2,
                                 ident.2 = c(0:1)) 

# Add gene symbols to the DE table
epithelial.markers = epithelial.markers %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
epithelial.markers$pct.diff = epithelia.markers$pct.1 - epithelia.markers$pct.2

# Reorder columns and sort by padj   
epithelia.markers = epithelia.markers[, c(1, 3:5,2,8,6:7)]

epithelia.markers = epithelia.markers %>%
  dplyr::arrange(p_val_adj)

write.csv(epithelia.markers,'annotation/epithelia.markers_cluster0.csv',row.names=FALSE)

FeaturePlot(epithelial, features = c('Esr1','Prlr','Fgfr2','Krt14','Krt5','Trp63','Krt17'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

epithelial.combined = JoinLayers(epithelial)
epithelial.markers = FindAllMarkers(epithelial.combined,
                                    only.pos=TRUE,
                                    logfc.threshold = 0.25,
                                    min.pct=0.1,
                                    min.diff.pct = -Inf) 

# Combine markers with gene descriptions 
annotations = read.csv('C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/annotation/annotationhub_mice_genes.csv')
epithelial.ann.markers = left_join(x = epithelial.markers, 
                                   y = annotations[, c("gene_name", "description")],
                                   by = c("gene" = "gene_name")) %>%
  unique()

# add pct.diff col
epithelial.ann.markers$pct.diff = epithelial.ann.markers$pct.1 - epithelial.ann.markers$pct.2

# Rearrange the columns to be more intuitive
epithelial.ann.markers = epithelial.ann.markers[ , c(6, 7, 2:4, 9, 1, 5,8)]

# Order the rows by p-adjusted values
epithelial.ann.markers = epithelial.ann.markers %>%
  dplyr::arrange(cluster, p_val_adj)

# save markers
write.csv(epithelial.ann.markers,paste0(annot_dir,'epithelial_annotated_pos_markers_res.0.8.csv'),row.names=FALSE)

# Extract top 10 markers per cluster
epithelial_top10 = epithelial.ann.markers %>% 
  group_by(cluster) %>%
  slice_max(order_by=tibble(avg_log2FC,pct.diff,p_val_adj),n=10)

# get conserved markers
epithelial_get_conserved = function(cluster){
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

epithelial.conserved.markers = map_dfr(c(0:2), epithelial_get_conserved) # 3 clusters at resolution 0.8

# add avg.pct.diff, max_adj_pval, and min.pct 1 cols
epithelial.conserved.markers = epithelial.conserved.markers %>%
  mutate(avg_log2FC = (ctrl_avg_log2FC + pb_avg_log2FC) / 2,
         max_adj_pval = ifelse(ctrl_p_val_adj > pb_p_val_adj,ctrl_p_val_adj,pb_p_val_adj),
         avg.pct.diff = ((pb_pct.1 - pb_pct.2) + (ctrl_pct.1 - ctrl_pct.2)) / 2,
         min.pct.1 = ifelse(ctrl_pct.1 < pb_pct.1, ctrl_pct.1,pb_pct.1))

# Rearrange the columns to be more intuitive
epithelial.conserved.markers = epithelial.conserved.markers[ , c(1:2,16:19,15,3:14)]

# Order the rows by p-adjusted values
epithelial.conserved.markers = epithelial.conserved.markers %>%
  dplyr::arrange(cluster_id, max_adj_pval)

# save markers
write.csv(conserved.pos.markers, paste0(annot_dir,'epithelial_conserved_pos_markers_res.0.8.csv'),row.names=FALSE)


# Extract top 10 conserved markers per cluster
epithelial.top10.conserved = epithelial.conserved.markers %>% 
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_log2FC,avg.pct.diff,max_adj_pval),n=10)






# FindMarkers() for clusters still in question ----
## Cd4 T cells
# cluster 4 vs 1 and 10
cd4_tcells = FindMarkers(seurat_umap_sct,
                         ident.1 = 4,
                         ident.2 = c(1,10)) 

# Add gene symbols to the DE table
cd4_tcells = cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
cd4_tcells$pct.diff = cd4_tcells$pct.1 - cd4_tcells$pct.2

# Reorder columns and sort by padj   
cd4_tcells = cd4_tcells[, c(1, 3:5,2,8,6:7)]

cd4_tcells = cd4_tcells %>%
  dplyr::arrange(p_val_adj)

write.csv(cd4_tcells,'annotation/cd4_tcells_4vs1and10.csv',row.names=FALSE)

# cluster 10 vs 1 and 4
cd4_tcells = FindMarkers(seurat_umap_sct,
                         ident.1 = 10,
                         ident.2 = c(1,4)) 

# Add gene symbols to the DE table
cd4_tcells = cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
cd4_tcells$pct.diff = cd4_tcells$pct.1 - cd4_tcells$pct.2

# Reorder columns and sort by padj   
cd4_tcells = cd4_tcells[, c(1, 3:5,2,8,6:7)]

cd4_tcells = cd4_tcells %>%
  dplyr::arrange(p_val_adj)

write.csv(cd4_tcells,'annotation/cd4_tcells_10vs1and4.csv',row.names=FALSE)

#B cells
# cluster 3 vs 0,11,and 16
bcells = FindMarkers(seurat_umap_sct,
                     ident.1 = 3,
                     ident.2 = c(0,11,16)) 

# Add gene symbols to the DE table
bcells = bcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
bcells$pct.diff = bcells$pct.1 - bcells$pct.2

# Reorder columns and sort by padj   
bcells = bcells[, c(1,3:5,2,8,6:7)]

bcells = bcells %>%
  dplyr::arrange(p_val_adj)

write.csv(bcells,'annotation/bcells_cluster3.csv',row.names=FALSE)


# cluster 11 vs 0,3,and 16
bcells = FindMarkers(seurat_umap_sct,
                     ident.1 = 11,
                     ident.2 = c(0,3,16)) 

# Add gene symbols to the DE table
bcells = bcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
bcells$pct.diff = bcells$pct.1 - bcells$pct.2

# Reorder columns and sort by padj   
bcells = bcells[, c(1,3:5,2,8,6:7)]

bcells = bcells %>%
  dplyr::arrange(p_val_adj)

write.csv(bcells,'annotation/bcells_cluster11.csv',row.names=FALSE)


# cluster 16 vs 0,3,and 11
bcells = FindMarkers(seurat_umap_sct,
                     ident.1 = 16,
                     ident.2 = c(0,3,11)) 

# Add gene symbols to the DE table
bcells = bcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
bcells$pct.diff = bcells$pct.1 - bcells$pct.2

# Reorder columns and sort by padj   
bcells = bcells[, c(1,3:5,2,8,6:7)]

bcells = bcells %>%
  dplyr::arrange(p_val_adj)

write.csv(bcells,'annotation/bcells_cluster16.csv',row.names=FALSE)



## adipocytes
# cluster 21 vs 7
adipocytes = FindMarkers(seurat_umap_sct,
                         ident.1 = 21,
                         ident.2 = c(7))

# Add gene symbols to the DE table
adipocytes = adipocytes %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
adipocytes$pct.diff = adipocytes$pct.1 - adipocytes$pct.2

# Reorder columns and sort by padj   
adipocytes = adipocytes[, c(1,3:5,2,8,6:7)]

adipocytes = adipocytes %>%
  dplyr::arrange(p_val_adj)

write.csv(adipocytes,'annotation/adipocytes_cluster21.csv',row.names=FALSE)


# cluster 7 vs 21
adipocytes = FindMarkers(seurat_umap_sct,
                         ident.1 = 7,
                         ident.2 = 21)

# Add gene symbols to the DE table
adipocytes = adipocytes %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
adipocytes$pct.diff = adipocytes$pct.1 - adipocytes$pct.2

# Reorder columns and sort by padj   
adipocytes = adipocytes[, c(1,3:5,2,8,6:7)]

adipocytes = adipocytes %>%
  dplyr::arrange(p_val_adj)

write.csv(adipocytes,'annotation/adipocytes_cluster7.csv',row.names=FALSE)


## fibroblasts
# clusters 17 vs 6 and 24
fibroblasts = FindMarkers(seurat_umap_sct,
                          ident.1 = 17,
                          ident.2 = c(6,24))

# Add gene symbols to the DE table
fibroblasts = fibroblasts %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
fibroblasts$pct.diff = fibroblasts$pct.1 - fibroblasts$pct.2

# Reorder columns and sort by padj   
fibroblasts = fibroblasts[, c(1,3:5,2,8,6:7)]

fibroblasts = fibroblasts %>%
  dplyr::arrange(p_val_adj)

write.csv(fibroblasts,'annotation/fibroblasts_cluster17.csv',row.names=FALSE)

# cluster 24 vs 6 and 17
fibroblasts = FindMarkers(seurat_umap_sct,
                          ident.1 = 24,
                          ident.2 = c(6,17))

# Add gene symbols to the DE table
fibroblasts = fibroblasts %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
fibroblasts$pct.diff = fibroblasts$pct.1 - fibroblasts$pct.2

# Reorder columns and sort by padj   
fibroblasts = fibroblasts[, c(1,3:5,2,8,6:7)]

fibroblasts = fibroblasts %>%
  dplyr::arrange(p_val_adj)

write.csv(fibroblasts,'annotation/schwann_cells_cluster24.csv',row.names=FALSE)


## epithelial
# clusters 23 vs 18 and 24
epithelial = FindMarkers(seurat_umap_sct,
                         ident.1 = 23,
                         ident.2 = c(18,24))

# Add gene symbols to the DE table
epithelial = epithelial %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
epithelial$pct.diff = epithelial$pct.1 - epithelial$pct.2

# Reorder columns and sort by padj   
epithelial = epithelial[, c(1,3:5,2,8,6:7)]

epithelial = epithelial %>%
  dplyr::arrange(p_val_adj)

write.csv(epithelial,'annotation/epithelial_cluster23.csv',row.names=FALSE)


# clusters 18 vs 23 and 24
epithelial = FindMarkers(seurat_umap_sct,
                         ident.1 = 18,
                         ident.2 = c(23,24))

# Add gene symbols to the DE table
epithelial = epithelial %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
epithelial$pct.diff = epithelial$pct.1 - epithelial$pct.2

# Reorder columns and sort by padj   
epithelial = epithelial[, c(1,3:5,2,8,6:7)]

epithelial = epithelial %>%
  dplyr::arrange(p_val_adj)

write.csv(epithelial,'annotation/epithelial_cluster18.csv',row.names=FALSE)



## endothelial

# myeloid

rownames(seurat_umap_sct)


# Rename cluster #'s to cell type annotations ----
Idents(seurat_umap_sct) = 'SCT_snn_res.0.8'
seurat_umap_sct = RenameIdents(seurat_umap_sct,
                               '0' = 'B cells',
                               '1' = 'CD4+ T cells',
                               '2' = 'CD8+ T cells',
                               '3' = 'Tregs',
                               '4' = 'B cells',
                               '5' = 'Endothelial',
                               '6' = 'B cells', # refined assessment says T cells. Naive/quiescent or memory. Refined uses top 10 markers from all and conserved
                               '7' = 'Fibroblast',
                               '8' = 'Proliferating cells',
                               '9' = 'Dendritic cell', # or DC-like
                               '10' = 'Adipocytes', # possibly white adipose
                               '11' = 'Activated CD4+ T Cells', # or immediate-early T cell cluster
                               '12' = 'Macrophage.Ma', # m2-skewed or tissue resident type
                               '13' = 'T helper 17 cells', #nk-like. an “innate-like” phenotype (sometimes called “NK‐like γδ T cells”)
                               '14' = 'Myeloid.Mb',
                               '15' = 'Myeloid.Mb', # refined says M2-like macrophages
                               '16' = 'Reticular cells', # reticular cells? FRC-like. Specialized stroma or transitional cell states
                               '17' = 'Endothelial', # refined says lymphatic associated. or transitional between blood and lymphatic. possibly high endothelial venule‐like or sinusoidal)
                               '18' = 'CD8+ T cells', # 
                               '19' = 'Basal-Myoepithelial', # or vascular smooth muscle. Subclustering shows small part is myoepithelial. large could be pericytes
                               '20' = 'B cells',
                               '21' = 'Adipocytes', # possibly beige/brown adipose tissue
                               '22' = 'Myeloid.Mb', # potentially inflammatory or specialized
                               '23' = 'Schwann cells', # or glial-like ells
                               '24' = 'Myeloid.Mb', #cDC1 lineage in mice. 
                               '25' = 'Luminal.HS', # hormone responsive
                               '26' = 'Muscle', # skeletal muscle
                               '27' = 'Erythroid cells'# mature RBCs
)

# seurat_umap_sct = RenameIdents(seurat_umap_sct,
#                                  '0' = 'B cells',
#                                  '1' = 'CD4+ T cells',
#                                  '2' = 'CD8+ T cells',
#                                  '3' = 'B cells',
#                                  '4' = 'Tregs',
#                                  '5' = 'Endothelial',
#                                  '6' = 'Fibroblast',
#                                  '7' = 'Adipocytes',
#                                  '8' = 'Proliferating cells',
#                                  '9' = 'Dendritic cell',
#                                  '10' = 'Activated CD4+ T cells',
#                                  '11' = 'B cells',
#                                  '12' = 'Macrophage.Ma',
#                                  '13' = 'T helper 17 cells',
#                                  '14' = 'Myeloid.Mb',
#                                  '15' = 'Myeloid.Mb',
#                                  '16' = 'B cells',
#                                  '17' = 'Reticular cells',
#                                  '18' = 'Basal-Myoepithelial',
#                                  '19' = 'Endothelial',
#                                  '20' = 'CD8+ T cells',
#                                  '21' = 'Adipocytes',
#                                  '22' = 'Myeloid.Mb',
#                                  '23' = 'Luminal.HS',
#                                  '24' = 'Schwann cells',
#                                  '25' = 'Endothelial',
#                                  '26' = 'Muscle',
#                                  '27' = 'Erythroid cells')

# assign cluster names to meta data
seurat_umap_sct$celltype = Idents(seurat_umap_sct)
Idents(seurat_umap_sct) = 'celltype'

if (!dir.exists(paste0(annot_dir,'figures/'))) { dir.create(paste0(annot_dir,'figures/'),
                                                            recursive=TRUE) }

DimPlot(seurat_umap_sct, 
        group.by = 'celltype',
        label=T,
        label.size=6,
        repel=T,
        shuffle=T)
ggsave(paste0(annot_dir,'figures/annotated_seurat_res.0.8.png'), dpi=320, width=20, height=20)

qs_save(seurat_umap_sct, paste0(data_dir,'annotated_seurat_res.0.8.qs'),nthreads=2)

# plot qc metrics for cell-annotated clusters
if (!dir.exists(paste0(fig_dir,'celltype/'))) { dir.create(paste0(fig_dir,'celltype/'),
                                                           recursive=TRUE) }

plot_umaps(seurat_umap_sct, 
           resolution='celltype',
           directory=fig_dir)

umap_clustering_qc(seurat_umap_sct, 
                   resolution='celltype',
                   directory=fig_dir)

# shrink seurat object in preparation of DE analysis
de_seurat = DietSeurat(seurat_umap_sct, assays='RNA', layers=c('counts'))
qs_save(de_seurat, paste0(data_dir,'de_seurat_annotated.qs'),nthreads=2)

