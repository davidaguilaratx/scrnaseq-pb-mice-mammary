library(Seurat)
library(SeuratObject)
library(sctransform)
library(harmony)
# library(patchwork)
# library(cowplot)
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
# library(pheatmap)
library(presto)
library(qs2) # for fast data saving/loading
library(zeallot) # for unpacking multiple values from functions
library(openai)
library(GPTCelltype)
library(CyteTypeR)


source(here::here("scripts/R/utils.R"))
source(here::here("scripts/R/plot_utils.R"))

# set seed
seed <- 2024
set.seed(seed)

c(nthreads, bp_param) %<-% configure_parallelism(rng_seed = seed)

# read in data ----
# directories for saving
fig_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis_old/integrated/cellbender_analysis/FPR_0.0_MAD/clustering"
data_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis_old/integrated/cellbender_analysis/FPR_0.0_MAD/data"
annot_dir = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis_old/integrated/cellbender_analysis/FPR_0.0_MAD/annotation"

# Create directories to save results if they don't already exist:
dirs_to_make = c(
  annot_dir, 
  data_dir, 
  fig_dir)

walk(dirs_to_make, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))



integrated_seurat = qs_read(file.path(data_dir,'integrated_seurat_annotated.qs'), nthreads=nthreads)





## test subclustering t cells ----
test = subset(DietSeurat(integrated_seurat, assays='RNA'), subset=(seurat_clusters == 1) | 
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
epithelial = subset(DietSeurat(integrated_seurat, assays='RNA'), subset=(seurat_clusters == 19) | 
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

# qs_save(epithelial,paste0(data_dir,'epithelial_19&25_umap_sct.qs'),nthreads=2)

epithelial = qs_read(file.path(data_dir,'epithelial.qs'), nthreads=2)

DimPlot(epithelial,
        reduction='tsne_harmony',
        group.by='subtype',
        label=TRUE,
        shuffle=TRUE)

DimPlot(subset(epithelial,subset=(timepoint=='wk10')&(condition=='pb')),
        reduction='umap_harmony',
        group.by='subtype',
        split.by = 'sample',
        label=TRUE,
        shuffle=TRUE)

FeaturePlot(epithelial,
        features = c('Muc1'),
        # split.by ='condition',
        split.by ='timepoint',
        slot='scale.data',
        # slot='data',
        reduction='umap_harmony',
        min.cutoff='q10',
        max.cutoff='q50',
        cols=viridis(256)) &
  DarkTheme()


# Check presence of mammary stem / progenitor / bipotent markers in epithelial ----
mammary_markers = list(
  stem_basal = c(
    'Krt5','Krt14','Krt17','Trp63','Acta2','Myh11','Cnn1','Oxtr',
    'Cd44','Lgr5','Lgr6','Procr','Bmi1','Itga6','Itgb1','Cd24a',
    'Sox9','Nrg1','Lrp5','Axin2','Wnt5a','Wnt10a','Tspan8','Snai2',
    'Vim','Zeb1','Zeb2','Twist1','Twist2', 'Id1', 'Sox2', 'Aldha1a3',
    'Aldh1a1','Spy1'
  ),
  basal_myoepithelial = c(
    'Acta2','Myh11','Mylk','Oxtr','Krt5','Krt14','Krt17','Trp63','Tagln',
    'Cnn1','Sparc','Myl9'
  ),
  luminal_progenitor = c(
    'Kit','Aldh1a3','Elf5','Krt8','Krt18','Krt19','Itgb3','Cd14',
    'Hey1','Notch1','Notch2','Sox10','Tspan8','Cd61','Csn2','Wap',
    'Lalba','Foxc1','Foxm1','Mki67','Top2a'
  ),
  mature_luminal.hs = c(
    'Esr1','Pgr','Foxa1','Prom1','Areg','Prlr','Gata3','Tbx3',
    'Cited1','Mfge8','Ltf','Muc1','Sca1','Wnt4','Tnfsf11'
  ),
  bipotent = c(
    'Procr','Bmi1','Lgr5','Axin2','Sox9','Lrp5','Tspan8','Cd24a',
    'Itga6','Itgb1','Krt5','Krt8','Krt14','Krt18','Ybx1','Eno1'
  )
)

# Use default assay (likely SCT after subclustering)
all_genes = rownames(epithelial)
cat('Default assay:', DefaultAssay(epithelial),
    '| n genes:', length(all_genes), '\n\n')

# Exact-match check per category
marker_presence = imap(mammary_markers, function(genes, category) {
  present = intersect(genes, all_genes)
  missing = setdiff(genes, all_genes)
  cat('---', category, '---\n')
  cat('  present (', length(present), '/', length(genes), '):',
      paste(present, collapse=', '), '\n')
  cat('  missing:', paste(missing, collapse=', '), '\n\n')
  tibble(category=category, gene=genes, present=genes %in% all_genes)
}) %>% bind_rows()

# Fuzzy grep — catch case differences, family members (e.g. Krt5a, Wnt5b)
all_markers_unique = unique(unlist(mammary_markers))
grep_hits = map(all_markers_unique, function(g) {
  hits = grep(paste0('^', g, '$|^', g, '[0-9a-z]*$'),
              all_genes, ignore.case=TRUE, value=TRUE)
  if (length(hits)) tibble(query=g, hit=hits) else NULL
}) %>% compact() %>% bind_rows()

print(grep_hits, n=Inf)

# Sanity-check across all assays in case markers live in RNA but not SCT
walk(Assays(epithelial), function(a) {
  hits = intersect(all_markers_unique, rownames(epithelial[[a]]))
  cat('assay', a, ':', length(hits), '/', length(all_markers_unique),
      'markers present\n')
})



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
cd4_tcells = FindMarkers(integrated_seurat,
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
cd4_tcells = FindMarkers(integrated_seurat,
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
bcells = FindMarkers(integrated_seurat,
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
bcells = FindMarkers(integrated_seurat,
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
bcells = FindMarkers(integrated_seurat,
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
adipocytes = FindMarkers(integrated_seurat,
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
adipocytes = FindMarkers(integrated_seurat,
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
fibroblasts = FindMarkers(integrated_seurat,
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
fibroblasts = FindMarkers(integrated_seurat,
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


## epithelial ----
# clusters 23 vs 18 and 24
epithelial = FindMarkers(integrated_seurat,
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
epithelial = FindMarkers(integrated_seurat,
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



## endothelial ----

# myeloid ----

myeloid = subcluster_by_celltype(
  integrated_seurat,
  celltype_col    = "major_celltype",
  celltype_values = "Myeloid",
  sct_mode        = "per_sample",
  split_by        = "sample",
  batch_var       = "batch",
  vars_to_regress = c("cc.difference", "percent.mt"),
  resolutions     = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
  seed            = seed
)

myeloid <- subset(integrated_seurat, subset = major_celltype == "Myeloid")
myeloid <- DietSeurat(myeloid, assays = "RNA")
gc()

options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB
tmp <- myeloid
tmp[['RNA']] <- split(
  tmp[['RNA']],
  f = tmp$sample
)

myeloid <- SCTransform(tmp,
                      vars.to.regress = c('percent.mt', 'cc.difference'), # https://satijalab.org/seurat/articles/cell_cycle_vignette.html
                      variable.features.n = 3000
                )
rm(tmp) # free up memory
gc()

myeloid <- RunPCA(myeloid, assay = "SCT", npcs = 50)
test <- check_pc_contamination(myeloid, species='mouse', top_n=30)
test$flagged_pcs
test$loadings
test$summary

elbow <- calc_pcs_elbow(myeloid)
plot_pcs_elbow(elbow)

DimPlot(myeloid, reduction = "pca", group.by = "timepoint", split.by = 'Phase')
FeaturePlot(myeloid, reduction = "pca", features = c('percent.mt','S.Score','G2M.Score'),
             ncol=3,
             # split.by = 'condition',
             # slot='scale.data',
             # slot='data',
            # min.cutoff='q10',
            # max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)

myeloid <- FindNeighbors(myeloid, reduction = "pca", dims = 1:50)
resolutions = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
myeloid <- FindClusters(myeloid, resolution = resolutions, algorithm = 4)
myeloid <- RunUMAP(myeloid, reduction = "pca", dims = 1:50)


myeloid <- RunHarmony(myeloid, 'batch')
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:50,
                         graph.name=c('harmony_nn','harmony_snn'))
myeloid <- FindClusters(myeloid, resolution = resolutions, algorithm = 4, graph.name = "harmony_snn")
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:50, reduction.name = "umap_harmony")

# qs_save(myeloid, file.path(data_dir, "myeloid.qs"), nthreads = nthreads)
myeloid <- qs_read(file.path(data_dir, "myeloid.qs"), nthreads = nthreads)

DimPlot(subset(myeloid, subset = timepoint == "wk3"),
        reduction = "umap_harmony",
        group.by  = "sample",
        # split.by = 'condition'
        # label     = TRUE,
        repel     = TRUE,
        shuffle = TRUE
      )

DimPlot(myeloid,
        reduction = "umap_harmony",
        group.by  = "condition",
        # split.by = 'condition'
        # label     = TRUE,
        repel     = TRUE,
        shuffle = TRUE
      )



DimPlot(myeloid,
        reduction = "umap_harmony",
        group.by  = "celltype",
        # split.by = 'timepoint',
        repel     = TRUE,
        shuffle = TRUE
      )

DimPlot(myeloid,
        reduction = "umap",
        group.by  = "timepoint",
        # split.by = 'timepoint',
        repel     = TRUE,
        shuffle = TRUE
      )



FeaturePlot(myeloid, features = c('percent.mt'),
            reduction='umap_harmony',
            # min.cutoff='q10',
            # max.cutoff='q50',
            ncol=3,
            cols=viridis(256)) &
  DarkTheme()


VlnPlot(myeloid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",''),
        group.by = "harmony_snn_res.0.3", ncol = 3) &
  theme(legend.position = "none")


clustree::clustree(myeloid, prefix = 'harmony_snn_res.') +
  guides(edge_colour = "none")


if (!dir.exists(file.path(fig_dir, "myeloid"))) {
  dir.create(file.path(fig_dir, "myeloid"), recursive = TRUE)
}

DimPlot(
  myeloid,
  reduction = "umap_harmony",
  group.by  = "harmony_snn_res.1",
  label     = TRUE,
  repel     = TRUE,
  shuffle   = TRUE
)



DimPlot(
  myeloid,
  reduction = "umap_harmony",
  group.by  = "harmony_snn_res.0.25",
  label     = TRUE,
  repel     = TRUE,
  shuffle   = TRUE
)


ggsave(
  file.path(fig_dir, "myeloid", "umap_harmony_res.0.4.png"),
  dpi = 320, width = 10, height = 8
)

myeloid = JoinLayers(myeloid,assay='RNA')
gc()










# Rename cluster #'s to cell type annotations ----
Idents(integrated_seurat) = 'SCT_snn_res.0.8'
integrated_seurat = RenameIdents(integrated_seurat,
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

# integrated_seurat = RenameIdents(integrated_seurat,
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
integrated_seurat$celltype = Idents(integrated_seurat)
Idents(integrated_seurat) = 'celltype'

if (!dir.exists(paste0(annot_dir,'figures/'))) { dir.create(paste0(annot_dir,'figures/'),
                                                            recursive=TRUE) }

DimPlot(integrated_seurat, 
        group.by = 'celltype',
        label=T,
        label.size=6,
        repel=T,
        shuffle=T)
ggsave(paste0(annot_dir,'figures/annotated_seurat_res.0.8.png'), dpi=320, width=20, height=20)

# qs_save(integrated_seurat, paste0(data_dir,'annotated_seurat_res.0.8.qs'),nthreads=2)
qs_save(integrated_seurat, paste0(data_dir,'integrated_seurat_annotated.qs'),nthreads=2)

# plot qc metrics for cell-annotated clusters
if (!dir.exists(paste0(fig_dir,'celltype/'))) { dir.create(paste0(fig_dir,'celltype/'),
                                                           recursive=TRUE) }

plot_umaps(integrated_seurat, 
           resolution='celltype',
           directory=fig_dir)

umap_clustering_qc(integrated_seurat, 
                   resolution='celltype',
                   directory=fig_dir)

# shrink seurat object in preparation of DE analysis
de_seurat = DietSeurat(integrated_seurat, assays='RNA', layers=c('counts'))
qs_save(de_seurat, paste0(data_dir,'de_seurat_annotated.qs'),nthreads=2)

