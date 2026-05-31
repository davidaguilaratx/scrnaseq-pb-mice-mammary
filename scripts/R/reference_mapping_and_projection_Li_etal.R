setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(SeuratObject)
library(sctransform)
library(ggpubr)
library(cowplot)

# load data
seurat_umap_sct = readRDS('data/annotated_seurat_scCATCH_0.8_df.rds')

## Create seurat object from reference data
# setwd('C:/Users/david/Documents/Research/PhD/scRNAseq/References work flows/Li_2020_age_associated_alterations/')
# 
# matrix = read.csv('GSE150580_Li_Brugge_10XscRNAseq_GeneCellMatrix_RNAcounts_GEO.csv',
#                   header=TRUE,
#                   row.names=1)
# metadata = read.csv('GSE150580_Li_Brugge_10XscRNAseq_Metadata_GEO.csv',
#                     header=TRUE,
#                     row.names=1)
# 
# setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")
# 
# ref = CreateSeuratObject(counts=matrix,meta.data=metadata)
# saveRDS(reference, 'data/Li2020_reference.rds')
#
# rm(matrix,metadata) # free RAM
# gc()

# load reference data ----
ref = readRDS('data/Li2020_reference.rds')
ref = subset(ref, m.AgeGroup=='Young')

# process reference data identically to sample data----
ref$cc.difference = ref$S.Score - ref$G2M.Score

# Check quarterly values of mitochondrial gene expression
ref$percent.mt = ref$mitoRatio*100 # turn proportions into percentages
summary(ref$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02353 0.03407 0.03704 0.04723 0.09995 

# Turn percent.mt into categorical factor vector based on quartile values
ref@meta.data$mito.level = cut(ref@meta.data$percent.mt, 
                                        breaks=c(-Inf, 2.353, 3.407, 4.723, Inf), 
                                        labels=c("1st","2nd","3rd", "4th"))

ref = NormalizeData(ref)
ref = FindVariableFeatures(ref)
ref = ScaleData(ref)
ref = RunPCA(ref, assay='RNA', npcs=50)
# 
# DimPlot(ref,
#         reduction='pca',
#         group.by='mito.level',
#         repel=TRUE,
#         split.by='mito.level')+
#   labs(title='PCA by mitochondrial expresion quartile',
#        color='mito.level')+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.title = element_text(hjust=0.5, face="bold"),
#         plot.subtitle = element_text(hjust=0.5),
#         axis.title = element_text(size=12),
#         axis.text = element_text(size=12),
#         title = element_text(size=18))
# 
# DimPlot(ref,
#         reduction='pca',
#         group.by='Phase',
#         repel=TRUE,
#         split.by='Phase')+
#   labs(title='PCA by cell cycle phase',
#        color='phase')+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.title = element_text(hjust=0.5, face="bold"),
#         plot.subtitle = element_text(hjust=0.5),
#         axis.title = element_text(size=12),
#         axis.text = element_text(size=12),
#         title = element_text(size=18))


# # normalize reference data with sct transform ----
# ref = SCTransform(ref, assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
#                 vars.to.regress=c('cc.difference','percent.mt'), 
#                 return.only.var.genes=FALSE) 


top25_sct = head(VariableFeatures(ref), 25)
p = VariableFeaturePlot(ref, selection.method='sct', assay='SCT')
LabelPoints(plot=p, points = top25_sct, repel = TRUE)+
labs(title='Top 25 variably expressed genes in Li et al. reference data')+
theme_classic2()


# # calculate number of principal componentes to use for reference data ----
# pct = ref[["pca"]]@stdev / sum(ref[["pca"]]@stdev) * 100
# pct
# #
# # # Calculate cumulative percents for each PC
# cumu = cumsum(pct)
# cumu
# #
# # # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
# co1 = which(cumu >= 90 & pct <= 5)[1]
# co1 # 41
# #
# # # Determine the difference between variation of PC and subsequent PC.
# # # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
# co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# co2 # 17
# #
# # # Minimum of the two calculations
# pcs = min(co1, co2)
# pcs 
# # [1] 17
# 
# # Create a dataframe to visualize percent variation captured by PCs
# plot_df = data.frame(pct = pct, 
#                      cumu = cumu, 
#                      rank = 1:length(pct))
# 
# ### Elbow plot to visualize percent variation captured by PCs ----
# ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
#   geom_text() + 
#   geom_vline(xintercept = 90, color = "grey") + 
#   geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
#   theme_classic2()+
#   scale_x_continuous(n.breaks=10)+
#   scale_y_continuous(n.breaks=10)+
#   labs(title = 'Eblow Plot of percent variation explained by PCs',
#        x='Cumulative percentage',
#        y='Percent of variation')+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         plot.title = element_text(hjust=0.5, face="bold"),
#         plot.subtitle = element_text(hjust=0.5),
#         axis.title = element_text(size=12), 
#         axis.text = element_text(size=12),
#         title = element_text(size=18))




# run PCA and UMAP on reference data ----
# ref = RunPCA(ref, assay='SCT', npcs=50)
ref = FindNeighbors(ref, dims=1:50)
ref = FindClusters(ref)

ref = RunUMAP(ref,
              dims=1:50,
              reduction='pca',
              return.model=TRUE) # need to return model to run MapQuery

# saveRDS(ref, 'data/Li2020_reference_umap_sct.rds')
ref = readRDS('data/Li2020_reference_umap_sct.rds')

Idents(ref) = 'SCT_snn_res.0.8'
DimPlot(ref,
        reduction = "umap",
        group.by='m.CellTypes',
        label.size=6,
        label=T)


# Find transfer anchors between reference and data ----
anchors = FindTransferAnchors(reference=ref, query=seurat_umap_sct,
                              normalization.method='LogNormalize',
                              reference.reduction='pca',
                              dims=1:50)

gc() # clean up RAM
# saveRDS(anchors,'data/anchors_Li_etal')
anchors = readRDS('data/anchors_Li_etal')

seurat_umap_sct = MapQuery(anchorset=anchors, reference=ref, query=seurat_umap_sct,
                           refdata=list(celltype="m.CellTypes"), reference.reduction="pca", 
                           reduction.model="umap")

gc() # free RAM
saveRDS(seurat_umap_sct, 'data/annotated_seurat_Li_etal_res.0.8_df.rds')

DimPlot(seurat_umap_sct, 
        reduction = "ref.umap", 
        group.by = "predicted.celltype", 
        label = TRUE, 
        label.size = 6, 
        repel = TRUE)

DimPlot(seurat_umap_sct, 
        reduction = "umap", 
        group.by = "predicted.celltype", 
        label = TRUE, 
        label.size = 6)

### "de novo visualization ----
ref = DietSeurat(ref, dimreducs='pca')
seurat_umap_sct = DietSeurat(seurat_umap_sct, dimreducs= c('pca','umap','ref.pca'))

#merge reference and query
ref$id = 'reference'
seurat_umap_sct$id = 'query'
refquery = merge(ref, seurat_umap_sct)
refquery[['pca']] =  merge(ref[['pca']], seurat_umap_sct[['ref.pca']])

# umap and plot by id
refquery = RunUMAP(refquery, reduction='pca', dims=1:50)
saveRDS(refquery,'data/refquery_Li_etal')
refquery = readRDS(refquery,'data/refquery_Li_etal')


p1 = DimPlot(refquery, 
        group.by='id', 
        shuffle=TRUE)

p2 = DimPlot(refquery, 
        group.by='m.CellTypes', 
        shuffle=TRUE)

p1+p2






### sessionInfo()
sessionInfo()
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
#   [1] cowplot_1.1.3      ggpubr_0.6.0       ggplot2_3.5.1      sctransform_0.4.1  Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                tibble_3.2.1                polyclip_1.10-7            
# [7] graph_1.80.0                XML_3.99-0.17               fastDummies_1.7.4           lifecycle_1.0.4             rstatix_0.7.2               globals_0.16.3             
# [13] lattice_0.22-6              MASS_7.3-60.0.1             backports_1.5.0             magrittr_2.0.3              plotly_4.10.4               httpuv_1.6.15              
# [19] spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0           pbapply_1.7-2               DBI_1.2.3                   RColorBrewer_1.1-3         
# [25] abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1        purrr_1.0.2                 BiocGenerics_0.48.1        
# [31] RCurl_1.98-1.16             GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6               irlba_2.3.5.1              
# [37] listenv_0.9.1               spatstat.utils_3.1-0        GSVA_1.50.5                 goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-2      
# [43] annotate_1.80.0             fitdistrplus_1.2-1          parallelly_1.38.0           DelayedMatrixStats_1.24.0   leiden_0.4.3.1              codetools_0.2-20           
# [49] DelayedArray_0.28.0         tidyselect_1.2.1            farver_2.1.2                ScaledMatrix_1.10.0         matrixStats_1.4.1           stats4_4.3.1               
# [55] spatstat.explore_3.3-2      jsonlite_1.8.9              Formula_1.2-5               progressr_0.14.0            ggridges_0.5.6              survival_3.7-0             
# [61] tools_4.3.1                 ica_1.0-3                   Rcpp_1.0.13                 glue_1.7.0                  gridExtra_2.3               SparseArray_1.2.4          
# [67] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.8         dplyr_1.1.4                 HDF5Array_1.30.1            withr_3.0.1                 fastmap_1.2.0              
# [73] rhdf5filters_1.14.1         fansi_1.0.6                 digest_0.6.37               rsvd_1.0.5                  R6_2.5.1                    mime_0.12                  
# [79] colorspace_2.1-1            scattermore_1.2             tensor_1.5                  spatstat.data_3.1-2         RSQLite_2.3.7               utf8_1.2.4                 
# [85] tidyr_1.3.1                 generics_0.1.3              data.table_1.16.0           httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.2.1             
# [91] uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.5                blob_1.2.4                  lmtest_0.9-40               SingleCellExperiment_1.24.0
# [97] XVector_0.42.0              htmltools_0.5.8.1           carData_3.0-5               dotCall64_1.1-1             GSEABase_1.64.0             scales_1.3.0               
# [103] Biobase_2.62.0              png_0.1-8                   spatstat.univar_3.0-1       rstudioapi_0.16.0           reshape2_1.4.4              nlme_3.1-166               
# [109] cachem_1.1.0                zoo_1.8-12                  rhdf5_2.46.1                stringr_1.5.1               KernSmooth_2.23-24          parallel_4.3.1             
# [115] miniUI_0.1.1.1              AnnotationDbi_1.64.1        pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                 RANN_2.6.2                 
# [121] promises_1.3.0              car_3.1-3                   BiocSingular_1.18.0         beachmat_2.18.1             xtable_1.8-4                cluster_2.1.6              
# [127] cli_3.6.3                   compiler_4.3.1              rlang_1.1.4                 crayon_1.5.3                ggsignif_0.6.4              future.apply_1.11.2        
# [133] plyr_1.8.9                  stringi_1.8.4               viridisLite_0.4.2           deldir_2.0-4                BiocParallel_1.36.0         munsell_0.5.1              
# [139] Biostrings_2.70.3           lazyeval_0.2.2              spatstat.geom_3.3-3         Matrix_1.6-5                RcppHNSW_0.6.0              patchwork_1.3.0            
# [145] sparseMatrixStats_1.14.0    bit64_4.5.2                 future_1.34.0               Rhdf5lib_1.24.2             KEGGREST_1.42.0             shiny_1.9.1                
# [151] SummarizedExperiment_1.32.0 ROCR_1.0-11                 broom_1.0.7                 igraph_2.0.3                memoise_2.0.1               bit_4.5.0