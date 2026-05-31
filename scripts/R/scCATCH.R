setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(scCATCH)

# sc_data is the scRNA-seq data matrix 
# sc_cluster is a character containing the cluster information


seurat_umap_sct = readRDS('data/annotated_seurat_Li_etal_sctype.rds')
seurat_umap_sct = JoinLayers(seurat_umap_sct, assay='RNA')

## resolution 1.4 ----
Idents(seurat_umap_sct) = 'SCT_snn_res.0.8' # clusters to annotate
seurat_umap_sct$seurat_clusters = Idents(seurat_umap_sct) # set # of clusters to 'SCT_snn_res.1.4'

obj = createscCATCH(data = seurat_umap_sct@assays[['RNA']]$data, 
                    cluster = as.character(seurat_umap_sct@meta.data$seurat_clusters))


# find marker gene for each cluster
obj = findmarkergene(obj, species='Mouse',marker=cellmatch,tissue='Mammary gland')

# find cell type for each cluster
obj = findcelltype(obj)

saveRDS(obj, 'annotation/scCATCH_annotations_res.0.8_df.rds')

obj@celltype



# # Rename cluster #'s to cell type annotations from scCATCH at 0.2 resolution ----
# seurat_umap_sct = RenameIdents(seurat_umap_sct,
#                                '6' = 'Stromal',
#                                '0' = 'B cell',
#                                '12' = 'Macrophage',
#                                '5' = 'Endo/Luminal/Prog/Stromal',
#                                '4' = 'Stem cell',
#                                '9' = 'Macrophage',
#                                '7' = 'Dividing',
#                                '8' = 'Dendritic',
#                                '11'= 'Macrophage',
#                                '13' = 'Stromal',
#                                '14' = 'Macrophage',
#                                '15' = 'Luminal')

# # Rename cluster #'s to cell type annotations from scCATCH at 1.4 resolution ----
# seurat_umap_sct = RenameIdents(seurat_umap_sct,
#                                '6' = 'Endothelial',
#                                '7' = 'Stromal',
#                                '9' = 'Stromal',
#                                '11' = 'Dendritic',
#                                '14' = 'Macrophage',
#                                '16' = 'Macrophage',
#                                '18' = 'Stromal',
#                                '20' = 'Myoepithelial',
#                                '22' = 'Luminal',
#                                '23' = 'Macrophage',
#                                '24' = 'Macrophage',
#                                '25' = 'Luminal Progenitor',
#                                '27' = 'Luminal Progenitor',
#                                '28' = 'Myoepithelial',
#                                '29' = 'Stem')

# Rename cluster #'s to cell type annotations from scCATCH at 0.8 resolution ----
seurat_umap_sct = RenameIdents(seurat_umap_sct,
                               '5' = 'Endothelial',
                               '6' = 'Stromal',
                               '7' = 'Stromal',
                               '9' = 'Dendritic',
                               '12' = 'Macrophage',
                               '14' = 'Macrophage',
                               '15' = 'Macrophage',
                               '17' = 'Stromal',
                               '18' = 'Myoepithelial',
                               '21' = 'Luminal',
                               '22' = 'Macrophage',
                               '23' = 'Luminal Progenitor',
                               '25' = 'Luminal Progenitor',
                               '26' = 'Myoepithelial',
                               '27' = 'Stem')
                               
seurat_umap_sct$scCATCH_clusters = Idents(seurat_umap_sct)

saveRDS(seurat_umap_sct, 'data/annotated_seurat_scCATCH_0.8_df.rds')

# # resolution 1.4
# Idents(seurat_umap_sct) = 'SCT_snn_res.1.4' # clusters to annotate
# seurat_umap_sct$seurat_clusters = Idents(seurat_umap_sct) # set # of clusters to 'SCT_snn_res.1.4'
# 
# obj = createscCATCH(data = seurat_umap_sct@assays[['RNA']]$data, 
#                      cluster = as.character(seurat_umap_sct@meta.data$seurat_clusters))
# 
# 
# # find marker gene for each cluster
# obj = findmarkergene(obj, species='Mouse',marker=cellmatch,tissue='Mammary gland')
# 
# # find cell type for each cluster
# obj = findcelltype(obj)
# 
# saveRDS(obj, 'annotation/scCATCH_annotations.rds')
# 
# obj@celltype
# 
# 
# # Rename cluster #'s to cell type annotations from scCATCH ----
# seurat_umap_sct = RenameIdents(seurat_umap_sct,
#                                '6' = 'Endothelial',
#                                '7' = 'Stromal',
#                                '9' = 'Stromal',
#                                '11' = 'Dendritic',
#                                '14' = 'Macrophage',
#                                '16' = 'Macrophage',
#                                '18' = 'Stromal',
#                                '20' = 'Myoepithelial',
#                                '22' = 'Luminal',
#                                '23' = 'Macrophage',
#                                '24' = 'Macrophage',
#                                '25' = 'Luminal Progenitor',
#                                '27' = 'Luminal Progenitor',
#                                '28' = 'Myoepithelial',
#                                '29' = 'Stem')
# 
# seurat_umap_sct$scCATCH_clusters = Idents(seurat_umap_sct)
# 
# saveRDS(seurat_umap_sct, 'data/annotated_seurat_scCATCH.rds')






### sessionInfo() ----
sessionInfo()
# version 4.3.1 (2023-06-16 ucrt)
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
#   [1] scCATCH_3.2.2      Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                tibble_3.2.1                polyclip_1.10-7            
# [7] graph_1.80.0                XML_3.99-0.17               fastDummies_1.7.4           lifecycle_1.0.4             globals_0.16.3              lattice_0.22-6             
# [13] MASS_7.3-60.0.1             magrittr_2.0.3              plotly_4.10.4               httpuv_1.6.15               sctransform_0.4.1           spam_2.10-0                
# [19] spatstat.sparse_3.1-0       reticulate_1.39.0           cowplot_1.1.3               pbapply_1.7-2               DBI_1.2.3                   RColorBrewer_1.1-3         
# [25] abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1        purrr_1.0.2                 BiocGenerics_0.48.1        
# [31] RCurl_1.98-1.16             GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6               irlba_2.3.5.1              
# [37] listenv_0.9.1               spatstat.utils_3.1-0        GSVA_1.50.5                 goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-2      
# [43] annotate_1.80.0             fitdistrplus_1.2-1          parallelly_1.38.0           DelayedMatrixStats_1.24.0   leiden_0.4.3.1              codetools_0.2-20           
# [49] DelayedArray_0.28.0         tidyselect_1.2.1            farver_2.1.2                ScaledMatrix_1.10.0         matrixStats_1.4.1           stats4_4.3.1               
# [55] spatstat.explore_3.3-2      jsonlite_1.8.9              progressr_0.14.0            ggridges_0.5.6              survival_3.7-0              progress_1.2.3             
# [61] tools_4.3.1                 ica_1.0-3                   Rcpp_1.0.13                 glue_1.7.0                  gridExtra_2.3               SparseArray_1.2.4          
# [67] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.8         dplyr_1.1.4                 HDF5Array_1.30.1            withr_3.0.1                 fastmap_1.2.0              
# [73] rhdf5filters_1.14.1         fansi_1.0.6                 digest_0.6.37               rsvd_1.0.5                  R6_2.5.1                    mime_0.12                  
# [79] colorspace_2.1-1            scattermore_1.2             tensor_1.5                  spatstat.data_3.1-2         RSQLite_2.3.7               utf8_1.2.4                 
# [85] tidyr_1.3.1                 generics_0.1.3              data.table_1.16.0           prettyunits_1.2.0           httr_1.4.7                  htmlwidgets_1.6.4          
# [91] S4Arrays_1.2.1              uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.5                blob_1.2.4                  lmtest_0.9-40              
# [97] SingleCellExperiment_1.24.0 XVector_0.42.0              htmltools_0.5.8.1           dotCall64_1.1-1             GSEABase_1.64.0             scales_1.3.0               
# [103] Biobase_2.62.0              png_0.1-8                   spatstat.univar_3.0-1       rstudioapi_0.16.0           reshape2_1.4.4              nlme_3.1-166               
# [109] cachem_1.1.0                zoo_1.8-12                  rhdf5_2.46.1                stringr_1.5.1               KernSmooth_2.23-24          parallel_4.3.1             
# [115] miniUI_0.1.1.1              AnnotationDbi_1.64.1        pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                 RANN_2.6.2                 
# [121] promises_1.3.0              BiocSingular_1.18.0         beachmat_2.18.1             xtable_1.8-4                cluster_2.1.6               cli_3.6.3                  
# [127] compiler_4.3.1              rlang_1.1.4                 crayon_1.5.3                future.apply_1.11.2         plyr_1.8.9                  stringi_1.8.4              
# [133] viridisLite_0.4.2           deldir_2.0-4                BiocParallel_1.36.0         munsell_0.5.1               Biostrings_2.70.3           lazyeval_0.2.2             
# [139] spatstat.geom_3.3-3         Matrix_1.6-5                RcppHNSW_0.6.0              hms_1.1.3                   patchwork_1.3.0             sparseMatrixStats_1.14.0   
# [145] bit64_4.5.2                 future_1.34.0               ggplot2_3.5.1               Rhdf5lib_1.24.2             KEGGREST_1.42.0             shiny_1.9.1                
# [151] SummarizedExperiment_1.32.0 ROCR_1.0-11                 igraph_2.0.3                memoise_2.0.1               bit_4.5.0 