# remotes::install_github("Winnie09/GPTCelltype")
# install.packages("openai")
library(Seurat)
library(SeuratObject)
library(presto)
library(dplyr)
library(openai)
library(GPTCelltype)

Sys.setenv(OPENAI_API_KEY = '') # to generate prompt for chatGPT

all.pos.markers = read.csv('C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/annotation/annotated_pos_markers_res.0.8_df.csv')
n_genes = 25
results = gptcelltype(all.pos.markers, 
                   tissuename = 'developing 3 week old mouse mammary tissue', 
                   model = 'gpt-4',
                   topgenenumber=n_genes)
results # fed into chatgpt manually to verify celltypes with top n_genes

conserved.pos.markers = read.csv('aC:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/annotation/conserved_pos_markers_res.0.8_df.csv')


conserved.results = gptcelltype(conserved.pos.markers, 
                      tissuename = 'developing 3 week old mouse mammary tissue', 
                      model = 'gpt-4',
                      topgenenumber=n_genes)
conserved.results # fed into chatgpt manually to verify celltypes with top n_genes


# cell_markers_indata = read.csv('annotation/aggregate_sources_celltype_markers_actually_in_data.csv')

# # res 0.1
# Idents(seurat_umap_sct) = 'SCT_snn_res.0.1'
# seurat_umap_sct = RenameIdents(seurat_umap_sct,
#                                '0' = 'T cells',
#                                '1' = 'B cells',
#                                '2' = 'Tregs',
#                                '3' = 'Endothelial',
#                                '4' = 'Adipocytes',
#                                '5' = 'Fibroblasts',
#                                '6' = 'Fibroblasts',
#                                '7' = 'Proliferating',
#                                '8' = 'Neurons',
#                                '9' = 'Myeloid',
#                                '10' = 'Macrophages',
#                                '11' = 'Dendritic',
#                                '12' = 'Muscle')
# 
# seurat_umap_sct$gptcelltype_clusters = Idents(seurat_umap_sct)
# 
# 
# # res 0.2
# Idents(seurat_umap_sct) = 'SCT_snn_res.0.2'
# seurat_umap_sct = RenameIdents(seurat_umap_sct,
#                                '0' = 'B cells',
#                                '1' = 'T cells',
#                                '2' = 'Cytotoxic T cells',
#                                '3' = 'Tregs',
#                                '4' = 'Endothelial',
#                                '5' = 'Adipocytes',
#                                '6' = 'Fibroblasts',
#                                '7' = 'Proliferating',
#                                '8' = 'Neurons',
#                                '9' = 'Macrophages',
#                                '10' = 'Smooth Muscle / Fibroblasts',
#                                '11' = 'Myeloid',
#                                '12' = 'Dendritic',
#                                '13' = 'Endothelial',
#                                '14' = 'Myeloid',
#                                '15' = 'Muscle')
# 
# seurat_umap_sct$gptcelltype_clusters = Idents(seurat_umap_sct)


# res 1.4
Idents(seurat_umap_sct) = seurat_umap_sct$SCT_snn_res.0.8
seurat_umap_sct = RenameIdents(seurat_umap_sct,
                               '0' = 'B cells',
                               '1' = 'CD4+ T cells',
                               '2' = 'CD8+ T cells',
                               '3' = 'Tregs',
                               '4' = 'B cells',
                               '5' = 'Endothelial cells',
                               '6' = 'B cells', # refined assessment says T cells. Naive/quiescent or memory. Refined uses top 10 markers from all and conserved
                               '7' = 'Fibroblast',
                               '8' = 'Proliferating',
                               '9' = 'Dendritic cells', # or DC-like
                               '10' = 'Adipocytes', # possibly white adipose
                               '11' = 'Activated T Cells', # or immediate-early T cell cluster
                               '12' = 'Macrophages', # m2-skewed or tissue resident type
                               '13' = 'Innate lymphoid cells or γδ T cells ', #nk-like. an “innate-like” phenotype (sometimes called “NK‐like γδ T cells”)
                               '14' = 'Plasmacytoid dendritic cells',
                               '15' = 'Monocytes', # refined says M2-like macrophages
                               '16' = 'Fibroblasts (stromal)', # reticular cells? FRC-like. Specialized stroma or transitional cell states
                               '17' = 'Endothelial cells', # refined says lymphatic associated. or transitional between blood and lymphatic. possibly high endothelial venule‐like or sinusoidal)
                               '18' = 'Interferon response T cells', # 
                               '19' = 'Pericytes or myoepithelial', # or vascular smooth muscle. Subclustering shows small part is myoepithelial. large could be pericytes
                               '20' = 'Interferon response B cells',
                               '21' = 'Adipocytes', # possibly beige/brown adipose tissue
                               '22' = 'Macrophages', # potentially inflammatory or specialized
                               '23' = 'Schwann cells', # or glial-like ells
                               '24' = 'Conventional dendritic cells', #cDC1 lineage in mice. 
                               '25' = 'Luminal epithelial', # hormone responsive
                               '26' = 'Muscle', # skeletal muscle
                               '27' = 'Erythrocytes'# mature RBCs
                               )

seurat_umap_sct$gptcelltype_clusters = Idents(seurat_umap_sct)



### sessionInfo() ----
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
#   [1] GPTCelltype_1.0.1  openai_0.4.1       dplyr_1.1.4        presto_1.0.0       data.table_1.16.0  Rcpp_1.0.13        Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3              deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.4            magrittr_2.0.3         RcppAnnoy_0.0.22       spatstat.geom_3.3-3   
# [9] matrixStats_1.4.1      ggridges_0.5.6         compiler_4.3.1         png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          pkgconfig_2.0.3       
# [17] fastmap_1.2.0          utf8_1.2.4             promises_1.3.0         purrr_1.0.2            jsonlite_1.8.9         goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-0  
# [25] irlba_2.3.5.1          parallel_4.3.1         cluster_2.1.6          R6_2.5.1               ica_1.0-3              spatstat.data_3.1-2    stringi_1.8.4          RColorBrewer_1.1-3    
# [33] reticulate_1.39.0      spatstat.univar_3.0-1  parallelly_1.38.0      lmtest_0.9-40          scattermore_1.2        tensor_1.5             future.apply_1.11.2    zoo_1.8-12            
# [41] IRanges_2.36.0         sctransform_0.4.1      httpuv_1.6.15          Matrix_1.6-5           splines_4.3.1          igraph_2.0.3           tidyselect_1.2.1       abind_1.4-8           
# [49] rstudioapi_0.16.0      spatstat.random_3.3-2  spatstat.explore_3.3-2 codetools_0.2-20       miniUI_0.1.1.1         listenv_0.9.1          lattice_0.22-6         tibble_3.2.1          
# [57] plyr_1.8.9             Biobase_2.62.0         shiny_1.9.1            withr_3.0.1            ROCR_1.0-11            Rtsne_0.17             future_1.34.0          fastDummies_1.7.4     
# [65] survival_3.7-0         polyclip_1.10-7        fitdistrplus_1.2-1     pillar_1.9.0           KernSmooth_2.23-24     stats4_4.3.1           plotly_4.10.4          generics_0.1.3        
# [73] RcppHNSW_0.6.0         S4Vectors_0.40.2       ggplot2_3.5.1          munsell_0.5.1          scales_1.3.0           globals_0.16.3         xtable_1.8-4           glue_1.7.0            
# [81] lazyeval_0.2.2         tools_4.3.1            RSpectra_0.16-2        RANN_2.6.2             leiden_0.4.3.1         dotCall64_1.1-1        XML_3.99-0.17          cowplot_1.1.3         
# [89] grid_4.3.1             tidyr_1.3.1            colorspace_2.1-1       nlme_3.1-166           patchwork_1.3.0        cli_3.6.3              spatstat.sparse_3.1-0  spam_2.10-0           
# [97] fansi_1.0.6            viridisLite_0.4.2      uwot_0.2.2             gtable_0.3.5           digest_0.6.37          progressr_0.14.0       BiocGenerics_0.48.1    ggrepel_0.9.6         
# [105] htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7             mime_0.12              MASS_7.3-60.0.1 
