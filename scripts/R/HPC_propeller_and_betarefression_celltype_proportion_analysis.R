library(Seurat)
library(speckle)
library(tidyverse)


nthreads = 4

analysis <- 'cellbender_analysis'
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

integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_sctype_annotated.qs'))


# # directories for saving
# output_annot_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/annotation/"
# 
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/clustering/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/data/"
# 
# output_stats_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/"
# 
# if (!dir.exists(paste0(output_data_dir))) { dir.create(paste0(output_data_dir),
#                                                        recursive=TRUE) }
# if (!dir.exists(paste0(output_fig_dir))) { dir.create(paste0(output_fig_dir),
#                                                       recursive=TRUE) }
# if (!dir.exists(paste0(output_stats_dir))) { dir.create(paste0(output_stats_dir),
#                                                         recursive=TRUE) }
# if (!dir.exists(paste0(output_annot_dir))) { dir.create(paste0(output_annot_dir),
#                                                        recursive=TRUE) }
# 
# # data directories
# data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/data/"
annotated_seurat = readRDS(paste0(data_dir_integrated,'de_seurat_annotated.rds')) # read in data

DefaultAssay(annotated_seurat)
Idents(annotated_seurat)

# remove clusters with only one condition
annotated_seurat = subset(annotated_seurat, 
                          idents=c('Muscle'),
                          invert=TRUE)

annotated_seurat$celltype = Idents(annotated_seurat) # reassign with clusters removed

  
p = propeller(clusters = annotated_seurat$celltype, sample = annotated_seurat$sample, 
          group = annotated_seurat$condition, transform='logit')
p
write.csv(p,pasteo(output_stats_dir,'propeller_props.csv'),row.names=FALSE)

plotCellTypeProps(clusters = annotated_seurat$celltype, sample = annotated_seurat$sample)
ggsave(paste0(output_fig_dir,'celltype/cell_type_proportions_bars_bysample.png'), dpi=320, width=20, height= 10)

plotCellTypeProps(clusters = annotated_seurat$celltype, sample = annotated_seurat$condition)
ggsave(paste0(output_fig_dir,'celltype/cell_type_proportions_bars_bycondition.png'), dpi=320, width=10, height= 10)













library(betareg)

## by cluster numbers 0:27
# Cluster Percentages by Pb
freq_table <- prop.table(table(annotated_seurat$seurat_clusters, annotated_seurat$condition), 2)
freq_table

### set up covariates
props <- prop.table(table(annotated_seurat$seurat_clusters,annotated_seurat$condition),2)
treat <- gsub('\\..*','',colnames(props))
props
treat

### beta regression
beta.res<-data.frame(mat.or.vec(28,3))
names(beta.res)<-c("cluster","coef","p-value")
props.zeros <- props
props.zeros[props.zeros==0] <- 1e-6

control_params <- betareg.control(hessian=T)

for (i in 0:27){
  #table(annotated_seurat$celltype==i, annotated_seurat$condition)
  test<-betareg(props.zeros[i+1,]~factor(treat),
                control=control_params)
  beta.res[i+1,1]<-i
  beta.res[i+1,2]<-round(test$coefficients$mean[2],3)
  beta.res[i+1,3]<-summary(test)$coefficients$mean[2,4]
}

beta.res


## by celltype annotation
# Cluster Percentages by Pb
freq_table <- prop.table(table(annotated_seurat$celltype, annotated_seurat$condition), 2)
freq_table

### set up covars
props <- prop.table(table(annotated_seurat$celltype,annotated_seurat$condition),2)
treat <- gsub('\\..*','',colnames(props))
props
treat

### beta regression
beta.res<-data.frame(mat.or.vec(18,3))
names(beta.res)<-c("cluster","coef","p-value")
props.zeros <- props
props.zeros[props.zeros==0] <- 1e-6

control_params <- betareg.control(hessian=T)

for (i in 0:17){
  #table(annotated_seurat$celltype==i, annotated_seurat$condition)
  test<-betareg(props.zeros[i+1,]~factor(treat),
                control=control_params)
  beta.res[i+1,1]<-i
  beta.res[i+1,2]<-round(test$coefficients$mean[2],3)
  beta.res[i+1,3]<-summary(test)$coefficients$mean[2,4]
}

beta.res

### combined clusters with same predicted cell type
# olig_combi <- props[3,] + props[11,]
# test_oli <- betareg(olig_combi~treat)
# summary(test_oli)
# predict(test_oli)
# 
# micg_combi <- props[2,] + props[10,]
# test_mcg <- betareg(micg_combi~treat)
# summary(test_mcg)
# predict(test_mcg)
# 
# peri_combi <- props[4,] + props[8,]
# test_peri <- betareg(peri_combi~treat)
# summary(test_peri)
# predict(test_peri)

# beta.res[2,3] <- summary(test_mcg)$coefficients$mean[2,4]
# beta.res[10,3] <- NA
# beta.res[3,3] <- summary(test_oli)$coefficients$mean[2,4]
# beta.res[11,3] <- NA
# beta.res[4,3] <- summary(test_peri)$coefficients$mean[2,4]
# beta.res[8,3] <- NA

### table of cluster % by Pb
freq_table <- data.frame(cluster=rownames(freq_table),Ctrl=freq_table[,1],Pb=freq_table[,2])
tab.clust_x_pb <- cbind(Ncells=as.numeric(table(annotated_seurat$celltype)),freq_table,beta.res[,2:3])
tab.clust_x_pb <- tab.clust_x_pb[,c(2,1,3:6)]
tab.clust_x_pb$Ctrl <- round(100*tab.clust_x_pb$Ctrl,2)
tab.clust_x_pb$Pb <- round(100*tab.clust_x_pb$Pb,2)
tab.clust_x_pb <- tab.clust_x_pb[,-5]
write.csv(tab.clust_x_pb, file=paste0(output_stats_dir,'beta_regression_cluster_props_by_pb.csv'), row.names = F)


###stacked bar chart

#set up data frame
props.m <- t(props)
props.m <- as.data.frame.matrix(props.m)
props.m$mouse <- rownames(props.m)

# #combine 1+9 and 2+10 and 3+7
# props.m$'1' <- props.m$'1' + props.m$'9'
# props.m$'2' <- props.m$'2' + props.m$'10'
# props.m$'3' <- props.m$'3' + props.m$'7'
# props.m <- props.m[,-c(8,10,11)]
# rowSums(props.m[,-10])
# names(props.m)[1:9] <- c('Endothelia','Microglia','Oligodendrocytes','Pericytes','Astrocytes',
#                          'Oligodendrocyte Progenitors','Choroid Plexus','Neurons','Fibroblasts')

#convert to long
long_props <- gather(props.m, cluster, pct, -mouse)
# long_props$cluster <- factor(long_props$cluster, 
#                              levels= c('Endothelia','Microglia','Oligodendrocytes','Pericytes','Astrocytes',
#                                        'Oligodendrocyte Progenitors','Choroid Plexus','Neurons','Fibroblasts'))


ggplot(data=long_props, aes(x=mouse, y=pct, fill=cluster)) +
  geom_bar(stat="identity")+ 
  scale_fill_discrete(name='') +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_text(hjust=-0.1, angle=-45),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12)) +
  xlab("") + ylab("Proportion of Cells")
ggsave(paste0(output_fig_dir,'celltype/cell_type_proportions_bars_bycondition_betareg.png'), dpi=320, width=10, height= 10)


#### sessionInfo() ----
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
#   [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1     
# [10] tidyverse_2.0.0    speckle_1.2.0      Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                polyclip_1.10-7             graph_1.80.0               
# [7] XML_3.99-0.17               fastDummies_1.7.4           lifecycle_1.0.4             edgeR_4.0.16                globals_0.16.3              lattice_0.22-6             
# [13] MASS_7.3-60.0.1             magrittr_2.0.3              limma_3.58.1                plotly_4.10.4               httpuv_1.6.15               sctransform_0.4.1          
# [19] spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0           cowplot_1.1.3               pbapply_1.7-2               DBI_1.2.3                  
# [25] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1        BiocGenerics_0.48.1        
# [31] RCurl_1.98-1.16             GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6               irlba_2.3.5.1              
# [37] listenv_0.9.1               spatstat.utils_3.1-0        GSVA_1.50.5                 goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-2      
# [43] annotate_1.80.0             fitdistrplus_1.2-1          parallelly_1.38.0           DelayedMatrixStats_1.24.0   leiden_0.4.3.1              codetools_0.2-20           
# [49] DelayedArray_0.28.0         tidyselect_1.2.1            farver_2.1.2                ScaledMatrix_1.10.0         matrixStats_1.4.1           stats4_4.3.1               
# [55] spatstat.explore_3.3-2      jsonlite_1.8.9              progressr_0.14.0            ggridges_0.5.6              survival_3.7-0              tools_4.3.1                
# [61] ica_1.0-3                   Rcpp_1.0.13                 glue_1.7.0                  gridExtra_2.3               SparseArray_1.2.4           DESeq2_1.42.1              
# [67] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.8         HDF5Array_1.30.1            withr_3.0.1                 fastmap_1.2.0               rhdf5filters_1.14.1        
# [73] fansi_1.0.6                 digest_0.6.37               rsvd_1.0.5                  timechange_0.3.0            R6_2.5.1                    mime_0.12                  
# [79] colorspace_2.1-1            scattermore_1.2             tensor_1.5                  spatstat.data_3.1-2         RSQLite_2.3.7               utf8_1.2.4                 
# [85] generics_0.1.3              data.table_1.16.0           httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.2.1              uwot_0.2.2                 
# [91] pkgconfig_2.0.3             gtable_0.3.5                blob_1.2.4                  lmtest_0.9-40               SingleCellExperiment_1.24.0 XVector_0.42.0             
# [97] htmltools_0.5.8.1           dotCall64_1.1-1             GSEABase_1.64.0             scales_1.3.0                Biobase_2.62.0              png_0.1-8                  
# [103] spatstat.univar_3.0-1       rstudioapi_0.16.0           tzdb_0.4.0                  reshape2_1.4.4              nlme_3.1-166                cachem_1.1.0               
# [109] zoo_1.8-12                  rhdf5_2.46.1                KernSmooth_2.23-24          parallel_4.3.1              miniUI_0.1.1.1              AnnotationDbi_1.64.1       
# [115] pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                 RANN_2.6.2                  promises_1.3.0              BiocSingular_1.18.0        
# [121] beachmat_2.18.1             xtable_1.8-4                cluster_2.1.6               cli_3.6.3                   locfit_1.5-9.10             compiler_4.3.1             
# [127] rlang_1.1.4                 crayon_1.5.3                future.apply_1.11.2         plyr_1.8.9                  stringi_1.8.4               viridisLite_0.4.2          
# [133] deldir_2.0-4                BiocParallel_1.36.0         munsell_0.5.1               Biostrings_2.70.3           lazyeval_0.2.2              spatstat.geom_3.3-3        
# [139] Matrix_1.6-5                RcppHNSW_0.6.0              hms_1.1.3                   patchwork_1.3.0             sparseMatrixStats_1.14.0    bit64_4.5.2                
# [145] future_1.34.0               Rhdf5lib_1.24.2             statmod_1.5.0               KEGGREST_1.42.0             shiny_1.9.1                 SummarizedExperiment_1.32.0
# [151] ROCR_1.0-11                 igraph_2.0.3                memoise_2.0.1               bit_4.5.0    