###Let's do the cell type identification, see what our best bets for cell types are using scType

###Let's test ScType for cell type ID for fun
# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T) # load packages
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

set.seed(2024)

# directories for saving
output_annot_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
# output_annot_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0_MAD/annotation/"
# output_annot_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_month/cellbender_analysis/FPR_0.0_MAD/annotation/"
# output_annot_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0_MAD/annotation/"

output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0_MAD/qc/"
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weekshs/cellbender_analysis/FPR_0.0_MAD/qc/"
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0_MAD/qc/"
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0_MAD/qc/"

output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0_MAD/data/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0_MAD/data/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0_MAD/data/"


# seurat_umap_sct = readRDS('data/seurat_umap_sct_rm29and32_regressout_cc.diffandmt_doubletfinder.rds') # read in data
seurat_umap_sct = readRDS(paste0(output_data_dir,'seurat_umap_sct_regressout_cc.diffandmt.rds'))

# Scale RNA data for scType ----
DefaultAssay(seurat_umap_sct) = "RNA"
seurat_umap_sct_rna_sc = ScaleData(seurat_umap_sct, features = rownames(seurat_umap_sct))

# build custom scType database ----
# cell_markers_indata = read.csv('annotation/aggregate_sources_celltype_markers_actually_in_data.csv')
# 
# # read in sctype immune cell excel file
# sctype = read_excel('annotation/ScTypeDB_immune.xlsx')
# 
# # make gene lists lowercase for proper mice gene names
# sctype$geneSymbolmore1 = lapply(sctype$geneSymbolmore1, function(x) {
#   genes = strsplit(x, ",")
#   genes_titlecase = sapply(genes, function(g) paste(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))), sep = ""))
#   return(paste(genes_titlecase, collapse = ","))
# })
# 
# sctype$geneSymbolmore2 = lapply(sctype$geneSymbolmore2, function(x) {
#   genes = strsplit(x, ",")
#   genes_titlecase = sapply(genes, function(g) paste(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))), sep = ""))
#   return(paste(genes_titlecase, collapse = ","))
# })
# 
# # convert NANA into NA
# sctype$geneSymbolmore2 = lapply(sctype$geneSymbolmore2, function(x) {
#   ifelse(x == "NANA", '', x)
#   })
# 
# write.xlsx(sctype, 'annotation/ScTypeDB_immune_processed.xlsx')
# sctype = read_excel('annotation/ScTypeDB_immune_processed.xlsx')
# 
# 
# # # combine custom markers with sctype file into new df ----
# cell_markers_sctype = lapply(cell_markers_indata, function (x) {
#   c(names(x),paste(x[x != ""], collapse = ","))
# }) %>%
#   bind_rows() %>%
#   pivot_longer(cols=everything(), names_to='cellName', values_to='geneSymbolmore1')
# 
# sctype_custom = rbind.fill(sctype,cell_markers_sctype) # combine dfs
# sctype_custom$tissueType = 'mammary' # fill in tissue type
# sctype_custom = sctype_custom[which(!is.na(sctype_custom$cellName)),] # get rid of blanks
# write.xlsx(sctype_custom, 'annotation/ScTypeDB_mammary.xlsx')

## build scType database from cellmarker 2.0 ----

cm2.0 = read.xlsx('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/marker_genes/Cell_marker_All.xlsx')
cm2.0 = subset(cm2.0, subset=(species == 'Mouse') &
                             ((tissue_class == 'Mammary gland') |
                                (tissue_class == 'Breast') |
                                (tissue_class == 'Blood') |
                                (tissue_class == 'Blood vessel') |
                                (tissue_class == 'Lymph') |
                                (tissue_class == 'Abdomen') |
                                (tissue_class == 'Adipose tissue') |
                                (tissue_class == 'Nerve'))
)
cm2.0$marker = str_to_title(cm2.0$marker)
cm2.0[(cm2.0 == 'c-Kit') | (cm2.0 == 'C-Kit')] = 'Kit'
cm2.0[cm2.0 == 'Vimentin'] = 'Vim'
cm2.0[cm2.0 == 'Sca-1'] = 'Sca1'
cm2.0[cm2.0 == 'Pecaml'] = 'Pecam1'
cm2.0[cm2.0 == 'Granzyme B'] = 'Gzmb'
cm2.0[cm2.0 == 'Gr-1'] = 'Gr1'
cm2.0[cm2.0 == 'Fc-Epsilon Ri-Alpha'] = 'Fcer1a'
cm2.0[cm2.0 == 'Alpha Actin'] = 'Acta2'
cm2.0[cm2.0 == 'Asc-1'] = 'Asc1'
cm2.0[cm2.0 == 'Alpha Actin'] = 'Acta2'
cm2.0[cm2.0 == 'F4/80'] = 'Adgre1'
cm2.0[cm2.0 == 'Mmcp-8'] = 'Mmcp8'
cm2.0[cm2.0 == 'Moma-2'] = 'Moma2'
cm2.0[cm2.0 == 'N-Cadherin'] = 'Cdh2'
cm2.0[cm2.0 == 'E-Cadherin'] = 'Cdh1'
cm2.0[(cm2.0 == 'A-Sma') | (cm2.0 == 'Asma')] = 'Acta2'
cm2.0[cm2.0 == 'Fabp-4'] = 'Fabp4'
cm2.0[cm2.0 == 'Cytokeratin 14'] = 'Krt14'
cm2.0[cm2.0 == 'Ck14'] = 'Krt14'
cm2.0[cm2.0 == 'Ck5'] = 'Krt5'
cm2.0[cm2.0 == 'Cox-2'] = 'Cox2'
cm2.0[cm2.0 == 'Ppargc1α'] = 'Ppargc1a'
cm2.0[cm2.0 == 'C-Myb'] = 'Myb'
cm2.0[cm2.0 == 'Pdgfrα'] = 'Pdgfra'
cm2.0[cm2.0 == 'Pdgfrβ'] = 'Pdgfrb'
cm2.0[cm2.0 == 'Cd8'] = 'Cd8a'

















## build scType database from PanglaoDB ----




# run scType ----

## load DB file ----
###I made a custom cell type list to use that includes immune and mammary cells
db_ = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/marker_genes/scType/ScTypeDB_mammary.xlsx";
tissue = c("mammary") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

gs_list = gene_sets_prepare(db_, tissue)

# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 = isFALSE('counts' %in% names(attributes(seurat_umap_sct[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))



## extract scaled scRNA-seq matrix ----
# scRNAseqData_scaled = if (seurat_package_v5) as.matrix(seurat_umap_sct_rna_sc[["RNA"]]$scale.data) else as.matrix(seurat_umap_sct_rna_sc[["RNA"]]@scale.data)
scRNAseqData_scaled = if (seurat_package_v5) as.matrix(seurat_umap_sct[["SCT"]]$data) else as.matrix(seurat_umap_sct_rna_sc[["RNA"]]@data)


rm(seurat_umap_sct_rna_sc) # free RAM
gc()

## run ScType ----
es.max = sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.8' # clusters to annotate
seurat_umap_sct$seurat_clusters = Idents(seurat_umap_sct) # set # of clusters to 'SCT_snn_res.1.4'

### merge ScType results by cluster ----
cL_results = do.call("rbind", lapply(unique(seurat_umap_sct@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_umap_sct@meta.data[seurat_umap_sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_umap_sct@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

### set low-confident (low ScType score) clusters to "unknown" 
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])  ###Oh, now we're in business, this definitely seems to have worked

sctype_scores = sctype_scores[order(sctype_scores$cluster),] # order by cluster id#

### assign scType calssifications to seurat object ----
seurat_umap_sct@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_umap_sct@meta.data$sctype_classification[seurat_umap_sct@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# plot sctype annotations
DimPlot(seurat_umap_sct, 
        reduction = "umap", 
        label = TRUE,
        shuffle=TRUE,
        group.by = 'sctype_classification')  


# replace cluster names
sctype_scores$type<-gsub('ISG expressing immune cells',"T helper cells",sctype_scores$type)

###Vascular endothelial cells is too long, replace with "Vascular Endo"
sctype_scores$type<-gsub("Vascular.endothelial","Vascular Endo",sctype_scores$type)

###"Memory CD4+ T cells ends up being too long, replace with T cells
sctype_scores$type<-gsub("γδ-T cells","GammaDelta-Treg",sctype_scores$type)

###Lymphatic endothelial cells is too long, replace with "Lymphatic Endo"
sctype_scores$type<-gsub("Myofibroblasts","Myofibroblasts or Pericyte",sctype_scores$type)

write.csv(sctype_scores, paste0(output_annot_dir,'sctype_scores_regressout_ccdiffandmt_SCT_snn_res.0.8.csv'), row.names=FALSE)

seurat_umap_sct@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_umap_sct@meta.data$sctype_classification[seurat_umap_sct@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


# plot sctype annotations
DimPlot(seurat_umap_sct, 
        reduction = "umap", 
        label = TRUE,
        shuffle = TRUE,
        repel=TRUE,
        group.by = 'sctype_classification',
        label.size = 6) +
  guides(color = guide_legend(override.aes = list(size=4), ncol=1) )
ggsave(paste0(output_fig_dir,'SCT_snn_res.0.8/ScType_annotations_res.0.8.png'), dpi=320, width=20, height=20)

rm(scRNAseqData_scaled) # clean RAM
gc()

saveRDS(seurat_umap_sct, paste0(output_data_dir,'annotated_seurat_Li_etal_sctype.rds'))








## from justin's script ----

paste0(sctype_scores$type,"_",sctype_scores$cluster)
ScType.cluster<-paste0(sctype_scores$type,"-",sctype_scores$cluster)

names(ScType.cluster) <- levels(flex)

###Set the idents back to cluster instead of treatment, to rename with our best guess cell types
Idents(flex) <- "seurat_clusters"
names(ScType.cluster) <- levels(flex)

flex <- RenameIdents(flex, ScType.cluster)
levels(flex) ###This looks good!

flex$ScType.cluster<-Idents(flex)

DimPlot(flex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(flex, reduction = "umap", label = TRUE, pt.size = 0.5,split.by="Treatment_Tumor") + NoLegend()

###I want to also split by individual to satisfy my curiosity
flex$Treatment_Tumor_Ind<-paste(flex$Treatment_Tumor,flex$MouseID,sep="_")

DimPlot(flex, reduction = "umap",split.by="Treatment_Tumor_Ind") ###Seems to be pretty consistent, although this plot didn't really work

###What's going on with Clusters 1 and 8 (pretty much only in the tumor)
VlnPlot(flex,features=c("Epcam","Cd19","Cd4"))

###Calculate the individual level cell type proportions, this will give me an idea if things worked (seems like it definitely did, brocoli_tumor has mega high pre-B cells cluster 3)
heatmap(prop.table(table(flex$Treatment_Tumor_Ind,flex$ScType.cluster),margin = 1))







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
#   [1] openxlsx_4.2.7.1   HGNChelper_0.8.14  Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4           dplyr_1.1.4       
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                tibble_3.2.1                polyclip_1.10-7            
# [7] graph_1.80.0                XML_3.99-0.17               fastDummies_1.7.4           lifecycle_1.0.4             globals_0.16.3              lattice_0.22-6             
# [13] MASS_7.3-60.0.1             magrittr_2.0.3              plotly_4.10.4               httpuv_1.6.15               sctransform_0.4.1           zip_2.3.1                  
# [19] spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0           cowplot_1.1.3               pbapply_1.7-2               DBI_1.2.3                  
# [25] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1        purrr_1.0.2                
# [31] BiocGenerics_0.48.1         RCurl_1.98-1.16             GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6              
# [37] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-0        GSVA_1.50.5                 goftest_1.2-3               RSpectra_0.16-2            
# [43] spatstat.random_3.3-2       annotate_1.80.0             fitdistrplus_1.2-1          parallelly_1.38.0           DelayedMatrixStats_1.24.0   leiden_0.4.3.1             
# [49] codetools_0.2-20            DelayedArray_0.28.0         tidyselect_1.2.1            farver_2.1.2                ScaledMatrix_1.10.0         matrixStats_1.4.1          
# [55] stats4_4.3.1                spatstat.explore_3.3-2      jsonlite_1.8.9              progressr_0.14.0            ggridges_0.5.6              survival_3.7-0             
# [61] tools_4.3.1                 ica_1.0-3                   Rcpp_1.0.13                 glue_1.7.0                  gridExtra_2.3               SparseArray_1.2.4          
# [67] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.8         HDF5Array_1.30.1            withr_3.0.1                 fastmap_1.2.0               rhdf5filters_1.14.1        
# [73] fansi_1.0.6                 digest_0.6.37               rsvd_1.0.5                  R6_2.5.1                    mime_0.12                   colorspace_2.1-1           
# [79] scattermore_1.2             tensor_1.5                  spatstat.data_3.1-2         RSQLite_2.3.7               utf8_1.2.4                  tidyr_1.3.1                
# [85] generics_0.1.3              data.table_1.16.0           httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.2.1              uwot_0.2.2                 
# [91] pkgconfig_2.0.3             gtable_0.3.5                blob_1.2.4                  lmtest_0.9-40               SingleCellExperiment_1.24.0 XVector_0.42.0             
# [97] htmltools_0.5.8.1           dotCall64_1.1-1             GSEABase_1.64.0             scales_1.3.0                Biobase_2.62.0              png_0.1-8                  
# [103] spatstat.univar_3.0-1       rstudioapi_0.16.0           reshape2_1.4.4              nlme_3.1-166                cachem_1.1.0                zoo_1.8-12                 
# [109] rhdf5_2.46.1                stringr_1.5.1               KernSmooth_2.23-24          parallel_4.3.1              miniUI_0.1.1.1              AnnotationDbi_1.64.1       
# [115] pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                 RANN_2.6.2                  promises_1.3.0              BiocSingular_1.18.0        
# [121] beachmat_2.18.1             xtable_1.8-4                cluster_2.1.6               cli_3.6.3                   compiler_4.3.1              rlang_1.1.4                
# [127] crayon_1.5.3                future.apply_1.11.2         plyr_1.8.9                  stringi_1.8.4               viridisLite_0.4.2           deldir_2.0-4               
# [133] BiocParallel_1.36.0         munsell_0.5.1               Biostrings_2.70.3           lazyeval_0.2.2              spatstat.geom_3.3-3         Matrix_1.6-5               
# [139] RcppHNSW_0.6.0              patchwork_1.3.0             sparseMatrixStats_1.14.0    bit64_4.5.2                 future_1.34.0               ggplot2_3.5.1              
# [145] Rhdf5lib_1.24.2             KEGGREST_1.42.0             shiny_1.9.1                 SummarizedExperiment_1.32.0 ROCR_1.0-11                 igraph_2.0.3               
# [151] memoise_2.0.1               bit_4.5.0                   splitstackshape_1.4.8