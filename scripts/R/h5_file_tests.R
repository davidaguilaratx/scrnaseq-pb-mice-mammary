raw = Read10X_h5("D:/PhD/scRNAseq/8651-JC-mammary/10x_analysis_8651-JC/8651-JC-21/count/sample_raw_feature_bc_matrix.h5")
raw_seurat = CreateSeuratObject(raw, assay='RNA')
dim(raw_seurat)
# 32344 109502   genes x cells


fil = Read10X_h5("D:/PhD/scRNAseq/8651-JC-mammary/10x_analysis_8651-JC/8651-JC-21/count/sample_filtered_feature_bc_matrix.h5")
fil_seurat = CreateSeuratObject(fil, assay='RNA')
dim(fil_seurat)
# 19059  2142

fil2 = Read10X("D:/PhD/scRNAseq/8651-JC-mammary/10x_analysis_8651-JC/8651-JC-21/count/sample_filtered_feature_bc_matrix/")
fil2_seurat = CreateSeuratObject(fil2, assay='RNA')
dim(fil2_seurat)
# 19059  2142

rownames(fil2_seurat)

cb_raw = Read10X_h5("C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs_for_seurat/FPR_0.0/output_cellbender_8651-JC-21_FPR_0.0_seurat.h5")
cb_raw_seurat = CreateSeuratObject(cb_raw, assay='RNA')
dim(cb_raw_seurat)
# 32344 109502

cb_fil = Read10X_h5("C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs_for_seurat/FPR_0.0/output_cellbender_8651-JC-21_FPR_0.0_filtered_seurat.h5")
cb_fil_seurat = CreateSeuratObject(cb_fil, assay='RNA')
dim(cb_fil_seurat)
# 32344  3965

joined = Create_CellBender_Merged_Seurat(raw_cell_bender_matrix = cb_fil,
                                         raw_counts_matrix = raw)

joined$log10genes = log10(joined$nFeature_RNA)
joined$log10counts = log10(joined$nCount_RNA)

# mean and standard deviation of log10 ngenes and ncount
mean_ngene = mean(joined$log10genes)
mean_ncount = mean(joined$log10counts)
std_ngene = sd(joined$log10genes)
std_ncount = sd(joined$log10counts)

# z-scores of log10 ngenes and ncounts
joined$log10_genes_zscored = (joined$log10genes - mean_ngene) / std_ngene 
joined$log10_counts_zscored = (joined$log10counts - mean_ncount) / std_ncount

# Might as well do mitochondrial expression percentage here
joined = PercentageFeatureSet(joined, pattern = "^mt-", col.name = "percent.mt") # percent mitochondrial genes



# Add number of log10 genes per log10 UMI for each cell to metadata
# acts as novelty score to measure complexity
joined$log10GenesPerUMI = log10(joined$nFeature_RNA) / log10(joined$nCount_RNA)

# create metadata dataframe
metadata = joined@meta.data

# add a cell IDs column to metadata
metadata$cells = colnames(joined)

# save metadata dataframe back to merged seurat data
joined@meta.data = metadata


## proposed QC cutoffs ----
percent.mt.cutoff = 10
gene_cutoff_low = 200
gene_cutoff_high = 5000
umi_cutoff_low = 150
umi_cutoff_high = 50000 
complexity_cutoff = 0.8
min_cells_per_gene = 20 # filtered genes must be expressed in at least this # of cells



filtered_joined = subset(x = joined, 
                         subset=(nCount_RNA >= umi_cutoff_low) & 
                           (nCount_RNA <= umi_cutoff_high) & 
                           (nFeature_RNA >= gene_cutoff_low) &
                           (nFeature_RNA <= gene_cutoff_high) &
                           (log10GenesPerUMI > 0.80) & 
                           (percent.mt < percent.mt.cutoff))


cat(paste0(length(joined$orig.ident),' cells before cell-level filtering\n',
           length(filtered_joined$orig.ident),' cells left after cell-level filtering\n',
           length(joined$orig.ident) - length(filtered_joined$orig.ident),' cells removed after cell-level filtering'))


joined