library(goat)


nthreads <- 2 

# set seed
set.seed(2024)

# directories for saving
# selected_timepoint <- 'week3'
# selected_timepoint <- 'week10'
# selected_timepoint <- 'month7'
selected_timepoint <- 'month18'

de_directory <- "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/results/"

# directory to save figures and results in ----

# analysis type, cellbender or cellranger
analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analys

data_dir =  "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/GOAT/data/"

directory_figs = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/GOAT/figures/"
directory_results = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/GOAT/results/"

dirs = c(data_dir, directory_figs, directory_results)

for (d in dirs) {
  if (!dir.exists(d)) { dir.create(d,
                                   recursive=TRUE)}
}

### grab celltypes from DE results ----


### grab de result file names and extract cell types 
de_results = list.files(path=de_directory, pattern='all_genes')

### List of cell types (replace with your actual cell types)
celltypes = sub(pattern='_condition.*',replacement='',x=de_results)
celltypes



# load in data ----
res_tbl_adip = read.csv(paste0(de_directory,'Adipocytes_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_endo = read.csv(paste0(de_directory,'Endothelial_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_macro = read.csv(paste0(de_directory,'Macrophage.Ma_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_lum = read.csv(paste0(de_directory,'Luminal.HS_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_bcell = read.csv(paste0(de_directory,'B cells_condition_pb_vs_ctrl_all_genes_deseq2.csv'))

# to map our restuling genes to ensemble and entrez ids
gene_id = read.csv('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/annotation/ensembl_and_entrez_gene_ids_and_human_homologs.csv')

res_tbl_endo = get_ids_and_homologs(res_tbl_endo)
genelist_endo = names(to_genelist(res_tbl_endo, id='entrezgene_id'))
