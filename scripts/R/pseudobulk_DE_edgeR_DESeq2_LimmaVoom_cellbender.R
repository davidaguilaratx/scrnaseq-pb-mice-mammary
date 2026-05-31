# to install Matrix.utils
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)

# Load libraries
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(DESeq2)
library(limma)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(SeuratObject)
library(pheatmap)
library(apeglm)
# library(zinbwave)
library(png)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(EnhancedVolcano)
library(VennDiagram)

# set seed for reproducibility
set.seed(2024)

## Set a color-blind friendly palette
heat_colors <- rev(brewer.pal(11, "PuOr"))

# read in data ----

# directories for saving. Make sure they do not end in /
output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/figures/"
output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0/stats/figures/"
output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0/stats/figures/"
output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0/stats/figures/"

output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/data/"
output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0/data/"
output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0/data/"
output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0/data/"

# save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2/"
save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/"
save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/"
save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/"
save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/"

save_edger_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/edger_fltrd/"
save_edger_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0/stats/edger_fltrd/"
save_edger_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0/stats/edger_fltrd/"
save_edger_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0/stats/edger_fltrd/"

save_limma_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/limma_fltrd/"
save_limma_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/10_weeks/cellbender_analysis/FPR_0.0/stats/limma_fltrd/"
save_limma_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/7_months/cellbender_analysis/FPR_0.0/stats/limma_fltrd/"
save_limma_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/18_months/cellbender_analysis/FPR_0.0/stats/limma_fltrd/"



# Create directories to save results if they don't already exist:
dirs = c(output_data_dir, output_fig_dir, save_deseq2_dir, save_edger_dir, save_limma_dir)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# load in data
seurat = readRDS(paste0(output_data_dir,'de_seurat_annotated.rds'))

# directories for saving
save_deseq2_dir = 'stats/DESeq2_gene_fltrd/'
save_edger_dir = 'stats/edgeR/'
save_limma_dir = 'stats/limma/'

# criteria for filtering later
n_pb = length(grep('pb',unique(seurat$sample)))
n_ctrl = length(grep('ctrl',unique(seurat$sample)))
smallestGroupSize = ifelse(n_pb < n_ctrl, n_pb, n_ctrl)
cat('Size of smaller group is',smallestGroupSize,'\n There are',n_pb,
    'Pb samples and',n_ctrl,'controls.')
# Size of smaller group is 4 
# There are 4 Pb samples and 6 controls.

## How many cells per condition per cell type?
# Looks like not enough muscle or erythroid cells for analysis
table(seurat$celltype, seurat$condition)
# B cells                6992 7366
# CD4+ T cells           4425 3799
# CD8+ T cells           3609 3072
# Tregs                   827  826
# Endothelial            1157  843
# Fibroblast              880  595
# Adipocytes              987  628
# Proliferating cells     584  650
# Dendritic cell          592  584
# Activated CD4+ T cells  610  445
# Macrophage.Ma           492  423
# T helper 17 cells       433  348
# Myeloid.Mb              793  727
# Reticular cells         213  243
# Basal-Myoepithelial     237  163
# Luminal.HS              148   69
# Schwann cells           127   68
# Muscle                   67    0
# Erythroid cells          17   33

## How many cells will each pseudobulk sample have?
# looks like there won't be enough cells to analyze muscle or erythroid cells
table(seurat$celltype, seurat$sample)
# ctrl21 ctrl22 ctrl23 ctrl24 ctrl25 ctrl26 pb27 pb28 pb30 pb31
# B cells                   520    976   1204    460   2008   1824 3000 1364 2415  587
# CD4+ T cells              409    666    844    361   1039   1106 1599  712 1200  288
# CD8+ T cells              437    564    661    319    780    848 1163  523  984  402
# Tregs                      85    110    130     70    184    248  336  123  280   87
# Endothelial                97     55    272     65    392    276  284  315  197   47
# Fibroblast                 45     61    346     36    273    119  171  252  129   43
# Adipocytes                135     88    248     65    251    200  149  229  186   64
# Proliferating cells        67     81     84     42    125    185  229  125  185  111
# Dendritic cell             66     85     57     46    157    181  182   94  207  101
# Activated CD4+ T cells     56     79    133     59    116    167  154   86  155   50
# Macrophage.Ma              29     29    158     37    170     69   93  192  114   24
# T helper 17 cells          20     57     84     21    127    124  139   54  115   40
# Myeloid.Mb                 67    109    136     59    253    169  255  155  243   74
# Reticular cells             9     40     29     15     69     51   91   41   68   43
# Basal-Myoepithelial        32     14     41     15     78     57   45   75   31   12
# Luminal.HS                 55      7     33      6     36     11    9   28   22   10
# Schwann cells               4      4     33      5     57     24   32   22   13    1
# Muscle                      1      1     62      0      3      0    0    0    0    0
# Erythroid cells             1      2      5      2      4      3   31    1    1    0

# remove muscle cells from analysis since it's only made up of ctrl samples. No condition to compare to
seurat = subset(seurat, idents=c('Muscle','Erythroid cells'), invert=TRUE) # removed 67 muscle cells from data
seurat$celltype = Idents(seurat) # reassign with muscle and erythroid cluster removed

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- seurat@assays$RNA$counts

metadata <- seurat@meta.data

# Set up metadata as desired for aggregation and DE analysis
metadata$cluster_id <- factor(seurat@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

## Check the assays present
assays(sce)
# List of length 1
# names(1): counts

## Check the counts matrix
dim(counts(sce))
# [1] 13249 43955

counts(sce)[1:6, 1:6]

# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(sce)$celltype)
cluster_names
# [1] "B cells"                "CD4+ T cells"           "CD8+ T cells"           "Tregs"                  "Endothelial"            "Fibroblast"             "Adipocytes"            
# [8] "Proliferating cells"    "Dendritic cell"         "Activated CD4+ T cells" "Macrophage.Ma"          "T helper 17 cells"      "Myeloid.Mb"             "Reticular cells"       
# [15] "Basal-Myoepithelial"    "Luminal.HS"             "Schwann cells"  

# Total number of clusters
length(cluster_names)
# [1] 17


# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- unique(colData(sce)$sample)
sample_names
# [1] "ctrl21" "ctrl22" "ctrl23" "ctrl24" "ctrl25" "ctrl26" "pb27"   "pb28"   "pb30"  
# [10] "pb31"  

# Total number of samples
length(sample_names)
# [1] 10

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("celltype", "sample")]
groups$sample = as.factor(groups$sample) # coerce sample to factor data type

dim(groups)
# [1] 43955     2

head(groups)


# Aggregate across cluster-sample groups ----
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 

# Explore output matrix
class(aggr_counts)
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"

dim(aggr_counts)
# [1]   170 13249

aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

dim(aggr_counts)
# [1] 13249   179

# Understanding tstrsplit()

# Exploring structure of function output (list) ----
tstrsplit(colnames(aggr_counts), "_") %>% str()

## Comparing the first 10 elements of our input and output strings
head(colnames(aggr_counts), n = 10)
head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 10)

# Using which() to look up tstrsplit() output to isolate celltype. This is just an example
# here we isolate B cells
b_cell_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == "B cells")
b_cell_idx

colnames(aggr_counts)[b_cell_idx]
aggr_counts[1:10, b_cell_idx]

# Loop over all cell types to extract corresponding counts, and store information in a list ----

## Initiate empty list
counts_ls <- list()

for (i in 1:length(cluster_names)) {
  
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
  
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
  
}

# Explore the different components of the list
str(counts_ls)

length(counts_ls)
# [1] 17

# Reminder: explore structure of metadata ----
head(colData(sce))

# Extract sample-level variables
metadata <- colData(sce) %>% 
  as.data.frame() %>% 
  dplyr::select(condition,  sample)

dim(metadata)
# [1] 44005     2
head(metadata)

rownames(metadata) # barcodes

# Exclude duplicated rows, reducing metadata to 10 rows, one for each sample
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
# [1] 10  2

head(metadata)

# Rename rows to sample name
rownames(metadata) <- metadata$sample
head(metadata)

# Createa table with the number of cells per sample and cluster
t <- table(colData(sce)$sample,
           colData(sce)$cluster_id)
t[1:6, 1:6]

# Creating metadata list by appending cell count info to metadata table ----
# generates one metadata data frame specific of each cell type

## Initiate empty list
metadata_ls <- list()

for (i in 1:length(counts_ls)) {
  
  ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  
  ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  df$sample  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  
  ## Retrieve cell count information for this cluster from global cell count table
  idx <- which(colnames(t) == unique(df$cluster_id))
  cell_counts <- t[, idx]
  
  ## Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]
  
  ## Match order of cell_counts and sample_ids
  sample_order <- match(df$sample, names(cell_counts))
  cell_counts <- cell_counts[sample_order]
  
  ## Append cell_counts to data frame
  df$cell_count <- cell_counts
  
  
  ## Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- plyr::join(df, metadata, 
                   by = intersect(names(df), names(metadata)))
  
  ## Update rownames of metadata to match colnames of count matrix, as needed later for DE using DeSeq2
  rownames(df) <- df$cluster_sample_id
  
  ## Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- unique(df$cluster_id)
  
}

# Explore the different components of the list
str(metadata_ls)

# Double-check that both lists have same names
all(names(counts_ls) == names(metadata_ls))


# DESeq2 analysis ----
## first run of deseq2 using B cells as practice----
# idx <- which(names(counts_ls) == "B cells")
# cluster_counts <- counts_ls[[idx]]
# cluster_metadata <- metadata_ls[[idx]]
# 
# # Check contents of extracted objects
# cluster_counts[1:6, 1:6]
# head(cluster_metadata)
# 
# 
# # Check matching of matrix columns and metadata rows
# all(colnames(cluster_counts) == rownames(cluster_metadata))
# 
# 
# ### Create DESeq2 object ----   
# dds <- DESeqDataSetFromMatrix(cluster_counts, 
#                               colData = cluster_metadata, 
#                               design = ~ condition)
# 
# # Transform counts for data visualization
# rld <- rlog(dds, blind=TRUE)
# # Plot PCA
# DESeq2::plotPCA(rld, ntop = 500, intgroup = "condition")
# ggsave("figures/overall_specific_PCAplot.png",dpi=320, height=10,width=10)
# 
# DESeq2::plotPCA(rld, ntop = 500, intgroup = "cell_count")
# 
# 
# # Extract the rlog matrix from the object and compute pairwise correlation values
# rld_mat <- assay(rld)
# rld_cor <- cor(rld_mat)
# 
# # Plot heatmap
# pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop=F])
# 
# 
# ### Run DESeq2 differential expression analysis ----
# dds <- DESeq(dds)
# 
# ### Plot dispersion estimates ----
# plotDispEsts(dds)
# 
# # Check the coefficients for the comparison
# resultsNames(dds)
# 
# # Generate results object
# res <- results(dds, 
#                name = "condition_pb_vs_ctrl",
#                alpha = 0.05)
# 
# # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
# res <- lfcShrink(dds, 
#                  coef = "condition_pb_vs_ctrl",
#                  res=res,
#                  type = "apeglm")
# 
# ### table DE results ----
# # Turn the DESeq2 results object into a tibble for use with tidyverse functions
# res_tbl <- res %>%
#   data.frame() %>%
#   rownames_to_column(var = "gene") %>%
#   as_tibble() %>%
#   arrange(padj)
# 
# # Check results output
# res_tbl 
# 
# # Write all results to file
# write.csv(res_tbl,
#           paste0("stats/DESeq2/results/", unique(cluster_metadata$cluster_id), "_", 
#                  unique(cluster_metadata$condition)[2], "_vs_", unique(cluster_metadata$condition)[1], "_all_genes_deseq2.csv"),
#           quote = FALSE,
#           row.names = FALSE)
# 
# 
# ### table DE results for significant genes only ----
# # Set thresholds
# padj_cutoff <- 0.05
# 
# # Subset the significant results
# sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
#   dplyr::arrange(padj)
# 
# # Check significant genes output
# sig_res
# 
# # Write significant results to file
# write.csv(sig_res,
#           paste0("stats/DESeq2/results/", unique(cluster_metadata$cluster_id), "_", 
#                  unique(cluster_metadata$condition)[2], "_vs_", unique(cluster_metadata$condition)[1], "_signif_genes_deseq2.csv"),
#           quote = FALSE, 
#           row.names = FALSE)
# 
# # Set thresholds
# log2fc_cutoff <- 0.58
# 
# # Count significantly up/down genes above threshold
# n_sig_up <- dplyr::filter(sig_res, log2FoldChange >= log2fc_cutoff) %>% 
#   nrow()
# n_sig_dn <- dplyr::filter(sig_res, log2FoldChange <= -log2fc_cutoff) %>% 
#   nrow()
# n_sig_up
# n_sig_dn
# 
# ### Scatterplot of top genes for celltype----
# ## Extract normalized counts from dds object
# normalized_counts <- counts(dds, normalized = TRUE)
# 
# ## Extract top 20 DEG from resLFC (make sure to order by padj)
# top20_sig_genes <- sig_res %>%
#   dplyr::arrange(padj) %>%
#   dplyr::pull(gene) %>%
#   head(n = min(20,nrow(sig_res)))
# 
# length(top20_sig_genes)
# head(top20_sig_genes)
# 
# ## Extract matching normalized count values from matrix
# top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
# top20_sig_counts
# 
# ## Convert wide matrix to long data frame for ggplot2
# top20_sig_df <- data.frame(top20_sig_counts)
# top20_sig_df$gene <- rownames(top20_sig_counts)
# 
# top20_sig_df <- melt(setDT(top20_sig_df), 
#                      id.vars = c("gene"),
#                      variable.name = "cluster_sample_id") %>% 
#   data.frame()
# 
# ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
# top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
# 
# ## Join counts data frame with metadata
# top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
#                            by = "cluster_sample_id")
# head(top20_sig_df)
# 
# ## Generate plot
# ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
#   geom_boxplot(outlier.shape=NA, width=0.2)+
#   geom_jitter(height = 0, width = 0.5)+
#   scale_y_continuous(trans = 'log10') +
#   ylab("log10 of normalized expression level") +
#   xlab("condition") +
#   ggtitle("Top 17 Significant DE Genes overall") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   facet_wrap(~ gene)
# ggsave('stats/DESeq2/figures/overall_top17DEgenes.png',dpi=320,width=20,height=10)
# 
# ### Heatmap of overall significant genes ----
# 
# ## Extract normalized counts for significant genes only
# sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
# 
# ## Run pheatmap using the metadata data frame for the annotation
# png(filename='stats/DESeq2/figures/overall_pheatmap_top17DEgene.png', units='in', width=10, height=10, res=320)
# pheatmap(sig_counts, 
#          color = heat_colors, 
#          cluster_rows = TRUE, 
#          show_rownames = TRUE,
#          annotation_col = cluster_metadata[,c('condition','sample')],
#          border_color = NA, 
#          fontsize = 10, 
#          scale = "row", 
#          fontsize_row = 10, 
#          height = 20)
# dev.off()
# graphics.off()
# 
# # Volcano plot
# res_table_thres <- res_tbl[!is.na(res_tbl$padj), ] %>% 
#   mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
# min(log10(res_table_thres$padj))
# 
# ## Generate volcano plot
# ggplot(res_table_thres) +
#   geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
#   ggtitle("Volcano plot of pb-exposed mammary gland cells relative to control") +
#   xlab("log2 fold change") +
#   xlim(-1.5, 1.5) +
#   ylab("-log10 adjusted p-value") +
#   scale_y_continuous(limits = c(0, 4)) +
#   scale_color_manual(values = c("grey60", "red3")) +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.3), hjust = 0.5),
#         axis.title = element_text(size = rel(1.15)))  
# 
# ## generate enhanced volcano plot
# EnhancedVolcano(res_tbl,
#                 lab=res_tbl$gene,
#                 title='Overall DE Genes',
#                 subtitle='Pb vs control',
#                 x='log2FoldChange',
#                 y='padj',
#                 xlim=c(-1.5,1.5),
#                 ylim=c(0,4),
#                 pCutoff=0.05,
#                 FCcutoff=0.58, # roughly 50% change
#                 pointSize=4,
#                 labSize=8,
#                 colAlpha=0.8,
#                 cutoffLineType='twodash',
#                 cutoffLineWidth=0.5,
#                 legendLabSize=16,
#                 legendIconSize = 5)
# ggsave('stats/DEseq2/figures/overall_DEgenes_evolcano.png',dpi=320,width=15,height=10)
# 





### run pseudobulk DESeq2 aggregating all cells to get overall differences by condition ----
# Subset metadata to include only the variables you want to 
# aggregate across (here, we want to aggregate by sample and by cluster)

# Create directories to save results if they don't already exist:
if (!dir.exists(paste0(save_deseq2_dir,'/figures'))) { dir.create(paste0(save_deseq2_dir,'/figures'),
                                                      recursive=TRUE) }
if (!dir.exists(paste0(save_deseq2_dir,'/results'))) { dir.create(paste0(save_deseq2_dir,'/results'),
                                                                  recursive=TRUE) }

groups_bulk <- colData(sce)[, c("sample")]
groups_bulk = as.factor(groups_bulk) # coerce sample to factor data type

length(groups_bulk)
# [1] 43955     

head(groups_bulk)

#### Aggregate across cluster-sample groups ----
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts_bulk <- aggregate.Matrix(t(counts(sce)), 
                                     groupings = groups_bulk, fun = "sum") 

# Explore output matrix
class(aggr_counts_bulk)
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"

dim(aggr_counts_bulk)
# [1]   10 13249

# How many reads per pseudobulk sample?
rowSums(aggr_counts_bulk)
# ctrl21  ctrl22  ctrl23  ctrl24  ctrl25  ctrl26    pb27    pb28    pb30    pb31 
# 1882676 2206272 3705188 1801452 3046964 4670048 4882792 2363042 4383457 1984282 

# How many reads per pseudobulk sample?

aggr_counts_bulk[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts_bulk <- t(aggr_counts_bulk)
aggr_counts_bulk[1:6, 1:6]

dim(aggr_counts_bulk)
# [1] 13249    10



# create copy of metadata and add cell count col
metadata_bulk = copy(metadata)
metadata_bulk$cell_count = as.numeric(table(seurat$sample))

# Check matching of matrix columns and metadata rows
all(colnames(aggr_counts_bulk) == rownames(metadata_bulk))

# zero-count proportion?
(sum(rowSums(aggr_counts_bulk == 0))) / (nrow(aggr_counts_bulk) * ncol(aggr_counts_bulk))
# [1] 0.01206129

# Create DESeq2 object        
dds_bulk <- DESeqDataSetFromMatrix(aggr_counts_bulk, 
                                   colData = metadata_bulk, 
                                   design = ~ condition)

# prefiltering low count genes. 
# Want genes that have at least 10 counts in at least 4 samples.
# following recommendation from
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
# size of smaller group. there are 4 Pb samples and 6 controls
keep = rowSums(counts(dds_bulk) >= 10) >= smallestGroupSize
table(keep)
# keep
# FALSE  TRUE 
# 1614 11635 discarding 1614 genes
dds_bulk = dds_bulk[keep,]

# zinb_bulk = zinbwave(dds_bulk, K=0, observationalWeights=TRUE,
#                      BPPARAM=BiocParallel::SerialParam(), epsilon=1e12)


# Transform counts for data visualization
rld_bulk <- rlog(dds_bulk, blind=TRUE)
# Plot PCA
DESeq2::plotPCA(rld_bulk, ntop = 500, intgroup = "condition")
ggsave(paste0(save_deseq2_dir,"/figures/overall_specific_PCAplot.png"),dpi=320,height=10,width=10)

DESeq2::plotPCA(rld_bulk, ntop = 500, intgroup = "cell_count")
ggsave(paste0(save_deseq2_dir,"/figures/overall_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat_bulk <- assay(rld_bulk)
rld_cor_bulk <- cor(rld_mat_bulk)

# Plot heatmap
png(paste0(save_deseq2_dir,'/figures/overall_specific_heatmap.png'),
    res=320, height = 10, width = 10, units = "in")
pheatmap(rld_cor_bulk, annotation = metadata_bulk[, c("condition"), drop=F])
dev.off()
dev.off()


# compute size factors using scran::computeSumFactors() as recommended for
# single cell analysis on deseq2 vignette here:
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
# scr_bulk = scran::computeSumFactors(dds_bulk) # scran size factors
# 
# sizeFactors(dds_bulk) = sizeFactors(scr_bulk)

#### Run DESeq2 differential expression analysis ----
dds_bulk <- DESeq(dds_bulk)

#### Plot dispersion estimates ----
png(paste0(save_deseq2_dir,'/figures/overall_dispersion_plot.png'),
    res=320, height = 5, width = 6, units = "in")
plotDispEsts(dds_bulk)
dev.off()

# Check the coefficients for the comparison
resultsNames(dds_bulk)

# Generate results object
res_bulk <- results(dds_bulk, 
                    name = "condition_pb_vs_ctrl",
                    alpha = 0.05)

# Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
res_bulk <- lfcShrink(dds_bulk, 
                      coef = "condition_pb_vs_ctrl",
                      res=res_bulk,
                      type = "apeglm")

#### table DE results ----
# Turn the DESeq2 results object into a tibble for use with tidyverse functions
res_tbl_bulk <- res_bulk %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Check results output
res_tbl_bulk 

# Write all results to file
write.csv(res_tbl_bulk,
          paste0(save_deseq2_dir,"/results/overall_condition_", 
                 unique(metadata_bulk$condition)[2], "_vs_", unique(metadata_bulk$condition)[1], "_all_genes_deseq2.csv"),
          quote = FALSE,
          row.names = FALSE)


#### table DE results for significant genes only ----
# Set thresholds
padj_cutoff <- 0.05

# Subset the significant results
sig_res_bulk <- dplyr::filter(res_tbl_bulk, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res_bulk

# Write significant results to file
write.csv(sig_res_bulk,
          paste0(save_deseq2_dir,"/results/overall_condition_", 
                 unique(metadata_bulk$condition)[2], "_vs_", unique(metadata_bulk$condition)[1], "_signif_genes_deseq2.csv"),
          quote = FALSE, 
          row.names = FALSE)

# Set thresholds
log2fc_cutoff <- 0.58

# Count significantly up/down genes above threshold
n_sig_up <- dplyr::filter(sig_res_bulk, log2FoldChange >= log2fc_cutoff) %>% 
  nrow()
n_sig_dn <- dplyr::filter(sig_res_bulk, log2FoldChange <= -log2fc_cutoff) %>% 
  nrow()


#### Scatterplot of top genes across all cells ----

## Extract normalized counts from dds object
normalized_counts <- counts(dds_bulk, normalized = TRUE)

## Extract top 20 DEG from resLFC (make sure to order by padj)
top20_sig_genes <- sig_res_bulk %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n = min(20,nrow(sig_res_bulk)))

## Extract matching normalized count values from matrix
top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
top20_sig_counts

## Convert wide matrix to long data frame for ggplot2
top20_sig_df <- data.frame(top20_sig_counts)
top20_sig_df$gene <- rownames(top20_sig_counts)

top20_sig_df <- melt(setDT(top20_sig_df), 
                     id.vars = c("gene"),
                     variable.name = "sample") %>% 
  data.frame()


## Join counts data frame with metadata
top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds_bulk)),
                           by = "sample")
head(top20_sig_df)

## Generate plot
ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
  geom_boxplot(outlier.shape=NA, width=0.2)+
  geom_jitter(height = 0, width = 0.5)+
  scale_y_continuous(trans = 'log10') +
  ylab("log10 of normalized expression level") +
  xlab("condition") +
  ggtitle(paste0("Top ", min(20,nrow(sig_res_bulk)), " Significant DE Genes overall")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene)
ggsave(paste0(save_deseq2_dir,'/figures/overall_top17DEgenes.png'),dpi=320,width=20,height=10)

#### Heatmap of overall significant genes ----

## Extract normalized counts for significant genes only
sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res_bulk$gene, ]

## Run pheatmap using the metadata data frame for the annotation
png(filename=paste0(save_deseq2_dir,'/figures/overall_pheatmap_top17DEgene.png'), units='in', width=10, height=10, res=320)
pheatmap(sig_counts, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         annotation_col = metadata_bulk[,c('condition','sample')],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)
dev.off()

# Volcano plot
res_table_thres <- res_tbl_bulk[!is.na(res_tbl_bulk$padj), ] %>% 
  mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
min(log10(res_table_thres$padj))

## Generate volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of pb-exposed mammary gland cells relative to control") +
  xlab("log2 fold change") +
  xlim(-1.5, 1.5) +
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0, 4)) +
  scale_color_manual(values = c("grey60", "red3")) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.3), hjust = 0.5),
        axis.title = element_text(size = rel(1.15)))  

## generate enhanced volcano plot
EnhancedVolcano(res_tbl_bulk,
                lab=res_tbl_bulk$gene,
                title='Overall DE Genes',
                subtitle='Pb vs control mammary gland',
                x='log2FoldChange',
                y='padj',
                xlim=c(-2,2),
                ylim=c(0,6),
                pCutoff=0.05,
                FCcutoff=0.5, 
                pointSize=4,
                labSize=8,
                colAlpha=0.8,
                cutoffLineType='twodash',
                cutoffLineWidth=0.5,
                legendLabSize=16,
                legendIconSize = 5,
                drawConnectors=TRUE)
ggsave(paste0(save_deseq2_dir,'/figures/overall_DEgenes_evolcano.png'),dpi=320,width=15,height=10)





### DE genes by celltype using DESeq2 ----

# Function to run DESeq2 Wald Test and get results for any cluster:
## clustx is the name of the cluster (cell type) on which to run the function
## A is the sample group to compare (e.g. stimulated condition)
## B is the sample group to compare against (base/control level)
## padj_cutoff defines the ajusted p-value cutoff for significance (set to 0.05 by default)

## This function assumes the counts matrices and metadata for all clusters have been prepared
## and arranged in matching named lists (as illustrated in tutorial above)
## This function assumes the contrast (e.g. stim vs. control) is stored in a variable named "condition"

#### keep track of zero count proportions for each cell type ----
zero_props = list()
zero_props_post_filtering = list()
change_in_zero_props = list()

get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
  
  cat('Beggining cluster ',clustx,'\n\n') # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  zero_prop = (sum(rowSums(cluster_counts == 0))) / (nrow(cluster_counts) * ncol(cluster_counts))
  cat('Proportion of zero counts in',clustx,'for 13249 genes across 10 samples:\n')
  print(zero_prop)
  
  zero_props[[clustx]] = zero_prop
  
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ condition)
  
  # prefiltering low count genes. 
  # Want genes that have at least 10 counts in at least 4 samples.
  # following recommendation from
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  smallestGroupSize = 4 # size of smaller group. there are 4 Pb samples
  keep = rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds = dds[keep,]
  
  
  zero_prop_post_filter = (sum(rowSums(dds@assays@data$counts == 0))) / 
    (nrow(dds@assays@data$counts) * ncol(dds@assays@data$counts))
  
  cat('Proportion of zero counts in',clustx,'for 13249 genes across 10 samples
      after low-count gene filtering:\n')
  print(zero_prop_post_filter)
  
  zero_props_post_filtering[[clustx]] = zero_prop_post_filter
  
  change_zero_prop = zero_prop - zero_prop_post_filter
  cat('Change in zero count proportions after filtering for',clustx,':')
  print(change_zero_prop)
  
  change_in_zero_props[[clustx]] = change_zero_prop
  
  # # zinb-wave used to intoduce weights for zero-count proportions
  # zinb = zinbwave(dds, K=0, observationalWeights=TRUE,
  #                 BPPARAM=BiocParallel::SerialParam(), epsilon=1e12)
  # 
  # dds = DESeqDataSet(zinb, design=~condition)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)
  
  # Generate QC plots
  
  ## Plot and save PCA plot
  DESeq2::plotPCA(rld, intgroup = "condition")
  ggsave(paste0(save_deseq2_dir,"/figures/", clustx, "_specific_PCAplot.png"))
  
  DESeq2::plotPCA(rld, intgroup = "cell_count")
  ggsave(paste0(save_deseq2_dir,"/figures/", clustx, "_specific_PCAplot_cellcounts.png"))
  
  ## Extract rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  ## Plot and save heatmap
  png(paste0(save_deseq2_dir,"/figures/", clustx, "_specific_heatmap.png"),
      height = 6, width = 7.5, units = "in", res = 300)
  pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop = FALSE])
  dev.off()
  
  # compute size factors using scran::computeSumFactors() as recommended for
  # single cell analysis on deseq2 vignette here:
  # https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
  scr = scran::computeSumFactors(dds) # scran size factors
  
  sizeFactors(dds) = sizeFactors(scr)
  
  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds,
               test='LRT',
               fitType='parametric', # the default
               reduced= ~ 1,
               minReplicatesForReplace=Inf,
               useT=TRUE,
               minmu=1e-6)
  
  cat('Finished fitting DESeq2 model for',clustx,'\n')
  
  ## Plot dispersion estimates
  png(paste0(save_deseq2_dir,"/figures/", clustx, "_dispersion_plot.png"),
      height = 5, width = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
  
  ## Output and shrink results of Wald test for contrast A vs B
  contrast <- paste(c("condition", A, "vs", B), collapse = "_")
  print(resultsNames(dds))
  
  res <- results(dds, name = contrast, alpha = 0.05, independentFiltering = FALSE)
  # res <- results(dds, name = 'conditionpb', alpha = 0.05)
  res <- lfcShrink(dds, coef = contrast, res = res, type = "apeglm")
  # res <- lfcShrink(dds, coef = 'conditionpb', res = res, type = "apeglm")
  
  ## Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  write.csv(res_tbl,
            paste0(save_deseq2_dir,"/results/", clustx, "_", contrast, "_all_genes_deseq2.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  ## Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            paste0(save_deseq2_dir,"/results/", clustx, "_", contrast, "_signif_genes_deseq2.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Generate results visualization plots
  
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds, normalized = TRUE)
  
  if ((nrow(sig_res) > 0) & (!is.null(nrow(sig_res)))) {
    ## Extract top 20 DEG from resLFC (make sure to order by padj)
    n = min(20, nrow(sig_res)) # top 20 or fewer DE genes
    top20_sig_genes <- sig_res %>%
      dplyr::arrange(padj) %>%
      dplyr::pull(gene) %>%
      head(n = n)
    
    print(head(top20_sig_genes))
    
    ## Convert wide matrix to long data frame for ggplot2
    
    ## Extract matching normalized count values from matrix
    top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
    
    # Convert the matrix to a data frame and add the 'gene' column
    top20_sig_df <- data.frame(top20_sig_counts)
    top20_sig_df$gene <- rownames(top20_sig_counts)
    
    # Convert to data.table and melt into long format
    top20_sig_df <- melt(setDT(top20_sig_df), 
                         id.vars = c("gene"),
                         variable.name = "cluster_sample_id")
    
    ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
    top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
    
    # Proceed with joining the metadata and plotting
    top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
                               by = "cluster_sample_id")
    
    ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
      geom_boxplot(outlier.shape=NA, width=0.2)+
      geom_jitter(height = 0, width = 0.5)+
      scale_y_continuous(trans = 'log10') +
      ylab("log10 of normalized expression level") +
      xlab("condition") +
      ggtitle(paste0("Top ", n, " Significant DE Genes")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ gene)
    ggsave(paste0(save_deseq2_dir,"/figures/", clustx, "_", contrast, "_top20_DE_genes.png"))
    
    ## Extract normalized counts for significant genes only
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
    
    ## generate pheatmap
    png(filename=paste0(save_deseq2_dir,'/figures/',clustx,'_',contrast,'_pheatmap_top17DEgene.png'), units='in', width=10, height=10, res=320)
    pheatmap(sig_counts, 
             color = heat_colors, 
             cluster_rows = TRUE, 
             show_rownames = FALSE,
             annotation_col = cluster_metadata[,c('condition','sample')],
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 10, 
             height = 20)
    dev.off()
    
  } else {
    warning(paste0("No significant genes found for ",clustx))
    cat(paste0("No significant genes found for ",clustx),'\n\n')
  }
  
  
  ## generate enhanced volcano plot
  EnhancedVolcano(res_tbl,
                  lab=res_tbl$gene,
                  title=paste0(clustx,' DE Genes'),
                  subtitle='Pb vs control',
                  x='log2FoldChange',
                  y='padj',
                  pCutoff=0.05,
                  FCcutoff=0.5, # roughly 50% change equates to 0.58 log2FC, or 1.5 FC
                  pointSize=4,
                  labSize=8,
                  colAlpha=0.8,
                  cutoffLineType='twodash',
                  cutoffLineWidth=0.5,
                  legendLabSize=16,
                  legendIconSize = 5,
                  drawConnectors=TRUE)
  ggsave(paste0(save_deseq2_dir,'/figures/',clustx,'_',contrast,'_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
  
  cat('\nFinished analyzing',clustx,'\n')
}

# Run the script on all clusters comparing stimulated condition relative to control condition
map(cluster_names, get_dds_resultsAvsB, A = "pb", B = "ctrl", padj_cutoff = 0.05)





# edgeR analysis ----
# Create DGEList object for use in edgeR. practice with B cells
dge = DGEList(counts=counts_ls$`B cells`, samples=metadata_ls$`B cells`)
head(dge)

# remove label-sample combinations with fewer than 10 cells
discarded = dge$samples$cell_count < 10
dge = dge[,!discarded]
summary(discarded)
# Mode   FALSE 
# logical      10 
# none discarded

unique(seurat[,seurat$celltype == 'Erythroid cells']$sample) # no erythroid cells in pb31

# remove genes that are not expressed above a log-CPM threshold in a minimum number
# of samples (determined from the size of the smallest treatment group in 
# the experimental design, or pb group in this case.)
keep = filterByExpr(dge, group=dge$samples$condition,
                    min.count=10, min.total.count=15, min.prop=0.4)
dge = dge[keep,]
summary(keep)
# Mode   FALSE    TRUE 
# logical    7770    5479 


# correct for composition biases by computing normalization factors with the 
# trimmed mean of M-values method (Robinson and Oshlack 2010), convertering raw library
# sizes in normalized effective library sizes
dge = calcNormFactors(dge)
head(dge$samples)

# Create directories to save results if they don't already exist:
if (!dir.exists(paste0(save_edger_dir,"/figures"))) { dir.create(paste0(save_edger_dir,"/figures"),
                                                                 recursive=TRUE) }
if (!dir.exists(paste0(save_edger_dir,"/results"))) { dir.create(paste0(save_edger_dir,"/results"),
                                                                 recursive=TRUE) }

# Plot mean-difference plot for each pseudobulk profile
png(filename=paste0(save_edger_dir,'/figures/',clustx,'_mean-difference_plots.png'), units='in', width=10, height=10, res=320)
par(mfrow=c(3,4))
for (i in seq_len(ncol(dge))) {
  plotMD(dge, column=i)
}
mtext(paste0('Mean-difference plots for pseudobulk samples - ', unique(dge$samples$cluster_id)),
      outer=T,cex=1.15, side='top', adj=0, padj=1)
dev.off()

# Plot multi-dimensional scaling plots for each pseudobulk profile
png(filename=paste0(save_edger_dir,'/figures/',unique(dge$samples$cluster_id),'_multidimensional_scaling_plot.png'), units='in', width=10, height=10, res=320)
plotMDS(edgeR::cpm(dge, log=TRUE), 
        col=ifelse(dge$samples$condition == 'pb', "red", "blue"))
title(paste0('MDS plot for pseudobulk samples - ', unique(dge$samples$cluster_id)))
dev.off()

# set up design matrix for edgeR
design = model.matrix(~factor(condition), dge$samples)

# estimate negative bibnomial dispersion
dge = estimateDisp(dge, design)
summary(dge$trended.dispersion)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.009666 0.011287 0.013003 0.015265 0.016352 0.032906 

# plot biological coefficient of variation for each gene as a function of average abundance
# The BCV is computed as the square root of the NB dispersion after 
# empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are 
# shown in blue and red, respectively.
png(filename=paste0(save_edger_dir,'/figures/',unique(dge$samples$cluster_id),'_BCV_plot.png'), units='in', width=10, height=10, res=320)
plotBCV(dge)
title(paste0('BCV for ',unique(dge$samples$cluster_id), ' genes'))
dev.off()

# Fit a GLM to the counts for each gene and estimates the QL dispersion from the 
# GLM deviance. We set robust=TRUE to avoid distortions from 
# highly variable clusters (Phipson et al. 2016). The QL dispersion models the 
# uncertainty and variability of the per-gene variance - which is not 
# well handled by the NB dispersions, so the two dispersion types 
# complement each other in the final analysis.
fit = glmQLFit(dge, design, robust=TRUE)
summary(fit$var.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5798  0.8138  0.8780  0.8610  0.9201  1.0786 

summary(fit$df.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2738 27.1197 27.1197 25.7818 27.1197 27.1197 

# Plot QL dispersion estimates for each gene as a function of abundance.
# Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).
png(filename=paste0(save_edger_dir,'/figures/',unique(dge$samples$cluster_id),'_QLdispersion_estimates.png'), units='in', width=10, height=10, res=320)
plotQLDisp(fit)
fittitle(paste0('QL dispersion estimates for ',unique(dge$samples$cluster_id), ' genes'))
dev.off()

# Test for differences in expression due to condition, Pb exposure
result = glmQLFTest(fit, coef=ncol(design))
summary(decideTests(result))
# factor(condition)pb
# Down                     2
# NotSig                8016
# Up                       0

topTags(result)
# Coefficient:  factor(condition)pb 
# logFC   logCPM        F       PValue         FDR
# Stip1   -0.7140595 7.318869 39.45676 3.246866e-07 0.002603337
# Frat2   -1.0573577 5.004866 27.45923 7.685404e-06 0.030810785
# Cfd     -1.6839362 3.817630 18.59192 1.248610e-04 0.244284705
# Adck2    0.8289027 4.970495 18.36890 1.347894e-04 0.244284705
# Slc14a1  0.7624198 5.610101 18.01437 1.523352e-04 0.244284705
# Picalm   0.4540081 7.016770 16.93063 2.227258e-04 0.297635950
# Narf    -0.5097320 7.608240 16.86363 2.861539e-04 0.309910485
# Tmem238 -1.1407662 4.323068 16.01478 3.092148e-04 0.309910485
# Gtf3c3   0.8315441 4.929941 15.58878 3.610260e-04 0.317912243
# St13    -0.6393512 6.403725 15.64738 4.248224e-04 0.317912243



## overall bulk analysys by condition with edgeR ----
# Create DGEList object for use in edgeR
# Create directories to save results if they don't already exist:
if (!dir.exists(paste0(save_edger_dir,"/figures"))) { dir.create(paste0(save_edger_dir,"/figures"),
                                                                 recursive=TRUE) }
if (!dir.exists(paste0(save_edger_dir,"/results"))) { dir.create(paste0(save_edger_dir,"/results"),
                                                                 recursive=TRUE) }

dge_bulk = DGEList(counts=aggr_counts_bulk, samples=metadata_bulk)
head(dge_bulk)


# remove label-sample combinations with fewer than 10 cells
discarded_bulk = dge_bulk$samples$cell_count <= 10
dge_bulk = dge_bulk[,!discarded_bulk]
summary(discarded_bulk)
# Mode   FALSE 
# logical      10 
# none discarded


# # remove genes that are not expressed above a log-CPM threshold in a minimum number
# # of samples (determined from the size of the smallest treatment group in 
# # the experimental design, or pb group in this case.)
# keep_bulk = filterByExpr(dge_bulk, group=dge_bulk$samples$condition,
#                          min.count=10, min.total.count=15, min.prop=0.7) # default
# dge_bulk = dge_bulk[keep_bulk,]
# summary(keep_bulk)
# Mode   FALSE    TRUE
# logical    1867   11382

keep_bulk = rowSums(dge_bulk$counts >= 10) >= smallestGroupSize # smallest group is 4
dge_bulk = dge_bulk[keep_bulk, ]
summary(keep_bulk)
# Mode   FALSE    TRUE 
# logical    1607   11642 

# correct for composition biases by computing normalization factors with the 
# trimmed mean of M-values method (Robinson and Oshlack 2010), convertering raw library
# sizes in normalized effective library sizes
dge_bulk = calcNormFactors(dge_bulk)
dge_bulk$samples

# Plot mean-difference plot for each pseudobulk profile
png(filename=paste0(save_edger_dir,'/figures/overall_mean-difference_plots.png'), units='in', width=10, height=10, res=320)
par(mfrow=c(3,4))
for (i in seq_len(ncol(dge_bulk))) {
  plotMD(dge_bulk, column=i)
}
mtext(paste0('Mean-difference plots for pseudobulk samples - overall'),
      outer=T,cex=1.15, side='top', adj=0, padj=1)
dev.off()

# Plot multi-dimensional scaling plots for each pseudobulk profile
png(filename=paste0(save_edger_dir,'/figures/overall_multidimensional_scaling_plot.png'), units='in', width=10, height=10, res=320)
plotMDS(edgeR::cpm(dge_bulk, log=TRUE), 
        col=ifelse(dge_bulk$samples$condition == 'pb', "red", "blue"))
title(paste0('MDS plot for pseudobulk samples - overall'))
dev.off()

# set up design matrix for edgeR
design_bulk = model.matrix(~factor(condition), dge_bulk$samples)

# estimate negative bibnomial dispersion
dge_bulk = estimateDisp(dge_bulk, design_bulk)
summary(dge_bulk$trended.dispersion)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01593 0.01808 0.02892 0.04591 0.06649 0.12761 

# plot biological coefficient of variation for each gene as a function of average abundance
# The BCV is computed as the square root of the NB dispersion after 
# empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are 
# shown in blue and red, respectively.
png(filename=paste0(save_edger_dir,'/figures/overall_BCV_plot.png'), units='in', width=10, height=10, res=320)
plotBCV(dge_bulk)
title(paste0('BCV for overall genes'))
dev.off()

# Fit a GLM to the counts for each gene and estimates the QL dispersion from the 
# GLM deviance. We set robust=TRUE to avoid distortions from 
# highly variable clusters (Phipson et al. 2016). The QL dispersion models the 
# uncertainty and variability of the per-gene variance - which is not 
# well handled by the NB dispersions, so the two dispersion types 
# complement each other in the final analysis.
fit_bulk = glmQLFit(dge_bulk, design_bulk, robust=TRUE)
summary(fit_bulk$var.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8389  0.8983  1.1135  1.2870  1.7131  2.0317  

summary(fit_bulk$df.prior)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.895   5.428   5.428   5.186   5.428   5.428 

# Plot QL dispersion estimates for each gene as a function of abundance.
# Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).
png(filename=paste0(save_edger_dir,'/figures/overall_QLdispersion_estimates.png'), units='in', width=10, height=10, res=320)
plotQLDisp(fit_bulk)
title(paste0('QL dispersion estimates for overall genes'))
dev.off()

# Test for differences in expression due to condition, Pb exposure
result_bulk = glmQLFTest(fit_bulk, coef=ncol(design_bulk))
summary(decideTests(result_bulk))
# factor(condition)pb
# Down                     0
# NotSig               11624
# Up                       0
# no significantly expressed genes after multiple testing correction

res_bulk = topTags(result_bulk, n=Inf)
head(res_bulk, 10)
# Coefficient:  factor(condition)pb 
# logFC   logCPM        F       PValue       FDR
# Slc27a2  -3.5089015 3.260035 39.37358 2.176507e-05 0.2073778
# Stip1    -0.7127911 7.391454 36.44987 3.646203e-05 0.2073778
# Apoe     -0.6275331 8.691676 29.24379 1.067695e-04 0.3849811
# Spry1    -1.2044394 3.706505 26.24331 1.780773e-04 0.3849811
# Hsph1    -1.1617826 6.480592 28.92426 1.939593e-04 0.3849811
# Ctso      0.4755548 6.872607 24.35576 2.478778e-04 0.3849811
# Haao      0.4534229 7.146709 23.38232 2.971393e-04 0.3849811
# Timeless  0.6040163 5.536668 22.78974 3.327378e-04 0.3849811
# Cd180     0.5680880 6.390895 22.56116 3.476418e-04 0.3849811
# Tatdn1    0.3855929 6.246507 21.92553 3.935629e-04 0.3849811

# format data for saving to csv file
res_tbl_edger.bulk = res_bulk$table[order(res_bulk$table$FDR),] %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as_tibble()

# save to csv
write.csv(res_tbl_edger.bulk,
          paste0(save_edger_dir,'/results/Overall_condition_pb_vs_ctrl_all_genes_edger.csv'),
          quote = FALSE, 
          row.names = FALSE)

padj_cutoff = 0.05

sig_res_tbl_edger.bulk = dplyr::filter(res_tbl_edger.bulk, FDR < padj_cutoff) %>%
  dplyr::arrange(FDR)

write.csv(sig_res_tbl_edger.bulk,
          paste0(save_edger_dir,'/results/Overall_condition_pb_vs_ctrl_signif_genes_edger.csv'),
          quote = FALSE, 
          row.names = FALSE)
 

## generate enhanced volcano plot
EnhancedVolcano(res_tbl_edger.bulk,
                lab=res_tbl_edger.bulk$gene,
                title='Overall DE Genes',
                subtitle='Pb vs control mammary gland',
                x='logFC',
                y='FDR',
                xlim=c(-2,2),
                ylim=c(0,6),
                pCutoff=0.05,
                FCcutoff=0.5, 
                pointSize=4,
                labSize=8,
                colAlpha=0.8,
                cutoffLineType='twodash',
                cutoffLineWidth=0.5,
                legendLabSize=16,
                legendIconSize = 5,
                drawConnectors=TRUE)
ggsave(paste0(save_edger_dir,'/figures/overall_pb_vs_ctrl_DEgenes_evolcano.png'),dpi=320,width=15,height=10)



## apply edgeR DE anlaysis to all cell types across samples and conditions ----
get_edger_resultsAvsB <- function(clustx, padj_cutoff = 0.05) {
  
  cat('Beggining cluster ',clustx,'\n\n') # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls_limma[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  # create DGEList object
  dge = DGEList(counts=cluster_counts, samples=cluster_metadata)
  
  # remove label-sample combinations with fewer than 10 cells
  discarded = dge$samples$cell_count < 10
  dge = dge[,!discarded]
  cat('# of samples being discarded for ',clustx,'?\n')
  print(summary(discarded))
  
  # remove genes that are not expressed above a log-CPM threshold in a minimum number
  # of samples (determined from the size of the smallest treatment group in 
  # the experimental design, or pb group in this case.)
  # keep = filterByExpr(dge, group=dge$samples$condition,
  #                     min.count=10, min.total.count=15, min.prop=0.4)
  # dge = dge[keep,]
  # cat('# of genes being kept for ',clustx,'?\n')
  # print(summary(keep))
  
  keep = rowSums(dge$counts >= 10) >= smallestGroupSize # smallest group is 4
  dge = dge[keep, ]
  cat('# of genes being kept for ',clustx,'?\n')
  print(summary(keep))
  
  # correct for composition biases by computing normalization factors with the 
  # trimmed mean of M-values method (Robinson and Oshlack 2010), converting raw library
  # sizes in normalized effective library sizes
  dge = calcNormFactors(dge)
  
  # Plot mean-difference plot for each pseudobulk profile
  png(filename=paste0(save_edger_dir,'/figures/',clustx,'_mean-difference_plots.png'), units='in', width=10, height=10, res=320)
  par(mfrow=c(3,4))
  for (i in seq_len(ncol(dge))) {
    try(
      plotMD(dge, column=i),
      silent = TRUE
    )
  }
  try(
    mtext(paste0('Mean-difference plots for pseudobulk samples - ', clustx),
          outer=T,cex=1.15, side='top', adj=0, padj=1),
    silent=TRUE)
  dev.off()
  
  # Plot multi-dimensional scaling plots for each pseudobulk profile
  png(filename=paste0(save_edger_dir,'/figures/',clustx,'_multidimensional_scaling_plot.png'), units='in', width=10, height=10, res=320)
  plotMDS(edgeR::cpm(dge, log=TRUE), 
          col=ifelse(dge$samples$condition == 'pb', "red", "blue"))
  title(paste0('MDS plot for pseudobulk samples - ', clustx))
  dev.off()
  
  # set up design matrix
  design = model.matrix(~factor(condition), dge$samples)
  
  # estimate negative bibnomial dispersion
  dge = estimateDisp(dge, design)
  
  # plot biological coefficient of variation for each gene as a function of average abundance
  # The BCV is computed as the square root of the NB dispersion after 
  # empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are 
  # shown in blue and red, respectively.
  png(filename=paste0(save_edger_dir,'/figures/',clustx,'_BCV_plot.png'), units='in', width=10, height=10, res=320)
  plotBCV(dge)
  title(paste0('Biological coefficient of variation for ',clustx,' genes'))
  dev.off()
  
  # Fit a GLM to the counts for each gene and estimates the QL dispersion from the 
  # GLM deviance. We set robust=TRUE to avoid distortions from 
  # highly variable clusters (Phipson et al. 2016). The QL dispersion models the 
  # uncertainty and variability of the per-gene variance - which is not 
  # well handled by the NB dispersions, so the two dispersion types 
  # complement each other in the final analysis.
  fit = glmQLFit(dge, design, robust=TRUE)
  
  # Plot QL dispersion estimates for each gene as a function of abundance.
  # Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).
  png(filename=paste0(save_edger_dir,'/figures/',clustx,'_QLdispersion_estimates.png'), units='in', width=10, height=10, res=320)
  plotQLDisp(fit)
  title(paste0('QL dispersion estimates for ',clustx,' genes'))
  dev.off()
  
  # Test for differences in expression due to condition, Pb exposure
  result = glmQLFTest(fit, coef=ncol(design))
  cat('Differential expression results',)
  print(summary(decideTests(result)))
  
  res_tbl_edger = result %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    as.tibble()
  
  # save results
  write.csv(res_tbl_edger,
            paste0(save_edger_dir,'/results/',clustx,'_condition_pb_vs_ctrl_all_genes_edger.csv'),
            quote = FALSE, 
            row.names = FALSE)
  
  # subset the significant results
  sig_res_tbl_edger = dplyr::filter(res_tbl_edger, adj.P.Val < padj_cutoff) %>%
    dplyr::arrange(adj.P.Val)
  
  write.csv(sig_res_tbl_edger,
            paste0(save_edger_dir,'/results/',clustx,'_condition_pb_vs_ctrl_all_signif_genes_edger.csv'),
            quote = FALSE, 
            row.names = FALSE)
  
  # Generate results visualization plots
  
  ## Extract normalized counts from dds object
  # normalized_counts <- counts(dds, normalized = TRUE)
  # 
  # if ((nrow(sig_res) > 0) & (!is.null(nrow(sig_res)))) {
  #   ## Extract top 20 DEG from resLFC (make sure to order by padj)
  #   n = min(20, nrow(sig_res)) # top 20 or fewer DE genes
  #   top20_sig_genes <- sig_res %>%
  #     dplyr::arrange(padj) %>%
  #     dplyr::pull(gene) %>%
  #     head(n = n)
  #   
  #   print(head(top20_sig_genes))
  #   
  #   ## Convert wide matrix to long data frame for ggplot2
  #   
  #   ## Extract matching normalized count values from matrix
  #   top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  #   
  #   # Convert the matrix to a data frame and add the 'gene' column
  #   top20_sig_df <- data.frame(top20_sig_counts)
  #   top20_sig_df$gene <- rownames(top20_sig_counts)
  #   
  #   # Convert to data.table and melt into long format
  #   top20_sig_df <- melt(setDT(top20_sig_df), 
  #                        id.vars = c("gene"),
  #                        variable.name = "cluster_sample_id")
  #   
  #   ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  #   top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
  #   
  #   # Proceed with joining the metadata and plotting
  #   top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
  #                              by = "cluster_sample_id")
  #   
  #   ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
  #     geom_boxplot(outlier.shape=NA, width=0.2)+
  #     geom_jitter(height = 0, width = 0.5)+
  #     scale_y_continuous(trans = 'log10') +
  #     ylab("log10 of normalized expression level") +
  #     xlab("condition") +
  #     ggtitle(paste0("Top ", n, " Significant DE Genes")) +
  #     theme(plot.title = element_text(hjust = 0.5)) +
  #     facet_wrap(~ gene)
  #   ggsave(paste0("figures/", clustx, "_", contrast, "_top20_DE_genes.png"))
  #   
  #   ## Extract normalized counts for significant genes only
  #   sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
  #   
  #   ## generate pheatmap
  #   png(filename=paste0('figures/',clustx,'_',contrast,'_pheatmap_top17DEgene.png'), units='in', width=10, height=10, res=320)
  #   pheatmap(sig_counts, 
  #            color = heat_colors, 
  #            cluster_rows = TRUE, 
  #            show_rownames = FALSE,
  #            annotation_col = cluster_metadata[,c('condition','sample')],
  #            border_color = NA, 
  #            fontsize = 10, 
  #            scale = "row", 
  #            fontsize_row = 10, 
  #            height = 20)
  #   dev.off()
  #   
  # } else {
  #   warning(paste0("No significant genes found for ",clustx))
  #   cat(paste0("No significant genes found for ",clustx),'\n\n')
  # }
  # 
  
  ## generate enhanced volcano plot
  EnhancedVolcano(res_tbl_edger,
                  lab=res_tbl_voom$gene,
                  title=paste0(clustx,' DE Genes'),
                  subtitle='Pb vs control',
                  x='logFC',
                  y='adj.P.Val',
                  pCutoff=0.05,
                  FCcutoff=0.5, # roughly 50% change equates to 0.58 log2FC, or 1.5 FC
                  pointSize=4,
                  labSize=8,
                  colAlpha=0.8,
                  cutoffLineType='twodash',
                  cutoffLineWidth=0.5,
                  legendLabSize=16,
                  legendIconSize = 5,
                  drawConnectors=TRUE)
  ggsave(paste0(save_edger_dir,'/figures/',clustx,'_pb_vs_ctrl_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
  
}



# Run the limma script on all clusters comparing stimulated condition relative to control condition
map(names(counts_ls), get_edger_resultsAvsB, padj_cutoff = 0.05)







# summed = scuttle::aggregateAcrossCells(sce, 
#                                        id=colData(sce)[,c('celltype', 'sample')])
# 
# summed.filt = summed[,summed$ncells >= 10] # remove pseuduobulk samples with less than 10 cells
# 
# # keep = filterByExpr(summed.filt, group=summed.filt$condition,
# #                     min.count=10, min.total.count=15)
# # table(keep)
# # keep
# # FALSE  TRUE 
# # 10156  3093
# 
# # summed.filt = summed.filt[keep,] scran::pseudoBulkDGE does this filtering already
# 
# de.results = scran::pseudoBulkDGE(summed.filt, 
#                                   label=summed.filt$celltype,
#                                   design=~factor(condition),
#                                   coef='factor(condition)pb',
#                                   condition=summed.filt$condition,
#                                   sorted=TRUE,
#                                   method='edgeR',
#                                   robust=TRUE
# )
# 
# 
# 
# 
# # Create directories to save results if they don't already exist:
# if (!dir.exists("stats/edgeR/results")) { dir.create("stats/edgeR/results",
#                                                      recursive=TRUE) }
# if (!dir.exists("stats/edgeR/figures")) { dir.create("stats/edgeR/figures",
#                                                      recursive=TRUE) }
# 
# # set threshold
# padj_cutoff <- 0.05
# 
# ### loop through edgeR results to generate figures and save data for each cell type ----
# for (cluster in unique(sce$celltype)) {
#   
#   cat(paste0('Started generating tables and graphs for ',cluster, '\n'))
#   
#   curr.result = de.results[[cluster]]
#   res_tbl = curr.result[order(curr.result$FDR),] %>%
#     data.frame() %>%
#     rownames_to_column(var = 'gene') %>%
#     as.tibble()
#   
#   curr.metadata = metadata(curr.result)
#   
#   write.csv(res_tbl,
#             paste0('stats/edgeR/results/',cluster,"_condition_pb_vs_ctrl_all_genes_edger.csv"),
#             quote = FALSE, 
#             row.names = FALSE)
#   
#   ## Subset the significant results
#   sig_res = dplyr::filter(res_tbl, FDR < padj_cutoff) %>%
#     dplyr::arrange(FDR)
#   
#   write.csv(sig_res,
#             paste0('stats/edgeR/results/',cluster,"condition_pb_vs_ctrl_signif_genes_edger.csv"),
#             quote = FALSE, 
#             row.names = FALSE)
#   
#   png(filename=paste0('stats/edgeR/figures/',cluster,'_mean-difference_plots.png'), units='in', width=10, height=10, res=320)
#   par(mfrow=c(3,4))
#   for (i in seq_len(ncol(curr.metadata$y))) {
#     plotMD(curr.metadata$y, column=i)
#   }
#   mtext(paste0('Mean-difference plots for pseudobulk samples - ',cluster,
#                outer=T,cex=1.15, side='top', adj=0, padj=1))
#   dev.off()
#   
#   # Plot multi-dimensional scaling plots for each pseudobulk profile
#   png(filename=paste0('stats/edgeR/figures/',cluster,'_multidimensional_scaling_plot.png'), units='in', width=10, height=10, res=320)
#   plotMDS(edgeR::cpm(curr.metadata$y, log=TRUE), 
#           col=ifelse(curr.metadata$y$samples$condition == 'pb', "red", "blue"))
#   title(paste0('MDS plot for pseudobulk samples - ',cluster))
#   dev.off()
#   
#   # plot biological coefficient of variation for each gene as a function of average abundance
#   # The BCV is computed as the square root of the NB dispersion after 
#   # empirical Bayes shrinkage towards the trend. Trended and common BCV estimates are 
#   # shown in blue and red, respectively.
#   png(filename=paste0('stats/edgeR/figures/',cluster,'_BCV_plot.png'), units='in', width=10, height=10, res=320)
#   plotBCV(curr.metadata$y)
#   title(paste0('BCV for ',cluster, ' genes'))
#   dev.off()
#   
#   
#   # Plot QL dispersion estimates for each gene as a function of abundance.
#   # Raw estimates (black) are shrunk towards the trend (blue) to yield squeezed estimates (red).
#   png(filename=paste0('stats/edgeR/figures/',cluster,'_QLdispersion_estimates.png'), units='in', width=10, height=10, res=320)
#   plotQLDisp(curr.metadata$fit)
#   title(paste0('QL dispersion estimates for ',cluster, ' genes'))
#   dev.off()
#   
#   ## generate enhanced volcano plot
#   EnhancedVolcano(res_tbl,
#                   lab=res_tbl$gene,
#                   title=paste0(cluster,' DE Genes'),
#                   subtitle='Pb vs control',
#                   x='logFC',
#                   y='FDR',
#                   pCutoff=0.05,
#                   FCcutoff=0.5, # roughly 50% change equates to 0.58 log2FC, or 1.5 FC
#                   pointSize=4,
#                   labSize=8,
#                   colAlpha=0.8,
#                   cutoffLineType='twodash',
#                   cutoffLineWidth=0.5,
#                   legendLabSize=16,
#                   legendIconSize = 5,
#                   drawConnectors=TRUE)
#   ggsave(paste0('stats/edgeR/figures/',cluster,'_pb_vs_ctrl','_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
#   
#   cat(paste0('Finished with ',cluster, '\n'))
# }













# Limma Voom analysis ----
# Create directories to save results if they don't already exist:
if (!dir.exists("stats/limma/results")) { dir.create("stats/limma/results",
                                                     recursive=TRUE) }
if (!dir.exists("stats/limma/figures")) { dir.create("stats/limma/figures",
                                                     recursive=TRUE) }

## overall bulk analysis by condition with limma ----

# will use same filtered dge_bulk and design_bulk from edgeR analysis
# to fit limma voom glm.

# The voom pipeline from the limma package allows us to use sample weights 
# to better account for the variation in the precision of 
# each pseudo-bulk profile.
png(filename=paste0(save_limma_dir,'/figures/overall_mean-variance_trend_and_sample_weights.png'), units='in', width=10, height=10, res=320)
voom.bulk = voomWithQualityWeights(dge_bulk, design_bulk, plot=TRUE) # uses same dge_bulk object as edgeR
mtext('Overall', outer=T,cex=1.15, side='top', adj=0, padj=1)
dev.off()

# voom.bulk = voom(dge_bulk, design_bulk)
fit.voom.bulk = lmFit(voom.bulk)
# fit.voom.bulk = lmFit(dge_bulk, design_bulk)
# png(filename=paste0('stats/limma/figures/overall_mean-variance_trend_and_sample_weights.png'), units='in', width=10, height=10, res=320)
# fit.voom.bulk = voomLmFit(dge_bulk, design_bulk,
#                           sample.weights=TRUE, plot=TRUE,
#                           keep.EList=TRUE)
# mtext('Overall', outer=T,cex=1.15, side='top', adj=0, padj=1)
# dev.off()

fit.voom.bulk = eBayes(fit.voom.bulk, trend = T, robust=T)

# MD plots are same as those generated for edgeR above so won't
# plot them again here

result.voom.bulk = topTable(fit.voom.bulk, sort.by='p', n=Inf,
                            coef='factor(condition)pb')
head(result.voom.bulk, 10)

res_tbl_voom.bulk = result.voom.bulk %>%
  data.frame() %>%
  rownames_to_column(var = 'gene') %>%
  as.tibble()

write.csv(res_tbl_voom.bulk,
          paste0(save_limma_dir,'/results/overall_condition_pb_vs_ctrl_all_genes_limma.csv'),
          quote = FALSE, 
          row.names = FALSE)

sig_res_tbl_voom.bulk = dplyr::filter(res_tbl_voom.bulk, adj.P.Val < padj_cutoff) %>%
  dplyr::arrange(adj.P.Val)

write.csv(sig_res_tbl_voom.bulk,
          paste0(save_limma_dir,'/results/overall_condition_pb_vs_ctrl_all_signif_genes_limma.csv'),
          quote = FALSE, 
          row.names = FALSE)

# set threshold
padj_cutoff = 0.05


## Run limma voom by looping through each cluster----

get_limma_resultsAvsB <- function(clustx, padj_cutoff = 0.05) {
  
  cat('Beggining cluster ',clustx,'\n\n') # useful for debugging
  
  # Extract counts matrix and metadata for cluster x
  idx <- which(names(counts_ls) == clustx)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls_limma[[idx]]
  
  # Print error message if sample names do not match
  if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
    print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
  }
  
  # create DGEList object
  dge = DGEList(counts=cluster_counts, samples=cluster_metadata)

  # remove label-sample combinations with fewer than 10 cells
  discarded = dge$samples$cell_count < 10
  dge = dge[,!discarded]
  cat('TRUE/FALSE for samples being discarded for ',clustx,'?\n')
  print(summary(discarded))

  # remove genes that are not expressed above a log-CPM threshold in a minimum number
  # of samples (determined from the size of the smallest treatment group in 
  # the experimental design, or pb group in this case.)
  keep = filterByExpr(dge, group=dge$samples$condition,
                      min.count=10, min.total.count=15)
  dge = dge[keep,]
  cat('TRUE/FALSE for gene being kept for ',clustx,'?\n')
  print(summary(keep))
  
  
  # correct for composition biases by computing normalization factors with the 
  # trimmed mean of M-values method (Robinson and Oshlack 2010), convertering raw library
  # sizes in normalized effective library sizes
  dge = calcNormFactors(dge)

  # Plot mean-difference plot for each pseudobulk profile
  png(filename=paste0(save_limma_dir,'/figures/',unique(dge$samples$cluster_id),'_mean-difference_plots.png'), units='in', width=10, height=10, res=320)
  par(mfrow=c(3,4))
  for (i in seq_len(ncol(dge))) {
    try(
    plotMD(dge, column=i),
    silent = TRUE
    )
  }
  try(
    mtext(paste0('Mean-difference plots for pseudobulk samples - ', unique(dge$samples$cluster_id)),
        outer=T,cex=1.15, side='top', adj=0, padj=1),
    silent=TRUE)
  dev.off()
  
  # Plot multi-dimensional scaling plots for each pseudobulk profile
  png(filename=paste0(save_limma_dir,'/figures/',unique(dge$samples$cluster_id),'_multidimensional_scaling_plot.png'), units='in', width=10, height=10, res=320)
  plotMDS(edgeR::cpm(dge, log=TRUE), 
          col=ifelse(dge$samples$condition == 'pb', "red", "blue"))
  title(paste0('MDS plot for pseudobulk samples - ', unique(dge$samples$cluster_id)))
  dev.off()
  
  # set up design matrix
  design = model.matrix(~factor(condition), dge$samples)
  
  
  # fit limma voom model with trend and sample weights
  png(filename=paste0(save_limma_dir,'/figures/',clustx,'mean-variance_trend_and_sample_weights.png'), units='in', width=10, height=10, res=320)
  lv = voomWithQualityWeights(dge, design, plot=TRUE) # Combine voom observational-level precision weights with sample-specific quality weights 
  mtext(clustx, outer=T,cex=1.15, side='top', adj=0, padj=1)
  dev.off()
  
  # voom.bulk = voom(dge_bulk, design_bulk)
  fit.voom = lmFit(lv)
  # fit.voom.bulk = lmFit(dge, design)
  # png(filename=paste0('stats/limma/figures/',clustx,'_mean-variance_trend_and_sample_weights.png'), units='in', width=10, height=10, res=320)
  # fit.voom = voomLmFit(dge, design,
  #                           sample.weights=TRUE, plot=TRUE,
  #                           keep.EList=TRUE)
  # mtext('Overall', outer=T,cex=1.15, side='top', adj=0, padj=1)
  # dev.off()
  
  fit.voom = eBayes(fit.voom, trend = T, robust=T)
  
  result.voom = topTable(fit.voom, sort.by='p', n=Inf,
                              coef='factor(condition)pb')
  cat('\nTop 10 DE genes for ',clustx,'\n')
  print(head(result.voom, 10))
  
  res_tbl_voom = result.voom %>%
    data.frame() %>%
    rownames_to_column(var = 'gene') %>%
    as.tibble()
  
  # save results
  write.csv(res_tbl_voom,
            paste0(save_limma_dir,'/results/',clustx,'_condition_pb_vs_ctrl_all_genes_limma.csv'),
            quote = FALSE, 
            row.names = FALSE)
  
  # subset the significant results
  sig_res_tbl_voom = dplyr::filter(res_tbl_voom, adj.P.Val < padj_cutoff) %>%
    dplyr::arrange(adj.P.Val)
  
  write.csv(sig_res_tbl_voom,
            paste0(save_limma_dir,'/results/',clustx,'_condition_pb_vs_ctrl_all_signif_genes_limma.csv'),
            quote = FALSE, 
            row.names = FALSE)

  # Generate results visualization plots
  
  ## Extract normalized counts from dds object
  # normalized_counts <- counts(dds, normalized = TRUE)
  # 
  # if ((nrow(sig_res) > 0) & (!is.null(nrow(sig_res)))) {
  #   ## Extract top 20 DEG from resLFC (make sure to order by padj)
  #   n = min(20, nrow(sig_res)) # top 20 or fewer DE genes
  #   top20_sig_genes <- sig_res %>%
  #     dplyr::arrange(padj) %>%
  #     dplyr::pull(gene) %>%
  #     head(n = n)
  #   
  #   print(head(top20_sig_genes))
  #   
  #   ## Convert wide matrix to long data frame for ggplot2
  #   
  #   ## Extract matching normalized count values from matrix
  #   top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
  #   
  #   # Convert the matrix to a data frame and add the 'gene' column
  #   top20_sig_df <- data.frame(top20_sig_counts)
  #   top20_sig_df$gene <- rownames(top20_sig_counts)
  #   
  #   # Convert to data.table and melt into long format
  #   top20_sig_df <- melt(setDT(top20_sig_df), 
  #                        id.vars = c("gene"),
  #                        variable.name = "cluster_sample_id")
  #   
  #   ## Replace "." by " " in cluster_sample_id variable (melt() introduced the ".")
  #   top20_sig_df$cluster_sample_id <- gsub("\\.", " ", top20_sig_df$cluster_sample_id)
  #   
  #   # Proceed with joining the metadata and plotting
  #   top20_sig_df <- plyr::join(top20_sig_df, as.data.frame(colData(dds)),
  #                              by = "cluster_sample_id")
  #   
  #   ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
  #     geom_boxplot(outlier.shape=NA, width=0.2)+
  #     geom_jitter(height = 0, width = 0.5)+
  #     scale_y_continuous(trans = 'log10') +
  #     ylab("log10 of normalized expression level") +
  #     xlab("condition") +
  #     ggtitle(paste0("Top ", n, " Significant DE Genes")) +
  #     theme(plot.title = element_text(hjust = 0.5)) +
  #     facet_wrap(~ gene)
  #   ggsave(paste0("figures/", clustx, "_", contrast, "_top20_DE_genes.png"))
  #   
  #   ## Extract normalized counts for significant genes only
  #   sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
  #   
  #   ## generate pheatmap
  #   png(filename=paste0('figures/',clustx,'_',contrast,'_pheatmap_top17DEgene.png'), units='in', width=10, height=10, res=320)
  #   pheatmap(sig_counts, 
  #            color = heat_colors, 
  #            cluster_rows = TRUE, 
  #            show_rownames = FALSE,
  #            annotation_col = cluster_metadata[,c('condition','sample')],
  #            border_color = NA, 
  #            fontsize = 10, 
  #            scale = "row", 
  #            fontsize_row = 10, 
  #            height = 20)
  #   dev.off()
  #   
  # } else {
  #   warning(paste0("No significant genes found for ",clustx))
  #   cat(paste0("No significant genes found for ",clustx),'\n\n')
  # }
  # 
  
  ## generate enhanced volcano plot
  EnhancedVolcano(res_tbl_voom,
                  lab=res_tbl_voom$gene,
                  title=paste0(clustx,' DE Genes'),
                  subtitle='Pb vs control',
                  x='logFC',
                  y='adj.P.Val',
                  pCutoff=0.05,
                  FCcutoff=0.5, # roughly 50% change equates to 0.58 log2FC, or 1.5 FC
                  pointSize=4,
                  labSize=8,
                  colAlpha=0.8,
                  cutoffLineType='twodash',
                  cutoffLineWidth=0.5,
                  legendLabSize=16,
                  legendIconSize = 5,
                  drawConnectors=TRUE)
  ggsave(paste0(save_limma_dir,'/figures/',clustx,'_pb_vs_ctrl_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
  
}



# Run the limma script on all clusters comparing stimulated condition relative to control condition
map(names(counts_ls), get_limma_resultsAvsB, padj_cutoff = 0.05)




# Compare DE results across DESeeq2, edgeR, and Limma ----

# Define paths for the directories containing the result files
deseq2rr_dir <- '3_weeks/stats_rerun/DESeq2/results/'
deseq2zinb_dir <- '3_weeks/stats/DESeq2/results/'
deseq2_dir <- '3_weeks/stats/DESeq2_unfiltered/results/'
edger_dir <- '3_weeks/stats/edgeR/results/'
limma_dir <-'3_weeks/stats/limma/results/'

# Function to read gene lists for a cell type from different methods,
# will return gene lists only be default but can return logFC too if indicated.
read_gene_list <- function(cell_type, return_logFC=FALSE) {
  if (return_logFC) {
    deseq2rr_file <- paste0(deseq2rr_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    deseq2zinb_file <- paste0(deseq2zinb_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    deseq2_file <- paste0(deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    edger_file <- paste0(edger_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_edger.csv")
    limma_file <- paste0(limma_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_limma.csv")
    
    # Read the full data for log fold changes
    deseq2rr_data <- read.csv(deseq2rr_file)[, c("gene", "log2FoldChange")]
    deseq2zinb_data <- read.csv(deseq2zinb_file)[, c("gene", "log2FoldChange")]
    deseq2_data <- read.csv(deseq2_file)[, c("gene", "log2FoldChange")]
    edger_data <- read.csv(edger_file)[, c("gene", "logFC")]
    limma_data <- read.csv(limma_file)[, c("gene", "logFC")]
    
    # Rename columns for clarity
    colnames(deseq2rr_data)[2] = "logFC_DESeq2rr_havard_low-count-filter"
    colnames(deseq2zinb_data)[2] = "logFC_DESeq2zinb_deseq2_vignette"
    colnames(deseq2_data)[2] = "logFC_DESeq2_harvard_vignette"
    colnames(edger_data)[2] = "logFC_edgeR"
    colnames(limma_data)[2] = "logFC_limma"
    
    # Merge all data frames by gene name to build dataframe of genes x logFC_method
    merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(deseq2rr_data, deseq2zinb_data, deseq2_data, edger_data, limma_data))
    
    return(merged_data)
  
  } else {
    deseq2rr_file <- paste0(deseq2rr_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    deseq2zinb_file <- paste0(deseq2zinb_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    deseq2_file <- paste0(deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    edger_file <- paste0(edger_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_edger.csv")
    limma_file <- paste0(limma_dir, cell_type, "_condition_pb_vs_ctrl_all_signif_genes_limma.csv")
    
    # Read the gene column
    deseq2rr_genes = read.csv(deseq2rr_file)$gene
    deseq2zinb_genes = read.csv(deseq2zinb_file)$gene
    deseq2_genes = read.csv(deseq2_file)$gene
    edger_genes = read.csv(edger_file)$gene
    limma_genes = read.csv(limma_file)$gene
    
    
    return(list(deseq2rr = deseq2rr_genes, deseq2zinb = deseq2zinb_genes, deseq2 = deseq2_genes, edger = edger_genes, limma = limma_genes))
  }
}

# Function to create a Venn diagram for a cell type
colors <- brewer.pal(3, "Pastel1")

create_venn_diagram <- function(gene_lists, cell_type) {

  # get intersecting gene names and format them
  intersections = calculate.overlap(gene_lists)
  formatted_intersections <- lapply(intersections, function(x) {
    count = length(x)
    if (count > 0) {
      paste0(count, "\n", paste(x, collapse = "\n"))
    } else {
      "0"
    }
  })
  
  
  vennplot <- venn.diagram(gene_lists,
                           filename = paste0('3_weeks/stats/venn_diagrams/',cell_type,'_venn_diagram.png'),
                           output=TRUE,
                           category.names = c("DESeq2_filtered", 'DESeq2_zinb','DESeq2',"edgeR", "limma"),
                           fill = colors,
                           lwd = 2,
                           lty = 'blank',
                           # cex = .4,
                           fontface = "bold",
                           fontfamily = "sans",
                           # cat.cex = 0.5,
                           cat.fontface = "bold",
                           cat.default.pos = "outer",
                           cat.pos = c(-27, 27, 135),
                           cat.dist = c(0.055, 0.055, 0.085),
                           cat.fontfamily = "sans",
                           rotation = 1,
                           label.col='black',
                           cex=0.3,
                           cat.cex=0.7)
  
  # Update labels with gene counts and names for each section
  for (i in 1:length(vennplot)) {
    if (!is.null(vennplot[[i]]$label)) {
      vennplot[[i]]$label <- formatted_intersections[[i]]
    }
  }
  
  return(vennplot)
}


### grab celltypes from results ----
### directory of differential expression results 
de_directory = '3_weeks/stats_rerun/DESeq2/results/'

### grab de result file names and extract cell types 
de_results = list.files(path=de_directory, pattern='all_genes')

### List of cell types (replace with your actual cell types)
celltypes = sub(pattern='_condition.*',replacement='',x=de_results)


## Venn diagrams of overlapping differential expressed genes (FDR <0.05) across methods ----
if (!dir.exists("3_weeks/stats/venn_diagrams/")) { dir.create("3_Weeks/stats/venn_diagrams/",
                                                          recursive=TRUE)}
# Loop over each cell type and create the Venn diagrams
for (celltype in celltypes) {
  gene_lists = read_gene_list(celltype)
  create_venn_diagram(gene_lists, celltype)
}

### try other ways. Want to include gene names in venn diagrams ----
library(ggVennDiagram)

# Loop over each cell type and create the Venn diagrams with ggVennDiagram
for (celltype in celltypes) {
  gene_lists <- read_gene_list(celltype)
  
  # Create the Venn diagram
  p <- ggVennDiagram(gene_lists, label = "both") + 
    scale_fill_gradient(low = "lightblue", high = "purple") +
    theme(legend.position = "none") +
    labs(title = paste("Venn Diagram for", cell_type))
  
  # Save the plot to file
  ggsave(paste0('3_weeks/stats/venn_diagrams/', cell_type, '_ggvenn_diagram.png'), plot = p)
}



# try another way
cell_type = 'Endothelial'
gene_lists <- read_gene_list(cell_type)
intersections = calculate.overlap(gene_lists[1:3])
formatted_intersections <- lapply(intersections, function(x) {
  count = length(x)
  if (count > 0) {
    paste0(count, "\n", paste(x, collapse = "\n"))
  } else {
    "0"
  }
})
vennplot <- venn.diagram(gene_lists[1:3],
                         # filename=NULL,
                         filename = paste0('3_weeks/stats/venn_diagrams/deseq2s/',cell_type,'_venn_diagram.png'),
                         output=TRUE,
                         category.names = c("DESeq2_filtered", 'DESeq2_zinb','DESeq2'),
                         fill = colors,
                         lwd = 2,
                         lty = 'blank',
                         # cex = .4,
                         fontface = "bold",
                         fontfamily = "sans",
                         # cat.cex = 0.5,
                         cat.fontface = "bold",
                         cat.default.pos = "outer",
                         cat.pos = c(-27, 27, 135),
                         cat.dist = c(0.055, 0.055, 0.085),
                         cat.fontfamily = "sans",
                         rotation = 1,
                         label.col='black',
                         cex=0.7,
                         cat.cex=0.7)

# run if set filename above to NULL
# grid.newpage()
# grid.draw(vennplot)

# Update labels with gene counts and names for each section
for (i in 1:length(vennplot)) {
  if (!is.null(vennplot[[i]]$label)) {
    vennplot[[i]]$label <- formatted_intersections[[i]]
  }
}



vennplot <- venn.diagram(gene_lists,
                         filename=NULL,
                         category.names = c("DESeq2", "edgeR", "limma"),
                         alpha=c(0.5,0.5,0.5),
                         fill = colors,
                         lwd = 2,
                         lty = 'blank',
                         cex = 1.5,
                         cat.cex = 1.5,
                         fontface = "bold",
                         fontfamily = "sans",
                         cat.fontface = "bold",
                         cat.default.pos = "outer",
                         cat.pos = c(-27, 27, 135),
                         cat.dist = c(0.055, 0.055, 0.085),
                         cat.fontfamily = "sans",
                         rotation = 1,
                         label.col='black',
                         print.mode='raw')
grid.newpage()
grid.draw(vennplot)

lapply(vennplot, function(i) i$label)


# Calculate overlap site
overlaps <- calculate.overlap(gene_lists)
overlaps <- overlaps[str_sort(names(overlaps), numeric = TRUE)] # sort base on numeric value


# Apply name to global variable
# Index of venn diagram start at calculate.overlaps + 8. You have to find the index value by yourself for (3,5,6,7,.. venn)
walk2(seq(overlaps) + 8, seq(overlaps),  
      function(x,y) {vennplot[[x]]$label <<- paste0(overlaps[[y]], collapse = "\n")})

# Draw plot
grid.draw(v)
grid.newpage()
grid.draw(vennplot)


library(ggVennDiagram)

ggVennDiagram(gene_lists,
              show_intersect = TRUE,
              set_color='black') + 
  theme(legend.position = "none") +
  labs(title = "Venn Diagram for Adipocytes")



## Check correlation of log-fold change between all three methods ----
celltype = 'Macrophage.ma'
merged_de_res = read_gene_list(celltype,return_logFC = TRUE)

# Calculate correlation matrix
# Use pairwise complete observations to handle missing values
correlation_matrix = cor(merged_de_res[, -1], use = "pairwise.complete.obs")

# Print correlation matrix
print(correlation_matrix)

# Optionally, visualize the correlation matrix
library(corrplot)
corrplot(correlation_matrix, method = "square",col=viridis(option='viridis',n=256))

pheatmap(correlation_matrix, color=viridis(option='viridis',n=256))


# sessionInfo() ----
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] EnhancedVolcano_1.20.0      ggrepel_0.9.6               ggpubr_0.6.0                data.table_1.16.0           RColorBrewer_1.1-3          png_0.1-8                  
# [7] apeglm_1.24.0               pheatmap_1.0.12             SeuratObject_5.0.2          sp_2.1-4                    SingleCellExperiment_1.24.0 reshape2_1.4.4             
# [13] DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0       matrixStats_1.4.1           GenomicRanges_1.54.1       
# [19] GenomeInfoDb_1.38.8         IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1         edgeR_4.0.16                limma_3.58.1               
# [25] Matrix.utils_0.9.8          Matrix_1.6-5                cowplot_1.1.3               lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
# [31] dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1              
# [37] tidyverse_2.0.0            
# 
# loaded via a namespace (and not attached):
#   [1] later_1.3.2                   bitops_1.0-8                  filelock_1.0.3                graph_1.80.0                  XML_3.99-0.17                 lifecycle_1.0.4              
# [7] rstatix_0.7.2                 globals_0.16.3                lattice_0.22-6                MASS_7.3-60.0.1               backports_1.5.0               magrittr_2.0.3               
# [13] yaml_2.3.10                   metapod_1.10.1                httpuv_1.6.15                 spam_2.10-0                   DBI_1.2.3                     abind_1.4-8                  
# [19] zlibbioc_1.48.2               RCurl_1.98-1.16               rappdirs_0.3.3                GenomeInfoDbData_1.2.11       irlba_2.3.5.1                 listenv_0.9.1                
# [25] GSVA_1.50.5                   annotate_1.80.0               dqrng_0.4.1                   parallelly_1.38.0             DelayedMatrixStats_1.24.0     codetools_0.2-20             
# [31] DelayedArray_0.28.0           scuttle_1.12.0                tidyselect_1.2.1              ScaledMatrix_1.10.0           BiocFileCache_2.10.2          BiocNeighbors_1.20.2         
# [37] progressr_0.14.0              Formula_1.2-5                 bbmle_1.0.25.1                tools_4.3.1                   Rcpp_1.0.13                   glue_1.7.0                   
# [43] SparseArray_1.2.4             HDF5Array_1.30.1              withr_3.0.1                   numDeriv_2016.8-1.1           BiocManager_1.30.25           fastmap_1.2.0                
# [49] rhdf5filters_1.14.1           bluster_1.12.0                fansi_1.0.6                   digest_0.6.37                 rsvd_1.0.5                    timechange_0.3.0             
# [55] R6_2.5.1                      mime_0.12                     colorspace_2.1-1              RSQLite_2.3.7                 utf8_1.2.4                    generics_0.1.3               
# [61] httr_1.4.7                    S4Arrays_1.2.1                pkgconfig_2.0.3               gtable_0.3.5                  blob_1.2.4                    XVector_0.42.0               
# [67] htmltools_0.5.8.1             carData_3.0-5                 dotCall64_1.1-1               GSEABase_1.64.0               scales_1.3.0                  scran_1.30.2                 
# [73] rstudioapi_0.16.0             tzdb_0.4.0                    coda_0.19-4.1                 curl_5.2.3                    bdsmatrix_1.3-7               cachem_1.1.0                 
# [79] rhdf5_2.46.1                  BiocVersion_3.18.1            parallel_4.3.1                AnnotationDbi_1.64.1          pillar_1.9.0                  grid_4.3.1                   
# [85] vctrs_0.6.5                   promises_1.3.0                BiocSingular_1.18.0           car_3.1-3                     dbplyr_2.5.0                  beachmat_2.18.1              
# [91] xtable_1.8-4                  cluster_2.1.6                 mvtnorm_1.3-1                 cli_3.6.3                     locfit_1.5-9.10               compiler_4.3.1               
# [97] rlang_1.1.4                   crayon_1.5.3                  grr_0.9.5                     future.apply_1.11.2           ggsignif_0.6.4                emdbook_1.3.13               
# [103] plyr_1.8.9                    stringi_1.8.4                 BiocParallel_1.36.0           munsell_0.5.1                 Biostrings_2.70.3             hms_1.1.3                    
# [109] sparseMatrixStats_1.14.0      bit64_4.5.2                   future_1.34.0                 Rhdf5lib_1.24.2               KEGGREST_1.42.0               statmod_1.5.0                
# [115] shiny_1.9.1                   interactiveDisplayBase_1.40.0 AnnotationHub_3.10.1          igraph_2.0.3                  broom_1.0.7                   memoise_2.0.1                
# [121] bit_4.5.0


