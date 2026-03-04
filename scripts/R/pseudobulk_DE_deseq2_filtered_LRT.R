# to install Matrix.utils
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)


# Load libraries
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(EnhancedVolcano)
library(qs2)

# set seed for reproducibility
set.seed(2024)

## Set a color-blind friendly palette
heat_colors <- rev(brewer.pal(11, "PuOr"))

# read in data ----
# 
# # directories for saving. Make sure they do not end in /
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/figures/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/data/"
# save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2/"
# # save_deseq2_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2_gene_fltrd/"
# 
# # Create directories to save results if they don't already exist:
# dirs = c(output_data_dir, output_fig_dir, save_deseq2_dir)
# 
# for (dir in dirs) {
#   if (!dir.exists(dir)) { dir.create(dir,
#                                      recursive=TRUE) }
# }


nthreads = 2

analysis <- 'cellbender_analysis'
## directories for integrated data ----
fig_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")
stats_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/")

# Create directories to save results if they don't already exist:
dirs = c(fig_dir_integrated, data_dir_integrated, annot_dir_integrated, stats_dir_integrated)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# timepoints. 3 weeks, 10 weeks, 7 months, 18 months
timepoints = c('wk3','wk10','mn7','mn18')



# seurat = readRDS(paste0(data_dir_integrated,'de_seurat_annotated.rds'))

seurat = qs_read(paste0(data_dir_integrated,'epithelial.qs'), nthreads = 2)

DefaultAssay(seurat) = 'RNA'
seurat = DietSeurat(seurat,assays = 'RNA')



# criteria for filtering later
n_pb = length(grep('pb',unique(seurat$sample)))
n_ctrl = length(grep('ctrl',unique(seurat$sample)))
smallestGroupSize = ifelse(n_pb < n_ctrl, n_pb, n_ctrl)
cat('Size of smaller group is',smallestGroupSize,'\n There are',n_pb,
    'Pb samples and',n_ctrl,'controls.')
# Size of smaller group is 4 
# There are 4 Pb samples and 6 controls.

# How many cells per condition per cell type? ----
# Looks like not enough muscle or erythroid cells for analysis
table(seurat$celltype, seurat$condition)
#                         ctrl    pb
# B cells                 8160 10218
# CD4+ T cells            5086  5816
# CD8+ T cells            4080  4753
# Tregs                    887  1247
# Endothelial             1255  1082
# Fibroblast               976   732
# Proliferating cells      618   935
# Dendritic cells          618   864
# Adipocytes              1044   865
# Activated CD4+ T Cells   558   583
# Macrophage.Ma            583   525
# T helper 17 cells        528   524
# Myeloid.Mb               886  1021
# Reticular cells          230   368
# Basal-Myoepithelial      254   229
# Schwann cells            204   131
# Luminal.HS                91    79
# Muscle                    76     0
# Erythroid cells           12    31

# How many cells will each pseudobulk sample have? ----
# looks like there won't be enough cells to analyze muscle or erythroid cells
table(seurat$celltype, seurat$sample)
#                         ctrl21 ctrl22 ctrl23 ctrl24 ctrl25 ctrl26 pb27 pb28 pb29 pb30 pb31 pb32
# B cells                   907   1410   1461    934   1443   2005 2343 1312 1722 2292 1502 1047
# CD4+ T cells              655    814    982    563    872    14712007  691 1084 1208  697  659
# CD8+ T cells              634    690    762    419    706    869 1034  520  977 1007  623  592
# Tregs                     118    135    135     90    161    248  310  114  273  254  159  137
# Endothelial               181     82    301    110    283    298  231  311  123  207  102  108
# Fibroblast                 57     93    424    114    157    131  123  221   87  120   89   92
# Proliferating cells        86     86     92     52    116    186  223  123  175  182  127  105
# Dendritic cells            76     98     61     56    142    185  176   89  162  213  120  104
# Adipocytes                171    102    258     79    226    208  145  231  112  193   74  110
# Activated CD4+ T Cells     60     83    110     66     87    152  110   71  123  125   85   69
# Macrophage.Ma              61     50    212     66    117     77   71  180   59  117   39   59
# T helper 17 cells          53     65    122     52     87    149  102   60  112  116   80   54
# Myeloid.Mb                 84    142    177     95    206    182  214  143  166  245  124  129
# Reticular cells            12     44     36     24     58     56   75   39   84   70   52   48
# Basal-Myoepithelial        45     16     48     26     57     62   42   73   25   34   23   32
# Schwann cells              30      5     50     15     63     41   25   39   10   23    7   27
# Luminal.HS                 44      6     19      4     13      5    6   14   17   17    8   17
# Muscle                      1      1     71      0      3      0    0    0    0    0    0    0
# Erythroid cells             1      1      4      2      2      2   24    1    2    3    0    1

# remove muscle cells from analysis since it's only made up of ctrl samples. No condition to compare to
seurat = subset(seurat, idents=c('Muscle','Erythroid cells'), invert=TRUE) # removed 67 muscle cells from data
seurat$celltype = Idents(seurat) # reassign with muscle and erythroid cluster removed

# Extract raw counts and metadata to create SingleCellExperiment object
counts <- JoinLayers(seurat)@assays$RNA$counts

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
# [1] 13568 56030 cellbender processed data

counts(sce)[1:6, 1:6]

# Extract unique names of clusters (= levels of cluster_id factor variable)
cluster_names <- levels(colData(sce)$celltype)
cluster_names
# [1] "B cells"                "CD4+ T cells"           "CD8+ T cells"           "Tregs"                  "Endothelial"           
# [6] "Fibroblast"             "Proliferating cells"    "Dendritic cells"        "Adipocytes"             "Activated CD4+ T Cells"
# [11] "Macrophage.Ma"          "T helper 17 cells"      "Myeloid.Mb"             "Reticular cells"        "Basal-Myoepithelial"   
# [16] "Schwann cells"          "Luminal.HS"

# Total number of clusters
length(cluster_names)
# [1] 17


# Extract unique names of samples (= levels of sample_id factor variable)
sample_names <- unique(colData(sce)$sample)
sample_names
# [1] "ctrl21" "ctrl22" "ctrl23" "ctrl24" "ctrl25" "ctrl26" "pb27"   "pb28"   "pb29"   "pb30"   "pb31"   "pb32"

# Total number of samples
length(sample_names)
# [1] 12

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("celltype", "sample")]
groups$sample = as.factor(groups$sample) # coerce sample to factor data type

dim(groups)
# [1] 56030     2

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
# [1]   204 13568

aggr_counts[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
aggr_counts[1:6, 1:6]

dim(aggr_counts)
# [1] 13249   179
# [1] 13568   204 cellbender processed data

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
# [1] 56030     2 cellbender processed data
head(metadata)

rownames(metadata) # barcodes

# Exclude duplicated rows, reducing metadata to 10 rows, one for each sample
metadata <- metadata[!duplicated(metadata), ]

dim(metadata)
# [1] 10  2

# [1] 12  2 cellbender processed data

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
# [1] 56030 cellbender processed data

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
# [1]    12 13568 cellbender processed data

# How many reads per pseudobulk sample?
rowSums(aggr_counts_bulk)
# ctrl21  ctrl22  ctrl23  ctrl24  ctrl25  ctrl26    pb27    pb28    pb30    pb31 
# 1882676 2206272 3705188 1801452 3046964 4670048 4882792 2363042 4383457 1984282 

# ctrl21  ctrl22  ctrl23  ctrl24  ctrl25  ctrl26    pb27    pb28    pb29    pb30    pb31    pb32  cellbender processed data
# 2150062 2392421 3859676 2055833 2422721 4424842 4143984 2164626 4271990 3820166 2650831 1725627

# How many reads per pseudobulk sample?

aggr_counts_bulk[1:6, 1:6]

# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts_bulk <- t(aggr_counts_bulk)
aggr_counts_bulk[1:6, 1:6]

dim(aggr_counts_bulk)
# [1] 13249    10
# [1] 13568    12



# create copy of metadata and add cell count col
metadata_bulk = copy(metadata)
metadata_bulk$cell_count = as.numeric(table(seurat$sample))

# Check matching of matrix columns and metadata rows
all(colnames(aggr_counts_bulk) == rownames(metadata_bulk))

# zero-count proportion?
(sum(rowSums(aggr_counts_bulk == 0))) / (nrow(aggr_counts_bulk) * ncol(aggr_counts_bulk))
# [1] 0.01206129
# [1] 0.01576626

# Create DESeq2 object        
dds_bulk <- DESeqDataSetFromMatrix(aggr_counts_bulk, 
                                   colData = metadata_bulk, 
                                   design = ~ condition)

# prefiltering low count genes.
# Want genes that have at least 10 counts in at least 4 samples. use 6 for cellbender processed data
# following recommendation from
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
# size of smaller group. there are 4 Pb samples and 6 controls
keep = rowSums(counts(dds_bulk) >= 10) >= smallestGroupSize
table(keep)
# keep
# FALSE  TRUE
# 2071 11497  discarding 2071 genes
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
scr_bulk = scran::computeSumFactors(dds_bulk) # scran size factors

sizeFactors(dds_bulk) = sizeFactors(scr_bulk)

#### Run DESeq2 differential expression analysis ----
dds_bulk <- DESeq(dds_bulk,
                  test='LRT',
                  fitType='parametric',
                  reduced= ~1,
                  minReplicatesForReplace=Inf,
                  useT=TRUE,
                  minmu=1e-6)

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
          paste0(save_deseq2_dir,"results/overall_condition_", 
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
ggsave(paste0(save_deseq2_dir,'figures/overall_DEgenes_evolcano.png'),dpi=320,width=15,height=10)





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
  
  # zero_prop = (sum(rowSums(cluster_counts == 0))) / (nrow(cluster_counts) * ncol(cluster_counts))
  # cat('Proportion of zero counts in',clustx,'for 13249 genes across 10 samples:\n')
  # print(zero_prop)
  # 
  # zero_props[[clustx]] = zero_prop
  
  
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ condition)
  
  # prefiltering low count genes. 
  # Want genes that have at least 10 counts in at least 4 samples. 6 for cellbender data since we haven't dropped any samples
  # following recommendation from
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  # smallestGroupSize = 6 # size of smaller group. there are 4 Pb samples in cellranger data. 6 per group in cellbender data
  # keep = rowSums(counts(dds) >= 10) >= smallestGroupSize
  # dds = dds[keep,]
  
  
  # zero_prop_post_filter = (sum(rowSums(dds@assays@data$counts == 0))) / 
  #   (nrow(dds@assays@data$counts) * ncol(dds@assays@data$counts))
  # 
  # cat('Proportion of zero counts in',clustx,'for 13249 genes across 10 samples
  #     after low-count gene filtering:\n')
  # print(zero_prop_post_filter)
  # 
  # zero_props_post_filtering[[clustx]] = zero_prop_post_filter
  # 
  # change_zero_prop = zero_prop - zero_prop_post_filter
  # cat('Change in zero count proportions after filtering for',clustx,':')
  # print(change_zero_prop)
  # 
  # change_in_zero_props[[clustx]] = change_zero_prop
  
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
    top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ,drop=FALSE]
    
    # Convert the matrix to a data frame and add the 'gene' column
    top20_sig_df <- data.frame(top20_sig_counts)
    top20_sig_df$gene <- rownames(top20_sig_counts)
    
    # Convert to data.table and melt into long format
    top20_sig_df <- melt(setDT(top20_sig_df), 
                         id.vars = c("gene"),
                         variable.name = "cluster_sample_id",
                         drop=FALSE)
    
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
    if ((nrow(sig_res) >= 2) & (!is.null(nrow(sig_res)))) {
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
      message("Skipping pheatmap: fewer than 2 genes found (nrow < 2).")
    }
    
    
    
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




# Compare DE results across DESeq2 variations----

# Define paths for the directories containing the result files
cr_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/'
deseq2rr_dir <- paste0(cr_deseq2_dir,'stats_rerun/DESeq2/results/') # filtered wald
deseq2LRT_dir <- paste0(cr_deseq2_dir,'stats/DESeq2/results/') # filtered LRT
deseq2_dir <- paste0(cr_deseq2_dir,'stats/DESeq2_unfiltered/results/') # unfiltered Wald
deseq2uLRT_dir <- paste0(cr_deseq2_dir,'stats/DESeq2_unfiltered_LRT/results/') # unfiltered LRT
cb_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2/results/' # cellbender unfiltered LRT
cb_deseq2_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/results/'
crv9_deseq2_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2_gene_fltrd/results/'
crv9_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2/results/'
# Function to read gene lists for a cell type from different methods,
# will return gene lists only be default but can return logFC too if indicated.
# read_gene_list <- function(cell_type, return_logFC=FALSE) {
#   if (return_logFC) {
#     deseq2rr_file <- paste0(deseq2rr_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#     deseq2_file <- paste0(deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
# 
#     # Read the full data for log fold changes
#     deseq2rr_data <- read.csv(deseq2rr_file)[, c("gene", "log2FoldChange")]
#     deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
#     deseq2_data <- read.csv(deseq2_file)[, c("gene", "log2FoldChange")]
#     deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
# 
# 
#     # Rename columns for clarity
#     colnames(deseq2rr_data)[2] = "logFC_DESeq2_Wald_filtered"
#     colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
#     colnames(deseq2_data)[2] = "logFC_DESeq2_Wald_unfiltered"
#     colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
# 
# 
#     # Merge all data frames by gene name to build dataframe of genes x logFC_method
#     merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(deseq2rr_data, deseq2LRT_data, deseq2_data, deseq2uLRT_data))
# 
#     return(merged_data)
# 
#   } else {
#     deseq2rr_file <- paste0(deseq2rr_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#     deseq2_file <- paste0(deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
# 
#     # Read the gene column
#     deseq2rr_genes = read.csv(deseq2rr_file)$gene
#     deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
#     deseq2_genes = read.csv(deseq2_file)$gene
#     deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
# 
# 
# 
#     return(list(Wald_filtered = deseq2rr_genes, LRT_filtered = deseq2LRT_genes, Wald_unfiltered = deseq2_genes, LRT_unfiltered = deseq2uLRT_genes))
#   }
# }

read_gene_list <- function(cell_type, return_logFC=FALSE) {
  if (return_logFC) {
    cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    
    # Read the full data for log fold changes
    deseq2rr_data <- read.csv(cb_deseq2_file)[, c("gene", "log2FoldChange")]
    deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
    cb_deseq2_fltrd_data <- read.csv(cb_deseq2_fltrd_file)[, c("gene", "log2FoldChange")]
    deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
    
    
    # Rename columns for clarity
    colnames(cb_deseq2_data)[2] = "logFC_DESeq2_LRT_filtered_Cellbender"
    colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
    colnames(cb_deseq2_fltrd_data)[2] = "logFC_DESeq2_LRT_unfiltered_Cellbender"
    colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
    
    
    # Merge all data frames by gene name to build dataframe of genes x logFC_method
    merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(cb_deseq2_data, deseq2LRT_data, cb_deseq2_fltrd_data, deseq2uLRT_data))
    
    return(merged_data)
    
  } else {
    cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    
    # Read the gene column
    cb_deseq2_genes = read.csv(cb_deseq2_file)$gene
    deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
    cb_deseq2_fltrd_genes = read.csv(cb_deseq2_fltrd_file)$gene
    deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
    
    
    
    return(list(CB_LRT_unfiltered = cb_deseq2_genes, LRT_filtered = deseq2LRT_genes, CB_LRT_filtered = cb_deseq2_fltrd_genes, LRT_unfiltered = deseq2uLRT_genes))
  }
}

# read_gene_list <- function(cell_type, return_logFC=FALSE) {
#   if (return_logFC) {
#     cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#     cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#
#     # Read the full data for log fold changes
#     deseq2rr_data <- read.csv(cb_deseq2_file)[, c("gene", "log2FoldChange")]
#     deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
#     cb_deseq2_fltrd_data <- read.csv(cb_deseq2_fltrd_file)[, c("gene", "log2FoldChange")]
#     deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
#
#
#     # Rename columns for clarity
#     colnames(cb_deseq2_data)[2] = "logFC_DESeq2_LRT_filtered_Cellbender"
#     colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
#     colnames(cb_deseq2_fltrd_data)[2] = "logFC_DESeq2_LRT_unfiltered_Cellbender"
#     colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
#
#
#     # Merge all data frames by gene name to build dataframe of genes x logFC_method
#     merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(cb_deseq2_data, deseq2LRT_data, cb_deseq2_fltrd_data, deseq2uLRT_data))
#
#     return(merged_data)
#
#   } else {
#     cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#     cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#
#     # Read the gene column
#     cb_deseq2_genes = read.csv(cb_deseq2_file)$gene
#     deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
#     cb_deseq2_fltrd_genes = read.csv(cb_deseq2_fltrd_file)$gene
#     deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
#
#
#
#     return(list(CB_LRT_unfiltered = cb_deseq2_genes, LRT_filtered = deseq2LRT_genes, CB_LRT_filtered = cb_deseq2_fltrd_genes, LRT_unfiltered = deseq2uLRT_genes))
#   }
# }




read_gene_list <- function(cell_type, return_logFC=FALSE) {
  if (return_logFC) {
    crv9_deseq2_file <- paste0(crv9_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    crv9_deseq2_fltrd_file <- paste0(crv9_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
    
    # Read the full data for log fold changes
    deseq2rr_data <- read.csv(cb_deseq2_file)[, c("gene", "log2FoldChange")]
    deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
    cb_deseq2_fltrd_data <- read.csv(cb_deseq2_fltrd_file)[, c("gene", "log2FoldChange")]
    deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
    
    
    # Rename columns for clarity
    colnames(cb_deseq2_data)[2] = "logFC_DESeq2_LRT_filtered_Cellbender"
    colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
    colnames(cb_deseq2_fltrd_data)[2] = "logFC_DESeq2_LRT_unfiltered_Cellbender"
    colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
    
    
    # Merge all data frames by gene name to build dataframe of genes x logFC_method
    merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(cb_deseq2_data, deseq2LRT_data, cb_deseq2_fltrd_data, deseq2uLRT_data))
    
    return(merged_data)
    
  } else {
    cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
    
    # Read the gene column
    cb_deseq2_genes = read.csv(cb_deseq2_file)$gene
    deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
    cb_deseq2_fltrd_genes = read.csv(cb_deseq2_fltrd_file)$gene
    deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
    
    
    
    return(list(CB_LRT_unfiltered = cb_deseq2_genes, LRT_filtered = deseq2LRT_genes, CB_LRT_filtered = cb_deseq2_fltrd_genes, LRT_unfiltered = deseq2uLRT_genes))
  }
}

# Function to create a Venn diagram for a cell type
# colors <- brewer.pal(3, "Pastel1")
colors <- brewer.pal(4, "Pastel1") # to compare 4 outputs


if (!dir.exists(paste0(output_fig_dir,'venn_diagrams/'))) { dir.create(paste0(output_fig_dir,'venn_diagrams/'),
                                                                       recursive=TRUE)}


library(VennDiagram)

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
                           filename = paste0(output_fig_dir,'venn_diagrams/',cell_type,'_venn_diagram_signif_genes.png'),
                           output=TRUE,
                           category.names = c('CB_unfiltered', 'CR_filtered','CB_filtered', 'CR_unfiltered'),
                           fill = colors,
                           lwd = 2,
                           lty = c('blank', 'blank', 'blank', 'blank'),
                           # cex = .4,
                           fontface = "bold",
                           fontfamily = "sans",
                           # cat.cex = 0.5,
                           cat.fontface = "bold",
                           cat.default.pos = "outer",
                           cat.pos = c(-27, 27, 135, 135),
                           cat.dist = c(0.055, 0.055, 0.085, 0.085),
                           cat.fontfamily = "sans",
                           # rotation = 1,
                           label.col='black',
                           # cex=c(0.3,15)
                           # cat.cex=c(0.7,4)
  )
  
  # Update labels with gene counts and names for each section
  # for (i in 1:length(vennplot)) {
  #   if (!is.null(vennplot[[i]]$label)) {
  #     vennplot[[i]]$label <- formatted_intersections[[i]]
  #   }
  # }
  
  return(vennplot)
}


### grab celltypes from results ----
### directory of differential expression results 
de_directory = paste0(save_deseq2_dir,'results/')

### grab de result file names and extract cell types 
de_results = list.files(path=de_directory, pattern='all_genes')

### List of cell types (replace with your actual cell types)
celltypes = sub(pattern='_condition.*',replacement='',x=de_results)

# Combine all significant genes across all 4 versions of DE analysis done on cell ranger data----
signif_genes_by_celltype = list()

for (celltype in celltypes) {
  gene_lists = read_gene_list(celltype)
  signif_genes_by_celltype[[celltype]] = unique(unlist(gene_lists))
}
signif_genes = unique(unlist(signif_genes_by_celltype))


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
    labs(title = paste("Venn Diagram for", celltype))
  
  # Save the plot to file
  ggsave(paste0('3_weeks/stats/DESeq2_unfiltered_LRT/venn_diagrams/', celltype, '_ggvenn_diagram.png'), plot = p)
}



# try another way
cell_type = 'Endothelial'
gene_lists <- read_gene_list(cell_type)
intersections = calculate.overlap(gene_lists)
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
if (!dir.exists("3_weeks/stats/DESeq2_unfiltered_LRT/logFC_correlation/")) { 
  dir.create("3_Weeks/stats/DESeq2_unfiltered_LRT/logFC_correlation/",
             recursive=TRUE)}

library(corrplot)
celltype = 'Overall'
merged_de_res = read_gene_list(celltype,return_logFC = TRUE)

# Calculate correlation matrix
# Use pairwise complete observations to handle missing values
correlation_matrix = cor(merged_de_res[, -1], use = "pairwise.complete.obs")

# Print correlation matrix
print(correlation_matrix)

png(paste0(save_deseq2_dir,'/logFC_correlation/Overall_logFC_corrplot.png'),
    res=320, height = 10, width = 10, units = "in")
corrplot(correlation_matrix, col=viridis(option='viridis',n=256), type='upper',
         tl.cex=0.8, cl.cex=0.8)
dev.off()


# Optionally, visualize the correlation matrix
png(paste0(save_deseq2_dir,'/logFC_correlation/Overall_logFC_correlation_matrix.png'),
    res=320, height = 10, width = 10, units = "in")
pheatmap(correlation_matrix, 
         color=viridis(option='viridis',n=256),
         fontsize=12)
dev.off()

for (i in 1:length(celltypes)) {
  print(i)
  celltype=celltypes[i]
  cat('\nCalculating logFC correlations for',celltype)
  
  merged_de_res = read_gene_list(celltype, return_logFC = TRUE)
  
  # Calculate correlation matrix
  # Use pairwise complete observations to handle missing values
  correlation_matrix = cor(merged_de_res[, -1], use = "pairwise.complete.obs")
  
  # Print correlation matrix
  cat('\nlogFC correlations for',celltype,'\n')
  print(correlation_matrix)
  
  
  # Optionally, visualize the correlation matrix
  png(paste0(save_deseq2_dir,'/logFC_correlation/',celltype,'_logFC_corrplot.png'),
      res=320, height = 10, width = 10, units = "in")
  corrplot(correlation_matrix, col=viridis(option='viridis',n=256),type='upper',
           tl.cex=0.8, cl.cex=0.8)
  dev.off()
  
  
  png(paste0(save_deseq2_dir,'/logFC_correlation/',celltype,'_logFC_correlation_matrix.png'),
      res=320, height = 10, width = 10, units = "in")
  pheatmap(correlation_matrix, color=viridis(option='viridis',n=256),
           fontsize=12)
  dev.off()
}

# sessionInfo() ----
sessionInfo()




