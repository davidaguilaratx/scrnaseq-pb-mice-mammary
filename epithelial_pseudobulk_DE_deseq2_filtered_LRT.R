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
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/epithelial/figures/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/data/"
# stats_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2/"
# # stats_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2_gene_fltrd/"
# 
# # Create directories to save results if they don't already exist:
# dirs = c(output_data_dir, output_fig_dir, stats_dir)
# 
# for (dir in dirs) {
#   if (!dir.exists(dir)) { dir.create(dir,
#                                      recursive=TRUE) }
# }


nthreads = 2

analysis <- 'cellbender_analysis'
## directories for integrated data ----
fig_dir = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")
stats_dir = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/DESeq2/")

# Create directories to save results if they don't already exist:
dirs = c(fig_dir, data_dir, annot_dir, stats_dir)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# timepoints. 3 weeks, 10 weeks, 7 months, 18 months
timepoints = c('wk3','wk10','mn7','mn18')



# seurat = readRDS(paste0(data_dir_integrated,'de_seurat_annotated.rds'))

epithelial = qs_read(paste0(data_dir,'epithelial.qs'), nthreads = 2)

DefaultAssay(epithelial) = 'RNA'
epithelial = DietSeurat(epithelial,assays = 'RNA')


epithelial$celltype = epithelial$subtype

# diagnostic How many samples per condition per timepoint?
for (time in timepoints) {
  test = subset(epithelial, subset=(timepoint == time))
  
  print(table(test$sample, test$subtype))
}

zero.count.props = list()
# for (time in timepoints[2:4]) {
for (time in timepoints) {
  
  
  if (!dir.exists(paste0(stats_dir,'epithelial/figures'))) { dir.create(paste0(stats_dir,'epithelial/figures'),
                                                              recursive=TRUE) }
  if (!dir.exists(paste0(stats_dir,'epithelial/results'))) { dir.create(paste0(stats_dir,'epithelial/results'),
                                                              recursive=TRUE) }
  
  seurat = subset(epithelial, subset=(timepoint == time))
  
  
  # How many cells per condition per cell type? ----
  cat('\n',time,'# of cell by celltype and condition\n')
  print(table(seurat$celltype, seurat$condition))
  
  
  # How many cells will each pseudobulk sample have? ----
  cat('\n',time,'# of cell by celltype and sample\n')
  print(table(seurat$celltype, seurat$sample))
  
  seurat$celltype = Idents(seurat) 
  
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
  
  ## Check the counts matrix
  dim(counts(sce))

  
  # counts(sce)[1:6, 1:6]
  
  # Extract unique names of clusters (= levels of cluster_id factor variable)
  cluster_names <- levels(colData(sce)$celltype)
  cluster_names
  
  
  # Total number of clusters
  length(cluster_names)
  
  
  # Extract unique names of samples (= levels of sample_id factor variable)
  sample_names <- unique(colData(sce)$sample)
  sample_names

  # Total number of samples
  length(sample_names)

  # Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
  groups <- colData(sce)[, c("celltype", "sample")]
  groups$sample = as.factor(groups$sample) # coerce sample to factor data type
  
  dim(groups)

  head(groups)
  
  
  # Aggregate across cluster-sample groups ----
  # transposing row/columns to have cell_ids as row names matching those of groups
  aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                  groupings = groups, fun = "sum") 
  
  # Explore output matrix
  class(aggr_counts)
  
  dim(aggr_counts)

  # aggr_counts[1:6, 1:6]
  
  # Transpose aggregated matrix to have genes as rows and samples as columns
  aggr_counts <- t(aggr_counts)
  # aggr_counts[1:6, 1:6]
  
  dim(aggr_counts)
  
  
  # Understanding tstrsplit()
  
  # Exploring structure of function output (list) ----
  # tstrsplit(colnames(aggr_counts), "_") %>% str()
  
  ## Comparing the first 10 elements of our input and output strings
  # head(colnames(aggr_counts), n = 20)
  # head(tstrsplit(colnames(aggr_counts), "_")[[1]], n = 20)
  # 
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
  cat('\nStrucuture of counts_ls for',time,'\n')
  str(counts_ls)
  cat('\n\n')
  
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
  
  head(metadata)
  
  # Rename rows to sample name
  rownames(metadata) <- metadata$sample
  head(metadata)
  
  # Createa table with the number of cells per sample and cluster
  t <- table(colData(sce)$sample,
             colData(sce)$cluster_id)
  # t[1:6, 1:6]
  
  # Creating metadata list by appending cell count info to metadata table ----
  # generates one metadata data frame specific of each cell type
  
  ## Initiate empty list
  metadata_ls <- list()
  
  for (i in 1:length(counts_ls)) {
    # Get current cluster name
    current_cluster <- names(counts_ls)[i]
    cat("Processing metadata for cluster:", current_cluster, "\n")
    
    # Skip if cluster data is not a proper matrix with column names
    if (!inherits(counts_ls[[i]], "dgCMatrix") && !is.matrix(counts_ls[[i]])) {
      cat("Skipping", current_cluster, "because it's not a proper matrix\n")
      next
    }
    
    # Skip if no column names
    col_names <- colnames(counts_ls[[i]])
    if (is.null(col_names) || length(col_names) == 0) {
      cat("Skipping", current_cluster, "because it has no column names\n")
      next
    }
    
    # Create metadata for this cluster
    df <- data.frame(cluster_sample_id = col_names, stringsAsFactors = FALSE)
    
    # Use base R function strsplit instead of tstrsplit
    split_cols <- strsplit(df$cluster_sample_id, "_")
    
    # Extract cluster_id and sample safely
    df$cluster_id <- sapply(split_cols, function(x) {
      if (length(x) >= 1) return(x[1]) else return(NA)
    })
    
    df$sample <- sapply(split_cols, function(x) {
      if (length(x) >= 3) return(paste0(x[2], "_", x[3])) else return(NA)
    })
    
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
    
    # Check if parsing worked correctly
    if (any(is.na(df$cluster_id)) || any(is.na(df$sample))) {
      cat("Skipping", current_cluster, "due to parsing issues\n")
      print(head(df))
      next
    }
    
    # The rest of your code for cell counts and joining metadata
    # ...
    
    # Join with the main metadata to get condition information
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df <- plyr::join(df, metadata,
                     by = intersect(names(df), names(metadata)))
    
    # Check if the merge worked correctly
    if (any(is.na(df$condition))) {
      cat("Skipping", current_cluster, "because some samples couldn't be matched to conditions\n")
      print(head(df))
      next
    }
    
    # Set row names to match column names in count matrix
    rownames(df) <- df$cluster_sample_id
    
    # Store in metadata list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- current_cluster
  }
  
  # Clean up the lists to remove any skipped clusters
  valid_clusters <- names(metadata_ls)[!sapply(metadata_ls, is.null)]
  counts_ls <- counts_ls[valid_clusters]
  metadata_ls <- metadata_ls[valid_clusters]
  
  # Update cluster_names to only include valid clusters
  cluster_names <- valid_clusters
  cat('\n',time,'has suitable clusters/celltypes:\n')
  cat(cluster_names,'\n',sep='    ')
  
    # for (i in 1:length(counts_ls)) {
  #   
  #   ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
  #   df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
  #   
  #   ## Use tstrsplit() to separate cluster (cell type) and sample IDs
  #   df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
  #   df$sample  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
  #   
  #   ## Retrieve cell count information for this cluster from global cell count table
  #   idx <- which(colnames(t) == unique(df$cluster_id))
  #   cell_counts <- t[, idx]
  #   
  #   ## Remove samples with zero cell contributing to the cluster
  #   cell_counts <- cell_counts[cell_counts > 0]
  #   
  #   ## Match order of cell_counts and sample_ids
  #   sample_order <- match(df$sample, names(cell_counts))
  #   cell_counts <- cell_counts[sample_order]
  #   
  #   ## Append cell_counts to data frame
  #   df$cell_count <- cell_counts
  #   
  #   
  #   ## Join data frame (capturing metadata specific to cluster) to generic metadata
  #   df <- plyr::join(df, metadata, 
  #                    by = intersect(names(df), names(metadata)))
  #   
  #   ## Update rownames of metadata to match colnames of count matrix, as needed later for DE using DeSeq2
  #   rownames(df) <- df$cluster_sample_id
  #   
  #   ## Store complete metadata for cluster i in list
  #   metadata_ls[[i]] <- df
  #   names(metadata_ls)[i] <- unique(df$cluster_id)
  #   
  # }
  
  # Explore the different components of the list
  cat('\nStrucuture of metadata_ls for',time,'\n')
  str(metadata_ls)
  cat('\n\n')
  
  # Double-check that both lists have same names
  all(names(counts_ls) == names(metadata_ls))
  
  
  # DESeq2 analysis ----        
  
  groups_bulk <- colData(sce)[, c("sample")]
  groups_bulk = as.factor(groups_bulk) # coerce sample to factor data type
  
  length(groups_bulk)
  
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
  cat('\n;',time,'zero-count proportion:\n')
  zeroprop = (sum(rowSums(aggr_counts_bulk == 0))) / (nrow(aggr_counts_bulk) * ncol(aggr_counts_bulk))
  print(zeroprop)
  zero.count.props[[paste0(time,'_overall')]] = zeroprop
  
  
  
  # Create DESeq2 object        
  dds_bulk <- DESeqDataSetFromMatrix(aggr_counts_bulk, 
                                     colData = metadata_bulk, 
                                     design = ~ condition)
  
  # prefiltering low count genes.
  # Want genes that have at least 5 counts in at least 4 samples.
  # following recommendation from
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  
  # smallestGroupSize = floor(0.5 * min(table(colData(dds_bulk)$condition)))
  smallestGroupSize = 4
  # cat('\nSize of smaller group is',smallestGroupSize)
  cat('\n# of samples per condition for',time,'overall\n')
  print(table(colData(dds_bulk)$condition))

  keep = rowSums(counts(dds_bulk) >= 5) >= smallestGroupSize
  
  cat('\n# of genes filtered and kept:\n')
  print(table(keep))

  dds_bulk = dds_bulk[keep,]
  
  # zinb_bulk = zinbwave(dds_bulk, K=0, observationalWeights=TRUE,
  #                      BPPARAM=BiocParallel::SerialParam(), epsilon=1e12)
  
  
  # Transform counts for data visualization
  rld_bulk <- rlog(dds_bulk, blind=TRUE)
  # Plot PCA
  
  # Generate QC plots - only use condition for PCA
  tryCatch({
    pca_data_condition = DESeq2::plotPCA(rld_bulk, intgroup = "condition", returnData= TRUE)
    ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = condition, label = name)) +
      geom_point() + 
      geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
      theme_classic2() +
      xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 
    ggsave(paste0(stats_dir,"epithelial/figures/",time,"_overall_specific_PCAplot.png"),dpi=320,height=10,width=10)
  }, error = function(e) {
    cat("Error in plotPCA across condition:", conditionMessage(e), "\n")
  })

  tryCatch({
    pca_data_n_cells = DESeq2::plotPCA(rld_bulk, intgroup = "cell_count", returnData=TRUE)
    ggplot(pca_data_n_cells, aes(x = PC1, y = PC2, color = cell_count, label = name)) +
      geom_point() + 
      geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
      theme_classic2() +
      xlab(paste0("PC1: ", round(attr(pca_data_n_cells, "percentVar")[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(attr(pca_data_n_cells, "percentVar")[2] * 100), "% variance"))
    ggsave(paste0(stats_dir,"epithelial/figures/",time,"_overall_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
  }, error = function(e) {
    cat("Error in plotPCA across condition, including cell count:", conditionMessage(e), "\n")
  })
  

  # DESeq2::plotPCA(rld_bulk, ntop = 500, intgroup = "cell_count")
  # ggsave(paste0(stats_dir,"epithelial/figures/",time,"_overall_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat_bulk <- assay(rld_bulk)
  rld_cor_bulk <- cor(rld_mat_bulk)
  
  # Plot heatmap
  tryCatch({
    png(paste0(stats_dir,'epithelial/figures/',time,'_overall_specific_heatmap.png'),
        res=320, height = 10, width = 10, units = "in")
    pheatmap(rld_cor_bulk, annotation = metadata_bulk[, c("condition"), drop=F])
    if (dev.cur() > 1) {dev.off()}
  }, error = function(e) {
    cat("Error in creating heatmap:", conditionMessage(e), "\n")
    if (dev.cur() > 1) dev.off()
  })
  
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
  tryCatch({
    png(paste0(stats_dir,'epithelial/figures/',time,'_overall_dispersion_plot.png'),
        res=320, height = 5, width = 6, units = "in")
    plotDispEsts(dds_bulk)
    if (dev.cur() > 1) {dev.off()}
  }, error = function(e) {
    cat("Error in creating dispersion plot:", conditionMessage(e), "\n")
    if (dev.cur() > 1) dev.off()
  })
  
  # Check the coefficients for the comparison
  resultsNames(dds_bulk)
  
  # Generate results object
  res_bulk <- results(dds_bulk, 
                      name = "condition_pb_vs_ctrl",
                      # name = "conditionpb",
                      alpha = 0.05,
                      independentFiltering = TRUE)
  
  # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
  res_bulk <- lfcShrink(dds_bulk, 
                        coef = "condition_pb_vs_ctrl",
                        # coef = "conditionpb",
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
            paste0(stats_dir,"epithelial/results/overall_condition_", 
                   unique(metadata_bulk$condition)[2], "_vs_", unique(metadata_bulk$condition)[1],'_',time, "_all_genes_deseq2.csv"),
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
            paste0(stats_dir,"epithelial/results/overall_condition_", 
                   unique(metadata_bulk$condition)[2], "_vs_", unique(metadata_bulk$condition)[1],'_',time,"_signif_genes_deseq2.csv"),
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
  
  # Look at distributions for up to top 20 significant genes if any
  if ((nrow(sig_res) > 0) & (!is.null(nrow(sig_res)))) {
    ## Extract top 20 DEG from resLFC (make sure to order by padj)
    n = min(20, nrow(sig_res)) # top 20 or fewer DE genes
    top20_sig_genes <- sig_res_bulk %>%
      dplyr::arrange(padj) %>%
      dplyr::pull(gene) %>%
      head(n = min(20,nrow(sig_res_bulk)))
  
    print(head(top20_sig_genes))
    
    ## Extract matching normalized count values from matrix
    top20_sig_counts <- normalized_counts[rownames(normalized_counts) %in% top20_sig_genes, ]
    top20_sig_counts
    
    ## Convert wide matrix to long data frame for ggplot2
    tryCatch({
      top20_sig_df <- data.frame(top20_sig_counts)
      top20_sig_df$gene <- rownames(top20_sig_counts)
      
      # Check if 'gene' is actually in the dataframe
      if (!"gene" %in% colnames(top20_sig_df)) {
        stop("gene' column not found in top20_sig_df\n")
      }
      
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
      ggsave(paste0(stats_dir,'epithelial/figures/',time,'_overall_topDEgenes.png'),dpi=320,width=20,height=10)
    
    }, error = function(e) {
      cat("Error in creating top DE genes plot:", conditionMessage(e), "\n")
    })
    
    #### Heatmap of overall significant genes ----
    
    ## Extract normalized counts for significant genes only
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res_bulk$gene, ]
    
    ## Run pheatmap using the metadata data frame for the annotation
    png(filename=paste0(stats_dir,'epithelial/figures/',time,'_overall_pheatmap_top17DEgene.png'), units='in', width=10, height=10, res=320)
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
    if (dev.cur() > 1) {dev.off()}
  
    # Volcano plot
    # res_table_thres <- res_tbl_bulk[!is.na(res_tbl_bulk$padj), ] %>% 
    #   mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= log2fc_cutoff)
    # min(log10(res_table_thres$padj))
    
    ## Generate volcano plot
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

  }
  ## generate enhanced volcano plot
  tryCatch({
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
  ggsave(paste0(stats_dir,'epithelial/figures/',time,'_overall_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
  }, error = function(e) {
    cat("Error in creating volcano plot:", conditionMessage(e), "\n")
  })
  
  cat('\nFinished analyzing',time, 'overall\n')
  
  
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
  # zero_props = list()
  # zero_props_post_filtering = list()
  # change_in_zero_props = list()
  
  #### DE analysis by cell type ----
  get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
    cat('Beginning cluster ', clustx, 'for timepoint:',time,'\n\n')
    
    # Check if cluster exists in our valid list
    if (!(clustx %in% names(counts_ls)) || !(clustx %in% names(metadata_ls))) {
      cat('Skipping cluster', clustx, 'because it was not properly processed\n')
      return(NULL)
    }
    
    # Extract counts matrix and metadata for cluster x
    idx <- which(names(counts_ls) == clustx)
    cluster_counts <- counts_ls[[idx]]
    cluster_metadata <- metadata_ls[[idx]]
    
    # Make sure condition is a factor
    cluster_metadata$condition <- as.factor(cluster_metadata$condition)
    
    # Check if we have samples from both conditions
    condition_counts <- table(cluster_metadata$condition)
    if (length(condition_counts) < 2 || any(condition_counts == 0)) {
      cat('Skipping cluster', clustx, 'because it does not have samples from both conditions\n')
      cat('Condition counts:', paste(names(condition_counts), condition_counts, sep="=", collapse=", "), '\n')
      return(NULL)
    }
    
    
    # Print error message if sample names do not match
    if ( all(colnames(cluster_counts) != rownames(cluster_metadata)) ) {
      print("ERROR: sample names in counts matrix columns and metadata rows do not match!")
    }
    
    # zero-count proportions?
    cat('\n;',time,clustx,'zero-count proportion:\n')
    zeroprop = (sum(rowSums(cluster_counts == 0))) / (nrow(cluster_counts) * ncol(cluster_counts))
    print(zeroprop)
    zero.count.props[[paste0(time,'_',clustx)]] = zeroprop

    
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata, 
                                  design = ~ condition)
    
    # prefiltering low count genes. 
    # Want genes that have at least 10 counts in at least 4 samples. 6 for cellbender data since we haven't dropped any samples
    # following recommendation from
    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
    # smallestGroupSize = floor(0.5 * min(table(colData(dds)$condition)))
    smallestGroupSize = 4
    # cat('\nSize of smaller group is',smallestGroupSize)
    cat('\n# of samples per condition for',time,clustx,'\n')
    print(table(colData(dds)$condition))
    
    keep = rowSums(counts(dds_bulk) >= 5) >= smallestGroupSize
    
    cat('\n# of genes filtered and kept:\n')
    print(table(keep))
   
    dds = dds[keep,]
    
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
    tryCatch({
      pca_data_condition = DESeq2::plotPCA(rld, intgroup = "condition", returnData= TRUE)
      ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = condition, label = name)) +
        geom_point() + 
        geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
        theme_classic2() +
        xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 
      ggsave(paste0(stats_dir,"epithelial/figures/",time,"_overall_specific_PCAplot.png"),dpi=320,height=10,width=10)
    }, error = function(e) {
      cat("Error in plotPCA across condition:", conditionMessage(e), "\n")
    })
    
    tryCatch({
      pca_data_n_cells = DESeq2::plotPCA(rld, intgroup = "cell_count", returnData=TRUE)
      ggplot(pca_data_n_cells, aes(x = PC1, y = PC2, color = cell_count, label = name)) +
        geom_point() + 
        geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
        theme_classic2() +
        xlab(paste0("PC1: ", round(attr(pca_data_n_cells, "percentVar")[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(attr(pca_data_n_cells, "percentVar")[2] * 100), "% variance"))
      ggsave(paste0(stats_dir,"epithelial/figures/",time,"_overall_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
    }, error = function(e) {
      cat("Error in plotPCA across condition, including cell count:", conditionMessage(e), "\n")
    })
    
    ## Extract rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    ## Plot and save heatmap
    tryCatch({
      png(paste0(stats_dir,"/epithelial/figures/", clustx,'_',time, "_specific_heatmap.png"),
          height = 6, width = 7.5, units = "in", res = 300)
      pheatmap(rld_cor, annotation = cluster_metadata[, c("condition"), drop = FALSE])
      if (dev.cur() > 1) {dev.off()}
    }, error = function(e) {
      cat("Error in creating heatmap:", conditionMessage(e), "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    # compute size factors using scran::computeSumFactors() as recommended for
    # single cell analysis on deseq2 vignette here:
    # https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#recommendations-for-single-cell-analysis
    scr = scran::computeSumFactors(dds) # scran size factors
    
    sizeFactors(dds) = sizeFactors(scr)
    
    # Run DESeq2 differential expression analysis
    dds <- DESeq(dds,
                 test='LRT',
                 fitType='parametric', # the default
                 reduced= ~1,
                 minReplicatesForReplace=Inf,
                 useT=TRUE,
                 minmu=1e-6)
    
    cat('Finished fitting DESeq2 model for',time,clustx,'\n')
    
    ## Plot dispersion estimates
    tryCatch({
    png(paste0(stats_dir,"/epithelial/figures/", clustx,'_',time, "_dispersion_plot.png"),
        height = 5, width = 6, units = "in", res = 300)
    plotDispEsts(dds)
    if (dev.cur() > 1) {dev.off()}
    }, error = function(e) {
      cat("Error in creating dispersion plot:", conditionMessage(e), "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    
    ## Output and shrink results of Wald test for contrast A vs B
    contrast <- paste(c("condition", A, "vs", B), collapse = "_")
    print(resultsNames(dds))
    
    res <- results(dds, name = contrast, alpha = 0.05, independentFiltering = TRUE)
    # res <- results(dds, name = 'conditionpb', alpha = 0.05)
    res <- lfcShrink(dds, coef = contrast, res = res, type = "apeglm")
    # res <- lfcShrink(dds, coef = 'conditionpb', res = res, type = "apeglm")
    
    ## Turn the results object into a tibble for use with tidyverse functions
    res_tbl <- res %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble()
    
    write.csv(res_tbl,
              paste0(stats_dir,"/epithelial/results/", clustx, "_", contrast,'_',time, "_all_genes_deseq2.csv"),
              quote = FALSE, 
              row.names = FALSE)
    
    ## Subset the significant results
    sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
      dplyr::arrange(padj)
    
    write.csv(sig_res,
              paste0(stats_dir,"/epithelial/results/", clustx, "_", contrast,'_',time, "_signif_genes_deseq2.csv"),
              quote = FALSE, 
              row.names = FALSE)
    
    # Generate results visualization plots
    
    ## Extract normalized counts from dds object
    normalized_counts <- counts(dds, normalized = TRUE)
    
    # Look at distributions for up to top 20 significant genes if any
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
      ggsave(paste0(stats_dir,"/epithelial/figures/", clustx, "_", contrast,'_',time, "_top20_DE_genes.png"))
      
      ## Extract normalized counts for significant genes only
      sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
      
      ## generate pheatmap
      if ((nrow(sig_res) >= 2) & (!is.null(nrow(sig_res)))) {
        png(filename=paste0(stats_dir,'/epithelial/figures/',clustx,'_',contrast,'_',time,'_pheatmap_top17DEgene.png'), units='in', width=10, height=10, res=320)
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
        if (dev.cur() > 1) {dev.off()}
      } else {
        message("Skipping pheatmap: fewer than 2 genes found (nrow < 2).")
      }
      
      
    } else {
      warning(paste0("No significant genes found for ",clustx))
      cat(paste0("No significant genes found for ",clustx),'\n\n')
    }
    
    
    ## generate enhanced volcano plot
    tryCatch({
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
      ggsave(paste0(stats_dir,'/epithelial/figures/',clustx,'_',contrast,'_',time,'_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
    }, error = function(e) {
      cat("Error in creating volcano plot:", conditionMessage(e), "\n")
    })
    
    cat('\nFinished analyzing',time,clustx,'\n')
  }
  
  # Run the script on all clusters comparing stimulated condition relative to control condition
  map(cluster_names, get_dds_resultsAvsB, A = "pb", B = "ctrl", padj_cutoff = 0.05)
}





  
  
#   # Compare DE results across DESeq2 variations----
#   
#   # Define paths for the directories containing the result files
#   cr_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/'
#   deseq2rr_dir <- paste0(cr_deseq2_dir,'stats_rerun/DESeq2/results/') # filtered wald
#   deseq2LRT_dir <- paste0(cr_deseq2_dir,'stats/DESeq2/results/') # filtered LRT
#   deseq2_dir <- paste0(cr_deseq2_dir,'stats/DESeq2_unfiltered/results/') # unfiltered Wald
#   deseq2uLRT_dir <- paste0(cr_deseq2_dir,'stats/DESeq2_unfiltered_LRT/results/') # unfiltered LRT
#   cb_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2/results/' # cellbender unfiltered LRT
#   cb_deseq2_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/results/'
#   crv9_deseq2_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2_gene_fltrd/results/'
#   crv9_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2/results/'
#   # Function to read gene lists for a cell type from different methods,
#   # will return gene lists only be default but can return logFC too if indicated.
#   # read_gene_list <- function(cell_type, return_logFC=FALSE) {
#   #   if (return_logFC) {
#   #     deseq2rr_file <- paste0(deseq2rr_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #     deseq2_file <- paste0(deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   # 
#   #     # Read the full data for log fold changes
#   #     deseq2rr_data <- read.csv(deseq2rr_file)[, c("gene", "log2FoldChange")]
#   #     deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
#   #     deseq2_data <- read.csv(deseq2_file)[, c("gene", "log2FoldChange")]
#   #     deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
#   # 
#   # 
#   #     # Rename columns for clarity
#   #     colnames(deseq2rr_data)[2] = "logFC_DESeq2_Wald_filtered"
#   #     colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
#   #     colnames(deseq2_data)[2] = "logFC_DESeq2_Wald_unfiltered"
#   #     colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
#   # 
#   # 
#   #     # Merge all data frames by gene name to build dataframe of genes x logFC_method
#   #     merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(deseq2rr_data, deseq2LRT_data, deseq2_data, deseq2uLRT_data))
#   # 
#   #     return(merged_data)
#   # 
#   #   } else {
#   #     deseq2rr_file <- paste0(deseq2rr_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #     deseq2_file <- paste0(deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   # 
#   #     # Read the gene column
#   #     deseq2rr_genes = read.csv(deseq2rr_file)$gene
#   #     deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
#   #     deseq2_genes = read.csv(deseq2_file)$gene
#   #     deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
#   # 
#   # 
#   # 
#   #     return(list(Wald_filtered = deseq2rr_genes, LRT_filtered = deseq2LRT_genes, Wald_unfiltered = deseq2_genes, LRT_unfiltered = deseq2uLRT_genes))
#   #   }
#   # }
#   
#   read_gene_list <- function(cell_type, return_logFC=FALSE) {
#     if (return_logFC) {
#       cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       
#       # Read the full data for log fold changes
#       deseq2rr_data <- read.csv(cb_deseq2_file)[, c("gene", "log2FoldChange")]
#       deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
#       cb_deseq2_fltrd_data <- read.csv(cb_deseq2_fltrd_file)[, c("gene", "log2FoldChange")]
#       deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
#       
#       
#       # Rename columns for clarity
#       colnames(cb_deseq2_data)[2] = "logFC_DESeq2_LRT_filtered_Cellbender"
#       colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
#       colnames(cb_deseq2_fltrd_data)[2] = "logFC_DESeq2_LRT_unfiltered_Cellbender"
#       colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
#       
#       
#       # Merge all data frames by gene name to build dataframe of genes x logFC_method
#       merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(cb_deseq2_data, deseq2LRT_data, cb_deseq2_fltrd_data, deseq2uLRT_data))
#       
#       return(merged_data)
#       
#     } else {
#       cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       
#       # Read the gene column
#       cb_deseq2_genes = read.csv(cb_deseq2_file)$gene
#       deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
#       cb_deseq2_fltrd_genes = read.csv(cb_deseq2_fltrd_file)$gene
#       deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
#       
#       
#       
#       return(list(CB_LRT_unfiltered = cb_deseq2_genes, LRT_filtered = deseq2LRT_genes, CB_LRT_filtered = cb_deseq2_fltrd_genes, LRT_unfiltered = deseq2uLRT_genes))
#     }
#   
#   # read_gene_list <- function(cell_type, return_logFC=FALSE) {
#   #   if (return_logFC) {
#   #     cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #     cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#   #
#   #     # Read the full data for log fold changes
#   #     deseq2rr_data <- read.csv(cb_deseq2_file)[, c("gene", "log2FoldChange")]
#   #     deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
#   #     cb_deseq2_fltrd_data <- read.csv(cb_deseq2_fltrd_file)[, c("gene", "log2FoldChange")]
#   #     deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
#   #
#   #
#   #     # Rename columns for clarity
#   #     colnames(cb_deseq2_data)[2] = "logFC_DESeq2_LRT_filtered_Cellbender"
#   #     colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
#   #     colnames(cb_deseq2_fltrd_data)[2] = "logFC_DESeq2_LRT_unfiltered_Cellbender"
#   #     colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
#   #
#   #
#   #     # Merge all data frames by gene name to build dataframe of genes x logFC_method
#   #     merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(cb_deseq2_data, deseq2LRT_data, cb_deseq2_fltrd_data, deseq2uLRT_data))
#   #
#   #     return(merged_data)
#   #
#   #   } else {
#   #     cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #     deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #     cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #     deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#   #
#   #     # Read the gene column
#   #     cb_deseq2_genes = read.csv(cb_deseq2_file)$gene
#   #     deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
#   #     cb_deseq2_fltrd_genes = read.csv(cb_deseq2_fltrd_file)$gene
#   #     deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
#   #
#   #
#   #
#   #     return(list(CB_LRT_unfiltered = cb_deseq2_genes, LRT_filtered = deseq2LRT_genes, CB_LRT_filtered = cb_deseq2_fltrd_genes, LRT_unfiltered = deseq2uLRT_genes))
#   #   }
#   # }
#   
#   
#   
#   
#   read_gene_list <- function(cell_type, return_logFC=FALSE) {
#     if (return_logFC) {
#       crv9_deseq2_file <- paste0(crv9_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       crv9_deseq2_fltrd_file <- paste0(crv9_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
#       
#       # Read the full data for log fold changes
#       deseq2rr_data <- read.csv(cb_deseq2_file)[, c("gene", "log2FoldChange")]
#       deseq2LRT_data <- read.csv(deseq2LRT_file)[, c("gene", "log2FoldChange")]
#       cb_deseq2_fltrd_data <- read.csv(cb_deseq2_fltrd_file)[, c("gene", "log2FoldChange")]
#       deseq2uLRT_data <- read.csv(deseq2uLRT_file)[, c("gene", "log2FoldChange")]
#       
#       
#       # Rename columns for clarity
#       colnames(cb_deseq2_data)[2] = "logFC_DESeq2_LRT_filtered_Cellbender"
#       colnames(deseq2LRT_data)[2] = "logFC_DESeq2_LRT_filtered"
#       colnames(cb_deseq2_fltrd_data)[2] = "logFC_DESeq2_LRT_unfiltered_Cellbender"
#       colnames(deseq2uLRT_data)[2] = "logFC_DESeq2_LRT_unfiltered"
#       
#       
#       # Merge all data frames by gene name to build dataframe of genes x logFC_method
#       merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(cb_deseq2_data, deseq2LRT_data, cb_deseq2_fltrd_data, deseq2uLRT_data))
#       
#       return(merged_data)
#       
#     } else {
#       cb_deseq2_file <- paste0(cb_deseq2_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       deseq2LRT_file <- paste0(deseq2LRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       cb_deseq2_fltrd_file <- paste0(cb_deseq2_fltrd_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       deseq2uLRT_file <- paste0(deseq2uLRT_dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
#       
#       # Read the gene column
#       cb_deseq2_genes = read.csv(cb_deseq2_file)$gene
#       deseq2LRT_genes = read.csv(deseq2LRT_file)$gene
#       cb_deseq2_fltrd_genes = read.csv(cb_deseq2_fltrd_file)$gene
#       deseq2uLRT_genes = read.csv(deseq2LRT_file)$gene
#       
#       
#       
#       return(list(CB_LRT_unfiltered = cb_deseq2_genes, LRT_filtered = deseq2LRT_genes, CB_LRT_filtered = cb_deseq2_fltrd_genes, LRT_unfiltered = deseq2uLRT_genes))
#     }
#   }
#   
#   # Function to create a Venn diagram for a cell type
#   # colors <- brewer.pal(3, "Pastel1")
#   colors <- brewer.pal(4, "Pastel1") # to compare 4 outputs
#   
#   
#   if (!dir.exists(paste0(output_fig_dir,'venn_diagrams/'))) { dir.create(paste0(output_fig_dir,'venn_diagrams/'),
#                                                                          recursive=TRUE)}
#   
#   
#   library(VennDiagram)
#   
#   create_venn_diagram <- function(gene_lists, cell_type) {
#     
#     # get intersecting gene names and format them
#     intersections = calculate.overlap(gene_lists)
#     formatted_intersections <- lapply(intersections, function(x) {
#       count = length(x)
#       if (count > 0) {
#         paste0(count, "\n", paste(x, collapse = "\n"))
#       } else {
#         "0"
#       }
#     })
#     
#     
#     vennplot <- venn.diagram(gene_lists,
#                              filename = paste0(output_fig_dir,'venn_diagrams/',cell_type,'_venn_diagram_signif_genes.png'),
#                              output=TRUE,
#                              category.names = c('CB_unfiltered', 'CR_filtered','CB_filtered', 'CR_unfiltered'),
#                              fill = colors,
#                              lwd = 2,
#                              lty = c('blank', 'blank', 'blank', 'blank'),
#                              # cex = .4,
#                              fontface = "bold",
#                              fontfamily = "sans",
#                              # cat.cex = 0.5,
#                              cat.fontface = "bold",
#                              cat.default.pos = "outer",
#                              cat.pos = c(-27, 27, 135, 135),
#                              cat.dist = c(0.055, 0.055, 0.085, 0.085),
#                              cat.fontfamily = "sans",
#                              # rotation = 1,
#                              label.col='black',
#                              # cex=c(0.3,15)
#                              # cat.cex=c(0.7,4)
#     )
#     
#     # Update labels with gene counts and names for each section
#     # for (i in 1:length(vennplot)) {
#     #   if (!is.null(vennplot[[i]]$label)) {
#     #     vennplot[[i]]$label <- formatted_intersections[[i]]
#     #   }
#     # }
#     
#     return(vennplot)
#   }
#   
#   
#   ### grab celltypes from results ----
#   ### directory of differential expression results 
#   de_directory = paste0(stats_dir,'results/')
#   
#   ### grab de result file names and extract cell types 
#   de_results = list.files(path=de_directory, pattern='all_genes')
#   
#   ### List of cell types (replace with your actual cell types)
#   celltypes = sub(pattern='_condition.*',replacement='',x=de_results)
#   
#   # Combine all significant genes across all 4 versions of DE analysis done on cell ranger data----
#   signif_genes_by_celltype = list()
#   
#   for (celltype in celltypes) {
#     gene_lists = read_gene_list(celltype)
#     signif_genes_by_celltype[[celltype]] = unique(unlist(gene_lists))
#   }
#   signif_genes = unique(unlist(signif_genes_by_celltype))
#   
#   
#   # Loop over each cell type and create the Venn diagrams
#   for (celltype in celltypes) {
#     gene_lists = read_gene_list(celltype)
#     create_venn_diagram(gene_lists, celltype)
#   }
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   ### try other ways. Want to include gene names in venn diagrams ----
#   library(ggVennDiagram)
#   
#   # Loop over each cell type and create the Venn diagrams with ggVennDiagram
#   for (celltype in celltypes) {
#     gene_lists <- read_gene_list(celltype)
#     
#     # Create the Venn diagram
#     p <- ggVennDiagram(gene_lists, label = "both") + 
#       scale_fill_gradient(low = "lightblue", high = "purple") +
#       theme(legend.position = "none") +
#       labs(title = paste("Venn Diagram for", celltype))
#     
#     # Save the plot to file
#     ggsave(paste0('3_weeks/stats/DESeq2_unfiltered_LRT/venn_diagrams/', celltype, '_ggvenn_diagram.png'), plot = p)
#   }
#   
#   
#   
#   # try another way
#   cell_type = 'Endothelial'
#   gene_lists <- read_gene_list(cell_type)
#   intersections = calculate.overlap(gene_lists)
#   formatted_intersections <- lapply(intersections, function(x) {
#     count = length(x)
#     if (count > 0) {
#       paste0(count, "\n", paste(x, collapse = "\n"))
#     } else {
#       "0"
#     }
#   })
#   vennplot <- venn.diagram(gene_lists[1:3],
#                            # filename=NULL,
#                            filename = paste0('3_weeks/stats/venn_diagrams/deseq2s/',cell_type,'_venn_diagram.png'),
#                            output=TRUE,
#                            category.names = c("DESeq2_filtered", 'DESeq2_zinb','DESeq2'),
#                            fill = colors,
#                            lwd = 2,
#                            lty = 'blank',
#                            # cex = .4,
#                            fontface = "bold",
#                            fontfamily = "sans",
#                            # cat.cex = 0.5,
#                            cat.fontface = "bold",
#                            cat.default.pos = "outer",
#                            cat.pos = c(-27, 27, 135),
#                            cat.dist = c(0.055, 0.055, 0.085),
#                            cat.fontfamily = "sans",
#                            rotation = 1,
#                            label.col='black',
#                            cex=0.7,
#                            cat.cex=0.7)
#   
#   # run if set filename above to NULL
#   # grid.newpage()
#   # grid.draw(vennplot)
#   
#   # Update labels with gene counts and names for each section
#   for (i in 1:length(vennplot)) {
#     if (!is.null(vennplot[[i]]$label)) {
#       vennplot[[i]]$label <- formatted_intersections[[i]]
#     }
#   }
#   
#   
#   
#   vennplot <- venn.diagram(gene_lists,
#                            filename=NULL,
#                            category.names = c("DESeq2", "edgeR", "limma"),
#                            alpha=c(0.5,0.5,0.5),
#                            fill = colors,
#                            lwd = 2,
#                            lty = 'blank',
#                            cex = 1.5,
#                            cat.cex = 1.5,
#                            fontface = "bold",
#                            fontfamily = "sans",
#                            cat.fontface = "bold",
#                            cat.default.pos = "outer",
#                            cat.pos = c(-27, 27, 135),
#                            cat.dist = c(0.055, 0.055, 0.085),
#                            cat.fontfamily = "sans",
#                            rotation = 1,
#                            label.col='black',
#                            print.mode='raw')
#   grid.newpage()
#   grid.draw(vennplot)
#   
#   lapply(vennplot, function(i) i$label)
#   
#   
#   # Calculate overlap site
#   overlaps <- calculate.overlap(gene_lists)
#   overlaps <- overlaps[str_sort(names(overlaps), numeric = TRUE)] # sort base on numeric value
#   
#   
#   # Apply name to global variable
#   # Index of venn diagram start at calculate.overlaps + 8. You have to find the index value by yourself for (3,5,6,7,.. venn)
#   walk2(seq(overlaps) + 8, seq(overlaps),  
#         function(x,y) {vennplot[[x]]$label <<- paste0(overlaps[[y]], collapse = "\n")})
#   
#   # Draw plot
#   grid.draw(v)
#   grid.newpage()
#   grid.draw(vennplot)
#   
#   
#   library(ggVennDiagram)
#   
#   ggVennDiagram(gene_lists,
#                 show_intersect = TRUE,
#                 set_color='black') + 
#     theme(legend.position = "none") +
#     labs(title = "Venn Diagram for Adipocytes")
#   
#   
#   
#   ## Check correlation of log-fold change between all three methods ----
#   if (!dir.exists("3_weeks/stats/DESeq2_unfiltered_LRT/logFC_correlation/")) { 
#     dir.create("3_Weeks/stats/DESeq2_unfiltered_LRT/logFC_correlation/",
#                recursive=TRUE)}
#   
#   library(corrplot)
#   celltype = 'Overall'
#   merged_de_res = read_gene_list(celltype,return_logFC = TRUE)
#   
#   # Calculate correlation matrix
#   # Use pairwise complete observations to handle missing values
#   correlation_matrix = cor(merged_de_res[, -1], use = "pairwise.complete.obs")
#   
#   # Print correlation matrix
#   print(correlation_matrix)
#   
#   png(paste0(stats_dir,'/logFC_correlation/Overall_logFC_corrplot.png'),
#       res=320, height = 10, width = 10, units = "in")
#   corrplot(correlation_matrix, col=viridis(option='viridis',n=256), type='upper',
#            tl.cex=0.8, cl.cex=0.8)
#   if (dev.cur() > 1) {dev.off()}
#   
#   
#   # Optionally, visualize the correlation matrix
#   png(paste0(stats_dir,'/logFC_correlation/Overall_logFC_correlation_matrix.png'),
#       res=320, height = 10, width = 10, units = "in")
#   pheatmap(correlation_matrix, 
#            color=viridis(option='viridis',n=256),
#            fontsize=12)
#   if (dev.cur() > 1) {dev.off()}
#   
#   for (i in 1:length(celltypes)) {
#     print(i)
#     celltype=celltypes[i]
#     cat('\nCalculating logFC correlations for',celltype)
#     
#     merged_de_res = read_gene_list(celltype, return_logFC = TRUE)
#     
#     # Calculate correlation matrix
#     # Use pairwise complete observations to handle missing values
#     correlation_matrix = cor(merged_de_res[, -1], use = "pairwise.complete.obs")
#     
#     # Print correlation matrix
#     cat('\nlogFC correlations for',celltype,'\n')
#     print(correlation_matrix)
#     
#     
#     # Optionally, visualize the correlation matrix
#     png(paste0(stats_dir,'/logFC_correlation/',celltype,'_logFC_corrplot.png'),
#         res=320, height = 10, width = 10, units = "in")
#     corrplot(correlation_matrix, col=viridis(option='viridis',n=256),type='upper',
#              tl.cex=0.8, cl.cex=0.8)
#     if (dev.cur() > 1) {dev.off()}
#     
#     
#     png(paste0(stats_dir,'/logFC_correlation/',celltype,'_logFC_correlation_matrix.png'),
#         res=320, height = 10, width = 10, units = "in")
#     pheatmap(correlation_matrix, color=viridis(option='viridis',n=256),
#              fontsize=12)
#     if (dev.cur() > 1) {dev.off()}
#   }
# }

# sessionInfo() ----
sessionInfo()




