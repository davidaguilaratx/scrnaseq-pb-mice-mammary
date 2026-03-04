# to install Matrix.utils
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", type = "source", repos = NULL)


# Load libraries
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
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


nthreads = 4

analysis <- 'cellbender_analysis'
## directories for integrated data ----
fig_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")
stats_dir = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/DESeq2_unfiltered/")

# Create directories to save results if they don't already exist:
# Create directories to save results if they don't already exist:
dirs = c(fig_dir, data_dir, annot_dir, stats_dir)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# timepoints. 3 weeks, 10 weeks, 7 months, 18 months
timepoints = c('wk3','wk10','mn7','mn18')



# seurat = readRDS(paste0(data_dir_integrated,'de_seurat_annotated.rds'))

de_seurat = qs_read(paste0(data_dir,'de_seurat_annotated.qs'),nthreads=nthreads)


# function for parallelization ----
get_bioc_param <- function(workers = nthreads) {
  # Check if running on Windows
  if (.Platform$OS.type == "windows") {
    message("Running on Windows, using SnowParam")
    return(SnowParam(workers = workers))
  } else {
    # On Unix-like systems
    message("Running on Unix-like system, using MulticoreParam")
    return(MulticoreParam(workers = workers))
  }
}


# diagnostic How many samples per condition per timepoint?
for (time in timepoints) {
  test = subset(de_seurat, subset=(timepoint == time))
  
  print(table(test$sample, test$celltype))
}

# Set thresholds
padj_cutoff = 0.05
log2fc_cutoff = 0.58 # 50% change


zero.count.props = list()
# for (time in timepoints[2:4]) {
for (time in timepoints) {
  
  
  if (!dir.exists(paste0(stats_dir,'figures'))) { dir.create(paste0(stats_dir,'figures'),
                                                                        recursive=TRUE) }
  if (!dir.exists(paste0(stats_dir,'results'))) { dir.create(paste0(stats_dir,'results'),
                                                                        recursive=TRUE) }
  
  seurat = subset(de_seurat, subset=(timepoint == time))
  
  # set cluster names (i.e. celltype, cluster number, etc.)
  cluster_names = as.character(sort(unique(seurat@meta.data[["celltype"]])))
  
  # How many cells per condition per cell type? ----
  cat('\n',time,'# of cell by celltype and condition\n')
  print(table(seurat$celltype, seurat$condition))
  
  # How many cells will each pseudobulk sample have? ----
  cat('\n',time,'# of cell by celltype and sample\n')
  print(table(seurat$celltype, seurat$sample))
  
  seurat$celltype = Idents(seurat) 
  
  # bulk expression analysis across conditions ----
  bulk <- AggregateExpression(
    seurat,
    return.seurat = TRUE,
    assays = "RNA",
    group.by = c("sample", "condition")
  )
  
  # Number of cells by sample
  n_cells <- seurat@meta.data %>% 
    dplyr::count(sample) %>% 
    rename("n"="n_cells")
  n_cells$sample <- str_replace(n_cells$sample, "_", "-")
  
  meta_bulk <- left_join(bulk@meta.data, n_cells)
  rownames(meta_bulk) <- meta_bulk$orig.ident
  bulk@meta.data <- meta_bulk
  
  # Turn condition into a factor
  bulk$condition <- factor(bulk$condition, levels=c("ctrl", "pb"))
  
  # bulk@meta.data %>% head()
  
  # plot number of cells per sample
  ggplot(bulk@meta.data, aes(x=sample, y=n_cells, fill=condition)) +
    geom_bar(stat="identity", color="black") +
    theme_classic2() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(title=paste0(time,' overall # of cells per sample'), x="Sample name", y="Number of cells") +
    geom_text(aes(label=n_cells), vjust=-0.5)
  ggsave(paste0(stats_dir,"figures/","Overall",time,"_ncells_per_sample.png"),dpi=320,height=10,width=10)
  
  # DESeq2 analysis ----        
  
  cluster_counts <- FetchData(bulk, layer="counts", vars=rownames(bulk))
  
  # zero-count proportion?
  cat('\n;',time,'zero-count proportion:\n')
  zeroprop = (sum(rowSums(cluster_counts == 0))) / (nrow(cluster_counts) * ncol(cluster_counts))
  print(zeroprop)
  zero.count.props[[paste0(time,'_overall')]] = zeroprop
  
  
  # Create DESeq2 object        
  dds_bulk <- DESeqDataSetFromMatrix(t(cluster_counts), 
                                     colData = bulk@meta.data, 
                                     design = ~ condition)
  
  # prefiltering low count genes.
  # Want genes that have at least 5 counts in at least as many samples as the smallest group size
  # following recommendation from
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering
  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
  
  # smallestGroupSize = min(table(colData(dds_bulk)$condition))
  # smallestGroupSize = 4
  # # cat('\nSize of smaller group is',smallestGroupSize)
  # cat('\n# of samples per condition for',time,'overall\n')
  # print(table(colData(dds_bulk)$condition))
  # 
  # keep = rowSums(counts(dds_bulk) >= 5) >= smallestGroupSize
  # 
  # cat('\n# of genes filtered and kept:\n')
  # print(table(keep))
  # 
  # dds_bulk = dds_bulk[keep,]
  
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
    ggsave(paste0(stats_dir,"figures/","Overall",time,"_specific_PCAplot.png"),dpi=320,height=10,width=10)
  }, error = function(e) {
    cat("Error in plotPCA across condition:", conditionMessage(e), "\n")
  })
  
  tryCatch({
    pca_data_n_cells = DESeq2::plotPCA(rld_bulk, intgroup = "n_cells", returnData=TRUE)
    ggplot(pca_data_n_cells, aes(x = PC1, y = PC2, color = n_cells, label = name)) +
      geom_point() + 
      geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
      theme_classic2() +
      xlab(paste0("PC1: ", round(attr(pca_data_n_cells, "percentVar")[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(attr(pca_data_n_cells, "percentVar")[2] * 100), "% variance"))
    ggsave(paste0(stats_dir,"figures/","Overall",time,"_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
  }, error = function(e) {
    cat("Error in plotPCA across condition, including cell count:", conditionMessage(e), "\n")
  })
  
  
  # DESeq2::plotPCA(rld_bulk, ntop = 500, intgroup = "cell_count")
  # ggsave(paste0(stats_dir,"figures/",time,"_overall_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
  
  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat_bulk <- assay(rld_bulk)
  rld_cor_bulk <- cor(rld_mat_bulk)
  
  # For nicer plots
  rename_samples <- bulk$sample
  colnames(rld_cor_bulk) <- str_replace_all(colnames(rld_cor_bulk), rename_samples)
  rownames(rld_cor_bulk) <- str_replace_all(rownames(rld_cor_bulk), rename_samples)
  
  anno <- bulk@meta.data %>%
    select(sample, condition) %>% 
    remove_rownames() %>% 
    column_to_rownames("sample")
  
  # Plot heatmap
  tryCatch({
    png(paste0(stats_dir,'figures/',"Overall",time,'_specific_heatmap.png'),
        res=320, height = 10, width = 10, units = "in")
    pheatmap(rld_cor_bulk, annotation_col=anno, annotation_row=anno,
             main = paste0(time,' Overall sample correlations'))
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
  # parallel = TRUE,
  # BPPARAM = bp)
  
  
  #### Plot dispersion estimates ----
  tryCatch({
    png(paste0(stats_dir,'figures/',"Overall",time,'_dispersion_plot.png'),
        res=320, height = 5, width = 6, units = "in")
    plotDispEsts(dds_bulk)
    if (dev.cur() > 1) {dev.off()}
  }, error = function(e) {
    cat("Error in creating dispersion plot:", conditionMessage(e), "\n")
    if (dev.cur() > 1) dev.off()
  })
  
  # Check the coefficients for the comparison
  print(resultsNames(dds_bulk))
  
  # Generate results object
  res_bulk <- results(dds_bulk, 
                      name = "condition_pb_vs_ctrl",
                      # name = "conditionpb",
                      alpha = 0.05,
                      independentFiltering = TRUE)
  
  res_bulk
  # Shrink the log2 fold changes to be more appropriate using the apeglm method - should cite [paper]() when using this method
  res_bulk <- lfcShrink(dds_bulk, 
                        coef = "condition_pb_vs_ctrl",
                        # coef = "conditionpb",
                        res=res_bulk,
                        type = "apeglm")
  # parallel = TRUE,
  # BPPARAM = bp)
  
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
            paste0(stats_dir,"results/overall_condition_", 
                   unique(bulk$condition)[2], "_vs_", unique(bulk$condition)[1],'_',time, "_all_genes_deseq2.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  
  #### table DE results for significant genes only ----
  
  
  # Subset the significant results
  sig_res_bulk <- dplyr::filter(res_tbl_bulk, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  # Check significant genes output
  sig_res_bulk
  
  # Write significant results to file
  write.csv(sig_res_bulk,
            paste0(stats_dir,"results/overall_condition_", 
                   unique(bulk$condition)[2], "_vs_", unique(bulk$condition)[1],'_',time,"_signif_genes_deseq2.csv"),
            quote = FALSE, 
            row.names = FALSE)
  
  # Count significantly up/down genes above threshold
  # n_sig_up <- dplyr::filter(sig_res_bulk, log2FoldChange >= log2fc_cutoff) %>% 
  #   nrow()
  # n_sig_dn <- dplyr::filter(sig_res_bulk, log2FoldChange <= -log2fc_cutoff) %>% 
  #   nrow()
  
  
  #### Scatterplot of top genes across all cells ----
  
  ## Extract normalized counts from dds object
  normalized_counts <- counts(dds_bulk, normalized = TRUE)
  
  # Look at distributions for up to top 20 significant genes if any
  if ((nrow(sig_res_bulk) > 0) & (!is.null(nrow(sig_res_bulk)))) {
    ## Extract top 20 DEG from resLFC (make sure to order by padj)
    n = min(20, nrow(sig_res_bulk)) # top 20 or fewer DE genes
    top20_sig_genes <- sig_res_bulk %>%
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
    
    # Check that gene column exists
    if (!"gene" %in% colnames(top20_sig_df)) {
      stop("Gene column not found in data frame")
    }
    
    # Convert to data.table and melt into long format
    top20_sig_df <- reshape2::melt(top20_sig_df, 
                                   id.vars = c("gene"),
                                   variable.name = "sample",
                                   value.name='value')
    
    ## Foramat sample name to be same as in dds object
    top20_sig_df$sample <- gsub(paste0(time,"\\."),paste0(time,'-'), top20_sig_df$sample)
    
    # Proceed with joining the metadata and plotting
    top20_sig_df <- merge(top20_sig_df, as.data.frame(colData(dds_bulk)) %>% rownames_to_column('rownames'),
                          by.x='sample',
                          by.y='rownames',
                          al.x=TRUE)
    
    ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
      geom_boxplot(outlier.shape=NA, width=0.2)+
      geom_jitter(height = 0, width = 0.5)+
      scale_y_continuous(trans = 'log10') +
      ylab("log10 of normalized expression level") +
      xlab("condition") +
      ggtitle(paste0("Top ", n, " Significant DE Genes")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ gene)
    ggsave(paste0(stats_dir,"figures/", "Overall_pb_vs_ctrl_",time,"_top_DE_genes.png"))
    
    ## Extract normalized counts for significant genes only
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res_bulk$gene, ]
    
    if (is.matrix(sig_counts) || is.data.frame(sig_counts)) {
      colnames(sig_counts) = bulk$sample
    }
    
    ## generate pheatmap
    if ((nrow(sig_res_bulk) >= 2) & (!is.null(nrow(sig_res_bulk)))) {
      png(filename=paste0(stats_dir,'figures/',"Overall_pb_vs_ctrl_",time,'_pheatmap_topDEgene.png'), units='in', width=10, height=10, res=320)
      pheatmap(sig_counts, 
               main = paste0(time,' Overall gene-sample correlations'),
               color = heat_colors, 
               cluster_rows = TRUE, 
               show_rownames = FALSE,
               annotation = anno,
               border_color = NA, 
               fontsize = 10, 
               scale = "row", 
               fontsize_row = 10, 
               height = 20)
      if (dev.cur() > 1) {dev.off()}
    } else {
      message("Skipping pheatmap: fewer than 2 genes found (nrow < 2).")
    }
    
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
                    title=paste0(time,' overall DE Genes'),
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
    ggsave(paste0(stats_dir,'figures/',"Overall",time,'_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
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
  
  ## create pseudobulk samples by celltype, sample, and condition
  pseudobulk = AggregateExpression(
    seurat,
    return.seurat = TRUE,
    assays = "RNA",
    group.by = c("celltype", "sample", "condition")
  )
  
  # Number of cells by sample and celltype
  n_cells <- seurat@meta.data %>% 
    dplyr::count(sample, celltype) %>% 
    rename("n"="n_cells")
  n_cells$sample <- str_replace(n_cells$sample, "_", "-")
  
  meta_pseudobulk <- left_join(pseudobulk@meta.data, n_cells)
  rownames(meta_pseudobulk) <- meta_pseudobulk$orig.ident
  pseudobulk@meta.data <- meta_pseudobulk
  
  # Turn condition into a factor
  pseudobulk$condition <- factor(pseudobulk$condition, levels=c("ctrl","pb"))
  
  # pseudobulk@meta.data %>% head()
  
  #### DE analysis by cell type ----
  get_dds_resultsAvsB <- function(clustx, A, B, padj_cutoff = 0.05) {
    cat('Beginning cluster ', clustx, 'for timepoint:',time,'\n\n')
    
    # Extract counts for cluster x
    bulk_clustx = subset(pseudobulk, subset= (celltype == clustx) & (condition %in% c("ctrl", "pb")))
    
    ggplot(bulk_clustx@meta.data, aes(x=sample, y=n_cells, fill=condition)) +
      geom_bar(stat="identity", color="black") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(title=paste0(time,' ',clustx,' # of cells per sample'), x="Sample name", y="Number of cells") +
      geom_text(aes(label=n_cells), vjust=-0.5)
    ggsave(paste0(stats_dir,"figures/",clustx,"_",time,"_ncells_per_sample.png"),dpi=320,height=10,width=10)
    
    # Check if we have samples from both conditions
    condition_counts <- table(bulk_clustx@meta.data$condition)
    if (length(condition_counts) < 2 || any(condition_counts == 0)) {
      cat('Skipping cluster', clustx, 'because it does not have samples from both conditions\n')
      cat('Condition counts:', paste(names(condition_counts), condition_counts, sep="=", collapse=", "), '\n')
      return(NULL)
    }
    
    cluster_counts <- FetchData(bulk_clustx, layer="counts", vars=rownames(bulk_clustx))
    
    # zero-count proportions?
    cat('\n;',time,clustx,'zero-count proportion:\n')
    zeroprop = (sum(rowSums(cluster_counts == 0))) / (nrow(cluster_counts) * ncol(cluster_counts))
    print(zeroprop)
    zero.count.props[[paste0(time,'_',clustx)]] = zeroprop
    
    
    dds <- DESeqDataSetFromMatrix(t(cluster_counts), 
                                  colData = bulk_clustx@meta.data, 
                                  design = ~ condition)
    
    # prefiltering low count genes. 
    # Want genes that have at least 5 counts in at least 4 samples. 
    # following recommendation from
    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
    # smallestGroupSize = min(table(colData(dds)$condition))
    # smallestGroupSize = 4
    # # cat('\nSize of smaller group is',smallestGroupSize)
    # cat('\n# of samples per condition for',time,clustx,'\n')
    # print(table(colData(dds)$condition))
    # 
    # keep = rowSums(counts(dds) >= 5) >= smallestGroupSize
    # 
    # cat('\n# of genes filtered and kept:\n')
    # print(table(keep))
    # 
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
    tryCatch({
      pca_data_condition = DESeq2::plotPCA(rld, intgroup = "condition", returnData= TRUE) %>%
        mutate(name_parsed = gsub(paste0(clustx,'_'),'',name))
      ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = condition, label = name_parsed)) +
        geom_point() + 
        geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
        theme_classic2() +
        ggtitle(paste0(time,'_',clustx)) +
        xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 
      ggsave(paste0(stats_dir,"figures/",clustx,'_',time,"_specific_PCAplot.png"),dpi=320,height=10,width=10)
    }, error = function(e) {
      cat("Error in plotPCA across condition:", conditionMessage(e), "\n")
    })
    
    tryCatch({
      pca_data_n_cells = DESeq2::plotPCA(rld, intgroup = "n_cells", returnData=TRUE) %>%
        mutate(name_parsed = gsub(paste0(clustx,'_'),'',name))
      ggplot(pca_data_n_cells, aes(x = PC1, y = PC2, color = n_cells, label = name_parsed)) +
        geom_point() + 
        geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
        theme_classic2() +
        ggtitle(paste0(time,'_',clustx)) +
        xlab(paste0("PC1: ", round(attr(pca_data_n_cells, "percentVar")[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(attr(pca_data_n_cells, "percentVar")[2] * 100), "% variance"))
      ggsave(paste0(stats_dir,"figures/",clustx,'_',time,"_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
    }, error = function(e) {
      cat("Error in plotPCA across condition, including cell count:", conditionMessage(e), "\n")
    })
    
    ## Extract rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    rename_samples <- bulk_clustx$sample
    colnames(rld_cor) <- rename_samples
    rownames(rld_cor) <- rename_samples
    
    anno <- bulk_clustx@meta.data %>%
      select(sample, condition) %>% 
      remove_rownames() %>% 
      column_to_rownames("sample")
    
    ## Plot and save heatmap
    tryCatch({
      png(paste0(stats_dir,"figures/", clustx,'_',time, "_specific_heatmap.png"),
          height = 6, width = 7.5, units = "in", res = 300)
      pheatmap(rld_cor, annotation_row=anno, annotation_col=anno,
               main = paste0(time,' ',clustx,' sample correlations'))
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
    # bp = get_bioc_param(workers=nthreads)
    dds <- DESeq(dds,
                 test='LRT',
                 fitType='parametric', # the default
                 reduced= ~1,
                 minReplicatesForReplace=Inf,
                 useT=TRUE,
                 minmu=1e-6)
                 # parallel = TRUE,
                 # BPPARAM = bp)
    
    cat('Finished fitting DESeq2 model for',time,clustx,'\n')
    
    ## Plot dispersion estimates
    tryCatch({
      png(paste0(stats_dir,"figures/", clustx,'_',time, "_dispersion_plot.png"),
          height = 5, width = 6, units = "in", res = 300)
      plotDispEsts(dds)
      if (dev.cur() > 1) {dev.off()}
    }, error = function(e) {
      cat("Error in creating dispersion plot:", conditionMessage(e), "\n")
      if (dev.cur() > 1) dev.off()
    })
    
    
    ## Output and shrink results of DE test for contrast A vs B
    contrast <- paste(c("condition", A, "vs", B), collapse = "_")
    print(resultsNames(dds))
    
    res <- results(dds, name = contrast, alpha = 0.05, independentFiltering = TRUE)
    # res <- results(dds, name = 'conditionpb', alpha = 0.05)
    res <- lfcShrink(dds, 
                     coef = contrast, 
                     res = res, 
                     type = "apeglm")
                     # parallel = TRUE,
                     # BPPARAM = bp)
    
    # res <- lfcShrink(dds, coef = 'conditionpb', res = res, type = "apeglm")
    
    ## Turn the results object into a tibble for use with tidyverse functions
    res_tbl <- res %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble()
    
    write.csv(res_tbl,
              paste0(stats_dir,"results/", clustx, "_", contrast,'_',time, "_all_genes_deseq2.csv"),
              quote = FALSE, 
              row.names = FALSE)
    
    ## Subset the significant results
    sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
      dplyr::arrange(padj)
    
    write.csv(sig_res,
              paste0(stats_dir,"results/", clustx, "_", contrast,'_',time, "_signif_genes_deseq2.csv"),
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
      
      # Check that gene column exists
      if (!"gene" %in% colnames(top20_sig_df)) {
        stop("Gene column not found in data frame")
      }
      
      # Convert to data.table and melt into long format
      top20_sig_df <- reshape2::melt(top20_sig_df, 
                                     id.vars = c("gene"),
                                     variable.name = "sample",
                                     value.name='value')
      
      ## Foramat sample name to be same as in dds object
      # top20_sig_df$sample <- grep(paste0(time, "-(ctrl|pb)\\d+"), top20_sig_df$sample)
      
      # top20_sig_df$sample <- gsub(paste0(time,"\\."),paste0(time,'-'), top20_sig_df$sample)
      
      
      top20_sig_df$sample <- gsub("\\.", "-", 
                                  regmatches(top20_sig_df$sample, 
                                             regexpr(paste0(time, "\\.(ctrl|pb)\\d+"), 
                                                     top20_sig_df$sample)))
      
      
      # Proceed with joining the metadata and plotting
      top20_sig_df <- merge(top20_sig_df, as.data.frame(colData(dds)) %>% rownames_to_column('rownames'),
                            by.x='sample',
                            by.y='sample',
                            al.x=TRUE)
      
      ggplot(top20_sig_df, aes(y = value, x = condition, col = condition)) +
        geom_boxplot(outlier.shape=NA, width=0.2)+
        geom_jitter(height = 0, width = 0.5)+
        scale_y_continuous(trans = 'log10') +
        ylab("log10 of normalized expression level") +
        xlab("condition") +
        ggtitle(paste0("Top ", n, " Significant DE Genes")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        facet_wrap(~ gene)
      ggsave(paste0(stats_dir,"figures/", clustx, "_", contrast,'_',time, "_top_DE_genes.png"))
      
      ## Extract normalized counts for significant genes only
      sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
      
      if (is.matrix(sig_counts) || is.data.frame(sig_counts)) {
        colnames(sig_counts) = bulk_clustx$sample
      }      
      ## generate pheatmap
      if ((nrow(sig_res) >= 2) & (!is.null(nrow(sig_res)))) {
        png(filename=paste0(stats_dir,'figures/',clustx,'_',contrast,'_',time,'_pheatmap_topDEgene.png'), units='in', width=10, height=10, res=320)
        pheatmap(sig_counts, 
                 main = paste0(time,' ',clustx,' gene-sample correlations'),
                 color = heat_colors, 
                 cluster_rows = TRUE, 
                 show_rownames = FALSE,
                 annotation_col = anno,
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
                      title=paste0(time,'_',clustx,' DE Genes'),
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
      ggsave(paste0(stats_dir,'figures/',clustx,'_',contrast,'_',time,'_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
    }, error = function(e) {
      cat("Error in creating volcano plot:", conditionMessage(e), "\n")
    })
    
    cat('\nFinished analyzing',time,clustx,'\n')
  }
  
  # Run the script on all clusters comparing stimulated condition relative to control condition
  # map(cluster_names, get_dds_resultsAvsB, A = "pb", B = "ctrl", padj_cutoff = 0.05)
  
  map(cluster_names, function(clustx) {
    tryCatch({
      get_dds_resultsAvsB(clustx, A = "pb", B = "ctrl", padj_cutoff = 0.05)
    }, error = function(e) {
      cat("Error processing",time, clustx, ":", conditionMessage(e), "\n")
      return(NULL)
    })
  })
}


saveRDS(zero.count.props,paste0(stats_dir,'zero_count_propotions.rds'))










# Timecourse analysis ----
cat('\nBeginning timecourse analysis\n')

# Create directories for timecourse results if they don't exist
tc_figs_dir <- paste0(stats_dir, "timecourse/figures/")
tc_results_dir <- paste0(stats_dir, "timecourse/results/")

dirs_tc <- c(tc_figs_dir, tc_results_dir)
for (dir in dirs_tc) {
  if (!dir.exists(dir)) { 
    dir.create(dir, recursive=TRUE) 
  }
}

# Ensure timepoint is a factor with correct ordering
de_seurat$timepoint <- factor(de_seurat$timepoint, levels = timepoints,ordered=TRUE)
de_seurat$timepoint <- relevel(de_seurat$timepoint, ref = "wk3")

# Create pseudobulk samples across all timepoints
tc_bulk <- AggregateExpression(
  de_seurat,
  return.seurat = TRUE,
  assays = "RNA",
  group.by = c("sample", "condition", "timepoint")
)

# Get cell counts
n_cells_tc <- de_seurat@meta.data %>% 
  dplyr::count(sample, timepoint) %>% 
  rename('n' = 'n_cells') 

# Update metadata
meta_tc_bulk <- left_join(tc_bulk@meta.data, n_cells_tc, by = c("sample", "timepoint"))
rownames(meta_tc_bulk) <- rownames(tc_bulk@meta.data)
tc_bulk@meta.data <- meta_tc_bulk

# Ensure condition is a factor
tc_bulk$condition <- factor(tc_bulk$condition, levels = c("ctrl", "pb"))

# Get counts for DESeq2
tc_counts <- FetchData(tc_bulk, layer = "counts", vars = rownames(tc_bulk))

# IMPORTANT: Set the factor levels in tc_bulk AFTER aggregation
tc_bulk$timepoint <- factor(tc_bulk$timepoint, levels = c('wk3','wk10','mn7','mn18'))
tc_bulk$timepoint <- relevel(tc_bulk$timepoint, ref = "wk3")

# Create DESeq2 object
dds_tc <- DESeqDataSetFromMatrix(t(tc_counts), 
                                 colData = tc_bulk@meta.data, 
                                 design = ~ condition + timepoint + condition:timepoint)

# Compute size factors 
scr_tc = scran::computeSumFactors(dds_tc)
sizeFactors(dds_tc) = sizeFactors(scr_tc)

# Run DESeq2 with likelihood ratio test
dds_tc <- DESeq(dds_tc, 
                test = "LRT", 
                reduced = ~ condition + timepoint,
                fitType = 'parametric',
                minReplicatesForReplace = Inf,
                useT = TRUE,
                minmu = 1e-6)

# Get results (genes with condition-specific time effects)
res_tc <- results(dds_tc, alpha = 0.05)

# Transform to data frame
res_tc_tbl <- res_tc %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(padj)

# Write results
write.csv(res_tc_tbl,
          paste0(tc_results_dir, "overall_condition_specific_time_effects_all_genes.csv"),
          quote = FALSE,
          row.names = FALSE)

# Filter for significant genes
sig_res_tc <- dplyr::filter(res_tc_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Write significant results
write.csv(sig_res_tc,
          paste0(tc_results_dir, "overall_condition_specific_time_effects_signif_genes.csv"),
          quote = FALSE,
          row.names = FALSE)

cat("Number of significant genes with condition-specific time effects:", nrow(sig_res_tc), "\n")

# Plot top genes if any significant results
if (nrow(sig_res_tc) > 0) {
  # Extract coefficients - these show the interaction terms
  betas <- coef(dds_tc)
  
  # Display coefficient names
  cat('\nCoefficient names:\n')
  print(colnames(betas))
  
  # Extract interaction terms
  interaction_cols <- grep("conditionpb.*timepoint", colnames(betas), value = TRUE)
  
  # Plot heatmap of interaction terms for top genes
  if (length(interaction_cols) > 0) {
    top_genes <- head(order(res_tc$padj), min(20, nrow(sig_res_tc)))
    mat <- betas[top_genes, interaction_cols, drop = FALSE]
    
    # Add gene names
    rownames(mat) <- res_tc_tbl$gene[top_genes]
    
    # Threshold values for better visualization
    thr <- 3
    mat[mat < -thr] <- -thr
    mat[mat > thr] <- thr
    
    # Create heatmap
    tryCatch({
      png(paste0(tc_figs_dir, "top_DE_genes_timecourse_heatmap.png"),
          res = 320, height = 10, width = 10, units = "in")
      pheatmap(mat, 
               breaks = seq(from = -thr, to = thr, length = 101),
               cluster_col = FALSE,
               main = "Condition-specific time effects")
      if (dev.cur() > 1) {dev.off()}
    }, error = function(e) {
      cat("Error in creating heatmap:", conditionMessage(e), "\n")
      if (dev.cur() > 1) dev.off()
    })
  }
  
  # Plot the top gene across time points
  top_gene <- sig_res_tc$gene[1]
  
  counts_data <- plotCounts(dds_tc, which.min(res_tc$padj),
                            intgroup = c("timepoint", "condition"), 
                            returnData = TRUE)
  
  # Plot
  tryCatch({
    # Extract numeric value from timepoint for plotting
    # This converts 'wk3', 'wk10', etc. to numerical values for proper ordering
    counts_data$time_order <- match(counts_data$timepoint, timepoints)
    
    ggplot(counts_data, aes(x = time_order, y = count, color = condition, group = condition)) +
      geom_point() +
      stat_summary(fun = mean, geom = "line") +
      scale_y_log10() +
      labs(title = paste("Expression of", top_gene, "over time"),
           x = "Time",
           y = "Normalized counts (log10)") +
      scale_x_continuous(breaks = 1:length(timepoints),
                         labels = timepoints) +
      theme_classic2()
    ggsave(paste0(tc_figs_dir, "top_gene_timecourse_plot.png"), dpi = 320, height = 8, width = 10)
  }, error = function(e) {
    cat("Error in creating top gene plot:", conditionMessage(e), "\n")
  })
  
  # Create a separate cell type analysis function if needed
  # ... (cell type specific code would go here)
}

cat('\nTimecourse analysis completed\n')







# sessionInfo() ----
sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: x86_64-pc-linux-gnu
# Running under: Red Hat Enterprise Linux 8.8 (Ootpa)
# 
# Matrix products: default
# BLAS:   /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRblas.so
# LAPACK: /sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.0/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
# 
# time zone: America/Detroit
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base
# 
# other attached packages:
#   [1] EnhancedVolcano_1.22.0      ggrepel_0.9.6               RColorBrewer_1.1-3
# [4] DESeq2_1.46.0               png_0.1-8                   apeglm_1.26.1
# [7] SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0 Biobase_2.64.0
# [10] GenomicRanges_1.56.2        GenomeInfoDb_1.40.1         IRanges_2.38.1
# [13] MatrixGenerics_1.16.0       matrixStats_1.5.0           S4Vectors_0.42.1
# [16] BiocGenerics_0.50.0         reshape2_1.4.4              edgeR_4.2.2
# [19] limma_3.60.6                Matrix.utils_0.9.8          Matrix_1.7-0
# [22] openai_0.4.1                GPTCelltype_1.0.1           UCell_2.8.0
# [25] qs2_0.1.5                   presto_1.0.0                data.table_1.17.0
# [28] pheatmap_1.0.12             Polychrome_1.5.1            glmGamPoi_1.16.0
# [31] scales_1.3.0                viridis_0.6.5               viridisLite_0.4.2
# [34] plyr_1.8.9                  ggpubr_0.6.0                ggthemes_5.1.0
# [37] lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1
# [40] dplyr_1.1.4                 purrr_1.0.4                 readr_2.1.5
# [43] tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.1
# [46] tidyverse_2.0.0             cowplot_1.1.3               patchwork_1.3.0
# [49] harmony_1.2.3               Rcpp_1.0.14                 sctransform_0.4.1
# [52] Seurat_5.2.1                SeuratObject_5.0.2          sp_2.2-0
# 
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0     httr_1.4.7                numDeriv_2016.8-1.1
# [4] tools_4.4.0               backports_1.5.0           utf8_1.2.4
# [7] R6_2.6.1                  lazyeval_0.2.2            uwot_0.2.3
# [10] sn_2.1.1                  withr_3.0.2               gridExtra_2.3
# [13] progressr_0.15.1          cli_3.6.4                 textshaping_1.0.0
# [16] spatstat.explore_3.3-4    fastDummies_1.7.5         sandwich_3.1-1
# [19] labeling_0.4.3            mvtnorm_1.3-3             spatstat.data_3.1-4
# [22] ggridges_0.5.6            pbapply_1.7-2             systemfonts_1.2.1
# [25] parallelly_1.42.0         plotrix_3.8-4             bbmle_1.0.25.1
# [28] rstudioapi_0.17.1         generics_0.1.3            ica_1.0-3
# [31] spatstat.random_3.3-2     car_3.1-3                 MetBrewer_0.2.0
# [34] abind_1.4-8               lifecycle_1.0.4           multcomp_1.4-28
# [37] scatterplot3d_0.3-44      carData_3.0-5             mathjaxr_1.6-0
# [40] SparseArray_1.4.8         Rtsne_0.17                grid_4.4.0
# [43] promises_1.3.2            dqrng_0.4.1               bdsmatrix_1.3-7
# [46] crayon_1.5.3              miniUI_0.1.1.1            lattice_0.22-6
# [49] beachmat_2.20.0           pillar_1.10.1             metapod_1.12.0
# [52] future.apply_1.11.3       codetools_0.2-20          mutoss_0.1-13
# [55] glue_1.8.0                leidenbase_0.1.32         spatstat.univar_3.1-2
# [58] vctrs_0.6.5               spam_2.11-1               Rdpack_2.6.2
# [61] gtable_0.3.6              emdbook_1.3.13            rbibutils_2.3
# [64] S4Arrays_1.4.1            mime_0.13                 coda_0.19-4.1
# [67] survival_3.5-8            statmod_1.5.0             bluster_1.14.0
# [70] TH.data_1.1-3             fitdistrplus_1.2-2        ROCR_1.0-11
# [73] nlme_3.1-164              RcppAnnoy_0.0.22          irlba_2.3.5.1
# [76] KernSmooth_2.23-22        colorspace_2.1-1          mnormt_2.1.1
# [79] tidyselect_1.2.1          compiler_4.4.0            BiocNeighbors_1.22.0
# [82] TFisher_0.2.0             DelayedArray_0.30.1       plotly_4.10.4
# [85] stringfish_0.16.0         lmtest_0.9-40             digest_0.6.37
# [88] goftest_1.2-3             spatstat.utils_3.1-2      RhpcBLASctl_0.23-42
# [91] XVector_0.44.0            htmltools_0.5.8.1         pkgconfig_2.0.3
# [94] sparseMatrixStats_1.16.0  fastmap_1.2.0             rlang_1.1.5
# [97] htmlwidgets_1.6.4         UCSC.utils_1.0.0          shiny_1.10.0
# [100] DelayedMatrixStats_1.26.0 farver_2.1.2              zoo_1.8-13
# [103] jsonlite_1.9.1            BiocParallel_1.38.0       BiocSingular_1.20.0
# [106] magrittr_2.0.3            Formula_1.2-5             scuttle_1.14.0
# [109] GenomeInfoDbData_1.2.12   dotCall64_1.2             munsell_0.5.1
# [112] reticulate_1.41.0.1       stringi_1.8.7             zlibbioc_1.50.0
# [115] MASS_7.3-60.2             parallel_4.4.0            listenv_0.9.1
# [118] deldir_2.0-4              splines_4.4.0             multtest_2.60.0
# [121] tensor_1.5                hms_1.1.3                 qqconf_1.3.2
# [124] locfit_1.5-9.12           igraph_2.1.4              spatstat.geom_3.3-5
# [127] ggsignif_0.6.4            RcppHNSW_0.6.0            ScaledMatrix_1.12.0
# [130] metap_1.12                RcppParallel_5.1.10       scran_1.32.0
# [133] tzdb_0.4.0                httpuv_1.6.15             grr_0.9.5
# [136] RANN_2.6.2                polyclip_1.10-7           future_1.34.0
# [139] scattermore_1.2           rsvd_1.0.5                broom_1.0.5
# [142] xtable_1.8-4              RSpectra_0.16-2           rstatix_0.7.2
# [145] later_1.4.1               ragg_1.3.3                cluster_2.1.6
# [148] timechange_0.3.0          globals_0.16.3



