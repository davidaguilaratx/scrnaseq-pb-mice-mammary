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
# library(BiocParallel)

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
  test = subset(epithelial, subset=(timepoint == time))
  
  print(table(test$sample, test$subtype))
}

# Set thresholds
padj_cutoff = 0.05
log2fc_cutoff = 0.58 # 50% change


zero.count.props = list()
# for (time in timepoints[2:4]) {
for (time in timepoints) {
  
  
  if (!dir.exists(paste0(stats_dir,'epithelial/figures'))) { dir.create(paste0(stats_dir,'epithelial/figures'),
                                                              recursive=TRUE) }
  if (!dir.exists(paste0(stats_dir,'epithelial/results'))) { dir.create(paste0(stats_dir,'epithelial/results'),
                                                              recursive=TRUE) }
  
  seurat = subset(epithelial, subset=(timepoint == time))
  
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
  ggsave(paste0(stats_dir,"epithelial/figures/","Overall",time,"_ncells_per_sample.png"),dpi=320,height=10,width=10)
  
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
    ggsave(paste0(stats_dir,"epithelial/figures/","Overall",time,"_specific_PCAplot.png"),dpi=320,height=10,width=10)
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
    ggsave(paste0(stats_dir,"epithelial/figures/","Overall",time,"_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
  }, error = function(e) {
    cat("Error in plotPCA across condition, including cell count:", conditionMessage(e), "\n")
  })
  

  # DESeq2::plotPCA(rld_bulk, ntop = 500, intgroup = "cell_count")
  # ggsave(paste0(stats_dir,"epithelial/figures/",time,"_overall_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
  
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
    png(paste0(stats_dir,'epithelial/figures/',"Overall",time,'_specific_heatmap.png'),
        res=320, height = 10, width = 10, units = "in")
    pheatmap(rld_cor_bulk, annotation_col=anno, annotation_row=anno)
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
    png(paste0(stats_dir,'epithelial/figures/',"Overall",time,'_dispersion_plot.png'),
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
            paste0(stats_dir,"epithelial/results/overall_condition_", 
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
            paste0(stats_dir,"epithelial/results/overall_condition_", 
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
    ggsave(paste0(stats_dir,"/epithelial/figures/", "Overall_pb_vs_ctrl_",time,"_top_DE_genes.png"))
    
    ## Extract normalized counts for significant genes only
    sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res_bulk$gene, ]
    
    if (is.matrix(sig_counts) || is.data.frame(sig_counts)) {
      colnames(sig_counts) = bulk$sample
    }
    
    ## generate pheatmap
    if ((nrow(sig_res_bulk) >= 2) & (!is.null(nrow(sig_res_bulk)))) {
      png(filename=paste0(stats_dir,'/epithelial/figures/',"Overall_pb_vs_ctrl_",time,'_pheatmap_topDEgene.png'), units='in', width=10, height=10, res=320)
      pheatmap(sig_counts, 
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
  ggsave(paste0(stats_dir,'epithelial/figures/',"Overall",time,'_DEgenes_evolcano.png'),dpi=320,width=15,height=10)
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
  
  ## create pseudobulk csamples by celltype, sample, and condition
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
    ggsave(paste0(stats_dir,"epithelial/figures/",clustx,"_",time,"_ncells_per_sample.png"),dpi=320,height=10,width=10)
    
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
    smallestGroupSize = 4
    # cat('\nSize of smaller group is',smallestGroupSize)
    cat('\n# of samples per condition for',time,clustx,'\n')
    print(table(colData(dds)$condition))
    
    keep = rowSums(counts(dds) >= 5) >= smallestGroupSize
    
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
      pca_data_condition = DESeq2::plotPCA(rld, intgroup = "condition", returnData= TRUE) %>%
        mutate(name_parsed = gsub(paste0(clustx,'_'),'',name))
      ggplot(pca_data_condition, aes(x = PC1, y = PC2, color = condition, label = name_parsed)) +
        geom_point() + 
        geom_text_repel(vjust = 1.5, hjust = 0.5, show.legend = FALSE) +
        theme_classic2() +
        ggtitle(paste0(time,'_',clustx)) +
        xlab(paste0("PC1: ", round(attr(pca_data_condition, "percentVar")[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(attr(pca_data_condition, "percentVar")[2] * 100), "% variance")) 
      ggsave(paste0(stats_dir,"epithelial/figures/",clustx,'_',time,"_specific_PCAplot.png"),dpi=320,height=10,width=10)
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
      ggsave(paste0(stats_dir,"epithelial/figures/",clustx,'_',time,"_specific_PCAplot_with_cellcounts.png"),dpi=320,height=10,width=10)
    }, error = function(e) {
      cat("Error in plotPCA across condition, including cell count:", conditionMessage(e), "\n")
    })
    
    ## Extract rlog matrix from the object and compute pairwise correlation values
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    
    rename_samples <- bulk_clustx$sample
    colnames(rld_cor) <- str_replace_all(colnames(rld_cor), rename_samples)
    rownames(rld_cor) <- str_replace_all(rownames(rld_cor), rename_samples)
    
    anno <- bulk_clustx@meta.data %>%
      select(sample, condition) %>% 
      remove_rownames() %>% 
      column_to_rownames("sample")
    
    ## Plot and save heatmap
    tryCatch({
      png(paste0(stats_dir,"/epithelial/figures/", clustx,'_',time, "_specific_heatmap.png"),
          height = 6, width = 7.5, units = "in", res = 300)
      pheatmap(rld_cor, annotation_row=anno, annotation_col=anno)
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
    bp = get_bioc_param(workers=nthreads)
    dds <- DESeq(dds,
                 test='LRT',
                 fitType='parametric', # the default
                 reduced= ~1,
                 minReplicatesForReplace=Inf,
                 useT=TRUE,
                 minmu=1e-6,
                 parallel = TRUE,
                 BPPARAM = bp)

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
    
    
    ## Output and shrink results of DE test for contrast A vs B
    contrast <- paste(c("condition", A, "vs", B), collapse = "_")
    print(resultsNames(dds))
    
    res <- results(dds, name = contrast, alpha = 0.05, independentFiltering = TRUE)
    # res <- results(dds, name = 'conditionpb', alpha = 0.05)
    res <- lfcShrink(dds, 
                     coef = contrast, 
                     res = res, 
                     type = "apeglm",
                     parallel = TRUE,
                     BPPARAM = bp)
    
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
      top20_sig_df <- merge(top20_sig_df, as.data.frame(colData(dds)) %>% rownames_to_column('rownames'),
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
      ggsave(paste0(stats_dir,"/epithelial/figures/", clustx, "_", contrast,'_',time, "_top_DE_genes.png"))
      
      ## Extract normalized counts for significant genes only
      sig_counts <- normalized_counts[rownames(normalized_counts) %in% sig_res$gene, ]
      
      if (is.matrix(sig_counts) || is.data.frame(sig_counts)) {
        colnames(sig_counts) = bulk_clustx$sample
      }      
      ## generate pheatmap
      if ((nrow(sig_res) >= 2) & (!is.null(nrow(sig_res)))) {
        png(filename=paste0(stats_dir,'/epithelial/figures/',clustx,'_',contrast,'_',time,'_pheatmap_topDEgene.png'), units='in', width=10, height=10, res=320)
        pheatmap(sig_counts, 
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




