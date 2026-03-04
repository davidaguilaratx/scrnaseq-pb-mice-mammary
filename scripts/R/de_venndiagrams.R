library(VennDiagram)
library(grid)

stats_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/"

### grab celltypes from results ----
### directory of differential expression results 
de_directory = paste0(stats_dir,'results/')

### grab de result file names and extract cell types 
de_results = list.files(path=de_directory, pattern='all_genes')

### List of cell types (replace with your actual cell types)
celltypes = sub(pattern='_condition.*',replacement='',x=de_results)
celltypes


# Function to create directory configuration
create_config <- function(dirs, labels) {
  if (length(dirs) < 2 || length(dirs) > 4) {
    stop("Please provide between 2 and 4 directories")
  }
  if (length(dirs) != length(labels)) {
    stop("Number of directories must match number of labels")
  }
  
  return(list(
    directories = dirs,
    labels = labels
  ))
}


# Add this to the read_gene_list function, right after the for loop starts:
read_gene_list <- function(cell_type, config, return_logFC=FALSE) {
  # Initialize lists to store data
  gene_lists <- list()
  all_data <- list()
  
  for(i in seq_along(config$directories)) {
    dir <- config$directories[i]
    label <- config$labels[i]
    
    if(return_logFC) {
      # Read full data for log fold changes
      file <- paste0(dir, cell_type, "_condition_pb_vs_ctrl_all_genes_deseq2.csv")
      print(paste("Attempting to read file:", file))  # Debug line
      if(!file.exists(file)) {
        warning(paste("File not found:", file))
        next
      }
      data <- read.csv(file)[, c("gene", "log2FoldChange")]
      colnames(data)[2] <- paste0("logFC_", label)
      all_data[[i]] <- data
    } else {
      # Read only significant genes
      file <- paste0(dir, cell_type, "_condition_pb_vs_ctrl_signif_genes_deseq2.csv")
      print(paste("Attempting to read file:", file))  # Debug line
      if(!file.exists(file)) {
        warning(paste("File not found:", file))
        next
      }
      gene_lists[[label]] <- read.csv(file)$gene
      print(paste("Found", length(gene_lists[[label]]), "genes for", label))  # Debug line
    }
  }
  
  if(return_logFC) {
    # Merge all data frames
    merged_data <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), all_data)
    return(merged_data)
  } else {
    return(gene_lists)
  }
}

# Helper function to format intersection labels with gene names
# Helper function to format intersection labels with gene names
format_venn_labels <- function(gene_lists) {
  # Use VennDiagram's calculate.overlap function
  intersections <- calculate.overlap(gene_lists)
  
  # Format each intersection with count and gene names
  labels <- sapply(intersections, function(genes) {
    if(length(genes) > 0) {
      paste0(length(genes), "\n", paste(genes, collapse="\n"))
    } else {
      "0"
    }
  })
  
  return(labels)
}

# make the venndiagram
create_venn_diagram <- function(gene_lists, cell_type, config, output_fig_dir, colors) {
  n_sets <- length(gene_lists)
  
  if (n_sets < 2 || n_sets > 4) {
    warning(paste("Cannot create Venn diagram for", cell_type, "- need 2-4 sets, got", n_sets))
    return(NULL)
  }
  
  # Format labels with gene names
  labels <- format_venn_labels(gene_lists)
  
  # Set up basic parameters
  venn_params <- list(
    x = gene_lists,
    filename = paste0(output_fig_dir, 'venn_diagrams/', cell_type, '_venn_diagram_signif_genes.png'),
    category.names = config$labels[1:n_sets],  # Only use the labels we need
    fill = colors[1:n_sets],
    lwd = 2,
    lty = rep('blank', n_sets),
    fontface = "bold",
    fontfamily = "sans",
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    label.col = 'black',
    cat.col = 'black',
    cat.cex = 1,
    cex = 1,
    euler.d = TRUE,
    scaled = TRUE,
    ext.text = TRUE,
    ext.percent = 0.15,
    margin = 0.1,
    cat.default.pos = "outer",
    output = TRUE,
    imagetype = "png",
    height = 3000,
    width = 3000,
    resolution = 300,
    compression = "lzw",
    sep.dist = 0.15,
    rotation.degree = 0,
    disable.logging = TRUE
  )
  
  # Set category positions based on number of sets
  if(n_sets == 2) {
    venn_params$cat.pos <- c(-50, 50)
    venn_params$cat.dist <- c(0.05, 0.05)
  } else if(n_sets == 3) {
    venn_params$cat.pos <- c(-40, 40, 180)
    venn_params$cat.dist <- c(0.1, 0.1, 0.1)
  } else if(n_sets == 4) {
    venn_params$cat.pos <- c(-27, 27, 135, 135)
    venn_params$cat.dist <- c(0.055, 0.055, 0.085, 0.085)
  }
  
  # Create the Venn diagram
  vennplot <- do.call(venn.diagram, venn_params)
  
  return(vennplot)
}

# Modify the run_venn_analysis function:
run_venn_analysis <- function(config, output_fig_dir, colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442")) {
  # Validate output directory
  venn_dir <- file.path(output_fig_dir, "venn_diagrams")
  dir.create(venn_dir, showWarnings = FALSE, recursive = TRUE)
  print(paste("Output directory:", venn_dir))
  
  # Get cell types from one of the directories and filter properly
  de_results <- list.files(path = config$directories[1], pattern = 'all_genes')
  # Filter out any files that don't end with the expected pattern
  de_results <- de_results[grep("_condition_pb_vs_ctrl_all_genes_deseq2.csv$", de_results)]
  print("Found files:")
  print(de_results)
  
  celltypes <- sub(pattern = '_condition.*', replacement = '', x = de_results)
  print("Cell types found:")
  print(celltypes)
  
  # Process each cell type
  signif_genes_by_celltype <- list()
  
  for (celltype in celltypes) {
    print(paste("\nProcessing cell type:", celltype))
    # Read gene lists
    gene_lists <- read_gene_list(celltype, config)
    if (length(gene_lists) > 0) {  # Only process if we have gene lists
      signif_genes_by_celltype[[celltype]] <- unique(unlist(gene_lists))
      # Create Venn diagram
      create_venn_diagram(gene_lists, celltype, config, output_fig_dir, colors)
    }
  }
  
  # Return all significant genes
  return(unique(unlist(signif_genes_by_celltype)))
}


# Example usage:
# config <- create_config(
#   dirs = c(
#     "path/to/first/directory/",
#     "path/to/second/directory/",
#     "path/to/third/directory/",
#     "path/to/fourth/directory/"
#   ),
#   labels = c(
#     "CB_LRT_unfiltered",
#     "LRT_filtered",
#     "CB_LRT_filtered",
#     "LRT_unfiltered"
#   )
# )
# 
# output_fig_dir <- "path/to/output/directory/"
# all_signif_genes <- run_venn_analysis(config, output_fig_dir)



cr_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/'
deseq2rr_dir <- paste0(cr_deseq2_dir,'stats_rerun/DESeq2/results/') # filtered wald
deseq2LRT_dir <- paste0(cr_deseq2_dir,'stats/DESeq2/results/') # filtered LRT
deseq2_dir <- paste0(cr_deseq2_dir,'stats/DESeq2_unfiltered/results/') # unfiltered Wald
deseq2uLRT_dir <- paste0(cr_deseq2_dir,'stats/DESeq2_unfiltered_LRT/results/') # unfiltered LRT
cb_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2/results/' # cellbender unfiltered LRT
cb_deseq2_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/DESeq2_gene_fltrd/results/'

cb_deseq2_MAD_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0_MAD/stats/DESeq2_gene_fltrd/results/'


crv9_deseq2_fltrd_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2_gene_fltrd/results/'
crv9_deseq2_dir = 'C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/DESeq2/results/'

# config <- create_config(
#   dirs = c(
#     crv9_deseq2_dir,
#     deseq2LRT_dir,
#     crv9_deseq2_fltrd_dir,
#     deseq2uLRT_dir
#   ),
#   labels = c(
#     "CRV9_LRT_unfiltered",
#     "LRT_filtered",
#     "CRv9_LRT_filtered",
#     "LRT_unfiltered"
#   )
# )

# config <- create_config(
#   dirs = c(
#     deseq2LRT_dir,
#     cb_deseq2_fltrd_dir,
#     crv9_deseq2_fltrd_dir
#   ),
#   labels = c(
#     "LRT_filtered",
#     "CRv9_LRT_filtered",
#     "CB_LRT_filtered"
#   )
# )


# config <- create_config(
#   dirs = c(
#     cb_deseq2_dir,
#     crv9_deseq2_dir,
#     cb_deseq2_fltrd_dir,
#     crv9_deseq2_fltrd_dir
#   ),
#   labels = c(
#     "CB_LRT_unfiltered",
#     "crv9_LRT_fltrd",
#     "CB_LRT_filtered",
#     "crv9_LRT_unfiltered"
#   )
# )

# config <- create_config(
#   dirs = c(
#     cb_deseq2_dir,
#     deseq2LRT_dir,
#     cb_deseq2_fltrd_dir,
#     deseq2uLRT_dir
#   ),
#   labels = c(
#     "CB_LRT_unfiltered",
#     "LRT_filtered",
#     "CB_LRT_filtered",
#     "LRT_unfiltered"
#   )
# )

config <- create_config(
  dirs = c(
    deseq2LRT_dir,
    crv9_deseq2_fltrd_dir,
    cb_deseq2_fltrd_dir,
    cb_deseq2_MAD_fltrd_dir
  ),
  labels = c(
    "og_cellranger_LRT_fltrd",
    "crv9_LRT_fltrd",
    "CB_LRT_filtered",
    "CB_MAD_LRT_filtered"
  )
)

# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellranger.v.9.0.0/stats/"
output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0_MAD/stats/"

if (!dir.exists(paste0(output_fig_dir,'venn_diagrams/'))) { dir.create(paste0(output_fig_dir,'venn_diagrams/'),
                                                                         recursive=TRUE)}

all_signif_genes <- run_venn_analysis(config, output_fig_dir)

