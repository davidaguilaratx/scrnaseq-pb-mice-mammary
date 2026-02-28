# Improved GSEA Analysis Script
# setwd("/path/to/your/project")  # Set your working directory here

# Load required libraries
library(tidyverse)
library(fgsea)
library(msigdb)
library(biomaRt)
library(ExperimentHub)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(ggthemes)
library(ggdark)



# Set seed for reproducibility
set.seed(2024)

###########################################
# Configuration Parameters
###########################################

config <- list(
  # Directories
  de_directory = "stats_rerun/DESeq2/results/",
  output_directory = "stats_rerun/GSEA/fGSEA/",
  annotation_file = "annotation/ensembl_and_entrez_gene_ids_and_human_homologs.csv",
  
  # GSEA parameters
  min_gene_set_size = 15,
  max_gene_set_size = 1000,
  padj_cutoff = 0.25,
  pvalue_threshold = 0.05,
  sample_size = 101,
  num_processors = 2,
  
  # Collections to analyze
  pathway_subcollections = c('h', 'GO:BP', 'GO:MF', 'GO:CC', 'CP:REACTOME', 'CP:KEGG', 'CP:WIKIPATHWAYS'),
  pathway_names = c('Hallmarks', 'GO: Biological Processes', 'GO: Molecular Function',
                    'GO: Cellular Compartment', 'Reactome', 'KEGG', 'Wikipathways'),
  
  # Species
  species = c("Mouse", "Human"),
  
  # Plot settings
  max_pathways = 30,
  plot_width = 10,
  plot_height = 15
)

###########################################
# Helper Functions
###########################################

# Function to load all DE results from a directory
load_all_de_results <- function(directory, annotation_file) {
  # Load gene IDs and human homologs
  gene_id <- read.csv(annotation_file)
  
  # Get all DE result files
  files <- list.files(path = directory, pattern = "all_genes_deseq2.csv", full.names = TRUE)
  
  # Extract cell types from filenames
  celltypes <- sub(pattern = '_condition.*', replacement = '', x = basename(files))
  
  # Create a named list to store results
  results_list <- setNames(vector("list", length(files)), celltypes)
  
  # Load each file and add gene IDs
  for (i in seq_along(files)) {
    message(paste0("Loading ", celltypes[i], " data..."))
    results_list[[celltypes[i]]] <- get_ids_and_homologs(read.csv(files[i]), gene_ids = gene_id)
  }
  
  return(results_list)
}

# Function to merge DE results with gene IDs and human homologs
get_ids_and_homologs <- function(res_tbl, gene_ids) {
  res_tbl <- merge(x = res_tbl, y = gene_ids, by.x = "gene", by.y = "external_gene_name")
  return(res_tbl)
}

# Function to convert DE results to a ranked gene list
to_genelist <- function(res_tbl, species_human = FALSE) {
  if (!species_human) {
    genelist <- res_tbl %>%
      as_tibble() %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
      dplyr::select(gene, ranking_metric) %>%
      na.omit() %>%
      arrange(desc(ranking_metric)) %>%
      distinct() %>%
      deframe()
  } else {
    genelist <- res_tbl %>%
      as_tibble() %>%
      group_by(hsapiens_homolog_associated_gene_name) %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
      dplyr::select(hsapiens_homolog_associated_gene_name, ranking_metric) %>%
      na.omit() %>%
      arrange(desc(ranking_metric)) %>%
      distinct() %>%
      deframe()
  }
  
  return(genelist)
}

# Additional helper functions are in the separate file
# Include them or source them before running this script

# Function to run fGSEA analysis for a given dataset and pathway collection
run_fgsea_analysis <- function(data, celltype, msigdb_obj, collection_type, collection_name, 
                               species, config) {
  # Create proper output directories
  results_dir <- file.path(config$output_directory, "results")
  figures_dir <- file.path(config$output_directory, "figures")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate a ranked gene list
  genelist <- to_genelist(data, species_human = (species == "Human"))
  
  # Get appropriate pathway collection
  if (collection_type == 'h') {
    pathways <- geneIds(subsetCollection(msigdb_obj, collection = collection_type))
  } else {
    pathways <- geneIds(subsetCollection(msigdb_obj, subcollection = collection_type))
  }
  
  # Run fGSEA
  message(paste0("Running fGSEA for ", celltype, " using ", species, " ", collection_name, " pathways..."))
  fgseaRes <- fgsea(pathways = pathways,
                    stats = genelist,
                    sampleSize = config$sample_size,
                    minSize = config$min_gene_set_size,
                    maxSize = config$max_gene_set_size,
                    eps = 1e-50,
                    nproc = config$num_processors)
  
  # Save results
  pathway_name <- paste0(species, " ", collection_name)
  results_file <- tidysave_fgsea(fgseaRes, 
                                 pathways = pathway_name, 
                                 celltype = celltype,
                                 directory = paste0(results_dir, "/"))
  
  # Collapse pathways to reduce redundancy (except for hallmarks)
  message(paste0("Collapsing pathways for ", celltype, " ", pathway_name))
  
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < config$padj_cutoff], 
                                        pathways, genelist,
                                        pval.threshold = config$pvalue_threshold,
                                        nperm = 10000)
  
  # Determine which pathways to use for plotting
  if (collection_type == 'h') {
    mainPathways <- names(pathways)  # Don't modify if hallmarks
  } else {
    mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
  }
  
  # Save collapsed pathways
  saveRDS(collapsedPathways, 
          file.path(results_dir, paste0(celltype, "_", pathway_name, "_collapsedPathways.rds")))
  
  # Create plots
  message(paste0("Creating plots for ", celltype, " using ", pathway_name, " pathways..."))
  
  # Determine which results to plot
  plot_data <- if(length(mainPathways) != 0) {
    fgseaRes[pathway %in% mainPathways,]
  } else {
    fgseaRes
  }
  
  # Create barplot
  barplot <- plot_fgsea_nes(fgsea_nes = plot_data,
                            pathways = pathway_name,
                            celltype = celltype,
                            padj_cutoff = config$padj_cutoff,
                            max_pathways = config$max_pathways,
                            plot_width = config$plot_width,
                            plot_height = config$plot_height,
                            directory = paste0(figures_dir, "/"))
  
  # Create dotplot  
  dotplot <- plot_fgsea_nes_dotplot(fgsea_nes = plot_data,
                                    pathways = pathway_name,
                                    celltype = celltype,
                                    padj_cutoff = config$padj_cutoff,
                                    max_pathways = config$max_pathways,
                                    plot_width = config$plot_width,
                                    plot_height = config$plot_height,
                                    directory = paste0(figures_dir, "/"))
  
  # Create enrichment plot
  update_geom_defaults("segment", list(color = "black"))
  
  plot_pathways <- if(length(pathways) > 50) {
    if(length(mainPathways) != 0) {
      pathways[mainPathways]
    } else {
      pathways
    }
  } else {
    pathways
  }
  
  p <- plotGseaTable(pathways = plot_pathways,
                     stats = genelist, 
                     fgseaRes = fgseaRes,
                     pathwayLabelStyle = list(color = 'black'),
                     headerLabelStyle = list(size = 20, color = 'black'),
                     valueStyle = list(color = 'black'),
                     axisLabelStyle = list(color = 'black'),
                     gseaParam = 0.5)
  
  enrichment_plot <- cowplot::ggdraw(p) +
    theme(plot.background = element_rect(fill = 'white', color = NA),
          panel.background = element_rect(fill = 'white', color = NA))
  
  # Save enrichment plot  
  ggsave(file.path(figures_dir, paste0(celltype, "_multilevel_fGSEA_NES_", pathway_name, ".png")),
         plot = enrichment_plot, dpi = 320, width = 10, height = 15, units = 'in')
  
  # Return results for further use if needed
  return(list(
    fgsea_results = fgseaRes,
    collapsed_pathways = collapsedPathways,
    main_pathways = mainPathways,
    barplot = barplot,
    dotplot = dotplot,
    enrichment_plot = enrichment_plot
  ))
}

# Main Execution Function
###########################################

run_gsea_workflow <- function(config) {
  # Load DE results
  de_results <- load_all_de_results(config$de_directory, config$annotation_file)
  
  # Import MSigDB databases
  message("Loading MSigDB databases...")
  msigdb.mm <- getMsigdb(org = 'mm', id = 'SYM')
  msigdb.mm <- appendKEGG(msigdb.mm)
  
  msigdb.hs <- getMsigdb(org = 'hs', id = 'SYM')
  msigdb.hs <- appendKEGG(msigdb.hs)
  
  # Create result storage
  all_results <- list()
  
  # Loop through each cell type
  for (celltype in names(de_results)) {
    message(paste0("\n\nProcessing ", celltype, "...\n"))
    all_results[[celltype]] <- list()
    
    # Loop through species (mouse, human)
    for (j in 1:2) {
      species <- config$species[j]
      msigdb_obj <- if(j == 1) msigdb.mm else msigdb.hs
      
      all_results[[celltype]][[species]] <- list()
      
      # Loop through pathway collections
      for (k in 1:length(config$pathway_subcollections)) {
        collection_type <- config$pathway_subcollections[k]
        collection_name <- config$pathway_names[k]
        
        # Run analysis
        analysis_results <- run_fgsea_analysis(
          data = de_results[[celltype]],
          celltype = celltype,
          msigdb_obj = msigdb_obj,
          collection_type = collection_type,
          collection_name = collection_name,
          species = species,
          config = config
        )
        
        # Store results
        all_results[[celltype]][[species]][[collection_name]] <- analysis_results
      }
    }
  }
  
  # Create comparison plots across cell types
  message("\n\nCreating comparison plots across cell types...")
  
  figures_dir <- file.path(config$output_directory, "figures")
  
  for (k in 1:length(config$pathway_names)) {
    collection_name <- config$pathway_names[k]
    
    # Create dotplot comparison
    plot_multicelltype_comparison(all_results, collection_name, 
                                  pathways_to_show = 10, 
                                  directory = paste0(figures_dir, "/"))
    
    # Create heatmap
    plot_pathway_heatmap(all_results, collection_name,
                         directory = paste0(figures_dir, "/"))
  }
  
  # Save complete results
  saveRDS(all_results, file.path(config$output_directory, "complete_gsea_results.rds"))
  return(all_results)
} 