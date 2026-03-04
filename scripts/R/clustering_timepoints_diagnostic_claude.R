diagnostic_integration_decision <- function(output_dir = "results/00_diagnostics") {
  
  # Detach plyr if loaded (causes conflicts with dplyr)
  if ("package:plyr" %in% search()) {
    detach("package:plyr", unload = TRUE)
    message("Detached plyr to avoid conflicts with dplyr")
  }
  
  # Load required packages
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(glue)
    library(Seurat)
    library(qs2)
  })
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load all timepoints (QC'd but not integrated)
  message("Loading timepoints...")
  timepoints <- list(
    wk3 = qs_read("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/merged_seurat_cb.qs"),
    wk10 = qs_read("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/10_weeks/cellbender_analysis/FPR_0.0_MAD/data/merged_seurat_cb.qs"),
    mn7 = qs_read("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/7_months/cellbender_analysis/FPR_0.0_MAD/data/merged_seurat_cb.qs"),
    mn18 = qs_read("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/cellbender_analysis/FPR_0.0_MAD/data/merged_seurat_cb.qs")
  )
  
  # Store results separately
  timepoints_clustered <- list()
  all_markers_list <- list()
  
  # Quick clustering of each
  for (tp in names(timepoints)) {
    message(glue("\nQuick clustering {tp}..."))
    
    obj <- timepoints[[tp]]
    
    # Normalize and cluster
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
    
    # Find markers
    message(glue("Finding markers for {tp}..."))
    markers <- FindAllMarkers(obj, 
                              only.pos = TRUE, 
                              max.cells.per.ident = 200,  # Subsample for speed
                              verbose = FALSE)
    
    # Store separately
    timepoints_clustered[[tp]] <- obj
    all_markers_list[[tp]] <- markers
    
    message(glue("✓ {tp} complete: {ncol(obj)} cells, {nlevels(Idents(obj))} clusters"))
  }
  
  # ==========================================================================
  # VISUALIZATION 1: Side-by-side UMAPs
  # ==========================================================================
  
  message("\nGenerating UMAP plots...")
  
  pdf(file.path(output_dir, "separate_timepoint_UMAPs.pdf"), 
      width = 16, height = 12)
  
  plots <- lapply(names(timepoints_clustered), function(tp) {
    DimPlot(timepoints_clustered[[tp]], label = TRUE, label.size = 5) +
      ggtitle(glue("{tp}: {ncol(timepoints_clustered[[tp]])} cells, {nlevels(Idents(timepoints_clustered[[tp]]))} clusters")) +
      NoLegend() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  })
  
  print(wrap_plots(plots, ncol = 2))
  dev.off()
  
  # ==========================================================================
  # VISUALIZATION 2: Cluster sizes comparison
  # ==========================================================================
  
  message("Comparing cluster sizes...")
  
  cluster_sizes <- lapply(names(timepoints_clustered), function(tp) {
    cluster_counts <- table(Idents(timepoints_clustered[[tp]]))
    data.frame(
      timepoint = tp,
      cluster = names(cluster_counts),
      n_cells = as.numeric(cluster_counts),
      stringsAsFactors = FALSE
    )
  }) %>% dplyr::bind_rows()
  
  pdf(file.path(output_dir, "cluster_size_comparison.pdf"), 
      width = 12, height = 8)
  
  p <- ggplot(cluster_sizes, aes(x = cluster, y = n_cells, fill = timepoint)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Cluster Sizes Across Timepoints",
         x = "Cluster",
         y = "Number of Cells") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  dev.off()
  
  # ==========================================================================
  # ANALYSIS 1: Top markers per timepoint
  # ==========================================================================
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("TOP MARKERS PER TIMEPOINT\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  for (tp in names(all_markers_list)) {
    cat(glue("\n{tp}:\n"))
    cat(rep("-", 70), "\n", sep = "")
    
    markers <- all_markers_list[[tp]]
    
    # Save full marker list
    write.csv(markers,
              file.path(output_dir, glue("{tp}_all_markers.csv")),
              row.names = FALSE)
    
    # Show top 3 per cluster using base R to avoid dplyr issues
    sig_markers <- markers[markers$p_val_adj < 0.05, ]
    
    if (nrow(sig_markers) > 0) {
      # Split by cluster
      marker_by_cluster <- split(sig_markers, sig_markers$cluster)
      
      # Get top 3 per cluster
      top_markers_list <- lapply(names(marker_by_cluster), function(cl) {
        cl_markers <- marker_by_cluster[[cl]]
        cl_markers <- cl_markers[order(-cl_markers$avg_log2FC), ]
        top_genes <- head(cl_markers$gene, 3)
        data.frame(
          cluster = cl,
          n_markers = nrow(cl_markers),
          top_genes = paste(top_genes, collapse = ", "),
          stringsAsFactors = FALSE
        )
      })
      
      top_markers_df <- do.call(rbind, top_markers_list)
      print(top_markers_df)
    } else {
      cat("No significant markers found\n")
    }
  }
  
  # ==========================================================================
  # ANALYSIS 2: Marker gene overlap (Jaccard similarity)
  # ==========================================================================
  
  cat("\n\n", rep("=", 70), "\n", sep = "")
  cat("MARKER GENE OVERLAP ANALYSIS\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Extract significant marker genes
  all_markers <- lapply(all_markers_list, function(markers) {
    sig_markers <- markers[markers$p_val_adj < 0.05, ]
    unique(sig_markers$gene)
  })
  
  # Calculate Jaccard similarity
  jaccard_matrix <- matrix(0, 
                           nrow = length(all_markers), 
                           ncol = length(all_markers),
                           dimnames = list(names(all_markers), names(all_markers)))
  
  for (i in 1:length(all_markers)) {
    for (j in 1:length(all_markers)) {
      if (i == j) {
        jaccard_matrix[i, j] <- 1
      } else {
        intersection <- length(intersect(all_markers[[i]], all_markers[[j]]))
        union <- length(union(all_markers[[i]], all_markers[[j]]))
        jaccard_matrix[i, j] <- intersection / union
      }
    }
  }
  
  cat("Jaccard Similarity Matrix:\n")
  print(round(jaccard_matrix, 3))
  cat("\n")
  
  # Save matrix
  write.csv(jaccard_matrix,
            file.path(output_dir, "jaccard_similarity_matrix.csv"))
  
  # Visualize as heatmap
  pdf(file.path(output_dir, "jaccard_similarity_heatmap.pdf"), 
      width = 8, height = 7)
  
  library(pheatmap)
  pheatmap(jaccard_matrix,
           display_numbers = TRUE,
           number_format = "%.3f",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "blue"))(50),
           main = "Marker Gene Overlap (Jaccard Similarity)")
  
  dev.off()
  
  # ==========================================================================
  # ANALYSIS 3: Unique vs shared markers
  # ==========================================================================
  
  cat("\nMarker Gene Statistics:\n")
  cat(rep("-", 70), "\n", sep = "")
  
  for (tp in names(all_markers)) {
    n_total <- length(all_markers[[tp]])
    
    # How many are unique to this timepoint?
    other_timepoints <- setdiff(names(all_markers), tp)
    other_markers <- unique(unlist(all_markers[other_timepoints]))
    n_unique <- length(setdiff(all_markers[[tp]], other_markers))
    
    cat(glue("{tp}: {n_total} total markers, {n_unique} unique ({round(100*n_unique/n_total, 1)}%)\n"))
  }
  
  # ==========================================================================
  # RECOMMENDATION
  # ==========================================================================
  
  cat("\n\n", rep("=", 70), "\n", sep = "")
  cat("INTEGRATION STRATEGY RECOMMENDATION\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  # Average pairwise Jaccard (excluding diagonal)
  avg_jaccard <- mean(jaccard_matrix[lower.tri(jaccard_matrix)])
  
  cat(glue("Average pairwise Jaccard similarity: {round(avg_jaccard, 3)}\n\n"))
  
  if (avg_jaccard > 0.6) {
    cat("✓ HIGH SIMILARITY (Jaccard > 0.6)\n")
    cat("→ RECOMMENDATION: Strategy A (Integrate → Cluster)\n")
    cat("  - Timepoints are biologically similar\n")
    cat("  - Integration will work well\n")
    cat("  - Can identify common cell types easily\n\n")
    recommendation <- "Strategy_A"
    
  } else if (avg_jaccard > 0.4) {
    cat("⚠ MODERATE SIMILARITY (0.4 < Jaccard < 0.6)\n")
    cat("→ RECOMMENDATION: Strategy B (Cluster → Transfer → Integrate)\n")
    cat("  - Some differences between timepoints\n")
    cat("  - Cluster separately first to capture stage-specific biology\n")
    cat("  - Use label transfer to harmonize annotations\n")
    cat("  - Then integrate for cross-timepoint comparisons\n\n")
    recommendation <- "Strategy_B"
    
  } else {
    cat("⚠⚠ LOW SIMILARITY (Jaccard < 0.4)\n")
    cat("→ RECOMMENDATION: Strategy C (Analyze Independently)\n")
    cat("  - Major biological differences between timepoints\n")
    cat("  - Integration may mask important biology\n")
    cat("  - Analyze each timepoint separately\n")
    cat("  - Compare results post-hoc\n\n")
    recommendation <- "Strategy_C"
  }
  
  # Check cluster numbers
  n_clusters <- sapply(timepoints_clustered, function(x) nlevels(Idents(x)))
  cat(glue("Number of clusters per timepoint: {paste(names(n_clusters), '=', n_clusters, collapse = ', ')}\n"))
  
  if (max(n_clusters) / min(n_clusters) > 2) {
    cat("⚠ Large variation in cluster numbers suggests different cell states\n")
    cat("  → Consider Strategy B or C\n\n")
  } else {
    cat("✓ Similar number of clusters across timepoints\n")
    cat("  → Integration may work well\n\n")
  }
  
  # Save recommendation
  writeLines(
    c(
      glue("Recommendation: {recommendation}"),
      glue("Average Jaccard: {round(avg_jaccard, 3)}"),
      glue("Cluster counts: {paste(names(n_clusters), '=', n_clusters, collapse = ', ')}")
    ),
    file.path(output_dir, "recommendation.txt")
  )
  
  # ==========================================================================
  # Return results
  # ==========================================================================
  
  message("\n✓ Diagnostic complete!")
  message(glue("Results saved to: {output_dir}"))
  message("\nReview:")
  message("  1. separate_timepoint_UMAPs.pdf - Visual comparison")
  message("  2. jaccard_similarity_heatmap.pdf - Marker overlap")
  message("  3. *_all_markers.csv - Full marker lists per timepoint")
  message("  4. recommendation.txt - Integration strategy")
  
  return(list(
    seurat_objects = timepoints_clustered,
    markers = all_markers_list,
    jaccard_matrix = jaccard_matrix,
    recommendation = recommendation,
    avg_jaccard = avg_jaccard
  ))
}

# Run the diagnostic
results <- diagnostic_integration_decision()

# View recommendation
cat("\n=== RECOMMENDATION ===\n")
cat("Strategy:", results$recommendation, "\n")
cat("Average Jaccard:", round(results$avg_jaccard, 3), "\n")

# Access results
# results$seurat_objects$wk3  # Clustered Seurat object for week 3
# results$markers$wk3          # Marker genes for week 3
# results$jaccard_matrix       # Similarity matrix
# results$recommendation       # Which strategy to use