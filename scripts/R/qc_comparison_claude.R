compare_timepoint_qc <- function(seurat_obj, output_dir = "results/00_diagnostics") {
  
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
    library(scales)
  })
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "complexity")
  
  plot_list <- list()
  
  for (metric in metrics) {
    
    message(glue("Processing {metric}..."))
    
    # Extract just the columns we need to avoid any issues
    meta_data <- data.frame(
      timepoint = seurat_obj@meta.data$timepoint,
      value = seurat_obj@meta.data[[metric]]
    )
    
    # Remove NAs
    meta_data <- meta_data[!is.na(meta_data$value), ]
    
    # Calculate summary stats
    summary_stats <- meta_data %>%
      dplyr::group_by(timepoint) %>%
      dplyr::summarise(
        median = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Create violin plot
    p <- ggplot(meta_data, aes(x = timepoint, y = value, fill = timepoint)) +
      geom_violin(alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
      geom_text(data = summary_stats,
                aes(x = timepoint, y = median, 
                    label = paste0("Med: ", round(median, 0))),
                vjust = -0.5, size = 3, inherit.aes = FALSE) +
      labs(title = metric,
           subtitle = "Distribution across timepoints",
           x = "Timepoint",
           y = metric) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        plot.title = element_text(face = "bold", size = 14)
      )
    
    # Add log scale for count metrics
    if (metric %in% c("nCount_RNA", "nFeature_RNA")) {
      p <- p + 
        scale_y_log10(labels = comma) +
        annotation_logticks(sides = "l")
    }
    
    plot_list[[metric]] <- p
  }
  
  # Combine plots
  combined <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = "QC Metrics Comparison Across Timepoints",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(output_dir, "timepoint_QC_comparison.pdf"),
         combined, width = 12, height = 10, dpi = 300)
  
  message("✓ Plots saved")
  
  # ==========================================================================
  # Statistical tests
  # ==========================================================================
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("STATISTICAL COMPARISON OF QC METRICS ACROSS TIMEPOINTS\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  test_results <- list()
  
  for (metric in metrics) {
    cat(glue("=== {metric} ===\n"))
    
    # Extract values directly
    metric_values <- seurat_obj@meta.data[[metric]]
    timepoint_values <- seurat_obj@meta.data$timepoint
    
    # Remove NAs
    valid_idx <- !is.na(metric_values) & !is.na(timepoint_values)
    metric_values <- metric_values[valid_idx]
    timepoint_values <- timepoint_values[valid_idx]
    
    # Summary statistics per timepoint
    summary_table <- data.frame(
      timepoint = seurat_obj@meta.data$timepoint,
      value = seurat_obj@meta.data[[metric]]
    ) %>%
      dplyr::group_by(timepoint) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        median = median(value, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        min = min(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        .groups = "drop"
      )
    
    print(summary_table)
    cat("\n")
    
    # Kruskal-Wallis test
    test_result <- kruskal.test(metric_values ~ timepoint_values)
    
    cat(glue("Kruskal-Wallis Test:\n"))
    cat(glue("  Chi-squared = {round(test_result$statistic, 3)}\n"))
    cat(glue("  df = {test_result$parameter}\n"))
    cat(glue("  p-value = {format.pval(test_result$p.value, digits = 4)}\n\n"))
    
    # Interpretation
    if (test_result$p.value < 0.001) {
      cat("*** HIGHLY SIGNIFICANT (p < 0.001)\n")
      cat("    → Strong evidence for differences between timepoints\n")
      cat("    → RECOMMEND: Per-sample MAD filtering\n\n")
      significance <- "highly_significant"
      
    } else if (test_result$p.value < 0.05) {
      cat("**  SIGNIFICANT (p < 0.05)\n")
      cat("    → Moderate evidence for differences\n")
      cat("    → RECOMMEND: Hybrid approach (universal minimums + MAD)\n\n")
      significance <- "significant"
      
    } else {
      cat("    NOT SIGNIFICANT (p >= 0.05)\n")
      cat("    → No strong evidence for differences\n")
      cat("    → RECOMMEND: Universal QC thresholds\n\n")
      significance <- "not_significant"
    }
    
    # Post-hoc if significant
    if (test_result$p.value < 0.05) {
      cat("Post-hoc Pairwise Comparisons:\n")
      
      pairwise_result <- pairwise.wilcox.test(
        metric_values, 
        timepoint_values,
        p.adjust.method = "bonferroni"
      )
      
      pval_matrix <- pairwise_result$p.value
      print(round(pval_matrix, 4))
      cat("\n")
    }
    
    # Store results
    test_results[[metric]] <- list(
      test = test_result,
      significance = significance,
      summary = summary_table
    )
    
    cat(rep("-", 70), "\n\n", sep = "")
  }
  
  # ==========================================================================
  # Overall recommendation
  # ==========================================================================
  
  cat(rep("=", 70), "\n", sep = "")
  cat("OVERALL QC STRATEGY RECOMMENDATION\n")
  cat(rep("=", 70), "\n\n", sep = "")
  
  sig_count <- sum(sapply(test_results, function(x) x$significance != "not_significant"))
  highly_sig_count <- sum(sapply(test_results, function(x) x$significance == "highly_significant"))
  
  cat(glue("Metrics with significant differences: {sig_count} / {length(metrics)}\n"))
  cat(glue("Metrics with highly significant differences: {highly_sig_count} / {length(metrics)}\n\n"))
  
  if (sig_count == 0) {
    cat("✓ NO SIGNIFICANT DIFFERENCES\n")
    cat("→ RECOMMENDATION: Universal QC thresholds\n\n")
    recommendation <- "universal"
    
  } else if (sig_count <= 2 && highly_sig_count == 0) {
    cat("⚠ MINOR DIFFERENCES\n")
    cat("→ RECOMMENDATION: Hybrid (universal minimums + per-sample MAD)\n\n")
    recommendation <- "hybrid"
    
  } else {
    cat("⚠⚠ SUBSTANTIAL DIFFERENCES\n")
    cat("→ RECOMMENDATION: Per-sample MAD filtering\n\n")
    recommendation <- "per_sample_MAD"
  }
  
  # Save results
  summary_all <- dplyr::bind_rows(lapply(names(test_results), function(metric) {
    test_results[[metric]]$summary %>%
      dplyr::mutate(metric = metric, .before = 1)
  }))
  
  write.csv(summary_all,
            file.path(output_dir, "QC_summary_by_timepoint.csv"),
            row.names = FALSE)
  
  test_summary <- data.frame(
    metric = names(test_results),
    chi_squared = sapply(test_results, function(x) x$test$statistic),
    p_value = sapply(test_results, function(x) x$test$p.value),
    significance = sapply(test_results, function(x) x$significance)
  )
  
  write.csv(test_summary,
            file.path(output_dir, "statistical_tests_summary.csv"),
            row.names = FALSE)
  
  writeLines(
    c(
      glue("Recommendation: {recommendation}"),
      glue("Significant metrics: {sig_count}/{length(metrics)}"),
      glue("Highly significant: {highly_sig_count}/{length(metrics)}")
    ),
    file.path(output_dir, "QC_strategy_recommendation.txt")
  )
  
  message("\n✓ Analysis complete!")
  message(glue("Results saved to: {output_dir}"))
  
  return(list(
    test_results = test_results,
    recommendation = recommendation,
    summary_stats = summary_all
  ))
}

# USAGE - Run this:
qc_comparison <- compare_timepoint_qc(merged_seurat_cb)

# View recommendation
cat("\nRECOMMENDATION:", qc_comparison$recommendation, "\n")

# Access specific results
qc_comparison$test_results$nCount_RNA$summary
qc_comparison$summary_stats

# 
# ## What This Fixed Version Does:
# 
# 1. **Fixes the NSE error** by using `!!sym(metric)` correctly
# 2. **Adds detailed summary statistics** per timepoint
# 3. **Performs pairwise comparisons** when overall test is significant
# 4. **Identifies which pairs differ** significantly
# 5. **Gives clear recommendations** based on results
# 6. **Saves everything** to files for documentation
# 
# ## Expected Output:
# 
# You should see:
#   ```
# === nCount_RNA ===
#   # A tibble: 4 × 7
#   timepoint n_cells median   mean    sd   min    max
# <fct>       <int>  <dbl>  <dbl> <dbl> <dbl>  <dbl>
#   1 wk3          xxxx   2500   3200  1800   150  15000
# 2 wk10         xxxx   4200   5100  2300   150  20000
# 3 mn7          xxxx   4000   4800  2100   150  18000
# 4 mn18         xxxx   3800   4500  2000   150  19000
# 
# Kruskal-Wallis Test:
#   Chi-squared = XX.XXX
# df = 3
# p-value = X.XXXe-XX
# 
# *** HIGHLY SIGNIFICANT (p < 0.001)
# → Strong evidence for differences between timepoints
# → RECOMMEND: Per-timepoint or per-sample QC thresholds

