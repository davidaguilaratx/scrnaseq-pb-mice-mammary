# Load required libraries
library(Seurat)
library(speckle)
library(limma)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(qs2)
library(RColorBrewer)

# Setting threads for parallel processing
nthreads = 2

# Define directory structure from your existing setup
analysis <- 'cellbender_analysis'
fig_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")
stats_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/DA_analysis/")

# Make sure directories exist
dirs = c(fig_dir_integrated, data_dir_integrated, annot_dir_integrated, stats_dir_integrated)
for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir, recursive=TRUE) }
}

# Load the annotated Seurat object
# This is your existing code, I'm keeping it as is
annotated_seurat = qs_read(paste0(data_dir_integrated,'de_seurat_annotated.qs'), nthreads=nthreads)
DefaultAssay(annotated_seurat) <- "RNA"
Idents(annotated_seurat) <- 'celltype'

table(annotated_seurat$timepoint, annotated_seurat$celltype)


# Remove clusters with only one condition 
annotated_seurat = subset(annotated_seurat, 
                          idents=c('Muscle'),
                          invert=TRUE)

table(annotated_seurat$sample, annotated_seurat$celltype)

annotated_seurat$celltype = Idents(annotated_seurat) # reassign with clusters removed

# Extract unique timepoints and conditions
timepoints <- unique(annotated_seurat$timepoint)
conditions <- unique(annotated_seurat$condition)

# Create directory for DA analysis results
da_dir <- paste0(stats_dir_integrated, "DA_res/")
if (!dir.exists(da_dir)) { dir.create(da_dir, recursive=TRUE) }

# Create directory for DA figures
da_fig_dir <- paste0(stats_dir_integrated, "DA_plots/")
if (!dir.exists(da_fig_dir)) { dir.create(da_fig_dir, recursive=TRUE) }

# 1. Overall differential abundance between conditions 
p_overall <- propeller(annotated_seurat,
                       clusters = annotated_seurat$celltype, 
                       sample = annotated_seurat$sample, 
                       group = annotated_seurat$condition, 
                       transform='logit')

# Save the overall results
write.csv(p_overall, paste0(da_dir, 'propeller_overall.csv'), row=TRUE)

# 2. Differential abundance at each timepoint
da_results_by_timepoint <- list()

for(tp in timepoints) {
  # Subset data for this timepoint
  seurat_subset <- subset(annotated_seurat, subset = timepoint == tp)
  
  # Check if we have multiple samples for each condition
  condition_counts <- table(seurat_subset$sample, seurat_subset$condition)
  if(ncol(condition_counts) < 2 || any(colSums(condition_counts > 0) < 2)) {
    cat("Skipping timepoint", tp, "- insufficient samples for comparison\n")
    next
  }
  
  # Run propeller analysis
  da_result <- propeller(seurat_subset,
    clusters = seurat_subset$celltype, 
    sample = seurat_subset$sample,
    group = seurat_subset$condition, 
    transform = 'logit'
  )
  
  # Add to results list
  da_results_by_timepoint[[tp]] <- da_result
  
  # Save to file
  write.csv(da_result, paste0(da_dir, 'propeller_', tp, '.csv'), row=TRUE)
  
  # Print summary of significant results
  cat("Differential abundance results for timepoint:", tp, "\n")
  print(da_result[da_result$FDR < 0.05, ])
  cat("\n")
}

# 3. Differential abundance across timepoints (for each condition)
da_results_across_time <- list()

for(cond in conditions) {
  # Subset data for this condition
  seurat_subset <- subset(annotated_seurat, subset = condition == cond)
  
  # Check if we have multiple timepoints with samples
  time_counts <- table(seurat_subset$sample, seurat_subset$timepoint)
  if(ncol(time_counts) < 2 || any(colSums(time_counts > 0) < 2)) {
    cat("Skipping condition", cond, "- insufficient timepoints for comparison\n")
    next
  }
  
  # Run propeller analysis
  da_result <- propeller(seurat_subset,
    clusters = seurat_subset$celltype, 
    sample = seurat_subset$sample,
    group = seurat_subset$timepoint, 
    transform = 'logit'
  )
  
  # Add to results list
  da_results_across_time[[cond]] <- da_result
  
  # Save to file
  write.csv(da_result, paste0(da_dir, 'propeller_', cond, '_across_time.csv'), row=TRUE)
  
  # Print summary of significant results
  cat("Differential abundance results across timepoints for condition:", cond, "\n")
  print(da_result[da_result$FDR < 0.05, ])
  cat("\n")
}

# 4. Global analysis across all condition-timepoint combinations
# Create a condition_timepoint variable
annotated_seurat$condition_timepoint <- paste0(annotated_seurat$condition, "_", annotated_seurat$timepoint)

# Check if we have enough samples for each combination
cond_time_counts <- table(annotated_seurat$sample, annotated_seurat$condition_timepoint)
if(ncol(cond_time_counts) >= 2 && all(colSums(cond_time_counts > 0) >= 1)) {
  # Run global propeller analysis
  da_results_global <- propeller(annotated_seurat,
    clusters = annotated_seurat$celltype, 
    sample = annotated_seurat$sample,
    group = annotated_seurat$condition_timepoint, 
    transform = 'logit'
  )
  
  # Save to file
  write.csv(da_results_global, paste0(da_dir, 'propeller_global.csv'), row=TRUE)
  
  # Print summary of significant results
  cat("Global differential abundance results:\n")
  print(da_results_global[da_results_global$FDR < 0.05, ])
  cat("\n")
} else {
  cat("Skipping global analysis - insufficient samples per condition-timepoint combination\n")
}

# 5. Visualization functions
# Function to extract cell proportions for plotting
get_cell_proportions <- function(seurat_obj) {
  # Create a metadata data frame
  meta_data <- seurat_obj@meta.data
  
  # Make sure we have condition, timepoint, and sample columns
  meta_data$condition_sample <- paste0(meta_data$condition, "_", meta_data$sample)
  meta_data$condition_timepoint <- paste0(meta_data$condition, "_", meta_data$timepoint)
  
  # Get counts per sample and cell type
  props_counts <- table(meta_data$sample, meta_data$celltype)
  
  # Convert to proportions
  props <- prop.table(props_counts, margin = 1)
  
  # Convert to data frame for ggplot
  props_df <- as.data.frame.matrix(props)
  props_df$sample <- rownames(props_df)
  
  # Add condition and timepoint information
  sample_info <- unique(meta_data[, c("sample", "condition", "timepoint")])
  rownames(sample_info) <- sample_info$sample
  
  props_df$condition <- sample_info[props_df$sample, "condition"]
  props_df$timepoint <- sample_info[props_df$sample, "timepoint"]
  
  # Reshape to long format for ggplot
  props_long <- pivot_longer(props_df, 
                             cols = -c(sample, condition, timepoint),
                             names_to = "celltype", 
                             values_to = "proportion")
  
  return(props_long)
}

# Get cell proportions data
cell_props <- get_cell_proportions(annotated_seurat)

# Get all unique cell types
celltypes <- unique(annotated_seurat$celltype)

# Create a color palette for cell types
n_colors <- length(celltypes)
color_palette <- colorRampPalette(brewer.pal(min(9, n_colors), "Set1"))(n_colors)
names(color_palette) <- celltypes

# Visualize cell type proportions across timepoints

# For each cell type, create a line plot showing changes over time
for(ct in celltypes) {
  # Filter data for this cell type
  ct_data <- filter(cell_props, celltype == ct)
  
  # Calculate mean and standard error per condition and timepoint
  ct_summary <- ct_data %>%
    group_by(condition, timepoint) %>%
    summarise(
      mean_prop = mean(proportion),
      se_prop = sd(proportion) / sqrt(n()),
      n = n(),
      .groups = 'drop'
    )
  
  # Create the line plot
  p <- ggplot(ct_summary, aes(x = timepoint, y = mean_prop, 
                              group = condition, color = condition)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_prop - se_prop, 
                      ymax = mean_prop + se_prop), 
                  width = 0.2) +
    labs(
      title = paste0("Proportion of ", ct, " cells across timepoints"),
      subtitle = "Error bars show standard error of the mean",
      x = "Timepoint",
      y = "Proportion of cells",
      color = "Condition"
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Save the plot
  ggsave(paste0(da_fig_dir, "proportion_", gsub("/", "_", ct), "_by_time.png"), 
         plot = p, width = 8, height = 6, dpi = 300)
}

# Create boxplots for each timepoint showing cell type proportions by condition
for(tp in timepoints) {
  # Filter data for this timepoint
  tp_data <- filter(cell_props, timepoint == tp)
  
  # Create the boxplot
  p <- ggplot(tp_data, aes(x = celltype, y = proportion, fill = condition)) +
    geom_boxplot(position = position_dodge(width = 0.9)) +
    labs(
      title = paste0("Cell type proportions at timepoint ", tp),
      x = "Cell Type",
      y = "Proportion",
      fill = "Condition"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
  
  # Save the plot
  ggsave(paste0(da_fig_dir, "proportion_boxplot_timepoint_", tp, ".png"), 
         plot = p, width = max(8, length(celltypes) * 0.5), height = 7, dpi = 300)
}

# Create a heatmap of cell type proportions across conditions and timepoints
# Calculate mean proportions per cell type, condition, and timepoint
heatmap_data <- cell_props %>%
  group_by(condition, timepoint, celltype) %>%
  summarise(mean_proportion = mean(proportion), .groups = 'drop')

# Create the heatmap
p_heatmap <- ggplot(heatmap_data, 
                    aes(x = timepoint, y = celltype, fill = mean_proportion)) +
  geom_tile() +
  facet_wrap(~condition) +
  scale_fill_viridis_c(name = "Mean\nProportion") +
  labs(
    title = "Cell Type Proportions Across Timepoints and Conditions",
    x = "Timepoint",
    y = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

# Save the heatmap
ggsave(paste0(da_fig_dir, "proportion_heatmap.png"), 
       plot = p_heatmap, width = 12, height = 10, dpi = 300)

# Create stacked bar plots to show overall composition at each timepoint
# Calculate mean proportions per cell type, condition, and timepoint
barplot_data <- cell_props %>%
  group_by(condition, timepoint, celltype) %>%
  summarise(mean_proportion = mean(proportion), .groups = 'drop')

# Create the stacked barplot
p_barplot <- ggplot(barplot_data, 
                    aes(x = timepoint, y = mean_proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~condition) +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Cell Type Composition at Each Timepoint",
    x = "Timepoint",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Save the stacked barplot
ggsave(paste0(da_fig_dir, "proportion_stacked_barplot.png"), 
       plot = p_barplot, width = 12, height = 8, dpi = 300)

# Create the standard proportion plots
p1 <- plotCellTypeProps(clusters = annotated_seurat$celltype, 
                        sample = annotated_seurat$sample)
ggsave(paste0(da_fig_dir, 'cell_type_proportions_by_sample.png'), 
       plot = p1, width = 20, height = 10, dpi = 320)

p2 <- plotCellTypeProps(clusters = annotated_seurat$celltype, 
                        sample = annotated_seurat$condition)
ggsave(paste0(da_fig_dir, 'cell_type_proportions_by_condition.png'), 
       plot = p2, width = 10, height = 10, dpi = 320)

# Create a plot showing proportions by timepoint
p3 <- plotCellTypeProps(clusters = annotated_seurat$celltype, 
                        sample = annotated_seurat$timepoint)
ggsave(paste0(da_fig_dir, 'cell_type_proportions_by_timepoint.png'), 
       plot = p3, width = 10, height = 10, dpi = 320)

# Create a plot showing proportions by condition and timepoint
p4 <- plotCellTypeProps(clusters = annotated_seurat$celltype, 
                        sample = annotated_seurat$condition_timepoint)
ggsave(paste0(da_fig_dir, 'cell_type_proportions_by_condition_timepoint.png'), 
       plot = p4, width = 20, height = 10, dpi = 320)

# Print a summary of the analysis
cat("\n=== Differential Abundance Analysis Summary ===\n")
cat("Overall DA between conditions: ", sum(p_overall$FDR < 0.05), " significant cell types\n")

# Summary for timepoint-specific analyses
for(tp in names(da_results_by_timepoint)) {
  sig_count <- sum(da_results_by_timepoint[[tp]]$FDR < 0.05)
  cat("Timepoint ", tp, ": ", sig_count, " differentially abundant cell types\n", sep="")
}

# Summary for condition-specific time analyses
for(cond in names(da_results_across_time)) {
  sig_count <- sum(da_results_across_time[[cond]]$FDR < 0.05)
  cat("Condition ", cond, " across time: ", sig_count, " differentially abundant cell types\n", sep="")
}

# Print completion message
cat("\nDifferential abundance analysis complete!\n")
cat("Results saved to:", da_dir, "\n")
cat("Plots saved to:", da_fig_dir, "\n")



# complex DA analysis ----


props <- getTransformedProps(annotated_seurat$celltype, annotated_seurat$sample, transform="logit")


plotCellTypeMeanVar(props$Counts)
plotCellTypePropsMeanVar(props$Counts)

# Create a sample info data frame for design matrix
# sample_info <- unique(annotated_seurat@meta.data[, c("sample", "condition", "timepoint", "batch")])
sample_info <- unique(annotated_seurat@meta.data[, c("sample", "condition", "timepoint")])
rownames(sample_info) <- sample_info$sample

# Make sure factors are properly set up
sample_info$condition <- factor(sample_info$condition)
sample_info$timepoint <- factor(sample_info$timepoint, ordered = F)  # Ordered factor for timepoint
# sample_info$batch <- factor(sample_info$batch)

# design <- model.matrix(~ 0 + condition + batch + timepoint + condition:timepoint, data = sample_info)
design <- model.matrix(~ 0 + condition + timepoint + condition:timepoint, data = sample_info)
colnames(design) <- make.names(colnames(design)) # Make sure column names are valid
print("Design matrix column names:")
print(colnames(design))
print("Design matrix dimensions:")
print(dim(design))

# mycontr <- makeContrasts(conditionctrl-conditionpb, levels=design)
# mycontr <- makeContrasts(conditionpb.timepointwk10 - conditionpb.timepointmn7,  levels=design)


propeller.ttest(props, design, contrasts = mycontr, robust=TRUE, trend=FALSE, 
                sort=TRUE)




library(sccomp)










library(betareg)









