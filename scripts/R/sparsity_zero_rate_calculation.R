library(future)
library(future.apply)

# Set up parallel processing using 5 cores
plan(multisession, workers = 4)

# Get the counts matrix
counts_matrix <- LayerData(seurat_umap_sct, assay='RNA', layer = "counts")

# Calculate zero rate
zero_rates <- rowSums(counts_matrix == 0) / ncol(counts_matrix)

# Identify genes with zero rate <= 0.95 (keep these)
genes_to_keep <- names(zero_rates)[zero_rates <= 0.95]
cat("Keeping", length(genes_to_keep), "out of", nrow(counts_matrix), "genes\n")

# Create a filtered matrix with only the genes to keep
filtered_counts <- counts_matrix[genes_to_keep, ]

# Process in parallel - this uses the matrix approach but splits columns
cells <- colnames(filtered_counts)
chunk_size <- ceiling(length(cells) / 4)  # Divide by number of workers
cell_chunks <- split(cells, ceiling(seq_along(cells) / chunk_size))

calc_chunk_avg <- function(cell_ids) {
  # Extract the subset of cells
  subset_idx <- match(cell_ids, colnames(filtered_counts))
  cell_matrix <- filtered_counts[, subset_idx, drop = FALSE]
  
  # Calculate column sums and number of non-zero entries
  col_sums <- Matrix::colSums(cell_matrix)
  nonzero_counts <- Matrix::colSums(cell_matrix > 0)
  
  # Calculate average non-zero value per cell
  result <- ifelse(nonzero_counts > 0, col_sums / nonzero_counts, 0)
  return(result)
}

# Process all chunks in parallel
chunk_results <- future_lapply(cell_chunks, calc_chunk_avg, future.seed = TRUE)

# Combine results
avg_nonzero_counts <- unlist(chunk_results)
names(avg_nonzero_counts) <- unlist(cell_chunks)

summary(avg_nonzero_counts)

# Add to metadata
seurat_umap_sct$avg_nonzero_expr <- avg_nonzero_counts[colnames(seurat_umap_sct)]
# 
# # View summary statistics
summary(seurat_umap_sct$avg_nonzero_expr)

# Reset to sequential processing when done
plan(sequential)


# Calculate overall sparsity (proportion of zeros in the entire matrix)
overall_sparsity <- sum(counts_matrix == 0) / (nrow(counts_matrix) * ncol(counts_matrix))
overall_sparsity_percentage <- overall_sparsity * 100

# Print the result
cat("Overall sparsity: ", round(overall_sparsity_percentage, 2), "%\n")

# Calculate zero rate for each gene
gene_zero_rates <- rowSums(counts_matrix == 0) / ncol(counts_matrix)
gene_zero_percentages <- gene_zero_rates * 100

# Summary of gene sparsity
summary(gene_zero_percentages)

# Histogram of gene sparsity
hist(gene_zero_percentages, 
     main = "Distribution of Gene Sparsity", 
     xlab = "Percentage of Zeros", 
     breaks = 50)



# Calculate zero rate for each cell
cell_zero_rates <- colSums(counts_matrix == 0) / nrow(counts_matrix)
cell_zero_percentages <- cell_zero_rates * 100

# Summary of cell sparsity
summary(cell_zero_percentages)

# Add to Seurat metadata
seurat_umap_sct$zero_percentage <- cell_zero_percentages

# Visualize on UMAP
FeaturePlot(seurat_umap_sct, features = "zero_percentage")
