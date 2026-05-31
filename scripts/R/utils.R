library(parallel)
library(BiocParallel)
library(Seurat)
library(yaml)
library(RhpcBLASctl)
library(scDblFinder)
library(SingleCellExperiment)
library(scCustomize)
library(harmony)

# ============================================================
# Parallelism Configuration ----
# ============================================================

#' Configure parallelism based on OS and compute environment
#'
#' Detects whether the session is running on Windows, Linux/macOS,
#' or an HPC (via SLURM environment variables) and returns an
#' appropriate BiocParallel BPPARAM object and thread count.
#'
#' @return A named list with:
#'   - `nthreads`: integer number of workers
#'   - `bpparam`:  a BiocParallelParam object
configure_parallelism <- function(rng_seed = NULL, max_workers = 8L) {
  os <- Sys.info()[["sysname"]]
  is_hpc <- nchar(Sys.getenv("SLURM_JOB_ID")) > 0
  max_hpc_threads <- max_workers # diminishing returns beyond ~8 for BiocParallel / qs2 I/O?

  if (is_hpc) {
    # Prefer SLURM_CPUS_PER_TASK (--cpus-per-task), fall back to SLURM_NTASKS
    slurm_cpus <- suppressWarnings(as.integer(Sys.getenv(
      "SLURM_CPUS_ON_NODE",
      unset = NA
    )))
    slurm_tasks <- suppressWarnings(as.integer(Sys.getenv(
      "SLURM_NTASKS",
      unset = NA
    )))
    detected <- Filter(Negate(is.na), c(slurm_cpus, slurm_tasks, 1L))[[1]]
    nthreads <- min(max(1L, detected - 2L), max_hpc_threads)
    bpparam <- BiocParallel::MulticoreParam(
      workers = nthreads,
      RNGseed = rng_seed,
      progressbar = TRUE
    )
  } else if (os == "Windows") {
    # Forking is unavailable on Windows; use SNOW-based parallelism
    physical_cores <- parallel::detectCores(logical = FALSE)
    nthreads <- max(1L, floor(physical_cores / 4L))
    bpparam <- BiocParallel::SnowParam(
      workers = nthreads,
      RNGseed = rng_seed,
      progressbar = TRUE
    )
  } else {
    # Linux / macOS — fork-based, leave 2 cores free for the OS
    physical_cores <- parallel::detectCores(logical = FALSE)
    nthreads <- max(1L, physical_cores - 2L)
    bpparam <- BiocParallel::MulticoreParam(
      workers = nthreads,
      RNGseed = rng_seed,
      progressbar = TRUE
    )
  }

  message(sprintf(
    "[configure_parallelism] OS: %s | HPC: %s | workers: %d | backend: %s",
    os,
    is_hpc,
    nthreads,
    class(bpparam)
  ))

  list(nthreads = nthreads, bpparam = bpparam)
}


# ============================================================
# Modes for thread pools ----
# ============================================================

threads_serial <- function(n = 4) {
  RhpcBLASctl::omp_set_num_threads(n)
  RhpcBLASctl::blas_set_num_threads(n)
  RcppParallel::setThreadOptions(numThreads = n)
  data.table::setDTthreads(n)
}

threads_capped <- function() {
  RhpcBLASctl::omp_set_num_threads(1)
  RhpcBLASctl::blas_set_num_threads(1)
  RcppParallel::setThreadOptions(numThreads = 1)
  data.table::setDTthreads(1)
}

# ============================================================
# Get yaml parameter ----
# ============================================================

#' Retrieve a parameter from a nested YAML config with timepoint-level override
#'
#' Looks up \code{param_name} first under the timepoint-specific block, then
#' falls back to the \code{default} block. Errors loudly if neither exists.
#'
#' @param configs    Named list parsed from a YAML config file via \code{yaml::read_yaml()}
#' @param param_name Character string naming the parameter to retrieve
#' @param timepoint  Character string matching a top-level key in \code{configs}
#'
#' @return The value of \code{param_name} for the given timepoint, or the default
#'
#' @examples
#' lower <- get_param(configs, "emptydrops", "lower_umi_bound", "7_months")
get_param <- function(
  configs,
  section,
  param_name,
  timepoint,
  fallback = NULL
) {
  value <- configs[[section]][[timepoint]][[param_name]] %||%
    configs[[section]]$default[[param_name]] %||%
    configs$default[[param_name]] %||%
    fallback
  if (is.null(value)) {
    stop(
      sprintf(
        "[get_param] Parameter '%s' not found in section '%s' ",
        param_name,
        section
      ),
      sprintf("for timepoint '%s' or in defaults.", timepoint)
    )
  }
  value
}


# ============================================================
# Quality metrics calculations ----
# ============================================================

calc_metrics <- function(seurat_object) {
  # ==========================================================
  # 1. Vectorized Metadata Calculations (MASSIVELY FASTER)
  # ==========================================================
  # Pull out metadata to avoid repeated S4 object copies
  meta <- seurat_object@meta.data
  original_cells <- rownames(meta)

  # Calculate Z-scores and complexity natively and vectorized
  meta <- meta |>
    dplyr::group_by(sample) |>
    dplyr::mutate(
      log10_genes = log10(nFeature_RNA),
      log10_counts = log10(nCount_RNA),
      log10_genes_zscored = (log10_genes - mean(log10_genes, na.rm = TRUE)) /
        sd(log10_genes, na.rm = TRUE),
      log10_counts_zscored = (log10_counts - mean(log10_counts, na.rm = TRUE)) /
        sd(log10_counts, na.rm = TRUE),
      complexity = log10_genes / log10_counts
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-log10_genes, -log10_counts) |>
    as.data.frame()

  # Restore rownames and inject back into the Seurat object
  rownames(meta) <- original_cells
  seurat_object@meta.data <- meta

  # Clean up memory
  rm(meta, original_cells)
  gc()

  # ==========================================================
  # 2. Percentage Feature Sets
  # ==========================================================
  seurat_object <- PercentageFeatureSet(
    seurat_object,
    pattern = "^mt-",
    col.name = "percent.mt"
  )
  seurat_object <- PercentageFeatureSet(
    seurat_object,
    pattern = "^Rp[sl]",
    col.name = "percent.rb"
  )
  seurat_object <- PercentageFeatureSet(
    seurat_object,
    pattern = "^Hb[^(P|E|S)]",
    col.name = "percent.hb"
  )

  # ==========================================================
  # 3. Cumulative proportions of top genes (Memory-Optimized)
  # ==========================================================
  layer_names <- Layers(seurat_object, assay = "RNA", search = "counts")
  qc_list <- vector("list", length(layer_names))
  names(qc_list) <- layer_names

  for (lyr in layer_names) {
    # Extract just this layer
    counts_matrix <- LayerData(seurat_object, assay = "RNA", layer = lyr)

    # Calculate stats
    qc_stats <- scuttle::perCellQCMetrics(
      counts_matrix,
      percent_top = c(20, 50, 100)
    )
    qc_list[[lyr]] <- as.data.frame(qc_stats)

    # Free memory immediately
    rm(counts_matrix, qc_stats)
    gc()
  }

  # Bind and assign efficiently to Seurat metadata
  all_qc_stats <- do.call(rbind, unname(qc_list))

  seurat_object$pct_counts_in_top_20_genes <- all_qc_stats$percent.top_20
  seurat_object$pct_counts_in_top_50_genes <- all_qc_stats$percent.top_50
  seurat_object$pct_counts_in_top_100_genes <- all_qc_stats$percent.top_100
  seurat_object$sum_counts <- all_qc_stats$sum
  seurat_object$n_genes_detected <- all_qc_stats$detected

  rm(all_qc_stats, qc_list)
  gc()

  # ==========================================================
  # 4. MALAT1 Thresholding
  # ==========================================================
  if (any(grepl("^Malat1", rownames(seurat_object), ignore.case = TRUE))) {
    has_norm <- "data" %in%
      slotNames(seurat_object[["RNA"]]) &&
      ncol(seurat_object[["RNA"]]@data) > 0

    if (has_norm) {
      message(
        "[calc_metrics] RNA normalized data already present; using it for MALAT1 threshold."
      )
    } else {
      message(
        "[calc_metrics] RNA normalized data missing; calling NormalizeData() before MALAT1 threshold."
      )
      seurat_object <- NormalizeData(seurat_object)
    }

    seurat_object <- Add_MALAT1_Threshold(
      object = seurat_object,
      species = "mouse",
      sample_col = "sample",
      save_plots = TRUE,
      save_plot_path = fig_dir,
      save_plot_name = "MALAT1_Threshold_Plots"
    )
  }

  return(seurat_object)
}


# ============================================================
# Remove deprecated genes ----
# ============================================================

# Remove genes with "DEPRECATED-" prefix in their names.
remove_deprecated_genes <- function(seurat_obj) {
  # Find deprecated genes
  deprecated_genes <- grep("^DEPRECATED-", rownames(seurat_obj), value = TRUE)
  
  # If deprecated genes exist, remove them
  if (length(deprecated_genes) > 0) {
    message(paste("Removing", length(deprecated_genes), "deprecated genes from the Seurat object"))
    # Keep all genes EXCEPT the deprecated ones
    seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), deprecated_genes))
    
    message(paste('Keeping', length(rownames(seurat_obj)),'genes'))
  } else {
    message("No deprecated genes found in this Seurat object")
  }
  
  return(seurat_obj)
}

# ============================================================
# Gene-level filtering ----
# ============================================================

# according to best practices, want to filter out genes that
# are not detected in at least 20 cells.
# https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

filter_genes <- function(seurat_obj, min_cells = 3L) {
  lyrs <- Layers(seurat_obj, assay = "RNA", search = "counts")
  gene_ncells <- Reduce(
    "+",
    lapply(lyrs, function(lyr) {
      m <- LayerData(seurat_obj, assay = "RNA", layer = lyr)
      tabulate(m@i + 1L, nbins = nrow(m))
    })
  )
  keep <- gene_ncells >= min_cells
  message(sprintf(
    "[filter_genes] %d genes before | %d after | %d removed (min_cells = %d)",
    nrow(seurat_obj),
    sum(keep),
    sum(!keep),
    min_cells
  ))
  seurat_obj[keep, ]
}


# ============================================================
# PC elbow selection ----
# ============================================================

#' Determine the number of PCs to retain via two elbow heuristics
#'
#' Computes two cutoffs on the PCA stdev spectrum and returns the larger so
#' that both criteria are satisfied:
#'   - \code{co1}: first PC where cumulative variance >= 90% AND per-PC
#'     variance <= 5% (the "variance-threshold" PC)
#'   - \code{co2}: one past the last PC where the drop in per-PC variance
#'     to the next PC exceeds 0.1% (the "elbow" PC)
#'
#' @param seurat_obj Seurat object with a PCA reduction already computed.
#'
#' @return A list with \code{pct}, \code{cumu}, \code{co1}, \code{co2}, and
#'   \code{pcs} (the recommended number of PCs, = max(co1, co2)).
calc_pcs_elbow <- function(seurat_obj) {
  stdev <- seurat_obj[["pca"]]@stdev
  pct <- stdev / sum(stdev) * 100
  cumu <- cumsum(pct)

  # cumulative >= 90% and this PC contributes <= 5%
  co1 <- which(cumu >= 90 & pct <= 5)[1]

  # last PC where the drop in per-PC variance to the next PC exceeds 0.1%
  diffs <- pct[1:(length(pct) - 1)] - pct[2:length(pct)]
  co2 <- sort(which(diffs > 0.1), decreasing = TRUE)[1] + 1

  # take the larger so both heuristics are satisfied
  pcs <- max(co1, co2)

  message(sprintf(
    "[calc_pcs_elbow] co1 (cumu>=90 & pct<=5): %d | co2 (last d(pct)>0.1): %d | pcs = max(co1, co2): %d",
    co1,
    co2,
    pcs
  ))

  list(pct = pct, cumu = cumu, co1 = co1, co2 = co2, pcs = pcs)
}


# ============================================================
# Doublet detection ----
# ============================================================

#' Two-round scDblFinder doublet detection
#'
#' Runs scDblFinder once on the full object, then re-runs it on the
#' round-1 singlets. Writes round-specific score/class columns and a
#' combined \code{doublet_final} label to the Seurat metadata.
#' This two-round approach following She et al. 2025
#'
#' @param seurat_obj  Seurat object (post QC, gene-filtered)
#' @param bpparam     BiocParallelParam used for scDblFinder's internal
#'                    PCA / KNN steps only.
#' @return the input Seurat object with new metadata columns:
#'   scDblFinder.score_r1, scDblFinder.class_r1,
#'   scDblFinder.score_r2, scDblFinder.class_r2,
#'   doublet_final ("singlet" | "doublet"),
#'   doublet_round ("singlet" | "round1_doublet" | "round2_doublet")
run_doublet_detection <- function(seurat_obj, bpparam, samples = NULL) {
  run_once <- function(so, label) {
    message(sprintf("[doublet] round %s: %d cells", label, ncol(so)))
    sce <- as.SingleCellExperiment(so)
    # Ensure a clean backend for this call — avoids stale worker sockets
    # from the previous round exhausting file descriptors.
    if (BiocParallel::bpisup(bpparam)) {
      BiocParallel::bpstop(bpparam)
    }
    BiocParallel::bpstart(bpparam)
    on.exit(
      if (BiocParallel::bpisup(bpparam)) BiocParallel::bpstop(bpparam),
      add = TRUE
    )
    sce <- scDblFinder(
      sce,
      clusters = TRUE,
      samples = samples,
      dbr.sd = 1,
      BPPARAM = bpparam
    )
    out <- data.frame(
      score = sce$scDblFinder.score,
      class = as.character(sce$scDblFinder.class),
      difficulty = sce$scDblFinder.difficulty,
      mostLikelyOrigin = as.character(sce$scDblFinder.mostLikelyOrigin),
      originAmbiguous = sce$scDblFinder.originAmbiguous,
      row.names = colnames(sce),
      stringsAsFactors = FALSE
    )
    rm(sce)
    gc()
    out
  }

  # round 1 — all cells
  r1 <- run_once(seurat_obj, "1")
  seurat_obj$scDblFinder.score_r1 <- r1[colnames(seurat_obj), "score"]
  seurat_obj$scDblFinder.class_r1 <- r1[colnames(seurat_obj), "class"]
  seurat_obj$scDblFinder.difficulty_r1 <- r1[colnames(seurat_obj), "difficulty"]
  seurat_obj$scDblfinder.mostLikelyOrigin_r1 <- r1[colnames(seurat_obj), "mostLikelyOrigin"]
  seurat_obj$scDblFinder.originAmbiguous_r1 <- r1[colnames(seurat_obj), "originAmbiguous"]

  # round 2 — only round-1 singlets
  singlet_cells <- rownames(r1)[r1$class == "singlet"]
  r2 <- run_once(seurat_obj[, singlet_cells], "2")

  seurat_obj$scDblFinder.score_r2 <- NA_real_
  seurat_obj$scDblFinder.class_r2 <- NA_character_
  seurat_obj$scDblFinder.difficulty_r2 <- NA_real_
  seurat_obj$scDblFinder.mostLikelyOrigin_r2 <- NA_character_
  seurat_obj$scDblFinder.originAmbiguous_r2 <- NA
  r2_rows <- match(rownames(r2), colnames(seurat_obj))
  seurat_obj$scDblFinder.score_r2[r2_rows] <- r2$score
  seurat_obj$scDblFinder.class_r2[r2_rows] <- r2$class
  seurat_obj$scDblFinder.difficulty_r2[r2_rows] <- r2$difficulty
  seurat_obj$scDblFinder.mostLikelyOrigin_r2[r2_rows] <- r2$mostLikelyOrigin
  seurat_obj$scDblFinder.originAmbiguous_r2[r2_rows] <- r2$originAmbiguous

  # final label: doublet in either round
  is_doublet <- seurat_obj$scDblFinder.class_r1 == "doublet" |
    (!is.na(seurat_obj$scDblFinder.class_r2) &
      seurat_obj$scDblFinder.class_r2 == "doublet")
  seurat_obj$doublet_final <- factor(
    ifelse(is_doublet, "doublet", "singlet"),
    levels = c("singlet", "doublet")
  )

  # per-round label: which round flagged the cell
  seurat_obj$doublet_round <- dplyr::case_when(
    seurat_obj$scDblFinder.class_r1 == "doublet" ~ "round1_doublet",
    seurat_obj$scDblFinder.class_r2 == "doublet" ~ "round2_doublet",
    TRUE ~ "singlet"
  )

  message(sprintf(
    "[doublet] summary: singlet=%d | doublet=%d (R1=%d, R2=%d)",
    sum(!is_doublet),
    sum(is_doublet),
    sum(seurat_obj$scDblFinder.class_r1 == "doublet"),
    sum(
      seurat_obj$scDblFinder.class_r1 == "singlet" &
        !is.na(seurat_obj$scDblFinder.class_r2) &
        seurat_obj$scDblFinder.class_r2 == "doublet"
    )
  ))

  seurat_obj
}


# ============================================================
# Subclustering on a celltype ----
# ============================================================

#' Subset a celltype, renormalize with SCT, integrate with Harmony,
#' and Leiden-cluster on the harmony graph.
#'
#' Subsets `seurat_obj` to cells whose `celltype_col` value is in
#' `celltype_values`, drops existing reductions / non-RNA assays, then
#' runs the SCT -> PCA -> Harmony -> Leiden -> UMAP pipeline used
#' elsewhere in this project.
#'
#' `sct_mode`:
#'   * "per_sample" — split by `split_by`, SCTransform each split,
#'     merge, then carry shared integration features. Canonical Seurat
#'     SCT-integration workflow. Preferred when each split has enough
#'     cells (~>100).
#'   * "merged"     — single SCTransform on the full subset. Faster and
#'     more robust when some splits have very few cells of the target
#'     celltype; relies entirely on Harmony for batch correction.
#'
#' @param seurat_obj        Seurat object with `celltype_col`, `split_by`,
#'                          and `batch_var` columns present, and an RNA
#'                          assay with raw counts.
#' @param celltype_col      Metadata column to filter on (e.g.
#'                          `"major_celltype"`).
#' @param celltype_values   Character vector of values in `celltype_col`
#'                          to keep.
#' @param sct_mode          "per_sample" or "merged".
#' @param split_by          Metadata column used to split for per-sample
#'                          SCT (ignored if `sct_mode == "merged"`).
#' @param batch_var         Metadata column(s) passed to `RunHarmony` as
#'                          the batch variable(s).
#' @param vars_to_regress   Covariates regressed in SCTransform.
#' @param per_sample_nfeatures Per-sample SCT `variable.features.n`
#'                          (only used in per_sample mode).
#' @param nfeatures         Integration feature count (per_sample) or
#'                          merged SCT `variable.features.n` (merged).
#' @param npcs              PCs computed by `RunPCA`.
#' @param dims              Dims used in `RunHarmony` / `FindNeighbors` /
#'                          `RunUMAP`.
#' @param resolutions       Resolutions passed to `FindClusters`.
#' @param reduction_suffix  Suffix appended to harmony / UMAP / graph
#'                          names so subset reductions don't collide if
#'                          merged back (default: "").
#' @param seed              RNG seed.
#'
#' @return Subsetted Seurat object with new SCT assay; `pca`,
#'   `harmony<suffix>`, `umap_harmony<suffix>` reductions; and
#'   `harmony<suffix>_snn_res.*` Leiden cluster columns.
subcluster_by_celltype <- function(
  seurat_obj,
  celltype_col,
  celltype_values,
  sct_mode = c("per_sample", "merged"),
  split_by = "sample",
  batch_var = "batch",
  vars_to_regress = c("cc.difference", "percent.mt"),
  per_sample_nfeatures = 3000,
  nfeatures = 5000,
  npcs = 50,
  dims = 1:50,
  resolutions = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
  reduction_suffix = "",
  seed = 2024
) {
  sct_mode <- match.arg(sct_mode)
  set.seed(seed)

  meta_cols <- colnames(seurat_obj@meta.data)
  if (!celltype_col %in% meta_cols) {
    stop(sprintf(
      "[subcluster_by_celltype] `%s` not in metadata.",
      celltype_col
    ))
  }
  missing_batch <- setdiff(batch_var, meta_cols)
  if (length(missing_batch)) {
    stop(sprintf(
      "[subcluster_by_celltype] batch_var(s) missing from metadata: %s",
      paste(missing_batch, collapse = ", ")
    ))
  }

  # ---- subset and strip prior reductions / assays ---------------------
  keep <- as.character(seurat_obj@meta.data[[celltype_col]]) %in%
    as.character(celltype_values)
  if (!any(keep)) {
    stop(sprintf(
      "[subcluster_by_celltype] no cells matched %s in {%s}",
      celltype_col,
      paste(celltype_values, collapse = ", ")
    ))
  }
  message(sprintf(
    "[subcluster_by_celltype] keeping %d / %d cells where %s in {%s}",
    sum(keep),
    length(keep),
    celltype_col,
    paste(celltype_values, collapse = ", ")
  ))

  sub <- DietSeurat(seurat_obj, assays = "RNA")
  sub <- sub[, keep]

  # ---- SCT normalization ----------------------------------------------
  if (sct_mode == "per_sample") {
    if (!split_by %in% meta_cols) {
      stop(sprintf(
        "[subcluster_by_celltype] split_by column `%s` not in metadata.",
        split_by
      ))
    }
    n_splits <- length(unique(sub@meta.data[[split_by]]))
    message(sprintf(
      "[subcluster_by_celltype] SCT per_sample on %d splits of `%s`",
      n_splits,
      split_by
    ))

    split_obj <- SplitObject(sub, split.by = split_by)
    split_obj <- lapply(split_obj, function(x) {
      SCTransform(
        x,
        assay = "RNA",
        vst.flavor = "v2",
        variable.features.n = per_sample_nfeatures,
        vars.to.regress = vars_to_regress,
        return.only.var.genes = FALSE,
        verbose = FALSE
      )
    })

    features <- SelectIntegrationFeatures(
      split_obj,
      nfeatures = nfeatures,
      fvf.nfeatures = nfeatures
    )

    sub <- merge(
      split_obj[[1]],
      split_obj[-1],
      merge.data = TRUE,
      merge.dr = FALSE
    )
    VariableFeatures(sub) <- features
    # SCT scale.data is per-split after merge; reset variable features
    # to the rows present in scale.data so PCA uses them all.
    # https://github.com/satijalab/seurat/issues/2814
    VariableFeatures(sub[["SCT"]]) <- rownames(sub[["SCT"]]@scale.data)
    rm(split_obj)
    gc()
  } else {
    message("[subcluster_by_celltype] SCT merged (single pass)")
    sub <- SCTransform(
      sub,
      assay = "RNA",
      vst.flavor = "v2",
      variable.features.n = nfeatures,
      vars.to.regress = vars_to_regress,
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
  }

  # ---- PCA / Harmony / Leiden / UMAP ----------------------------------
  DefaultAssay(sub) <- "SCT"

  harmony_name <- paste0("harmony", reduction_suffix)
  umap_name <- paste0("umap_harmony", reduction_suffix)
  nn_name <- paste0(harmony_name, "_nn")
  snn_name <- paste0(harmony_name, "_snn")

  sub <- RunPCA(sub, assay = "SCT", npcs = npcs, verbose = FALSE)

  sub <- RunHarmony(
    sub,
    group.by.vars = batch_var,
    reduction = "pca",
    reduction.save = harmony_name,
    assay.use = "SCT",
    verbose = FALSE
  )

  sub <- FindNeighbors(
    sub,
    reduction = harmony_name,
    dims = dims,
    graph.name = c(nn_name, snn_name),
    verbose = FALSE
  )

  sub <- FindClusters(
    sub,
    resolution = resolutions,
    algorithm = 4, # Leiden
    graph.name = snn_name,
    random.seed = seed,
    verbose = FALSE
  )

  sub <- RunUMAP(
    sub,
    reduction = harmony_name,
    dims = dims,
    reduction.name = umap_name,
    verbose = FALSE
  )

  sub
}



# ============================================================
# Check PCs for cell cycle and mt genes ----
# ============================================================

#' Diagnose which PCs are driven by mitochondrial or cell-cycle genes
#'
#' Inspects the top contributing genes per PC (by signed loading, separately
#' for the positive and negative pole) and flags PCs whose top loadings
#' include mitochondrial genes (`^mt-` for mouse, `^MT-` for human) or
#' Seurat's `cc.genes` S / G2M phase sets. Useful for deciding which PCs to
#' drop, or which covariates to regress out before re-running PCA.
#'
#' For mouse, the human-uppercase `cc.genes` lists are converted with
#' `stringr::str_to_title()` to match MGI symbol case before searching.
#'
#' @param seurat_obj  Seurat object with a PCA reduction already computed.
#' @param species     "mouse" or "human". Controls MT regex case and whether
#'                    `cc.genes` are title-cased before matching.
#' @param reduction   Name of the PCA reduction to inspect. Default "pca".
#' @param top_n       Number of top genes per pole (positive / negative) to
#'                    examine per PC. Default 30.
#' @param pcs         Optional integer vector of PCs to scan. Default: all
#'                    PCs present in the reduction.
#' @param verbose     If TRUE (default), prints a per-PC summary with
#'                    `[MT]` / `[S]` / `[G2M]` tags next to flagged genes.
#'
#' @return A list with:
#'   - `flagged_pcs`: integer vector of PCs containing any MT or CC gene
#'     in their top loadings, ordered by total hit count (desc).
#'   - `loadings`: tibble with cols `PC`, `direction` ("pos"/"neg"),
#'     `rank`, `gene`, `loading`, `category`
#'     ("mt" | "s_phase" | "g2m" | "other") — top `top_n` per pole for
#'     every scanned PC.
#'   - `summary`: tibble of flagged PCs only, with per-PC hit counts
#'     (`mt`, `s_phase`, `g2m`, `total_hits`), sorted by `PC` ascending.
check_pc_contamination <- function(
  seurat_obj,
  species = c("mouse", "human"),
  reduction = "pca",
  top_n = 30,
  pcs = NULL,
  verbose = TRUE
) {
  species <- match.arg(species)

  if (!reduction %in% Reductions(seurat_obj)) {
    stop(sprintf(
      "[check_pc_contamination] reduction `%s` not found.",
      reduction
    ))
  }

  # ---- gene sets ------------------------------------------------------
  mt_pattern <- if (species == "mouse") "^mt-" else "^MT-"
  if (species == "mouse") {
    s_set <- stringr::str_to_title(Seurat::cc.genes$s.genes)
    g2m_set <- stringr::str_to_title(Seurat::cc.genes$g2m.genes)
  } else {
    s_set <- Seurat::cc.genes$s.genes
    g2m_set <- Seurat::cc.genes$g2m.genes
  }

  # ---- loadings + one-shot gene classification ------------------------
  loadings_mat <- Loadings(seurat_obj, reduction = reduction)
  all_genes <- rownames(loadings_mat)
  # Classify the full feature space once; lookups per-PC are O(1).
  gene_cat <- dplyr::case_when(
    grepl(mt_pattern, all_genes) ~ "mt",
    all_genes %in% s_set ~ "s_phase",
    all_genes %in% g2m_set ~ "g2m",
    TRUE ~ "other"
  )

  available_pcs <- seq_len(ncol(loadings_mat))
  pcs <- if (is.null(pcs)) {
    available_pcs
  } else {
    intersect(as.integer(pcs), available_pcs)
  }
  if (length(pcs) == 0) {
    stop("[check_pc_contamination] no requested PCs are in the reduction.")
  }
  top_n <- min(top_n, nrow(loadings_mat))

  # One sort per PC; head/tail of ascending order gives both poles.
  per_pc <- lapply(pcs, function(i) {
    v <- loadings_mat[, i]
    ord <- order(v)
    neg_idx <- ord[seq_len(top_n)]
    pos_idx <- rev(ord[(length(ord) - top_n + 1L):length(ord)])
    idx <- c(pos_idx, neg_idx)
    tibble::tibble(
      PC = i,
      direction = factor(
        rep(c("pos", "neg"), each = top_n),
        levels = c("pos", "neg")
      ),
      rank = rep(seq_len(top_n), 2L),
      gene = all_genes[idx],
      loading = v[idx],
      category = unname(gene_cat[idx])
    )
  })
  loadings_tbl <- dplyr::bind_rows(per_pc)

  summary_tbl <- loadings_tbl |>
    dplyr::group_by(PC) |>
    dplyr::summarise(
      mt = sum(category == "mt"),
      s_phase = sum(category == "s_phase"),
      g2m = sum(category == "g2m"),
      total_hits = mt + s_phase + g2m,
      .groups = "drop"
    ) |>
    dplyr::filter(total_hits > 0) |>
    dplyr::arrange(PC)

  flagged <- summary_tbl$PC

  # ---- verbose printout -----------------------------------------------
  if (verbose) {
    if (length(flagged) == 0) {
      message(sprintf(
        "[check_pc_contamination] no MT or cell-cycle genes in top %d loadings of %d PC(s) scanned (species=%s).",
        top_n,
        length(pcs),
        species
      ))
    } else {
      message(sprintf(
        "[check_pc_contamination] %d / %d PC(s) flagged (top %d per pole | species=%s):",
        length(flagged),
        length(pcs),
        top_n,
        species
      ))
      tag_map <- c(mt = "[MT]", s_phase = "[S]", g2m = "[G2M]")
      pc_chunks <- split(loadings_tbl, loadings_tbl$PC)
      for (i in flagged) {
        ct <- summary_tbl[summary_tbl$PC == i, ]
        message(sprintf(
          "  PC%-3d  mt=%d  s_phase=%d  g2m=%d",
          i, ct$mt, ct$s_phase, ct$g2m
        ))
        pc_tbl <- pc_chunks[[as.character(i)]]

        # Flagged genes, listed separately from the loadings dump.
        flagged_genes <- pc_tbl[pc_tbl$category != "other", ]
        if (nrow(flagged_genes) > 0) {
          flagged_genes <- flagged_genes[
            order(-abs(flagged_genes$loading)),
          ]
          lines <- sprintf(
            "      %-6s %-15s (%s, loading=%+.3f)",
            tag_map[flagged_genes$category],
            flagged_genes$gene,
            flagged_genes$direction,
            flagged_genes$loading
          )
          message("    flagged genes:\n", paste(lines, collapse = "\n"))
        }

        # Clean top-N loadings without inline tags.
        pos_g <- pc_tbl$gene[pc_tbl$direction == "pos"]
        neg_g <- pc_tbl$gene[pc_tbl$direction == "neg"]
        message(sprintf(
          "    top loadings:\n      pos: %s\n      neg: %s",
          paste(pos_g, collapse = ", "),
          paste(neg_g, collapse = ", ")
        ))
      }
    }
  }

  list(
    flagged_pcs = flagged,
    loadings = loadings_tbl,
    summary = summary_tbl
  )
}

