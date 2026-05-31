library(Seurat)
library(tidyverse)
library(qs2)
library(sccomp)

# ---- point sccomp's Stan model cache into the scRNAseq project -------------
# sccomp defaults its compiled-model cache to ~/.sccomp_models (issue #218: the
# Bioconductor binary even froze this to the build server's home - fixed here
# by reinstalling with type = "source"). We relocate it into the scRNAseq
# project so the compiled Stan model lives with the rest of the analysis.
# Run this once per session, before any sccomp call.

local({
  ver   <- as.character(packageVersion("sccomp"))
  cache <- file.path("c:/Users/David/Documents/Research/PhD/scRNAseq",
                     ".sccomp_models", ver)
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  ns <- asNamespace("sccomp")
  if (bindingIsLocked("sccomp_stan_models_cache_dir", ns))
    unlockBinding("sccomp_stan_models_cache_dir", ns)
  assign("sccomp_stan_models_cache_dir", cache, envir = ns)
  message("sccomp model cache: ", cache)
})


seed = 2024
set.seed(seed)

nthreads = 2

analysis <- 'cellbender_analysis'
## directories for integrated data ----
fig_dir = paste0("c:/Users/David/Documents/Research/PhD/scRNAseq/analysis_old/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir = paste0("c:/Users/David/Documents/Research/PhD/scRNAseq/analysis_old/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir = paste0("c:/Users/David/Documents/Research/PhD/scRNAseq/analysis_old/integrated/",analysis,"/FPR_0.0_MAD/annotation/")
stats_dir = paste0("c:/Users/David/Documents/Research/PhD/scRNAseq/analysis_old/integrated/",analysis,"/FPR_0.0_MAD/stats/DA_analysis/")

sccomp_dir <- file.path(stats_dir, "sccomp")
voomCLR_dir <- file.path(stats_dir, "voomCLR")

dirs_to_make <- c(
  file.path(sccomp_dir, "DA_results"),
  file.path(sccomp_dir, "DA_plots"),
  file.path(voomCLR_dir, "DA_results"),
  file.path(voomCLR_dir, "DA_plots"),
  file.path(sccomp_dir, 'epithelial', "DA_results"),
  file.path(sccomp_dir, 'epithelial', "DA_plots"),
  file.path(voomCLR_dir, 'epithelial', "DA_results"),
  file.path(voomCLR_dir, 'epithelial', "DA_plots")
)

walk(dirs_to_make, ~dir.create(.x, recursive = TRUE, showWarnings = FALSE))


annotated_seurat = qs_read(paste0(data_dir,'de_seurat_annotated.qs'), nthreads=nthreads)
DefaultAssay(annotated_seurat) <- "RNA"
Idents(annotated_seurat) <- 'celltype'

epithelial = qs_read(paste0(data_dir,'epithelial.qs'), nthreads=nthreads)
unique(epithelial$subtype)

DimPlot(epithelial,
  reduction = 'umap_harmony',
  group.by = "subtype",
  split.by = "timepoint",
)

unique(annotated_seurat$celltype)

DimPlot(annotated_seurat,
  group.by = "celltype",
  split.by = "timepoint",
)

# port epithelial subtype labels into the full annotated object ----
# start from the celltype label for every cell in the annotated object
subcluster <- setNames(as.character(annotated_seurat$celltype), colnames(annotated_seurat))

# named vector of subtype labels from the epithelial object
epi_subtype <- setNames(as.character(epithelial$subtype), colnames(epithelial))

# cells shared by both objects (i.e. the epithelial cells)
shared <- intersect(names(subcluster), names(epi_subtype))

# sanity check: every epithelial cell should match a cell in annotated_seurat
length(shared) == ncol(epithelial)  # should be TRUE

# override the shared cells with their epithelial subtype label;
# all other cells keep their celltype label
subcluster[shared] <- epi_subtype[shared]

# add back to the annotated object (reordered to match the object's cells)
annotated_seurat$subtype <- subcluster[colnames(annotated_seurat)]

# check the result
table(annotated_seurat$subtype, useNA = "ifany")




# ============================================================
# Differential abundance with sccomp - per timepoint, pb vs ctrl ----
# ============================================================


# make sure the grouping variables are well-defined
annotated_seurat$condition <- factor(as.character(annotated_seurat$condition),
                                     levels = c("ctrl", "pb"))
annotated_seurat$subtype   <- as.character(annotated_seurat$subtype)

# collapse all epithelial subtypes (every label that came from the epithelial
# object) into a single "Epithelial" cell group for the sccomp test. The fine
# labels are preserved as `subtype_fine` for plotting / later analyses; the
# loop below keeps using `subtype` and is otherwise unchanged.
epi_labels <- unique(as.character(epithelial$subtype))
annotated_seurat$subtype_epi <- annotated_seurat$subtype
annotated_seurat$subtype_epi <- ifelse(annotated_seurat$subtype %in% epi_labels,
                                   "Epithelial", annotated_seurat$subtype)
cat("subtype after epithelial collapse:\n")
print(table(annotated_seurat$subtype_epi, useNA = "ifany"))

# cells with a missing subtype label cannot be modelled - drop them
keep <- !is.na(annotated_seurat$subtype_epi)
if (any(!keep)) {
  cat("Dropping", sum(!keep), "cells with NA subtype label\n")
  annotated_seurat <- annotated_seurat[, keep]
}

# ---- per-(timepoint, subtype) cell-count filter ---------------------------
# Drop (timepoint, subtype) combinations too sparse for either compositional
# method to estimate reliably: CLR becomes unstable (geometric mean polluted
# by near-zero entries; voom mean-variance trend rotated by high-variance
# rare features; BC mode density widened), and sccomp's prior cannot rescue
# a subtype that's effectively absent. Rule mirrors edgeR::filterByExpr's
# "at least one full group" logic: keep subtype within a timepoint iff
# >= min_cells in at least min(n_ctrl, n_pb) samples for that timepoint.
# Computed once here as `keep_rule` and applied inside each loop's timepoint
# subset (sccomp + voomCLR) - this drops the subtype from the *analysis* for
# the failing timepoint(s) without mutating annotated_seurat, so the object
# stays whole for the combined-timepoint model and any downstream uses.
# The filter is per-(timepoint, subtype), so a subtype rare at one timepoint
# but common at another stays in for the timepoints where it's detectable.

min_cells <- 10

subtype_counts <- annotated_seurat@meta.data |>
  dplyr::count(timepoint, condition, sample, subtype, name = "n") |>
  tidyr::complete(tidyr::nesting(timepoint, condition, sample), subtype,
                  fill = list(n = 0))

samples_per_tp <- subtype_counts |>
  dplyr::distinct(timepoint, condition, sample) |>
  dplyr::count(timepoint, condition, name = "n_samples") |>
  dplyr::group_by(timepoint) |>
  dplyr::summarise(n_min = min(n_samples), .groups = "drop")

keep_rule <- subtype_counts |>
  dplyr::mutate(passes = n >= min_cells) |>
  dplyr::group_by(timepoint, subtype) |>
  dplyr::summarise(
    n_passing = sum(passes),
    n_ctrl    = sum(passes & condition == "ctrl"),
    n_pb      = sum(passes & condition == "pb"),
    .groups   = "drop"
  ) |>
  dplyr::left_join(samples_per_tp, by = "timepoint") |>
  dplyr::mutate(keep = n_passing >= n_min)

cat(sprintf(
  "Subtype filter (>= %d cells in at least min(n_ctrl, n_pb) samples / timepoint):\n",
  min_cells))
print(keep_rule |> dplyr::filter(!keep) |> dplyr::arrange(timepoint, subtype),
      n = Inf)

readr::write_csv(keep_rule,
                 file.path(stats_dir, "subtype_filter_decisions.csv"))

timepoints     <- levels(annotated_seurat$timepoint)   # wk3, wk10, mn7, mn18
sccomp_results <- list()

# ---- helper: save every plot demoed at github.com/MangiolaLaboratory/sccomp -
# Uses sccomp's native plotting. Requires sccomp's plotting to be compatible
# with the installed ggplot2 (sccomp 2.x, or sccomp 1.10.0 with ggplot2 3.5.x);
# on sccomp 1.10.0 + ggplot2 >= 4.0 every plot here will report FAILED.
# plot() returns a named list of ggplots ($credible_intervals_1D / _2D); each
# plot is saved independently in tryCatch so one failure does not abort the rest.
save_sccomp_plots <- function(res, tp, out_dir) {

  save_one <- function(name, plt, w = 12, h = 10) {
    tryCatch({
      ggsave(file.path(out_dir, sprintf("sccomp_%s_%s.png", tp, name)),
             plot = plt, width = w, height = h, dpi = 320, bg = "white")
      cat("  saved  :", name, "\n")
    }, error = function(e)
      cat("  FAILED :", name, "-", conditionMessage(e), "\n"))
  }

  # plot() -> list with $credible_intervals_1D and $credible_intervals_2D
  plots <- tryCatch(plot(res), error = function(e) {
    cat("  plot() failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(plots)) {
    save_one("credible_intervals_1D", plots$credible_intervals_1D, 12, 14)
    save_one("credible_intervals_2D", plots$credible_intervals_2D, 12, 12)
  }

  # standalone interval plots + proportion boxplot. The plotting calls are
  # lazy args, so errors from them are caught inside save_one's tryCatch.
  save_one("1D_intervals", plot_1D_intervals(res),                   12, 14)
  save_one("2D_intervals", plot_2D_intervals(res),                   12, 12)
  save_one("boxplot",      sccomp_boxplot(res, factor = "condition"), 14, 11)
}

# sccomp loop ----
for (tp in timepoints) {

  cat("\n==================  sccomp:", tp, "==================\n")

  # subset to one timepoint via cell-name indexing (robust, no NSE);
  # additionally restrict to subtypes that passed the filter at this tp.
  # annotated_seurat itself is not mutated - we only narrow the cells fed
  # to the model so the count matrix omits the failing subtype rows.
  
  # cells_tp <- colnames(annotated_seurat)[annotated_seurat$timepoint == tp]
  # obj_tp   <- annotated_seurat[, cells_tp]
  # obj_tp$condition <- droplevels(obj_tp$condition)

  cells_tp <- colnames(epithelial)[epithelial$timepoint == tp]
  obj_tp   <- epithelial[, cells_tp]
  # obj_tp$condition <- droplevels(obj_tp$condition)

  dropped_at_tp <- keep_rule |>
    dplyr::filter(timepoint == tp, !keep) |> dplyr::pull(subtype)
  if (length(dropped_at_tp))
    cat("  excluded subtypes at", tp, ":",
        paste(dropped_at_tp, collapse = ", "), "\n")

  # sanity check: cells per subtype per condition (post-filter)
  cat("Cells per subtype x condition:\n")
  print(table(obj_tp$subtype, obj_tp$condition))

  # if (nlevels(obj_tp$condition) < 2) {
  #   cat("Skipping", tp, "- only one condition present\n")
  #   next
  # }

  res <- tryCatch({

    obj_tp |>
      sccomp_estimate(
        formula_composition = ~ condition,
        formula_variability = ~ condition, # default
        sample             = 'sample',
        cell_group         = 'subtype',
        bimodal_mean_variability_association = TRUE, # recommended for scRNAseq data
        inference_method = 'hmc', # default is "pathfinder" (fast, approximate), switched to "hmc" (slower, more accurate) for final results.
        cores               = nthreads,
        mcmc_seed         = seed,
        verbose             = FALSE
      ) |>
      # sccomp_remove_outliers(cores = nthreads) |>
      sccomp_test()

  }, error = function(e) {
    cat("sccomp FAILED for", tp, ":", conditionMessage(e), "\n")
    NULL
  })

  if (is.null(res)) next
  sccomp_results[[tp]] <- res

  # save full result object (keeps list-columns: count_data, draws, ...)
  # saveRDS(res, file.path(sccomp_dir, "DA_results", paste0("sccomp_", tp, ".rds")))
  saveRDS(res, file.path(sccomp_dir, 'epithelial', "DA_results", paste0("sccomp_", tp, ".rds")))

  # save a flat results table (drop list-columns so it writes to csv),
  # appending proportion_fold_change / average_uncertainty / statement
  # from sccomp_proportional_fold_change() right after c_ess_tail.
  fc_df <- tryCatch({
    pfc <- sccomp_proportional_fold_change(
      res,
      formula_composition = ~ condition,
      from = "ctrl",
      to   = "pb"
    )
    dplyr::select(pfc, subtype, proportion_fold_change,
                  average_uncertainty, statement)
  }, error = function(e) {
    cat("sccomp_proportional_fold_change failed for", tp, ":",
        conditionMessage(e), "\n")
    NULL
  })

  flat <- res |> dplyr::select(dplyr::where(~ !is.list(.x)))
  if (!is.null(fc_df))
    flat <- flat |>
      dplyr::left_join(fc_df, by = "subtype") |>
      dplyr::relocate(proportion_fold_change, average_uncertainty, statement,
                      .after = c_FDR)
  # readr::write_csv(flat, file.path(sccomp_dir, "DA_results",
  #                                  paste0("sccomp_", tp, "_pb_vs_ctrl.csv")))
  readr::write_csv(flat, file.path(sccomp_dir, 'epithelial',"DA_results",
                                   paste0("sccomp_", tp, "_pb_vs_ctrl.csv")))
  
  # ---- plots ----
  # (a) custom ggplot2 plots (forest + observed proportions)
  # (b) sccomp's native plots via save_sccomp_plots()

  # (a1) effect-size forest plot: pb-vs-ctrl composition effect per subtype
  tryCatch({
    fp_df <- res |>
      dplyr::filter(parameter == "conditionpb") |>
      dplyr::mutate(significant = c_FDR < 0.05,
                    subtype     = forcats::fct_reorder(subtype, c_effect))

    p_forest <- ggplot(fp_df, aes(x = c_effect, y = subtype, colour = significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_linerange(aes(xmin = c_lower, xmax = c_upper)) +
      geom_point(size = 2.5) +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                          name = "FDR < 0.05") +
      labs(title    = paste0("sccomp composition effect - ", tp),
           subtitle = "pb vs ctrl (ctrl = reference); bars = 95% credible interval",
           x = "Composition effect (log-fold change)", y = NULL) +
      theme_bw()

    # ggsave(file.path(sccomp_dir, "DA_plots", paste0("sccomp_", tp, "_forest.png")),
    #        plot = p_forest, width = 10, height = 8, dpi = 320, bg = "white")
    ggsave(file.path(sccomp_dir, 'epithelial', "DA_plots", paste0("sccomp_", tp, "_forest.png")),
           plot = p_forest, width = 10, height = 8, dpi = 320, bg = "white")
  }, error = function(e) cat("forest plot failed for", tp, ":",
                             conditionMessage(e), "\n"))

  # (a2) observed per-sample proportions per subtype, pb vs ctrl
  #      (file named _observed_proportions to avoid clashing with sccomp's
  #       own _boxplot.png written by save_sccomp_plots() below)
  tryCatch({
    prop_df <- obj_tp@meta.data |>
      dplyr::count(sample, condition, subtype, name = "n") |>
      dplyr::group_by(sample) |>
      dplyr::mutate(proportion = n / sum(n)) |>
      dplyr::ungroup()

    p_box <- ggplot(prop_df, aes(x = condition, y = proportion, fill = condition)) +
      geom_boxplot(outliers = FALSE, alpha = 0.6) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
      facet_wrap(~ subtype, scales = "free_y") +
      labs(title = paste0("Cell-type proportions per sample - ", tp),
           x = NULL, y = "Proportion of cells in sample") +
      ggpubr::theme_classic2()

    # ggsave(file.path(sccomp_dir, "DA_plots", paste0("sccomp_", tp, "_observed_proportions.png")),
    #        plot = p_box, width = 14, height = 11, dpi = 320, bg = "white")
    ggsave(file.path(sccomp_dir, 'epithelial', "DA_plots", paste0("sccomp_", tp, "_observed_proportions.png")),
           plot = p_box, width = 14, height = 11, dpi = 320, bg = "white")
  }, error = function(e) cat("observed proportions plot failed for", tp, ":",
                             conditionMessage(e), "\n"))

  # (b) sccomp's native plots (credible intervals 1D/2D, boxplot)
  cat("Saving sccomp native plots for", tp, ":\n")
  # save_sccomp_plots(res, tp, file.path(sccomp_dir, "DA_plots"))
  save_sccomp_plots(res, tp, file.path(sccomp_dir, 'epithelial', "DA_plots"))

  gc()  # free memory before next iteration
}


# ---- combined summary across timepoints ----
da_summary <- imap_dfr(sccomp_results, ~ .x |>
  dplyr::select(dplyr::where(~ !is.list(.x))) |>
  dplyr::mutate(timepoint = .y, .before = 1))

write_csv(da_summary,
          file.path(sccomp_dir, "DA_results", "sccomp_all_timepoints_pb_vs_ctrl.csv"))

# significant differentially abundant subtypes (condition effect, FDR < 0.05)
da_summary |>
  dplyr::filter(parameter != "(Intercept)", c_FDR < 0.05) |>
  dplyr::arrange(timepoint, c_FDR) |>
  dplyr::select(timepoint, subtype, parameter,
                c_effect, c_lower, c_upper, c_FDR) |>
  print(n = Inf)




# ============================================================
# Across-timepoint composition figure (combined sccomp model) ----
# ============================================================
# A single sccomp model over ALL timepoints, so the pb-vs-ctrl
# difference can be compared across time. We use the cell-means form
#   ~ 0 + group        (group = condition x timepoint)
# which is statistically the same model as  ~ condition * timepoint,
# but every pb-vs-ctrl comparison is then a clean two-term contrast
# with no ':' interaction names to escape (matches sccomp's documented
# contrast syntax, e.g. "typecancer - typebenign").
#
#   sccomp_predict()        -> estimated proportion + 95% CI per group
#   sccomp_test(contrasts=) -> pb-vs-ctrl test WITHIN each timepoint
# One publication figure per cell type -> DA_plots/across_timepoints/.

tp_order <- c("wk3", "wk10", "mn7", "mn18")

# condition x timepoint group factor (ctrl levels first, wk3 first)
annotated_seurat$timepoint <- factor(as.character(annotated_seurat$timepoint),
                                     levels = tp_order)
annotated_seurat$group <- factor(
  paste(annotated_seurat$condition, annotated_seurat$timepoint, sep = "_"),
  levels = paste(rep(c("ctrl", "pb"), each = length(tp_order)), tp_order, sep = "_")
)

# fit the combined model (one fit, all timepoints)
combined_estimate <- annotated_seurat |>
  sccomp_estimate(
    formula_composition = ~ 0 + group,
    formula_variability = ~ 1,
    .sample             = 'sample',
    .cell_group         = 'subtype',
    bimodal_mean_variability_association = TRUE,
    cores               = nthreads,
    mcmc_seed         = seed,
    verbose             = FALSE
  )

saveRDS(combined_estimate,
        file.path(sccomp_dir, "DA_results", "sccomp_combined_group_model.rds"))


# ---- pb-vs-ctrl test within each timepoint (contrasts) ----
tp_contrasts <- setNames(
  sprintf("grouppb_%s - groupctrl_%s", tp_order, tp_order),
  paste0("pb_vs_ctrl_", tp_order)
)
combined_test <- sccomp_test(combined_estimate, contrasts = tp_contrasts)
cat("combined_test columns:\n"); print(colnames(combined_test))

# defensively standardise the cell-group column name to "subtype"
if ("cell_group" %in% names(combined_test) && !"subtype" %in% names(combined_test))
  combined_test <- dplyr::rename(combined_test, subtype = cell_group)

# the contrast id lands in `parameter`; recover the timepoint from it
# (works whether parameter holds the name 'pb_vs_ctrl_wk10' or the
#  expression 'grouppb_wk10 - groupctrl_wk10' - both contain the tp)
sig_df <- combined_test |>
  dplyr::mutate(timepoint = stringr::str_extract(parameter,
                                                 "wk10|wk3|mn18|mn7")) |>
  dplyr::filter(!is.na(timepoint)) |>
  dplyr::transmute(
    subtype, c_effect, c_FDR,
    timepoint = factor(timepoint, levels = tp_order),
    sig_label = dplyr::case_when(c_FDR < 0.001 ~ "***",
                                 c_FDR < 0.01  ~ "**",
                                 c_FDR < 0.05  ~ "*",
                                 TRUE          ~ ""))

write_csv(sig_df, file.path(sccomp_dir, "DA_results",
                            "sccomp_combined_pb_vs_ctrl_per_timepoint.csv"))


# ---- model-estimated proportions (+ 95% CI) per group ----
pred_new_data <- tibble::tibble(
  sample = paste0("pred_", levels(annotated_seurat$group)),
  group  = factor(levels(annotated_seurat$group),
                  levels = levels(annotated_seurat$group))
)
combined_pred <- sccomp_predict(combined_estimate, new_data = pred_new_data)
cat("combined_pred columns:\n"); print(colnames(combined_pred))

if ("cell_group" %in% names(combined_pred) && !"subtype" %in% names(combined_pred))
  combined_pred <- dplyr::rename(combined_pred, subtype = cell_group)

# attach condition + timepoint
plot_df <- combined_pred |>
  dplyr::left_join(dplyr::select(pred_new_data, sample, group), by = "sample") |>
  tidyr::separate_wider_delim(group, delim = "_",
                              names = c("condition", "timepoint")) |>
  dplyr::mutate(condition = factor(condition, levels = c("ctrl", "pb")),
                timepoint = factor(timepoint, levels = tp_order))


# ---- one publication figure per cell type ----
fig_out <- file.path(sccomp_dir, "DA_plots", "across_timepoints")
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)

# y position for the significance star: just above the higher CI per timepoint
sig_pos <- plot_df |>
  dplyr::group_by(subtype, timepoint) |>
  dplyr::summarise(y = max(proportion_upper), .groups = "drop") |>
  dplyr::left_join(sig_df, by = c("subtype", "timepoint")) |>
  dplyr::filter(!is.na(sig_label), sig_label != "")

cond_cols <- c(ctrl = "#F8766D", pb = "#00BFC4")

for (st in sort(unique(plot_df$subtype))) {

  d <- dplyr::filter(plot_df, subtype == st)
  s <- dplyr::filter(sig_pos, subtype == st)

  p <- ggplot(d, aes(timepoint, proportion_mean,
                      colour = condition, group = condition)) +
    geom_line(linewidth = 1) +
    geom_errorbar(aes(ymin = proportion_lower, ymax = proportion_upper),
                  width = 0.12, linewidth = 0.7) +
    geom_point(size = 3) +
    geom_text(data = s, aes(timepoint, y, label = sig_label),
              inherit.aes = FALSE, vjust = -0.3, size = 7, colour = "grey15") +
    scale_colour_manual(values = cond_cols, name = "Condition") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(title    = paste0("Proportion of ", st, " across timepoints"),
         subtitle = "sccomp-estimated proportion; bars = 95% credible interval;  * FDR<0.05  ** <0.01  *** <0.001",
         x = "Timepoint", y = "Proportion of cells") +
    theme_minimal(base_size = 14) +
    theme(panel.grid      = element_blank(),       # no background gridlines
          axis.line       = element_line(colour = "grey30"),
          axis.ticks      = element_line(colour = "grey30"),
          plot.title      = element_text(face = "bold"),
          plot.subtitle   = element_text(colour = "grey35", size = 9),
          legend.position = "top")

  fname <- paste0("proportion_across_timepoints_",
                  gsub("[^A-Za-z0-9]+", "_", st), ".png")
  ggsave(file.path(fig_out, fname), plot = p,
         width = 7, height = 5, dpi = 320, bg = "white")
}
cat("Saved", length(unique(plot_df$subtype)), "figures to:", fig_out, "\n")




# ============================================================
#  voomCLR - per timepoint, pb vs ctrl ----
# ============================================================
# Method: voomCLR (Van den Berge, Bioinformatics 2025; doi 10.1093/bioinformatics/btaf637).
# Pipeline per timepoint, mirroring the sccomp loop above:
#   counts (subtypes x samples) -> voomCLR()  (CLR transform with +0.5,
#       observation-level weights via lowess on the CLR-transformed counts;
#       varCalc = "empirical" is the default and is fine here since we have
#       ~10-20 subtypes per timepoint - the paper recommends "analytical"
#       with varDistribution = "NB" only when subtypes are few)
#   -> lmFit -> eBayes -> topTableBC(..., bootstrap = "parametric",
#       voomWeights = v$weights, confint = TRUE)
# topTableBC (not standard topTable) applies the paper's bias correction:
# regression coefficients on CLR are biased by the mean shift alpha_bar;
# topTableBC subtracts the mode across cell types (the "most cell types are
# not DA" assumption) and the parametric bootstrap accounts for the variance
# of that bias estimate - critical for FDR control with few cell types and
# small n (paper Figs 3-4). plotBeta() saves the per-coefficient density
# with the bias mode marked as a diagnostic for that assumption.
# Outputs per timepoint: rds of the fit, a flat topTableBC CSV (logFC, CI.L,
# CI.R, P.Value, adj.P.Val), a forest + volcano + observed proportions plot,
# and a plotBeta diagnostic. A combined summary CSV is written at the end.

library(limma)
library(voomCLR)

voomCLR_results <- list()

for (tp in timepoints) {

  cat("\n==================  voomCLR:", tp, "==================\n")

  # restrict to (timepoint, subtype) combos that passed the filter; this
  # drops the failing subtype rows from the count matrix without touching
  # annotated_seurat.
  cells_tp <- colnames(annotated_seurat)[annotated_seurat$timepoint == tp]
  obj_tp   <- annotated_seurat[, cells_tp]
  obj_tp$condition <- droplevels(obj_tp$condition)

  # (subtypes x samples) count matrix
  cnt_mat <- obj_tp@meta.data |>
    dplyr::count(subtype, sample, name = "n") |>
    tidyr::pivot_wider(names_from = sample, values_from = n, values_fill = 0) |>
    tibble::column_to_rownames("subtype") |>
    as.matrix()

  # one row per sample, in the same order as the count matrix columns
  smp_meta <- obj_tp@meta.data |>
    dplyr::distinct(sample, condition) |>
    dplyr::mutate(sample = as.character(sample)) |>
    (\(d) d[match(colnames(cnt_mat), d$sample), ])()
  stopifnot(identical(smp_meta$sample, colnames(cnt_mat)))
  smp_meta$condition <- factor(as.character(smp_meta$condition),
                               levels = c("ctrl", "pb"))

  cat("Cells per subtype x condition:\n")
  print(table(obj_tp$subtype, obj_tp$condition))

  design   <- model.matrix(~ condition, data = smp_meta)
  coef_nm  <- grep("^condition", colnames(design), value = TRUE)

  res <- tryCatch({
    v   <- voomCLR(counts = cnt_mat, 
                   design = design, plot = TRUE, 
                   varCalc = "analytical", 
                   varDistribution = "NB")
    fit <- eBayes(lmFit(v, design))
    list(v = v, fit = fit, counts = cnt_mat, meta = smp_meta, coef = coef_nm)
  }, error = function(e) {
    cat("voomCLR FAILED for", tp, ":", conditionMessage(e), "\n"); NULL
  })

  if (is.null(res)) next
  voomCLR_results[[tp]] <- res

  saveRDS(res, file.path(voomCLR_dir, "DA_results", paste0("voomCLR_", tp, ".rds")))

  # flat topTableBC: bias-corrected logFC + CI + p-values for pb vs ctrl.
  # bootstrap = "parametric" accounts for the variance of the bias estimate;
  # voomWeights = v$weights feeds the parametric bootstrap (required for it).
  tt <- topTableBC(res$fit, coef = coef_nm, number = Inf,
                   sort.by = "p", confint = TRUE,
                   bootstrap = "parametric",
                   voomWeights = res$v$weights) |>
    tibble::rownames_to_column("subtype")
  readr::write_csv(tt, file.path(voomCLR_dir, "DA_results",
                                 paste0("voomCLR_", tp, "_pb_vs_ctrl.csv")))

  # plotBeta diagnostic: density of coefficients across subtypes with the
  # mode (the bias correction subtracts this) marked as an orange line.
  # Base R graphics -> wrap in a png device.
  tryCatch({
    png(file.path(voomCLR_dir, "DA_plots",
                  paste0("voomCLR_", tp, "_plotBeta.png")),
        width = 1600, height = 1200, res = 200)
    plotBeta(res$fit)
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat("plotBeta failed for", tp, ":", conditionMessage(e), "\n")
  })

  # (a) effect-size forest plot
  tryCatch({
    fp_df <- tt |>
      dplyr::mutate(significant = adj.P.Val < 0.05,
                    subtype     = forcats::fct_reorder(subtype, logFC))

    p_forest <- ggplot(fp_df, aes(x = logFC, y = subtype, colour = significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_linerange(aes(xmin = CI.L, xmax = CI.R)) +
      geom_point(size = 2.5) +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                          name = "FDR < 0.05") +
      labs(title    = paste0("voomCLR composition effect - ", tp),
           subtitle = "pb vs ctrl (ctrl = reference); bars = 95% confidence interval",
           x = "Composition effect (logFC, CLR scale)", y = NULL) +
      theme_bw()

    ggsave(file.path(voomCLR_dir, "DA_plots", paste0("voomCLR_", tp, "_forest.png")),
           plot = p_forest, width = 10, height = 8, dpi = 320, bg = "white")
  }, error = function(e) cat("forest plot failed for", tp, ":",
                             conditionMessage(e), "\n"))

  # (b) observed per-sample proportions per subtype, pb vs ctrl
  tryCatch({
    prop_df <- obj_tp@meta.data |>
      dplyr::count(sample, condition, subtype, name = "n") |>
      dplyr::group_by(sample) |>
      dplyr::mutate(proportion = n / sum(n)) |>
      dplyr::ungroup()

    p_box <- ggplot(prop_df, aes(x = condition, y = proportion, fill = condition)) +
      geom_boxplot(outliers = FALSE, alpha = 0.6) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
      facet_wrap(~ subtype, scales = "free_y") +
      labs(title = paste0("Cell-type proportions per sample - ", tp),
           x = NULL, y = "Proportion of cells in sample") +
      ggpubr::theme_classic2()

    ggsave(file.path(voomCLR_dir, "DA_plots", paste0("voomCLR_", tp, "_observed_proportions.png")),
           plot = p_box, width = 14, height = 11, dpi = 320, bg = "white")
  }, error = function(e) cat("observed proportions plot failed for", tp, ":",
                             conditionMessage(e), "\n"))

  # (c) volcano plot
  tryCatch({
    v_df <- tt |>
      dplyr::mutate(significant = adj.P.Val < 0.05,
                    nlp         = -log10(P.Value))

    p_volc <- ggplot(v_df, aes(logFC, nlp, colour = significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_point(size = 2) +
      ggrepel::geom_text_repel(
        data = dplyr::filter(v_df, significant),
        aes(label = subtype), size = 3, max.overlaps = Inf, show.legend = FALSE) +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                          name = "FDR < 0.05") +
      labs(title    = paste0("voomCLR volcano - ", tp),
           subtitle = "pb vs ctrl (ctrl = reference)",
           x = "logFC (CLR scale)", y = expression(-log[10]~P)) +
      theme_bw()

    ggsave(file.path(voomCLR_dir, "DA_plots", paste0("voomCLR_", tp, "_volcano.png")),
           plot = p_volc, width = 9, height = 7, dpi = 320, bg = "white")
  }, error = function(e) cat("volcano plot failed for", tp, ":",
                             conditionMessage(e), "\n"))

  gc()
}


# ---- combined summary across timepoints ----
voomCLR_summary <- imap_dfr(voomCLR_results, function(r, tp) {
  topTableBC(r$fit, coef = r$coef, number = Inf,
             sort.by = "none", confint = TRUE,
             bootstrap = "parametric",
             voomWeights = r$v$weights) |>
    tibble::rownames_to_column("subtype") |>
    dplyr::mutate(timepoint = tp, .before = 1)
})

write_csv(voomCLR_summary,
          file.path(voomCLR_dir, "DA_results",
                    "voomCLR_all_timepoints_pb_vs_ctrl.csv"))

# significant differentially abundant subtypes (FDR < 0.05)
voomCLR_summary |>
  dplyr::filter(adj.P.Val < 0.05) |>
  dplyr::arrange(timepoint, adj.P.Val) |>
  dplyr::select(timepoint, subtype, logFC, CI.L, CI.R, P.Value, adj.P.Val) |>
  print(n = Inf)
