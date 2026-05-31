# fGSEA_utils.R
# Shared helpers for fGSEA analysis. Sourced by fGSEA.R (interactive driver)
# and fgsea_array_task.R (Slurm per-task script).

library(tidyverse)
library(fgsea)
library(ggplot2)
library(viridis)
library(scales)

# ============================================================
# Shared colorbar config ----
# ============================================================

# Fixed log10 ticks across the meaningful FDR range. Pure decades give
# perceptually-uniform color steps under transform = "log10". The 5e-3 and
# 5e-2 ticks add the conventional biological thresholds; they sit closer to
# their decade neighbors, which is the correct log-space placement.
# Values below GSEA_PADJ_LIMITS[1] are clamped to the deepest color via
# oob = scales::squish at the call site (rather than being dropped or getting
# their own arbitrary tick label like "6.144e-6").
GSEA_PADJ_BREAKS <- c(1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
GSEA_PADJ_LABELS <- c("1e-4", "0.001", "0.005", "0.01", "0.05", "0.1")
GSEA_PADJ_LIMITS <- c(1e-4, 1e-1)

# ============================================================
# Annotation join + ranked gene list ----
# ============================================================

get_ids_and_homologs <- function(res_tbl, gene_ids = gene_id) {
  merge(x = res_tbl, y = gene_ids,
        by.x = "gene", by.y = "external_gene_name")
}

to_genelist <- function(res_tbl, species_human = FALSE) {
  if (!species_human) {
    res_tbl %>%
      as_tibble() %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
      dplyr::select(gene, ranking_metric) %>%
      na.omit() %>%
      arrange(desc(ranking_metric)) %>%
      distinct() %>%
      deframe()
  } else {
    res_tbl %>%
      as_tibble() %>%
      group_by(hsapiens_homolog_associated_gene_name) %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
      summarize(ranking_metric = max(sign(log2FoldChange) * -log10(pvalue),
                                     na.rm = TRUE)) %>%
      dplyr::select(hsapiens_homolog_associated_gene_name, ranking_metric) %>%
      na.omit() %>%
      arrange(desc(ranking_metric)) %>%
      distinct() %>%
      deframe()
  }
}

# ============================================================
# Tidy + save fGSEA results to CSV ----
# ============================================================

# Output filename: <celltype>_<timepoint>_<pathway_name>_fgseaRes.csv when
# timepoint is supplied; falls back to the legacy multilevel/simple naming
# otherwise so older calls without timepoint still work.
tidysave_fgsea <- function(fgseaRes, pathways, pathway_name, celltype,
                           timepoint = NULL, simple_tf = FALSE,
                           directory = directory_results) {
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)

  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    mutate(leadingEdge2 = leadingEdge) %>%
    rowwise() %>%
    mutate_at(c("leadingEdge2"), ~ paste(unlist(.), collapse = ",")) %>%
    dplyr::select(!leadingEdge) %>%
    rename(leadingEdge = leadingEdge2)

  if (!is.null(timepoint)) {
    out_csv <- paste0(directory, celltype, "_", timepoint, "_", pathway_name,
                      "_fgseaRes.csv")
  } else {
    out_csv <- paste0(directory, celltype,
                      ifelse(simple_tf, "_simple", "_multilevel"),
                      "_fGSEA_NES_", pathway_name, ".csv")
  }

  write.csv(x = fgseaResTidy, file = out_csv, row.names = FALSE)
  fgseaResTidy
}

# ============================================================
# Wrap long pathway names ----
# ============================================================

wrap_pathway_names <- function(pathways, width = 50) {
  vapply(pathways, function(p) {
    if (nchar(p) >= 60) {
      stringr::str_wrap(gsub("_", " ", p), width)
    } else {
      p
    }
  }, character(1))
}

# ============================================================
# NES horizontal barplot ----
# ============================================================

plot_fgsea_nes <- function(fgsea_nes, celltype, pathways, pathway_name,
                           timepoint = NULL, padj_cutoff = 0.25,
                           simple_tf = FALSE, directory = directory_figs) {

  if (nrow(fgsea_nes) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    return(ggplot() +
             geom_text(aes(x = 0, y = 0,
                           label = paste("No significant pathways at padj <=",
                                         padj_cutoff))) +
             theme_minimal())
  }

  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)

  filtered_data <- fgsea_nes %>%
    filter(padj <= padj_cutoff) %>%
    mutate(
      pathway_wrapped = wrap_pathway_names(pathway, width = 50),
      marker  = ifelse(padj <= 0.05, "*", ""),
      marker2 = ifelse((padj <= 0.1) & (padj > 0.05), "#", ""),
      marker3 = ifelse(padj <= 0.005, "**", "")
    ) %>%
    filter(!is.na(NES), !is.na(pathway_wrapped))

  if (nrow(filtered_data) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    return(ggplot() +
             geom_text(aes(x = 0, y = 0,
                           label = paste("No significant pathways at padj <=",
                                         padj_cutoff))) +
             theme_minimal())
  }

  n_pathways <- nrow(filtered_data)
  nes_values <- filtered_data$NES[is.finite(filtered_data$NES)]
  if (length(nes_values) > 0) {
    x_min <- floor(min(nes_values) - 0.5)
    x_max <- ceiling(max(nes_values) + 0.5)
  } else {
    x_min <- -3; x_max <- 3
  }
  if (x_min == x_max) { x_min <- x_min - 1; x_max <- x_max + 1 }
  step <- if (x_max - x_min > 20) ceiling((x_max - x_min) / 10) else 1

  p <- ggplot(filtered_data,
              aes(x = NES, y = reorder(pathway_wrapped, NES))) +
    geom_col(aes(fill = padj)) +
    labs(x = "Normalized Enrichment Score", y = "Pathway",
         title = paste0("GSEA - ", pathway_name, " Pathways: ", celltype),
         subtitle = "** denotes padj < 0.005, * denotes padj < 0.05, # denotes padj < 0.1",
         fill = paste0("FDR up to ", padj_cutoff)) +
    theme_light() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(color = "grey95"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey95"),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      text = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text.y = element_text(size = 12, margin = margin(r = 5)),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(size = 16, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, margin = margin(b = 5)),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      legend.key.height = unit(0.8, "lines"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
    ) +
    scale_fill_viridis(
      option = "viridis", direction = -1,
      transform = "log10",
      limits = GSEA_PADJ_LIMITS,
      breaks = GSEA_PADJ_BREAKS,
      labels = GSEA_PADJ_LABELS,
      oob    = scales::squish,
      na.value = "grey90",
      guide = guide_colorbar(barwidth = 1, barheight = 10, ticks.linewidth = 1)
    ) +
    geom_text(aes(label = marker,
                  x = ifelse(NES > 0, NES + 0.25, NES - 0.4)),
              hjust = 0, vjust = 0.5, size = 6, color = "black") +
    geom_text(aes(label = marker2,
                  x = ifelse(NES > 0, NES + 0.25, NES - 0.4)),
              hjust = 0, vjust = 0.5, size = 4, color = "black") +
    geom_text(aes(label = marker3,
                  x = ifelse(NES > 0, NES + 0.25, NES - 0.4)),
              hjust = 0, vjust = 0.5, size = 6, color = "red")

  width  <- max(12, 5 + (x_max - x_min) * 0.5)
  height <- max(8, min(25, 4 + n_pathways * 0.7))
  if (max(nchar(filtered_data$pathway)) > 30) {
    p <- p + theme(axis.text.y = element_text(margin = margin(r = 15)))
    width <- width + 2
  }

  tp_tag <- if (!is.null(timepoint)) paste0("_", timepoint) else ""
  ggsave(paste0(directory, celltype, tp_tag,
                ifelse(simple_tf, "_simple", "_multilevel"),
                "_fGSEA_NES_", pathway_name, "_hbarplot.png"),
         plot = p, dpi = 320, width = width, height = height)

  p
}

# ============================================================
# NES dot plot ----
# ============================================================

plot_fgsea_nes_dotplot <- function(fgsea_nes, celltype, pathways, pathway_name,
                                   timepoint = NULL, padj_cutoff = 0.25,
                                   simple_tf = FALSE, directory = directory_figs) {

  if (nrow(fgsea_nes) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    return(ggplot() +
             geom_text(aes(x = 0, y = 0,
                           label = paste("No significant pathways at padj <=",
                                         padj_cutoff))) +
             theme_minimal())
  }

  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)

  filtered_data <- fgsea_nes %>%
    filter(padj <= padj_cutoff) %>%
    mutate(
      pathway_wrapped = wrap_pathway_names(pathway, width = 50),
      marker  = ifelse(padj <= 0.05, "*", ""),
      marker2 = ifelse((padj <= 0.1) & (padj > 0.05), "#", ""),
      marker3 = ifelse(padj <= 0.005, "**", "")
    ) %>%
    filter(!is.na(NES), !is.na(pathway_wrapped))

  if (nrow(filtered_data) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    return(ggplot() +
             geom_text(aes(x = 0, y = 0,
                           label = paste("No significant pathways at padj <=",
                                         padj_cutoff))) +
             theme_minimal())
  }

  n_pathways <- nrow(filtered_data)
  nes_values <- filtered_data$NES[is.finite(filtered_data$NES)]
  if (length(nes_values) > 0) {
    x_min <- floor(min(nes_values) - 0.5)
    x_max <- ceiling(max(nes_values) + 0.5)
  } else {
    x_min <- -3; x_max <- 3
  }
  if (x_min == x_max) { x_min <- x_min - 1; x_max <- x_max + 1 }
  step <- if (x_max - x_min > 20) ceiling((x_max - x_min) / 10) else 1

  p <- ggplot(filtered_data,
              aes(x = NES, y = reorder(pathway_wrapped, NES),
                  size = size, color = padj)) +
    geom_point() +
    labs(x = "Normalized Enrichment Score", y = "Pathway",
         title = paste0("GSEA - ", pathway_name, " Pathways: ", celltype),
         subtitle = "** denotes padj < 0.005, * denotes padj < 0.05, # denotes padj < 0.1",
         size = "Geneset Size",
         color = paste0("FDR up to ", padj_cutoff)) +
    theme_light() +
    theme(
      axis.line = element_line(colour = "black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_line(color = "grey95"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey95"),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      legend.box = "vertical",
      text = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text.y = element_text(size = 12, margin = margin(r = 5)),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(size = 16, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, margin = margin(b = 5)),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      legend.key.height = unit(0.8, "lines"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
    ) +
    scale_color_viridis(
      option = "viridis", direction = -1,
      transform = "log10",
      limits = GSEA_PADJ_LIMITS,
      breaks = GSEA_PADJ_BREAKS,
      labels = GSEA_PADJ_LABELS,
      oob    = scales::squish,
      na.value = "grey90",
      guide = guide_colorbar(barwidth = 1, barheight = 10, ticks.linewidth = 1)
    ) +
    scale_size_continuous(range = c(1, 6), transform = "identity") +
    guides(size = guide_legend(title = "Geneset Size")) +
    geom_text(aes(label = marker,
                  x = ifelse(NES > 0, NES + 0.25, NES - 0.4)),
              hjust = 0, vjust = 0.5, size = 6, color = "black") +
    geom_text(aes(label = marker2,
                  x = ifelse(NES > 0, NES + 0.25, NES - 0.4)),
              hjust = 0, vjust = 0.5, size = 4, color = "black") +
    geom_text(aes(label = marker3,
                  x = ifelse(NES > 0, NES + 0.25, NES - 0.4)),
              hjust = 0, vjust = 0.5, size = 6, color = "red") +
    scale_x_continuous(limits = c(x_min, x_max),
                       breaks = seq(x_min, x_max, step)) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.05)))

  width  <- max(12, 5 + (x_max - x_min) * 0.5)
  height <- max(8, min(25, 4 + n_pathways * 0.7))
  if (max(nchar(filtered_data$pathway)) > 30) {
    p <- p + theme(axis.text.y = element_text(margin = margin(r = 15)))
    width <- width + 2
  }

  tp_tag <- if (!is.null(timepoint)) paste0("_", timepoint) else ""
  ggsave(paste0(directory, celltype, tp_tag,
                ifelse(simple_tf, "_simple", "_multilevel"),
                "_fGSEA_NES_", pathway_name, "_dotplot.png"),
         plot = p, dpi = 320, width = width, height = height)

  p
}
