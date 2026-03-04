library(tidyverse)
library(fgsea)
# library(msigdbr)
library(msigdb)
library(biomaRt)
# library(org.Mm.eg.db)
# library(org.Hs.eg.db)
library(ExperimentHub)
# library(enrichplot)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(ggthemes)
library(ggdark)
library(igraph)
library(ggraph)
# library(ggrepel)
library(qs2)
library(jsonlite)


# set number of threads to use for saving qs formatted data
nthreads <- 2 

# set seed
set.seed(2024)

# directories for saving
# selected_timepoint <- 'week3'
# selected_timepoint <- 'week10'
# selected_timepoint <- 'month7'
# selected_timepoint <- 'month18'

de_directory <- "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/results/"

# directory to save figures and results in ----

# analysis type, cellbender or cellranger
analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analys

data_dir =  "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/data/"

directory_figs = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/fGSEA/figures/"
directory_results = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/fGSEA/results/"

dirs = c(data_dir, directory_figs, directory_results)

for (d in dirs) {
  if (!dir.exists(d)) { dir.create(d,
                                    recursive=TRUE)}
}

### grab celltypes from DE results ----


### grab de result file names and extract cell types 
de_results = list.files(path=de_directory, pattern='all_genes')

### List of cell types (replace with your actual cell types)
celltypes = sub(pattern='_condition.*',replacement='',x=de_results)
celltypes



# load in data ----
res_tbl_adip = read.csv(paste0(de_directory,'Adipocytes_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_endo = read.csv(paste0(de_directory,'Endothelial_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_macro = read.csv(paste0(de_directory,'Macrophage.Ma_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_lum = read.csv(paste0(de_directory,'Luminal.HS_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_bcell = read.csv(paste0(de_directory,'B cells_condition_pb_vs_ctrl_all_genes_deseq2.csv'))

# # get gene annotations, ensembl and entrez IDs and human gene homologs ----

# ensembl = useMart('ensembl')
# listDatasets(ensembl)
# ensembl=useDataset("mmusculus_gene_ensembl",mart=ensembl)
# 
# listFilters(ensebml)
# grep('sapien',listAttributes(ensembl,page='homologs'), value=TRUE)
# 
# test = listAttributes(ensembl)


# mice_id = getBM(attributes= c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
#                  mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))) %>%
#   mutate(external_gene_name = na_if(external_gene_name, ""))
# 
# 
# human_id = getBM(attributes= c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name",
#                                "hsapiens_homolog_ensembl_gene"),
#                  mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))) %>%
#   mutate(hsapiens_homolog_associated_gene_name = na_if(hsapiens_homolog_associated_gene_name,"")) %>%
#   mutate(hsapiens_homolog_ensembl_gene = na_if(hsapiens_homolog_ensembl_gene,""))
# 
# 
# gene_id = merge(mice_id, human_id, by='ensembl_gene_id') %>%
#   distinct()
# 
# write.csv(gene_id, 'annotation/ensembl_and_entrez_gene_ids_and_human_homologs.csv')

# load fetched gene IDs and human homologs ----
gene_id = read.csv('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/annotation/ensembl_and_entrez_gene_ids_and_human_homologs.csv')

# test = left_join(res_tbl_adip, msigdbr::msigdbr(species='Mus musculus',category='H'),
                                               # by=join_by('gene'=='gene_symbol'))



# define functions ----
get_ids_and_homologs= function(res_tbl, gene_ids=gene_id) {
  res_tbl = merge(x=res_tbl, y=gene_ids, by.x ="gene", by.y="external_gene_name")
}

to_genelist = function(res_tbl, species_human=F) {
  if(!species_human) {
    genelist = res_tbl %>%
      as_tibble() %>%
      # group_by(hsapiens_homolog_associated_gene_name)
      # summarize(pvalue=mean(pvalue)) %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue))  %>%
      dplyr::select(gene, ranking_metric) %>%
      na.omit() %>%
      # group_by(hsapiens_homolog_associated_gene_name) %>%
      arrange(desc(ranking_metric)) %>%
      distinct() %>%
      deframe()
  } else {
    genelist = res_tbl %>%
      as_tibble() %>%
      group_by(hsapiens_homolog_associated_gene_name) %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
      summarize(ranking_metric = max(sign(log2FoldChange) * -log10(pvalue), na.rm = TRUE)) %>%
      dplyr::select(hsapiens_homolog_associated_gene_name, ranking_metric) %>%
      na.omit() %>%
      arrange(desc(ranking_metric)) %>%
      distinct() %>%
      deframe()
  }
  
  return(genelist)
}

## tidy up and sort fGSEA results by NES and save to csv
tidysave_fgsea = function(fgseaRes, pathways, pathway_name, celltype, simple_tf=F, 
                          directory=directory_results) {
  
  if (!dir.exists(directory)) { dir.create(directory,
                                           recursive=TRUE)}
  
  fgseaResTidy = fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    mutate(leadingEdge2 = leadingEdge) %>% # for saving to csv
    rowwise() %>% 
    mutate_at(c('leadingEdge2'), ~paste(unlist(.), collapse = ',')) %>%
    dplyr::select(!leadingEdge) %>%
    rename(leadingEdge = leadingEdge2)
  
  # save to csv, omit leadingEdge col, located at index 8,
  # from saved file since it contains lists. 
  # write.csv(x=fgseaResTidy,
  #           file = paste0(directory,celltype,ifelse(simple_tf,'_simple','_multilevel'),
  #                         '_fGSEA_NES_', deparse(substitute(pathways)),'.csv'),
  #           row.names=FALSE)
  
  write.csv(x=fgseaResTidy,
            file = paste0(directory,celltype,ifelse(simple_tf,'_simple','_multilevel'),
                          '_fGSEA_NES_', pathway_name,'.csv'),
            row.names=FALSE)
  
  return(fgseaResTidy)
}


# Function to wrap pathway names
wrap_pathway_names <- function(pathways, width = 50) {
  wrapped <- character(length(pathways))
  
  for (i in seq_along(pathways)) {
    path <- pathways[i]
    
    # Only wrap if longer than 60 characters
    if (nchar(path) >= 60) {
      # Replace underscores with spaces first
      path_readable <- gsub("_", " ", path)
      
      # Then wrap the text
      wrapped[i] <- stringr::str_wrap(path_readable, width)
    } else {
      # Keep original if not long
      wrapped[i] <- path
    }
  }
  
  return(wrapped)
}

## function to plot and save normalized enrichment scores (NES) horizontal barplot ----
plot_fgsea_nes = function(fgsea_nes, celltype, pathways, pathway_name,
                          padj_cutoff=0.25, simple_tf=F, directory=directory_figs) {
  
  
  
  # stop if no data
  if(nrow(fgsea_nes) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    # Create an empty plot or return NULL
    return(ggplot() + 
             geom_text(aes(x = 0, y = 0, label = paste("No significant pathways at padj <=", padj_cutoff))) +
             theme_minimal())
  }
  
  
  if (!dir.exists(directory)) { dir.create(directory,
                                           recursive=TRUE)}
  
  
  # Filter data for significant results
  filtered_data <- fgsea_nes %>% 
  filter(padj <= padj_cutoff) %>%
  mutate(
    pathway_wrapped = wrap_pathway_names(pathway, width = 50),
    marker = ifelse(padj <= 0.05, '*', ''),
    marker2 = ifelse((padj <= 0.1) & (padj > 0.05), '#', ''),
    marker3 = ifelse(padj <= 0.005, '**', '')
  ) %>%
    filter(!is.na(NES), !is.na(pathway_wrapped))
  
  if(nrow(filtered_data) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    # Create an empty plot or return NULL
    return(ggplot() + 
             geom_text(aes(x = 0, y = 0, label = paste("No significant pathways at padj <=", padj_cutoff))) +
             theme_minimal())
  }
  
  cat("Number of rows in fgsea_nes:", nrow(fgsea_nes), "\n")
  cat("Number of rows after filtering for padj <=", padj_cutoff, ":", nrow(filtered_data), "\n")
  cat("NA values in NES:", sum(is.na(filtered_data$NES)), "\n")
  cat("NA values in pathway_wrapped:", sum(is.na(filtered_data$pathway_wrapped)), "\n")
  
  # Count number of pathways to adjust plot dimensions later
  n_pathways <- nrow(filtered_data)
  
  # First check if filtered_data has rows
  if (nrow(filtered_data) > 0) {
    # Remove non-finite values
    nes_values <- filtered_data$NES[is.finite(filtered_data$NES)]
    
    # Check if there are any valid NES values left
    if (length(nes_values) > 0) {
      x_min <- floor(min(nes_values) - 0.5)
      x_max <- ceiling(max(nes_values) + 0.5)
    } else {
      # Default values if no valid NES values
      x_min <- -3
      x_max <- 3
    }
  } else {
    # Default values if no filtered data
    x_min <- -3
    x_max <- 3
  }
  
  # Make sure x_min and x_max are not equal (would cause seq to fail)
  if (x_min == x_max) {
    x_min <- x_min - 1
    x_max <- x_max + 1
  }
  
  # Make sure the range isn't too large (to avoid creating too many breaks)
  if (x_max - x_min > 20) {
    # Use a larger step size
    step <- ceiling((x_max - x_min) / 10)
  } else {
    step <- 1
  }
  
  p = ggplot(data = filtered_data,
             aes(x=NES, y=reorder(pathway_wrapped, NES))) +
    geom_col(aes(fill=padj)) +
    labs(x="Normalized Enrichment Score", y="Pathway",
         title=paste0('GSEA - ', pathway_name,
                      ' Pathways: ',celltype),
         subtitle="** denotes padj < 0.005, * denotes padj < 0.05, # denotes padj < 0.1",
         fill=paste0('FDR up to ', padj_cutoff)) + 
    theme_light()+
    theme(axis.line = element_line(colour = "black"),
          panel.border=element_blank(),
          
          # Other theme elements
          panel.grid.major.x = element_line(color = "grey95"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey95"),
          panel.grid.minor.y = element_blank(),
          legend.position = "right",
          legend.box = "vertical",
          
          # Increase text size for all elements
          text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12, margin = margin(r = 5)),
          axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 16, margin = margin(b = 5)),
          plot.subtitle = element_text(size = 12, margin = margin(b = 5)),
          
          # Minimize margins around the plot
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          
          # Reduce space between legend items
          legend.key.height = unit(0.8, "lines"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
          
    )+
    
    scale_fill_viridis(option='viridis',
                        direction=-1,
                        transform='log',
                        n.breaks=10,
                        # breaks=c(0,0.001,0.005,0.01,0.025,0.05,0.1,0.25),
                        # labels=c(0,0.001,0.005,0.01,0.025,0.05,0.1,0.25),
                       # breaks=c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.07,0.1,0.17,0.25),
                       # labels=c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.07,0.1,0.17,0.25),
                       guide = guide_colorbar(
                         barwidth = 1,      # Make the bar wider
                         barheight = 10,    # Make the bar taller
                         ticks.linewidth = 1
                       ))+
    
    geom_text(aes(label = marker, x = ifelse(NES >0, NES + 0.25, NES - 0.4)), # Adjust x to position the marker next to the bars
              hjust = 0, # Adjust text alignment
              vjust = 0.5, # Vertical alignment
              size = 6, # Adjust size of the asterisk
              color = "black") +
    geom_text(aes(label = marker2, x = ifelse(NES >0, NES + 0.25, NES - 0.4)), # Adjust x to position the marker next to the bars
              hjust = 0, # Adjust text alignment
              vjust = 0.5, # Vertical alignment
              size = 4, # Adjust size of the pount sign
              color = "black") +
    geom_text(aes(label = marker3, x = ifelse(NES >0, NES + 0.25, NES - 0.4)), # Adjust x to position the marker next to the bars
              hjust = 0, # Adjust text alignment
              vjust = 0.5, # Vertical alignment
              size = 6, # Adjust size of the double asterisk
              color = "red")
  
  # Calculate dimensions - define before using
  width <- max(12, 5 + (x_max - x_min) * 0.5)  # Increased from 10 to 12
  height <- max(8, min(25, 4 + n_pathways * 0.7))  # Increased height per pathway
  
  if(max(nchar(filtered_data$pathway)) > 30) {
    p <- p + theme(axis.text.y = element_text(margin = margin(r = 15)))
    width <- width + 2  # More width for larger text
  }
  
  # ggsave(paste0(directory,celltype,ifelse(simple_tf,'_simple','_multilevel'),
  #               '_fGSEA_NES_',ifelse(is.symbol(substitute(pathways)),
  #                                              deparse(substitute(pathways)),
  #                                              pathways),'_hbarplot.png'),
  #        plot=p,dpi=320,width=10,height=20)
  
  ggsave(paste0(directory,celltype,ifelse(simple_tf,'_simple','_multilevel'),
                '_fGSEA_NES_', pathway_name,'_hbarplot.png'),
         plot=p,dpi=320,width=width,height=height)
  
  return(p)
}  

## function to plot and save normalized enrichment scores (NES) as dot plot----
plot_fgsea_nes_dotplot = function(fgsea_nes, celltype, pathways,pathway_name,
                                  padj_cutoff=0.25, simple_tf=F, directory=directory_figs) {

  
  
  # stop if no data
  if(nrow(fgsea_nes) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    # Create an empty plot or return NULL
    return(ggplot() + 
             geom_text(aes(x = 0, y = 0, label = paste("No significant pathways at padj <=", padj_cutoff))) +
             theme_minimal())
  }
  
  if (!dir.exists(directory)) { dir.create(directory,
                                           recursive=TRUE)}
  
  # Filter data for significant results
  filtered_data <- fgsea_nes %>% 
    filter(padj <= padj_cutoff) %>%
    mutate(
      pathway_wrapped = wrap_pathway_names(pathway, width = 50),
      marker = ifelse(padj <= 0.05, '*', ''),
      marker2 = ifelse((padj <= 0.1) & (padj > 0.05), '#', ''),
      marker3 = ifelse(padj <= 0.005, '**', '')
    ) %>%
    filter(!is.na(NES), !is.na(pathway_wrapped))
  
  if(nrow(filtered_data) == 0) {
    warning("No pathways meet the significance threshold of padj <= ", padj_cutoff)
    # Create an empty plot or return NULL
    return(ggplot() + 
             geom_text(aes(x = 0, y = 0, label = paste("No significant pathways at padj <=", padj_cutoff))) +
             theme_minimal())
  }
    
    
    cat("Number of rows in fgsea_nes:", nrow(fgsea_nes), "\n")
    cat("Number of rows after filtering for padj <=", padj_cutoff, ":", nrow(filtered_data), "\n")
    cat("NA values in NES:", sum(is.na(filtered_data$NES)), "\n")
    cat("NA values in pathway_wrapped:", sum(is.na(filtered_data$pathway_wrapped)), "\n")
  
  # Count number of pathways to adjust plot dimensions later
  n_pathways <- nrow(filtered_data)
  
  # First check if filtered_data has rows
  if (nrow(filtered_data) > 0) {
    # Remove non-finite values
    nes_values <- filtered_data$NES[is.finite(filtered_data$NES)]
    
    # Check if there are any valid NES values left
    if (length(nes_values) > 0) {
      x_min <- floor(min(nes_values) - 0.5)
      x_max <- ceiling(max(nes_values) + 0.5)
    } else {
      # Default values if no valid NES values
      x_min <- -3
      x_max <- 3
    }
  } else {
    # Default values if no filtered data
    x_min <- -3
    x_max <- 3
  }
  
  # Make sure x_min and x_max are not equal (would cause seq to fail)
  if (x_min == x_max) {
    x_min <- x_min - 1
    x_max <- x_max + 1
  }
  
  # Make sure the range isn't too large (to avoid creating too many breaks)
  if (x_max - x_min > 20) {
    # Use a larger step size
    step <- ceiling((x_max - x_min) / 10)
  } else {
    step <- 1
  }
  
  p = ggplot(data = filtered_data,
             aes(x=NES, y=reorder(pathway_wrapped, NES),
                 size=size, color=padj)) +
    geom_point() +
    labs(x="Normalized Enrichment Score", y="Pathway",
         title=paste0('GSEA - ', pathway_name,
                      ' Pathways: ',celltype),
         subtitle="** denotes padj < 0.005, * denotes padj < 0.05, # denotes padj < 0.1",
         size='Geneset Size',
         color=paste0('FDR up to ', padj_cutoff)) + 
    theme_light()+
    theme(axis.line = element_line(colour = "black"),
          panel.border=element_blank(),
          
          # Other theme elements
          panel.grid.major.x = element_line(color = "grey95"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey95"),
          panel.grid.minor.y = element_blank(),
          legend.position = "right",
          legend.box = "vertical",
          
          # Increase text size for all elements
          text = element_text(size = 12),
          axis.title = element_text(size = 14),
          axis.text.y = element_text(size = 12, margin = margin(r = 5)),
          axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 16, margin = margin(b = 5)),
          plot.subtitle = element_text(size = 12, margin = margin(b = 5)),
          
          # Minimize margins around the plot
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
          
          # Reduce space between legend items
          legend.key.height = unit(0.8, "lines"),
          legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
          
          )+
    
    scale_color_viridis(option='viridis',
                       direction=-1,
                       transform='log',
                       n.breaks=10,
                       # breaks=c(0,0.001,0.005,0.01,0.025,0.05,0.1,0.25),
                       # labels=c(0,0.001,0.005,0.01,0.025,0.05,0.1,0.25),
                       # breaks=c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.07,0.1,0.17,0.25),
                       # labels=c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.07,0.1,0.17,0.25),
                       guide = guide_colorbar(
                         barwidth = 1,      # Make the bar wider
                         barheight = 10,    # Make the bar taller
                         ticks.linewidth = 1
                       ))+
                       # breaks=c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.07,0.1,0.17,0.25),
                       # labels=c(0,0.001,0.002,0.005,0.01,0.02,0.05,0.07,0.1,0.17,0.25))+
    
    # guides(size=guide_legend(override.aes=list(color='black',
    #                                            fill='black')))+
    scale_size_continuous(range = c(1,6),
                          # breaks = c(10,,20,30,40,50, 100, 250, 500, 1000),
                          # labels= c(10,,20,30,40,50, 100, 250, 500, 1000),
                          transform = "identity") + 
    
    guides(size = guide_legend(title = "Geneset Size")) +
    
    geom_text(aes(label = marker, x = ifelse(NES >0, NES + 0.25, NES - 0.4)), # Adjust x to position the marker next to the bars
               hjust = 0, # Adjust text alignment
               vjust = 0.5, # Vertical alignment
               size = 6, # Adjust size of the asterisk
               color = "black") +
    geom_text(aes(label = marker2, x = ifelse(NES >0, NES + 0.25, NES - 0.4)), # Adjust x to position the marker next to the bars
              hjust = 0, # Adjust text alignment
              vjust = 0.5, # Vertical alignment
              size = 4, # Adjust size of the asterisk
              color = "black") +
    geom_text(aes(label = marker3, x = ifelse(NES >0, NES + 0.25, NES - 0.4)), # Adjust x to position the marker next to the bars
              hjust = 0, # Adjust text alignment
              vjust = 0.5, # Vertical alignment
              size = 6, # Adjust size of the asterisk
              color = "red") +
    
    scale_x_continuous(limits = c(x_min, x_max),
                       breaks = seq(x_min, x_max, step))+
    # Add this line to reduce spacing between y-axis elements
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.05)))
  
  
  # Calculate dimensions - define before using
  width <- max(12, 5 + (x_max - x_min) * 0.5)  # Increased from 10 to 12
  height <- max(8, min(25, 4 + n_pathways * 0.7))  # Increased height per pathway
  
  if(max(nchar(filtered_data$pathway)) > 30) {
    p <- p + theme(axis.text.y = element_text(margin = margin(r = 15)))
    width <- width + 2  # More width for larger text
  }
  
  # ggsave(paste0(directory,celltype,ifelse(simple_tf,'_simple','_multilevel'),
  #               '_fGSEA_NES_',ifelse(is.symbol(substitute(pathways)),
  #                                              deparse(substitute(pathways)),
  #                                              pathways),'_dotplot.png'),
  #        plot=p,dpi=320,width=10,height=20)
  
  ggsave(paste0(directory,celltype,ifelse(simple_tf,'_simple','_multilevel'),
                '_fGSEA_NES_', pathway_name,'_dotplot.png'),
         plot=p,dpi=320,width=width,height=height)

  
  return(p)
}



# deparse(substitute(variable_name)) # to getvariable name as character string


### merge cell-specific deseq2 data with gene IDs and human homologs ----
res_tbl_adip = get_ids_and_homologs(res_tbl_adip)
res_tbl_endo = get_ids_and_homologs(res_tbl_endo)
res_tbl_macro = get_ids_and_homologs(res_tbl_macro)
res_tbl_lum = get_ids_and_homologs(res_tbl_lum)
res_tbl_bcell = get_ids_and_homologs(res_tbl_bcell)


# adipo_genelist = res_tbl_adip %>%
#   as_tibble() %>%
#   # group_by(hsapiens_homolog_associated_gene_name)
#   # summarize(pvalue=mean(pvalue)) %>%
#   mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
#   dplyr::select(hsapiens_homolog_associated_gene_name, ranking_metric) %>%
#   na.omit() %>%
#   group_by(hsapiens_homolog_associated_gene_name) %>%
#   summarize(ranking_metric = mean(ranking_metric)) %>%
#   arrange(desc(ranking_metric)) $>$
#   distinct()

# adipo_genelist = res_tbl_adip %>%
#   as_tibble() %>%
#   group_by(hsapiens_homolog_associated_gene_name) %>%
#   summarize(pvalue=mean(pvalue)) %>%
#   mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
#   dplyr::select(gene, ranking_metric) %>%
#   na.omit() %>%
#   # group_by(hsapiens_homolog_associated_gene_name) %>%
#   arrange(desc(ranking_metric)) %>%
#   distinct() %>%
#   deframe()


### Import mice gene sets for analysis ----
# msigdb.mm = getMsigdb(org = 'mm', id = 'SYM') # snapshotDate(): 2024-04-29 # snapshotDate(): 2024-10-24
# msigdb.mm = appendKEGG(msigdb.mm,version = '2023.1')
# qs_save(msigdb.mm, 'c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/msigdb.mm.qs',nthreads=2)
msigdb.mm = qs_read('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/msigdb.mm.qs',nthreads=2)

listCollections(msigdb.mm)
# [1] "c1" "c3" "c2" "c8" "c6" "c7" "c4" "c5" "h" 
listSubCollections(msigdb.mm)
# [1] "MIR:MIR_LEGACY"  "TFT:TFT_LEGACY"  "CGP"             "TFT:GTRD"       
# [5] "VAX"             "CP:BIOCARTA"     "CGN"             "GO:BP"          
# [9] "GO:CC"           "IMMUNESIGDB"     "GO:MF"           "HPO"            
# [13] "MIR:MIRDB"       "CM"              "CP"              "CP:PID"         
# [17] "CP:REACTOME"     "CP:WIKIPATHWAYS" "CP:KEGG"

### Import human gene sets for analaysis ----
# msigdb.hs = getMsigdb(org='hs', id ='SYM') # snapshotDate(): 2024-04-29 # snapshotDate(): 2024-10-24
# msigdb.hs = appendKEGG(msigdb.hs, version='2023.1')
# qs_save(msigdb.hs, 'c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/msigdb.hs.qs',nthreads=2)
msigdb.hs = qs_read('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/msigdb.hs.qs',nthreads=2)


listCollections(msigdb.hs)
# [1] "c1" "c2" "c3" "c4" "c6" "c7" "c8" "h"  "c5"
listSubCollections(msigdb.hs)
# [1] "CGP"             "CP:BIOCARTA"     "CP:PID"          "CP"             
# [5] "MIR:MIRDB"       "MIR:MIR_LEGACY"  "TFT:TFT_LEGACY"  "CGN"            
# [9] "CM"              "IMMUNESIGDB"     "VAX"             "TFT:GTRD"       
# [13] "HPO"             "CP:REACTOME"     "CP:WIKIPATHWAYS" "GO:BP"          
# [17] "GO:CC"           "GO:MF"           "CP:KEGG"   


### list GSEA pathways of interest
pathway_subs = c('h','GO:BP','GO:MF','GO:CC','CP:REACTOME','CP:KEGG','CP:WIKIPATHWAYS','ALL')
pathway_names = c('Hallmarks','GO: Biological Processes','GO: Molecular Function',
                  'GO: Cellular Compartment','Reactome','KEGG','Wikipathways', 'ALL')

all(length(pathway_names) == length(pathway_subs))
# [1] TRUE if take away ALL from pathway_names

# list species database to be used
msigdb_species = c('Mouse','Human')


### gene set collections for mice ----
hallmarks_mm = geneIds(subsetCollection(msigdb.mm, collection='h'))
gobp_mm = geneIds(subsetCollection(msigdb.mm, subcollection='GO:BP'))
gomf_mm = geneIds(subsetCollection(msigdb.mm, subcollection='GO:MF'))
gocc_mm = geneIds(subsetCollection(msigdb.mm, subcollection='GO:CC'))
reactome_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:REACTOME'))
kegg_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:KEGG'))
cgp_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CGP'))
ccca_mm = geneIds(subsetCollection(msigdb.mm, collection='3CA'))
c6_mm = geneIds(subsetCollection(msigdb.mm, collection = 'c6'))
hpo_mm = geneIds(subsetCollection(msigdb.mm, subcollection='HPO'))
wiki_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:WIKIPATHWAYS'))
# all_mm = c(hallmarks_mm, gobp_mm, gomf_mm, reactome_mm, kegg_mm, wiki_mm)


mm_collections = list(
  'hallmarks_mm' = hallmarks_mm,
  'gobp_mm' = gobp_mm,
  'gomf_mm' = gomf_mm,
  'gocc_mm' = gocc_mm,
  'reactome_mm' = reactome_mm,
  'kegg_mm' = kegg_mm,
  'wiki_mm' = wiki_mm
  # 'all_mm' = all_mm
)


# de_results = c('Adipocytes','Endothelial','Macrophage.Ma','Luminal.HS','B cells')

# fGSEA on all celltype differential expression results ----
# for (i in 1:length(de_results)) {
# for (i in c(2,3,8,10,11)) {
for (i in c(10,11)) {  
  res_tbl = read.csv(paste0(de_directory,de_results[i]))
  res_tbl = get_ids_and_homologs(res_tbl)
  
  genelist = to_genelist(res_tbl)
    
  # for (j in 1:2) { # to go through mouse and human versions of msigdb
  for (j in 1) {
    for (k in 1:length(pathway_names)) {
      gc()
      cat('\n\n','Setting up pathways for ',celltypes[i],' using ',pathway_names[k],' pathways...\n\n')
      # select gsea pathway
      
      pathways = mm_collections[[k]]
      # pathways = ifelse(j == 1,
      #                   ifelse(pathway_subs[k] == 'h',
      #                   geneIds(subsetCollection(msigdb.mm, collection=pathway_subs[k])),
      #                   geneIds(subsetCollection(msigdb.mm, subcollection=pathway_subs[k]))
      #                   ),
      #                   ifelse(pathway_subs[k] == 'h',
      #                          geneIds(subsetCollection(msigdb.hs, collection=pathway_subs[k])),
      #                          geneIds(subsetCollection(msigdb.hs, subcollection=pathway_subs[k]))
      #                   )
      # )
      
      # pathway database name. Species followed by collection name
      pathway_name = ifelse(j==1, paste0(msigdb_species[1],'_',mm_collections[k]),
                            paste0(msigdb_species[2],'_',names(mm_collections[k])))
      
      cat('\n\n','Running fGSEA for ',celltypes[i],' using ',pathway_name,' pathways...\n\n')
      
      fgseaRes = fgsea(pathways=pathways,
                       stats=genelist,
                       sampleSize=101, # default 101. change to nperm for simple-fGSEA
                       minSize=10,
                       maxSize=1000,
                       eps=1e-50, # default 1e-50
                       nproc=6)
      gc()
      
      ### tidy fGSEA results, sort by NES, and save to csv
      tidysave_fgsea(fgseaRes,
                     pathways=pathway_name,
                     pathway_name = pathway_name,
                     celltype=celltypes[i])
      
    
      
      # select independent pathways to get rid of redundancy
      cat('collapsing pathways for ',celltypes[i],' ',pathway_name,
          'if not part of 50 hallmarks...\n\n')
      
      collapsedPathways = collapsePathways(fgseaRes[order(pval)][padj < 0.1], 
                                           pathways[names(pathways) %in% fgseaRes$pathway], 
                                           genelist,
                                           pval.threshold=0.05,
                                           nperm=1000)
      
      qs_save(collapsedPathways, paste0(directory_results,celltypes[i],'_',pathway_name,'_collapsedpathways.qs'),nthreads=nthreads)
      
      if (length(collapsedPathways$mainPathways) > 0) {
        # save collapedPathways as csv by transforming to dataframe
        main_df <- data.frame(
          pathway = collapsedPathways$mainPathways,
          parent = NA  # No parent for main pathways
        )
        
        # For parentPathways - create a data frame with pathway names, their parents, and type
        parent_df <- data.frame(
          pathway = names(collapsedPathways$parentPathways),
          parent = collapsedPathways$parentPathways
        )
        
        # Combine both data frames
        combined_df <- rbind(main_df, parent_df)
        
        write.csv(combined_df,paste0(directory_results,celltypes[i],'_',pathway_name,'_collapsedpathways.csv'))
        
        # create graph of collapsedPathways
        # Create edge list from your data
        edges <- data.frame(
          from = collapsedPathways$parentPathways,
          to = names(collapsedPathways$parentPathways)
        )
        
        # Remove NA parents (these are root nodes)
        edges <- edges[!is.na(edges$from), ]
        
        # Create the graph
        g <- graph_from_data_frame(edges, directed = TRUE)
        
        png(paste0(directory_figs,'_',celltypes[i],'_',pathway_name,'_collapsedPathways_hierarchy.png'),
            unit='in',height = 15,width = 15,res = 180)
        # Plot with non-overlapping labels
        plot(g, 
             layout = layout_with_fr(g),  # Force-directed layout spreads nodes
             vertex.size = 6,             # Smaller nodes
             vertex.label.dist = 0.7,     # Move labels away from nodes
             vertex.label.cex = 0.7,      # Smaller label text
             vertex.label.degree = -pi/4, # Position labels at an angle
             edge.arrow.size = 0.5,
             margin = c(0.2, 0.2, 0.2, 0.2))  # Add margin around the plot
        dev.off()
        
        
        if (pathway_subs[k] == 'h') {
          mainPathways = names(pathways)
        } else {
          mainPathways = fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
        }
        
        
        saveRDS(collapsedPathways,paste0(directory_results,celltypes[i],'_',pathway_name,'_collapsedPathways.rds'))
      } else {
        mainPathways = NULL
      }
      
      padj_cutoff = ifelse(sum(fgseaRes$padj < 0.25) > 70, 0.1, 0.25)
      
      
      ### plot NES results
      cat('plotting NES barplot for',celltypes[i],'using',pathway_name,'pathways p.adj cuttoff of:',padj_cutoff,'...\n\n')
      if((length(mainPathways) != 0) & (sum(fgseaRes$padj < 0.25) > 50)) {
        plot_fgsea_nes(fgsea_nes=fgseaRes[pathway %in% mainPathways,],
                       pathways=pathways,
                       pathway_name = pathway_name,
                       celltype=celltypes[i],
                       padj_cutoff=padj_cutoff,
                       directory = directory_figs)
      } else {
        plot_fgsea_nes(fgsea_nes=fgseaRes,
                       pathways=pathways,
                       pathway_name = pathway_name,
                       celltype=celltypes[i],
                       padj_cutoff=padj_cutoff,
                       directory = directory_figs)
      }
      
      cat('plotting NES dotplot for ',celltypes[i],' using ',pathway_name,' pathways p.adj cuttoff of: ',padj_cutoff,'...\n\n')
      if((length(mainPathways) != 0) & (sum(fgseaRes$padj < 0.25) > 50)) {
        plot_fgsea_nes_dotplot(fgsea_nes=fgseaRes[pathway %in% mainPathways,],
                               pathways=pathways,
                               pathway_name = pathway_name,
                               celltype=celltypes[i],
                               padj_cutoff=padj_cutoff)
      } else {
        plot_fgsea_nes_dotplot(fgsea_nes=fgseaRes,
                               pathways=pathways,
                               pathway_name = pathway_name,
                               celltype=celltypes[i],
                               padj_cutoff=padj_cutoff)
      }
      
      # cutoff = 0.25
      # i = names(hallmarks_mm) %in% adipo.fgsea.hallmarks.mm$pathway[adipo.fgsea.hallmarks.mm$padj < cutoff]
      
      # update_geom_defaults("segment", list(color = "black"))
      
      # cat('plotting enrichment graphs for ',celltypes[i],' using ',pathway_name,' pathways...\n\n')
      # p = plotGseaTable(pathways=ifelse(length(pathways) > 50,
      #                                   ifelse(length(mainPathways) != 0,
      #                                          pathways[mainPathways],
      #                                          pathways),
      #                                   pathways),
      #                   stats=genelist, 
      #                   fgseaRes=fgseaRes,
      #                   pathwayLabelStyle=list(color='black'),
      #                   headerLabelStyle =list(size=20, color='black'),
      #                   valueStyle=list(color='black'),
      #                   axisLabelStyle=list(color='black'),
      #                   gseaParam=0.5)
      # cowplot::ggdraw(p)+
      #   theme(plot.background = element_rect(fill='white', color = NA),
      #         panel.background = element_rect(fill='white', color = NA)
      #   )
      # ggsave(paste0(directory_figs,celltypes[i],'_multilevel',
      #               '_fGSEA_NES_',pathway_name,'.png'),
      #        plot=p,dpi=320,width=10,height=15,units='in')
      rm(fgseaRes,padj_cutoff,collapsePathway, mainPathways)
      
      cat('\nfinished ',pathway_name,' for ',celltypes[i],'\n')
    }
  }
}







## Adipocytes fGSEA ----
## hallmarks ----
### retrieve hallmarks collections ----
hallmarks_mm = geneIds(subsetCollection(msigdb.mm, collection='h'))
hallmarks_hs = geneIds(subsetCollection(msigdb.hs, collection='h'))

hallmarks_mm %>%
  head() %>%
  lapply(head)

hallmarks_hs %>%
  head() %>%
  lapply(head)

### convert adip0cyte table to ranked gene list ----
adipo_genelist_mm = to_genelist(res_tbl_adip)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nperm=10000 if simple is used.
adipo.fgsea.hallmarks.mm = fgsea(pathways=hallmarks_mm,
                                stats=adipo_genelist_mm,
                                sampleSize=101, # default 101. change to nproc for simple-fGSEA
                                minSize=10,
                                maxSize=1000,
                                eps=1e-50, # default 1e-50
                                nproc=4) # default 0

# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.hallmarks.mm[order(pval)][padj < 0.1], 
                                     hallmarks_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.hallmarks.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.hallmarks.mm,
               pathways=hallmarks_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.hallmarks.mm,
               pathways=hallmarks_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.hallmarks.mm,
               pathways=hallmarks_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)


### plot enrichment scores with fGSEA for adipocytes ----

# get indices where pathways are below an adjusted p-value threshold
update_geom_defaults("segment", list(color = "black"))

p = plotGseaTable(pathways=hallmarks_mm, 
              stats=adipo_genelist_mm, 
              fgseaRes=adipo.fgsea.hallmarks.mm,
              pathwayLabelStyle=list(color='black'),
              headerLabelStyle =list(size=20, color='black'),
              valueStyle=list(color='black'),
              axisLabelStyle=list(color='black'),
              gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs, 'adipocytes_enrichment_graphs_hallmarks_mm.png'), units='in', width=10, height=15, dpi=320)



# # Show in a nice table:
# adipo.fgsea.hallmarks.mm %>% 
#   dplyr::select(-leadingEdge, -ES) %>% 
#   arrange(padj) %>% 
#   DT::datatable()


# view differential expression of genes present in hallmark pathways
View(hallmarks_mm %>% 
  enframe("pathway", "gene") %>% 
  unnest(cols=c(gene)) %>% 
  inner_join(res_tbl_adip, by='gene') %>%
  filter(padj < 0.1))


## GO Biological Processes ----
### retrieve GO:BP collections ----
gobp_mm = geneIds(subsetCollection(msigdb.mm, subcollection='GO:BP'))
gobp_hs = geneIds(subsetCollection(msigdb.hs, subcollection='GO:BP'))

gobp_mm %>%
  head() %>%
  lapply(head)

gobp_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.gobp.mm = fgsea(pathways=gobp_mm,
                                 stats=adipo_genelist_mm,
                                 sampleSize=101, # default 101. change to nproc for simple-fGSEA
                                 minSize=10,
                                 maxSize=1000,
                                 eps=1e-50, # default 1e-50
                                 nproc=4) # default 0

# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.gobp.mm[order(pval)][padj < 0.1], 
                                      gobp_mm, adipo_genelist_mm,
                                      pval.threshold=0.05,
                                      nperm=1000)

mainPathways = adipo.fgsea.gobp.mm[pathway %in% collapsedPathways$mainPathway][
  order(-NES), pathway]


### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.gobp.mm,
               pathways=gobp_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.gobp.mm[pathway %in% mainPathways,],
               pathways=gobp_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.gobp.mm[pathway %in% mainPathways,],
               pathways=gobp_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)



### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))

# # get top upregulated and downregulated pathways
# topPathwaysUp <- adipo.fgsea.gobp.mm[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- adipo.fgsea.gobp.mm[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# plotGseaTable(msigdb_kegg_ids[topPathwaysUp], Day2.rnk, kegg.gsea.day2,
#               gseaParam=0.5)

p = plotGseaTable(pathways=gobp_mm[mainPathways], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.gobp.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_gobp_mm.png'), units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in gobp pathways
View(gobp_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))


## GO Molecular Functions ----
### retrieve GO:MF collections ----
gomf_mm = geneIds(subsetCollection(msigdb.mm, subcollection='GO:MF'))
gomf_hs = geneIds(subsetCollection(msigdb.hs, subcollection='GO:MF'))

gomf_mm %>%
  head() %>%
  lapply(head)

gomf_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.gomf.mm = fgsea(pathways=gomf_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=10,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.gomf.mm[order(pval)][padj < 0.1], 
                                     gomf_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.gomf.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.gomf.mm,
               pathways=gomf_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.gomf.mm[pathway %in% mainPathways,],
               pathways=gomf_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.gomf.mm[pathway %in% mainPathways,],
               pathways=gomf_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)


### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))

# # get top upregulated and downregulated pathways
# topPathwaysUp <- adipo.fgsea.gomf.mm[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- adipo.fgsea.gomf.mm[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# plotGseaTable(msigdb_kegg_ids[topPathwaysUp], Day2.rnk, kegg.gsea.day2,
#               gseaParam=0.5)

p = plotGseaTable(pathways=gomf_mm[mainPathways], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.gomf.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_gomf_mm.png'), units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in gomf pathways
View(gomf_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))




## GO Cellular Compartment ----
### retrieve GO:CC collections ----
gocc_mm = geneIds(subsetCollection(msigdb.mm, subcollection='GO:CC'))
gocc_hs = geneIds(subsetCollection(msigdb.hs, subcollection='GO:CC'))

gocc_mm %>%
  head() %>%
  lapply(head)

gocc_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.gocc.mm = fgsea(pathways=gocc_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=10,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.gocc.mm[order(pval)][padj < 0.1], 
                                     gocc_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.gocc.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.gocc.mm,
               pathways=gocc_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.gocc.mm[pathway %in% mainPathways,],
               pathways=gocc_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.gocc.mm[pathway %in% mainPathways,],
                       pathways=gocc_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.25)


### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))

# # get top upregulated and downregulated pathways
# topPathwaysUp <- adipo.fgsea.gocc.mm[ES > 0][head(order(pval), n=10), pathway]
# topPathwaysDown <- adipo.fgsea.gocc.mm[ES < 0][head(order(pval), n=10), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
# plotGseaTable(msigdb_kegg_ids[topPathwaysUp], Day2.rnk, kegg.gsea.day2,
#               gseaParam=0.5)

p = plotGseaTable(pathways=gocc_mm[mainPathways], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.gocc.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_gocc_mm.png'), units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in gocc pathways
View(gocc_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))






## try runnning fGSEA for human homologs for adipocytes ----
adipo_genelist_hs = to_genelist(res_tbl_adip,species_human=TRUE)

### run fGSEA on adipocytes ----

# check out multilevel fGSEA

# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.hallmarks.hs = fgsea(pathways=hallmarks_hs,
                                 stats=adipo_genelist_hs,
                                 sampleSize=101, # default 101. change to nproc for simple-fGSEA
                                 minSize=10,
                                 maxSize=1000,
                                 eps=1e-50, # default 1e-50
                                 nproc=4) # default 0


### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.hallmarks.hs,
               pathways=hallmarks_hs,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.hallmarks.hs,
               pathways=hallmarks_hs,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.hallmarks.hs,
               pathways=hallmarks_hs,
               celltype='Adipocytes',
               padj_cutoff=0.25)




### plot enrichment scores with fGSEA for adipocytes ----

# get indices where pathways are below an adjusted p-value threshold
cutoff = 0.25
i = names(hallmarks_hs) %in% adipo.fgsea.hallmarks.hs$pathway[adipo.fgsea.hallmarks.hs$padj < cutoff]

update_geom_defaults("segment", list(color = "black"))

p = plotGseaTable(pathways=hallmarks_hs[i], 
                  stats=adipo_genelist_hs, 
                  fgseaRes=adipo.fgsea.hallmarks.hs,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_hallmarks_hs.png'), units='in', width=10, height=15, dpi=320)



# Show in a nice table:
adipo.fgsea.hallmarks.mm %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()


# view differential expression of genes present in hallmark pathways
View(hallmarks_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))


## Reactome ----
### retrieve Reactome collections ----
reactome_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:REACTOME'))
reactome_hs = geneIds(subsetCollection(msigdb.hs, subcollection='CP:REACTOME'))

reactome_mm %>%
  head() %>%
  lapply(head)

reactome_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.reactome.mm = fgsea(pathways=reactome_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=10,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.reactome.mm[order(pval)][padj < 0.1],
                                     reactome_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.reactome.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.reactome.mm,
               pathways=reactome_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.reactome.mm,
               pathways=reactome_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.reactome.mm,
               pathways=reactome_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)



if (!dir.exists("stats/GSEA/fGSEA/figures")) { dir.create("stats/GSEA/fGSEA/figures",
                                                          recursive=TRUE)}
### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))

p = plotGseaTable(pathways=reactome_mm, 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.reactome.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_reactome_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in reactome pathways
View(reactome_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))




## KEGG ----
### retrieve Reactome collections ----
kegg_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:KEGG'))
kegg_hs = geneIds(subsetCollection(msigdb.hs, subcollection='CP:KEGG'))

kegg_mm %>%
  head() %>%
  lapply(head)

kegg_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.kegg.mm = fgsea(pathways=kegg_mm,
                                stats=adipo_genelist_mm,
                                sampleSize=101, # default 101. change to nproc for simple-fGSEA
                                minSize=10,
                                maxSize=1000,
                                eps=1e-50, # default 1e-50
                                nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.kegg.mm[order(pval)][padj < 0.1],
                                     kegg_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.kegg.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.kegg.mm,
               pathways=kegg_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.kegg.mm,
               pathways=kegg_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.kegg.mm,
                       pathways=kegg_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.25)



if (!dir.exists("stats/GSEA/fGSEA/figures")) { dir.create("stats/GSEA/fGSEA/figures",
                                                          recursive=TRUE)}
### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=kegg_mm, 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.kegg.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_kegg_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in kegg pathways
View(kegg_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))






## ALL Msigdb ----
### retrieve Reactome collections ----
all_mm = c(hallmarks_mm, gobp_mm, gobp_mm, reactome_mm,)
all_hs = geneIds(msigdb.hs)

all_mm %>%
  head() %>%
  lapply(head)

all_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.all.mm = fgsea(pathways=all_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=15,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.all.mm[order(pval)][padj < 0.25],
                                     all_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.all.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

saveRDS(collapsedPathways,'stats/GSEA/fGSEA/results/Adipocytes_allpathways_collapsedPathways.rds')
### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.all.mm,
               pathways=all_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.all.mm[pathway %in% mainPathways,],
               pathways=all_mm,
               celltype='Adipocytes',
               padj_cutoff=0.1)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.all.mm[pathway %in% mainPathways,],
                       pathways=all_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.05)



if (!dir.exists("stats/GSEA/fGSEA/figures")) { dir.create("stats/GSEA/fGSEA/figures",
                                                          recursive=TRUE)}
### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=all_mm[pathway %in% mainPathways,], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.all.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_all_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in all pathways
View(all_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))




##  CGP, chemical and genetic perturbations ----
### retrieve Reactome collections ----
cgp_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CGP'))
cgp_hs = geneIds(subsetCollection(msigdb.hs, subcollection='CGP'))

cgp_mm %>%
  head() %>%
  lapply(head)

cgp_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.cgp.mm = fgsea(pathways=cgp_mm,
                           stats=adipo_genelist_mm,
                           sampleSize=101, # default 101. change to nproc for simple-fGSEA
                           minSize=10,
                           maxSize=1000,
                           eps=1e-50, # default 1e-50
                           nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.cgp.mm[order(pval)][padj < 0.1],
                                     cgp_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.cgp.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.cgp.mm,
               pathways=cgp_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.cgp.mm[pathway %in% mainPathways,],
               pathways=cgp_mm,
               celltype='Adipocytes',
               padj_cutoff=0.05)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.cgp.mm[pathway %in% mainPathways,],
                       pathways=cgp_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.05)




### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=cgp_mm[pathway %in% mainPathways,], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.cgp.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_cgp_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in cgp pathways
View(cgp_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))




## c4, computational sets, includes 3CA? curated cancer cell atlas Msigdb ----
ccca_mm = geneIds(subsetCollection(msigdb.mm, collection='c4'))
ccca_hs = geneIds(subsetCollection(msigdb.hs, collection='c4'))

ccca_mm %>%
  head() %>%
  lapply(head)

ccca_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.ccca.mm = fgsea(pathways=ccca_mm,
                           stats=adipo_genelist_mm,
                           sampleSize=101, # default 101. change to nproc for simple-fGSEA
                           minSize=10,
                           maxSize=1000,
                           eps=1e-50, # default 1e-50
                           nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.ccca.mm[order(pval)][padj < 0.1],
                                     ccca_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.ccca.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

saveRDS(collapsedPathways,paste0(directory_results, 'Adipocytes_cccapathways_collapsedPathways.rds'))
### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.ccca.mm,
               pathways=ccca_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.ccca.mm[pathway %in% mainPathways,],
               pathways=ccca_mm,
               celltype='Adipocytes',
               padj_cutoff=0.1)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.ccca.mm[pathway %in% mainPathways,],
                       pathways=ccca_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.1)



### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=ccca_mm[pathway %in% mainPathways,], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.ccca.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_ccca_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in ccca pathways
View(ccca_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))





## C6, oncogene gene signatures ----
c6_mm = geneIds(subsetCollection(msigdb.mm, collection = 'c6'))
c6_hs = geneIds(subsetCollection(msigdb.hs, collection = 'c6'))

c6_mm %>%
  head() %>%
  lapply(head)

c6_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.c6.mm = fgsea(pathways=c6_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=10,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.c6.mm[order(pval)][padj < 0.1],
                                     c6_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.c6.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

saveRDS(collapsedPathways,paste0(directory_results,'Adipocytes_c6pathways_collapsedPathways.rds'))
### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.c6.mm,
               pathways=c6_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.c6.mm[pathway %in% mainPathways,],
               pathways=c6_mm,
               celltype='Adipocytes',
               padj_cutoff=0.1)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.c6.mm[pathway %in% mainPathways,],
                       pathways=c6_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.05)



### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=c6_mm[pathway %in% mainPathways,], 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.c6.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_c6_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in c6 pathways
View(c6_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))






## HPO ----
hpo_mm = geneIds(subsetCollection(msigdb.mm, subcollection='HPO'))
hpo_hs = geneIds(subsetCollection(msigdb.hs, subcollection='HPO'))

hpo_mm %>%
  head() %>%
  lapply(head)

hpo_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.hpo.mm = fgsea(pathways=hpo_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=10,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.hpo.mm[order(pval)][padj < 0.1],
                                     hpo_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.hpo.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.hpo.mm,
               pathways=hpo_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.hpo.mm,
               pathways=hpo_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.hpo.mm,
                       pathways=hpo_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.25)



if (!dir.exists("stats/GSEA/fGSEA/figures")) { dir.create("stats/GSEA/fGSEA/figures",
                                                          recursive=TRUE)}
### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=hpo_mm, 
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.hpo.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave('stats/GSEA/fGSEA/figures/adipocytes_enrichment_graphs_hpo_mm.png', units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in hpo pathways
View(hpo_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))


## wikipathways ----
wiki_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:WIKIPATHWAYS'))
wiki_hs = geneIds(subsetCollection(msigdb.hs, subcollection='CP:WIKIPATHWAYS'))

wiki_mm %>%
  head() %>%
  lapply(head)

wiki_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.wiki.mm = fgsea(pathways=wiki_mm,
                           stats=adipo_genelist_mm,
                           sampleSize=101, # default 101. change to nproc for simple-fGSEA
                           minSize=10,
                           maxSize=1000,
                           eps=1e-50, # default 1e-50
                           nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.wiki.mm[order(pval)][padj < 0.1],
                                     wiki_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.wiki.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.wiki.mm,
               pathways=wiki_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.wiki.mm,
               pathways=wiki_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.wiki.mm,
                       pathways=wiki_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.25)



if (!dir.exists("stats/GSEA/fGSEA/figures")) { dir.create("stats/GSEA/fGSEA/figures",
                                                          recursive=TRUE)}
### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=wiki_mm[mainPathways],
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.wiki.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_wiki_mm.png'), units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in wiki pathways
View(wiki_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))



## c8, celltypes ----
c8_celltype_mm = geneIds(subsetCollection(msigdb.mm, collection='c8'))
c8_celltype_hs = geneIds(subsetCollection(msigdb.hs, collection='c8'))

c8_celltype_mm %>%
  head() %>%
  lapply(head)

c8_celltype_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.c8_celltype.mm = fgsea(pathways=c8_celltype_mm,
                            stats=adipo_genelist_mm,
                            sampleSize=101, # default 101. change to nproc for simple-fGSEA
                            minSize=10,
                            maxSize=1000,
                            eps=1e-50, # default 1e-50
                            nproc=4) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.c8_celltype.mm[order(pval)][padj < 0.1],
                                     c8_celltype_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.c8_celltype.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.c8_celltype.mm,
               pathways=c8_celltype_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.c8_celltype.mm[pathway %in% mainPathways,],
               pathways=c8_celltype_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.c8_celltype.mm[pathway %in% mainPathways,],
                       pathways=c8_celltype_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.25)




### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=c8_celltype_mm[mainPathways],
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.c8_celltype.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_c8_celltype_mm.png'), units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in c8_celltype pathways
View(c8_celltype_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))


# ALL, almost ----
all_mm = c(hallmarks_mm, gobp_mm, gomf_mm, reactome_mm, kegg_mm)


all_mm %>%
  head() %>%
  lapply(head)

all_hs %>%
  head() %>%
  lapply(head)

### run fGSEA on adipocytes ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nproc=10000 if simple is used.
adipo.fgsea.all.mm = fgsea(pathways=all_mm,
                                   stats=adipo_genelist_mm,
                                   sampleSize=101, # default 101. change to nproc for simple-fGSEA
                                   minSize=10,
                                   maxSize=1000,
                                   eps=1e-50, # default 1e-50
                                   nproc=6) # default 0


# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(adipo.fgsea.all.mm[order(pval)][padj < 0.1],
                                     all_mm, adipo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = adipo.fgsea.all.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy adipocytes fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(adipo.fgsea.all.mm,
               pathways=all_mm,
               celltype='Adipocytes')

### plot NES for adipocytes ----
plot_fgsea_nes(fgsea_nes=adipo.fgsea.all.mm[pathway %in% mainPathways,],
               pathways=all_mm,
               celltype='Adipocytes',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=adipo.fgsea.all.mm[pathway %in% mainPathways,],
                       pathways=all_mm,
                       celltype='Adipocytes',
                       padj_cutoff=0.25)




### plot enrichment scores with fGSEA for adipocytes ----

update_geom_defaults("segment", list(color = "black"))


p = plotGseaTable(pathways=all_mm[mainPathways],
                  stats=adipo_genelist_mm, 
                  fgseaRes=adipo.fgsea.all.mm,
                  pathwayLabelStyle=list(color='black'),
                  headerLabelStyle =list(size=20, color='black'),
                  valueStyle=list(color='black'),
                  axisLabelStyle=list(color='black'),
                  gseaParam=0.5)
cowplot::ggdraw(p)+
  theme(plot.background = element_rect(fill='white', color = NA),
        panel.background = element_rect(fill='white', color = NA)
  )
ggsave(paste0(directory_figs,'adipocytes_enrichment_graphs_all_mm.png'), units='in', width=10, height=15, dpi=320)



# view differential expression of genes present in all pathways
View(all_mm %>% 
       enframe("pathway", "gene") %>% 
       unnest(cols=c(gene)) %>% 
       inner_join(res_tbl_adip, by='gene') %>%
       filter(padj < 0.25))



# Endothelial cells fGSEA ----
## hallmarks ----
### convert adip0cyte table to ranked gene list ----
endo_genelist_mm = to_genelist(res_tbl_endo)

### run fGSEA on Endothelial ----
# Can set fGSEA method to multi-level or simple depending on
# whether sampleSize or nproc paramter is used. Set nperm=10000 if simple is used.
endo.fgsea.hallmarks.mm = fgsea(pathways=hallmarks_mm,
                                 stats=endo_genelist_mm,
                                 sampleSize=101, # default 101. change to nproc for simple-fGSEA
                                 minSize=10,
                                 maxSize=1000,
                                 eps=1e-50, # default 1e-50
                                 nproc=4) # default 0

# select independent pathways to get rid of redundancy
collapsedPathways = collapsePathways(endo.fgsea.hallmarks.mm[order(pval)][padj < 0.1], 
                                     hallmarks_mm, endo_genelist_mm,
                                     pval.threshold=0.05,
                                     nperm=1000)

mainPathways = endo.fgsea.hallmarks.mm[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

### tidy Endothelial fGSEA results, sort by NES, and save to csv ----
tidysave_fgsea(endo.fgsea.hallmarks.mm,
               pathways=hallmarks_mm,
               celltype='Endothelial')

### plot NES for Endothelial ----
plot_fgsea_nes(fgsea_nes=endo.fgsea.hallmarks.mm,
               pathways=hallmarks_mm,
               celltype='Endothelial',
               padj_cutoff=0.25)

plot_fgsea_nes_dotplot(fgsea_nes=endo.fgsea.hallmarks.mm,
                       pathways=hallmarks_mm,
                       celltype='Endothelial',
                       padj_cutoff=0.25)

