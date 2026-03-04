setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(SeuratObject)
library(scCustomize)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(glmGamPoi)
library(metap)
library(RCurl)

sessionInfo()

## read in data ----
# initialize empty list to store seurat objects
mammary.list = list()

# Loop over numbers in file directory names and read in data into seurat objects
for (i in 21:32) {
  # Construct the file name
  file_name = paste0('10x_summary_8651-JC/8651-JC-', i, '_sample_filtered_feature_bc_matrix/')
  
  # Read the data into a Seurat object
  data = Read10X(file_name) 
  
  # Create variable name for data
  var_name = paste0('count', i)
  
  # Create a Seurat object and assign to var_name
  assign(var_name, CreateSeuratObject(counts = data, 
                                      assay='RNA'), envir=.GlobalEnv)
  
  # Add the Seurat object to the list
  mammary.list[i-20] = get(var_name, envir=.GlobalEnv)
}
# assign names to count data in list
names(mammary.list) = sprintf('count%s', seq(21,32)) 

rm(data, var_name, file_name, i) # clean up RAM
gc()
# count21@assays[['RNA']]$counts # to view counts matrix
# count21@assays[['RNA']]$data # to view normalized matrix
# count21@assays[['RNA']]$c # to view scaled

# which mt genes are present?
# rownames(count21)[which(startsWith(rownames(count21),'mt'))]
#  [1] "mt-Nd1"  "mt-Nd2"  "mt-Co1"  "mt-Co2"  "mt-Atp8" "mt-Atp6" "mt-Co3"  "mt-Nd3" 
# [9] "mt-Nd4l" "mt-Nd4"  "mt-Nd5"  "mt-Nd6"  "mt-Cytb"


# "Complexity
# We can evaluate each cell in terms of how complex the RNA species are by 
# using a measure called the novelty score. The novelty score is computed by 
# taking the ratio of nGenes over nUMI. If there are many captured 
# transcripts (high nUMI) and a low number of genes detected in a cell, 
# this likely means that you only captured a low number of genes and simply 
# sequenced transcripts from those lower number of genes over and over again. 
# These low complexity (low novelty) cells could represent a specific 
# cell type (i.e. red blood cells which lack a typical transcriptome), 
# or could be due to an artifact or contamination. Generally, we expect 
# the novelty score to be above 0.80 for good quality cells."
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

## calculate percent mitochondrial genes and novelty and log10 z-scores score ----
mammary.list = lapply(mammary.list, function(x) {
  log10genes = log10(x$nFeature_RNA)
  log10counts = log10(x$nCount_RNA)
  x$log10GenesPerUMI = log10genes / log10counts # novelty score to measure complexity
  
  # mean and standard deviation of log10 ngenes and ncount
  mean_ngene = mean(log10genes)
  mean_ncount = mean(log10counts)
  std_ngene = sd(log10genes)
  std_ncount = sd(log10counts)
  
  # z-scores of log10 ngenes and ncounts
  x$log10_genes_zscored = (log10genes - mean_ngene) / std_ngene 
  x$log10_counts_zscored = (log10counts - mean_ncount) / std_ncount
  x = PercentageFeatureSet(x, pattern = "^mt-", col.name = "percent.mt") # percent mitochondrial genes
})


# calculate log transformed median, median absolute deviation (MAD), and +- 3x MAD
mads_df = data.frame(median_gene=numeric(12), MAD_gene=numeric(12), plus_3MAD_gene=numeric(12), minus_3MAD_gene=numeric(12), 
                     median_count=numeric(12), MAD_count=numeric(12), plus_3MAD_count=numeric(12), minus_3MAD_count=numeric(12))
for (i in 1:12) {
  s = unique(metadata$sample)[i]
  MAD_gene = mad(log(metadata[metadata$sample == s,]$nFeature_RNA))
  MAD_count = mad(log(metadata[metadata$sample == s,]$nCount_RNA))
  median_gene = median(log(metadata[metadata$sample == s,]$nFeature_RNA))
  median_count = median(log(metadata[metadata$sample == s,]$nCount_RNA))
  plus_3MAD_gene = median_gene + 3*MAD_gene
  minus_3MAD_gene = median_gene - 3*MAD_gene
  plus_3MAD_count = median_count + 3*MAD_count
  minus_3MAD_count = median_count - 3*MAD_count
  
  mads_df[i,] = list(median_gene, MAD_gene, plus_3MAD_gene, minus_3MAD_gene, median_count, MAD_count, plus_3MAD_count, minus_3MAD_count)
}

  
## plot proportion of mitochondrial genes per cell ----
plots = lapply(names(mammary.list), function(name) {
  VlnPlot(mammary.list[[name]], features = "percent.mt", assay='RNA', pt.size = 0.05)+
    NoLegend()+
    labs(title=name)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank())
})

plot_grid = wrap_plots(plots, ncol=4)+
  plot_annotation(title='Proportion of mitochondrial genes per cell',
                  theme = theme(plot.title = element_text(hjust = 0.5)))+
  plot_layout(guides='collect')

plot_grid = wrap_elements(plot_grid) +
  labs(tag = "percent.mt (%)") +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )
plot_grid
ggsave('percent_mt.png', plot=plot_grid, dpi=320, width=10, height= 10)


## plot number of reads per cell ----
plots = lapply(names(mammary.list), function(name) {
  VlnPlot(mammary.list[[name]], features = "nCount_RNA", assay='RNA', pt.size = 0.05)+
    NoLegend()+
    labs(title=name)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank())+
    scale_y_log10(n.breaks=8, labels = comma_format())+
    annotation_logticks(sides = "l")
})

plot_grid = wrap_plots(plots, ncol=4)+
  plot_annotation(title='Number of transcript reads per cell',
                  theme = theme(plot.title = element_text(hjust = 0.5)))+
  plot_layout(guides='collect')

plot_grid = wrap_elements(plot_grid) +
  labs(tag = "# of reads (log scaled axis)") +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )
plot_grid
ggsave('figures/read_counts_vln.png', plot=plot_grid, dpi=320, width=10, height= 10)


## plot number of genes per cell ----
plots = lapply(names(mammary.list), function(name) {
  VlnPlot(mammary.list[[name]], features = "nFeature_RNA", assay='RNA', pt.size = 0.05)+
    NoLegend()+
    labs(title=name)+
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.x=element_blank())+
    scale_y_log10(n.breaks=8, labels = comma_format())+
    annotation_logticks(sides = "l")
})

plot_grid = wrap_plots(plots, ncol=4)+
  plot_annotation(title='Number of genes per cell',
                  theme = theme(plot.title = element_text(hjust = 0.5)))+
  plot_layout(guides='collect')

plot_grid = wrap_elements(plot_grid) +
  labs(tag = "# of genes (log scaled axis)") +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )
plot_grid
ggsave('figures/gene_counts_vln.png', plot=plot_grid, dpi=320, width=10, height= 10)

## plot count depth (number of reads) vs number of genes, colored percent mitochondrial gene counts ---- 
plots = lapply(names(mammary.list), function(name) {
  data = mammary.list[[name]]
  p = ggplot(data@meta.data,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
    geom_point(size=0.5)+
    scale_color_viridis(option='viridis')+
    DarkTheme()+
    labs(title=paste0(name,' - count depth vs # of genes'),
         subtitle=paste0('colored by percent mitochondrial (mt) counts. # of cells = ',
                         length(data$nCount_RNA)),
         x='# of UMI counts (log scale)',
         y='# of genes (log scale)',
         color='% mt counts')+
    scale_x_log10(n.breaks=8, labels=comma_format())+
    scale_y_log10(n.breaks=8, labels=comma_format())+
    annotation_logticks(sides = "bl", colour='white', size=1.5)+
    theme(plot.title=element_text(colour='white'), 
          axis.title = element_text(size=12), 
          axis.text = element_text(size=12),
          title = element_text(size=18),
          legend.text=element_text(size=18))+
    NoGrid()+
    stat_smooth(method=lm,size=0.5, colour='#61e4e8')
  ggsave(paste0('figures/',name, '_ncounts_vs_ngenes_coloredby_percentmt.png'), 
                plot=p, dpi=320, width=10, height= 10)
})

## plot z-scores for count depth (number of reads) vs number of genes, colored percent mitochondrial gene counts ---- 
plots = lapply(names(mammary.list), function(name) {
  data = mammary.list[[name]]
  p = ggplot(data@meta.data,aes(x=log10_counts_zscored,y=log10_genes_zscored,color=percent.mt))+
    geom_point(size=0.5)+
    scale_color_viridis(option='viridis')+
    DarkTheme()+
    labs(title=paste0(name,' - count depth vs # of genes (z-scores) '),
         subtitle=paste0('colored by percent mitochondrial (mt) counts. # of cells = ',
                         length(data$nCount_RNA)),
         x='# of UMI counts (z-scores)',
         y='# of genes (z-scores)',
         color='% mt counts')+
    scale_x_continuous(n.breaks=10)+
    scale_y_continuous(n.breaks=10)+
    theme(plot.title=element_text(colour='white'), 
          axis.title = element_text(size=12), 
          axis.text = element_text(size=12),
          title = element_text(size=18),
          legend.text=element_text(size=18),
          axis.ticks=)+
    NoGrid()+
    stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')
  ggsave(paste0('figures/',name, '_z-scores_ncounts_vs_ngenes_coloredby_percentmt.png'), 
         plot=p, dpi=320, width=10, height= 10)
})

## plot log10genes / log10UMI distirbutions ----
plots = lapply(names(mammary.list), function(name) {
  data=mammary.list[[name]]@meta.data
  ggplot(data, aes(x=log10GenesPerUMI))+
    geom_density()+
    theme_classic2()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())+
    labs(title=name)
})

plot_grid = wrap_plots(plots, ncol=4)+
  plot_annotation(title='Genes detected per UMI',
                  theme = theme(plot.title = element_text(hjust = 0.5)))+
  plot_layout(guides='collect')

plot_grid = wrap_elements(plot_grid) +
  labs(tag = "density") +
  theme(
    plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left"
  )
plot_grid
ggsave('complexity_novelty_scores.png', plot=plot_grid, dpi=320, width=10, height= 10)


# show plots
# for (p in plots) {
#   print(p)
# }


## plot histogram of reads per cell ----
# plots = lapply(names(mammary.list), function(name) {
#   ggplot(mammary.list[[name]]@meta.data, aes(x=nCount_RNA))+
#     geom_histogram(bins=50)+
#     labs(title=name)+
#     theme(axis.ticks.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.title.x=element_blank())+
#     xlim(0,2000)+
#     theme_classic2()
# })
# 
# plot_grid = wrap_plots(plots, ncol=4)+
#   plot_annotation(title='Number of transcript reads per cell',
#                   theme = theme(plot.title = element_text(hjust = 0.5)))+
#   plot_layout(guides='collect')
# 
# plot_grid


lapply(names(mammary.list), function(name) {
  x = mammary.list[[name]]
  pre = length(x$nCount_RNA)
  print(paste0('Dataset: ', name))
  print(paste0('# of cells = ', pre))
  x = subset(x, subset=(log10_genes_zscored <= 3) & (log10_counts_zscored <= 3) &
               (percent.mt <= 15))
  post = length(x$nCount_RNA)
  print(paste0(pre - post, ' cells filtered out'))
  print(paste0('# of cells after filtering = ', post))
  cat('\n')
})

# cell cycle genes, formatted to match genes in matrix with first letter uppercased
s.genes = stringr::str_to_title(cc.genes$s.genes)
g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)



## ensebmle and annotation hub to get cell cycle and ribosomal genes



## filter and normalize data and calculate cell cycle scores ----
ncells.list = vector(length=12)
ncells_filtered.list = vector(length=12)
nvariable_genes.list = vector(length=12)
i=1
mammary.list = lapply(seq_along(mammary.list), function(i) {
  x = mammary.list[[i]]
  ncells = length(x$nCount_RNA)
  ncells.list[[i]] = ncells
  x = subset(x, )
  x = NormalizeData(x)
  x = FindVariableFeatures(x)
  x = ScaleData(x)
  # x = SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
  x = CellCycleScoring(x, s.features=s.genes, g2m.features=g2m.genes)
})
length(count21$nCount_RNA)

# rownames(count32)[which(s.genes %in% rownames(count24@assays[["RNA"]]$counts))]\

# save.image('sctransformed_data.RData')
# saveRDS(mammary.list, 'mammary_counts_list.rds')
mammary.list = readRDS('data/mammary_counts_list.rds')

# update individual variables to reflect changes to variables in list
for (name in names(mammary.list)) {
  assign(name, mammary.list[[name]], envir = .GlobalEnv)
}

# integrate control cells
ctrl.list = c(count21, count22, count23, count24, count25, count26)



# merge lead exposed cells


## integrate data
mammary.anchors = FindIntegrationAnchors(object.list = mammary.list, dims = 1:30)


# calculate log transformed median, median absolute deviation (MAD), and +- 3x MAD ----
mads_df = data.frame(median_gene=numeric(12), MAD_gene=numeric(12), plus_3MAD_gene=numeric(12), minus_3MAD_gene=numeric(12), 
                     median_count=numeric(12), MAD_count=numeric(12), plus_3MAD_count=numeric(12), minus_3MAD_count=numeric(12))
for (i in 1:12) {
  s = unique(metadata$sample)[i]
  MAD_gene = mad(log(metadata[metadata$sample == s,]$nFeature_RNA))
  MAD_count = mad(log(metadata[metadata$sample == s,]$nCount_RNA))
  median_gene = median(log(metadata[metadata$sample == s,]$nFeature_RNA))
  median_count = median(log(metadata[metadata$sample == s,]$nCount_RNA))
  plus_3MAD_gene = median_gene + 3*MAD_gene
  minus_3MAD_gene = median_gene - 3*MAD_gene
  plus_3MAD_count = median_count + 3*MAD_count
  minus_3MAD_count = median_count - 3*MAD_count
  
  mads_df[i,] = list(median_gene, MAD_gene, plus_3MAD_gene, minus_3MAD_gene, median_count, MAD_count, plus_3MAD_count, minus_3MAD_count)
}
saveRDS(mads_df, 'data/median_absolute_deviations_of_log_transformed_genesandcounts.rds')




