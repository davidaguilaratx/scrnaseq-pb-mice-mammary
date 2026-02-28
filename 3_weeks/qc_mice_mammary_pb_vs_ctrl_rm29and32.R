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
library(Polychrome)



# set seed
set.seed(2024)

# create polychrome color palette ----

p12 = createPalette(12, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
swatch(p12)
names(p12) = NULL

## read in data ----
# initialize empty list to store seurat objects
mammary.list = list()

# Loop over numbers in file directory names and read in data into seurat objects.
for (i in 21:32) {
  # Construct the file name
  file_name = paste0('10x_summary_8651-JC/8651-JC-', i, '_sample_filtered_feature_bc_matrix/')
  
  # Read the data into a Seurat object
  data = Read10X(file_name) 
  
  # Create variable name for data
  var_name = paste0('count', i)
  
  # Create a Seurat object and assign to var_name\
  assign(var_name, CreateSeuratObject(counts = data, 
                                        assay='RNA', 
                                        project=paste0('count',i)), envir=.GlobalEnv)
 
  # Add the Seurat object to the list
  mammary.list[i-20] = get(var_name, envir=.GlobalEnv)
}

# assign names to count data in list
names(mammary.list) = sprintf('count%s', seq(21,32)) 

rm(data, var_name, file_name, i) # clean up RAM
gc()

## calculate log10 z-scores of genes and counts ----
mammary.list = lapply(mammary.list, function(x) {
  log10genes = log10(x$nFeature_RNA)
  log10counts = log10(x$nCount_RNA)

  # mean and standard deviation of log10 ngenes and ncount
  mean_ngene = mean(log10genes)
  mean_ncount = mean(log10counts)
  std_ngene = sd(log10genes)
  std_ncount = sd(log10counts)
  
  # z-scores of log10 ngenes and ncounts
  x$log10_genes_zscored = (log10genes - mean_ngene) / std_ngene 
  x$log10_counts_zscored = (log10counts - mean_ncount) / std_ncount
  
  # Might as well do mitochondrial expression percentage here
  x = PercentageFeatureSet(x, pattern = "^mt-", col.name = "percent.mt") # percent mitochondrial genes
})

# create list of sample names, separating count data into ctrl and pb.
cell_ids = unlist(lapply(names(mammary.list), function(name) {
  n = as.numeric(gsub('[^0-9]', '', name))
  n = ifelse(n < 27, paste0('ctrl',n), paste0('pb',n))
}))

# merge data into one object ----
merged_seurat = merge(mammary.list[[1]], mammary.list[-1], add.cell.ids=cell_ids)

# Add number of log10 genes per log10 UMI for each cell to metadata
# acts as novelty score to measure complexity
merged_seurat$log10GenesPerUMI = log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# # Compute percent mito percentage
# merged_seurat = PercentageFeatureSet(object = merged_seurat, 
#                                                 pattern = "^mt-", 
#                                                 col.name = "percent.mt")
rm(mammary.list) # clean up RAM
gc()

# creat metadata dataframe
metadata = merged_seurat@meta.data

# add a cell IDs column to metadata
metadata$cells = rownames(metadata)

# add sample type column to metadata
metadata$sample = str_extract(metadata$cells, '(\\D+)\\d{2}')

# add condition column to metadata
metadata$condition = str_extract(metadata$cells, '(\\D+)')

metadata$exposed = ifelse(metadata$condition == 'pb', 1, 0)

# save metadata dataframe back to merged seurat data
merged_seurat@meta.data = metadata


saveRDS(merged_seurat, 'data/merged_seurat.rds')
saveRDS(metadata, 'data/metadata.rds')

# number of cells per dataset
table(Idents(merged_seurat))
# count21 count22 count23 count24 count25 count26 count27 count28 count29 count30 count31 
# 2142    3035    4567    1720    6151    5671    7976    4396    1249    6560    2096 
# count32 
# 529 

## plot number of cells per dataset ----
ggplot(metadata, aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic2() +
  NoLegend()+
  labs(title="Number of cells per experiment",
       subtitle = paste0('Total cells = ',dim(metadata)[1],', # ctrl cells = ',
                         sum(metadata$condition =='ctrl'),', # pb cells = ',
                         sum(metadata$condition == 'pb')))+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  scale_y_continuous(n.breaks=8)
ggsave('figures/qc/ncells.png', dpi=320, width=10, height= 10)


## proposed QC cutoffs ----
percent.mt.cutoff = 10
gene_cutoff_low = 150
gene_cutoff_high = 5000
umi_cutoff_low = 150
umi_cutoff_high = 50000 
complexity_cutoff = 0.8
min_cells_per_gene = 20 # filtered genes must be expressed in at least this # of cells



## plot proportion of mitochondrial genes per cell ----
VlnPlot(merged_seurat, features = "percent.mt", group.by='sample', pt.size = 0.05)+
  NoLegend()+
  labs(title='Percent mitochondrial genes',
       subtitle=paste0('cutoff above ', percent.mt.cutoff,'%'),
       x='sample',
       y='percent.mt')+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  geom_hline(yintercept=percent.mt.cutoff, linetype='dashed',color='red')+
  annotate('text', x=7.5, y=10.5, label='cutoff', size=5, color='red')
ggsave('figures/qc/percent_mt.png', dpi=320, width=10, height= 10)


## plot novelty score (complexity = log10genes / log10UMI) ----
ggplot(metadata, aes(color=sample, x=log10GenesPerUMI, fill=sample)) + 
  geom_density(alpha = 0.1)+
  theme_classic2()+
  labs(title='Overal complexity via genes detected per UMI (novelty score)',
       x='log10(nGenes) / log10(nUMI)',
       y='Cell density')+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  geom_vline(xintercept = complexity_cutoff, linetype='dashed', color='red')+
  annotate('text', x=complexity_cutoff-0.02, y=40, label='cutoff', size=5, color='red')
ggsave('figures/qc/complexity_novelty_scores.png', dpi=320, width=10, height= 10)



## plot number of reads per cell as vln plots ----
VlnPlot(merged_seurat, features = "nCount_RNA", group.by='sample', pt.size = 0.05)+
  theme_classic2()+
  NoLegend()+
  labs(title='Number of transcript reads per cell',
       subtitle=paste0('cutoffs below ',umi_cutoff_low,' and above ',umi_cutoff_high, ' UMIs'),
       x='sample',
       y='nUMI')+
  scale_y_log10(n.breaks=10, limits=c(100,NA),labels = comma_format())+
  annotation_logticks(sides = "l")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  geom_hline(yintercept=umi_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=umi_cutoff_high, linetype='dashed',color='red')
ggsave('figures/qc/read_counts.png', dpi=320, width=10, height= 10)


# plot number of reads per cell as density plots ----
ggplot(metadata, aes(color=sample, x=nCount_RNA, fill=sample))+ 
  geom_density(alpha = 0.1)+
  labs(title='Number of transcript reads per cell',
       subtitle=paste0('cutoffs below ',umi_cutoff_low,' and above ',umi_cutoff_high, ' UMIs'),
       x='nUMI',
       y='Cell density')+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  scale_x_log10(n.breaks=10, labels=comma_format())+
  annotation_logticks(sides = "b")+
  geom_vline(xintercept = umi_cutoff_low, linetype='dashed', color='red')+
  geom_vline(xintercept = umi_cutoff_high, linetype='dashed', color='red')
ggsave('figures/qc/read_counts_density.png', dpi=320, width=10, height= 10)


## plot number of genes per cell as vln plots ----
VlnPlot(merged_seurat, features = "nFeature_RNA", group.by='sample', pt.size = 0.05)+
  theme_classic2()+
  NoLegend()+
  labs(title='Number of genes detected per cell',
       subtitle=paste0('cutoffs below ',gene_cutoff_low,' and above ',gene_cutoff_high, ' genes'),
       x='sample',
       y='nGenes')+
  scale_y_log10(n.breaks=10, limits=c(100,NA),labels = comma_format())+
  annotation_logticks(sides = "l")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')
ggsave('figures/qc/gene_counts.png', dpi=320, width=10, height= 10)

## plot number of genes per cell as density plots ----
ggplot(metadata, aes(color=sample, x=nFeature_RNA, fill=sample))+
  geom_density(alpha = 0.1)+
  labs(title='Number of genes detected per cell',
       subtitle=paste0('cutoffs below ',gene_cutoff_low,' and above ',gene_cutoff_high, ' genes'),
       x='nGene',
       y='Cell density')+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))+
  scale_x_log10(n.breaks=10, labels=comma_format())+
  annotation_logticks(sides = "b")+
  geom_vline(xintercept = gene_cutoff_low, linetype='dashed', color='red')+
  geom_vline(xintercept = gene_cutoff_high, linetype='dashed', color='red')
ggsave('figures/qc/gene_counts_density.png', dpi=320, width=10, height= 10)


## plot count depth (number of reads) vs number of genes, colored by sample name (left) and percent mitochondrial gene counts (right) ----
p1 = FeatureScatter(merged_seurat, feature1='nCount_RNA', feature2='nFeature_RNA', group.by='sample', pt.size=0.5, shuffle=TRUE)+
  DarkTheme()+
  labs(title='Count depth vs # of genes',
       subtitle='colored by sample',
       x='# of UMI counts (log scale)',
       y='# of genes (log scale)',
       color='sample')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm, linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
  geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='blue')+
  geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='blue')

p2 = ggplot(metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
  geom_point(size=0.5)+
  DarkTheme()+
  labs(title='Count depth vs # of genes',
       subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
       x='# of UMI counts (log scale)',
       y='# of genes (log scale)',
       color='% mt counts')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
  geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='blue')+
  geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='blue')+
  scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#292a2d')

plot_grid(p1, p2)
ggsave('figures/qc/ncounts_vs_ngenes_coloredby_percentmt.png', dpi=320, width=20, height= 10)


## plot count depth (number of reads) vs number of genes, colored by exposure condition and mitochondrial gene percentage ----
p1 = FeatureScatter(merged_seurat, feature1='nCount_RNA', feature2='nFeature_RNA', group.by='condition', pt.size=0.5, shuffle=TRUE)+
  DarkTheme()+
  labs(title='Count depth vs # of genes',
       subtitle='colored by condition',
       x='# of UMI counts (log scale)',
       y='# of genes (log scale)',
       color='condition')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm, linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
  geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='blue')+
  geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='blue')

p2 = ggplot(metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
  geom_point(size=0.5)+
  DarkTheme()+
  labs(title='Count depth vs # of genes',
       subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
       x='# of UMI counts (log scale)',
       y='# of genes (log scale)',
       color='% mt counts')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
  geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='blue')+
  geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='blue')+
  scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#292a2d')
  

plot_grid(p1, p2)
ggsave('figures/qc/ncounts_vs_ngenes_coloredby_percentmt_and_exposure_condition.png', dpi=320, width=20, height= 10)

# > length(metadata$orig.ident)
# [1] 46092

# > length(metadata[metadata$percent.mt < 10,]$orig.ident)
# [1] 46009

# > length(metadata[metadata$percent.mt < 5,]$orig.ident)
# [1] 45530. at 5% cutoff ~550 cells filtered out

length(metadata[(metadata$percent.mt < 5) & ((metadata$nCount_RNA > 150) & (metadata$nCount_RNA < 50000)) & 
  ((metadata$nFeature_RNA > 150) & (metadata$nFeature_RNA < 6000)) & (metadata$log10GenesPerUMI > 0.8),]$orig.ident)

## more vln plots for UMIs, genes, and percent.mt by exposure condition ----            
QC_Plots_Genes(merged_seurat, group.by = 'condition', low_cutoff = gene_cutoff_low, high_cutoff = gene_cutoff_high)
ggsave('figures/qc/genes_per_cell_by_exposure_condition.png', dpi=320, width=10, height= 10)

QC_Plots_UMIs(merged_seurat, group.by = 'condition', low_cutoff = umi_cutoff_low, high_cutoff = umi_cutoff_high, y_axis_log = TRUE)+
  scale_y_log10(n.breaks=10, labels=comma_format())+
  annotation_logticks(sides = "l", size=.5)
ggsave('figures/qc/UMIs_per_cell_by_exposure_condition.png', dpi=320, width=10, height= 10)

QC_Plots_Mito(merged_seurat, mito_name='percent.mt',group.by = 'condition', low_cutoff = percent.mt.cutoff)
ggsave('figures/qc/percent.mt_per_cell_by_exposure_condition.png', dpi=320, width=10, height= 10)



# # calculate log transformed median, median absolute deviation (MAD), and +- 5x MAD
mads_df = data.frame(median_gene=numeric(12), MAD_gene=numeric(12), plus_5MAD_gene=numeric(12), minus_5MAD_gene=numeric(12),
                     median_count=numeric(12), MAD_count=numeric(12), plus_5MAD_count=numeric(12), minus_5MAD_count=numeric(12))
for (i in 1:12) {
  s = unique(metadata$sample)[i]
  MAD_gene = mad(metadata[metadata$sample == s,]$nFeature_RNA)
  MAD_count = mad(metadata[metadata$sample == s,]$nCount_RNA)
  median_gene = median(metadata[metadata$sample == s,]$nFeature_RNA)
  median_count = median(metadata[metadata$sample == s,]$nCount_RNA)
  plus_5MAD_gene = median_gene + 5*MAD_gene
  minus_5MAD_gene = median_gene - 5*MAD_gene
  plus_5MAD_count = median_count + 5*MAD_count
  minus_5MAD_count = median_count - 5*MAD_count

  mads_df[i,] = list(median_gene, MAD_gene, plus_5MAD_gene, minus_5MAD_gene, median_count, MAD_count, plus_5MAD_count, minus_5MAD_count)
}
saveRDS(mads_df, 'data/mads_df.rds')

#### Plot log transformed coutns vs genes with average 3x MAD cutoff ----
plus_5MAD_gene = mean(mads_df$plus_5MAD_gene)
minus_5MAD_gene = mean(mads_df$minus_5MAD_gene)
plus_5MAD_count = mean(mads_df$plus_5MAD_count)
minus_5MAD_count = mean(mads_df$minus_5MAD_count)

ggplot(metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
  geom_point(size=0.5)+
  DarkTheme()+
  labs(title='Count depth vs # of genes',
       subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
       x='log(# of UMI counts)',
       y='log(# of genes)',
       color='% mt counts')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=minus_5MAD_gene, linetype='dashed',color='red')+
  geom_hline(yintercept=plus_5MAD_gene, linetype='dashed',color='red')+
  annotate('text', x=200, y=1700, label=round(plus_5MAD_gene,0), size=5, color='red')+
  geom_vline(xintercept=minus_5MAD_count, linetype='dashed',color='blue')+
  geom_vline(xintercept=plus_5MAD_count, linetype='dashed',color='blue')+
  annotate('text', x=2500, y=200, label=round(plus_5MAD_count,0), size=5, color='blue')+
  scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(10, NA), na.value ='#292a2d')
ggsave('figures/qc/ncounts_vs_ngenes_using_median_absolute_deviation_cutoffs.png', dpi=320, width=10, height= 10)


count93 = quantile(merged_seurat@meta.data$nCount_RNA, 0.93) # calculate value in the 93rd percentile -- drawn from 10xGenomics doublet rate estimates

## Remove poor quality samples 29 and 32 ----
filtered_seurat = subset(merged_seurat, idents=c('count29','count32'),invert=TRUE)

length(filtered_seurat$nCount_RNA[filtered_seurat$sample == 'ctrl24'])
# 1720

length(filtered_seurat$nCount_RNA[filtered_seurat$sample == 'pb31'])
# 2096


filtered_seurat
# An object of class Seurat 
# 19059 features across 44314 samples within 1 assay 
# Active assay: RNA (19059 features, 0 variable features)

## Filter samples 24 and 31 individually ----
## since they still contain some background signal according to the
## cellranger summary files
filtered_seurat = subset(filtered_seurat, 
                        subset=(orig.ident == 'count24' &
                                  nCount_RNA >= 558) | 
                          (orig.ident != 'count24'))

length(filtered_seurat$nCount_RNA[filtered_seurat$sample == 'ctrl24'])
# 1694

filtered_seurat
# An object of class Seurat 
# 19059 features across 44288 samples within 1 assay 
# Active assay: RNA (19059 features, 0 variable features)

filtered_seurat = subset(filtered_seurat, 
                         subset=(orig.ident == 'count31' &
                                   nCount_RNA >= 568) | 
                           (orig.ident != 'count31'))


length(filtered_seurat$nCount_RNA[filtered_seurat$sample == 'pb31'])
# 1990

filtered_seurat
# An object of class Seurat 
# 19059 features across 44182 samples within 1 assay 
# Active assay: RNA (19059 features, 0 variable features)


# number of cell before and after sample-level filtering
ncells = length(merged_seurat$orig.ident)
ncells_post_cell_filter = length(filtered_seurat$orig.ident)
cat(paste0(ncells,' cells before sample-level filtering\n',
           ncells_post_cell_filter,' cells left after sample-level filtering\n',
           ncells - ncells_post_cell_filter,' cells removed after sample-level filtering'))
# 46092 cells before sample-level filtering
# 44182 cells left after sample-level filtering
# 1910 cells removed after sample-level filtering



## cell-level filtering of all remaining samples ----
ncells = length(filtered_seurat$orig.ident) # for tracking # cells before filtering

filtered_seurat = subset(x = filtered_seurat, 
                          subset=(nCount_RNA >= umi_cutoff_low) & 
                            (nCount_RNA <= umi_cutoff_high) & 
                            (nFeature_RNA >= gene_cutoff_low) &
                            (nFeature_RNA <= gene_cutoff_high) &
                            (log10GenesPerUMI > 0.80) & 
                            (percent.mt < percent.mt.cutoff))

# number of cell before and after cell-level filtering
ncells_post_cell_filter = length(filtered_seurat$orig.ident)
cat(paste0(ncells,' cells before cell-level filtering\n',
           ncells_post_cell_filter,' cells left after cell-level filtering\n',
           ncells - ncells_post_cell_filter,' cells removed after cell-level filtering'))
# 44182 cells before cell-level filtering
# 44102 cells left after cell-level filtering
# 80 cells removed after cell-level filtering

## gene-level filtering ----
counts = filtered_seurat[['RNA']]$counts

# extract counts from joint layers
jointcounts = JoinLayers(filtered_seurat)
counts = LayerData(jointcounts, assay='RNA', layer='counts')

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero = counts > 0

# Sums all TRUE values and returns TRUE if more than 20 TRUE values per gene
keep_genes = Matrix::rowSums(nonzero) >= 20

# Only keeping those genes expressed in more than 20 cells
filtered_counts = counts[keep_genes, ]

# Reassign to filtered Seurat object and filtered metadata
filtered_seurat = CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
saveRDS(filtered_seurat, 'data/filtered_seurat_rm29and32.rds')

filtered_metadata = filtered_seurat@meta.data
saveRDS(filtered_metadata, 'data/filtered_metadata_rm29and32.rds')

rm(counts, filtered_counts, jointcounts, nonzero) # clean up RAM
gc()

## number of cells after filtering ----
ncells_post_gene_filter = length(filtered_seurat$orig.ident)
cat(paste0(ncells,' cells before cell-level filtering\n',
           ncells - ncells_post_cell_filter,' cells removed after cell-level filtering\n',
           ncells_post_cell_filter - ncells_post_gene_filter,' more cells removed after gene-level filtering\n',
           ncells - ncells_post_gene_filter,' total cells removed after cell and gene-level filtering\n',
           ncells_post_gene_filter, ' cells after filtering'))
# 44182 cells before cell-level filtering
# 80 cells removed after cell-level filtering
# 0 more cells removed after gene-level filtering
# 80 total cells removed after cell and gene-level filtering
# 44102 cells after filtering


## plot number of cells before and after filtering dataset ----
p1 = ggplot(metadata, aes(x=sample, fill=sample)) + 
      geom_bar() +
      theme_classic2() +
      NoLegend()+
      labs(title="Number of cells per experiment",
           subtitle = paste0('Total cells = ',dim(metadata)[1],', # ctrl cells = ',
                             sum(metadata$condition =='ctrl'),', # pb cells = ',
                             sum(metadata$condition == 'pb')))+
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      scale_y_continuous(n.breaks=8)
    
p2 = ggplot(filtered_metadata, aes(x=sample, fill=sample)) + 
      geom_bar() +
      theme_classic2() +
      NoLegend()+
      labs(title="Number of cells per experiment after filtering",
           subtitle = paste0('Total cells = ',dim(filtered_metadata)[1],', # ctrl cells = ',
                             sum(filtered_metadata$condition =='ctrl'),', # pb cells = ',
                             sum(filtered_metadata$condition == 'pb')))+
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            plot.subtitle = element_text(hjust=0.5),
            axis.title = element_text(size=12), 
            axis.text = element_text(size=12),
            title = element_text(size=18))+
      scale_y_continuous(n.breaks=8)

plot_grid(p1, p2)
ggsave('figures/qc/ncells_before_and_after_filtering_rm29and32.png', dpi=320, width=20, height= 10)


## plot count depth (number of UMIs) vs number of genes, colored by percent mitochondrial gene counts, before and after filtering ----
p1 = ggplot(metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
  geom_point(size=0.5)+
  DarkTheme()+
  labs(title='Count depth vs # of genes',
       subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
       x='# of UMI counts (log scale)',
       y='# of genes (log scale)',
       color='% mt counts')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
  geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='blue')+
  geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='blue')+
  scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#4b4d53')

p2 = ggplot(filtered_metadata,aes(x=nCount_RNA,y=nFeature_RNA,color=percent.mt))+
  geom_point(size=0.5)+
  DarkTheme()+
  labs(title='Count depth vs # of genes after filtering',
       subtitle=paste0('colored by percent mitochondrial (mt) counts >', percent.mt.cutoff),
       x='# of UMI counts (log scale)',
       y='# of genes (log scale)',
       color='% mt counts')+
  scale_x_log10(n.breaks=8, labels=comma_format())+
  scale_y_log10(n.breaks=8, limits=c(NA,10000), labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  stat_smooth(method=lm,linewidth=0.5, colour='#61e4e8')+
  geom_hline(yintercept=gene_cutoff_low, linetype='dashed',color='red')+
  geom_hline(yintercept=gene_cutoff_high, linetype='dashed',color='red')+
  geom_vline(xintercept=umi_cutoff_low, linetype='dashed',color='blue')+
  geom_vline(xintercept=umi_cutoff_high, linetype='dashed',color='blue')+
  scale_color_gradientn(colors=viridis_plasma_light_high, limits = c(percent.mt.cutoff, NA), na.value ='#4b4d53')

plot_grid(p1, p2)
ggsave('figures/qc/ncounts_vs_ngenes_coloredby_percentmt_before_and_after_filtering_rm29and32.png', dpi=320, width=20, height= 10)










#### misc stuff like MAD and ave count of genes




median(merged_seurat$nCount_RNA)

mad(merged_seurat$nCount_RNA) * 5


c = median(merged_seurat$nCount_RNA) + mad(merged_seurat$nCount_RNA) * 5
c # 1921

g = median(merged_seurat$nFeature_RNA) + mad(merged_seurat$nFeature_RNA) * 5
g # 1529

m = median(merged_seurat$percent.mt) + mad(merged_seurat$percent.mt) * 5
m # 4.47

i = WhichCells(merged_seurat, expression=(nCount_RNA >= 3034) |
             (nFeature_RNA >= 2946) | (percent.mt >= 4.72) | (log10GenesPerUMI <= 0.8))
i


samples = unique(metadata$orig.ident)
samples
# [1] "count21" "count22" "count23" "count24" "count25" "count26" "count27" "count28" "count29" "count30" "count31" "count32"

metadata = readRDS('data/metadata.rds')
metadata = metadata[!(metadata$orig.ident == 'count29' | metadata$orig.ident == 'count32'),]

samples = unique(metadata$orig.ident)
samples
# [1] "count21" "count22" "count23" "count24" "count25" "count26" "count27" "count28" "count29" "count30" "count31" "count32"


total_cells_filtered = vector(mode='numeric',length=12) # assuming use all 12 samples
cells_filtered_by_sample = vector(mode='numeric',length=12)
for (i in seq_along(samples)) {
  data = metadata[metadata$orig.ident == samples[i],]
  
  c = median(data$nCount_RNA) + 5*mad(data$nCount_RNA)
  
  g = median(data$nFeature_RNA) + 5*mad(data$nFeature_RNA)
  
  m = median(data$percent.mt) + 5*mad(data$percent.mt)

  t = WhichCells(subset(merged_seurat,idents=samples[i]), expression=(nCount_RNA >= c) |
                   (nFeature_RNA >= g) | (percent.mt >= 10) | (log10GenesPerUMI <= 0.8))
  
  cat(paste0(i, '\n',
             'MADx5 nCount_RNA: ',c,'\n',
             'MADx5 nFeature_RNA: ',g,'\n',
             'MADx5 percent.mt: ',m,'\n',
             '# of cells filtered: ',length(t),'\n',
             '\n'))
  total_cells_filtered[i] = length(t)
  cells_filtered_by_sample[i] = c(t)
  names(cells_filtered_by_sample[i]) = samples[i]
}
cat(paste0('total # of cells filtered using 5xMAD rule for each sample,
           not including samples 29 and 32: ', sum(total_cells_filtered)))
# 1724


cells_filtered_by_sample

data = subset(merged_seurat,idents='count21')
c = median(data$nCount_RNA) + 5*mad(data$nCount_RNA)

g = median(data$nFeature_RNA) + 5*mad(data$nFeature_RNA)

m = median(data$percent.mt) + 5*mad(data$percent.mt)

t = WhichCells(subset(merged_seurat,idents=samples[i]), expression=(nCount_RNA >= c) |
                 (nFeature_RNA >= g) | (percent.mt >= 10) | (log10GenesPerUMI <= 0.8))



ggplot(metadata, aes(x=condition,y=nFeature_RNA,fill=sample)) +
  geom_violin(draw_quantiles=0.5)+
  scale_y_log10() +
  ylab("Total number of genes detected") +
  ggtitle("Number of Genes detected") +
  theme_bw() +
  scale_fill_manual(values=p12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


View(JoinLayers(merged_seurat))
AverageExpression(JoinLayers(merged_seurat), layer='counts')

joined = JoinLayers(merged_seurat)

joined@assays$RNA$counts

ave.counts = rowMeans(joined@assays$RNA$counts)

hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)



#### sessionInfo() ----
sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22631)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Polychrome_1.5.1   RCurl_1.98-1.16    metap_1.11         glmGamPoi_1.14.3   scales_1.3.0       viridis_0.6.5      viridisLite_0.4.2  plyr_1.8.9         ggpubr_0.6.0      
# [10] ggthemes_5.1.0     Matrix_1.6-5       lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.2        readr_2.1.5        tidyr_1.3.1       
# [19] tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0    cowplot_1.1.3      patchwork_1.3.0    scCustomize_2.1.2  Seurat_5.1.0       SeuratObject_5.0.2 sp_2.1-4          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                polyclip_1.10-7             janitor_2.2.0              
# [7] fastDummies_1.7.4           lifecycle_1.0.4             Rdpack_2.6.1                rstatix_0.7.2               globals_0.16.3              lattice_0.22-6             
# [13] MASS_7.3-60.0.1             backports_1.5.0             magrittr_2.0.3              plotly_4.10.4               plotrix_3.8-4               qqconf_1.3.2               
# [19] httpuv_1.6.15               sn_2.1.1                    sctransform_0.4.1           spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0          
# [25] pbapply_1.7-2               RColorBrewer_1.1-3          multcomp_1.4-26             abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                 
# [31] GenomicRanges_1.54.1        BiocGenerics_0.48.1         TH.data_1.1-2               sandwich_3.1-1              circlize_0.4.16             GenomeInfoDbData_1.2.11    
# [37] IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6               irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-0       
# [43] TFisher_0.2.0               goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-2       fitdistrplus_1.2-1          parallelly_1.38.0          
# [49] DelayedArray_0.28.0         leiden_0.4.3.1              codetools_0.2-20            tidyselect_1.2.1            shape_1.4.6.1               farver_2.1.2               
# [55] matrixStats_1.4.1           stats4_4.3.1                spatstat.explore_3.3-2      mathjaxr_1.6-0              jsonlite_1.8.9              multtest_2.58.0            
# [61] progressr_0.14.0            Formula_1.2-5               ggridges_0.5.6              survival_3.7-0              tools_4.3.1                 ica_1.0-3                  
# [67] Rcpp_1.0.13                 glue_1.7.0                  mnormt_2.1.1                SparseArray_1.2.4           gridExtra_2.3               MatrixGenerics_1.14.0      
# [73] GenomeInfoDb_1.38.8         numDeriv_2016.8-1.1         withr_3.0.1                 fastmap_1.2.0               fansi_1.0.6                 digest_0.6.37              
# [79] timechange_0.3.0            R6_2.5.1                    mime_0.12                   ggprism_1.0.5               colorspace_2.1-1            scattermore_1.2            
# [85] tensor_1.5                  spatstat.data_3.1-2         utf8_1.2.4                  generics_0.1.3              data.table_1.16.0           httr_1.4.7                 
# [91] htmlwidgets_1.6.4           S4Arrays_1.2.1              scatterplot3d_0.3-44        uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.5               
# [97] lmtest_0.9-40               XVector_0.42.0              htmltools_0.5.8.1           carData_3.0-5               dotCall64_1.1-1             Biobase_2.62.0             
# [103] png_0.1-8                   spatstat.univar_3.0-1       snakecase_0.11.1            rstudioapi_0.16.0           tzdb_0.4.0                  reshape2_1.4.4             
# [109] nlme_3.1-166                zoo_1.8-12                  GlobalOptions_0.1.2         KernSmooth_2.23-24          parallel_4.3.1              miniUI_0.1.1.1             
# [115] vipor_0.4.7                 ggrastr_1.0.2               pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                 RANN_2.6.2                 
# [121] promises_1.3.0              car_3.1-3                   xtable_1.8-4                cluster_2.1.6               beeswarm_0.4.0              paletteer_1.6.0            
# [127] mvtnorm_1.3-1               cli_3.6.3                   compiler_4.3.1              crayon_1.5.3                rlang_1.1.4                 mutoss_0.1-13              
# [133] future.apply_1.11.2         ggsignif_0.6.4              rematch2_2.1.2              ggbeeswarm_0.7.2            stringi_1.8.4               deldir_2.0-4               
# [139] munsell_0.5.1               lazyeval_0.2.2              spatstat.geom_3.3-3         RcppHNSW_0.6.0              hms_1.1.3                   future_1.34.0              
# [145] shiny_1.9.1                 SummarizedExperiment_1.32.0 rbibutils_2.2.16            ROCR_1.0-11                 igraph_2.0.3                broom_1.0.7