setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(SeuratObject)
library(sctransform)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(glmGamPoi)
library(metap)
library(RCurl)
library(Polychrome)
library(pheatmap)
library(presto)
library(DoubletFinder)

# set seed
set.seed(2024)

# sessionInfo()

# load in filtered data
filtered_seurat = readRDS('data/filtered_seurat_rm29and32.rds')


### create polychrome color palette ----
p12 = createPalette(12, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
swatch(p12)
names(p12) = NULL

### exploring sources of variation, cell cycle and mitochondrial gene expression----

# cell cycle genes, formatted to match genes in matrix with first letter upper-cased
s.genes = stringr::str_to_title(cc.genes$s.genes)
g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)


# Normalize the counts
seurat_phase = NormalizeData(filtered_seurat)
rm(filtered_seurat)# free RAM
gc()

# Score cells for cell cycle
seurat_phase = CellCycleScoring(seurat_phase, 
                                g2m.features = g2m.genes, 
                                s.features = s.genes,
                                set.ident=T)

# From https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Seurat suggests regressing out the difference between S and G2M phase scores
# in tissues with differentiating cells, like our developing mammary gland.
seurat_phase$cc.difference = seurat_phase$S.Score - seurat_phase$G2M.Score


# Check quartile values of mitochondrial gene expression
summary(seurat_phase@meta.data$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3536  0.7692  1.0381  1.3996  9.8712 

# Turn percent.mt into categorical factor vector based on quartile values
seurat_phase@meta.data$mito.level = cut(seurat_phase@meta.data$percent.mt, 
                                        breaks=c(-Inf, 0.3610, 0.7853, 1.0592, Inf), 
                                        labels=c("1st","2nd","3rd", "4th"))

# Identify the 25 most highly variable genes
top25 = head(VariableFeatures(seurat_pca), 25)
# plot variable features with and without labels
p = VariableFeaturePlot(seurat_pca)
LabelPoints(plot=p, points = top25, repel = TRUE)+
  labs(title='Top 25 variably expressed genes')+
  theme_classic2()
ggsave('figures/clustering_rm29and32/top25_variable_genes_seurat_phase.png', dpi=320, width=20, height= 10)


# Identify the most variable genes
seurat_pca = FindVariableFeatures(seurat_phase, 
                                  selection.method = "vst", # default
                                  nfeatures = 2000, # default
                                  verbose = FALSE)

# Scale the counts
seurat_pca = ScaleData(seurat_pca)

### Perform PCA on merged non-integrated data using cell cycle genes ----
seurat_pca = RunPCA(seurat_pca, features=c(s.genes,g2m.genes), npcs=50) # default params

#### Plot the PCA colored by cell cycle phase using cell cycle genes ----
DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "Phase",)+
  labs(title='PCA by cell cycle phase',
       color='phase')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32/PCA_by_cellcycle_phase_genes.png', dpi=320, width=20, height= 10)
# there is a lot of variation between S and G2M phase. 
# Perhaps better to regress this out in our developing mammary tissue


# From https://satijalab.org/seurat/articles/cell_cycle_vignette.html
# Regress out the difference between S and G2M phase scores
# takes a long while to regress out. Pushing up on 20-30 mins
seurat_pca = ScaleData(seurat_pca, vars.to.regress = "cc.difference", features = rownames(seurat_pca))
saveRDS(seurat_pca, 'data/seurat_phase_rm29and32_regressout_cc.diff.rds')

# run PCA again w/ cc.diff now regressed out
seurat_pca = RunPCA(seurat_pca, features=c(s.genes,g2m.genes), npcs=50) # default params

## What affect did regressing out the cc.diff have?
DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "Phase",
        split.by='Phase')+
  labs(title='PCA by cell cycle phase',
       color='phase')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32/PCA_by_cellcycle_phase_genes_ccdiff_regressed_out.png', dpi=320, width=20, height= 10)
# still some variability between phases

# run PCA on all genes with cc.diff regressed out
seurat_pca = RunPCA(seurat_pca, npcs=50) # default params

DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "Phase",
        split.by='Phase',
        shuffle=T)+
  labs(title='PCA by cell cycle phase',
       color='phase')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32/PCA_by_cellcycle_phase_allgenes_.png', dpi=320, width=20, height= 10)


### Plot the PCA colored by mitochondrial expression quartile ----
DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "mito.level",
        split.by = "mito.level")+
  labs(title='PCA by mitochondrial expression quartile',
       color='quartile')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32/PCA_by_mitochondrial_expression_quartile.png', dpi=320, width=20, height= 10)

saveRDS(seurat_pca, 'data/seurat_phase_rm29and32_regressout_cc.diff.rds')# save seurat_phase object
##looking at the plots there doesn't seem to be any significant source of variation from cell cycle phase or mitochondrial expression

### plot the PCA colored by condition ----
DimPlot(seurat_pca,
        reduction = "pca",
        group.by= "condition",
        split.by='condition')+
  labs(title='PCA by condition',
       color='condition')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32/PCA_by_condition.png', dpi=320, width=20, height= 10)


rm(seurat_pca) # free RAM. Don't need any further

### Split seurat_phase by sample and apply sctransform to each individually ----
# Split seurat object by sample to perform cell cycle scoring and SCTransform on all samples
## will normalize seurat_phase with SCTransform as don't need PCA exploration from seurat_pca object.
split_seurat = SplitObject(seurat_phase, split.by = "sample")

rm(seurat_phase) # save RAM
gc()

#### doublet finder ----
# adjust object size limits of variables/objects in R
options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB

# initiate vectors to store optimal parameters for each sample
pks = c() 
nexps = c()
# pks = c(0.05,0.01,0.03,0.27,0.01,0.01,0.05,0.01,0.01,0.03)

#### apply SCTransform to each of the 10 samples, then merge normalized samples and save object ----
# takes a few hours to run. Do for each sample separately according to https://github.com/satijalab/seurat/issues/6003
# return all genes/features according to https://github.com/satijalab/seurat/issues/5205
# also run doublet finder for each sample, setting sct = TRUE.

for(i in seq_along(split_seurat)) { # Takes a few hours to run
  cat('\n Beginning sample:', names(split_seurat[i]), '\n\n')
  # split_seurat[i] = NormalizeData(split_seurat[i])
  # split_seurat[i] = FindVariableFeatures(split_seurat[i])
  # split_seurat[i] = ScaleData(split_seurat[i])
  split_seurat[i] = SCTransform(split_seurat[[i]], assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                  vars.to.regress=c('cc.difference','percent.mt'),
                  return.only.var.genes=FALSE)
  cat('/n running PCA for,', names(split_seurat[i]), '\n\n')
  split_seurat[i] = RunPCA(split_seurat[[i]], assay='SCT', npcs=50)
  split_seurat[i] = FindNeighbors(split_seurat[[i]], dims=1:50)
  split_seurat[i] = FindClusters(object = split_seurat[[i]])
  split_seurat[i] = RunUMAP(split_seurat[[i]], dims=1:50, reduction='pca')
  
  # estimate pK, the neighborhood size used to compute
  # the number of artificial doublets
  cat('\n Running doublet finder parameter sweep for',names(split_seurat[i]), '\n\n')
  sweep.res.list = paramSweep(split_seurat[[i]], PCs = 1:50, sct = TRUE)
  sweep.stats = summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn = find.pK(sweep.stats)

  ggplot(bcmvn, aes(pK, BCmetric, group=1))+
    geom_point()+
    geom_line()+
    labs(title=paste0(names(split_seurat[i]),' pK vs bcmvn parameter sweep'))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          plot.title = element_text(hjust=0.5, face="bold"),
          axis.title = element_text(size=12),
          axis.text = element_text(size=12),
          title = element_text(size=18))+
    theme_classic2()
  ggsave(paste0('figures/qc/doublet_finder/',names(split_seurat[i]),'_doublet_finder_pKvsbcmvn.png'),dpi=320, width=20, height=20)

  pK = bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK = as.numeric(as.character(pK[[1]]))

  pks[i] = pK
  cat(pK)
  cat(pks)
  pK = pks[i]
  
  ## Homotypic Doublet Proportion Estimate
  cat('\n estimating homotypic doublet proportions for', names(split_seurat[i]), '\n\n')
  annotations = split_seurat[[i]]@meta.data$seurat_clusters
  homotypic.prop = modelHomotypic(annotations) #ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi = round(0.0008*nrow(split_seurat[[i]]@meta.data)) # Assuming 7.5% doublet formation rate - tailor for your dataset. 10x said 0.8 per thousand for flex method
  nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
  
  nexps[i] = nExp_poi.adj
  cat(nExp_poi.adj)
  cat(nexps)
  
  ## Run DoubletFinder on sample
  cat('\n Running doublet finder with optimized parameters for', names(split_seurat[i]), '\n ')
  split_seurat[i] = doubletFinder(split_seurat[[i]], PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)

}

# save doublet finder parameters used for each sample
names(pks) = unique(seurat_merged_sct$sample)
pks
# ctrl21 ctrl22 ctrl23 ctrl24 ctrl25 ctrl26   pb27   pb28   pb30   pb31 
# 0.05   0.01   0.03   0.27   0.01   0.01   0.05   0.01   0.01   0.03 

nexps
# [1] 2 2 3 1 4 4 5 3 4 2

pns = rep(0.25, 10)

params_doubletfinder = data.frame(pns,pks,nexps)
write.csv(params_doubletfinder, 'data/doublet_finder_params_3week_mice_mammarygland.csv')

saveRDS(split_seurat, 'data/split_seurat__Sct_rm29and32_regressout_cc.diffandmt_doubletfinder.rds')


# adjust object size limits of variables/objects in R
# options(future.globals.maxSize = 4 * 1024^3) # set limit to 4 GB

#### apply SCTransform to each of the 10 samples, then merge normalized samples and save object ----
# takes a while to run. Do for each sample separately according to https://github.com/satijalab/seurat/issues/6003
# return all genes/features according to https://github.com/satijalab/seurat/issues/5205
# split_seurat = lapply(split_seurat, function(x) {
#   x = SCTransform(x, assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
#                   vars.to.regress=c('cc.difference','percent.mt'),
#                   return.only.var.genes=FALSE)
# })

# saveRDS(split_seurat, 'data/split_seurat_rm29and32_regressout_cc.diffandmt.rds')


# merge sctransform normalized samples. Takes a little while to run
seurat_merged_sct = merge(split_seurat[[1]], split_seurat[-1], merge.data=TRUE,merge.dr=FALSE)

### merge doublet classification scores and classification metadata fields ----
# into one since they were saved separately for each set of parameters
seurat_merged_sct$doublet_class = apply(seurat_merged_sct@meta.data[, grep('DF.classifications_', colnames(seurat_merged_sct@meta.data))], MARGIN=1, FUN=na.omit)
seurat_merged_sct$doublet_score = apply(seurat_merged_sct@meta.data[, grep('pANN_', colnames(seurat_merged_sct@meta.data))], MARGIN=1, FUN=na.omit)

# remove individual DF.classification and pANN cols
seurat_merged_sct@meta.data[, grep('DF.classifications_', colnames(seurat_merged_sct@meta.data))] = NULL
seurat_merged_sct@meta.data[, grep('pANN_', colnames(seurat_merged_sct@meta.data))] = NULL

### subset doublets out of dataset ----

# tabulate doublets by sample
table(seurat_merged_sct$sample, seurat_merged_sct$doublet_class)
#         Doublet Singlet
# ctrl21       2    2135
# ctrl22       2    3028
# ctrl23       3    4560
# ctrl24       1    1683
# ctrl25       4    6122
# ctrl26       4    5662
# pb27         5    7962
# pb28         3    4391
# pb30         4    6545
# pb31         2    1984

# tabulate total doublets
table(seurat_merged_sct$doublet_class)
# Doublet Singlet 
# 30   44072    ## 44102 cells total

x = split_seurat$pb27
p1 = DimPlot(x,
        group.by=grep('DF.classifications_', colnames(x@meta.data), value=TRUE),
        reduction='umap',
        shuffle=TRUE)
p2 = DimPlot(x,
             reduction='umap',
             shuffle=TRUE)
p1+p2

# subset only singlets
seurat_merged_sct = subset(seurat_merged_sct, subset=doublet_class == 'Singlet')

saveRDS(seurat_merged_sct, 'data/seurat_merged_sct_rm29and32_regressout_cc.diffandmt_doubletfinder.rds')

gc() #free up RAM

### continue analysis without integration using sctransform normalization ----
# select integration features according to https://github.com/satijalab/seurat/issues/
# select 5000 features to account for the possible different rankings of highly variable genes
# between samples. Since we're normalizing with SCTtransform this # shouldn't change things much.
features = SelectIntegrationFeatures(split_seurat, nfeatures=5000, fvf.nfeatures=5000)
saveRDS(features,'data/variable_features.rds')
VariableFeatures(seurat_merged_sct) = features

sct_results = SCTResults(seurat_merged_sct, slot = "feature.attributes")
saveRDS(sct_results, 'figures/clustering_rm29and32/sct_results_rm29and32_regressout_cc.diffandmt_doubletfinder.rds')
rm(sct_results)

# Variable features lost/unset after merging transformed objects. Need to set again.
# should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(seurat_merged_sct[["SCT"]]) = rownames(seurat_merged_sct[["SCT"]]@scale.data)

# can plot variable genes by plotting gmean vs residual variance
# according to https://github.com/satijalab/seurat/issues/8321
# z = seurat_merged_sct@assays$SCT@SCTModel.list[[3]]@feature.attributes
# qplot(z$gmean, z$residual_variance) + scale_y_log10() + scale_x_log10()

# Identify the 25 most highly variable genes
top25_sct = head(VariableFeatures(seurat_merged_sct), 25)
top25_sct

sum(top25 %in% top25_sct) # = 0. completely different top 25 variable features between sctransform and log norm

# VariableFeatures() won't work for merged object after SCTtransform.
# Need to run for each sample individually
lapply(seq_along(split_seurat), function(x,n,i) {
  top25_sct = head(VariableFeatures(x[[i]]), 25)
  p = VariableFeaturePlot(x[[i]], selection.method='sct', assay='SCT') 
  LabelPoints(plot=p, points = top25_sct, repel = TRUE)+
    labs(title=paste0('Top 25 variably expressed genes for ', n[[i]]))+
    theme_classic2()
  ggsave(paste0('figures/clustering_rm29and32_regress_ccdiffandmt/top25_variable_genes_',n[[i]],'.png'), dpi=320, width=20, height= 10)
  
}, x=split_seurat, n=names(split_seurat))

rm(split_seurat) # free RAM
gc()

### Run PCA on non-integrated merged SCTransformed data ---- 
seurat_merged_sct = RunPCA(seurat_merged_sct, assay='SCT', npcs=50) # default params
gc() # free up RAM

saveRDS(seurat_merged_sct, 'data/seurat_merged_sct_rm29and32_regressout_cc.diffandmt_doubletfinder.rds')


### Visualize expression of most highly weighted genes for PCs 1-50 (loadings) ----
VizDimLoadings(seurat_merged_sct, dims = 1:12, ncol=4)
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/PCA_loadings_sct_1-12.png', dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 13:24, ncol=4)
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/PCA_loadings_sct_13-24.png', dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 25:36, ncol=4)
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/PCA_loadings_sct_25-36.png', dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 37:43, ncol=4)
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/PCA_loadings_sct_37-43.png', dpi=320, width=20, height= 15)

VizDimLoadings(seurat_merged_sct, dims = 44:50, ncol=4)
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/PCA_loadings_sct_44-50.png', dpi=320, width=20, height= 15)

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC1-6_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 1:6,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC7-12_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 7:12,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC13-18_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 13:18,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC19-24_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 19:24,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC25-30_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 25:30,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC31-36_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 31:36,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC37-42_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 37:42,
           ncol=3,
           cells=500,
           balanced = TRUE)
dev.off()

png(filename = 'figures/clustering_rm29and32_regress_ccdiffandmt/sct/high_score_genes_per_PC43-50_heatmap_sct.png', units='in', width=20, height=10, res=320)
DimHeatmap(object = seurat_merged_sct, 
           reduction = "pca", 
           dims = 43:50,
           ncol=4,
           cells=500,
           balanced = TRUE)
dev.off()

# ranked list (PC score) of genes for PC1
names(sort(Loadings(seurat_merged_sct, reduction="pca")[,1], decreasing=TRUE))


# Printing out the top 10 most variable genes driving PCs
print(x = seurat_merged_sct[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

### Plot the elbow plot to see variation captured by PCs ----
ElbowPlot(object = seurat_merged_sct, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/elbow_plot_sct.png', dpi=320, width=20, height= 10)

# We can calculate where the principal components start to elbow by taking the smaller value of:
#
# 1)The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# 2)The point where the percent change in variation between the consecutive PCs is less than 0.1%.

### Determine percent of variation associated with each PC ----
pct = seurat_merged_sct[["pca"]]@stdev / sum(seurat_merged_sct[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 42
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 12
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 
# [1] 12

# Create a dataframe to visualize percent variation captured by PCs
plot_df = data.frame(pct = pct, 
                     cumu = cumu, 
                     rank = 1:length(pct))

### Elbow plot to visualize percent variation captured by PCs ----
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_classic2()+
  scale_x_continuous(n.breaks=10)+
  scale_y_continuous(n.breaks=10)+
  labs(title = 'Eblow Plot of percent variation explained by PCs',
       x='Cumulative percentage',
       y='Percent of variation')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/quantitative_elbow_plot_sct.png', dpi=320, width=20, height= 10)
rm(plot_df)

### cluster sctransformed non-integrated cells----
seurat_umap_sct = FindNeighbors(seurat_merged_sct, dims=1:50)
seurat_umap_sct = FindClusters(object = seurat_umap_sct,
                           resolution = c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 2.0))

## clusters after running with dims=1:50 filtered nUMI<10000, removed low quality samples 29 and 31
# Resolution 0.1: 14 clusters
# Resolution 0.2: 18 clusters
# Resolution 0.4: 22 clusters
# Resolution 0.6: 26 clusters
# Resolution 0.8: 28 clusters
# Resolution 1.0: 30 clusters
# Resolution 1.2: 31 clusters
# Resolution 1.4: 32 clusters
# Resolution 2.0: 39 clusters


## clusters after running with dims=1:50. filtered nUMI<50000
# Resolution 0.1: 14 clusters
# Resolution 0.2: 17 clusters
# Resolution 0.4: 21 clusters
# Resolution 0.6: 26 clusters
# Resolution 0.8: 29 clusters
# Resolution 1.0: 30 clusters
# Resolution 1.2: 31 clusters
# Resolution 1.4: 33 clusters
# Resolution 2.0: 39 clusters

## clusters after running with dims=1:50. filtered nUMI<50000. filtered 30 doublets
# Resolution 0.1: 13 clusters
# Resolution 0.2: 17 clusters
# Resolution 0.4: 21 clusters
# Resolution 0.6: 25 clusters
# Resolution 0.8: 28 clusters
# Resolution 1.0: 29 clusters
# Resolution 1.2: 30 clusters
# Resolution 1.4: 30 clusters
# Resolution 2.0: 38 clusters



saveRDS(seurat_umap_sct, 'data/seurat_umap_sct_rm29and32_regressout_cc.diffandmt_doubletfinder.rds')
rm(seurat_merged_sct) # free RAM
gc()


### run UMAP on non-integrated, merged, sct-normalized data ----
seurat_umap_sct = RunUMAP(seurat_umap_sct,
                      dims=1:50,
                      reduction='pca') 

saveRDS(seurat_umap_sct, 'data/seurat_umap_sct_rm29and32_regressout_cc.diffandmt_doubletfinder.rds')
# seurat_umap_sct = readRDS('data/seurat_umap_sct.rds')

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'Phase',
        label=F,
        label.size=6,
        shuffle=TRUE )

DimPlot(seurat_umap_sct,
        reduction = "umap",group.by = 'condition',
        label=F,
        label.size=6,
        shuffle=TRUE)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.1'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.2'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.4'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.6'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.0.8'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1.2'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

Idents(seurat_umap_sct) = 'SCT_snn_res.2'
DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=T,
        label.size=6)

# start at 0.2
# 0.4 maybe
# 0.8 looks like best to start with
# 1.2 is the highest res to use probably
# 1.4 has same clusters as 1.2, but slightly different shape
# between cluster 19 and 12, which are the same cluster comparataively speaking.
# 1.2 has a better defined cluster shape
# 1. just has one less B cell cluster than 1.2 and 1.4

# 
# Can use 0.1, 0.2, 0.4, and 0.8


## visualize clustering results with clustree ----
library(clustree)

clustree(seurat_umap_sct, prefix='SCT_snn_res.')
ggsave('figures/clustering_rm29and32_regress_ccdiffandmt/sct/clustree.png', dpi=320, width=20, height= 10)

clustree(seurat_umap_sct, prefix = 'SCT_snn_res.',
         node_colour = 'Pgr', node_colour_aggr = 'mean')


## function to plot initial umaps ----
plot_umaps = function (seurat_umap, resolution, directory) {
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # Assign umap resolution
  Idents(seurat_umap) = resolution
  
  # Plot the UMAP of non-integrated data
  DimPlot(seurat_umap,
          reduction = "umap",
          label=TRUE,
          repel=TRUE,
          label.size=6)
  ggsave(paste0(directory,resolution,'/UMAP.png'), dpi=320, width=10, height=10)
  
  
  # plot UMAP by sample and condition
  p1 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='sample',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='condition',
               shuffle=T)
  
  png(filename=paste0(directory,resolution,'/UMAP_by_sample_and_condition.png'), units='in', width=20, height=10, res=320)
  print(plot_grid(p1, p2))
  dev.off()
  
  
  # plot UMAP by Phase and mitochondrial expression quartile
  p1 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='Phase',
               shuffle=T)
  p2 = DimPlot(seurat_umap,
               reduction = "umap",
               group.by='mito.level',
               shuffle=T)
  
  png(filename=paste0(directory,resolution,'/UMAP_by_cellcycle_and_mitoqrtl.png'), units='in', width=20, height=10, res=320)
  print(plot_grid(p1, p2))
  graphics.off()
}

## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.1 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')


#### UMAP clustering QC ----
metrics =  c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "percent.mt") # metrics to overlay on umap

## create function to plot all QC plots and umap overlays ----
umap_clustering_qc = function(seurat_umap, resolution, directory) {
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # set umap resolution
  Idents(seurat_umap) = resolution
  seurat_umap@meta.data$resolution = Idents(seurat_umap)
  
  # Explore whether clusters segregate by cell cycle phase
  DimPlot(seurat_umap,
          label = TRUE, 
          split.by = "Phase")  + NoLegend()
  ggsave(paste0(directory,resolution,'/cell_cycle_split_umap.png'), dpi=320, width=30, height=20)
  
  # Metrics to plot on UMAP
  FeaturePlot(seurat_umap, 
              reduction = "umap", 
              features = metrics,
              pt.size = 0.4, 
              order = TRUE,
              min.cutoff = 'q10',
              label = TRUE)
  ggsave(paste0(directory,resolution,'/counts_cycle_mt_metrics.png'), dpi=320, width=20, height=20)
  
  # Boxplot of nUMIs per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nCount_RNA, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nUMIs_per_cluster.png'), dpi=320, width=20, height= 10)
  
  # Boxplot of nGenes per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=nFeature_RNA, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nGenes_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Boxplot of percent.mt per cluster. 
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=percent.mt, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/percent.mt_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Boxplot of S.Score per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=S.Score, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/S.score_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Boxplot of G2M.Score per cluster
  ggplot(seurat_umap@meta.data) +
    geom_boxplot(aes(x=resolution, y=G2M.Score, fill=resolution))+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    NoLegend()+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/G2M.score_per_cluster.png'), dpi=320, width=20, height= 10)
  
  
  # Extract identity and sample information from seurat object to determine 
  # the number of cells per cluster per sample
  n_cells = FetchData(seurat_umap_sct, 
                       vars = c("ident", "sample")) %>%
    dplyr::count(ident, sample)
  
  # Barplot of number of cells per cluster by sample
  ggplot(n_cells, aes(x=ident, y=n, fill=sample))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
    
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_sample.png'), dpi=320, width=20, height= 10)
  
  
  # Extract identity and condition information from seurat object to determine 
  # the number of cells per cluster per condition
  n_cells = FetchData(seurat_umap, 
                      vars = c("ident", "condition")) %>%
    dplyr::count(ident, condition)
  
  # Barplot of number of cells per cluster by condition
  ggplot(n_cells, aes(x=ident, y=n, fill=condition))+
    geom_bar(position=position_dodge(), stat="identity")+
    geom_text(aes(label=n), vjust = -.2, position=position_dodge(1))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5))
  ggsave(paste0(directory,resolution,'/nCells_per_cluster_by_condition.png'), dpi=320, width=20, height= 10)
  
  # Barplot of proportion of cells in each cluster by sample
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=sample), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5),
          legend.key.size = unit(1,'cm'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20))
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_sample.png'), dpi=320, width=20, height= 10)
  
  
  # Barplot of proportion of cells in each cluster by condition
  ggplot(seurat_umap@meta.data) +
    geom_bar(aes(x=resolution, fill=condition), position=position_fill())+
    scale_x_discrete(name = 'Ident') +
    scale_fill_discrete(name = paste0("Resolution: ", resolution))+
    theme_classic2()+
    theme(axis.text.x = element_text(angle=-15,vjust=.5),
          legend.key.size = unit(1,'cm'),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20))
  ggsave(paste0(directory,resolution,'/proportion_of_cells_per_cluster_by_condition.png'), dpi=320, width=20, height= 10)
}


# clustering QC at SCT_snn_res.0.1 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')


## Which PCs drive different clusters? ----

# determine max number of PCs present in seurat object
max_pcs = length(grep("^PC_", colnames(seurat_umap_sct@reductions$pca@cell.embeddings)))

# Create PC names based on the max PCs present
columns = c(paste0("PC_", 1:max_pcs),
            "ident",
            "umap_1", "umap_2")

## overlay PCs on umap ----
# takes a while to run so will use parallel processing
# library(future.apply)
# library(parallel)
# 
# options(future.globals.maxSize = 8 * 1024^3)  # Set limit to 8 GB

## function to overlay PCs on umap
# takes a while to run so will use parallel processing with 3 workers
pcs_umap_overlay = function(seurat_umap, resolution, directory) {
  
  # plan(list(sequential, tweak(multisession, workers=2)))
  
  DefaultAssay(seurat_umap) = 'SCT'
  
  # set resolution
  Idents(seurat_umap) = resolution
  
  # Extract columns from the seurat object. FetchData() can pull any data
  # from expression matrices, cell embeddings, or metadata
  pc_data = FetchData(seurat_umap, vars = columns)
  
  # Adding cluster label to center of cluster on UMAP
  umap_label = FetchData(seurat_umap, 
                          vars = c("ident", "umap_1", "umap_2"))  %>%
    group_by(ident) %>%
    dplyr::summarise(x=mean(umap_1), y=mean(umap_2))
  
  # Chunk PCs into blocks of 12
  pc_chunks = split(paste0("PC_", 1:max_pcs), ceiling(seq_along(1:max_pcs)/12))
  
  # Plotting a UMAP plot for each of the PCs
  lapply(seq_along(pc_chunks), function(chunk){
    
    pc_range = pc_chunks[[chunk]]
    first_pc = pc_range[1]
    last_pc = pc_range[length(pc_range)]
    
    lapply(pc_chunks[[chunk]], function(pc) {
      ggplot(pc_data, aes(umap_1, umap_2))+
        geom_point(aes(color=!!sym(pc)), alpha = 0.7)+
        scale_color_gradient(guide = 'none', low = "grey90", high = "blue")+
        geom_text(data=umap_label, aes(x=x, y=y, label=ident))+
        ggtitle(pc)+
        theme(plot.title = element_text(hjust = 0.5, size=24, face='bold'))+
        theme_classic2()
    }) %>% 
    plot_grid(plotlist = .)
    ggsave(paste0(directory,resolution,'/UMAP_by_PC_components_',first_pc,'-',last_pc,'.png'), dpi=320, width=20, height=15)
    
  })
  # on.exit(plan(sequential)) # restore to default non-parallel settings
}

## plot PCs on umap for resolution SCT_snn_res.0.1 ----
# takes a while to run
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# PC_7 looks like it's driving cell cycle genes
# PC_13 looks like it's driving myoepithelial cells
# PC_14 marcophage mb
# PC_17 mystery cluster by endothelial cells
# PC_43 and PC_46 maybe cells transitioning between fibroblast and epithelial 


#### exploring cell marker genes ----
library(readxl)

# read in cell markers from excel file
cell_markers = read_excel('annotation/aggregate_sources_celltype_markers.xlsx')
cell_markers = lapply(cell_markers, na.omit) # convert to list of lists, omit NAs

# set RNA counts slot to be default and normalize for visualizing genes
# SCTransform normalization was performed only on the 3000 most variable genes, 
# so many of our genes of interest may not be present in the SCT data. Hence,
# we switch over to the original RNA count to explore cell type marker genes.
DefaultAssay(seurat_umap_sct) = 'RNA'
# seurat_umap_sct = NormalizeData(seurat_umap_sct, verbose = FALSE)
# DefaultAssay(seurat_umap_sct) = 'SCT'
# saveRDS(seurat_umap_sct, 'data/seurat_umap_sct.rds')

# Want to only consider cell markers that are actually present in the dataset
cell_markers = lapply(cell_markers, function(x){
  x = unique(x[x %in% rownames(seurat_umap_sct)])
})

# remove any cell types with zero markers in the dataset
cell_markers = Filter(length, cell_markers)
# Erythrocyte cell type removed since no gene markers present in this dataset

# process and save the gene markers actually present in the data into a new file
max_length = max(sapply(cell_markers, length)) # max length lists of gene markers

# Pad shorter lists with NA to match max length so all gene lists have same length
padded_markers = lapply(cell_markers, function(x) {
  length(x) = max_length
  return(x)
})

# save gene markers actually present in the dataset to csv
write.csv(as.data.frame(padded_markers), 
          file='annotation/aggregate_sources_celltype_markers_actually_in_data.csv',
          row.names=FALSE,
          na='')


# function to plot all gene markers on umap and dotplots for each cell type in batches of 16 genes ----
# takes a while to run
plot_gene_markers = function(seurat_umap, resolution, directory) {
  
  DefaultAssay(seurat_umap) = 'RNA'
  
  # set resolution
  Idents(seurat_umap) = resolution
  
  lapply(seq_along(cell_markers), function(i) {
    gene_chunks = split(cell_markers[[i]], ceiling(seq_along(cell_markers[[i]])/16)) # plot 16 gene markers at a time
    lapply(seq_along(gene_chunks), function(j) {
      FeaturePlot(seurat_umap, # overlay gene marker expression on umap
                  reduction='umap',
                  features=gene_chunks[[j]],
                  order=T,
                  min.cutoff='q10',
                  label=T)+
        plot_annotation(
          title=paste0(names(cell_markers[i]),'-',j),
          theme=theme(plot.title = element_text(hjust = 0.5,size=24,face='bold')))
      ggsave(paste0(directory,resolution,'/annotation/',names(cell_markers[i]),'-',j,'.png'), 
                    dpi=320, width=20, height=20)
      
      png(filename=paste0(directory,resolution,'/annotation/',names(cell_markers[i]),'-',j,'_dotplot.png'), 
          units='in', width=20, height=10, res=320)
      print(
        DotPlot(seurat_umap, features=gene_chunks[[j]])+ # construct dotplot of gene markers for each cell type
          RotatedAxis()+
          ggtitle(paste0(names(cell_markers[i]),'-',j))+
          theme(plot.title = element_text(hjust = 0.5, size=24, face='bold'))
      )
      graphics.off()
    })
  })
}
  
  
## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.1. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.1', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# looks like we could benefit from increasing resolution
# resolutions available. For reference
resolutions = sprintf('SCT_snn_res.%.1f',c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 2.0))

# Assign identity of clusters to desired resolution
Idents(seurat_umap_sct) = 'SCT_snn_res.1.4'

DimPlot(seurat_umap_sct,
        reduction = "umap",
        label=TRUE,
        label.size=6)


#### Resolution SCT_snn_res.0.2 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.2 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# clustering QC at SCT_snn_res.0.2 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## plot PCs on umap for resolution SCT_snn_res.0.2 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.2. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')


#### Resolution SCT_snn_res.0.4 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.4 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# clustering QC at SCT_snn_res.0.4 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## plot PCs on umap for resolution SCT_snn_res.0.4 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.4. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')


#### Resolution SCT_snn_res.0.8 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.0.8 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# clustering QC at SCT_snn_res.0.8 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## plot PCs on umap for resolution SCT_snn_res.0.8 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')


## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.0.8. ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.0.8', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')



#### Resolution SCT_snn_res.1.2 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.1.2 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# clustering QC at SCT_snn_res.1.2 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## plot PCs on umap for resolution SCT_snn_res.1.2 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.1.2 ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.2', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')


#### Resolution SCT_snn_res.1.4 ----
## plot regular umap, as well as by sample, condition, cell cycle, and mitochondrial expression quartile for resolution SCT_snn_res.1.4 ----
plot_umaps(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

# clustering QC at SCT_snn_res.1.4 resolution ----
umap_clustering_qc(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## plot PCs on umap for resolution SCT_snn_res.1.4 ----
pcs_umap_overlay(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')

## overlay umap and create dotplot with gene markers for resolution SCT_snn_res.1.4 ----
# takes a while to run
plot_gene_markers(seurat_umap=seurat_umap_sct, resolution='SCT_snn_res.1.4', directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/')



#### annotate clusters ----
# will annotate at a resolution of 0.8 and subcluster as needed
Idents(seurat_umap_sct) = 'SCT_snn_res.0.8' # set cluster #s as ident

DefaultAssay(seurat_umap_sct) = "RNA"

#### Find marker genes ----

## Find all positive markers
# join layers to be able to run FindAllMarkers()
combined_seurat = JoinLayers(seurat_umap_sct)

all.pos.markers = FindAllMarkers(combined_seurat,
                                 only.pos=TRUE,
                                 logfc.threshold = 0.25,
                                 min.pct=0.1,
                                 min.diff.pct = -Inf) 

# Combine markers with gene descriptions 
annotations = read.csv('annotation/annotationhub_mice_genes.csv')
ann.pos.markers = left_join(x = all.pos.markers, 
                          y = annotations[, c("gene_name", "description")],
                          by = c("gene" = "gene_name")) %>%
  unique()

# add pct.diff col
ann.pos.markers$pct.diff = ann.pos.markers$pct.1 - ann.pos.markers$pct.2

# Rearrange the columns to be more intuitive
ann.pos.markers = ann.pos.markers[ , c(6, 7, 2:4, 9, 1, 5,8)]

# Order the rows by p-adjusted values
ann.pos.markers = ann.pos.markers %>%
  dplyr::arrange(cluster, p_val_adj)

# save markers
write.csv(ann.pos.markers,'annotation/annotated_pos_markers_res.0.8_df.csv',row.names=FALSE)

# Extract top 10 markers per cluster
top10 = ann.pos.markers %>% 
  group_by(cluster) %>%
  slice_max(order_by=tibble(avg_log2FC,pct.diff,p_val_adj),n=10)


## Find all positive conserved markers between conditions ----
# Create function to get conserved markers for any given cluster
get_conserved = function(cluster){
  FindConservedMarkers(combined_seurat, # need to specify exact seurat object here
                       ident.1 = cluster,
                       grouping.var = "condition", # condition variable
                       only.pos=TRUE,
                       logfc.threshold = 0.25,
                       min.pct=0.1,
                       min.diff.pct = -Inf) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}


# get conserved markers
conserved.pos.markers = map_dfr(c(0:27), get_conserved) # 28 clusters at resolution 0.8

# add avg.pct.diff, max_adj_pval, and min.pct 1 cols
conserved.pos.markers = conserved.pos.markers %>%
  mutate(avg_log2FC = (ctrl_avg_log2FC + pb_avg_log2FC) / 2,
         max_adj_pval = ifelse(ctrl_p_val_adj > pb_p_val_adj,ctrl_p_val_adj,pb_p_val_adj),
         avg.pct.diff = ((pb_pct.1 - pb_pct.2) + (ctrl_pct.1 - ctrl_pct.2)) / 2,
         min.pct.1 = ifelse(ctrl_pct.1 < pb_pct.1, ctrl_pct.1,pb_pct.1))

# Rearrange the columns to be more intuitive
conserved.pos.markers = conserved.pos.markers[ , c(1:2,16:19,15,3:14)]

# Order the rows by p-adjusted values
conserved.pos.markers = conserved.pos.markers %>%
  dplyr::arrange(cluster_id, max_adj_pval)

# save markers
write.csv(conserved.pos.markers,'annotation/conserved_pos_markers_res.0.8_df.csv',row.names=FALSE)


# Extract top 10 conserved markers per cluster
top10.conserved = conserved.pos.markers %>% 
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)

# top10.conserved = top10.conserved[!top10.conserved$cluster_id == 15,]
top10.conserved = top10.conserved[!top10.conserved$cluster_id == 26,]
# top10.conserved = top10.conserved[!top10.conserved$cluster_id == 28,]
# top10.conserved = top10.conserved[!(top10.conserved$cluster_id == 29 |
#                                            top10.conserved$cluster_id == 30),]


# # clusters 15 is made up of one condition. 
# # Will treat individual values as avg and max to easily bind rows to top10.conserved
# topc15 = conserved.pos.markers[conserved.pos.markers$cluster_id == 15,] %>%
#   mutate(avg_fc = ctrl_avg_log2FC,
#          max_adj_pval = ctrl_p_val_adj,
#          avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
#   group_by(cluster_id) %>%
#   slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)
# 
# top10.conserved = rbind(top10.conserved, topc15)


# clusters 26 is made up of one condition.
# Will treat individual values as avg and max to easily bind rows to top10.conserved
topc26 = conserved.pos.markers[conserved.pos.markers$cluster_id == 26,] %>%
  mutate(avg_fc = ctrl_avg_log2FC,
         max_adj_pval = ctrl_p_val_adj,
         avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
  group_by(cluster_id) %>%
  slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)

top10.conserved = rbind(top10.conserved, topc26)

# # clusters 28 is made up of one condition.
# # Will treat individual values as avg and max to easily bind rows to top10.conserved
# topc28 = conserved.pos.markers[conserved.pos.markers$cluster_id == 15,] %>%
#   mutate(avg_fc = ctrl_avg_log2FC,
#          max_adj_pval = ctrl_p_val_adj,
#          avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
#   group_by(cluster_id) %>%
#   slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)
# 
# top10.conserved = rbind(top10.conserved, topc28)


# # clusters 29 and 30 are made up of one condition. 
# # Will treat individual values as avg and max to easily bind rows to top10.conserved
# topc29 = conserved.pos.markers[conserved.pos.markers$cluster_id == 29,] %>%
#   mutate(avg_fc = pb_avg_log2FC,
#          max_adj_pval = pb_p_val_adj,
#          avg.pct.diff = pb_pct.1 - pb_pct.2) %>%
#   group_by(cluster_id) %>%
#   slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)
# 
# topc30 = conserved.pos.markers[conserved.pos.markers$cluster_id == 30,] %>%
#   mutate(avg_fc = ctrl_avg_log2FC,
#          max_adj_pval = ctrl_p_val_adj,
#          avg.pct.diff = ctrl_pct.1 - ctrl_pct.2) %>%
#   group_by(cluster_id) %>%
#   slice_max(order_by=tibble(avg_fc,avg.pct.diff,max_adj_pval),n=10)
# 
# top10.conserved = rbind(top10.conserved, topc29, topc30)

## test subclustering t cells ----
test = subset(DietSeurat(seurat_umap_sct, assays='RNA'), subset=(seurat_clusters == 1) | 
                (seurat_clusters == 2) | (seurat_clusters == 4) | 
                (seurat_clusters == 10) | (seurat_clusters == 13) | 
                (seurat_clusters == 20))

test = SplitObject(test, split.by='sample')

test = lapply(test, function(x) {
  x = SCTransform(x, assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                  vars.to.regress=c('cc.difference','percent.mt'),
                  return.only.var.genes=FALSE)
})

features = SelectIntegrationFeatures(test, nfeatures=5000, fvf.nfeatures=5000)

# merge sctransform normalized samples. Takes a little while to run
test = merge(test[[1]], test[-1], merge.data=TRUE,merge.dr=FALSE)

VariableFeatures(test) = features

sct_results = SCTResults(test, slot = "feature.attributes")
saveRDS(sct_results, 'figures/clustering_rm29and32_regress_ccdiffandmt/sct_results_test.rds')
rm(sct_results)

# Variable features lost/unset after merging transformed objects. Need to set again.
# should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(test[["SCT"]]) = rownames(test[["SCT"]]@scale.data)

test = RunPCA(test, assay='SCT', npcs=50)
test = FindNeighbors(test, dims=1:50)
test = FindClusters(object=test, resolution=c(0.1,0.4,0.8,1.2), graph.name='SCT_snn')
test = RunUMAP(test, dims=1:50, reduction='pca')

Idents(test) = 'SCT_snn_res.0.4'
DimPlot(test,
        reduction='umap',
        label=TRUE,
        shuffle=TRUE)

table(test$condition, test$SCT_snn_res.0.4)




## test subclustering epithelial cells ----
# normalization
test = subset(DietSeurat(seurat_umap_sct, assays='RNA'), subset=(seurat_clusters == 18) | 
                (seurat_clusters == 23))

test = SplitObject(test, split.by='sample')

test = lapply(test, function(x) {
  x = NormalizeData(x)
  x = FindVariableFeatures(x)
  x = ScaleData(x, vars.to.regress=c('cc.difference','percent.mt'))
})

test = merge(test[[1]], test[-1], merge.data=TRUE,merge.dr=FALSE)
test = JoinLayers(test)


test = RunPCA(test, npcs=50)

ElbowPlot(object = test, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))

### Determine percent of variation associated with each PC ----
pct = test[["pca"]]@stdev / sum(test[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 42
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 12
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 


test = FindNeighbors(test, dims=1:20)
test = FindClusters(object=test, resolution=c(0.1,0.4,0.8,1.0,1.2))
test = RunUMAP(test, dims=1:20, reduction='pca')

Idents(test) = 'RNA_snn_res.0.4'
DimPlot(test,
        reduction='umap',
        label=TRUE,
        shuffle=TRUE)

table(test$condition, test$RNA_snn_res.0.4)
# 0   1   2   3   4
# ctrl 182 174 117  37   2
# pb   114  77  71  37   1






# SCTransform
test = subset(DietSeurat(seurat_umap_sct, assays='RNA'), subset=(seurat_clusters == 18) | 
                (seurat_clusters == 23))

test = SplitObject(test, split.by='sample')

test = lapply(test, function(x) {
  x = SCTransform(x, assay='RNA', vst.flavor ="v2", variable.features.n=3000,# default params
                  vars.to.regress=c('cc.difference','percent.mt'),
                  return.only.var.genes=FALSE)
})

features = SelectIntegrationFeatures(test, nfeatures=3000, fvf.nfeatures=3000)

# merge sctransform normalized samples. Takes a little while to run
test = merge(test[[1]], test[-1], merge.data=TRUE,merge.dr=FALSE)

VariableFeatures(test) = features

sct_results = SCTResults(test, slot = "feature.attributes")
saveRDS(sct_results, 'figures/clustering_rm29and32_regress_ccdiffandmt/sct_results_test_epithelial.rds')
rm(sct_results)

# # Variable features lost/unset after merging transformed objects. Need to set again.
# # should do this according to https://github.com/satijalab/seurat/issues/2814
VariableFeatures(test[["SCT"]]) = rownames(test[["SCT"]]@scale.data)

test = RunPCA(test, assay='SCT', npcs=50)

ElbowPlot(object = test, 
          ndims = 50)+
  ggtitle('Elbow plot of PCs')+
  theme_classic2()+
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18))

### Determine percent of variation associated with each PC ----
pct = test[["pca"]]@stdev / sum(test[["pca"]]@stdev) * 100
pct
#
# # Calculate cumulative percents for each PC
cumu = cumsum(pct)
cumu
#
# # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 = which(cumu >= 90 & pct <= 5)[1]
co1 # 42
#
# # Determine the difference between variation of PC and subsequent PC.
# # This is the last point where change of % of variation is more than 0.1%. Afterwards, % of variation change is < 0.1%
co2 = sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # 12
#
# # Minimum of the two calculations
pcs = min(co1, co2)
pcs 



test = FindNeighbors(test, dims=1:30)
test = FindClusters(object=test, resolution=c(0.1,0.4,0.8,1.0,1.2), graph.name='SCT_snn')
test = RunUMAP(test, dims=1:30, reduction='pca')

Idents(test) = 'SCT_snn_res.0.4'
DimPlot(test,
        reduction='umap',
        label=TRUE,
        shuffle=TRUE)

table(test$condition, test$SCT_snn_res.1.2)
# 0   1   2   3   4
# ctrl 182 174 117  37   2
# pb   114  77  71  37   1

DoHeatmap(test, features[1:100],size=4,angle=90)+NoLegend()


## test epithelial cells ----
# cluster 0
test_epithelial = FindMarkers(test,
                         ident.1 = 0,
                         ident.2 = c(1:4)) 

# Add gene symbols to the DE table
test_epithelial = test_epithelial %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
test_epithelial$pct.diff = test_epithelial$pct.1 - test_epithelial$pct.2

# Reorder columns and sort by padj   
test_epithelial = test_epithelial[, c(1, 3:5,2,8,6:7)]

test_epithelial = test_epithelial %>%
  dplyr::arrange(p_val_adj)

write.csv(test_epithelial,'annotation/test_epithelial_cluster0.csv',row.names=FALSE)

FeaturePlot(test, features = c('Esr1','Prlr','Fgfr2','Areg','Egfr','Igf1','Gh'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()






### FindMarkers() for clusters still in question ----
## Cd4 T cells
# cluster 4 vs 1 and 10
cd4_tcells = FindMarkers(seurat_umap_sct,
                         ident.1 = 4,
                         ident.2 = c(1,10)) 

# Add gene symbols to the DE table
cd4_tcells = cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
cd4_tcells$pct.diff = cd4_tcells$pct.1 - cd4_tcells$pct.2

# Reorder columns and sort by padj   
cd4_tcells = cd4_tcells[, c(1, 3:5,2,8,6:7)]

cd4_tcells = cd4_tcells %>%
  dplyr::arrange(p_val_adj)

write.csv(cd4_tcells,'annotation/cd4_tcells_4vs1and10.csv',row.names=FALSE)

# cluster 10 vs 1 and 4
cd4_tcells = FindMarkers(seurat_umap_sct,
                         ident.1 = 10,
                         ident.2 = c(1,4)) 

# Add gene symbols to the DE table
cd4_tcells = cd4_tcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
cd4_tcells$pct.diff = cd4_tcells$pct.1 - cd4_tcells$pct.2

# Reorder columns and sort by padj   
cd4_tcells = cd4_tcells[, c(1, 3:5,2,8,6:7)]

cd4_tcells = cd4_tcells %>%
  dplyr::arrange(p_val_adj)

write.csv(cd4_tcells,'annotation/cd4_tcells_10vs1and4.csv',row.names=FALSE)

#B cells
# cluster 3 vs 0,11,and 16
bcells = FindMarkers(seurat_umap_sct,
                         ident.1 = 3,
                         ident.2 = c(0,11,16)) 

# Add gene symbols to the DE table
bcells = bcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
bcells$pct.diff = bcells$pct.1 - bcells$pct.2

# Reorder columns and sort by padj   
bcells = bcells[, c(1,3:5,2,8,6:7)]

bcells = bcells %>%
  dplyr::arrange(p_val_adj)

write.csv(bcells,'annotation/bcells_cluster3.csv',row.names=FALSE)


# cluster 11 vs 0,3,and 16
bcells = FindMarkers(seurat_umap_sct,
                     ident.1 = 11,
                     ident.2 = c(0,3,16)) 

# Add gene symbols to the DE table
bcells = bcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
bcells$pct.diff = bcells$pct.1 - bcells$pct.2

# Reorder columns and sort by padj   
bcells = bcells[, c(1,3:5,2,8,6:7)]

bcells = bcells %>%
  dplyr::arrange(p_val_adj)

write.csv(bcells,'annotation/bcells_cluster11.csv',row.names=FALSE)


# cluster 16 vs 0,3,and 11
bcells = FindMarkers(seurat_umap_sct,
                     ident.1 = 16,
                     ident.2 = c(0,3,11)) 

# Add gene symbols to the DE table
bcells = bcells %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
bcells$pct.diff = bcells$pct.1 - bcells$pct.2

# Reorder columns and sort by padj   
bcells = bcells[, c(1,3:5,2,8,6:7)]

bcells = bcells %>%
  dplyr::arrange(p_val_adj)

write.csv(bcells,'annotation/bcells_cluster16.csv',row.names=FALSE)



## adipocytes
# cluster 21 vs 7
adipocytes = FindMarkers(seurat_umap_sct,
                     ident.1 = 21,
                     ident.2 = c(7))

# Add gene symbols to the DE table
adipocytes = adipocytes %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
adipocytes$pct.diff = adipocytes$pct.1 - adipocytes$pct.2

# Reorder columns and sort by padj   
adipocytes = adipocytes[, c(1,3:5,2,8,6:7)]

adipocytes = adipocytes %>%
  dplyr::arrange(p_val_adj)

write.csv(adipocytes,'annotation/adipocytes_cluster21.csv',row.names=FALSE)


# cluster 7 vs 21
adipocytes = FindMarkers(seurat_umap_sct,
                         ident.1 = 7,
                         ident.2 = 21)

# Add gene symbols to the DE table
adipocytes = adipocytes %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
adipocytes$pct.diff = adipocytes$pct.1 - adipocytes$pct.2

# Reorder columns and sort by padj   
adipocytes = adipocytes[, c(1,3:5,2,8,6:7)]

adipocytes = adipocytes %>%
  dplyr::arrange(p_val_adj)

write.csv(adipocytes,'annotation/adipocytes_cluster7.csv',row.names=FALSE)


## fibroblasts
# clusters 17 vs 6 and 24
fibroblasts = FindMarkers(seurat_umap_sct,
                         ident.1 = 17,
                         ident.2 = c(6,24))

# Add gene symbols to the DE table
fibroblasts = fibroblasts %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
fibroblasts$pct.diff = fibroblasts$pct.1 - fibroblasts$pct.2

# Reorder columns and sort by padj   
fibroblasts = fibroblasts[, c(1,3:5,2,8,6:7)]

fibroblasts = fibroblasts %>%
  dplyr::arrange(p_val_adj)

write.csv(fibroblasts,'annotation/fibroblasts_cluster17.csv',row.names=FALSE)

# cluster 24 vs 6 and 17
fibroblasts = FindMarkers(seurat_umap_sct,
                          ident.1 = 24,
                          ident.2 = c(6,17))

# Add gene symbols to the DE table
fibroblasts = fibroblasts %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
fibroblasts$pct.diff = fibroblasts$pct.1 - fibroblasts$pct.2

# Reorder columns and sort by padj   
fibroblasts = fibroblasts[, c(1,3:5,2,8,6:7)]

fibroblasts = fibroblasts %>%
  dplyr::arrange(p_val_adj)

write.csv(fibroblasts,'annotation/schwann_cells_cluster24.csv',row.names=FALSE)


## epithelial
# clusters 23 vs 18 and 24
epithelial = FindMarkers(seurat_umap_sct,
                          ident.1 = 23,
                          ident.2 = c(18,24))

# Add gene symbols to the DE table
epithelial = epithelial %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
epithelial$pct.diff = epithelial$pct.1 - epithelial$pct.2

# Reorder columns and sort by padj   
epithelial = epithelial[, c(1,3:5,2,8,6:7)]

epithelial = epithelial %>%
  dplyr::arrange(p_val_adj)

write.csv(epithelial,'annotation/epithelial_cluster23.csv',row.names=FALSE)


# clusters 18 vs 23 and 24
epithelial = FindMarkers(seurat_umap_sct,
                         ident.1 = 18,
                         ident.2 = c(23,24))

# Add gene symbols to the DE table
epithelial = epithelial %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# add pct.diff col
epithelial$pct.diff = epithelial$pct.1 - epithelial$pct.2

# Reorder columns and sort by padj   
epithelial = epithelial[, c(1,3:5,2,8,6:7)]

epithelial = epithelial %>%
  dplyr::arrange(p_val_adj)

write.csv(epithelial,'annotation/epithelial_cluster18.csv',row.names=FALSE)



## endothelial

# myeloid

rownames(seurat_umap_sct)


# Rename cluster #'s to cell type annotations ----
Idents(seurat_umap_sct) = 'SCT_snn_res.0.8'
seurat_umap_sct = RenameIdents(seurat_umap_sct,
                                 '0' = 'B cells',
                                 '1' = 'CD4+ T cells',
                                 '2' = 'CD8+ T cells',
                                 '3' = 'B cells',
                                 '4' = 'Tregs',
                                 '5' = 'Endothelial',
                                 '6' = 'Fibroblast',
                                 '7' = 'Adipocytes',
                                 '8' = 'Proliferating cells',
                                 '9' = 'Dendritic cell',
                                 '10' = 'Activated CD4+ T cells',
                                 '11' = 'B cells',
                                 '12' = 'Macrophage.Ma',
                                 '13' = 'T helper 17 cells',
                                 '14' = 'Myeloid.Mb',
                                 '15' = 'Myeloid.Mb',
                                 '16' = 'B cells',
                                 '17' = 'Reticular cells',
                                 '18' = 'Basal-Myoepithelial',
                                 '19' = 'Endothelial',
                                 '20' = 'CD8+ T cells',
                                 '21' = 'Adipocytes',
                                 '22' = 'Myeloid.Mb',
                                 '23' = 'Luminal.HS',
                                 '24' = 'Schwann cells',
                                 '25' = 'Endothelial',
                                 '26' = 'Muscle',
                                 '27' = 'Erythroid cells')

# assign cluster names to meta data
seurat_umap_sct$celltype = Idents(seurat_umap_sct)
Idents(seurat_umap_sct) = 'celltype'

DimPlot(seurat_umap_sct, 
        group.by = 'celltype',
        label=T,
        label.size=6,
        repel=T,
        shuffle=T)
ggsave('annotation/figures/annotated_seurat_res.0.8_df.png', dpi=320, width=20, height=20)

saveRDS(seurat_umap_sct, 'data/annotated_seurat_res.0.8_df.rds')

# plot qc metrics for cell-annotated clusters
plot_umaps(seurat_umap_sct, 
           resolution='celltype',
           directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/SCT_snn_res.0.8/')

umap_clustering_qc(seurat_umap_sct, 
                   resolution='celltype',
                   directory='figures/clustering_rm29and32_regress_ccdiffandmt/sct/SCT_snn_res.0.8/')

# shrink seurat object in preparation of DE analysis
de_seurat = DietSeurat(seurat_umap_sct, assays='RNA', layers=c('counts'))
saveRDS(de_seurat, 'data/de_seurat_annotated_19clusters.rds')




#### sessionInfo ----
sessionInfo()
# 
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
#   [1] DoubletFinder_2.0.4 presto_1.0.0        data.table_1.16.0   Rcpp_1.0.13         pheatmap_1.0.12     Polychrome_1.5.1    RCurl_1.98-1.16     metap_1.11          glmGamPoi_1.14.3   
# [10] scales_1.3.0        viridis_0.6.5       viridisLite_0.4.2   plyr_1.8.9          ggpubr_0.6.0        ggthemes_5.1.0      Matrix_1.6-5        lubridate_1.9.3     forcats_1.0.0      
# [19] stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2         readr_2.1.5         tidyr_1.3.1         tibble_3.2.1        ggplot2_3.5.1       tidyverse_2.0.0     cowplot_1.1.3      
# [28] patchwork_1.3.0     sctransform_0.4.1   Seurat_5.1.0        SeuratObject_5.0.2  sp_2.1-4           
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22            splines_4.3.1               later_1.3.2                 bitops_1.0-8                polyclip_1.10-7             fastDummies_1.7.4          
# [7] lifecycle_1.0.4             Rdpack_2.6.1                rstatix_0.7.2               globals_0.16.3              lattice_0.22-6              MASS_7.3-60.0.1            
# [13] backports_1.5.0             magrittr_2.0.3              plotly_4.10.4               plotrix_3.8-4               qqconf_1.3.2                httpuv_1.6.15              
# [19] sn_2.1.1                    spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0           pbapply_1.7-2               RColorBrewer_1.1-3         
# [25] multcomp_1.4-26             abind_1.4-8                 zlibbioc_1.48.2             Rtsne_0.17                  GenomicRanges_1.54.1        BiocGenerics_0.48.1        
# [31] TH.data_1.1-2               sandwich_3.1-1              GenomeInfoDbData_1.2.11     IRanges_2.36.0              S4Vectors_0.40.2            ggrepel_0.9.6              
# [37] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-0        TFisher_0.2.0               goftest_1.2-3               RSpectra_0.16-2            
# [43] spatstat.random_3.3-2       fitdistrplus_1.2-1          parallelly_1.38.0           leiden_0.4.3.1              codetools_0.2-20            DelayedArray_0.28.0        
# [49] tidyselect_1.2.1            farver_2.1.2                matrixStats_1.4.1           stats4_4.3.1                spatstat.explore_3.3-2      mathjaxr_1.6-0             
# [55] jsonlite_1.8.9              multtest_2.58.0             progressr_0.14.0            Formula_1.2-5               ggridges_0.5.6              survival_3.7-0             
# [61] tools_4.3.1                 ica_1.0-3                   glue_1.7.0                  mnormt_2.1.1                gridExtra_2.3               SparseArray_1.2.4          
# [67] MatrixGenerics_1.14.0       GenomeInfoDb_1.38.8         numDeriv_2016.8-1.1         withr_3.0.1                 fastmap_1.2.0               fansi_1.0.6                
# [73] digest_0.6.37               timechange_0.3.0            R6_2.5.1                    mime_0.12                   colorspace_2.1-1            scattermore_1.2            
# [79] tensor_1.5                  spatstat.data_3.1-2         utf8_1.2.4                  generics_0.1.3              httr_1.4.7                  htmlwidgets_1.6.4          
# [85] S4Arrays_1.2.1              scatterplot3d_0.3-44        uwot_0.2.2                  pkgconfig_2.0.3             gtable_0.3.5                lmtest_0.9-40              
# [91] XVector_0.42.0              htmltools_0.5.8.1           carData_3.0-5               dotCall64_1.1-1             Biobase_2.62.0              png_0.1-8                  
# [97] spatstat.univar_3.0-1       rstudioapi_0.16.0           tzdb_0.4.0                  reshape2_1.4.4              nlme_3.1-166                zoo_1.8-12                 
# [103] KernSmooth_2.23-24          parallel_4.3.1              miniUI_0.1.1.1              pillar_1.9.0                grid_4.3.1                  vctrs_0.6.5                
# [109] RANN_2.6.2                  promises_1.3.0              car_3.1-3                   xtable_1.8-4                cluster_2.1.6               mvtnorm_1.3-1              
# [115] cli_3.6.3                   compiler_4.3.1              rlang_1.1.4                 crayon_1.5.3                mutoss_0.1-13               future.apply_1.11.2        
# [121] ggsignif_0.6.4              stringi_1.8.4               deldir_2.0-4                munsell_0.5.1               lazyeval_0.2.2              spatstat.geom_3.3-3        
# [127] RcppHNSW_0.6.0              hms_1.1.3                   future_1.34.0               shiny_1.9.1                 SummarizedExperiment_1.32.0 rbibutils_2.2.16  