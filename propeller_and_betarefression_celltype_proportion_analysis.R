library(Seurat)
library(speckle)
library(tidyverse)
library(qs2)


nthreads = 2

analysis <- 'cellbender_analysis'
## directories for integrated data ----
fig_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/clustering/")
data_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/data/")
annot_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/annotation/")
stats_dir_integrated = paste0("c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/DA_analysis/")


# Create directories to save results if they don't already exist:
dirs = c(fig_dir_integrated, data_dir_integrated, annot_dir_integrated)

for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir,
                                     recursive=TRUE) }
}

# integrated_seurat = qs_read(paste0(data_dir_integrated,'integrated_seurat_sctype_annotated.qs'))


# # directories for saving
# output_annot_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/annotation/"
# 
# output_fig_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/clustering/"
# output_data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/data/"
# 
# output_stats_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/stats/"
# 
# if (!dir.exists(paste0(output_data_dir))) { dir.create(paste0(output_data_dir),
#                                                        recursive=TRUE) }
# if (!dir.exists(paste0(output_fig_dir))) { dir.create(paste0(output_fig_dir),
#                                                       recursive=TRUE) }
# if (!dir.exists(paste0(output_stats_dir))) { dir.create(paste0(output_stats_dir),
#                                                         recursive=TRUE) }
# if (!dir.exists(paste0(output_annot_dir))) { dir.create(paste0(output_annot_dir),
#                                                        recursive=TRUE) }
# 
# # data directories
# data_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary/3_weeks/cellbender_analysis/FPR_0.0/data/"


# Extract unique timepoints and conditions
timepoints <- unique(annotated_seurat$timepoint)
conditions <- unique(annotated_seurat$condition)


annotated_seurat = qs_read(paste0(data_dir_integrated,'de_seurat_annotated.qs'),nthreads=nthreads) # read in data

DefaultAssay(annotated_seurat)
Idents(annotated_seurat) = 'celltype'

table(annotated_seurat$sample, annotated_seurat$celltype)


# remove clusters with only one condition
annotated_seurat = subset(annotated_seurat, 
                          idents=c('Muscle'),
                          invert=TRUE)


table(annotated_seurat$sample, annotated_seurat$celltype)

annotated_seurat$celltype = Idents(annotated_seurat) # reassign with clusters removed


# Run differential proportion analysis

p = propeller(clusters = annotated_seurat$celltype, sample = annotated_seurat$sample, 
              group = annotated_seurat$condition, transform='logit')
p

p = propeller(clusters = annotated_seurat, sample = annotated_seurat$sample, 
              group = annotated_seurat$condition, transform='logit')
p
write.csv(p,pasteo(output_stats_dir,'propeller_props.csv'),row.names=FALSE)

plotCellTypeProps(clusters = annotated_seurat$celltype, sample = annotated_seurat$sample)


ggsave(paste0(output_fig_dir,'celltype/cell_type_proportions_bars_bysample.png'), dpi=320, width=20, height= 10)

plotCellTypeProps(clusters = annotated_seurat$celltype, sample = annotated_seurat$condition)
ggsave(paste0(output_fig_dir,'celltype/cell_type_proportions_bars_bycondition.png'), dpi=320, width=10, height= 10)



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




# epithelial subcluster props ----
epithelial = qs_read(paste0(data_dir,'epithelial.qs'), nthreads = 2)
epithelial$celltype = epithelial$subtype
table(epithelial$sample, epithelial$celltype)

# Create directory for DA analysis results
da_dir <- paste0(stats_dir_integrated, "epithelial/DA_res/")
if (!dir.exists(da_dir)) { dir.create(da_dir, recursive=TRUE) }

# Create directory for DA figures
da_fig_dir <- paste0(stats_dir_integrated, "epithelial/DA_plots/")
if (!dir.exists(da_fig_dir)) { dir.create(da_fig_dir, recursive=TRUE) }


propres <- propeller(epithelial, 
                     clusters = epithelial$subtype,
                     sample=epithelial$sample,
                     group = epithelial$condition,
                     transform='logit')
propres %>%  View()


# Get cell proportions data
cell_props <- get_cell_proportions(epithelial)

# Get all unique cell types
celltypes <- unique(epithelial$subtype)

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

























library(betareg)

## by cluster numbers 0:27
# Cluster Percentages by Pb
freq_table <- prop.table(table(annotated_seurat$seurat_clusters, annotated_seurat$condition), 2)
freq_table

### set up covariates
props <- prop.table(table(annotated_seurat$seurat_clusters,annotated_seurat$condition),2)
treat <- gsub('\\..*','',colnames(props))
props
treat

### beta regression
beta.res<-data.frame(mat.or.vec(22,3))
names(beta.res)<-c("cluster","coef","p-value")
props.zeros <- props
props.zeros[props.zeros==0] <- 1e-6

control_params <- betareg.control(hessian=T)

for (i in 1:22){
  #table(annotated_seurat$celltype==i, annotated_seurat$condition)
  test<-betareg(props.zeros[i,]~factor(treat),
                control=control_params)
  beta.res[i,1]<-i
  beta.res[i,2]<-round(test$coefficients$mean[2],3)
  beta.res[i,3]<-summary(test)$coefficients$mean[2,4]
}

beta.res


## by celltype annotation
# Cluster Percentages by Pb
freq_table <- prop.table(table(annotated_seurat$celltype, annotated_seurat$condition), 2)
freq_table

### set up covars
props <- prop.table(table(annotated_seurat$celltype,annotated_seurat$condition),2)
treat <- gsub('\\..*','',colnames(props))
props
treat

### beta regression
beta.res<-data.frame(mat.or.vec(22,3))
names(beta.res)<-c("cluster","coef","p-value")
props.zeros <- props
props.zeros[props.zeros==0] <- 1e-6

control_params <- betareg.control(hessian=T)

for (i in 1:22){
  #table(annotated_seurat$celltype==i, annotated_seurat$condition)
  test<-betareg(props.zeros[i,]~factor(treat),
                control=control_params)
  beta.res[i,1]<-i
  beta.res[i,2]<-round(test$coefficients$mean[2],3)
  beta.res[i,3]<-summary(test)$coefficients$mean[2,4]
}

beta.res

### combined clusters with same predicted cell type
# olig_combi <- props[3,] + props[11,]
# test_oli <- betareg(olig_combi~treat)
# summary(test_oli)
# predict(test_oli)
# 
# micg_combi <- props[2,] + props[10,]
# test_mcg <- betareg(micg_combi~treat)
# summary(test_mcg)
# predict(test_mcg)
# 
# peri_combi <- props[4,] + props[8,]
# test_peri <- betareg(peri_combi~treat)
# summary(test_peri)
# predict(test_peri)

# beta.res[2,3] <- summary(test_mcg)$coefficients$mean[2,4]
# beta.res[10,3] <- NA
# beta.res[3,3] <- summary(test_oli)$coefficients$mean[2,4]
# beta.res[11,3] <- NA
# beta.res[4,3] <- summary(test_peri)$coefficients$mean[2,4]
# beta.res[8,3] <- NA

### table of cluster % by Pb
freq_table <- data.frame(cluster=rownames(freq_table),Ctrl=freq_table[,1],Pb=freq_table[,2])
tab.clust_x_pb <- cbind(Ncells=as.numeric(table(annotated_seurat$celltype)),freq_table,beta.res[,2:3])
tab.clust_x_pb <- tab.clust_x_pb[,c(2,1,3:6)]
tab.clust_x_pb$Ctrl <- round(100*tab.clust_x_pb$Ctrl,2)
tab.clust_x_pb$Pb <- round(100*tab.clust_x_pb$Pb,2)
tab.clust_x_pb <- tab.clust_x_pb[,-5]
write.csv(tab.clust_x_pb, file=paste0(output_stats_dir,'beta_regression_cluster_props_by_pb.csv'), row.names = F)


###stacked bar chart

#set up data frame
props.m <- t(props)
props.m <- as.data.frame.matrix(props.m)
props.m$mouse <- rownames(props.m)

# #combine 1+9 and 2+10 and 3+7
# props.m$'1' <- props.m$'1' + props.m$'9'
# props.m$'2' <- props.m$'2' + props.m$'10'
# props.m$'3' <- props.m$'3' + props.m$'7'
# props.m <- props.m[,-c(8,10,11)]
# rowSums(props.m[,-10])
# names(props.m)[1:9] <- c('Endothelia','Microglia','Oligodendrocytes','Pericytes','Astrocytes',
#                          'Oligodendrocyte Progenitors','Choroid Plexus','Neurons','Fibroblasts')

#convert to long
long_props <- gather(props.m, cluster, pct, -mouse)
# long_props$cluster <- factor(long_props$cluster, 
#                              levels= c('Endothelia','Microglia','Oligodendrocytes','Pericytes','Astrocytes',
#                                        'Oligodendrocyte Progenitors','Choroid Plexus','Neurons','Fibroblasts'))


ggplot(data=long_props, aes(x=mouse, y=pct, fill=cluster)) +
  geom_bar(stat="identity")+ 
  scale_fill_discrete(name='') +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_text(hjust=-0.1, angle=-45),
        axis.title=element_text(size=18),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12)) +
  xlab("") + ylab("Proportion of Cells")
ggsave(paste0(output_fig_dir,'celltype/cell_type_proportions_bars_bycondition_betareg.png'), dpi=320, width=10, height= 10)


#### sessionInfo() ----
sessionInfo()