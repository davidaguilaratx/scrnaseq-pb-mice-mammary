library(Seurat)
library(SeuratObject)
library(scCustomize)
library(patchwork)
library(cowplot)
library(tidyverse)
library(Matrix)
# library(scDblFinder)
# library(DropletUtils)
library(rhdf5)
library(SoupX)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(plyr)
library(viridis)
library(scales)
library(qs2) # for fast data saving/loading
library(BiocParallel)
# library(metap)
# library(RCurl)
# library(Polychrome)

# set number of threads to use for saving qs formatted data
nthreads <- 2 

# set seed
set.seed(2024)

# create polychrome color palette ----

# p12 = createPalette(12, c("#FF0000", "#00FF00", "#0000FF"), range = c(20, 80))
# swatch(p12)
# names(p12) = NULL

## read in data ----

# data directories

# cellbender data directories ----
cellbender_fulloutput_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs/FPR_0.0/"
cellbender_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs_for_seurat/FPR_0.0/"
# cellbender_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs_for_seurat/FPR_0.01/"

# cell ranger data directory
# cellranger_dir <- "//sph-colacino-win.turbo.storage.umich.edu/sph-colacino/aguilada/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/"





# analysis type, cellbender or cellranger
analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analysis' # we set emptydrops-minimum-umis, 100 and default values for all else

# directories for saving
# selected_timepoint <- 'week3'
# selected_timepoint <- 'week10'
# selected_timepoint <- 'month7'
selected_timepoint <- 'month18'

i = 25
cellranger_dir <- ifelse(selected_timepoint == 'week3', 
                         "D:/PhD/scRNAseq/8651-JC-mammary/10x_analysis_JC-8651_cellranger_v.9.0.0/",
                         "D:/PhD/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0/")



raw_file_path <- paste0(cellranger_dir,
                        "11389-DA-", i, "/count/sample_raw_feature_bc_matrix.h5")

fltrd_file_path <- paste0(cellranger_dir,
                          "11389-DA-", i, "/count/sample_filtered_feature_bc_matrix.h5")


raw_matrix <- Read10X_h5(filename=raw_file_path)
fltrd_matrix <- Read10X_h5(filename=fltrd_file_path)
raw_matrix_subset <- raw_matrix[rownames(raw_matrix) %in% rownames(fltrd_matrix), ]


fltrd_seurat_obj <- CreateSeuratObject(
  counts  = fltrd_matrix,
  assay   = "RNA",
  project = paste0("sample", i)
)


# perform clustering ----
fltrd_seurat_obj = SCTransform(fltrd_seurat_obj)
fltrd_seurat_obj = RunPCA(fltrd_seurat_obj, assay='SCT')
# fltrd_seurat_obj = RunTSNE(fltrd_seurat_obj, dims=1:50,
#                            perplexity=100)
fltrd_seurat_obj = FindNeighbors(fltrd_seurat_obj, dims=1:50)
fltrd_seurat_obj = FindClusters(fltrd_seurat_obj, verbose=T)
fltrd_seurat_obj = RunUMAP(fltrd_seurat_obj, dims=1:50)


DimPlot(fltrd_seurat_obj,
        reduction='umap',
        # reduction='tsne',
        # group.by = 'seurat_clusters',
        # group.by = 'SCT_snn_res.1.4',
        group.by = 'SCT_snn_res.0.8',
        label=T,
        repel=T)

features = c('Adipoq')

FeaturePlot(fltrd_seurat_obj,features = features,
            reduction='umap',
            label=T,
            repel=T)

library(readxl)
# Function to convert Excel data to gene signatures list
excel_to_signatures <- function(excel_file) {
  # Load required libraries
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Please install the 'readxl' package: install.packages('readxl')")
  }
  
  # Read the Excel file
  cell_data <- readxl::read_excel(excel_file)
  
  # Initialize an empty list to store the signatures
  signatures <- list()
  
  # Process each row in the data frame
  for (i in 1:nrow(cell_data)) {
    # Get the cell name and replace spaces with underscores
    cell_name <- gsub(" ", "_", cell_data$cellName[i])
    
    # Process geneSymbolmore1 - split comma-separated genes and trim whitespace
    genes1 <- character(0)
    if (!is.na(cell_data$geneSymbolmore1[i]) && cell_data$geneSymbolmore1[i] != "") {
      genes1 <- strsplit(cell_data$geneSymbolmore1[i], ",")[[1]]
      genes1 <- trimws(genes1)  # Remove leading/trailing whitespace
    }
    
    # Process geneSymbolmore2 (if it exists) - add "-" suffix to each gene
    genes2 <- character(0)
    if (!is.na(cell_data$geneSymbolmore2[i]) && cell_data$geneSymbolmore2[i] != "") {
      genes2 <- strsplit(cell_data$geneSymbolmore2[i], ",")[[1]]
      genes2 <- paste0(trimws(genes2), "-")  # Add "-" suffix and trim whitespace
    }
    
    # Combine both gene lists
    all_genes <- c(genes1, genes2)
    
    # Add to signatures list if there are genes
    if (length(all_genes) > 0) {
      signatures[[cell_name]] <- all_genes
    }
  }
  
  # Return the signatures list
  return(signatures)
}

signatures = excel_to_signatures('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/marker_genes/scType/ScTypeDB_mammary.xlsx')

fltrd_seurat_obj <- UCell::AddModuleScore_UCell(fltrd_seurat_obj, 
                                      features=signatures, 
                                      name=NULL,
                                      BPPARAM=SnowParam(4))
# 
# fltrd_seurat_obj <- UCell::SmoothKNN(fltrd_seurat_obj,
#                            signature.names = names(signatures),
#                            reduction="pca",
#                            BPPARAM = SnowParam(4))

DefaultAssay(fltrd_seurat_obj) = 'SCT'
FeaturePlot(fltrd_seurat_obj, 
            slot = 'data',
            # min.cutoff = 'q10',
            # max.cutoff = 'q80',
            reduction='umap',
            features = c('Krt8'), # Ccl21a is lymphatic endothelial cell marker. Krt6a, luminal stem/progenitor?
            )+
  scale_color_viridis(option='inferno')+
  DarkTheme()

# assign cluster labels with scType ----
DefaultAssay(fltrd_seurat_obj) = "SCT"
# fltrd_seurat_obj = NormalizeData(fltrd_seurat_obj)
# seurat_umap_sct_rna_sc = ScaleData(fltrd_seurat_obj, features = rownames(fltrd_seurat_obj))
seurat_package_v5 = isFALSE('counts' %in% names(attributes(fltrd_seurat_obj[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))


# fltrd_seurat_obj$seurat_clusters = fltrd_seurat_obj$SCT_snn_res.1.4
fltrd_seurat_obj$seurat_clusters = fltrd_seurat_obj$SCT_snn_res.0.2



scRNAseqData_scaled = if (seurat_package_v5) as.matrix(fltrd_seurat_obj[["SCT"]]$scale.data) else as.matrix(fltrd_seurat_obj[["SCT"]]@scale.data)

# run ScType ----
es.max = sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_results = do.call("rbind", lapply(unique(fltrd_seurat_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(fltrd_seurat_obj@meta.data[fltrd_seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(fltrd_seurat_obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

sctype_scores = sctype_scores[order(sctype_scores$cluster),] # order by cluster id#

fltrd_seurat_obj@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  fltrd_seurat_obj@meta.data$sctype_classification[fltrd_seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# plot sctype annotations
DimPlot(fltrd_seurat_obj, 
        reduction = "umap", 
        label = TRUE,
        shuffle=TRUE,
        repel=TRUE,
        group.by = 'sctype_classification') 

Idents(fltrd_seurat_obj) = 'sctype_classification'

fltrd_seurat_obj$seurat_clusters = Idents(fltrd_seurat_obj)

# soup x ----
soup.channel = SoupChannel(raw_matrix_subset, fltrd_matrix)
soup.channel = setClusters(soup.channel,
                           setNames(fltrd_seurat_obj@meta.data$seurat_clusters,
                                    rownames(fltrd_seurat_obj@meta.data)))
soup.channel = setDR(soup.channel, fltrd_seurat_obj@reductions$umap@cell.embeddings)
soup.channel = autoEstCont(soup.channel)
soup = soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ]
soup$gene = rownames(soup)
soup = soup[,c(3,1,2)]


dd = fltrd_seurat_obj@meta.data[colnames(soup.channel$toc), ]
dd$RD1 = fltrd_seurat_obj@reductions$umap@cell.embeddings[,1]
dd$RD2 = fltrd_seurat_obj@reductions$umap@cell.embeddings[,2]

mids = aggregate(cbind(RD1, RD2) ~ seurat_clusters, data = dd, FUN = mean)
gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = seurat_clusters), size = 0.2) + 
  geom_label(data = mids, aes(label = seurat_clusters)) + ggtitle("PBMC 4k Annotation") + 
  guides(colour = guide_legend(override.aes = list(size = 1)))
plot(gg)


dd$krt19 = soup.channel$toc["Krt19", ]
gg = ggplot(dd, aes(RD1, RD2)) + geom_point(aes(colour = krt19 > 0))
plot(gg)

gg = plotMarkerMap(soup.channel, 'Krt19')
plot(gg)


plotMarkerDistribution(soup.channel)

# krtGenes = c('Krt7','Krt8','Krt18','Krt19')
# 
# useToEst = estimateNonExpressingCells(soup.channel, nonExpressedGeneList = list(KRT = krtGenes), 
#                                       clusters = FALSE)
# 
# plotMarkerMap(soup.channel, geneSet = krtGenes, useToEst = useToEst)
# 
# useToEst = estimateNonExpressingCells(soup.channel, nonExpressedGeneList = list(KRT = krtGenes))
# 
# plotMarkerMap(soup.channel, geneSet = krtGenes, useToEst = useToEst)


# soup.channel = calculateContaminationFraction(soup.channel, list(KRT = krtGenes), useToEst = useToEst)
# head(soup.channel$metaData)

# correcting expression profile ----
out = adjustCounts(soup.channel, roundToInt = F)

cntSoggy = rowSums(soup.channel$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(rowSums(soup.channel$toc > out)/rowSums(soup.channel$toc > 0)), n = 20)

plotChangeMap(soup.channel, out, "Krt19")
