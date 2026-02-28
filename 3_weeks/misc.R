setwd("C:/Users/david/Documents/Research/PhD/scRNAseq/8651-JC-mammary")

library(Seurat)
library(SeuratObject)
library(ggplot2)
library(viridis)
library(ggthemes)
library(ggpubr)

seurat_umap_sct = readRDS('data/annotated_seurat_Li_etal_res.0.8_df.rds')


# how many cells removed if filter threshold according to 100% cells (dark blue color) from web summary
length(subset(x=seurat_umap_sct,
              subset=(sample == 'ctrl24') &
                (nCount_RNA >= 558))$orig.ident)
# 1683. Down from 1709. 36 cells filtered out.

length(subset(x=seurat_umap_sct,
              subset=(sample == 'pb29') &
                (nCount_RNA >= 2019))$orig.ident)
# 324. Down from 1226. 902 filtered out

length(subset(x=seurat_umap_sct,
              subset=(sample == 'pb31') &
                (nCount_RNA >= 568))$orig.ident)
# 1984. Down from 2091. 107 cells filtered out.

length(subset(x=seurat_umap_sct,
              subset=(sample == 'pb32') &
                (nCount_RNA >= 2431))$orig.ident)
# 24. down from 520. 496 cells filtered out



# easily get list of cluster genes to input in online databases
cell_markers_indata = read.csv('annotation/aggregate_sources_celltype_markers_actually_in_data.csv')

# Aggregating gene symbols into a comma-separated string for each cell type
aggregated_df = pivot_longer(cell_markers_indata, cols = everything(), names_to = "Cell", values_to = "GeneSymbol") %>%
  filter(GeneSymbol != "", !is.na(GeneSymbol)) %>%
  group_by(Cell) %>%
  summarise_all(paste,collapse=',') # %>%
  # mutate(GeneSymbol = toupper(GeneSymbol))

# Aggregating gene symbols into a comma-separated string for each cell type
aggregated_df2 = top10[,1:2] %>%
  group_by(cluster) %>%
  summarise_all(paste,collapse=',') %>%
  mutate(gene = toupper(gene))

aggregated_df3 = top10.conserved[,1:2] %>%
  group_by(cluster_id) %>%
  summarise_all(paste,collapse=',') %>%
  mutate(gene = toupper(gene))

aggregated_df3$gene[aggregated_df3$cluster_id == 17]


Idents(seurat_umap_sct) = 'SCT_snn_res.0.8'
DefaultAssay(seurat_umap_sct) = "RNA"



# seurat_umap_sct = JoinLayers(seurat_umap_sct, assay='RNA')


# Heatmap


# Proliferating cells ----
# cluster 8 and 27
FeaturePlot(seurat_umap_sct, features = c('Top2a','Rrm2','Slbp','Hjurp','Maz'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Top2a','Rrm2','Slbp','Hjurp','Maz'),
        assay='RNA')

# Erythroid-like and erythroid precursor cells ----
# cluster 27
FeaturePlot(seurat_umap_sct, features = c('Hbb-bs','Epb42','Alad','Ppox','Spta1','Slc4a1','Car2','Slc25a37','Samd11','Dmtn'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Top2a','Rrm2','Slbp','Hjurp','Maz'),
        assay='RNA')


##Hmmm... no epithelial cells ----
FeaturePlot(seurat_umap_sct, # overlay gene marker expression on umap. 
            reduction='umap',
            features=c("Trf", "Cited1", "Krt14", "Rgs5", "Csf1r", "Cd4", "Mmp12","Prox1","Epcam","Cdh5","Col1a1"),
            order=T,
            min.cutoff='q10',
            label=T)

FeaturePlot(seurat_umap_sct, features = c("Trf", "Cited1", "Krt14", "Rgs5", "Csf1r", "Cd4", "Mmp12","Prox1","Epcam","Cdh5","Col1a1")) # hmm no epithelial cells


# luminal progenitor ----
# cluster 23
FeaturePlot(seurat_umap_sct, 
            reduction='umap',
            features = c('Tspan8','Elf5','Aldh1a3'),
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()
# Some aldh1a3 expression in 6-fibroblast and 23-luminal.hs

DotPlot(seurat_umap_sct, features = c('Tspan8','Elf5','Aldh1a3'),
        assay='RNA')

# from Pal et al 2021
FeaturePlot(seurat_umap_sct, features = c('Ets2','Basp1','Epcam','Jund','Cldn4','Ly6e','Krt8','Krt18','Clu','Trps1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Btg2','Fosb','Egr1','Klf6','Jun','Tspan8','Igfbp5','Ehf','Aldh1a3',
                                          'Plet1','Ehf','Wfdc18','Kit','Cldn7','Ccnd1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()



# epithelial, luminal HS cells ----
# cluster 23
# Gata3 marker of T cells too. Cdh1?. Tspan8. Cldn3, Cldn4, cldn7. Maybe Kit but also expressed elsewhere
FeaturePlot(seurat_umap_sct, features = c('Krt8','Krt18','Elf5','Epcam','Prlr','Esr1','Pgr'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Krt8','Krt18','Elf5','Epcam','Prlr','Esr1','Pgr'),
        assay='RNA')

# From Saeki et al 2021
# luminal hormone sensing
FeaturePlot(seurat_umap_sct, features = c('Areg','Cited1','Ly6d','Prlr','Esr1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# From Saeki et al 2021
# luminal alveolar
FeaturePlot(seurat_umap_sct, features = c('Elf5'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# adipocytes ----
# clusters 7 and 21
FeaturePlot(seurat_umap_sct, features = c('Plin1','Lipe','Fabp4','Fasn','Lpl'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Plin1','Lipe','Fabp4'),
        assay='RNA')

# white adipose tissue ----
# cluster 7
FeaturePlot(seurat_umap_sct, features = c('Slc7a10','Slc25a10','Retn','Car3','Gsn','Igfbp4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# more metabolically active adipocyte, brown adipocyte ----
# cluster 21
FeaturePlot(seurat_umap_sct, features = c('Cpt1b','Pdk4','Ucp1','Ppara','Gk','Ntrk3','Coq8a'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Plin1','Lipe','Fabp4'),
        assay='RNA')

# Adipose progenitor? ----
FeaturePlot(seurat_umap_sct, features = c('Cidea','Ucp1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c(),
        assay='RNA')

# Preadipocytes ----
# not in adipose cluster, maybe in fibroblast cluster if the are indeed precurusos to fibroblasts?
FeaturePlot(seurat_umap_sct, features = c('Dpp4','Cd34','Pdgfra','Dlk1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Dpp4','Cd34','Pdgfra','Cd29','Dlk1','Cd24','Fabp4','Pparg','Cd9','Egfl6'),
        assay='RNA')


# basal ----
# cluster 18
FeaturePlot(seurat_umap_sct, features = c('Trp63','Krt14','Oxtr','Acta2','Sparc','Krt5','Myh11','Itgb1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Trp63','Krt14','Oxtr','Acta2','Sparc'),
        assay='RNA')

# basal/myoepithelial ----
# cluster 18
# Tagln? Il34? Notch3?
FeaturePlot(seurat_umap_sct, features = c('Krt5','Krt14','Cnn1','Cdh3','Oxtr','Acta2','Snai2','Myh11'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Krt5','Krt14','Cnn1','Cdh3','Oxtr','Acta2','Snai2','Myh11'),
        assay='RNA')

# basal from Pal et al 2021
FeaturePlot(seurat_umap_sct, features = c('Cxcl14','Krt14','Lgals7','Gas1','Vim','Postn','Mylk','Tpm2','Tagln',
                                          'Acta2','Sparc','Myl9','Tpm1','Cald1','Col4a1','Mgp','Igfbp3','Id4','Igfbp7',
                                          'Apoe','Fst','Dcn','Lmo1','Egr3','Apoc1','Pdpn'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Krt5','Krt14','Cnn1','Cdh3','Oxtr','Acta2','Snai2','Myh11'),
        assay='RNA')

# From Saeki et al 2021
# basal
FeaturePlot(seurat_umap_sct, features = c('Krt14','Acta2','Myl9', 'Tagln'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

Idents(seurat_umap_sct) = 'SCT_snn_res.1'
DotPlot(seurat_umap_sct, features = c('Krt14','Acta2','Myl9', 'Tagln'),
        assay='RNA')

Idents(seurat_umap_sct) = 'SCT_snn_res.1'
FeaturePlot(seurat_umap_sct, features = c('Tagln'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

 

# B cells ----
# clusters 0,3,11, and 16
FeaturePlot(seurat_umap_sct, features = c('Cd19','Cd74','Syk','Pax5','Cd79a','Cd79b','Ighd','Igkc','Ighm','Blnk','Ebf1','Tnfrsf13c'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd19','Cd74','Syk','Pax5','Cd79a','Cd79b','Ighd','Igkc','Ighm','Blnk','Ebf1','Tnfrsf13c'),
        assay='RNA')

# activated, transitional or regulatory, or B-T hybrids
# cluster 3. Bcar3 invovled in anti-estrogen resistance
FeaturePlot(seurat_umap_sct, features = c('Egr3','Bcar3'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# memory traits
FeaturePlot(seurat_umap_sct, features = c('Tnfrsf13c','Sell','Pax5'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# T-Bet cells
FeaturePlot(seurat_umap_sct, features = c('Tbx21'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# B cells with type 1 interferon response, viral immunity. maybe aka activated B cells ----
# cluster 16
FeaturePlot(seurat_umap_sct, features = c('Cd19','Cd74','Cd79a','Cd79b','Irf7','Irf9','
                                          Ifit3','Stat1','Stat2','Stat3','Tor3a','Ifi209','Trim30a','Zbp1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Irf7', 'Ifit3'),
        assay='RNA')

# interferon response genes 
FeaturePlot(seurat_umap_sct, features = c('Il15','Socs1','Usp18','Isg15','Csf1','Tnf','Irf7','Irf9'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# viral response interferon
FeaturePlot(seurat_umap_sct, features = c('Stat1','Irf7','Irf9','Usp18'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


FeaturePlot(seurat_umap_sct, features = c('Gata6','Mef2c','Egr2','Zeb1','Ehf'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()



# T cells ----
# clusters 1,2,4,10,13,and 20
FeaturePlot(seurat_umap_sct, features = c('Cd247','Lck','Lat'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd24a','Itgb1','Itga6','Sca1','Cd14','Prom1'),
        assay='RNA')

# CD4+ T cells ----
# cluster 1
FeaturePlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Cd4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Cd4'),
        assay='RNA')

# CD8+ T cells ----
# cytotoxic T cells ready to fight virus
# cluster 2
FeaturePlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Cd8a','Cd8b1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Cd8a','Cd8b1'),
        assay='RNA')

# T regulator cells (Tregs) ----
# cluster 4
FeaturePlot(seurat_umap_sct, features = c('Foxp3','Ikzf2','Ctla4','Il2ra','Il2rb','Lck','Lat','Fyn'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Foxp3','Ikzf2','Ctla4','Il2ra','Lck','Lat','Fyn'),
        assay='RNA')

#  Activated T cells ----
# cluster 10
FeaturePlot(seurat_umap_sct, features = c('Cd4','Cd5','Zap70','Lck','Lat','Cd69',
                                          'Nr4a1','Nr4a3', 'Nfkbid'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd4','Cd5','Zap70','Lck','Lat','Cd69','Nr41a'),
        assay='RNA')

# Th17 helper T cells (T helper 17 cells) ----
# cluster 13
FeaturePlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Rorc','Il17rb'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Rorc','Il17rb'),
        assay='RNA')

# cluster 20
FeaturePlot(seurat_umap_sct, features = c('Usp18','Irf7','Stat1','Xaf1','Zbp1','Trim30a','Ifi35'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd247','Lck','Lat','Rorc','Il17rb'),
        assay='RNA')

FeaturePlot(seurat_umap_sct, features = c('Klrb1c','Klrk1','Ncr1','Il2rb','Itgam','Cd27','Cd19'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# fractionation ----
FeaturePlot(seurat_umap_sct, features = c('Cd24a','Itgb1','Itga6','Sca1','Cd14','Prom1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Cd24a','Itgb1','Itga6','Sca1','Cd14','Prom1'),
        assay='RNA')

# epithelial progenitor? Basal? ----
# Il34 potential marker?
# cluster 24?
FeaturePlot(seurat_umap_sct, features = c('Sox10','Matn2','Cldn19'),
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Sox10','Matn2', 'Cldn19', 'Aatk'),
        assay='RNA')


# Endothelial ----
# clusters 19,25, and 5
FeaturePlot(seurat_umap_sct, features = c('Eng','S1pr1','Emcn'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# cluster 5
FeaturePlot(seurat_umap_sct, features = c('Pecam1','Cav1','Esam','Angptl4','Cxcl12','Hspg2'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# cluster 19
FeaturePlot(seurat_umap_sct, features = c('Pecam1','Egfl7','Eln','Tbx1','Myzap','Dock6',
                                          'Npr1','Il6st','Lyve1','Ccl21','Prox1','Flt4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# cluster 25
FeaturePlot(seurat_umap_sct, features = c('Pecam1','Vwf','Angptl4','Robo4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()



# Important for TEBs and ductal growth/branching ----
FeaturePlot(seurat_umap_sct, features = c('Msx2','Areg','Wnt4', 'Gata3','Foxm1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Msx2','Areg','Wnt4', 'Gata3','Foxm1'),
        assay='RNA')

# mesenchyme ----
FeaturePlot(seurat_umap_sct, features = c('Pth1r','Tbx3','Bmp4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Pth1r'),
        assay='RNA')

# macrophage ----
# clusters 14,15,12, and 22
FeaturePlot(seurat_umap_sct, features = c('Csf1r','Mertk','C1qa','Adgre1', 'Trem2'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Notch1'),
        assay='RNA')

# cluster 12
# Ma macrophages according to scTYPE and Li et al. 2020 genes
FeaturePlot(seurat_umap_sct, features = c('Mrc1','Cd163','Cd209d','Cd209f','Cd209g','Csf1r',
                                          'Csf2ra','C1qb','Clec10a','Stab1','Pltp','Apoe'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# cluster 14
# Mb macrophages according to scTYPE and Li et al. 2020 genes
# might just go with myeloid cell label for now. GPT4 suggests Clec9a and Xcr1 dendritic. These may be active dendritic cells
FeaturePlot(seurat_umap_sct, features = c('Tyrobp','Clec9a','H2-Ab1','H2-Eb1','Cd74','Lyz2',
                                          'Adam8','Fes','Hck','Grk3','Prkcd','Csf2ra'),
reduction='umap',
min.cutoff='q10',
max.cutoff='q50',
cols=viridis(256)) &
  DarkTheme()


# cluster 15
# Mb macrophages according to scTYPE and Li et al. 2020 genes
# might just go with myeloid cell label for now. GPT4 suggests osteoimmune cells or cells
# interacting wth physciology. Maybe those that are interacting with the developing mammary gland?
FeaturePlot(seurat_umap_sct, features = c('Siglech','Irf8','Lair1','Mpeg1','Grn','Cd300c',
                                          'Runx2','Sema4b','Adam11','Fgr','Fyb','Plxdc1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# cluster 22
# Mb macrophages according to scTYPE and Li et al. 2020 genes
# might just go with myeloid cell label for now.
# gpt4 suggests myeloid with both immune modulation, inflammatory and tissue repair contexts
FeaturePlot(seurat_umap_sct, features = c('Axl','Cd33','Mafb','Mertk','Csf1r','C1qb',
                                          'Sirpa','Clec10a','Fcgr2b','Cyth4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# muscle ----
# cluster 26
FeaturePlot(seurat_umap_sct, features = c('Obscn','Ttn','Ryr1','Tpm2','Mylk2','Tnnt3','Tnni3','Pygm','Myh4'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Obscn','Ttn','Ryr1','Tpm2','Mylk2','Tnnt3','Tnni3','Pygm','Myh4'),
        assay='RNA')

# natural killer cells ----
# luster 13?
# Looks like no NK cells....
FeaturePlot(seurat_umap_sct, features = c('Nkp46', 'Nkp30', 'Cd56', 'Klrb1', 'Nkg2d'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c(),
        assay='RNA')

# NK-like
FeaturePlot(seurat_umap_sct, features = c('Tbx21'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# natural killer cells, markers from scType ----
FeaturePlot(seurat_umap_sct, features = c('Nkg7','Gzmb','Ctsw','Gzma','Cd7','Klrk1','Pfn1','Cd2','Cd69','Gzmm','Cox6a2','Ncr'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c(),
        assay='RNA')

# looks like there are NK cells within the upper right portion of cluster 2


# dendritic ----
# cluster 9
FeaturePlot(seurat_umap_sct, features = c('Irf8','Siglech','Clec9a','Csf2ra'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c('Irf8','Siglech','Clec9a','Csf2ra'),
        assay='RNA')


FeaturePlot(seurat_umap_sct, features = c('Fcer1a','Batf3','Relb','H2-M2','H2-Q6','Ccr7','Il15ra'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# fibroblasts ----
# clusters 6, 17, and 24
FeaturePlot(seurat_umap_sct, features = c('Col1a1','Col1a2','Col3a1','Col4a1','Sparc','Hspg2','Mbp'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

DotPlot(seurat_umap_sct, features = c(),
        assay='RNA')

# fibroblast
# cluster 6
FeaturePlot(seurat_umap_sct, features = c('Col1a1', 'Col1a2', 'Col3a1', 'Col4a1', 'Col5a1', 'Col5a3',
                                          'Col6a1', 'Col6a2', 'Col6a3','Fn1','Fgfr1','Serpinh1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# Reticular cells
# cluster 17. activate remodeling. The Ccl19, Cxcl12, and Ccl21a contribute to tertiary lymphoid structures
FeaturePlot(seurat_umap_sct, features = c('Mfge8','Ccl19','Cxcl12','Ccl21a','Clu','Csf2rb','Des','Tpm1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Vegfa','Pdgfrb'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# Fibroblastic reticular cells
FeaturePlot(seurat_umap_sct, features = c('Vcam1','Pdpn','Pdgfra','Lyve1','Relb','Ccl19',
                                          'Cxcl12','Cd34','Bst1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# marginal reticular cells. from Katakai 2012
FeaturePlot(seurat_umap_sct, features = c('Cxcl13','Vcam1'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# nerves? oligodendrocytes? schwann cells? # likely nervous system related due to Kcna1, votlage gated ion channel protein
# cluster 24
# leaning towards schwann cells. Could also be mesnchymal stem/proginitor cells?
FeaturePlot(seurat_umap_sct, features = c('Mpz','Mbp','Plp1','Cnp','Mag','Sox10',
                                          'Aatk','Pmp22','Mog'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Mpz','Sox10','Pmp22','Lgi4','Kcna1','Kcna2','Kcna6','Mbp','Mag',
                                          'Plp1','Adamts20'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Sox10','Notch3','Igf2','Pdgfb','Lifr'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Sox9','Sox10'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# c(basal, alveolar, luminal HS)
FeaturePlot(seurat_umap_sct, features = c('Krt14','Wfdc18','Krt18'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# Mammary stem cell according to Wang etal 2015
FeaturePlot(seurat_umap_sct, features = c('Procr'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# important for development
# according to https://www.youtube.com/watch?v=s7s9LlPHnx4
FeaturePlot(test, features = c('Esr1','Prlr','Fgfr2','Areg','Egfr','Igf1','Gh'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

# also important. Guide axon in underdevelopment and also guide ductal morphogenesis
FeaturePlot(seurat_umap_sct, features = c('Tgfb1','Ntn1','Slit2','Reln'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Cd69','Cd4','Cd247'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()

FeaturePlot(seurat_umap_sct, features = c('Krt5','Krt8'),
            reduction='umap',
            min.cutoff='q10',
            max.cutoff='q50',
            cols=viridis(256)) &
  DarkTheme()


# how many cells per cluster by condition?
table(seurat_umap_sct$SCT_snn_res.0.8, seurat_umap_sct$condition)

# by sample and condition
table(seurat_umap_sct$sample, seurat_umap_sct$condition)

# by phase and condition
table(seurat_umap_sct$Phase, seurat_umap_sct$condition)
