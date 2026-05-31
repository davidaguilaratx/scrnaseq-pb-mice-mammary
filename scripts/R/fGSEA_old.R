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

source(here::here('scripts/R/fGSEA_utils.R'))

# set number of threads to use for saving qs formatted data
nthreads = 3

# set seed
set.seed(2024)

# directories for saving
# selected_timepoint <- 'week3'
# selected_timepoint <- 'week10'
# selected_timepoint <- 'month7'
# selected_timepoint <- 'month18'

# analysis type, cellbender or cellranger
analysis = 'cellbender_analysis'
# analysis = 'cellranger.v.9.0.0_analys

de_directory <- paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis_old/integrated/",analysis,"/FPR_0.0_MAD/stats/DESeq2/epithelial/results/")
# directory to save figures and results in ----

major_celltype = 'epithelial'
# directory_figs = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/DESeq2/GSEA/fGSEA/figures/")
# directory_results = paste0("/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/integrated/",analysis,"/FPR_0.0_MAD/stats/DESeq2/GSEA/fGSEA/results/")
directory_fgsea = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis_old/integrated/",analysis,"/FPR_0.0_MAD/stats/DESeq2/",major_celltype,"/GSEA/fGSEA/")

dirs = c(directory_fgsea)

for (d in dirs) {
  if (!dir.exists(d)) { dir.create(d,
                                    recursive=TRUE)}
}

### grab celltypes from DE results ----
### grab de result file names and extract cell types 
de_results = list.files(path=de_directory, pattern='all_genes')

### List of cell types (replace with your actual cell types)
celltypes = unique(sub(pattern='_condition.*',replacement='',x=de_results))
celltypes

# timepoints. 3 weeks, 10 weeks, 7 months, 18 months
timepoints = c('wk3','wk10','mn7','mn18')



# load in data ----
# res_tbl_adip = read.csv(paste0(de_directory,'Adipocytes_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
# res_tbl_endo = read.csv(paste0(de_directory,'Endothelial_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
# res_tbl_macro = read.csv(paste0(de_directory,'Macrophage.Ma_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
# res_tbl_lum = read.csv(paste0(de_directory,'Luminal.HS_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
# res_tbl_bcell = read.csv(paste0(de_directory,'B cells_condition_pb_vs_ctrl_all_genes_deseq2.csv'))

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
gene_id = read.csv('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/marker_genes/ensembl_and_entrez_gene_ids_and_human_homologs.csv')

# test = left_join(res_tbl_adip, msigdbr::msigdbr(species='Mus musculus',category='H'),
#                                                by=join_by('gene'=='gene_symbol'))


# deparse(substitute(variable_name)) # to getvariable name as character string


### merge cell-specific deseq2 data with gene IDs and human homologs ----
# res_tbl_adip = get_ids_and_homologs(res_tbl_adip)
# res_tbl_endo = get_ids_and_homologs(res_tbl_endo)
# res_tbl_macro = get_ids_and_homologs(res_tbl_macro)
# res_tbl_lum = get_ids_and_homologs(res_tbl_lum)
# res_tbl_bcell = get_ids_and_homologs(res_tbl_bcell)


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
# msigdb.mm = getMsigdb(org = 'mm', id = 'SYM') snapshotDate(): 2024-10-24
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
# msigdb.hs = qs_read('/nfs/turbo/sph-colacino/aguilada/scRNAseq/analysis/msigdb.hs.qs',nthreads=2)


# listCollections(msigdb.hs)
# [1] "c1" "c2" "c3" "c4" "c6" "c7" "c8" "h"  "c5"
# listSubCollections(msigdb.hs)
# [1] "CGP"             "CP:BIOCARTA"     "CP:PID"          "CP"             
# [5] "MIR:MIRDB"       "MIR:MIR_LEGACY"  "TFT:TFT_LEGACY"  "CGN"            
# [9] "CM"              "IMMUNESIGDB"     "VAX"             "TFT:GTRD"       
# [13] "HPO"             "CP:REACTOME"     "CP:WIKIPATHWAYS" "GO:BP"          
# [17] "GO:CC"           "GO:MF"           "CP:KEGG"   


# get M5: MPT: tumor phenotype ontology gene sets
mpt_mm = gson::read.gmt('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/m5.mpt.v2024.1.Mm.symbols.gmt')
mpt_mm = split(mpt_mm$gene, mpt_mm$term)
# ccca_hs = gson::read.gmt('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/c4.3ca.v2024.1.Hs.symbols.gmt')
# ccca_hs$gene = str_to_title(ccca_hs$gene)
# ccca_hs = split(ccca_hs$gene, ccca_hs$term)
# kegg_med_hs = gson::read.gmt('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt')
# kegg_med_hs$gene = str_to_title(kegg_med_hs$gene)
# kegg_med_hs = split(kegg_med_hs$gene, kegg_med_hs$term)

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
# c6_mm = geneIds(subsetCollection(msigdb.mm, collection = 'c6'))
hpo_mm = geneIds(subsetCollection(msigdb.mm, subcollection='HPO'))
wiki_mm = geneIds(subsetCollection(msigdb.mm, subcollection='CP:WIKIPATHWAYS'))
immune_mm = geneIds(subsetCollection(msigdb.mm, subcollection="IMMUNESIGDB"))
c6_mm = geneIds(subsetCollection(msigdb.mm, collection='c6'))
c3_mm = geneIds(subsetCollection(msigdb.mm, collection='c3'))
# all_mm = c(hallmarks_mm, gobp_mm, gomf_mm, reactome_mm, kegg_mm, wiki_mm, immune_mm)


# mm_collections = list(
#   'hallmarks_mm' = hallmarks_mm,
#   'gobp_mm' = gobp_mm,
#   'gomf_mm' = gomf_mm,
#   'gocc_mm' = gocc_mm,
#   'reactome_mm' = reactome_mm,
#   'kegg_mm' = kegg_mm,
#   'wiki_mm' = wiki_mm,
#   'immune_mm' = immune_mm,
#   'cgp_mm' = cgp_mm,
#   'hpo_mm'=hpo_mm,
#   'c6oncogenic_mm' = c6_mm,
#   'mpt_mm' = mpt_mm,
#   'c3regulatory_mm' = c3_mm,
#   'ccca_hs' = ccca_hs,
#   'kegg_med_hs' = kegg_med_hs
  
  
#   # 'all_mm' = all_mm
# )

mm_collections = list(
  'hallmarks_mm' = hallmarks_mm,
  'gobp_mm' = gobp_mm,
  'gomf_mm' = gomf_mm,
  'gocc_mm' = gocc_mm,
  'reactome_mm' = reactome_mm,
  'kegg_mm' = kegg_mm,
  'wiki_mm' = wiki_mm,
  'cgp_mm' = cgp_mm
  # 'all_mm' = all_mm
)


### list GSEA pathways of interest
# pathway_subs = c('h','GO:BP','GO:MF','GO:CC','CP:REACTOME', 'CP:KEGG',
#                  'CP:WIKIPATHWAYS',"IMMUNEDB",'CGP','HPO','c6','MPT', 'c3',
#                  '3CA','CP:KEGG_MEDICUS')
# pathway_names = c('Hallmarks','GO: Biological Processes','GO: Molecular Function',
#                   'GO: Cellular Compartment','Reactome','KEGG','Wikipathways', 
#                   "IMMUNEDB","CGP:Chemical and Genetic Perturbations",
#                   'HPO:Human Phenotype Ontology', 'C6:Oncogenic Sigs',
#                   'Mammalian Phenotype Ontology','C3:Regulatory Target Genes','3CA:Curated Cancer Cell Atlas',
#                   'KEGG_MEDICUS')

pathway_subs = c('h','GO:BP','GO:MF','GO:CC','CP:REACTOME','CP:KEGG','CP:WIKIPATHWAYS',"CGP:Chemical and Genetic Perturbations")
pathway_names = c('Hallmarks','GO: Biological Processes','GO: Molecular Function',
                  'GO: Cellular Compartment','Reactome','KEGG','Wikipathways', "CGP:Chemical and Genetic Perturbations")

all(length(pathway_names) == length(mm_collections))

all(length(pathway_names) == length(pathway_subs))

# test = reactome_mm
# test_list = c()
# for (i in seq_along(test)) {
#   # print(length(test[[i]]))
#   test_list = c(test_list, length(test[[i]]))
# }
# summary(test_list)

# t = which(lengths(test) > 5 & lengths(test) < 10)
# t

# de_results
# now = c(1:4,10,23,25:40,67,72)
# de_results = c('Adipocytes','Endothelial','Macrophage.Ma','Luminal.HS','B cells')

# fGSEA on all celltype differential expression results ----
# for (i in 1:12) {
# for (i in c(2,3,8,10,11)) {
# for (i in now) {
for (i in 1:length(de_results)) {
  de_res = de_results[i]
  res_tbl = read.csv(paste0(de_directory,de_res))
  res_tbl = get_ids_and_homologs(res_tbl)
  
  genelist = to_genelist(res_tbl)
  
  time = intersect(timepoints,str_split(de_res,'_', simplify = T))
  
  celltype = intersect(celltypes,str_split(de_results[i],'_', simplify = T))
  
  # create directories for saving figures and results
  save_dir_res = paste0(directory_fgsea,celltype,'/results/')
  if (!dir.exists(save_dir_res)) { dir.create(save_dir_res,
                                          recursive=TRUE)}
  
  save_dir_figs = paste0(directory_fgsea,celltype,'/figures/')
  if (!dir.exists(save_dir_figs)) { dir.create(save_dir_figs,
                                          recursive=TRUE)}
      # for (k in 9) {
  for (k in 1:length(pathway_names)) {
    invisible(gc())
    cat('\n\nSetting up pathways for ',time,celltype,' using ',pathway_names[k],' pathways...\n\n')
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
    
    # pathway database name. Species followed by collection name.
    # assumes first 13 collections are mouse and 14 and beyond are human, adjust threshold from k < 14 as needed
    pathway_name = ifelse(k < 14, paste0(msigdb_species[1],'_',names(mm_collections[k])),
                          paste0(msigdb_species[2],'_',names(mm_collections[k])))
    
    cat('\n\nRunning fGSEA for ',time, celltype,' using ',pathway_name,' pathways...\n\n')
    
    fgseaRes = fgsea(pathways=pathways,
                     stats=genelist,
                     sampleSize=101, # default 101. change to nperm for simple-fGSEA
                     minSize=5,
                     maxSize=1000,
                     nPermSimple = 10000,
                     eps=1e-50, # default 1e-50
                     nproc=nthreads)
    invisible(gc())
    
    if (length(fgseaRes$padj) > 0) {
  
      qs_save(fgseaRes, paste0(save_dir_res,celltype,'_',time,'_',pathway_name,'fgseaRes.qs'),nthreads=nthreads)
      
      ### tidy fGSEA results, sort by NES, and save to csv
      tidysave_fgsea(fgseaRes,
                     pathway_name = pathway_name,
                     celltype=celltype,
                     time=time,
                     directory = save_dir_res)
      
      # select independent pathways to get rid of redundancy
      cat('\ncollapsing pathways for',time,celltype,pathway_name,
          'if not part of 50 hallmarks...\n\n')
      
      collapsedPathways = collapsePathways(fgseaRes=fgseaRes[order(pval)][padj < 0.1], 
                                           pathways=pathways[names(pathways) %in% fgseaRes$pathway], 
                                           stats=genelist,
                                           pval.threshold=0.05, # default 0.05
                                           nperm=10000) # defaults to 200
      
      cat('\nsuccessfully collapsed pathways for',time,celltype,pathway_name,'\n\n')
      
      if (!is.null(collapsedPathways) && length(collapsedPathways$mainPathways) > 1) {
        # save collapedPathways as csv by transforming to dataframe
        qs_save(collapsedPathways, paste0(save_dir_res,celltype,'_',time,'_',pathway_name,'_collapsedpathways.qs'),nthreads=nthreads)
        
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
        
        write.csv(combined_df,paste0(save_dir_res,celltype,'_',time,'_',pathway_name,'_collapsedpathways.csv'))
        
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
        
        png(paste0(save_dir_figs,'_',celltype,'_',time,'_',pathway_name,'_collapsedPathways_hierarchy.png'),
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
        
        
      } else {
        mainPathways = NULL # there are no hierarchical pathways
      }
      
      if ((sum(fgseaRes$padj < 0.1, na.rm=TRUE)) > 0) {
        # padj_cutoff = ifelse(sum(fgseaRes$padj < 0.25, na.rm=TRUE) > 70, 0.1, 0.25)
        padj_cutoff = 0.1
        
        
        ### plot NES results
        cat('\nplotting NES barplot for',time,celltype,'using',pathway_name,'pathways p.adj cuttoff of:',padj_cutoff,'...\n\n')
        # if there are more than 50 significant pathways, plot only the main pathways to avoid overcrowding the plot. Otherwise, plot all significant pathways.
        if(!is.null(mainPathways) && length(mainPathways) > 0 && sum(fgseaRes$padj < 0.1, na.rm = TRUE) > 50) {
          plot_fgsea_nes(fgsea_nes=fgseaRes[pathway %in% mainPathways,],
                         pathways=pathways,
                         time=time,
                         pathway_name = pathway_name,
                         celltype=celltype,
                         padj_cutoff=padj_cutoff,
                         directory = save_dir_figs)
        } else {
          plot_fgsea_nes(fgsea_nes=fgseaRes,
                         pathways=pathways,
                         time=time,
                         pathway_name = pathway_name,
                         celltype=celltype,
                         padj_cutoff=padj_cutoff,
                         directory = save_dir_figs)
        }
        
        cat('\nplotting NES dotplot for',time,celltype,' using ',pathway_name,' pathways p.adj cuttoff of: ',padj_cutoff,'...\n\n')
        if(!is.null(mainPathways) && length(mainPathways) > 0 && sum(fgseaRes$padj < 0.1, na.rm = TRUE) > 50) {
          plot_fgsea_nes_dotplot(fgsea_nes=fgseaRes[pathway %in% mainPathways,],
                                 pathways=pathways,
                                 time=time,
                                 pathway_name = pathway_name,
                                 celltype=celltype,
                                 padj_cutoff=padj_cutoff,
                                 directory = save_dir_figs)
        } else {
          plot_fgsea_nes_dotplot(fgsea_nes=fgseaRes,
                                 pathways=pathways,
                                 time=time,
                                 pathway_name = pathway_name,
                                 celltype=celltype,
                                 padj_cutoff=padj_cutoff,
                                 directory = save_dir_figs)
        }
      } else {
        cat('\nNo pathways with FDR < 0.1 for',time, celltype,pathway_name,'\n')
      }
      
      # cutoff = 0.1
      # i = names(hallmarks_mm) %in% adipo.fgsea.hallmarks.mm$pathway[adipo.fgsea.hallmarks.mm$padj < cutoff]
      
      # update_geom_defaults("segment", list(color = "black"))
      
      # cat('plotting enrichment graphs for ',celltype,' using ',pathway_name,' pathways...\n\n')
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
      # ggsave(paste0(save_dir_figs,celltype,'_multilevel',
      #               '_fGSEA_NES_',pathway_name,'.png'),
      #        plot=p,dpi=320,width=10,height=15,units='in')
    } else {
    cat('\nNo enriched pathways for',time,celltype,'using',pathway_name,'\n')
      }
    }
      
      rm(fgseaRes,collapsedPathways,mainPathways)
      
      cat('\nfinished fGSEA for',time, celltype,'using',pathway_name,'\n')
      invisible(gc())
}






# end of loop -----
fgseaRes <- qs_read(file.path(directory_fgsea,'Basal.Luminal','results','Basal.Luminal_wk10_Mouse_hallmarks_mmfgseaRes.qs'),nthreads=2)


de_res = de_results[3]
res_tbl = read.csv(paste0(de_directory,de_res))
res_tbl = get_ids_and_homologs(res_tbl)

genelist = to_genelist(res_tbl)


plotEnrichment(hallmarks_mm[['HALLMARK_WNT_BETA_CATENIN_SIGNALING']],
               genelist) + labs(title="Glycolysis Enrichment in Macrophages at Week 10") + theme_light()