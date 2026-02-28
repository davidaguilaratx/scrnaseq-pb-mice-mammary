library(clusterProfiler)
library(enrichplot)
library(DOSE)

set.seed(2024)

de_directory <- "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/results/"

fgsea_directory = "C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/stats_rerun/DESeq2/GSEA/fGSEA/results/"


res_tbl_adip = read.csv(paste0(de_directory,'Adipocytes_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_endo = read.csv(paste0(de_directory,'Endothelial_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_macro = read.csv(paste0(de_directory,'Macrophage.Ma_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_lum = read.csv(paste0(de_directory,'Luminal.HS_condition_pb_vs_ctrl_all_genes_deseq2.csv'))
res_tbl_bcell = read.csv(paste0(de_directory,'B cells_condition_pb_vs_ctrl_all_genes_deseq2.csv'))


fgseaRes = read.csv(paste0(fgsea_directory,'Endothelial_multilevel_fGSEA_NES_Mouse_hallmarks_mm.csv'))


# define some helper functions we'll use ----
get_ids_and_homologs= function(res_tbl, gene_ids=gene_id) {
  res_tbl = merge(x=res_tbl, y=gene_ids, by.x ="gene", by.y="external_gene_name")
}

to_genelist = function(res_tbl, id='gene', species_human=F) { # id should be from c('gene','ensembl_gene_id','entrezgene_id')
  if(!species_human) {
    genelist = res_tbl %>%
      as_tibble() %>%
      # group_by(hsapiens_homolog_associated_gene_name)
      # summarize(pvalue=mean(pvalue)) %>%
      mutate(ranking_metric = sign(log2FoldChange) * -log10(pvalue))  %>%
      dplyr::select(id, ranking_metric) %>%
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




# to map our restuling genes to ensemble and entrez ids
gene_id = read.csv('c:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/annotation/ensembl_and_entrez_gene_ids_and_human_homologs.csv')

res_tbl_endo = get_ids_and_homologs(res_tbl_endo)
genelist_endo = names(to_genelist(res_tbl_endo, id='entrezgene_id'))

ego = enrichGO(genelist_endo, 
               OrgDb = 'org.Mm.eg.db',
               ont='BP',
               qvalueCutoff=0.05,
               minGSSize=15,
               maxGSSize=500)


mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore", showCategory=20)

? enrichGO


