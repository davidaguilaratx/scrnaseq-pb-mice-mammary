library(rhdf5)
library(Seurat)
library(qs2)

cellbender_output_dir <- "C:/Users/david/Documents/Research/PhD/scRNAseq/cellbender/cellbender_h5_outputs/FPR_0.0/"

# cellbender_original_h5_path <- paste0(cellbender_output_dir,
#                        "output_cellbender_11389-JC-28_FPR_0.0.h5")

cellbender_original_h5_path <- paste0(cellbender_output_dir,
                                      "output_cellbender_8651-JC-32_FPR_0.0.h5")

analysis = 'cellbender_analysis'
# selected_timepoint <- 'month18'
data_dir = paste0("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/18_months/",analysis,"/FPR_0.0_MAD/data/")

ctrl28_cb = qs_read(paste0(data_dir, 'merged_seurat_cb.qs'),nthreads = 2)
ctrl28_cb = subset(ctrl28_cb, subset = (sample == 'ctrl28')) 
ctrl28_cr = qs_read(paste0(data_dir, 'merged_seurat_fltrd.qs'),nthreads = 2)
ctrl28_cr = subset(ctrl28_cr, subset = (sample == 'ctrl28'))


h5ls(cellbender_original_h5_path, all=T)
# group                               name         ltype cset       otype num_attrs  dclass          dtype  stype rank        dim     maxdim
# 0                 /                    droplet_latents H5L_TYPE_HARD    0   H5I_GROUP         3                                  0                      
# 1  /droplet_latents                background_fraction H5L_TYPE_HARD    0 H5I_DATASET         3   FLOAT H5T_IEEE_F64LE SIMPLE    1      30000      30000
# 2  /droplet_latents        barcode_indices_for_latents H5L_TYPE_HARD    0 H5I_DATASET         3 INTEGER  H5T_STD_I64LE SIMPLE    1      30000      30000
# 3  /droplet_latents                   cell_probability H5L_TYPE_HARD    0 H5I_DATASET         3   FLOAT H5T_IEEE_F64LE SIMPLE    1      30000      30000
# 4  /droplet_latents                          cell_size H5L_TYPE_HARD    0 H5I_DATASET         3   FLOAT H5T_IEEE_F64LE SIMPLE    1      30000      30000
# 5  /droplet_latents                 droplet_efficiency H5L_TYPE_HARD    0 H5I_DATASET         3   FLOAT H5T_IEEE_F64LE SIMPLE    1      30000      30000
# 6  /droplet_latents           gene_expression_encoding H5L_TYPE_HARD    0 H5I_DATASET         3   FLOAT H5T_IEEE_F64LE SIMPLE    2 64 x 30000 64 x 30000
# 7                 /                     global_latents H5L_TYPE_HARD    0   H5I_GROUP         3                                  0                      
# 8   /global_latents                 ambient_expression H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1      32344      32344
# 9   /global_latents            cell_size_lognormal_std H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SCALAR    0      ( 0 )      ( 0 )
# 10  /global_latents   empty_droplet_size_lognormal_loc H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SCALAR    0      ( 0 )      ( 0 )
# 11  /global_latents empty_droplet_size_lognormal_scale H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SCALAR    0      ( 0 )      ( 0 )
# 12  /global_latents      swapping_fraction_dist_params H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1          2          2
# 13                /                             matrix H5L_TYPE_HARD    0   H5I_GROUP         3                                  0                      
# 14          /matrix                           barcodes H5L_TYPE_HARD    0 H5I_DATASET         3  STRING     H5T_STRING SIMPLE    1     269455     269455
# 15          /matrix                               data H5L_TYPE_HARD    0 H5I_DATASET         3 INTEGER  H5T_STD_I32LE SIMPLE    1    7396752    7396752
# 16          /matrix                           features H5L_TYPE_HARD    0   H5I_GROUP         3                                  0                      
# 17 /matrix/features                       feature_type H5L_TYPE_HARD    0 H5I_DATASET         3  STRING     H5T_STRING SIMPLE    1      32344      32344
# 18 /matrix/features                             genome H5L_TYPE_HARD    0 H5I_DATASET         3  STRING     H5T_STRING SIMPLE    1      32344      32344
# 19 /matrix/features                                 id H5L_TYPE_HARD    0 H5I_DATASET         3  STRING     H5T_STRING SIMPLE    1      32344      32344
# 20 /matrix/features                               name H5L_TYPE_HARD    0 H5I_DATASET         3  STRING     H5T_STRING SIMPLE    1      32344      32344
# 21          /matrix                            indices H5L_TYPE_HARD    0 H5I_DATASET         3 INTEGER  H5T_STD_I32LE SIMPLE    1    7396752    7396752
# 22          /matrix                             indptr H5L_TYPE_HARD    0 H5I_DATASET         3 INTEGER  H5T_STD_I32LE SIMPLE    1     269456     269456
# 23          /matrix                              shape H5L_TYPE_HARD    0 H5I_DATASET         3 INTEGER  H5T_STD_I32LE SIMPLE    1          2      16384
# 24                /                           metadata H5L_TYPE_HARD    0   H5I_GROUP         3                                  0                      
# 25        /metadata                  barcodes_analyzed H5L_TYPE_HARD    0 H5I_DATASET         4  STRING     H5T_STRING SIMPLE    1      30000      30000
# 26        /metadata             barcodes_analyzed_inds H5L_TYPE_HARD    0 H5I_DATASET         4 INTEGER  H5T_STD_I64LE SIMPLE    1      30000      30000
# 27        /metadata                          estimator H5L_TYPE_HARD    0 H5I_DATASET         4  STRING     H5T_STRING SIMPLE    1          1          1
# 28        /metadata             features_analyzed_inds H5L_TYPE_HARD    0 H5I_DATASET         4 INTEGER  H5T_STD_I64LE SIMPLE    1      15586      15586
# 29        /metadata     fraction_data_used_for_testing H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1          1          1
# 30        /metadata learning_curve_learning_rate_epoch H5L_TYPE_HARD    0 H5I_DATASET         4 INTEGER  H5T_STD_I64LE SIMPLE    1        150        150
# 31        /metadata learning_curve_learning_rate_value H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1        150        150
# 32        /metadata           learning_curve_test_elbo H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1         30         30
# 33        /metadata          learning_curve_test_epoch H5L_TYPE_HARD    0 H5I_DATASET         4 INTEGER  H5T_STD_I64LE SIMPLE    1         30         30
# 34        /metadata          learning_curve_train_elbo H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1        150        150
# 35        /metadata         learning_curve_train_epoch H5L_TYPE_HARD    0 H5I_DATASET         4 INTEGER  H5T_STD_I64LE SIMPLE    1        150        150
# 36        /metadata         target_false_positive_rate H5L_TYPE_HARD    0 H5I_DATASET         4   FLOAT H5T_IEEE_F64LE SIMPLE    1          1          1

metadata = h5read(cellbender_original_h5_path, 'metadata')
droplet_latents = h5read(cellbender_original_h5_path, 'droplet_latents')
droplet_latents = as.data.frame(droplet_latents[-c(2,6)]) # don't want gene_expression encoding or barcode indices
rownames(droplet_latents) = metadata$barcodes_analyzed

matrices = h5read(cellbender_original_h5_path, 'matrix')
global_latents = h5read(cellbender_original_h5_path, 'global_latents')
global_latents = as.data.frame(global_latents)[1]
# rownames(global_latents) = matrices$features$id
global_latents$gene = matrices$features$name
global_latents$ensemble_ids = matrices$features$id
global_latents = global_latents[,c(2,1,3)]
# global_latents$ensemble_ids = rownames(global_latents)

all(droplet_latents$barcode_indices_for_latents == metadata$barcodes_analyzed_inds)

test_cb$cells = sub("^[^_]+_[^_]+_", "",test_cb$cells)
test$cells = sub("^[^_]+_[^_]+_", "",test$cells)
colnames(test_cb) = test_cb$cells
colnames(test) = test$cells

test_cb = AddMetaData(test_cb, metadata = droplet_latents)
test = AddMetaData(test, metadata = droplet_latents)


FeatureScatter(test, feature1='cell_probability', feature2='nCount_RNA',pt.size=0.5, shuffle=TRUE)+
  DarkTheme()+
  scale_x_log10(n.breaks=8,limits=c(NA,NA), labels=comma_format())+
  scale_y_log10(n.breaks=8, labels=comma_format())+
  annotation_logticks(sides = "bl", colour='white', size=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        axis.title = element_text(size=12), 
        axis.text = element_text(size=12),
        title = element_text(size=18),
        plot.background=element_rect(color='black'))+
  NoGrid()+
  geom_hline(yintercept=150, linetype='dashed',color='red')

VlnPlot(test_cb, features = "cell_probability")

hist(test_cb$cell_probability,breaks=1000,xlim = c(0.98,1))

dim(test_cb)
sum(test_cb$cell_probability >0.90)

summary(test_cb$cell_probability)
summary(test$cell_probability)
t.test(test_cb$cell_probability, test$cell_probability,conf.level = 0.99)
