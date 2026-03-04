library(Seurat)
library(DropletUtils)
library(BiocParallel)
library(ggplot2)
library(ggpubr)
library(scales)
library(qs2)

set.seed(2024)


pb_cb = qs_read("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/merged_seurat_cb.qs",nthreads=2)
pb_raw = qs_read("C:/Users/david/Documents/Research/PhD/scRNAseq/analysis/3_weeks/cellbender_analysis/FPR_0.0_MAD/data/merged_seurat_raw.qs",nthreads=2)

sum(pb_raw$nCount_RNA < 100)

pb_raw = GetAssayData(subset(pb_raw, subset=(sample == 'wk3_pb32')), assay='RNA',layer='counts')

gc() # free RAM
br.out <- barcodeRanks(pb_raw)

# Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))



gc()


e.out <- emptyDrops(pb_raw, 
                    lower=100, 
                    test.ambient=TRUE,
                    niters=10000,
                    # ignore=1,
                    retain=Inf,
                    # retain=1000,
                    BPPARAM = SnowParam(workers = 4))


# As per the documentation, if the null hypothesis were true, p-values for 
# low-count barcodes should have a uniform distribution.Any strong peaks 
# in the p-values near zero indicate that emptyDrops is not controlling the FDR 
# correctly.
hist(e.out[which(e.out$Total < 100),]$PValue,breaks = 100,
     main = 'EmptyDrops null hypothesis p-values',
     xlab='p-value')

is.cell <- e.out$FDR <= 0.01
e.out$is.cell = is.cell
sum(is.cell, na.rm = TRUE)
table(Limited=e.out$Limited, Significant=is.cell)

# test_cb = subset(merged_seurat_cb, subset=(sample == 'wk3_pb32'))
test_cb = subset(pb_cb, subset=(sample == 'wk3_pb32'))

sum(test_cb$nCount_RNA <= 100) 
# 571

sum(test_cb$nCount_RNA > 100) 
# 3930


test = subset(test_cb, cells = WhichCells(test_cb, cells = rownames(e.out[which(e.out$is.cell==T),])))
dim(test)

sum(test$nCount_RNA <= 100)
# 241

sum(test$nCount_RNA > 100)
# 3654

# 32344  3902 we have 3902 cells after cellbender and emptyDrops processing. cellranger returned ~ 1100. v7.1 returned ~ 500

# test = subset(pb_cb, cells = WhichCells(pb_cb, cells = rownames(e.out[which(e.out$is.cell==T),])))

dim(test)
# WhichCells(pb_cb, expression = nCount_RNA > 100)

# plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     # xlab="Total UMI count", ylab="-Log Probability",log='xy')


# Create a more advanced plot
ggplot(e.out, aes(x = Total, y = -LogProb)) +
  # Add points with color based on cell classification
  geom_jitter(aes(color = is.cell), alpha = 0.6, size = 1.2) +
  
  # Use log10 scale for x-axis with comma formatting for readability
  scale_x_log10(labels = comma, n.breaks = 6) +
  scale_y_log10(labels = comma, n.breaks = 6) +
  annotation_logticks(color='black') +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black"),
    labels = c("TRUE" = "Cell", "FALSE" = "Background"),
    name = "Classification"
  ) +
  labs(
    title = "Cell Classification by UMI Count and Probability",
    subtitle = "Cells with < 100 UMI counts considerd empty",
    x = "Total UMI count (log scale)",
    y = "-Log Probability"
  ) +
  
  # Add a vertical reference line at a threshold value
  geom_vline(xintercept = 100, 
             linetype = "dashed", color = "blue", alpha = 0.7) +
  theme_classic2() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )


