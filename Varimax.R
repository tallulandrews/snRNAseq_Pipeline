merged_obj <- readRDS("All_merged_universal_genes_harmony_integrated_v2.rds")
require(Seurat)
merged_obj@reductions

harmony_loading <- merged_obj@reductions$harmony@feature.loadings
pca_loading <- merged_obj@reductions$pca@feature.loadings

dim(pca_loading)
dim(harmony_loading)

harmony_loading <- harmony_loading[match(rownames(pca_loading), rownames(harmony_loading)),]

harmony_loading <- harmony_loading/rowSums(harmony_loading)
pca_loading <- pca_loading/rowSums(pca_loading)





set.seed(1029)
test_pca <- varimax(pca_loading)
dim(test_pca$rotmat)

source("../scripts/LiverMap2.0/Colour_Scheme.R")

new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
merged_obj@meta.data$consistent_labs <- map_cell_types(merged_obj@meta.data[,"marker_labs"])
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% merged_obj@meta.data[,"marker_labs"],]

cell_colours <- new_colour_scheme[factor(merged_obj@meta.data[,"marker_labs"]),2]
assay_colours <- c("salmon", "dodgerblue")[factor(merged_obj@meta.data[,"assay_type"])]

require(scales)
sample_col_scheme <- hue_pal()(length(unique(merged_obj@meta.data$sample)))
sample_colours <- sample_col_scheme[factor(merged_obj@meta.data$sample)]

rotated_pca <- merged_obj@reductions$pca@cell.embeddings %*% test_pca$rotmat

par(mfrow=c(1,3))
dim1 = rotated_pca[,19]
dim2 = rotated_pca[,20]
plot(dim1, dim2, pch=16, cex=0.1, col=cell_colours)
plot(dim1, dim2, pch=16, cex=0.1, col=assay_colours);
legend("topright", levels(factor(merged_obj@meta.data[,"assay_type"])), fill=c("salmon", "dodgerblue"))
plot(dim1, dim2, pch=16, cex=0.1, col=sample_colours); 
legend("topright", levels(factor(merged_obj@meta.data$sample)), fill=sample_col_scheme, bty="n")


#Dim associations
# 1=hepatocytes?
# 2=LSECs / cholangio?
# 3=Macrophage
# 4=Stellate
# 5=portal hep
# 6=Cholangiocyte
# 7=?
# 8=high var in sn (Pericentral vs periportal hep?)
# 9=T cell
# 10=LSECs
# 11=hep
# 12=Mac (infl?)
# 13=?
# 14=LSEC.
# 15=?
# 16=cholangiocyte?
# 17=Endothelial vs Hepatocyte
# 18=Macrophage
# 19=?
# 20=?

source("../scripts/GSEA_code.R")

require(fgsea)
require(ggplot2)
require(fastmatch)
require(dplyr)
gene_set <- read_gmt("../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")

for (d in 1:ncol(test_pca$loadings)) {
	scores <- signif(test_pca$loadings[,d], digits=1)
	pval=2
	out <- fgsea::fgsea(pathways = gene_set, 
                           stats = scores,
                           minSize=10,
                           maxSize=1000,
                           nperm=10000) %>% 
                  as.data.frame() %>% 
                  dplyr::filter(padj < !!pval) %>% print
}













set.seed(1029)
test_harmony <- varimax(harmony_loading)
dim(test_harmony$rotmat)

source("../scripts/LiverMap2.0/Colour_Scheme.R")

new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
merged_obj@meta.data$consistent_labs <- map_cell_types(merged_obj@meta.data[,"marker_labs"])
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% merged_obj@meta.data[,"marker_labs"],]

cell_colours <- new_colour_scheme[factor(merged_obj@meta.data[,"marker_labs"]),2]
assay_colours <- c("salmon", "dodgerblue")[factor(merged_obj@meta.data[,"assay_type"])]

require(scales)
sample_col_scheme <- hue_pal()(length(unique(merged_obj@meta.data$sample)))
sample_colours <- sample_col_scheme[factor(merged_obj@meta.data$sample)]

rotated_harmony <- merged_obj@reductions$harmony@cell.embeddings %*% test_harmony$rotmat

par(mfrow=c(1,3))
dim1 = rotated_harmony[,1]
dim2 = rotated_harmony[,9]
plot(dim1, dim2, pch=16, cex=0.1, col=cell_colours)
plot(dim1, dim2, pch=16, cex=0.1, col=assay_colours);
legend("topright", levels(factor(merged_obj@meta.data[,"assay_type"])), fill=c("salmon", "dodgerblue"))
plot(dim1, dim2, pch=16, cex=0.1, col=sample_colours); 
legend("topright", levels(factor(merged_obj@meta.data$sample)), fill=sample_col_scheme, bty="n")


#Dim associations
# 1=hepatocytes?
# 2=LSECs
# 3=Macrophage
# 4=Stellate
# 5=?
# 6=Cholangiocyte
# 7=?
# 8=high var in sn (Pericentral?)
# 9=T cell
# 10=LSECs?
# 11=?
# 12=Mac diversity?
# 13=?
# 14=LSEC diversity?
# 15=?
# 16=?
# 17=mac diver?
# 18=mac diver?
# 19=LSEC diver?
# 20=?

source("../scripts/GSEA_code.R")

require(fgsea)
require(ggplot2)
gene_set <- read_gmt("../ExternalData/BaderLab25Aug2020/Human_MSigdb_August_01_2020_symbol.gmt.txt")

for (d in 1:ncol(test_harmony$loadings)) {
	scores <- signif(test_harmony$loadings[,d], digits=1)
	pval=2
	out <- fgsea::fgsea(pathways = gene_set, 
                           stats = scores,
                           minSize=10,
                           maxSize=1000,
                           nperm=10000) %>% 
                  as.data.frame() %>% 
                  dplyr::filter(padj < !!pval)
}
