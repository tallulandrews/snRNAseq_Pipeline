obj <- readRDS("Sc_vs_sn_integrated_harmony_plus_analysis_Anno.rds")
OPTS <- list(out_prefix="All_integrated_harmony")

source("~/scripts/LiverMap2.0/Colour_Scheme.R")
require(Seurat);

png(paste(OPTS$out_prefix, "All_Manual_cell_type_umap.png", sep="_"), width=10, height=8, units="in", res=300)
print(Type_DimPlot(obj, type_col="All_Integrated_Manual", reduction="umap", cluster_col="Coarse_clusters"))
dev.off()

pdf(paste(OPTS$out_prefix, "Integrated_sample_umap.pdf", sep="_"), width=10, height=8)
print(DimPlot(obj, group.by="sample", reduction="umap"))
dev.off()

pdf(paste(OPTS$out_prefix, "Integrated_assay_umap.pdf", sep="_"), width=10, height=8)
print(DimPlot(obj, group.by="assay_type", reduction="umap"))
dev.off()

tmp <- par()
pdf(paste(OPTS$out_prefix, "QC_Umi_vs_Gene.pdf", sep="_"), width=8, height=8)
par(mar=c(4,4,1,1))
plot(obj$nCount_RNA, obj$nFeature_RNA, col=c("dodgerblue", "darkmagenta")[factor(obj$assay_type)], pch=16, xlab="UMI per cell", ylab="Genes per cell")
dev.off()
par(mar=tmp$mar)

png(paste(OPTS$out_prefix, "All_Piecewise_cell_type_umap.png", sep="_"), width=10, height=8, units="in", res=300)
Type_DimPlot(obj, type_col="Piecewise_ManualAnno", reduction="umap", cluster_col="Coarse_clusters")
dev.off()


