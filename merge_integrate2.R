require("Seurat")
set.seed(3921)
files_healthy <- Sys.glob("*.rds")
files_psc<- "../PSC/PSC004S2\ SeurObj.rds"
seurfiles <- c(files_healthy, files_psc)
samp_names <- unlist(lapply(strsplit(seurfiles, " "), function(x){x[[1]]}))
samp_names[5] <- "PSC004S2"
obj_list <- list()
common_genes <- c();
for (i in 1:length(seurfiles)) {
        n <- samp_names[i];
        obj <- readRDS(seurfiles[i]);
        obj@meta.data$cell_barcode <- colnames(obj);
        if (length(common_genes) == 0) {
                common_genes <- rownames(obj);
        } else {
                common_genes <- intersect(common_genes, rownames(obj));
        }
        obj@meta.data$donor <- rep(n, ncol(obj));
        obj_list[[n]] <- obj
}
common_genes <- sort(common_genes);
common_genes <- common_genes[!grepl("^MT-", common_genes)]
for (i in 1:length(obj_list)) {
        obj_list[[i]] <- obj_list[[i]][match(common_genes, rownames(obj_list[[i]])),]
}
merged_obj <- NULL;
common_genes <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
        if (i == 1) {
                merged_obj <- obj_list[[i]]
        } else {
                merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
        }
}
merged_obj <- Seurat::NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, pc.genes = VariableGenes(merged_obj), npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj@reductions$umap_no_harmony <- merged_obj@reductions$umap
merged_obj@reductions$tsne_no_harmony <- merged_obj@reductions$tsne

png("merged_not_integrated_tsne.png", width=7, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="tsne_no_harmony", group.by="donor", pt.size=0.1)
dev.off();
png("merged_not_integrated_umap.png", width=7, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="umap_no_harmony", group.by="donor", pt.size=0.1)
dev.off();

require("ggplot2")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
tmp_colours <- new_colour_scheme[new_colour_scheme[,1] %in% merged_obj@meta.data$short_marker_anno,]
png("merged_not_integrated_umap_anno.png", width=9, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="umap_no_harmony", group.by="short_marker_anno", pt.size=.1)+scale_color_manual(values=tmp_colours[,2])
dev.off();
png("merged_not_integrated_tsne_anno.png", width=9, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="tsne_no_harmony", group.by="short_marker_anno", pt.size=.1)+scale_color_manual(values=tmp_colours[,2])
dev.off();

require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, "donor", plot_convergence = TRUE)
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj <- RunTSNE(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj@reductions$umap_harmony <- merged_obj@reductions$umap
merged_obj@reductions$tsne_harmony <- merged_obj@reductions$tsne

png("merged_harmony_integrated_umap.png", width=7, height =6, units="in", res=50)
DimPlot(merged_obj, reduction = "umap_harmony", group.by = "donor", pt.size = .1)
dev.off();
png("merged_harmony_integrated.png", width=7, height =6, units="in", res=50)
DimPlot(merged_obj, reduction = "tsne_harmony", group.by = "donor", pt.size = .1)
dev.off();
png("merged_harmony_integrated_umap_anno.png", width=9, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="umap_harmony", group.by="short_marker_anno", pt.size=.1)+scale_color_manual(values=tmp_colours[,2])
dev.off();
png("merged_harmony_integrated_tsne_anno.png", width=9, height =6, units="in", res=50)
DimPlot(merged_obj, reduction="tsne_harmony", group.by="short_marker_anno", pt.size=.1)+scale_color_manual(values=tmp_colours[,2])
dev.off();
saveRDS(merged_obj, "harmony_integrated.rds");
