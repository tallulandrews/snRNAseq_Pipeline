require(Seurat)

obj <- readRDS("Sn_vs_sc_merged_integrated_harmony_plus_clusters.rds")

# Compare all these clusterings
require(igraph)
require(gplots)
clust_table <- obj@meta.data[, grepl("^knn_", names(obj@meta.data))]
clust_table <- as.matrix(apply(clust_table,2,as.numeric))
require("proxy")
clust_dists <- proxy::dist(clust_table, method=function(x,y){igraph::compare(x,y,method="vi")}, by_rows=FALSE)
clust_similr1 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="nmi")}, by_rows=FALSE)
clust_similr2 <- proxy::simil(clust_table, method=function(x,y){igraph::compare(x,y,method="adjusted.rand")}, by_rows=FALSE)


# Find robust exemplar clustering(s)
require("apcluster")
require("gplots")
set.seed(18371)

res1 <- apcluster(-1*as.matrix(clust_dists), p=-2)

coarse_lvl <- names(res1@exemplars)[2]
fine_lvl <- names(res1@exemplars)[3]

#manually select which exemplar to use
obj@meta.data$Coarse_clusters <- obj@meta.data[[coarse_lvl]] 
obj@meta.data$Fine_clusters <- obj@meta.data[[fine_lvl]]

#apcluster::heatmap(res1, -1*as.matrix(clust_dists))

png(paste(proj_name,"compare_clusterings_heatmap.png",sep="_"), width=6, height=6, units="in", res=300)
lab <- matrix("", ncol=ncol(clust_table), nrow=ncol(clust_table))
lab[colnames(clust_table)==fine_lvl, colnames(clust_table)==fine_lvl] <- "F"
lab[colnames(clust_table)==coarse_lvl, colnames(clust_table)==coarse_lvl] <- "C"
heatmap.2(as.matrix(clust_dists), trace="none", distfun=function(x){return(as.dist(clust_dists))}, cellnote=lab)
dev.off()


# Visualize the Chosen clusterings
print("plotting1")
png(paste(proj_name,"coarse_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(obj, reduction = "umap", group.by="Coarse_clusters")
dev.off()

png(paste(proj_name,"coarse_tsne.png",sep="_"), width=6, height=6, units="in", res=150)
DimPlot(obj, reduction = "tsne", group.by="Coarse_clusters")
dev.off()

png(paste(proj_name,"fine_umap.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(obj, reduction = "umap", group.by="Fine_clusters")
dev.off()

png(paste(proj_name,"fine_tsne.png", sep="_"), width=6, height=6, units="in", res=150)
DimPlot(obj, reduction = "tsne", group.by="Fine_clusters")
dev.off()

saveRDS(obj, paste(proj_name, "harmony_plus_analysis.rds", sep="_"))

print("cluster anno")
#Use autoannotation to label clusters
cluster_assign <- function(scAnno) {
	scAnno <- factor(scAnno)
	freqs <- table(scAnno)/length(scAnno);
	# specific labels - majority rule
	if (max(freqs) > 0.5) {
		return(levels(scAnno)[freqs == max(freqs)])
	}
	# general labels - majority rule
	scAnno <- as.character(scAnno);
	scAnno[grepl("Tcell", scAnno)] <- "Tcell"
	scAnno[grepl("Bcell", scAnno)] <- "Bcell"
	scAnno[grepl("LSEC", scAnno)] <- "LSEC"
	scAnno[grepl("Macrophage", scAnno)] <- "Macrophage"
	scAnno[grepl("Hep", scAnno)] <- "Hepatocyte"
	scAnno <- factor(scAnno)
	freqs <- table(scAnno)/length(scAnno);
	if (max(freqs) > 0.5) {
		return(levels(scAnno)[freqs == max(freqs)])
	} else {
		return("ambiguous")
	}
}

consistent_cluster_labs <- sapply(split(obj@meta.data$consistent_labs, obj@meta.data$Fine_clusters), cluster_assign)
general_cluster_labs <- sapply(split(obj@meta.data$general_labs, obj@meta.data$Fine_clusters), cluster_assign)

obj@meta.data$cluster_quickanno <- consistent_cluster_labs[obj@meta.data$Fine_clusters];

# Cluster ID labelled figures
agg_coord_by_cluster <- function(coords, clusters) {
	x <- split(seq(nrow(coords)), clusters)
	result <- sapply(x, function(a) apply(coords[a,],2,median))
	return(result)
}

tsne_lab_pos <- agg_coord_by_cluster(obj@reductions$tsne@cell.embeddings, obj@meta.data$Fine_clusters)
umap_lab_pos <- agg_coord_by_cluster(obj@reductions$umap@cell.embeddings, obj@meta.data$Fine_clusters)
lab_id <- colnames(tsne_lab_pos)


require("ggplot2")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]
obj@meta.data$short_cluster_anno <- factor(map_cell_types(obj@meta.data$cluster_quickanno), levels=new_colour_scheme[,1]);
new_colour_scheme <- new_colour_scheme[new_colour_scheme[,1] %in% obj@meta.data$shoart_cluster_anno,]

print("plotting")
png("Autoanno_label_harmony_umap.png", width=7.5, height=6, units="in", res=300)
DimPlot(obj, reduction="umap", group.by="short_cluster_anno", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=lab_id, colour="grey35")
dev.off()


png("Autoanno_label_harmony_tsne.png", width=7.5, height=6, units="in", res=300)
DimPlot(obj, reduction="tsne", group.by="short_cluster_anno", pt.size=.1)+scale_color_manual(values=new_colour_scheme[,2])+annotate("text", x=tsne_lab_pos[1,], y=tsne_lab_pos[2,], label=lab_id, colour="grey35")
dev.off()

# T/NK = 10, 5, + 20
# Endo = 15, 12, 8, 19, 21
# Mac = 9, 11, 13, 25
# Hep Traj = 4, 1, 2, 24, 17, 3, 6, 7, 26
# Doublets? 20
# Hep / Erythroid = 0, 14, 27
