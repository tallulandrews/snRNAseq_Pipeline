require("Seurat")

set.seed(3921)

# Which do we include in the integrated map?
dir <- "/cluster/projects/macparland/TA/LiverMap2.0/Cleaned"
seurfiles <- c("C41_EmptyOnly.rds", 
	"C41_CST_EmptyOnly.rds",
	"C41_NST_EmptyOnly.rds",
	"C41_TST_EmptyOnly.rds",
	"C58_TST_EmptyOnly.rds",
	"C58_RESEQ_EmptyOnly.rds"
#	"C70_TST_EmptyOnly.rds",
#	"C70_RESEQ_EmptyOnly.rds",
#	"C72_TST_EmptyOnly.rds",
#	"C72_RESEQ_EmptyOnly.rds"
	);

samp_names <- unlist(lapply(strsplit(seurfiles, "_"), function(x){x <- x[c(-length(x))]; return(paste(x, collapse="_"))}))

obj_list <- list()
for (i in 1:length(seurfiles)) {
	n <- samp_names[i];
	obj <- readRDS(paste(dir,seurfiles[i], sep="/"));

	#Fix sample ID, and Donor ID
	obj@meta.data$sample <- obj@meta.data$orig.ident
	obj@meta.data$donor <- sapply(strsplit(as.character(obj@meta.data$sample), "_"), function(x){x[[1]]})

	# save sample specific clusters
	obj@meta.data$sample_specific_clusters <- paste(n, obj@meta.data$seurat_clusters, sep="_")

	# get rid of factors
	metadata_classes <- sapply(1:ncol(obj@meta.data), function(i){class(obj@meta.data[,i])})
	for (j in which(metadata_classes == "factor")) {
		obj@meta.data[,j] <- as.character(obj@meta.data[,j]);
	}

	
	obj <- Seurat::NormalizeData(obj, verbose = FALSE, normalization.method="LogNormalize", scale.factor=10000) 
	obj <- Seurat::ScaleData(obj);
	obj@meta.data$cell_barcode <- colnames(obj);
	obj@meta.data[[7]] <- obj@meta.data[[7]][[1]]
	obj@meta.data$sample <- rep(n, ncol(obj));
	obj@meta.data$cell_ID <- paste(obj@meta.data$sample, obj@meta.data$cell_barcode, sep="_")
	obj_list[[n]] <- obj
}

# Merge Datasets
#### Merging does not merge individually scaled datasets!!

merged_obj <- NULL;
universal_genes <- c(-1)
hvgs <- c();
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	if (i == 1) {
		merged_obj <- obj_list[[i]]
		universal_genes <- as.character(rownames(obj_list[[i]]))
		hvgs <- VariableFeatures(obj_list[[i]]);
	} else {
		merged_obj <- merge(merged_obj, y=obj_list[[i]], add.cell.ids=c("", n), project="LiverMap")
		universal_genes <- intersect(universal_genes, as.character(rownames(obj_list[[i]])))
		hvgs <- c(hvgs, VariableFeatures(obj_list[[i]]));
	}
}

# Merge Scaled Matrices.
all_Scaled <- c();
scaled_cell_ids <- c()
for (i in 1:length(obj_list)) {
        n <- samp_names[i];
	scaled <- obj_list[[i]]@assays$RNA@scale.data;
	scaled_cell_ids <- c(scaled_cell_ids, obj_list[[i]]@meta.data$cell_ID);

	if (i == 1) {
		all_Scaled <- scaled;
	} else {
		all_Scaled <- cbind(all_Scaled, scaled);
	}
}

merged_obj@assays$RNA@scale.data <- all_Scaled

# Keep HVGs seen in at least 2 datasets
hvgs <- unique(hvgs[duplicated(hvgs)])

merged_obj@misc$universal_genes <- universal_genes;
merged_obj@misc$repeated_hvgs <- hvgs;
merged_obj@misc$creation_date <- date();
VariableFeatures(merged_obj) <- hvgs;
saveRDS(merged_obj, "Merged_obj_Healthy_sn_sc.rds")

set.seed(9428)

merged_obj <- merged_obj[rownames(merged_obj) %in% universal_genes,]
merged_obj <- RunPCA(merged_obj, pc.genes = VariableFeatures(merged_obj), 
			npcs = 20, verbose = FALSE)
merged_obj <- RunTSNE(merged_obj, dims = 1:10, verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:10, verbose = FALSE)

png("SN_SC_merged_not_integrated_tsne.png", width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="tsne", group.by="sample", pt.size=0.1)
dev.off();
png("SN_SC_merged_not_integrated_umap.png", width=9, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="sample", pt.size=0.1)
dev.off();
png("SN_SC_merged_not_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="tsne", group.by="consistent_labs", pt.size=0.1)
dev.off();
png("SN_SC_merged_not_integrated_umap_autoanno.png", width=12, height =6, units="in", res=300)
DimPlot(merged_obj, reduction="umap", group.by="consistent_labs", pt.size=0.1)
dev.off();

require("harmony")
set.seed(10131)
merged_obj <- RunHarmony(merged_obj, "donor", plot_convergence = TRUE)

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)
merged_obj <- RunTSNE(merged_obj, reduction = "harmony", dims = 1:20)

png("SN_SC_merged_harmony_integrated_umap.png", width=9, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "umap", group.by = "sample", pt.size = .1)
dev.off();
png("SN_SC_merged_harmony_integrated_tsne.png", width=9, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "tsne", group.by = "sample", pt.size = .1)
dev.off();
png("SN_SC_merged_harmony_integrated_umap_autoanno.png", width=12, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "umap", group.by = "consistent_labs", pt.size = .1)
dev.off();
png("SN_SC_merged_harmony_integrated_tsne_autoanno.png", width=12, height =6, units="in", res=100)
DimPlot(merged_obj, reduction = "tsne", group.by = "consistent_labs", pt.size = .1)
dev.off();
saveRDS(merged_obj, "All_merged_universal_genes_harmony_integrated_v2.rds");

# add harmony dimensions to integrated object?
