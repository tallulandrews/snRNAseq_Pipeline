# Read in
auto_anno_dir = "/cluster/projects/macparland/TA/AutoAnnotation"
script_dir <- "~/scripts/LiverMap2.0/"
source(paste(script_dir, "Setup_autoannotation.R", sep="/"))
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
set.seed(8291)

args <- commandArgs(trailingOnly=TRUE)
# Read in EmptyDrops output.
# EmptyDrops output, sample name, folder with raw data.
NAME=args[2]
dat <- ReadRDS(args[1])
#raw_folder=args[3]


#Make mt vs n genes and umap wih mt perc
NPCS=20
N_GENE_per_CELL <- 100
N_CELL_per_GENE <- 10
# folder of CR output
# prefix for output

require("Seurat")

raw <- dat$raw_mat[,Matrix::colSums(dat$raw_mat) > N_GENE_per_CELL]
dat <- dat$raw_mat[,colnames(dat$raw_mat) %in% colnames(dat$corrected)];
dat <- dat[Matrix::rowSums(dat>0) > N_CELL_per_GENE, Matrix::colSums(dat>0)>N_GENE_per_CELL]

myseur <- CreateSeuratObject(counts = dat, project = NAME, min.cells = N_CELL_per_GENE, min.features = N_GENE_per_CELL)
#perc mito
myseur[["percent.mt"]] <- PercentageFeatureSet(myseur, pattern = "^MT-")

# Norm
myseur <- NormalizeData(myseur);

# Scale
myseur <- ScaleData(myseur);

#CC
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

myseur <- CellCycleScoring(myseur, s.features = s.genes, g2m.features=g2m.genes, set.ident=TRUE)

# HVG
myseur <- FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)

# PCA
myseur <- RunPCA(myseur, features = VariableFeatures(object = myseur))

# Clustering
myseur <- FindNeighbors(myseur, dims = 1:NPCS)
myseur <- FindClusters(myseur, resolution = 0.5)
# Visualization with UMAP
myseur <- RunUMAP(myseur, dims = 1:NPCS, parallel=FALSE)

png(paste(NAME, "perMT.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "percent.mt")
dev.off()
png(paste(NAME, "nFeature.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "nFeature_RNA")
dev.off()
png(paste(NAME, "nUMI.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "nCount_RNA")
dev.off()
png(paste(NAME, "CCphase.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(myseur, group.by="Phase")
dev.off()
png(paste(NAME, "Cluster.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(myseur, group.by="seurat_clusters")
dev.off()

png(paste(NAME, "MTscatter.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeatureScatter(myseur, feature1="nCount_RNA", feature2="percent.mt")
dev.off()
png(paste(NAME, "nFeatscatter.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeatureScatter(myseur, feature1="nCount_RNA", feature2="nFeature_RNA")
dev.off()


# Autoannotation
myseur <- run_scmap_seurat(myseur, scmap_ref=map1_ref);
anno_out <- Use_markers_for_anno(myseur@assays$RNA@counts, myseur$seurat_clusters)
myseur@meta.data$marker_anno <- anno_out$cell_assign

simplify_annotations <- function(annotations) {
        simplified <- as.character(annotations)
        simplified[simplified %in% c(
                "AntibodysecretingBcells",
                "MatureBcells")] <- "Bcells"
        simplified[simplified %in% c(
                "CD3abTcells", "gdTcells1", "gdTcells2")] <- "Tcells"
        simplified[simplified %in% c(
                "PericentralHep", "UnidentifiedHep", "PeriportalHep",
                "interzonalHep")] <- "Hepatocyte"
        return(simplified);
}


general_labs <- as.character(simplify_annotations(myseur@meta.data$scmap_cell_anno))
general_labs2 <- as.character(simplify_annotations(myseur@meta.data$scmap_cluster_anno))
myseur@meta.data$marker_simplfied_anno <- simplify_annotations(myseur@meta.data$marker_anno)
general_labs[general_labs != general_labs2] <- "ambiguous"
general_labs[general_labs == "unassigned"] <- "ambiguous"
myseur@meta.data$general_labs <- general_labs;
inconsistent <- myseur@meta.data$scmap_cell_anno != myseur@meta.data$scmap_cluster_anno;
myseur@meta.data$consistent_labs <- as.character(myseur@meta.data$scmap_cell_anno)
myseur@meta.data$consistent_labs[inconsistent] <- as.character(general_labs[inconsistent])

require("ggplot2")
source("~/scripts/LiverMap2.0/Colour_Scheme.R")
new_colour_scheme <- Cell_type_colours[order(Cell_type_colours[,1]),]

myseur@meta.data$short_cluster_anno <- factor(map_cell_types(myseur@meta.data$consistent_labs), levels=new_colour_scheme[,1]);

myseur@meta.data$short_marker_anno <- factor(map_cell_types(myseur@meta.data$marker_simplfied_anno), levels=new_colour_scheme[,1]);

print("plotting")

tmp_colours <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data$short_cluster_anno,]
png(paste(NAME,"Autoanno.png", sep="_"), width=7.5, height=6, units="in", res=300)
DimPlot(myseur, reduction="umap", group.by="short_cluster_anno", pt.size=.1)+scale_color_manual(values=tmp_colours[,2])#+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=lab_id, colour="grey35")
dev.off()

tmp_colours <- new_colour_scheme[new_colour_scheme[,1] %in% myseur@meta.data$short_marker_anno,]
png(paste(NAME,"Autoanno_marker.png", sep="_"), width=7.5, height=6, units="in", res=300)
DimPlot(myseur, reduction="umap", group.by="short_marker_anno", pt.size=.1)+scale_color_manual(values=tmp_colours[,2])#+annotate("text", x=umap_lab_pos[1,], y=umap_lab_pos[2,], label=lab_id, colour="grey35")
dev.off()

# SoupX
require("SoupX")
require("Seurat")
set.seed(3918)
# Create SoupX object
#raw <- Read10X(raw_folder)
raw <- raw[rownames(raw) %in% rownames(myseur),]
keep_cells <- colnames(myseur);
tot_umi <- Matrix::colSums(raw);
tot_is_cell <- min(tot_umi[colnames(raw) %in% colnames(myseur)])
raw_is_cell <- raw[,keep_cells]
myseur <- myseur[,match(colnames(raw_is_cell), colnames(myseur))]

sc <- SoupChannel(raw, raw_is_cell, metaData=myseur@meta.data, soupRange=c(10, tot_is_cell*0.9))

anno_cluster <- cell_anno_to_cluster_anno(myseur@meta.data$consistent_labs, myseur@meta.data$seurat_clusters)
anno_cluster <- anno_cluster$lab[match(myseur@meta.data$seurat_clusters, anno_cluster$cluster)]

sc <- setClusters(sc, anno_cluster);
sc <- setDR(sc, cbind(myseur@reductions$umap@cell.embeddings[,1], myseur@reductions$umap@cell.embeddings[,2]));

# Use known genes to estimate contamination fraction per cell
non_expressed_gene_list = list(HB =c("HBB", "HBA1", "HBA"), TCR=c("TRAC", "TRBC1", "TRBC2", "TRDC",
"TRGC1", "TRGC2"), Drug=c("CYP2E1", "CYP3A7", "CYP2A7", "CYP2A6"), BCR=c("IGKC", "IGHE", "IGHM", "IGLC1", "IGLC3", "IGLC2"), Clot=c("F10", "F5", "F9", "F7", "FGB"), Bile=c("SLC10A1", "SLC01B1", "SLCO1B3", "SLC22A1", "SLC22A7", "ABCB1", "ABCB11", "ABCB4", "ABCC3", "ABCC6"));
non_expressed_gene_list <- lapply(non_expressed_gene_list, function(x){return(x[x %in% rownames(myseur)])})


useToEst = estimateNonExpressingCells(sc, nonExpressedGeneList = non_expressed_gene_list, clusters=NULL)
sc = calculateContaminationFraction(sc, non_expressed_gene_list, useToEst = useToEst, cellSpecificEstimates=TRUE)
quantile(sc$metaData$rho)

# Correcting the expression matrix.
out = adjustCounts(sc)
out <- out[,match(colnames(myseur), colnames(out))]

myseur@assays$SoupCorrected <- out
saveRDS(myseur, paste(NAME,"SoupX_SeurObj.rds", sep="-"));
