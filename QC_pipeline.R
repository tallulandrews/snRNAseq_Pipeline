# Read in

#Make mt vs n genes and umap wih mt perc
NPCS=20

require("Seurat")
dat <- Read10X("Healthy/SnRNA-seq_Data_C41/TST/")
name="C41_TST"

dat <- dat[Matrix::rowSums(dat>0) > 10, Matrix::colSums(dat>0)>250]

myseur <- FindVariableFeatures(myseur, selection.method = "vst", nfeatures = 2000)

myseur <- CreateSeuratObject(counts = dat, project = name, min.cells = 10, min.features = 250)
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
# Visualization with UMAP
myseur <- RunUMAP(myseur, dims = 1:NPCS, parallel=FALSE)

png(paste(name, "perMT.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "percent.mt")
dev.off()
png(paste(name, "nFeature.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "nFeature_RNA")
dev.off()
png(paste(name, "nUMI.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeaturePlot(myseur, "nCount_RNA")
dev.off()
png(paste(name, "CCphase.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::DimPlot(myseur, group.by="Phase")
dev.off()

png(paste(name, "MTscatter.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeatureScatter(myseur, feature1="nCount_RNA", feature2="percent.mt")
dev.off()
png(paste(name, "nFeatscatter.png", sep="_"), width=6, height=6, units="in", res=150)
Seurat::FeatureScatter(myseur, feature1="nCount_RNA", feature2="nFeature_RNA")
dev.off()
