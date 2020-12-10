rename_genes <- function(myseur, new_names, old_names) {
        new_names[new_names==""] <- old_names[new_names==""]
        new_names[duplicated(new_names)] <- old_names[duplicated(new_names)]
        rownames(myseur@assays$RNA@counts) <- new_names
        if (nrow(myseur@assays$RNA@data)==length(new_names)) {
                rownames(myseur@assays$RNA@data) <- new_names
        }
        if (nrow(myseur@assays$RNA@scale.data)==length(new_names)) {
                rownames(myseur@assays$RNA@scale.data) <- new_names
        }
        return(myseur)
}

require(Seurat)

obj <- readRDS("Amby_6_integrated_harmony_plus_analysis.rds")

source("~/R-Scripts/Ensembl_Stuff.R")

oldnames <- rownames(obj)
newnames <- General_Map(oldnames, in.org="Hsap", out.org="Rat", in.name="symbol", out.name="symbol")

obj <- rename_genes(obj, newnames, oldnames)

genes <- read.table("~/scripts/snRNAseq_Pipeline/AmbyGenes.txt", header=F)

obj@meta.data$State <- rep("Norm", nrow(obj@meta.data));
obj@meta.data$State[obj@meta.data$donor %in% c("RAT11-1", "RAT9-3", "RAT9-5")] <- "Dis"

obj@meta.data$C_by_State <- paste(obj@meta.data$Core_clusters, obj@meta.data$State, sep="_")

png("Rat_Amby_Genes_DotPlot.png", width=10, height=8, units="in", res=300)
Seurat::DotPlot(obj, features=genes[,1], group.by="C_by_State")
dev.off()


png("Rat_Amby_Genes_umap.png", width=10, height=8, units="in", res=300)
Seurat::FeaturePlot(obj, features=as.character(genes[,1]), reduction="umap")
dev.off()










