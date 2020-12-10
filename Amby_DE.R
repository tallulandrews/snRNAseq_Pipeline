obj <- readRDS("Amby_6_integrated_harmony_plus_analysis.rds")

hepatocyte_DE <- Seurat::FindMarkers(obj, slot="counts", test.use="wilcox",
                                        group.by="Core_clusters", ident.1="3", 
					ident.2=c("0","1","4","7"), logfc.threshold=0)

saveRDS(hepatocyte_DE, "Amby_6_Empty_integrated_HepDE.rds")

Tcell_DE <- Seurat::FindMarkers(obj, slot="counts", test.use="wilcox",
                                        group.by="Core_clusters", ident.1="14", 
					ident.2=c("10"), logfc.threshold=0)



obj@meta.data$State <- rep("Normal", ncol(obj));
obj@meta.data$State[obj@meta.data$sample %in% c("RAT11-1_SN","RAT9-3","RAT9-5_SN")] <- "Disease"

macro <- obj[,obj@meta.data$Core_clusters %in% c("11", "6")]

stellate <- obj[,obj@meta.data$Core_clusters %in% c("5")]
chol <- obj[,obj@meta.data$Core_clusters %in% c("8")]


subset <- chol

DE <- Seurat::FindMarkers(subset, slot="counts", test.use="wilcox", group.by="State", ident.1="Disease", ident.2="Normal")

saveRDS(DE, "Amby_6_Empty_integrated_CholDE.rds")


#### Pathway Enrichment ####

require(fgsea)
react_gmt <- fgsea::gmtPathways("/cluster/projects/macparland/TA/ExternalData/GeneSets/Human_Reactome_August_01_2020_symbol.gmt")

GO_gmt <- fgsea::gmtPathways("/cluster/projects/macparland/TA/ExternalData/GeneSets/Human_GO_AllPathways_no_GO_iea_February_01_2020_symbol.gmt")

DE[,2] <- round(DE[,2], digits=2)
DE <- DE[order(DE$avg_logFC, decreasing=TRUE),]
DE <- DE[abs(DE[,2]) > 0.3,]

require(gprofiler)

cat(rownames(DE)[DE$avg_logFC < 0])





all_pathgene <- unlist(GO_gmt)
DE <- DE[rownames(DE) %in% all_pathgene,]

score <- DE[,2]
names(score) <- rownames(DE)

out_react <- fgsea(pathways=react_gmt, stats=score, minSize=15, maxSize=500)
out_go <- fgsea(pathways=GO_gmt, stats=score, minSize=50, maxSize=5000, nperm=10000)



