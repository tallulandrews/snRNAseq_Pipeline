source("../../scripts/LiverMap2.0/Colour_Scheme.R")

obj <- readRDS("Sc_vs_sn_integrated_harmony_plus_analysis_Anno.rds")


tmp <- table(obj@meta.data$All_Integrated_Manual, obj@meta.data$Coarse_clusters)
 
cluster2colour <- cbind(colnames(tmp), apply(tmp, 2, 
				function(x){rownames(tmp)[which(x == max(x))]}))



cluster2colour <- cbind(cluster2colour, Cell_type_colours[match(map_cell_types(cluster2colour[,2]), Cell_type_colours[,1]),2])
#factor_level <- c(PH,PH,PH,CH,IH,IH,CH,NKT,CH,LS,LS,IH,ST,B,PH,CH,IM,NM,PE,CH,CHOL,ST,NM,LS,PH,PH,CHOL,PH,CH,NM,CH,LS,PH,ER)
factor_level <- c(1,2,3,12,9,10,13,19,14,20,21,11,25,33,4,15,29,30,24,16,27,26,31,22,5,6,28,7,17,32,18,23,8,34)
cluster2colour <- cbind(cluster2colour, factor_level)

barplot_data <- table(factor(obj@meta.data$Coarse_clusters, level=cluster2colour[order(as.numeric(cluster2colour[,4])),1]), 
				obj@meta.data$sample)


# need to refactor clusters to group by cell-type
png("sn_vs_sc_celltype_by_coarse_cluster.png", width=8, height=6, units="in", res=300)
par(mar=c(4,6,1,1))
barplot((barplot_data), col=cluster2colour[order(as.numeric(cluster2colour[,4])),3], las=1, horiz=T, xlab="N Cells/Nuclei")
dev.off()


barplot_data <- barplot_data[,c("C41", "C58_RESEQ", "C70_RESEQ", "C72_RESEQ", 
			"C41_CST", "C41_NST", "C41_TST", "C58_TST", "C70_TST", "C72_TST")]

seurat_colours <- get_seurat_colours(obj, "assay_type")

donor_colours <- c(rep(seurat_colours[1], 4), rep(seurat_colours[2], 6))

donor_colours <- c("#F8766D", "#F8766DBB", "#F8766D88", "#F8766D44", 
				"#00BFC4", "#00BFC4AA", "#00BFC488", "#00BFC466",
				"#00BFC444", "#00BFC422")    

png("sn_vs_sc_assay_by_coarse_cluster.png", width=14, height=6, units="in", res=300)
par(mar=c(4,6,1,1))
barplot(t(barplot_data), col= donor_colours, ylab="N Cells/Nuclei")
dev.off()

pdf("sn_vs_sc_assay_by_coarse_cluster.pdf", width=14, height=6)
par(mar=c(4,6,1,1))
barplot(t(barplot_data), col= donor_colours, ylab="N Cells/Nuclei")
dev.off()