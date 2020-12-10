source("/cluster/home/tandrews/R-Scripts/GSEA_code.R")

files <- Sys.glob("*de_scoresonly.rds")


all_rank <- c()
for (f in files) {

	out <- readRDS(f)
	ranking <- round(rowMeans(out), digits=2)
	if (length(all_rank)==0) {
		all_rank <- cbind(all_rank, ranking)
		rownames(all_rank) <- names(ranking)
	} else {
		all_genes <- sort(unique(rownames(all_rank), names(ranking)))
		all_rank <- all_rank[match(all_genes, rownames(all_rank)),]
		ranking <- ranking[match(all_genes, names(ranking))]
		all_rank <- cbind(all_rank, ranking);
	}
}

colnames(all_rank) <- files
colnames(all_rank) <- c("Hep_0", "Chol_10", "Hep_11", "PortalEndo_12", "LSEC_13", "Hep_1", "Hep_2", "Hep_3", "LSEC_4", "Hep_5", "Mac_6", "NKT_7", "Hep_8", "Stellate_9")

overall <- rowMeans(all_rank, na.rm=T)

#GSEA or other Pathway tools.

#gsea_input <- data.frame(NAME=rownames(scaled_score_mat), score=scaled_score_mat[,4])
#GSEA(gsea_input, c(1,2), gs.db="../../ExternalData/ReactomePathways.gmt")

require(fgsea)
react_gmt <- fgsea::gmtPathways("/cluster/projects/macparland/TA/ExternalData/GeneSets/Human_Reactome_August_01_2020_symbol.gmt")

GO_gmt <- fgsea::gmtPathways("/cluster/projects/macparland/TA/ExternalData/GeneSets/Human_Reactome_August_01_2020_symbol.gmt")


expr_score <- scaled_score_mat[,4]
expr_score <- expr_score[expr_score != 0]
out_react <- fgsea(pathways=react_gmt, stats=expr_score, minSize=15, maxSize=500)
out_go <- fgsea(pathways=go_gmt, stats=expr_score, minSize=15, maxSize=500)




#out <- apply(react_gmt, 1, function(x){name=paste(x[2], gsub(" ", "_", x[1]), sep="_"); 
#				  genes <- x[3:length(x)]; genes <- genes[genes != ""]; 
#				  l <- list(); l[[name]] <- genes; return(l)})
#pathways <- list()
#for (p in out) {
#	tmp <- unlist(p)
#	names(tmp) <- NULL
#	pathways[[names(p)]] <- tmp
#}

