# too look into: Monocle DE, diffcyt : https://bioconductor.org/packages/release/bioc/vignettes/diffcyt/inst/doc/diffcyt_workflow.html

tidy_Seurat_output <- function(out) {
	out$score <- out$avg_logFC*(out$pct.1-out$pct.2)
	out <- out[order(out$score, decreasing=T),]
	return(out)
}

args <- commandArgs(trailingOnly=TRUE)
# file
# clustername
# projname
# cluster_col_name (Use_clusters)
# confounder_col_name (donor)
# assay_type

require(Seurat)
obj <- readRDS(args[1])

if (length(args) >= 4) {
	obj@meta.data$Use_clusters <- obj@meta.data[,args[4]]
}
if (length(args) >= 5) {
	obj@meta.data$donor <- obj@meta.data[,args[5]]
}

### Addition to enable use on merged then clustered object 
# Subset
if (length(args) >= 6) {
	obj <- obj[,obj$assay_type==args[6]]
}
###

exclude1 <- is.na(obj@meta.data$Use_clusters) | obj@meta.data$Use_clusters == "";
exclude2 <- is.na(obj@meta.data$donor) | obj@meta.data$donor == "";
obj <- obj[,!(exclude1 | exclude2)];

out_seur_wilcox <- FindMarkers(obj, ident.1=args[2], group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="wilcox")

out_seur_MAST <- FindMarkers(obj, ident.1=args[2], group.by="Use_clusters", logfc.threshold=-Inf, min.pct=0.1, min.diff.pct = -Inf, test.use="MAST", latent.vars=obj$donor)

all_DE <- list(seur_wilcox=out_seur_wilcox, seur_mast=out_seur_MAST)
saveRDS(all_DE, paste(args[3], args[2], "DE.rds", sep="_"))

## Pseudobulk ##
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

bulks <- get_pseudobulk(obj@assays$RNA@counts, 
			factor(obj@meta.data$Use_clusters), 
			factor(obj@meta.data$donor))
labs <- strsplit(colnames(bulks), "_")
c_labs <- sapply(labs, function(x){unlist(x[[1]])})
d_labs <- sapply(labs, function(x){unlist(x[[2]])})

bulks_file <- paste(args[3], "_pseudobulk_mat.rds", sep="")
saveRDS(list(mat=bulks, cluster=c_labs, donor=d_labs), bulks_file);


if(!identical(c_labs, d_labs)) {
require("edgeR")
edger_obj <- DGEList(bulks, samples=data.frame(cluster=c_labs, donor=d_labs), group=c_labs)
edger_obj <- calcNormFactors(edger_obj)
design <- model.matrix(~cluster+donor, data=edger_obj$samples)
design <- design[,colSums(design) > 0];

coef_names <- colnames(design);

edger_file <- paste(args[3], "_pseudobulk_edger.rds", sep="")
if (!file.exists(edger_file)) {
        edger_obj <- estimateDisp(edger_obj, design)
        fit <- glmQLFit(edger_obj, design)
        saveRDS(fit, edger_file);
} else {
        fit <- readRDS(edger_file)
}

contrast_vec <-rep(0, length(coef_names));
this_cluster <- paste("cluster", args[2], sep="");
contrast_vec[grepl("cluster", coef_names)] <- -1;
contrast_vec[1] <- -1
if (this_cluster %in% coef_names) {
        contrast_vec[coef_names==this_cluster]<- 1;
} else {
        contrast_vec[1] <- 1;
}

this_vs_all <- glmQLFTest(fit, contrast=contrast_vec)
out <- topTags(this_vs_all, nrow(obj));
all_DE <- list(seur_wilcox=out_seur_wilcox, seur_mast=out_seur_MAST, edger=out)
saveRDS(all_DE, paste(args[3], args[2], "_pseudobulk_DE.rds", sep="_"))
}

#######   END ########


