
all_scale <- c()

for (i in sort(unique(b@meta.data$orig.ident))) {
	tmp <- b[, b@meta.data$orig.ident == i]
	tmp <- NormalizeData(tmp, scale.factor=1500)
	tmp <- Seurat::ScaleData(tmp, features=rownames(tmp));
	all_scale <- cbind(all_scale, tmp@assays$RNA@scale.data)
}


