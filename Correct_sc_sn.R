
# Cluster specific DE in each assay-type
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

merged_obj@meta.data$unique_group_id <- paste(merged_obj@meta.data$assay_type, merged_obj@meta.data$donor, merged_obj@meta.data$seurat_clusters, sep="_")

group_means <- group_rowmeans(merged_obj@assays$RNA@data, merged_obj@meta.data$unique_group_id);


# remove groups with fewer than 15 cells.
# model: assay + donor + cluster

# 


# edgeR 




# Cluster specific DE in each assay-type
source("~/scripts/LiverMap2.0/My_R_Scripts.R")

merged_obj@meta.data$unique_group_id <- paste(merged_obj@meta.data$assay_type, merged_obj@meta.data$donor, merged_obj@meta.data$seurat_clusters, sep="_")

group_means <- group_rowmeans(merged_obj@assays$RNA@data, merged_obj@meta.data$unique_group_id);


# remove groups with fewer than 15 cells.
# model: assay + donor + cluster

# 


# edgeR 

