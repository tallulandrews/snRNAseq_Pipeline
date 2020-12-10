require("biomaRt")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
res <- getBM(c("ensembl_gene_id", "ensembl_transcript_id", "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end"), mart=ensembl, filter="with_hgnc", values=TRUE)
res <- res[!(is.na(res[,3]) & is.na(res[,4]) & is.na(res[,5]) & is.na(res[,6])),]
res$utr_5_len <- res[,4]-res[,3]
res$utr_3_len <- res[,6]-res[,5]
write.table(res, "UTRs_ensembl_hsap_7_July2020.txt", col.names=TRUE, row.names=FALSE)

