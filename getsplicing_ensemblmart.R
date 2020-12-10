library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
biomaRt::listAttributes(ensembl)
out <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'),
      mart = ensembl)

out <- getBM(attributes = c('ensembl_gene_id', 
		'ensembl_transcript_id', 
		'transcript_exon_intron', "percentage_gene_gc_content")
	, mart = ensembl)