##### NOTE ####
#
# This script is specifically designed for datasets where an excessive number of cells were 
# called using cellranger's pipeline. DO NOT USE on datasets where cell ranger identified 
# fewer than 20,000 cells!!!
#
##############

require(methods)
require(Matrix)
require("DropletUtils")
require(dplyr)
require(Seurat)
script_dir = "/cluster/home/tandrews/scripts/LiverMap2.0"

args <- as.character(commandArgs(trailingOnly=TRUE))
# FOLDER of raw data
# NAME of project
# max cells
# min cells
FOLDER <- args[1]
NAME <- args[2];
if (length(args) < 4) {
	MIN_CELLS <- 100
	if (length(args) < 3){
		MAX_CELLS <- 20000
	} else {
		MAX_CELLS <- args[3];
	}
} else {
	MAX_CELLS <- args[3];
	MIN_CELLS <- args[4];
}

FDR=0.01
TRIM=10 # change this for dynamic trimming to avoid memory overflow. - no don't

print(paste(FOLDER, NAME, MIN_CELLS, MAX_CELLS, FDR, TRIM, sep="\n"));

# Read the raw matrix
mydata <- Read10X(data.dir = FOLDER)
# Remove rows & columns that are completely zero
mydata <- mydata[,Matrix::colSums(mydata) > 0]
mydata <- mydata[Matrix::rowSums(mydata) > 0,]
print(paste("Stats out:", dim(mydata), "input",c("gene","cells"), "of", NAME))

# Rank droplet barcodes by total UMIs
br.out <- barcodeRanks(mydata)

# Knee and Inflection point plot
#plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
#o <- order(br.out$rank)
#lines(br.out$rank[o], br.out$fitted[o], col="red")

#abline(h=br.out$knee, col="dodgerblue", lty=2)
#abline(h=br.out$inflection, col="forestgreen", lty=2)
#legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#    legend=c("knee", "inflection"))

# Get total UMIs that correspond to the min & max thresholds.
plausible_cell_threshold <- max(br.out$total[br.out$rank > MAX_CELLS]);
mandatory_cell_threshold <- max(br.out$total[br.out$rank < MIN_CELLS]);

# My calling threshold
n_umi_sorted <- br.out$total[order(br.out$total, decreasing = T)]
rank_sorted <- 1:length(n_umi_sorted);

slope <- diff(log(n_umi_sorted))/diff(log(rank_sorted))
smooth_slope <- smooth.spline(slope, spar=0.5)
inflection <- which(smooth_slope$y == min(smooth_slope$y[MIN_CELLS:MAX_CELLS]))
my_inflection <- n_umi_sorted[inflection];

#TRIM <- min(plausible_cell_threshold, mandatory_cell_threshold, my_inflection, br.out$total[br.out$inflection], br.out$total[br.out$knee])/10;


if (br.out$knee > MIN_CELLS*10 || br.out$inflection > MIN_CELLS*10) {
        is_empty_threshold = max(plausible_cell_threshold,
                                br.out$total[br.out$knee],
                                br.out$total[br.out$inflection], my_inflection)
} else {
        is_empty_threshold = plausible_cell_threshold
}


set.seed(100)
# Run EmptyDrops
e.out <- emptyDrops(mydata, lower=is_empty_threshold, niters=100000, ignore=max(TRIM, is_empty_threshold/50), retain=mandatory_cell_threshold, test.ambient=TRUE)

# Clean up results
e.out <- e.out[!is.na(e.out$PValue),] # Remove NAs (too few UMIs to estimate PValue)
e.out$q.value <- e.out$PValue;
e.out$q.value[e.out$Limited] <- 0; # Set those with p < 1/n_simulations to p = 0 to ensure they are kept after multiple testing correction
e.out$q.value <- p.adjust(e.out$q.value, method="fdr"); # apply FDR multiple testing correction.

# Get number of detected genes/droplet
mydata <- mydata[,match(rownames(e.out),colnames(mydata))]
e.out$ngenes <- Matrix::colSums(mydata>0);

# Significance threshold
is.cell <- e.out$q.value <= FDR
sum(is.cell, na.rm=TRUE)

# Plot of selected genes
png(paste(NAME, "emptyDrops.png", sep="_"), width=6, height=6, units="in", res=150)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()

# Subset the matrix to the selected droplets and save it to a file.
outdata <- mydata[,is.cell]
print(paste("Stats out:", dim(outdata), "output",c("gene","cells"), "of", NAME))
#saveRDS(outdata, paste(NAME, "emptyDrops_table.rds", sep="_"));

# Estimate Background expression profile
umi.tot <- Matrix::colSums(mydata);
my_background <- Matrix::rowSums(mydata[, !is.cell])

# Estimate Background in cells.
# Get HVG, Calculate amount of background for each of those HVG
# Use lowest X percentile of those HVG as the actual estimate of background
# Subtract the background.

# Rational: HVGs are DE between cell-types
# If there exists genes specific to a cell-type they should be in HVGs 
# thus no cell-type should express all the HVGs.
# Thus the expression of the lowest relative expression among the HVGs is likely background.
# problem due to ton of zeros most HVGs won't be expressed at all.


my_background_correction <- function(sample_mat, is.cell) {
        #my_background <- Matrix::rowSums(sample_mat[,!is.cell]);
        zeros <- my_background==0;
        # don't want zeros in the background.
        # set them to 1s
        norm_factor <- (sum(my_background)+sum(zeros))/sum(my_background)
        my_background <- my_background*norm_factor;
        my_background[zeros] <- 1;
        my_background <- my_background/sum(my_background);

        correct_cell <- function(c) {
                # Ignore zeros since we assume all genes expressed in the background to some
                # extent, thus the background distribution is far more distributed than the cell
                # due to sparsity.
                # Assuming a cell is strictly made up of contamination + cell then this means
                # the zeros cannot tell us anything about the contamination.
                # Expectation if all background
                #expected <- my_background;
                #expected[c==0] <- 0;
                #expected <- expected/sum(expected)*sum(c);

                # if not much background theen there should be many big difference
                # however there will also be differences due to chance, and the differences
                # due to chance are proportional to the total.

                # But we also expect those with many counts to have more non-background.
                # Because of the scaling always expect about half to be -vs and half to be +vs
                # if same as background.
                #diff <- c-expected;

                # top most expressed genes in the background should be the guide to fitting the
                # amount.
                expected <- my_background;
                expected[c==0] <- 0;
                top_genes <- quantile(expected[expected>0], 0.80)

                top_bg_genes <- expected[expected >= top_genes]
                c_bg_genes <- c[names(c) %in% names(top_bg_genes)]
                bg_factor <-sum(c_bg_genes)/sum(top_bg_genes)
                expected <- my_background*bg_factor;

                diff <- c-expected;
                diff[diff < 0] <- 0
                return(diff)
        }
        corrected_mat <- apply(sample_mat, 2, correct_cell);
        return(corrected_mat);
}

print(paste("My correction matrix dimensions:", dim(mydata), "input",c("gene","cells"), "of", NAME))

saveRDS(mydata[,is.cell], paste(NAME, "emptyDrops_table.rds", sep="_"))
saveRDS(my_background_correction(mydata[,is.cell], is.cell), paste(NAME, "emptyDrops_correctedtable.rds", sep="_"));
saveRDS(my_background, paste(NAME, "background.rds", sep="_"))
