source("http://bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("pvca")
library(Biobase)
library(pvca)
library(data.table)

png(filename="/data/users/tbrunetti/batch-effects-UIC-UC-RNA-seq.png")

# read in gene expression table and metadata
total_raw_counts <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/express_effective_counts_matix_rounded", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/synapse-meta-clinical-technical-data-BrainGVEX-RNAseq_variables_to_regress_numerical_and_nonNumerical.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

# remove samples now that you don't want to use for analysis
total_raw_read_counts_dataframe <- within(total_raw_read_counts_dataframe, rm('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'))
metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'), ]

# sorts dataframes so columns and rows are in same BID order for count data and metadata
counts_sorted <- total_raw_read_counts_dataframe[,order(colnames(total_raw_read_counts_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]

# below statement should result in TRUE, else the two dataframes are not properly sorted
all(rownames(metadata_sorted)==colnames(counts_sorted))

# convert back to matrix to create object ExpressionSet
mat_counts_sorted <- as.matrix(counts_sorted)
metadata_extra_annotations <- data.frame(labelDescription = seq(1, nrow(metadata_sorted)), row.names = colnames(metadata_sorted))
phenoData <- new("AnnotatedDataFrame", data = metadata_sorted, varMetadata = metadata_extra_annotations)

# creates experssionSet object properly formatted to be used by pvca R package
full_experiment_obj <- ExpressionSet(assayData = mat_counts_sorted, phenoData = phenoData)

# illustrate the batch effects that account for x percentage of the variability in the dataset
pct_threshold <- 0.6
# column names in metadata to be considered for batch effects
batch.factors <- colnames(metadata_sorted)
# calculate pvca
pvcaObj <- pvcaBatchAssess(full_experiment_obj, batch.factors = batch.factors, pct_threshold)

# plot the barplot of variances from batch effects
variance_plot <- barplot(pvcaObj$dat, xlab="Effects", ylab = "Weighted average proportion variance",
                         ylim = c(0, 1.1),col=c("blue"), las=2,
                         main="PVCA estimation bar chart")
axis(1, at=variance_plot, labels = pvcaObj$label, xlab= "Effects", cex.axis=0.5, las=2)
values = pvcaObj$dat
new_values = round(values, 3)
text(variance_plot, pvcaObj$dat, labels=new_values, pos=3, cex=0.8)

dev.off()

print(sessionInfo())

