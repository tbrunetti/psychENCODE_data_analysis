library(ggfortify)
library(ggplot2)
library(data.table)
pdf(file="/data/users/tbrunetti/PCA_RNA_seq_data_BIDs_missing_data_removed.pdf")

# read in gene expression table and metadata

total_raw_counts <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/express_effective_counts_matix_rounded", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/synapse-meta-clinical-technical-data-BrainGVEX-RNAseq_variables_to_regress_numerical_and_nonNumerical.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

#remove sum zero row
total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe[rowSums(total_raw_read_counts_dataframe)!=0, ]
#add pseudo-count
total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe + 1
# log and transpose of dataframemetadata <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/synapse-meta-clinical-technical-data-BrainGVEX-RNAseq_numerically_converted_noKey.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")

tranformed_dataframe <- log(total_raw_read_counts_dataframe)

# remove BIDs that have missing or strange data
total_raw_read_counts_dataframe <- within(total_raw_read_counts_dataframe, rm('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'))
metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'), ]

# sort dataframes
counts_sorted <- total_raw_read_counts_dataframe[,order(colnames(total_raw_read_counts_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]

# transpose of expression dataframe
tranformed_dataframe <- t(tranformed_dataframe)
# make sure same order
all(rownames(metadata_sorted)==rownames(tranformed_dataframe))

pca_matrix <- prcomp(tranformed_dataframe, center=TRUE, scale. = TRUE)
#print(pca_matrix)
plot(pca_matrix, type ="l")

autoplot(pca_matrix, data=metadata_dataframe, colour = 'Sex')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'BrainBank')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'Hemisphere')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'TissueState')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'ERCC_Added')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'FlowcellBatch')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'SequencingPlatform')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'Multiplex_Oligo_Number')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'Final_Bead_Size_Selection_Ratio')
autoplot(pca_matrix, data= metadata_dataframe, colour = 'RNAIsolationBatch')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'LibraryBatch')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'Ethnicity')
autoplot(pca_matrix, data = metadata_dataframe, colour = 'YearAutopsy')
autoplot(pca_matrix, data=metadata_dataframe, colour = 'Diagnosis' )

dev.off()