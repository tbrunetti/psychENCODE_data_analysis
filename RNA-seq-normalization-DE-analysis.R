# if using half-rack VM, enviromental proxy setting is required
#Sys.setenv(http_proxy="http://cloud-proxy:3128")
#Sys.setenv(https_proxy="http://cloud-proxy:3128")
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
library(limma)
library("BiocParallel")
# number refers to the number or threads/processors to run analysis on (parallelization)
register(MulticoreParam(7))


#-----------------gene expression analysis portion (DE)----------------------
total_raw_counts <- read.table("/home/tonya/Desktop/express_effective_counts_matix", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.table("synapse-meta-clinical-technical-data-BrainGVEX-RNAseq.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
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


deseq_obj <- DESeqDataSetFromMatrix(countData= counts_sorted, colData=metadata_sorted , design = ~ LibraryBatch + FlowcellBatch + BrainBank + UF_MEDIAN_5PRIME_BIAS + UF_MEDIAN_5PRIME_TO_3PRIME_BIAS + UF_MEDIAN_3PRIME_BIAS + RIN + Average_bp_size_of_Library + ERCC_Added + TissueState + Final_Bead_Size_Selection_Ratio + Diagnosis )
#set comparison reference/control
print("Releveling reference to Control")
deseq_obj$Diagnosis <- relevel(factor(deseq_obj$Diagnosis), ref="Control")

# pre-filtering, remove rows in count data with 0 or 1 reads
deseq_obj <- deseq_obj[rowSums(counts(deseq_obj)) > 1, ]

# perform DE analysis (size factors, dispersion, negative binomial distribution)
deseq_obj <- DESeq(deseq_obj)

de_results <- results(deseq_obj)
resOrdered <- de_results[order(de_results$padj),]
summary(results)
write.csv(as.data.frame(resOrdered), file="/data/users/tbrunetti/UC-UIC-covariates-added-to-model-testing-condition.csv")

#those adj p-value of <0.05
res05 <- results(deseq_obj, alpha=0.05)
summary(res05)
adj05 <- subset(resOrdered, padj < 0.05)
write.csv(as.data.frame(adj05), file="/data/users/tbrunetti/UC-UIC-covariates-added-to-model-testing-condition-less-pval-0.05.csv")


#those adj p-value of <0.01
res01 <- results(deseq_obj, alpha=0.01)
summary(res01)
adj01 <- subset(resOrdered, padj < 0.01)
write.csv(as.data.frame(adj01), file="/data/users/tbrunetti/UC-UIC-covariates-added-to-model-testing-condition-less-pval-0.01.csv")



# output normalized count data
norm_counts <- counts(deseq_obj, normalized=TRUE)
write.csv(norm_counts, file='/data/users/tbrunetti/UC-UIC-deseq2-normalized-counts-added-covariates-to-model.csv')
