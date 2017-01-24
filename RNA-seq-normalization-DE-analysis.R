# if using half-rack VM, enviromental proxy setting is required
#Sys.setenv(http_proxy="http://cloud-proxy:3128")
#Sys.setenv(https_proxy="http://cloud-proxy:3128")
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
library(limma)
library("BiocParallel")
register(MulticoreParam(10))

#-----------------gene expression analysis portion (DE)----------------------
total_raw_counts <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/express_effective_counts_matrix_rounded_gene_level.csv", header=TRUE, row.names=432, check.names=FALSE, sep=",")
total_raw_counts <- total_raw_counts[ , -which(names(total_raw_counts) %in% c("nothing"))]
metadata <- read.table("/data/users/tbrunetti/pyschENCODE-DE-analysis/synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-MISSING-DATA-FILLED-WITH-AVERAGES.csv", header=TRUE, row.names=1, check.names=FALSE, sep=",")
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

# remove samples now that you don't want to use for analysis
#total_raw_read_counts_dataframe <- within(total_raw_read_counts_dataframe, rm('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'))
#metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'), ]

# sorts dataframes so columns and rows are in same BID order for count data and metadata
counts_sorted <- total_raw_read_counts_dataframe[,order(colnames(total_raw_read_counts_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]

# below statement should result in TRUE, else the two dataframes are not properly sorted
all(rownames(metadata_sorted)==colnames(counts_sorted))


deseq_obj <- DESeqDataSetFromMatrix(countData= counts_sorted, colData=metadata_sorted , design= ~UF_MEDIAN_5PRIME_TO_3PRIME_BIAS + LibraryBatch+ PMI + Sex + BrainBank + RIN +Diagnosis )
#set comparison reference/control
print("Releveling reference to Control")
deseq_obj$Diagnosis <- relevel(factor(deseq_obj$Diagnosis), ref="Control")

# pre-filtering, remove rows in count data with 0 or 1 reads
deseq_obj <- deseq_obj[rowSums(counts(deseq_obj)) > 1, ]

# perform DE analysis (size factors, dispersion, negative binomial distribution)
deseq_obj <- DESeq(deseq_obj, parallel=TRUE)

de_results_control_bp <- results(deseq_obj, parallel=TRUE, contrast=c("Diagnosis", "BP", "Control"))
de_results_control_scz <- results(deseq_obj, parallel=TRUE, contrast=c("Diagnosis", "SCZ", "Control"))
de_results_bp_scz <- results(deseq_obj, parallel=TRUE, contrast=c("Diagnosis", "BP", "SCZ"))

resOrdered_c_bp <- de_results_control_bp[order(de_results_control_bp$padj),]
summary(de_results_control_bp)
#write.csv(as.data.frame(resOrdered_c_bp), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-control_vs_bp.csv")


resOrdered_c_scz <- de_results_control_scz[order(de_results_control_scz$padj),]
summary(de_results_control_scz)
#write.csv(as.data.frame(resOrdered_c_scz), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-control_vs_scz.csv")


resOrdered_bp_scz <- de_results_bp_scz[order(de_results_bp_scz$padj),]
summary(de_results_bp_scz)
#write.csv(as.data.frame(resOrdered_bp_scz), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-bp_vs_scz.csv")

#those adj p-value of <0.05 control vs BP
res05_c_bp <- results(deseq_obj, alpha=0.05, contrast=c("Diagnosis", "BP", "Control"))
summary(res05_c_bp)
adj05_c_bp <- subset(resOrdered_c_bp, padj < 0.05)
#write.csv(as.data.frame(adj05_c_bp), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-control_vs_bp_less-pval-0.05.csv")


#those adj p-value of <0.05 control vs scz
res05_c_scz <- results(deseq_obj, alpha=0.05, contrast=c("Diagnosis", "SCZ", "Control"))
summary(res05_c_scz)
adj05_c_scz <- subset(resOrdered_c_scz, padj < 0.05)
#write.csv(as.data.frame(adj05_c_scz), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-control_vs_scz_less-pval-0.05.csv")

#those adj p-value of <0.05 BP vs SCZ
res05_bp_scz <- results(deseq_obj, alpha=0.05, contrast=c("Diagnosis", "BP", "SCZ"))
summary(res05_bp_scz)
adj05_bp_scz <- subset(resOrdered_bp_scz, padj < 0.05)
#write.csv(as.data.frame(adj05_bp_scz), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-bp_vs_scz_less-pval-0.05.csv")


#those adj p-value of <0.01 control vs BP
res01_c_bp <- results(deseq_obj, alpha=0.01, contrast=c("Diagnosis", "BP", "Control"))
summary(res01_c_bp)
adj01_c_bp <- subset(resOrdered_c_bp, padj < 0.01)
#write.csv(as.data.frame(adj01_c_bp), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-control_vs_bp_less-pval-0.01.csv")


#those adj p-value of <0.01 control vs scz
res01_c_scz <- results(deseq_obj, alpha=0.01, contrast=c("Diagnosis", "SCZ", "Control"))
summary(res01_c_scz)
adj05_c_scz <- subset(resOrdered_c_scz, padj < 0.01)
#write.csv(as.data.frame(adj05_c_scz), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-control_vs_scz_less-pval-0.01.csv")

#those adj p-value of <0.01 BP vs SCZ
res01_bp_scz <- results(deseq_obj, alpha=0.01, contrast=c("Diagnosis", "BP", "SCZ"))
summary(res01_bp_scz)
adj01_bp_scz <- subset(resOrdered_bp_scz, padj < 0.01)
#write.csv(as.data.frame(adj01_bp_scz), file="/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-covariates-added-to-model-testing-condition-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-bp_vs_scz_less-pval-0.01.csv")


# output normalized count data
norm_counts <- counts(deseq_obj, normalized=TRUE)
write.csv(norm_counts, file='/data/users/tbrunetti/pyschENCODE-DE-analysis/results/UC-UIC-deseq2-normalized-counts-all-430-samples-PMI-sex-brainbank-53bias-libraryBatch-RIN-Diagnosis-only-in-model.csv')
