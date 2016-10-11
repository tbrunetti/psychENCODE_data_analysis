# if using half-rack VM, enviromental proxy setting is required
Sys.setenv(http_proxy="http://cloud-proxy:3128")
Sys.setenv(https_proxy="http://cloud-proxy:3128")
library(data.table)
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")

total_raw_counts <- read.table("/home/tonya/Desktop/express_effective_counts_matix", header=TRUE, row.names=1, check.names=FALSE)
metadata <- read.table("synapse-meta-clinical-technical-data-BrainGVEX-RNAseq.csv", header=TRUE, row.names=1, check.names=FALSE)
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

# remove genes where all samples have zero reads
remove_zeros <- total_raw_read_counts_dataframe[rowSums(total_raw_read_counts_dataframe)!=0, ]

#transpose dataframe so genes names across top, samples names named in rows
reformatted_total_raw_reads_dataframe <- t(remove_zeros)
# prevent log(0), add a read count to every raw count
added_pseudocount <- reformatted_total_raw_reads_dataframe + 1
# log transform data
log_data <- log(added_pseudocount)

# calculate principal components
prin_comps <- prcomp(log_data, center=TRUE, scale.=TRUE)
# print out most important components
summary(prin_comps)