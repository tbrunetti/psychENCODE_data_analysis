library(corrplot)
library(ggplot2)
library(data.table)

geneExpMatrix <- read.table(file = "/home/tonya/psych_ENCODE/First-pass-ATAC-seq-results/first89marked.jaccard.csv", header = TRUE, row.names = 1, check.names = FALSE)
sample_ids <- colnames(geneExpMatrix)
pairwise_corrs <-as.data.frame(matrix(0, ncol(geneExpMatrix), ncol(geneExpMatrix)))
colnames(pairwise_corrs) <- sample_ids
rownames(pairwise_corrs) <- sample_ids

for (row in seq(1, length(sample_ids))){
  print(row)
  for (col in seq(1, length(sample_ids))){
    pairwise_corrs[sample_ids[row], sample_ids[col]] <- cor(geneExpMatrix[,sample_ids[row]], geneExpMatrix[,sample_ids[col]], method = "pearson")
    pairwise_corrs[sample_ids[col], sample_ids[row]] <- pairwise_corrs[sample_ids[row], sample_ids[col]]
  }
}

corrplot(as.matrix(pairwise_corrs), method="color", type="upper", order = "hclust", tl.col = "black", cl.lim=c(-1 ,1), tl.cex = 0.5)
