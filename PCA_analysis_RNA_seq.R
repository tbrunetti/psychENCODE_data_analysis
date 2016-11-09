#prcomp is used instread of princomp here to perform PCA since
# 1) princomp is limited to exps where observations > variables
# 2) prcomp can be used to easily visualize clusters using PCs and metadata

library(ggfortify)
library(ggplot2)
library(data.table)
pdf(file="/data/users/tbrunetti/testing_new_PCA.pdf")

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
transformed_dataframe <- log(total_raw_read_counts_dataframe)

# remove BIDs that have missing or strange data
transformed_dataframe <- within(transformed_dataframe, rm('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'))
metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('2015-1', '2015-739', '2015-2883', '2016-1531', '2016-1532'), ]

#remove sum zero row
transformed_dataframe <- transformed_dataframe[rowSums(transformed_dataframe)!=0, ]

# sort dataframes
counts_sorted <- transformed_dataframe[,order(colnames(transformed_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]

# make sure same order
all(rownames(metadata_sorted)==colnames(counts_sorted))


# scale.=TRUE means PCs will be baed on corrlation matrix and not covariance matrix
# corr is better if scales of variables are very different
pca_matrix <- prcomp(t(counts_sorted), center=TRUE, scale. = TRUE)

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

#linear_model <- lm(pca_matrix$x ~ metadata_sorted[,'Average_bp_size_of_Library'])
factors_on_pcs=list()
for (pc in seq(1,dim(pca_matrix$rotation)[2])){
  #print(pc)
  # correlation PCs on factor sex
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Sex']))
  factors_on_pcs[['Sex']][[as.character(pc)]]=list()
  factors_on_pcs[['Sex']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['Sex']][[as.character(pc)]][['-log10Pval']]=-log10(summary(linear_model)$coefficients[,4]) # -log10(pval)

  # correlation PCs on factor BrainBank
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'BrainBank']))
  factors_on_pcs[['BrainBank']][[as.character(pc)]]=list()
  factors_on_pcs[['BrainBank']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['BrainBank']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor TissueState
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'TissueState']))
  factors_on_pcs[['TissueState']][[as.character(pc)]]=list()
  factors_on_pcs[['TissueState']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['TissueState']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor ERCC_Added
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'ERCC_Added']))
  factors_on_pcs[['ERCC_Added']][[as.character(pc)]]=list()
  factors_on_pcs[['ERCC_Added']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['ERCC_Added']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor Hemisphere
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Hemisphere']))
  factors_on_pcs[['Hemisphere']][[as.character(pc)]]=list()
  factors_on_pcs[['Hemisphere']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['Hemisphere']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor LibraryBatch
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'LibraryBatch']))
  factors_on_pcs[['LibraryBatch']][[as.character(pc)]]=list()
  factors_on_pcs[['LibraryBatch']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['LibraryBatch']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor MEDIAN_5PRIME_TO_3PRIME_BIAS
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'MEDIAN_5PRIME_TO_3PRIME_BIAS']))
  factors_on_pcs[['MEDIAN_5PRIME_TO_3PRIME_BIAS']][[as.character(pc)]]=list()
  factors_on_pcs[['MEDIAN_5PRIME_TO_3PRIME_BIAS']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['MEDIAN_5PRIME_TO_3PRIME_BIAS']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor volume_RNA_starting_material_used_uL
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'volume_RNA_starting_material_used_uL']))
  factors_on_pcs[['volume_RNA_starting_material_used_uL']][[as.character(pc)]]=list()
  factors_on_pcs[['volume_RNA_starting_material_used_uL']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['volume_RNA_starting_material_used_uL']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor MEDIAN_CV_COVERAGE
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'MEDIAN_CV_COVERAGE']))
  factors_on_pcs[['MEDIAN_CV_COVERAGE']][[as.character(pc)]]=list()
  factors_on_pcs[['MEDIAN_CV_COVERAGE']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['MEDIAN_CV_COVERAGE']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor Multiplex_Oligo_Number
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Multiplex_Oligo_Number']))
  factors_on_pcs[['Multiplex_Oligo_Number']][[as.character(pc)]]=list()
  factors_on_pcs[['Multiplex_Oligo_Number']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['Multiplex_Oligo_Number']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor Final_Bead_Size_Selection_Ratio
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Final_Bead_Size_Selection_Ratio']))
  factors_on_pcs[['Final_Bead_Size_Selection_Ratio']][[as.character(pc)]]=list()
  factors_on_pcs[['Final_Bead_Size_Selection_Ratio']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['Final_Bead_Size_Selection_Ratio']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor Final_Library_Bioanalyzer_Concentration_ng_per_uL
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Final_Library_Bioanalyzer_Concentration_ng_per_uL']))
  factors_on_pcs[['Final_Library_Bioanalyzer_Concentration_ng_per_uL']][[as.character(pc)]]=list()
  factors_on_pcs[['Final_Library_Bioanalyzer_Concentration_ng_per_uL']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['Final_Library_Bioanalyzer_Concentration_ng_per_uL']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  
  # correlation PCs on factor Diagnosis
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Diagnosis']))
  factors_on_pcs[['Diagnosis']][[as.character(pc)]]=list()
  factors_on_pcs[['Diagnosis']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
  factors_on_pcs[['Diagnosis']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
}

pvalues <- data.frame()
adjRsq <- data.frame()

# look at first 10 PCs
for (total_factors in seq(1,length(factors_on_pcs))){
  for (pc in seq(1,10)){
    pvalues[total_factors, pc] <- unlist(factors_on_pcs[total_factors][[1]][[pc]][2])
  }
}
colnames(pvalues) <- names(factors_on_pcs)
rownames(pvalues) <- unlist(lapply(seq(1,10),function(x) paste(c('PC',x),collapse='')))

dev.off()


