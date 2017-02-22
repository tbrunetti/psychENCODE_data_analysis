#prcomp is used instread of princomp here to perform PCA since:
# 1) princomp is limited to experiments where observations >> variables
# 2) prcomp can be used to easily visualize clusters using PCs and metadata

# ggfortify and ggplot are required for plotting and clustering PCs to metadatas
library(ggfortify)
library(ggplot2)
library(data.table)
library(gplots)

# output_file is the name of the PDF the user would like to create and write to
# the pdf() function actually creates the pdf
output_file = "/data/users/tbrunetti/psychENCODE-differential-new/test-UC-UIC-RNAseq-BP-SCZ-noControl-raw-gene-counts-data-only.pdf"
pdf(file = output_file)

# counts_file is the FULL path to the counts matrix in CSV format--MUST HAVE HEADER AND ROW NAMES!!!
#             counts must be rounded to nearest integer
#             row names should be gene/transcript names
#             headers should be BIDs or sample ID names
# meta_file is the FULL path to the metadata matrix in CSV format--MUST HAVE HEADER AND ROW NAMES!!!
#             row names should be BIDs or sample ID names
#             headers should be the name of the metadata measurement (i.e. sex, diagnosis, chrM contamination, etc...)
# to keep it straight, it makes sense to have n x m counts files and a m x u because the matrices can be multiplied since
# their inner components are the same dimensions ( nxm matrix [dotproduct] mxu matrix = nxu matrix)
counts_file = "/data/users/tbrunetti/psychENCODE-differential-new/express_effective_counts_matrix_rounded_gene_level.csv"
meta_file = "/data/users/tbrunetti/psychENCODE-differential-new/synapse-meta-clinical-technical-data-BrainGVEX-RNAseq-TonyaCurrated-V2-MISSING-DATA-FILLED-WITH-AVERAGES.csv"

# read in gene expression table and metadata row.names needs to be set to the column number where the row names are listed
# it is important to set check.names=FALSE so that R doesn't try to change any numerical names into odd characters
total_raw_counts <- read.table(counts_file, header=TRUE, row.names=432, check.names=FALSE, sep=",")
metadata <- read.table(meta_file, header=TRUE, row.names=1, check.names=FALSE, sep=",")

# converts each table into a data frame
total_raw_read_counts_dataframe <- as.data.frame(total_raw_counts)
metadata_dataframe <- as.data.frame(metadata)

#pre-filter step, remove any ROWS that have zero counts in the count matrix, it does not mean that a sample cannot
# have zero counts, we are just removing the genes where no counts exist across all samples
total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe[rowSums(total_raw_read_counts_dataframe)!=0, ]

#add pseudo-count, why??  Because the log(0) = Does not exist, so to address this issue add 1 to all counts
total_raw_read_counts_dataframe <- total_raw_read_counts_dataframe + 1
# take the log of the counts data, this will help normalize the data
transformed_dataframe <- log(total_raw_read_counts_dataframe)

# remove BIDs that have missing or strange data
transformed_dataframe <- within(transformed_dataframe, rm('nothing'))
metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% c('nothing'), ]

# The remaining block of code can be uncommented out if you are going to use a SUBSET of the data for PCA
# in the metadata_to_be_removed variable you can indicate in which() statment the column with the specified variable
# of all the samples to remove for analysis.  For example, in this case we go to the column 'Diagnosis' and remove
# all samples that are listed as 'Control'.  This can be changed to anything.
#metadata_to_be_removed <- metadata_dataframe[which(metadata_dataframe$Diagnosis=='Control'), ]
#sample_removal <- levels(droplevels(metadata_to_be_removed$BID))
#transformed_dataframe <-transformed_dataframe[,-which(names(transformed_dataframe) %in% sample_removal)]
#metadata_dataframe <- metadata_dataframe[!rownames(metadata_dataframe) %in% sample_removal, ]
#print("Samples to be removed from counts matrix and metadata")
#print(sample_removal)


#remove rows that have a sum of zero
transformed_dataframe <- transformed_dataframe[rowSums(transformed_dataframe)!=0, ]

# sort dataframes.  Dataframes MUST be sorted by colnames in the counts matrix
# and sorted by the rownames in the metadata matrix.  This ensures that the sample names properly match
# each other between counts matix and metadata matrix.  Note, you will see the word 'TRUE' printed
# if they are properly match, else you will see 'FALSE' in which case you need to sort
counts_sorted <- transformed_dataframe[,order(colnames(transformed_dataframe))]
metadata_sorted <- metadata_dataframe[order(rownames(metadata_dataframe)),]
all(rownames(metadata_sorted)==colnames(counts_sorted))

# This part is more important if you are removing data from analysis.  It has to extract and collapse missing
# levels from the dataframe of data that was removed
print("getting levels")
sapply(metadata_dataframe,levels)

# Calculation of the prinipal components using prcomp
# the t() function means to take the transpose of the counts matrix
# scale.=TRUE means PCs will be baed on corrlation matrix and not covariance matrix
# corr is better if scales of variables are very different
# center=TRUE uses the centers the data, use centroid/mean
pca_matrix <- prcomp(t(counts_sorted), center=TRUE, scale. = TRUE)

# plot the PCs againts how much variance each PC contributes to the data set
# looking to incorporate the min number of PCs before "elbowing-effect" into the model unless
# a PC is strongly correlated with a variable in the metadata set, in which case,
# just regress out the variable rather than the PC
plot(pca_matrix, type ="l")


# ***************************Function to correlate PCs with metadata******************************
model_pcs <- function(pca_matrix){
  # this block of code will need to be changed depending on the metadata columns available 
  factors_affecting_pcs=list()
  # iterate through all PCs and get -log10(p-value) and adjusted R-squared value to every PC correlated to each
  # metadata column listed
  for (pc in seq(1,dim(pca_matrix$rotation)[2])){
    # correlation PCs on factor sex
    # lm() is the Limma library function to build a linear model, in this case each PC is being tested against 'Sex'
    # to determine if any PC is strongly correlated to Sex.  If yes, regress sex out in model
    # Note: in the linear model as.factor() should be placed around categorical variables, while it can be omitted
    # in continuous non-discrete variables
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Sex'])))
    factors_affecting_pcs[['Sex']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Sex']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Sex']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor BrainBank
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'BrainBank'])))
    factors_affecting_pcs[['BrainBank']][[as.character(pc)]]=list()
    factors_affecting_pcs[['BrainBank']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['BrainBank']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Hemisphere
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Hemisphere'])))
    factors_affecting_pcs[['Hemisphere']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Hemisphere']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Hemisphere']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor LibraryBatch
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'LibraryBatch'])))
    factors_affecting_pcs[['LibraryBatch']][[as.character(pc)]]=list()
    factors_affecting_pcs[['LibraryBatch']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['LibraryBatch']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_MEDIAN_5PRIME_TO_3PRIME_BIAS
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_MEDIAN_5PRIME_TO_3PRIME_BIAS']))
    factors_affecting_pcs[['UF_MEDIAN_5PRIME_TO_3PRIME_BIAS']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_MEDIAN_5PRIME_TO_3PRIME_BIAS']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_MEDIAN_5PRIME_TO_3PRIME_BIAS']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Ethnicity
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Ethnicity'])))
    factors_affecting_pcs[['Ethnicity']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Ethnicity']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Ethnicity']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    
    # correlation PCs on factor volume_RNA_starting_material_used_uL
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'volume_RNA_starting_material_used_uL']))
    factors_affecting_pcs[['volume_RNA_starting_material_used_uL']][[as.character(pc)]]=list()
    factors_affecting_pcs[['volume_RNA_starting_material_used_uL']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['volume_RNA_starting_material_used_uL']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Multiplex_Oligo_Number
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Multiplex_Oligo_Number'])))
    factors_affecting_pcs[['Multiplex_Oligo_Number']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Multiplex_Oligo_Number']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Multiplex_Oligo_Number']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Final_Bead_Size_Selection_Ratio
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Final_Bead_Size_Selection_Ratio'])))
    factors_affecting_pcs[['Final_Bead_Size_Selection_Ratio']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Final_Bead_Size_Selection_Ratio']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Final_Bead_Size_Selection_Ratio']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Final_Library_Bioanalyzer_Concentration_ng_per_uL
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Final_Library_Bioanalyzer_Concentration_ng_per_uL']))
    factors_affecting_pcs[['Final_Library_Bioanalyzer_Concentration_ng_per_uL']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Final_Library_Bioanalyzer_Concentration_ng_per_uL']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Final_Library_Bioanalyzer_Concentration_ng_per_uL']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor RIN
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'RIN']))
    factors_affecting_pcs[['RIN']][[as.character(pc)]]=list()
    factors_affecting_pcs[['RIN']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['RIN']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor PMI
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'PMI']))
    factors_affecting_pcs[['PMI']][[as.character(pc)]]=list()
    factors_affecting_pcs[['PMI']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['PMI']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor BrainWeight
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'BrainWeight']))
    factors_affecting_pcs[['BrainWeight']][[as.character(pc)]]=list()
    factors_affecting_pcs[['BrainWeight']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['BrainWeight']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor TIN_median
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'TIN_median']))
    factors_affecting_pcs[['TIN_median']][[as.character(pc)]]=list()
    factors_affecting_pcs[['TIN_median']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['TIN_median']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Average_bp_size_of_Library
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Average_bp_size_of_Library']))
    factors_affecting_pcs[['Average_bp_size_of_Library']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Average_bp_size_of_Library']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Average_bp_size_of_Library']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_MEDIAN_3PRIME_BIAS
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_MEDIAN_3PRIME_BIAS']))
    factors_affecting_pcs[['UF_MEDIAN_3PRIME_BIAS']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_MEDIAN_3PRIME_BIAS']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_MEDIAN_3PRIME_BIAS']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_MEDIAN_5PRIME_BIAS
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_MEDIAN_5PRIME_BIAS']))
    factors_affecting_pcs[['UF_MEDIAN_5PRIME_BIAS']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_MEDIAN_5PRIME_BIAS']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_MEDIAN_5PRIME_BIAS']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Total_Yield_ng
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'Total_Yield_ng']))
    factors_affecting_pcs[['Total_Yield_ng']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Total_Yield_ng']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Total_Yield_ng']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor SequencingPlatform
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'SequencingPlatform'])))
    factors_affecting_pcs[['SequencingPlatform']][[as.character(pc)]]=list()
    factors_affecting_pcs[['SequencingPlatform']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['SequencingPlatform']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    
    # correlation PCs on factor FlowcellBatch
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'FlowcellBatch'])))
    factors_affecting_pcs[['FlowcellBatch']][[as.character(pc)]]=list()
    factors_affecting_pcs[['FlowcellBatch']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['FlowcellBatch']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor TotalReads
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'TotalReads']))
    factors_affecting_pcs[['TotalReads']][[as.character(pc)]]=list()
    factors_affecting_pcs[['TotalReads']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['TotalReads']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor MappedReads_Primary
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'MappedReads_Primary']))
    factors_affecting_pcs[['MappedReads_Primary']][[as.character(pc)]]=list()
    factors_affecting_pcs[['MappedReads_Primary']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['MappedReads_Primary']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor MappedReads_Multimapped
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'MappedReads_Multimapped']))
    factors_affecting_pcs[['MappedReads_Multimapped']][[as.character(pc)]]=list()
    factors_affecting_pcs[['MappedReads_Multimapped']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['MappedReads_Multimapped']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_MEDIAN_CV_COVERAGE
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_MEDIAN_CV_COVERAGE']))
    factors_affecting_pcs[['UF_MEDIAN_CV_COVERAGE']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_MEDIAN_CV_COVERAGE']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_MEDIAN_CV_COVERAGE']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_PCT_RIBOSOMAL_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_RIBOSOMAL_BASES']))
    factors_affecting_pcs[['UF_PCT_RIBOSOMAL_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_RIBOSOMAL_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_RIBOSOMAL_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_PCT_CODING_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_CODING_BASES']))
    factors_affecting_pcs[['UF_PCT_CODING_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_CODING_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_CODING_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_PCT_UTR_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_UTR_BASES']))
    factors_affecting_pcs[['UF_PCT_UTR_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_UTR_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_UTR_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_PCT_INTRONIC_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_INTRONIC_BASES']))
    factors_affecting_pcs[['UF_PCT_INTRONIC_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_INTRONIC_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_INTRONIC_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_PCT_INTERGENIC_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_INTERGENIC_BASES']))
    factors_affecting_pcs[['UF_PCT_INTERGENIC_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_INTERGENIC_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_INTERGENIC_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor UF_PCT_MRNA_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_MRNA_BASES']))
    factors_affecting_pcs[['UF_PCT_MRNA_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_MRNA_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_MRNA_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    
    # correlation PCs on factor UF_PCT_USABLE_BASES
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_PCT_USABLE_BASES']))
    factors_affecting_pcs[['UF_PCT_USABLE_BASES']][[as.character(pc)]]=list()
    factors_affecting_pcs[['UF_PCT_USABLE_BASES']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['UF_PCT_USABLE_BASES']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor AgeDeath
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'AgeDeath']))
    factors_affecting_pcs[['AgeDeath']][[as.character(pc)]]=list()
    factors_affecting_pcs[['AgeDeath']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['AgeDeath']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
    
    # correlation PCs on factor Diagnosis
    linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Diagnosis'])))
    factors_affecting_pcs[['Diagnosis']][[as.character(pc)]]=list()
    factors_affecting_pcs[['Diagnosis']][[as.character(pc)]][['adj.r.squared']]=summary(linear_model)$adj.r.squared
    factors_affecting_pcs[['Diagnosis']][[as.character(pc)]][['-log10Pval']]=-log10(anova(linear_model)$Pr[1]) # -log10(pval)
  }
  
  # create heatmeap to visualize PC and metadata correlations
  # create a dataframe to store all the -log10 p-values and adjusted R squared vals for the visualization of PCA data
  pvalues <- data.frame()
  adjRsq <- data.frame()
  
  # iterate only through the look first 10 PCs and extract their -log10(pval) and adj R-sq values
  for (all_factors in seq(1,length(factors_affecting_pcs))){
    for (pc in seq(1,10)){
      pvalues[all_factors, pc] <- unlist(factors_affecting_pcs[all_factors][[1]][[pc]][[2]])
      adjRsq[all_factors, pc] <- unlist(factors_affecting_pcs[all_factors][[1]][[pc]][[1]])
    }
  }
  
  # get the row and column names match the p-values
  rownames(pvalues) <- names(factors_affecting_pcs)
  colnames(pvalues) <- unlist(lapply(seq(1,10),function(x) paste(c('PC',x),collapse='')))
  
  # get the row and column names matching the adj R-sq values
  rownames(adjRsq) <- names(factors_affecting_pcs)
  colnames(adjRsq) <- unlist(lapply(seq(1,10),function(x) paste(c('PC',x),collapse='')))
  
  # round all -log10(pvalue) in the dataframe to three decimal places
  is.num <- sapply(pvalues, is.numeric)
  pvalues[is.num] <- lapply(pvalues[is.num], round, 3)
  
  # create a heatmap of these values, value is -log10(p-val) and color is the adj R-sq value
  heatmap.2(as.matrix(adjRsq), cellnote=pvalues, notecol = "black", notecex = 0.5, cexRow = 0.3, dendrogram = "none", col=colorRampPalette(c("white", "yellow", "red"))(10))
  print("heatmap completed")
}
#*********************************************End of function*********************************************

model_pcs(pca_matrix)

#---------------------------regress out variables of interest-------------------------------------------------
# for each variable that is to be regressed a linear model must be made and the residuals of the linear model
# must be extract and replace the old PC matrix
# regress out FlowcellBatch
for (pc in seq(1,dim(pca_matrix$rotation)[2])){
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'FlowcellBatch'])))
  pca_matrix$x[,pc]  <- linear_model$residuals
}
# regress out Sex
for (pc in seq(1,dim(pca_matrix$rotation)[2])){
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(as.factor(metadata_sorted[,'Sex'])))
  pca_matrix$x[,pc]  <- linear_model$residuals
}
# regress out 5 prime to 3 prime bias
for (pc in seq(1,dim(pca_matrix$rotation)[2])){
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'UF_MEDIAN_5PRIME_TO_3PRIME_BIAS']))
  pca_matrix$x[,pc]  <- linear_model$residuals
}
# regress out RIN
for (pc in seq(1,dim(pca_matrix$rotation)[2])){
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'RIN']))
  pca_matrix$x[,pc]  <- linear_model$residuals
}
# regress out PMI
for (pc in seq(1,dim(pca_matrix$rotation)[2])){
  linear_model <- lm(pca_matrix$x[,pc] ~ na.omit(metadata_sorted[,'PMI']))
  pca_matrix$x[,pc]  <- linear_model$residuals
}

# call function again on regressed out variables
model_pcs(pca_matrix)
