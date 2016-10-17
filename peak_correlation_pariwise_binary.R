library(corrplot)
example(pairs)

imported_matrix <- read.table("/home/tonya/dumb_set.csv", row.names = 1, sep=',')
transposed_data <- as.data.frame(t(imported_matrix))

# gets number of samples to create corr-matrix data-frame
total_samples <- ncol(transposed_data)

# create custom corr_matrix and changes name of row and columns to sample names (based on transposed_data)
custom_corr_matrix <- as.data.frame(matrix(0, nrow=total_samples, ncol = total_samples), row.names = colnames(transposed_data), )
colnames(custom_corr_matrix) <- colnames(transposed_data)

for (i in 1:nrow(custom_corr_matrix)){
  for (j in 1:ncol(custom_corr_matrix)){
    message(sprintf("Progress is at row %s and col %s out of the %s by %s matrix", i, j, nrow(custom_corr_matrix), ncol(custom_corr_matrix)))
    corr_list_first <- list()
    corr_list_second <- list()
    # extracts column and row names
    first_sample <- rownames(custom_corr_matrix)[i]
    second_sample <- colnames(custom_corr_matrix)[j]
    # if both are equal to 0, then remove the peak, else add peak to list   
    for (x in 1:nrow(transposed_data)){
      if ((transposed_data[first_sample][x,]==0) & (transposed_data[second_sample][x,]==0)){
        next;
      }
      else{
        corr_list_first <- c(corr_list_first, unlist(transposed_data[first_sample][x,]))
        corr_list_second <- c(corr_list_second, unlist(transposed_data[second_sample][x,]))
      }
    }
    
    # calculate the pairwise correlation without the both zero regions and populate custom corrlation matrix
    pairwise_corr <- cor(unlist(corr_list_first), unlist(corr_list_second), method="spearman")
    custom_corr_matrix[i,j]<-pairwise_corr
  }
}

# take absolute value to convert all to positive values (removing double zeros gets all neg values, no std-dev??)
custom_corr_matrix[is.na(custom_corr_matrix)] <- 1
corrplot(abs(as.matrix(custom_corr_matrix)), method="color", type="upper", order = "hclust")

