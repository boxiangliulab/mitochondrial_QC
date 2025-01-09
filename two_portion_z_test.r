#####################
####generated by Chenyu Yang: chenyuy@mail.smu.edu ###
####Huang, Y., Wan, Z., Tang, Y. et al. Pantothenate kinase 2 interacts with PINK1 to regulate mitochondrial quality control via acetyl-CoA metabolism. Nat Commun 13, 2412 (2022). https://doi.org/10.1038/s41467-022-30178-x 
##### modified by Yinglu Tang: yinglut@smu.edu#####
#####2025-01-08####
# Set the working directory
setwd("D:/R_codes/colab_project")
rm(list = ls())

###########################################
# Part 1: Read in the dataset
###########################################
library(readxl)

# Define the path to the Excel workbook
filename<- "mito_data.xls"
my_data <- read_csv(filename)

# Assign names as 'Group_*' to the columns of features
my_list <- lapply(1:ncol(my_data), function(i) {
  col_data <- drop_na(my_data[[i]])
  assign(paste("Group_", LETTERS[i], sep = ""), as.vector(col_data), envir = .GlobalEnv)
  col_data
})

###########################################
# Part 2: Conduct Z-test
###########################################
Z_test <- c()
control <- my_list[[1]]

for (i in 2:length(my_list)) {
  group <- my_list[[i]]
  for (threshold in 2:10) { # Calculate the p-value for various thresholds
    p1 <- mean(control > threshold)
    p2 <- mean(group > threshold)
    all_data <- c(control, group)
    p0 <- mean(all_data > threshold)
    Z_stat <- (p1 - p2) / sqrt(p0 * (1 - p0) * (1 / length(control) + 1 / length(group)))
    Z_test <- c(Z_test, pnorm(abs(Z_stat), lower.tail = FALSE) * 2)
  }
}

###########################################
# Part 3: Save results to a CSV file
###########################################
# Reshape the p-values into a matrix
PVALUESZ <- matrix(Z_test, ncol = (ncol(my_data) - 1), byrow = FALSE)

# Save the matrix to a CSV file
output_path <- "Z_test_results.csv"
write.csv(PVALUESZ, file = output_path, row.names = FALSE)

# Inform the user of the file location
cat("Results saved to:", output_path, "\n")
