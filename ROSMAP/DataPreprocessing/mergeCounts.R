mergeCounts <- function(file1 = "Data/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv", 
                        file2 = "Data/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv", 
                        num_to_remove_beginning = 0L, num_to_remove_end = 2L,
                        num_col_with_name = c(1, 2)){
  # ---- Description ----
  # It merges a data set divided into two files and give format to the final data set 
  # removing some characters at the beginning and end of the column names, setting 
  # the row names, and deleting the columns with data used to the row names 
  #
  # ---- Parameters ----
  # file1: character 
  #     Path to the file 1
  # file 2: character
  #     Path to the file 2
  # num_to_remove_beginning: integer
  #     Number of characters to remove at the beginning of the column names 
  # num_to_remove_end: integer
  #     Number of characters to remove at the beginning of the column names 
  # num_col_with_name: vector
  #     Vector with indices of the columns with the row names 
  # 
  # ---- Returns ----
  # fullData: data.frame
  #     A complete data frame with the column names modified
  library(data.table)
  # read the files
  data1 <- fread(file = file1, sep = '\t', header = TRUE) 
  data2 <- fread(file = file2, sep = '\t', header = TRUE) 
  
  
  # merge the files 
  fullData <- as.data.frame(cbind(data1, data2[, !colnames(data2) %in% colnames(data1), with = FALSE]))  
  
  # Remove characters at the beginning and end 
  colnamesToChange <- (colnames(fullData) %in% c(colnames(fullData)[1],colnames(fullData)[2]) )==FALSE
  
  colnames(fullData)[  colnamesToChange ] <- substr(colnames(fullData)[  colnamesToChange ], 
                                                    num_to_remove_beginning + 1, 
                                                    nchar(colnames(fullData)[  colnamesToChange ])-num_to_remove_end)
  
  # set the row names of fullData 
  rownames(fullData) <- fullData[, num_col_with_name[1]]
  fullData <- fullData[, -num_col_with_name]
  
  fullData
}