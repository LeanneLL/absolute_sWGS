library(rascal)
library(tidyverse)

input_dir <- '/path/to/_output.csv/'

#input_dir <- '/scratchb/stlab-icgc/users/wu04/project/sWGS/results/qdnaseq/sarah_cohort/train/'

# Set the ploidy and cellularity values
# ploidy_list <- seq(1, 4, by=0.5)
# cellularity_list <- seq(0.05, 1, by=0.05)
ploidy_list <- seq(2, 3, by=0.5)
cellularity_list <- seq(0.05, 0.15, by=0.05)

# List all _output.csv files in the input directory
file_list <- list.files(input_dir, pattern = "_output\\.csv$", full.names = TRUE, recursive = TRUE)

# Define the function to convert relative CNA to absolute using given ploidy and cellularity values
relative_to_absolute_given_para <- function(input_file, ploidy_list, cellularity_list){
  result_each_para <- list()
  for (i in seq_along(ploidy_list)){
    for (j in seq_along(cellularity_list)){
      ploidy_i <- ploidy_list[i]
      cellularity_j <- cellularity_list[j]
      #print(ploidy_i)
      #print(cellularity_j)
      data <- read.csv(input_file)
      data$ploidy <- ploidy_i
      data$cellularity <- cellularity_j
      data$absolute_segmented <- relative_to_absolute_copy_number(data$segmented, ploidy_i, cellularity_j)
      data$absolute_copy_number <- relative_to_absolute_copy_number(data$copy_number, ploidy_i, cellularity_j)
      
      result_each_para <- append(result_each_para, list(data))
    }
  }
  results_df <- do.call(rbind, result_each_para)

  #Save the results in bin-size
  output_path_copy_number <- paste0(dirname(input_file), '/absolute_', basename(input_file))
  write.csv(results_df, output_path_copy_number, row.names = FALSE)
  #Save the results in segments
  results_df_segments <- copy_number_segments(results_df)
  output_path_segmented <- paste0(dirname(input_file), '/absolute_segments_', basename(input_file))
  write.csv(results_df_segments, output_path_segmented, row.names = FALSE)
}


# Apply the function to each file
lapply(file_list, function(file_path) {
  relative_to_absolute_given_para(file_path, ploidy_list, cellularity_list)
})
