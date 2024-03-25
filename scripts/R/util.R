source('relative_to_absolute_given_para.R')
source('auto_qc.R')
source('find_acceptable_solutions.R')

#input_dir <- '/path/to/_output.csv/'
input_dir <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/github_rascal/results"

# Output the solutions from RASCAL-------------------------------------------------------------------------------------------------------
file_list <- list.files(input_dir, pattern = "^binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)
## Acceptable solutions from RASCAL
lapply(file_list, function(file_path) {
  relative_copy_numbers <- fread(file_path)
  acceptable_solutions <- find_acceptable_solutions(relative_copy_numbers$segmented)
  write.table(acceptable_solutions, paste0(dirname(file_path), '/acceptable_solutions.txt'))
})

# Convert relative cna from QDNASeq to absolute using manually set solutions------------------------------------------------------------------------------

# Set the ploidy and cellularity values
#ploidy_list <- seq(1, 8, by=0.2)
#cellularity_list <- seq(0.05, 1, by=0.02)
ploidy_list <- seq(2, 3, by=0.5)
cellularity_list <- seq(0.05, 0.15, by=0.05)

# Apply the function to each file
results_list <- lapply(file_list, function(file_path) {
  relative_to_absolute_given_para(file_path, ploidy_list, cellularity_list)
})

#Save the results in bin-size and in segments
for (i in seq_along(results_list)) {
  results <- results_list[[i]]
  output_path_copy_number <- paste0(dirname(file_list[i]), '/absolute_', basename(file_list[i]))
  output_path_segmented <- paste0(dirname(file_list[i]), '/absolute_segments_', basename(file_list[i]))
  write.csv(results$results_df, output_path_copy_number, row.names = FALSE)
  write.csv(results$results_df_segments, output_path_segmented, row.names = FALSE)
}

# Select the best solution  --------------------------------------------------------------------------------------------------
seg_file_list <- list.files(input_dir, pattern = "^absolute_segments_binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)

lapply(seg_file_list, function(file_path) {
  cn_qc <- select_best_solution(file_path)

  # Write out
  data.table::fwrite(cn_qc, paste0(dirname(file_path), "/auto_solution_choice_auto_qc.tsv"), sep = "\t")
  print('============')
})


### END ----------------------------------------------------------------------------------------------------------------------------------
