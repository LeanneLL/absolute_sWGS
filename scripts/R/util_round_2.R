setwd('C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/back_up/scripts/R/')

source('relative_to_absolute_given_para.R')
source('auto_qc.R')
source('find_acceptable_solutions.R')


#input_dir <- '/path/to/_output.csv/'
input_dir <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/back_up/results/ndbe_np"

solution_list <- list.files(input_dir, pattern = "auto_solution_choice_auto_qc\\.tsv$", full.names = TRUE, recursive = TRUE)
file_list <- list.files(input_dir, pattern = "^binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)

# Read the solutions in the first round, and generate absolute CNA using closer ploidy and cellularity------------------------------------------------
generate_range <- function(x) {
  start <- x - 0.1
  stop <- x + 0.1
  if (start <= 0) {
    start <- 0.01
  }
  
  return(seq(start, stop, by = 0.01))
}

for (i in seq_along(solution_list)) {
  print(paste("Processing iteration for reading solution and reapply the au_qc:", i))
  solution_df <- read_tsv(solution_list[i], show_col_types = FALSE)
  selected_solution <- filter(solution_df, selected_solution_auto == TRUE)
  selected_ploidy <- selected_solution$ploidy
  selected_purity <- selected_solution$purity
  
  sample_name <- dirname(solution_list[i])
  cat("selected ploidy for sample", sample_name, selected_ploidy, selected_purity, "")
  
  ploidy_list <- generate_range(selected_ploidy)
  cellularity_list <- generate_range(selected_purity)
  #cat(ploidy_list, cellularity_list)
  
  # Apply the function to each file
  result_list <- relative_to_absolute_given_para(file_list[i], ploidy_list, cellularity_list)
  
  output_path_copy_number <- paste0(dirname(file_list[i]), '/absolute_round2_', basename(file_list[i]))
  output_path_segmented <- paste0(dirname(file_list[i]), '/absolute_round2_segments_', basename(file_list[i]))
  write.csv(result_list$results_df, output_path_copy_number, row.names = FALSE)
  write.csv(result_list$results_df_segments, output_path_segmented, row.names = FALSE)
  
}

# Select the best solution from closer ploidy and cellularity --------------------------------------------------------------------------------------------------
seg_file_list <- list.files(input_dir, pattern = "^absolute_round2_segments_binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)
bin_file_list <- list.files(input_dir, pattern = "^absolute_round2_binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)

lapply(1:length(seg_file_list), function(i) {
  print(paste("Processing iteration for reapplying qc functions:", i))
  cn_qc <- select_best_solution(seg_file_list[i])
  data.table::fwrite(cn_qc, paste0(dirname(seg_file_list[i]), "/auto_solution_round2_choice_auto_qc.tsv"), sep = "\t")
  
  print('======generate the aboslute CNA based on the optimal solution======')
  selected_solution <- filter(cn_qc, selected_solution_auto == TRUE)
  selected_ploidy <- selected_solution$ploidy
  selected_purity <- selected_solution$purity
  
  sample_name <- dirname(seg_file_list[i])
  cat("selected ploidy for sample", sample_name, selected_ploidy, selected_purity)
  
  relative_CNA_path <- paste0(dirname(seg_file_list[i]), "/binsize_500.copy_number_segmented_output.csv")
  result_list <- relative_to_absolute_given_para(relative_CNA_path, selected_ploidy, selected_purity)
  bin_df <- result_list$results_df
  bin_df <- bin_df %>% filter(!chromosome %in% c("X", "Y")) %>%
    mutate(chromosome = as.numeric(chromosome))
  
  seg_df <- result_list$results_df_segments
  seg_df <- seg_df %>% filter(!chromosome %in% c("X", "Y")) %>%
    mutate(chromosome = as.numeric(chromosome))
    
  output_path_best_sol_bin <- paste0(dirname(seg_file_list[i]), '/best_sol_bin_', i, basename(seg_file_list[i]))
  output_path_best_sol_segmented <- paste0(dirname(seg_file_list[i]), '/best_sol_segments_',i, basename(seg_file_list[i]))
  write.csv(bin_df, output_path_best_sol_bin, row.names = FALSE)
  write.csv(seg_df, output_path_best_sol_segmented, row.names = FALSE)
  
  absolute_CNA_plot <- genome_copy_number_plot (bin_df, seg_df, min_copy_number = 0, max_copy_number = 2)
  
  # Generate the lines indicating absolute copy number
  for (j in 0:9) {
    relative_j <- absolute_to_relative_copy_number(j, selected_ploidy, selected_purity)
    absolute_CNA_plot <- absolute_CNA_plot +
      geom_hline(yintercept = relative_j, color = "blue", size = 0.3, alpha = 0.6) +
      annotate("text", x = Inf, y = relative_j+0.03, label = j, hjust = 1.1, color = "black", size = 3, alpha = 0.7)
  }

  plot(absolute_CNA_plot)
  ggsave(paste0(dirname(seg_file_list[i]), "/absolute_CNA_plot.jpg"), plot = absolute_CNA_plot, width = 10, height = 4, dpi = 300)
})

# Draw the absolute copy number using the best solution


### END ----------------------------------------------------------------------------------------------------------------------------------
