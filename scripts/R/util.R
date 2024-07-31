setwd('C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/back_up/scripts/R/')

source('relative_to_absolute_given_para.R')
source('auto_qc.R')
source('find_acceptable_solutions.R')

#input_dir <- '/path/to/_output.csv/'
input_dir <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/back_up/results/imc"
file_list <- list.files(input_dir, pattern = "^binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)

# Output the solutions from RASCAL-------------------------------------------------------------------------------------------------------
## Acceptable solutions from RASCAL
#lapply(file_list, function(file_path) {
#  relative_copy_numbers <- fread(file_path)
#  acceptable_solutions <- find_acceptable_solutions(relative_copy_numbers$segmented)
#  write.table(acceptable_solutions, paste0(dirname(file_path), '/acceptable_solutions.txt'))
#})

# Convert relative cna from QDNASeq to absolute using manually set solutions------------------------------------------------------------------------------

# Set the ploidy and cellularity values
ploidy_list <- seq(1, 2, by=1)
cellularity_list <- seq(0.05, 1, by=0.05)
#ploidy_list <- seq( 1.5, 6, by = 0.1 )
#cellularity_list <- c(0.025, 0.05, 0.075, seq(0.1, 1, by=0.05))

# Apply the function to each file
results_list <- lapply(file_list, function(file_path) {
  relative_to_absolute_given_para(file_path, ploidy_list, cellularity_list)
})

#Save the results in bin-size and in segments
for (i in seq_along(results_list)) {
  print(i)
  results <- results_list[[i]]
  output_path_copy_number <- paste0(dirname(file_list[i]), '/absolute_', basename(file_list[i]))
  output_path_segmented <- paste0(dirname(file_list[i]), '/absolute_segments_', basename(file_list[i]))
  write.csv(results$results_df, output_path_copy_number, row.names = FALSE)
  write.csv(results$results_df_segments, output_path_segmented, row.names = FALSE)
}

# Bind segmented values of all samples 
#segslist <- rbindlist( lapply(results_list, function(x) x[[2]]) )

# Select the best solution  --------------------------------------------------------------------------------------------------
seg_file_list <- list.files(input_dir, pattern = "^absolute_segments_binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)
bin_file_list <- list.files(input_dir, pattern = "^absolute_binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)

i = 19
lapply(1:length(seg_file_list), function(i) {
  print(paste("Processing iteration:", i))
  cn_qc <- select_best_solution(seg_file_list[i])
  #cn_qc <- fread(paste0(dirname(seg_file_list[i]), "/auto_solution_choice_auto_qc.tsv"))

  ### Check if the sample is so low purity that we need to set a max possible purity ###
  # Check to see if there is any segment with a copy number change giving us signal to determine the purity
  # don't use small segments - can be noise
  bins <- fread(bin_file_list[i])
  
  # add segment names
  bins[ !is.na(copy_number),
        segname := paste(chromosome, min(start), max(end), sep = ':'),
        by = .(chromosome, segmented) ]

  # only take first solution (using relative CN so all the same)
  bins <- bins[ paste(ploidy, cellularity) == paste(ploidy, cellularity)[1] ]
  
  # For each segment - is it significantly different from 1
  bins[ !is.na(copy_number), p_value := t.test(copy_number, mu = 1, alternative = 'two.sided')$p.value, by = .(segname, sample) ]
  bins_segs <- unique( bins[ !is.na(segmented), .(sample, p_value, segmented, size = .N * 500000), by = segname ] )
  bins_segs[, q_value := p.adjust(p_value, method = 'bonferroni', n = .N), by = sample]
  bins_segs[, diff := abs(1 - segmented)]
  #bins_segs[, use_seg := size > 20000000 & diff > 0.1 ]
  bins_segs[, use_seg := diff > 0.1 ]
  low_purity <- bins_segs[, !any( q_value[(use_seg)] < 0.01, rm.na=T ) | !any(use_seg), sample]
  low_purity_value <- low_purity$V1[1]
  
  print(low_purity)

  # modify the autoQC sheet and add this solution to the output files
  cn_qc[, low_purity := low_purity_value]

  if( low_purity_value ){
    # Run a new solution forced at an estimated highest purity of 5% and a ploidy of 2
    ploidy_ <- 2
    cellularity_ <- 0.05
    
    results_list_low <- relative_to_absolute_given_para(file_list[i], ploidy_, cellularity_)

    output_path_forced_sol_bin <- paste0(dirname(file_list[i]), '/forced_sol_bin_', i, basename(file_list[i]))
    output_path_forced_sol_segmented <- paste0(dirname(file_list[i]), '/forced_sol_segments_',i, basename(file_list[i]))
    write.csv(results_list_low$results_df, output_path_forced_sol_bin, row.names = FALSE)
    write.csv(results_list_low$results_df_segments, output_path_forced_sol_segmented, row.names = FALSE)
    
    
    }

  # Write out
  data.table::fwrite(cn_qc, paste0(dirname(seg_file_list[i]), "/auto_solution_choice_auto_qc.tsv"), sep = "\t")

  print('============')

})


### END ----------------------------------------------------------------------------------------------------------------------------------
