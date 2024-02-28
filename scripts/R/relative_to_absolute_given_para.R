library(rascal)
library(tidyverse)

input_dir <- '/path/to/_output.csv/'

#input_dir <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/github_rascal/results/SLX-10722_SLX-11823.D701_D504"

# Set the ploidy and cellularity values
# ploidy_list <- seq(1, 4, by=0.5)
# cellularity_list <- seq(0.05, 1, by=0.05)
ploidy_list <- seq(2, 3, by=0.5)
cellularity_list <- seq(0.05, 0.15, by=0.05)

# List all _output.csv files in the input directory
file_list <- list.files(input_dir, pattern = "_output\\.csv$", full.names = TRUE, recursive = TRUE)

# Define the function to process results in segments
copy_number_segments_modified_from_rascal <- function (copy_number) 
{
  stopifnot(is.data.frame(copy_number))
  stopifnot("sample" %in% names(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("start" %in% names(copy_number), is.numeric(copy_number$start))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))
  stopifnot("segmented" %in% names(copy_number), is.numeric(copy_number$segmented))
  copy_number %>% filter(!is.na(segmented)) %>% mutate(length = end - 
                                                         start + 1) %>% arrange(ploidy, cellularity, sample, chromosome, start) %>% 
    mutate(new_segment = row_number() == 1 | !(sample == 
                                                 lag(sample) & chromosome == lag(chromosome) & segmented == 
                                                 lag(segmented) & ploidy == lag(ploidy) & cellularity == 
                                                 lag(cellularity))) %>% mutate(segment = cumsum(new_segment)) %>% 
    group_by(segment) %>% summarize(sample = first(sample), 
                                    chromosome = first(chromosome), start = first(start), 
                                    end = last(end), copy_number = first(segmented), absolute_copy_number = first(absolute_segmented), bin_count = n(), 
                                    sum_of_bin_lengths = sum(length), weight = sum(length)/median(length),
                                    ploidy = first(ploidy), cellularity = first(cellularity))
}


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
  results_df_segments <- copy_number_segments_modified_from_rascal(results_df)
  output_path_segmented <- paste0(dirname(input_file), '/absolute_segments_', basename(input_file))
  write.csv(results_df_segments, output_path_segmented, row.names = FALSE)
}


# Apply the function to each file
lapply(file_list, function(file_path) {
  relative_to_absolute_given_para(file_path, ploidy_list, cellularity_list)
})
