library(rascal)
library(tidyverse)
library(dplyr)

# Define the function to process results in segments
copy_number_segments_modified_from_rascal <- function (copy_number) 
{
  stopifnot(is.data.frame(copy_number))
  stopifnot("sample" %in% names(copy_number))
  stopifnot("chromosome" %in% names(copy_number))
  stopifnot("start" %in% names(copy_number), is.numeric(copy_number$start))
  stopifnot("end" %in% names(copy_number), is.numeric(copy_number$end))
  stopifnot("segmented" %in% names(copy_number), is.numeric(copy_number$segmented))
  segmented_CN <- copy_number %>% filter(!is.na(segmented)) %>% mutate(length = end - 
                                                         start + 1) %>% arrange(ploidy, cellularity, sample, chromosome, start) %>% 
    mutate(new_segment = row_number() == 1 | !(sample == 
                                                 lag(sample) & chromosome == lag(chromosome) & segmented == 
                                                 lag(segmented) & ploidy == lag(ploidy) & cellularity == 
                                                 lag(cellularity))) %>% mutate(segment = cumsum(new_segment)) %>% 
    group_by(segment) %>% summarize(sample = dplyr::first(sample), 
                                    chromosome = dplyr::first(chromosome), start = dplyr::first(start), 
                                    end = last(end), copy_number = dplyr::first(segmented), absolute_copy_number = dplyr::first(absolute_segmented), bin_count = n(), 
                                    sum_of_bin_lengths = sum(length), weight = sum(length)/median(length),
                                    ploidy = dplyr::first(ploidy), cellularity = dplyr::first(cellularity))
  return(segmented_CN)
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
  results_df_segments <- copy_number_segments_modified_from_rascal(results_df)
  
  return(list(results_df = results_df, results_df_segments = results_df_segments))
}

