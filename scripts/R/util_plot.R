library(rascal)
library(data.table)
source('relative_to_absolute_given_para.R')

input_dir <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/github_rascal/results"
file_list <- list.files(input_dir, pattern = "^binsize_500\\.copy_number_segmented_output\\.csv$", full.names = TRUE, recursive = TRUE)

# Draw cna plot based on relative copy number-------------------------------------------------------------------------------------------

lapply(file_list, function(file_path) {
  cna_df <- fread(file_path)
  cna_df <- cna_df[ !(cna_df$chromosome) %in% c(23, 24, 'chrX', 'chrY', 'X', 'Y') ]
  cna_df$copy_number <- log2(cna_df$copy_number)
  cna_df$segmented <- log2(cna_df$segmented)
  
  segmented_cna_df <- copy_number_segments(cna_df)
  cna_df$chromosome <- as.numeric(cna_df$chromosome)
  segmented_cna_df$chromosome <- as.numeric(segmented_cna_df$chromosome)
  
  genome_plot <- genome_copy_number_plot (cna_df, segmented_cna_df)
  ylims <- c(-2, 2)
  genome_plot <- genome_plot + coord_cartesian(ylim = ylims) + labs(y = "log2(Relative copy number)")
  ggsave(filename = paste0(dirname(file_path), "/genome_copy_number_plot.png"), plot = genome_plot, 
         width = 10, 
         height = 6)
})

# Draw cna density plot---------------------------------------------------------------------------------------
lapply(file_list, function(file_path) {
  cna_df <- fread(file_path)
  density_plot <- copy_number_density_plot(cna_df$segmented) + coord_cartesian(xlim = c(0, 2.5))
  ggsave(filename = paste0(dirname(file_path), "/Segmented_copy_number_density.png"), plot = density_plot, 
         width = 10, 
         height = 6)
})

# Draw heatmaps -------------------------------------------------------------------------------------------
lapply(file_list, function(file_path) {
  cna_df <- fread(file_path)
  segmented_cna_df <- copy_number_segments(cna_df)
  distances <- absolute_copy_number_distance_grid(segmented_cna_df$copy_number)
  heat_map_cna <- distance_heatmap (distances)
  ggsave(filename = paste0(dirname(file_path), "/heatmap.png"), plot = heat_map_cna, 
         width = 10, 
         height = 6)
})