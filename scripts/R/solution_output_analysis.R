library(readr)
library(dplyr)
library(ggplot2)
library(hrbrthemes)

root_path <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/back_up/results/"

groups <- c('imc', 'hgd', 'ndbe_np')

solutions_df <- lapply (1:length(groups), function(i) {
  result_path <- paste0(root_path, groups[i])
  
  # List all 'auto_solution_choice_auto_qc.tsv' files in the subdirectories
  file_paths <- list.files(path = result_path, pattern = "v2_auto_solution_choice_auto_qc.tsv", full.names = TRUE, recursive = TRUE)
  
  print(dirname(file_paths))
  # Read and filter each file, then combine them
  filtered_data <- lapply(file_paths, function(file) {
    data <- read_tsv(file)
    filtered <- filter(data, selected_solution_auto == TRUE)
    return(filtered)
  }) %>% bind_rows()
  filtered_data <- filtered_data %>% mutate(groups = groups[i]) 
  
  return(filtered_data)
}) %>% bind_rows()

solutions_filtered_df <- solutions_df %>% select(sample, ploidy, purity, num_gds, low_purity, groups)

solutions_single_group <- solutions_filtered_df %>% filter(groups=='ndbe_np')

ggplot(solutions_single_group, aes(x=purity, y=ploidy, color=groups)) + 
  geom_point(size=3, alpha=0.5) +
  ggtitle("Ploidy vs. Purity by Group") +
  theme_ipsum()+
  scale_color_manual(values = c("imc" = "#ff5252", "hgd" = "orange", "ndbe_np" = "green")) +
  xlim(0, 1)+
  ylim(1.5, 4.5)


# generate the data after low_purity force
solutions_df_with_low_purity <- solutions_df %>%
  mutate(
    ploidy = if_else(low_purity, 2, ploidy),  
    purity = if_else(low_purity, 0.05, purity) 
  )

solutions_single_group_with_low_purity <- solutions_df_with_low_purity %>% filter(groups=='hgd')

ggplot(solutions_single_group_with_low_purity, aes(x=purity, y=ploidy, color=groups)) + 
  geom_point(size=3, alpha=0.5) +
  ggtitle("Ploidy vs. Purity by Group") +
  theme_ipsum()+
  scale_color_manual(values = c("imc" = "#ff5252", "hgd" = "orange", "ndbe_np" = "green")) +
  xlim(0, 1)+
  ylim(1.5, 4.5)

#write_csv(filtered_data, paste0(result_path, "/combined_output.csv"))


# Read the CSV files
data1 <- read.csv(paste0(root_path, "/imc/combined_output_imc.csv"))
data2 <- read.csv(paste0(root_path, "/hgd/combined_output.csv"))
data3 <- read.csv(paste0(root_path, "/ndbe_np/combined_output_np_ndbe.csv"))

# Extract the 'ploidy' column and add a group identifier
data1 <- data.frame(purity = data1$purity, group = "imc")
data2 <- data.frame(purity = data2$purity, group = "hgd")
data3 <- data.frame(purity = data3$purity, group = "np_ndbe")

# Combine the data
combined_data <- rbind(data1, data2, data3)

# Create the violin plot
ggplot(combined_data, aes(x = group, y = purity, fill = group)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin plot of Ploidy across Three Groups",
       x = "Group",
       y = "Ploidy") +
  theme_minimal()

