
dir_path <- "C:/Users/lw748/OneDrive - University of Cambridge/sWGS/QDNASeq_and_model_training/RASCAL/github_rascal/results"

old_name <- "absolute_CNA_plot1.jpg"
new_name <- "v1_absolute_CNA_plot_round_1.jpg" # The absolute CNA based on the previous version of au_qc.R

old_name <- "absolute_CNA_plot2.jpg"
new_name <- "v1_absolute_CNA_plot_round_2.jpg" # The absolute CNA after rerunning the previous version of au_qc.R using solutions close to the selected one in the first run

old_name <- "v2_absolute_CNA_plot1.jpg"
new_name <- "v2_absolute_CNA_plot_round_1.jpg" # The absolute CNA based on the modified version of au_qc.R

old_name <- "v2_absolute_CNA_plot2.jpg"
new_name <- "v2_absolute_CNA_plot_round_2.jpg" # The absolute CNA after rerunning the modified version of au_qc.R using solutions close to the selected one in the first run

# List all files matching the pattern recursively
file_paths <- list.files(dir_path, pattern = paste0("^", old_name, "$"), full.names = TRUE, recursive = TRUE)

# Function to construct new filename based on the old one
construct_new_filename <- function(file_path) {
  dir_part <- dirname(file_path)
  file.path(dir_part, new_name)
}

# Rename each file
for (file_path in file_paths) {
  print(file_path)
  new_file_path <- construct_new_filename(file_path)
  print(new_file_path)
  file.rename(file_path, new_file_path)  # Rename the file
}

