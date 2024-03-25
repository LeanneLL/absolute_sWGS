# absolute_sWGS

Currently, this repository contains scripts designed for two primary purposes:

1. Generating solutions for diploid and cellularity using Rascal.
2. Draw the Copy Number Alterations (CNA) plots, heatmaps and CNA density plots.
3. Converting relative CNA to absolute values by manually setting diploid and cellularity parameters, and select the best one.

## Directory Structure

- **Data:** The input `.rds` data files generated from QDNAseq are located in the `./data` folder.
- **Scripts:** 
  - Bash scripts for processing are found in `./scripts/bash`.
  - Corresponding R scripts are located in `./scripts/R`.

## Getting Started

### Preparation

Before beginning, ensure you have a `.txt` file with the sequence IDs of the samples you plan to process. The input data should be placed in the `./data` directory.

### Processing Steps

1. **Extract QDNAseq CNVs:**
   - To extract bin-size and segmented copy number data, run the script `extract_qdna_cnv_givenID.sh` found in `./scripts/bash`. This script combines the extracted data into a CSV file named `*_output.csv`.
   - The corresponding R script used in this step is located at `./scripts/R`.

2. **Generate Solutions for Diploid and Cellularity:**
   - Use `fit_absolute_copy_numbers.sh` from `./scripts/bash` with `*_output.csv` to generate solutions for diploid and cellularity using Rascal.
   - The R script for this process can be found in `./scripts/R`.

3. **Draw plots:**
   - The R script for CNA plots, heatmaps and CNA density plots can be found in `./scripts/R/util_plot.R`.
     
4. **Convert Relative CNA to Absolute Values by forcing diploid and cellularity parameters, and select the best one:**
   - The R script for this process can be found in `./scripts/R/util.R`.

## Results

You can find examples of `*_output.csv` and the results generated at each step in the `./results` folder.


