# absolute_sWGS

Currently, this repository contains scripts designed for two primary purposes:

1. Generating solutions for diploid and cellularity using Rascal.
2. Converting relative Copy Number Alterations (CNA) to absolute values by manually setting diploid and cellularity parameters.

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

3. **Convert Relative CNA to Absolute Values:**
   - To manually set diploid and cellularity values and convert relative CNA to absolute, run `relative_to_absolute_given_para.sh` from `./scripts/bash` with `*_output.csv`.
   - This step also utilizes an R script located in `./scripts/R`.

## Results

You can find examples of `*_output.csv` and the results generated at each step in the `./results` folder.

