# absolute_sWGS

Currently, this repository contains scripts designed for two primary purposes:

1. Generating solutions for diploid and cellularity using Rascal.
2. Converting relative Copy Number Alterations (CNA) to absolute values by manually setting diploid and cellularity values.

## Getting Started

### Preparation

- The input `.rds` data files, directly generated from QDNAseq, are located in the `./data` folder.
- Prepare a `.txt` file containing the sequence IDs of the samples you wish to process.

### Processing Steps

1. **Extract QDNAseq Copy Number Variations (CNVs):**  
   Run `extract_qdna_cnv_givenID.sh` to extract bin-size and segmented copy number, combining them into a CSV file named `*_output.csv`.

2. **Generate Solutions for Diploid and Cellularity:**
   - To use Rascal for generating solutions for diploid and cellularity, run `fit_absolute_copy_numbers.sh` using `*_output.csv` as input.

3. **Convert Relative CAN to Absolute Values:**
   - If you wish to convert relative CAN to absolute by manually setting diploid and cellularity values, use `relative_to_absolute_given_para.sh` with `*_output.csv`.

## Results

- Examples of `*_output.csv` and the sequential results can be found in the `./results` folder.

