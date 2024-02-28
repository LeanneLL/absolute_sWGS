#!/bin/bash

INPUT_dir="/dir/path/to/_output.csv/" #the csv is the output of extract_qdnaseq_copy_number_data.sh
script_path="/path/to/fit_absolute_copy_numbers.R/script"

#modify the basename of _output.csv file
find "$INPUT_dir" -type f -name "binsize_500.copy_number_segmented_output.csv" | while read csv_file; do
    dir_name=$(dirname "$csv_file")
    echo "$dir_name"
    base_name=$(basename "$csv_file" .copy_number_segmented_output.csv)
    output_file="$dir_name/${base_name}_solutions.txt"
    Rscript "$script_path" -i "$csv_file" -o "$output_file"

    echo "Processed $csv_file -> $output_file"
done
