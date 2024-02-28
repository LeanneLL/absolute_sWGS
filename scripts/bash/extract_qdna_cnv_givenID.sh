#!/bin/bash

base_dir="/dir/to/.rds/file"
script_path="/path/to/extract_qdnaseq_copy_number_data.R/script"
sequence_ids_file="/path/to/sequence_id.txt"

while IFS= read -r sequence_id; do
    find "$base_dir/$sequence_id" -type f -name "*segmented.rds" | while read -r input_file; do
        echo "processing $input_file..."
        output_file="${input_file%.*}_output.csv"
#for input_file in "$input_directory"/*segmented.rds; do
    # Construct the output filename by replacing '.rds' with '_output.rds'
    #output_file="$input_directory/$(basename "$input_file" .rds)_output.csv"

    # Call the R script for each file
        Rscript "$script_path" --input "$input_file" --output "$output_file"
        echo "output saved to $output_file"
    done
done < "$sequence_ids_file"
