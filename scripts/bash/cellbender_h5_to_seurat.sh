#!/bin/bash

# Initialize the shell for conda
eval "$(conda shell.bash hook)"

# Activate the cellbender conda environment
conda activate cellbender


# Define input and output base directories
input_base_dir="/nfs/turbo/sph-colacino/aguilada/scRNAseq/cellbender/cellbender_h5_outputs_run2"
output_base_dir="/nfs/turbo/sph-colacino/aguilada/scRNAseq/cellbender/cellbender_h5_outputs_for_seurat_run2"

# Create the output base directory if it doesn't exist
mkdir -p "$output_base_dir"

# Iterate over all subdirectories in the input base directory
for sub_dir in "$input_base_dir"/*; do
    if [[ -d "$sub_dir" ]]; then
        # Extract the subdirectory name
        sub_dir_name=$(basename "$sub_dir")
        
        # Only process directories that start with FPR_
        if [[ "$sub_dir_name" != FPR_* ]]; then
            echo "Skipping non-FPR directory: $sub_dir_name"
            continue
        fi
        
        # Create a corresponding subdirectory in the output directory
        output_sub_dir="$output_base_dir/$sub_dir_name"
        mkdir -p "$output_sub_dir"
        
        echo "Processing directory: $sub_dir_name"
        
        # Process each file in the current subdirectory
        for file in "$sub_dir"/*; do
            if [[ -f "$file" ]]; then
                # Generate the output file path
                output_file="$output_sub_dir/$(basename "${file%.h5}_seurat.h5")"
                
                # Run the ptrepack command
                echo "Processing $file -> $output_file"
                ptrepack --complevel 5 "$file:/matrix" "$output_file:/matrix"

                # Check if the command was successful
                if [[ $? -ne 0 ]]; then
                    echo "Error: Failed to process $file"
                    exit 1
                fi
            fi
        done
    fi

done

echo "All files processed successfully."

conda deactivate