#!/bin/bash
#SBATCH --job-name=cellbender_4rm_rawcrv9_job_array
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL,ARRAY_TASKS
#SBATCH --mem-per-gpu=32g
#SBATCH --cpus-per-gpu=4
#SBATCH --gres=gpu:1
#SBATCH --output=/home/aguilada/slurm_output/slurm-%A_%a.out
#SBATCH --account=colacino0
#SBATCH --partition=gpu

# Runs Cellbender on our 60 Flex scRNAseq mammary gland samples.

# Script to submit SLURM job array for CellBender on multiple samples under the 2-week GPU walltime constraint
# to submit to cluster use: sbatch --array=0-59%5 cellbender_jobs_crv9_raw.sh 
# or, if not already in the correct directory, use: sbatch --array=0-59%5 /nfs/turbo/sph-colacino/aguilada/scRNAseq/slurm_jobs/cellbender_jobs_crv9_raw.sh
# 0-59 are array indices, of which there are 60. One for each of our samples
# %5 means 5 jobs will run at a time. This is the limit on our cluster. Can change the 5 to something else or leave out to run all at same time, though,
# limited to how many GPU nodes available on cluster.

# time needed per job depends on the number of cells in sample. 
# From observations, time to completion goes up exponentially with the number of cells in the sample.
# On a gpu partition, with 2 cores per gpu and 32gb of RAM, for these 60 samples,
# <5000 cells takes ~5-10 minutes, ~10000 cells takes 1hr-1.5hr. ~20000 cells takes 2hr-2.5hr, and so on.

# installed cellbender in own conda environment using
# pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@4334e8966217c3591bf7c545f31ab979cdc6590d
# then forked repository and change os.move() and os.replace() to shutil.move() to allow html report to be generated and saved across file systems.
# This is essential for html report compatability with distrubuted computing systems, like HPCs.
# Then installed using cloned fork with 
# git clone https://github.com/davidaguilaratx/CellBender.git
# pip install -e CellBender


# Initialize the shell for conda
eval "$(conda shell.bash hook)"

# Activate the cellbender conda environment
conda activate cellbender

# set up directories, make sure to use raw data from cell ranger v9.0.0 for consistency
directory_3w="/nfs/turbo/sph-colacino/aguilada/scRNAseq/8651-JC-mammary/10x_analysis_JC-8651_cellranger_v.9.0.0/8651-JC_P02/outs/per_sample_outs"
directory_rest="/nfs/turbo/sph-colacino/aguilada/scRNAseq/10x_analysis_11389-DA/10x_analysis_DA-11389_cellranger_v.9.0.0"
output_dir="/nfs/turbo/sph-colacino/aguilada/scRNAseq/cellbender/cellbender_h5_outputs_run2"

# Create the output directory, recursively
mkdir -p $output_dir

# cell bender parameters, expected-cells and total-droplets-included, for each sample
csv_file="/nfs/turbo/sph-colacino/aguilada/scRNAseq/cellbender_parameters.csv"

# Extract parameters from cellbender_parameters.csv
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" $csv_file)
echo "Line "$((SLURM_ARRAY_TASK_ID + 2))" of csv file: $line"

# Parse extracted csv line into variables. IFS means "Internal Field Seperator" and is
# a variable used to define how the shell splits input into fields. Without it we
# would save the output into each variable as something like 0,3000,6000 instead of individual numbers for each paramter like we want.
# Array index is in the first column, expected-cells is in the 6th column, total-droplets-included
# is in the 7th column of the csv file, and learning-rate is the 8th column. Adjust accodringly.
IFS=',' read -r array_index expected_cells total_droplets_included epochs learning_rate <<< "$(echo "$line" | awk -F ',' '{print $1 "," $6 "," $7 "," $8 "," $9}')"

# Validate the parameters
if [[ -z "$expected_cells" || -z "$total_droplets_included" || -z "$epochs" || -z "$learning_rate" ]]; then
    echo "Error: Missing expected_cells, total_droplets_included, epochs, or learning_rate for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID."
    exit 1
fi

# Determine which sample to process based on the job array index
if [ $SLURM_ARRAY_TASK_ID -le 11 ]; then
    sample_number=$((SLURM_ARRAY_TASK_ID + 21))
    input_file="$directory_3w/8651-JC-${sample_number}/count/sample_raw_feature_bc_matrix.h5"
    output_file="$output_dir/output_cellbender_8651-JC-${sample_number}.h5"
else
    sample_number=$((SLURM_ARRAY_TASK_ID - 11))
    if [ $sample_number -le 16 ]; then
        pool_number=1
    elif [ $sample_number -le 32 ]; then
        pool_number=2
    else
        pool_number=3
    fi
    input_file="$directory_rest/11389-DA_Pool0${pool_number}/outs/per_sample_outs/11389-DA-${sample_number}/count/sample_raw_feature_bc_matrix.h5"
    output_file="$output_dir/output_cellbender_11389-DA-${sample_number}.h5"
fi

# Check GPU allocation
echo "Checking GPU availability and usage..."
nvidia-smi

# Verify CUDA availability
python -c 'import torch; assert torch.cuda.is_available(), "CUDA not available"; print("CUDA is available")'

# need unique temporary directory for each job to avoid collisions when running multiple instances of cellbender at once
export TMPDIR=/scratch/colacino_root/colacino0/aguilada/tmp/cellbender_tmp_arrayid_${SLURM_ARRAY_TASK_ID}
mkdir -p "$TMPDIR"

cd $TMPDIR # to make sure all tmp files are writtin in TMPDIR and not the home directory

# Echo parameters for troubleshooting
echo "Troubleshooting Information:"
echo "--------------------------------------------"
echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Array Index (from CSV): $array_index"
echo "Sample Number: $sample_number"
echo "Pool Number: ${pool_number:-N/A}"  # Only defined for specific cases
echo "Expected Cells: $expected_cells"
echo "Total Droplets Included: $total_droplets_included"
echo "Input File Path: $input_file"
echo "Output Directory Path: $output_dir"
echo "Output File Path: $output_file"
echo "Tmp files directory: $TMPDIR"
echo "current working directory:"; pwd
echo "--------------------------------------------"

# Run CellBender

python /nfs/turbo/sph-colacino/aguilada/scRNAseq/slurm_jobs/run_cellbender_seeded.py \
	--cuda \
	--input "$input_file" \
	--output "$output_file" \
	--fpr 0.0 0.01 0.05 0.1 \
	--expected-cells $expected_cells \
	--total-droplets-included $total_droplets_included \
	--epochs $epochs \
	--learning-rate $learning_rate # setting learning rate to 1/2 or 1/4 default as per https://github.com/broadinstitute/CellBender/issues/276. Default learning rate is 0.0001
                                   # It's also suggested to decrease learning rate by factors of 2 and/or set epochs to 300 (default is 150) to improve convergence and PCA structure
# Clean up $TMPDIR
echo "Cleaning up temporary directory: $TMPDIR"
rm -rf "$TMPDIR"

# Confirm removal of $TMPDIR
if [ ! -d "$TMPDIR" ]; then
    echo "$TMPDIR successfully removed."
else
    echo "Failed to remove $TMPDIR."
fi
