#!/bin/bash

#BSUB -J ecaviar[1-100]
#BSUB -P project_code
#BSUB -q normal
#BSUB -n 1
#BSUB -R "rusage[mem=8GB]"
#BSUB -W 4:00
#BSUB -o logs/ecaviar_%J_%I.out
#BSUB -e logs/ecaviar_%J_%I.err

# Usage:
# bsub < run_ecaviar_jobarray.lsf
#
# Note: Update the job array range [1-100] to match the number of genes/isoforms
# Update the project_code to your specific computing cluster project code

# Create logs directory if it doesn't exist
mkdir -p logs

# Load necessary modules (customize based on your computing environment)
module load R/4.1.0
module load tabix/0.2.6
module load ecaviar/2.2

# Set paths to files and directories
CONFIG_FILE="/path/to/config.yaml"
ANALYSIS_TYPE=$1  # "twas" or "isotwas"
TRAIT=$2          # Trait name

# Check command line arguments
if [ -z "$ANALYSIS_TYPE" ] || [ -z "$TRAIT" ]; then
    echo "Usage: bsub -J ecaviar[1-N] < run_ecaviar_jobarray.lsf twas|isotwas trait_name"
    exit 1
fi

# Get index from LSF job array
INDEX=${LSB_JOBINDEX}

echo "Running eCAVIAR colocalization for ${ANALYSIS_TYPE} analysis on trait ${TRAIT}, index ${INDEX}"

# Run the colocalization analysis script
Rscript run_ecaviar_colocalization.R \
    --config ${CONFIG_FILE} \
    --type ${ANALYSIS_TYPE} \
    --trait ${TRAIT} \
    --index ${INDEX}

echo "Job completed."
