# Fine-mapping Pipeline for TWAS and isoTWAS Results

This repository contains a pipeline for statistical fine-mapping of Transcriptome-Wide Association Studies (TWAS) and isoform-level TWAS (isoTWAS) results to identify likely causal genes and isoforms in trait-associated genomic regions.

## Overview

The pipeline consists of three main stages:
1. Running TWAS and/or isoTWAS to identify gene/isoform-trait associations
2. Performing screening and confirmation testing (for isoTWAS only) to control false discoveries
3. Applying statistical fine-mapping to prioritize causal genes/isoforms in regions with multiple associations

## Prerequisites

The pipeline requires the following R packages:
```r
install.packages(c("data.table", "GenomicRanges", "Matrix", "bigsnpr", "yaml", "vroom", "isotwas"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```

## Directory Structure

Before running the pipeline, set up your directory structure as follows:
```
/path/to/base_dir/
├── models/              # Directory containing prediction models
│   └── [tissue]/       
│       ├── TWAS/       # TWAS models
│       └── isoTWAS/    # isoTWAS models
├── results/             # Directory for association results
├── reference/           # LD reference files
│   └── reference.bim
├── temp/                # Temporary files directory
└── output/              # Final output directory
```

## Configuration

Create a configuration file based on the provided `config_example.txt`:

```yaml
# Directory paths
results_dir: "/path/to/results_dir"
model_dir: "/path/to/model_dir"
reference_file: "/path/to/reference_file"
temp_dir: "/path/to/temp_dir"
output_dir: "/path/to/output_dir"

# Input files
isotwas_results_file: "AllResults_isoTWAS_PanCan.tsv"
twas_results_file: "AllResults_TWAS_PanCan.tsv"

# Configuration options
run_isotwas: true
run_twas: true

# Optional: Specify traits to process (leave empty to process all)
selected_traits:
  - "0_Breast_cancer.AUTO_Finalcleaned"
  - "0_Cimba_Ovarian_serous.AUTO_Finalcleaned"
  - "0_ovarian.European_overall.AUTO_Finalcleaned"
```

## Pipeline Workflow

### Step 1: Running TWAS and isoTWAS Analysis

Use the provided utility scripts to perform basic TWAS/isoTWAS analysis:

```bash
# Run TWAS for a specific cancer-tissue pair
Rscript run_twas.R --index 1

# Run isoTWAS for a specific cancer-tissue pair
Rscript run_isotwas.R --index 1
```

These scripts generate association results for each gene/isoform-trait pair.

### Step 2: Summarizing isoTWAS Results with Screening/Confirmation

For isoTWAS results, use the summarize_isotwas_results.R script to apply the screening/confirmation testing framework:

```bash
Rscript summarize_isotwas_results.R
```

This script:
1. Screens for genes with at least one significant isoform association
2. Performs confirmation testing on individual isoforms within significant genes
3. Generates summary files for significant associations

### Step 3: Running Fine-mapping

Use the run_finemapping.R script to perform statistical fine-mapping:

```bash
Rscript run_finemapping.R /path/to/config.yaml
```

The fine-mapping procedure:
1. Identifies overlapping loci with multiple associations
2. For each locus, calculates the posterior inclusion probability (PIP) for each gene/isoform
3. Determines the credible set of likely causal genes/isoforms
4. Produces consolidated results files

## Output Files

The pipeline generates several output files:

1. Trait-specific fine-mapping results:
   - `[trait]_isoTWAS_FOCUS_Pancan.tsv`
   - `[trait]_TWAS_FOCUS_Pancan.tsv`

2. Consolidated results across all traits:
   - `isoTWAS_FineMap_PanCan.tsv`
   - `TWAS_FineMap_PanCan.tsv`

3. Final timestamped results:
   - `isoTWAS_FineMap_PanCan_[date].tsv`
   - `TWAS_FineMap_PanCan_[date].tsv`

## Key Files and Functions

- `run_finemapping.R`: Main script for running the fine-mapping pipeline
- `finemapping_utilities.R`: Utility functions for fine-mapping operations
- `isotwas_functions.R`: Functions for running isoTWAS analysis
- `twas_functions.R`: Functions for running TWAS analysis
- `utility_functions.R`: Shared utility functions for data preparation and processing
- `summarize_isotwas_results.R`: Script to process and summarize isoTWAS results

## Example Usage

1. Set up your configuration file with appropriate paths
2. Run TWAS/isoTWAS analysis for all cancer-tissue pairs
3. Summarize isoTWAS results with screening/confirmation
4. Run fine-mapping:

```bash
# Example complete workflow
for i in {1..10}; do
  Rscript run_twas.R --index $i
  Rscript run_isotwas.R --index $i
done

Rscript summarize_isotwas_results.R

Rscript run_finemapping.R config.yaml
```

## Fine-mapping Details

The fine-mapping procedure:
1. Uses Bayesian inference to compute posterior inclusion probabilities
2. Handles overlapping loci by considering the correlation structure among genes/isoforms
3. Creates 90% credible sets of genes/isoforms most likely to be causal
4. Provides a well-calibrated statistical framework for multiple testing across multiple tissues and traits

## Citation

If you use this pipeline, please cite the appropriate references for TWAS, isoTWAS, and fine-mapping methodologies.

## Contact

For questions or issues, please open an issue in the GitHub repository.
