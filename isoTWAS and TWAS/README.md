# Fine-mapping Pipeline for isoTWAS Sample Data

This repository contains a simplified pipeline for statistical fine-mapping of isoform-level Transcriptome-Wide Association Studies (isoTWAS) results to identify likely causal isoforms in trait-associated genomic regions.

## Overview

The pipeline consists of three main stages:
1. Running isoTWAS to identify isoform-trait associations
2. Performing screening and confirmation testing to control false discoveries
3. Applying statistical fine-mapping to prioritize causal isoforms in regions with multiple associations

## Required Files

### Input Files
- `SampleGWASData.tsv.gz` - GWAS summary statistics file with columns for SNP, CHR, BP, A1, A2, BETA, SE, P, Z
- isoTWAS model files: 
  - `BABAM1_isoTWAS.RDS`
  - `KLF2_isotwas.RDS`
- LD matrix files:
  - `BABAM1_LD_matrix.rds`
  - `KLF2_LD_matrix.rds`

### R Scripts
- `run_sample_isotwas.R` - First tier: Performs isoTWAS analysis on sample data
- `run_sample_screening.R` - Second tier: Applies screening and confirmation testing
- `run_sample_finemapping.R` - Third tier: Runs fine-mapping on overlapping loci

### Configuration
- `config_example.txt` - Configuration file with file paths and settings

## Prerequisites

The pipeline requires the following R packages:
```r
install.packages(c("data.table", "Matrix", "yaml", "dplyr", "boot"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "biomaRt"))
# Install isotwas from GitHub
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("bhattacharya-a-bt/isotwas")
```

## Directory Structure

Before running the pipeline, set up your directory structure as follows:
```
/path/to/working_dir/
├── config_example.txt         # Configuration file
├── SampleGWASData.tsv.gz      # GWAS summary statistics
├── BABAM1_isoTWAS.RDS         # isoTWAS model file
├── KLF2_isotwas.RDS           # isoTWAS model file
├── BABAM1_LD_matrix.rds       # LD matrix file
├── KLF2_LD_matrix.rds         # LD matrix file
├── run_sample_isotwas.R       # Script for tier 1
├── run_sample_screening.R     # Script for tier 2
├── run_sample_finemapping.R   # Script for tier 3
├── results/                   # Directory for intermediate results (will be created)
└── output/                    # Directory for final output (will be created)
```

## Running the Pipeline

### Step 1: Run isoTWAS Analysis

First, run the isoTWAS analysis to generate Z-scores and P-values for each isoform:

```bash
Rscript run_sample_isotwas.R
```

This script:
1. Loads the GWAS summary statistics and isoTWAS models
2. Aligns SNPs between GWAS data and models
3. Calculates Z-scores using the weighted burden test
4. Performs permutation testing to assess significance
5. Outputs results to `results/isotwas_results.tsv`

### Step 2: Run Screening and Confirmation Testing

Next, apply the screening/confirmation testing framework:

```bash
Rscript run_sample_screening.R
```

This script:
1. Performs gene-level screening to identify genes with at least one significant isoform
2. Adjusts for multiple testing using FDR
3. Performs confirmation testing on individual isoforms within significant genes
4. Outputs results to `results/SignificantAssociations_isoTWAS_noFineMap.tsv`

### Step 3: Run Fine-mapping

Finally, perform statistical fine-mapping:

```bash
Rscript run_sample_finemapping.R
```

The fine-mapping procedure:
1. Identifies overlapping loci with multiple associations
2. Calculates the posterior inclusion probability (PIP) for each isoform
3. Determines the credible set of likely causal isoforms
4. Outputs final results to `output/isoTWAS_FineMap_PanCan_[date].tsv`

## Output Files

The pipeline generates several output files:

1. Tier 1 output:
   - `results/isotwas_results.tsv` - isoTWAS association results

2. Tier 2 output:
   - `results/ScreenConfirm_isoTWAS.tsv` - All isoforms with screening/confirmation p-values
   - `results/SignificantAssociations_isoTWAS_noFineMap.tsv` - Significant isoforms that pass both tests

3. Tier 3 output:
   - `output/sample_trait_isoTWAS_FOCUS_Pancan.tsv` - Fine-mapping results for each isoform
   - `output/isoTWAS_FineMap_PanCan.tsv` - Consolidated fine-mapping results
   - `output/isoTWAS_FineMap_PanCan_[date].tsv` - Final timestamped results

## Example Command Sequence

```bash
# Run the complete pipeline
Rscript run_sample_isotwas.R
Rscript run_sample_screening.R
Rscript run_sample_finemapping.R
```

## Notes

- This pipeline is simplified to work with the provided sample data
- The pipeline assumes all files are in the same directory
- Make sure the R scripts and configuration file are in the working directory
- You may need to adjust the configuration file paths depending on your setup

## Interpreting Results

The final output file contains the following columns:
- `Gene`, `Transcript`: Gene name and transcript ID
- `Z`: Z-score for the association
- `pip`: Posterior inclusion probability
- `in_cred_set`: Whether the isoform is in the 90% credible set
- Additional columns with chromosome location, p-values, etc.

Isoforms with high PIPs and inclusion in the credible set are considered the most likely causal variants.
