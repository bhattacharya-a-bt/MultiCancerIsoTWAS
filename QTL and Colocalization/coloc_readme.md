# eCAVIAR Colocalization Pipeline

This repository contains utilities for running colocalization analysis using eCAVIAR on TWAS and isoTWAS results to identify genetic variants that may be causal for both GWAS traits and gene/isoform expression.

## Overview

Colocalization analysis is used to determine whether a genetic variant influences both gene/isoform expression and disease risk. The eCAVIAR method estimates the Colocalization Posterior Probability (CLPP) for each variant in a genomic region, helping identify the most likely shared causal variants.

This pipeline:
1. Takes results from TWAS or isoTWAS fine-mapping
2. For each significant gene/isoform, extracts eQTL information from prediction models
3. Obtains GWAS summary statistics for the corresponding region
4. Runs eCAVIAR to compute CLPP scores
5. Collates results across all genes/isoforms

## Prerequisites

The pipeline requires the following:
- R (≥ 4.0) with packages: data.table, yaml, optparse, dplyr
- eCAVIAR executable (http://genetics.cs.ucla.edu/caviar/index.html)
- Tabix for accessing GWAS summary statistics
- LSF job scheduler (for running job arrays)

## Directory Structure

Prepare your directory structure as follows:
```
/path/to/project/
├── config.yaml               # Configuration file
├── run_ecaviar_colocalization.R  # Main script
├── run_ecaviar_jobarray.lsf      # LSF job array script
├── data/
│   ├── results/              # TWAS/isoTWAS results
│   ├── models/               # Gene/isoform prediction models
│   └── gwas/                 # GWAS summary statistics (tabix-indexed)
├── temp/                     # Temporary files
└── output/                   # Output directory
    └── colocalization/       # Colocalization results
```

## Configuration

Create a configuration file `config.yaml` with the following structure:

```yaml
# Directory paths
results_dir: "/path/to/results_dir"
model_dir: "/path/to/model_dir"
gwas_dir: "/path/to/gwas_dir"
temp_dir: "/path/to/temp_dir"
output_dir: "/path/to/output_dir"

# Input files
isotwas_results_file: "AllResults_isoTWAS_PanCan.tsv"
twas_results_file: "AllResults_TWAS_PanCan.tsv"
```

## Running the Pipeline

### Option 1: Run for a single gene/isoform

```bash
Rscript run_ecaviar_colocalization.R \
    --config config.yaml \
    --type isotwas \
    --trait "0_Breast_cancer.AUTO_Finalcleaned" \
    --index 1
```

### Option 2: Run for all genes/isoforms associated with a trait

```bash
Rscript run_ecaviar_colocalization.R \
    --config config.yaml \
    --type isotwas \
    --trait "0_Breast_cancer.AUTO_Finalcleaned"
```

### Option 3: Use LSF job array for parallel processing

```bash
# First determine the number of genes/isoforms N and update the job array range [1-N]
# in the run_ecaviar_jobarray.lsf file

# Submit job array
bsub < run_ecaviar_jobarray.lsf isotwas "0_Breast_cancer.AUTO_Finalcleaned"
```

## Output Files

For each gene/isoform, the pipeline generates a colocalization results file:
- `[gene]_[transcript]_[tissue]_[trait]_coloc.tsv` (for isoTWAS)
- `[gene]_[tissue]_[trait]_coloc.tsv` (for TWAS)

The final consolidated results are saved in:
- `[trait]_isotwas_all_coloc.tsv` (for isoTWAS)
- `[trait]_twas_all_coloc.tsv` (for TWAS)

These files contain the following columns:
- SNP: Variant ID
- CLPP: Colocalization posterior probability
- Chromosome, Position: Genomic coordinates
- Gene, HGNC: Gene identifiers
- Transcript: Transcript ID (for isoTWAS only)
- Tissue: Tissue name
- Indication: Trait/disease
- Additional columns with GWAS and eQTL statistics
