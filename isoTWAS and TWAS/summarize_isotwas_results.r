#########################################################################
# summarize_isotwas_results.R
#
# Purpose: Processes and summarizes isoTWAS results across different cancer types
# and tissues, performs screening and confirmation testing to identify significant 
# isoform-trait associations, and generates summary files for significant findings.
#
# The script implements a two-stage statistical testing approach:
# 1. Screening: Identifies genes with at least one significant isoform
# 2. Confirmation: Tests individual isoforms within significant genes
#########################################################################

# Set working directory to the base results folder
# This should be changed to your specific results directory
base_dir <- "PATH_TO_RESULTS_DIR"
setwd(base_dir)

# Load required packages
require(data.table)
require(isotwas)
require(tidyverse)
require(vroom)

# Connect to Ensembl database to get gene annotations
require(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                     dataset = "hsapiens_gene_ensembl")

# Read manifest file that contains cancer-tissue pairings
manifest <- as.data.frame(
  fread('PATH_TO_MANIFEST/twas_cancer_tissue_full.tsv', header=TRUE)
)

# Clean cancer names by removing MAF suffix if present
manifest$Cancer <- sapply(
  strsplit(manifest$Cancer, '_MAF'),
  function(x) x[1]
)

# Load gnomAD pLI (probability of loss-of-function intolerance) scores for genes
# pLI scores indicate how tolerant a gene is to loss of function mutations
gnomad <- fread('PATH_TO_GNOMAD/gnomad_plof.txt')
gnomad$HGNC <- gnomad$gene  # Ensure HGNC column exists for later merging

# Process each cancer type separately
for (cancer in unique(manifest$Cancer)) {
  # Navigate to the results directory for this cancer
  results_dir <- file.path(base_dir, "isoTWAS", cancer)
  setwd(results_dir)
  print(paste("Processing cancer type:", cancer))
  
  # Clean up directories: keep only directories for tissues associated with this cancer
  # This removes tissues not relevant to this cancer type
  all_dirs <- list.dirs()
  manifest_this_cancer <- subset(manifest, Cancer == cancer)
  dirs_to_remove <- all_dirs[!all_dirs %in% paste0('./', manifest_this_cancer$Tissue)]
  unlink(dirs_to_remove, recursive = TRUE)
  
  # Process each tissue directory for this cancer
  tissue_dirs <- list.dirs()
  tissue_dirs <- tissue_dirs[tissue_dirs != '.']  # Exclude current directory
  
  for (tissue_dir in tissue_dirs) {
    setwd(tissue_dir)
    tissue_name <- basename(tissue_dir)
    print(paste("Processing tissue:", tissue_name))
    
    # Read isoTWAS results file (assumes first file in directory is the results file)
    isotwas_res <- vroom::vroom(list.files()[1], show_col_types = FALSE)
    
    # Query Ensembl for gene metadata
    bm <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol',
                             'chromosome_name', 'start_position', 'end_position',
                             'gene_biotype'),
               filters = 'ensembl_gene_id',
               values = unique(isotwas_res$Gene), 
               mart = ensembl)
    
    # Rename columns for clarity
    colnames(bm) <- c('Gene', 'HGNC', 'Chromosome', 'Start', 'End', 'Biotype')
    
    # Clean up results: remove incomplete rows, sort by Z-score magnitude, remove duplicates
    isotwas_res <- isotwas_res[complete.cases(isotwas_res), ]
    isotwas_res <- isotwas_res[order(abs(isotwas_res$Z), decreasing = TRUE), ]
    isotwas_res <- isotwas_res[!duplicated(isotwas_res$Transcript) & abs(isotwas_res$Z) < Inf, ]
    
    # Rename R2 column for clarity
    colnames(isotwas_res)[8] <- 'R2'
    
    # Merge gene metadata with results
    isotwas_res <- merge(bm, isotwas_res, by='Gene')
    
    # Remove MHC region genes (chr 6: 27Mb-35Mb) which can have spurious associations
    # due to the complex LD structure in this region
    mhc_indices <- which(isotwas_res$Chromosome == 6 &
                         isotwas_res$Start < 35e6 &
                         isotwas_res$End > 27e6)
    if (length(mhc_indices) > 0) {
      isotwas_res <- isotwas_res[-mhc_indices, ]
    }
    
    # Perform gene-level screening using the screening p-value
    # This aggregates isoform p-values to get a gene-level p-value
    gene_summary <- isotwas_res %>%
      group_by(Gene) %>%
      summarise(HGNC = unique(HGNC),
                Chromosome = unique(Chromosome),
                Start = unique(Start),
                End = unique(End),
                Biotype = unique(Biotype),
                Screen.P = isotwas::p_screen(P))  # This implements Simes procedure
    
    # Adjust for multiple testing at gene level using FDR
    alpha1 <- 0.05  # Screening significance threshold
    G <- nrow(gene_summary)  # Total number of genes
    gene_summary$Screen.P.Adjusted <- p.adjust(gene_summary$Screen.P, method = 'fdr')
    
    # Determine significant genes and calculate confirmation threshold
    R <- length(unique(gene_summary$Gene[gene_summary$Screen.P.Adjusted < alpha1]))
    alpha2 <- (R * alpha1) / G  # Dynamic threshold for confirmation stage
    
    # Prepare data frame for storing results
    isoform_new <- as.data.frame(matrix(nrow = 0, ncol = ncol(isotwas_res) + 2))
    colnames(isoform_new) <- c(colnames(isotwas_res), 'Screen.P', 'Confirmation.P')
    
    # Sort genes by screening p-value and merge with isoform results
    gene_summary <- gene_summary[order(gene_summary$Screen.P), ]
    ttt <- merge(isotwas_res,
                gene_summary[, c('Gene', 'Screen.P', 'Screen.P.Adjusted')])
    
    # Perform confirmation testing for each isoform within significant genes
    isoform_new <- ttt %>%
      group_by(Gene) %>%
      summarise(Transcript = Transcript,
                Confirmation.P = isotwas::p_confirm(P, alpha = alpha2))
    
    # Merge confirmation p-values with original results
    isoform_new <- merge(isoform_new, ttt, by=c('Gene', 'Transcript'))
    
    # Set confirmation p-value to 1 for isoforms from non-significant genes
    isoform_new$Confirmation.P <- ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                      isoform_new$Confirmation.P,
                                      1)
    
    # Add tissue information
    isoform_new$Tissue <- tissue_name
    
    # Remove any duplicates
    isoform_new <- isoform_new[!duplicated(isoform_new), ]
    
    # Filter for significant isoforms: must pass screening, confirmation, and permutation tests
    isoform_sig <- subset(isoform_new, 
                        Screen.P < alpha1 &
                          Confirmation.P < alpha2 &
                          permute.P < 0.05)
    
    # Return to cancer directory to write results
    setwd(results_dir)
    
    # Write all results with screening and confirmation p-values
    fwrite(isoform_new,
          'ScreenConfirm_isoTWAS.tsv',
          append=TRUE,
          row.names=FALSE,
          quote=FALSE)
    
    # Process and write significant associations if any exist
    if (nrow(isoform_sig) > 1) {
      # Ensure chromosome is numeric for proper sorting
      isoform_sig$Chromosome <- as.numeric(isoform_sig$Chromosome)
      
      # Sort by genomic position
      isoform_sig <- isoform_sig[order(isoform_sig$Chromosome,
                                     isoform_sig$Start), ]
      
      # Add cancer indication information
      isoform_sig$Indication <- manifest$Indication[manifest$Cancer == cancer][1]
      
      # Select and rename columns for the final output
      isoform_sig <- isoform_sig[, c('Indication', 'Gene', 'HGNC', 'Tissue', 'Chromosome',
                                   'Start', 'End', 'Biotype', 'Transcript',
                                   'Z', 'P', 'permute.P', 'topSNP', 'topSNP.P',
                                   'Screen.P.Adjusted', 'Confirmation.P', 'R2')]
      
      colnames(isoform_sig) <- c('Indication', 'Gene', 'HGNC', 'Tissue', 'Chromosome',
                               'Start', 'End', 'Biotype', 'Transcript',
                               'Z', 'P', 'Permutation P', 'Top GWAS SNP', 'Top GWAS P',
                               'Screening Adjusted P', 'Confirmation P', 'R2')
      
      # Merge with gnomAD pLI scores
      isoform_sig <- merge(isoform_sig,
                         gnomad[, c('HGNC', 'pLI')],
                         by = 'HGNC',
                         all.x=TRUE)
      
      # Write significant associations to file
      fwrite(isoform_sig,
            'SignificantAssociations_isoTWAS_noFineMap.tsv',
            append=TRUE,
            row.names=FALSE,
            quote=FALSE)
    }
  }
  
  # Return to base directory before processing next cancer
  setwd(base_dir)
}
