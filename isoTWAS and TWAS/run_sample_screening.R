#########################################################################
# run_sample_screening.R
#
# Purpose: Run screening and confirmation on isoTWAS analysis results
#
# Usage: Rscript run_sample_screening.R
#########################################################################

# Load required packages
library(data.table)
library(isotwas)
library(dplyr)

# Read configuration
config <- yaml::read_yaml("config_example.txt")

# Create directories if they don't exist
for (dir in c(config$results_dir, config$output_dir)) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Load isotwas results
results_file <- file.path(config$results_dir, config$isotwas_results_file)
if (!file.exists(results_file)) {
  stop("isoTWAS results file not found. Run run_sample_isotwas.R first.")
}

isotwas_res <- fread(results_file)

# Clean up results
isotwas_res <- isotwas_res[complete.cases(isotwas_res), ]
isotwas_res <- isotwas_res[order(abs(isotwas_res$Z), decreasing = TRUE), ]
isotwas_res <- isotwas_res[!duplicated(isotwas_res$Transcript) & abs(isotwas_res$Z) < Inf, ]

# Perform gene-level screening
cat("Performing gene-level screening...\n")
gene_summary <- isotwas_res %>%
  group_by(Gene) %>%
  summarise(
    Chromosome = first(Chromosome),
    Start = first(Start),
    End = first(End),
    Screen.P = isotwas::p_screen(P)
  )

# Adjust for multiple testing at gene level using FDR
alpha1 <- 0.05  # Screening significance threshold
G <- nrow(gene_summary)  # Total number of genes
gene_summary$Screen.P.Adjusted <- p.adjust(gene_summary$Screen.P, method = 'fdr')

# Determine significant genes and calculate confirmation threshold
R <- sum(gene_summary$Screen.P.Adjusted < alpha1)
alpha2 <- (R * alpha1) / G  # Dynamic threshold for confirmation stage

cat("Found", R, "significant genes at screening level\n")
cat("Confirmation threshold:", alpha2, "\n")

# Merge with original results to get transcript-level data
isoform_new <- merge(isotwas_res, gene_summary[, c('Gene', 'Screen.P', 'Screen.P.Adjusted')])

# Perform confirmation testing for each transcript
isoform_confirm <- data.frame()
for (gene in unique(isoform_new$Gene)) {
  gene_data <- subset(isoform_new, Gene == gene)
  
  # Only do confirmation for genes that pass screening
  if (gene_data$Screen.P.Adjusted[1] < alpha1) {
    # Calculate confirmation p-values
    confirm_p <- isotwas::p_confirm(gene_data$P, alpha = alpha2)
    
    # Add to results
    gene_data$Confirmation.P <- confirm_p
  } else {
    # For genes that don't pass screening, set confirmation p-value to 1
    gene_data$Confirmation.P <- 1
  }
  
  isoform_confirm <- rbind(isoform_confirm, gene_data)
}

# Filter for significant isoforms
isoform_sig <- subset(isoform_confirm, 
                      Screen.P.Adjusted < alpha1 &
                        Confirmation.P < alpha2 &
                        permute.P < 0.05)

cat("Found", nrow(isoform_sig), "significant transcripts after screening and confirmation\n")

# Add Indication for compatibility with finemapping
isoform_confirm$Indication <- "sample_trait"
isoform_sig$Indication <- "sample_trait"

# Write all results with screening and confirmation p-values
screen_confirm_file <- file.path(config$results_dir, "ScreenConfirm_isoTWAS.tsv")
fwrite(isoform_confirm, screen_confirm_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Write significant associations
if (nrow(isoform_sig) > 0) {
  sig_file <- file.path(config$results_dir, "SignificantAssociations_isoTWAS_noFineMap.tsv")
  
  # Rename columns for consistency with the pipeline
  colnames(isoform_sig)[colnames(isoform_sig) == "Screen.P.Adjusted"] <- "Screening Adjusted P"
  colnames(isoform_sig)[colnames(isoform_sig) == "Confirmation.P"] <- "Confirmation P"
  colnames(isoform_sig)[colnames(isoform_sig) == "permute.P"] <- "Permutation P"
  colnames(isoform_sig)[colnames(isoform_sig) == "topSNP"] <- "Top GWAS SNP"
  colnames(isoform_sig)[colnames(isoform_sig) == "topSNP.P"] <- "Top GWAS P"
  
  fwrite(isoform_sig, sig_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Significant associations written to", sig_file, "\n")
} else {
  cat("No significant associations found after screening and confirmation\n")
}

cat("Screening and confirmation completed\n")