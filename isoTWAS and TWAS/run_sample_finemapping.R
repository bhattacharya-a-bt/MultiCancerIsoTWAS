#########################################################################
# run_sample_finemapping.R
#
# Purpose: Run fine-mapping on significant isoTWAS results
#
# Usage: Rscript run_sample_finemapping.R
#########################################################################

# Load required packages
library(data.table)
library(Matrix)
library(isotwas)

# Read configuration
config <- yaml::read_yaml("config_example.txt")

# Load the finemapping_utilities.R file for necessary functions
source("finemapping_utilities.R")

# Create directories if they don't exist
for (dir in c(config$results_dir, config$output_dir)) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Check for significant associations file
sig_file <- file.path(config$results_dir, "SignificantAssociations_isoTWAS_noFineMap.tsv")
if (!file.exists(sig_file)) {
  stop("No significant associations file found. Run run_sample_screening.R first.")
}

# Load significant associations
sig_results <- fread(sig_file)

# Prepare results for fine-mapping
results <- sig_results
results$Trait <- "sample_trait"
results$GTA <- paste(results$Gene, results$Tissue, sep = ":")
results$TTA <- paste(results$Transcript, results$Tissue, sep = ":")

# Identify overlapping loci
chr_unique <- unique(results$Chromosome)
keep_gta <- c()
keep_tta <- c()

for (chr in chr_unique) {
  res_chr <- subset(results, Chromosome == chr)
  res_chr <- res_chr[order(res_chr$Start),]
  
  if (nrow(res_chr) > 1) {
    for (i in 1:(nrow(res_chr) - 1)) {
      if (res_chr$End[i] > res_chr$Start[i+1] - 1e6) {
        keep_gta <- unique(c(keep_gta, c(res_chr$GTA[c(i, i+1)])))
        keep_tta <- unique(c(keep_tta, c(res_chr$TTA[c(i, i+1)])))
      }
    }
  }
}

# Define output file for this trait
output_file <- file.path(config$output_dir, paste0("sample_trait_isoTWAS_FOCUS_Pancan.tsv"))

# Remove existing output file if any
if (file.exists(output_file)) {
  file.remove(output_file)
}

# Process non-overlapping associations
new_res <- subset(results, !TTA %in% keep_tta)
if (nrow(new_res) > 0) {
  new_res$pip <- 1
  new_res$in_cred_set <- TRUE
  new_res$Group <- paste0('ind', 1:nrow(new_res))
  new_res$Overlap <- 'No'
  
  # Ensure column names match expected format
  colnames(new_res)[colnames(new_res) == "Screening Adjusted P"] <- "Screening.Adjusted.P"
  colnames(new_res)[colnames(new_res) == "Confirmation P"] <- "Confirmation.P"
  colnames(new_res)[colnames(new_res) == "Permutation P"] <- "Permutation.P"
  colnames(new_res)[colnames(new_res) == "Top GWAS SNP"] <- "Top.GWAS.SNP"
  colnames(new_res)[colnames(new_res) == "Top GWAS P"] <- "Top.GWAS.P"
  
  # Select columns for output
  cols <- c('Indication', 'Gene', 'Transcript', 'Chromosome',
            'Start', 'End', 'Z', 'Screening.Adjusted.P',
            'Confirmation.P', 'Permutation.P', 'Top.GWAS.SNP', 'Top.GWAS.P',
            'pip', 'in_cred_set', 'Group', 'Overlap')
  
  if (!all(cols %in% colnames(new_res))) {
    missing_cols <- cols[!cols %in% colnames(new_res)]
    cat("Missing columns in results:", paste(missing_cols, collapse = ", "), "\n")
    
    # Add any missing columns
    for (col in missing_cols) {
      new_res[[col]] <- NA
    }
  }
  
  fwrite(new_res[, cols], output_file, append = TRUE, row.names = FALSE, 
         quote = FALSE, sep = '\t', na = "NA")
  
  cat("Processed", nrow(new_res), "non-overlapping associations\n")
}

# Extract overlapping results
overlap_results <- subset(results, TTA %in% keep_tta)

# Process overlapping associations
if (nrow(overlap_results) > 0) {
  cat("Processing", nrow(overlap_results), "overlapping associations\n")
  
  # Get list of unique tissues
  tissues_list <- unique(overlap_results$Tissue)
  
  for (tissue in tissues_list) {
    cat("Processing tissue:", tissue, "\n")
    
    # Get results for current tissue
    res <- subset(overlap_results, Tissue == tissue)
    chr_list <- as.numeric(unique(res$Chromosome))
    
    if (length(chr_list) >= 1) {
      for (chr in chr_list) {
        cat("Processing chromosome:", chr, "\n")
        
        # Group by chromosome
        res_tot <- subset(res, Chromosome == chr)
        res_tot$Group <- 1
        group_id <- 1
        
        if (nrow(res_tot) > 1) {
          for (i in 1:(nrow(res_tot) - 1)) {
            if (res_tot$End[i] <= res_tot$Start[i+1] - 1e6) {
              group_id <- group_id + 1
              res_tot$Group[(i+1):nrow(res_tot)] <- group_id
            }
          }
        }
        
        # Process each group
        for (g in unique(res_tot$Group)) {
          res_group <- subset(res_tot, Group == g)
          
          if (nrow(res_group) == 1) {
            # For single transcript, set PIP to 1
            res_group$pip <- 1
            res_group$in_cred_set <- TRUE
            res_group$Group <- paste0('ind', 1:nrow(res_group))
            res_group$Overlap <- 'No'
            
            # Ensure column names match expected format
            colnames(res_group)[colnames(res_group) == "Screening Adjusted P"] <- "Screening.Adjusted.P"
            colnames(res_group)[colnames(res_group) == "Confirmation P"] <- "Confirmation.P"
            colnames(res_group)[colnames(res_group) == "Permutation P"] <- "Permutation.P"
            colnames(res_group)[colnames(res_group) == "Top GWAS SNP"] <- "Top.GWAS.SNP"
            colnames(res_group)[colnames(res_group) == "Top GWAS P"] <- "Top.GWAS.P"
            
            # Select columns for output
            cols <- c('Indication', 'Gene', 'Transcript', 'Chromosome',
                      'Start', 'End', 'Z', 'Screening.Adjusted.P',
                      'Confirmation.P', 'Permutation.P', 'Top.GWAS.SNP', 'Top.GWAS.P',
                      'pip', 'in_cred_set', 'Group', 'Overlap')
            
            if (!all(cols %in% colnames(res_group))) {
              missing_cols <- cols[!cols %in% colnames(res_group)]
              cat("Missing columns in results:", paste(missing_cols, collapse = ", "), "\n")
              
              # Add any missing columns
              for (col in missing_cols) {
                res_group[[col]] <- NA
              }
            }
            
            fwrite(res_group[, cols], output_file, append = TRUE, row.names = FALSE, 
                   quote = FALSE, sep = '\t', na = "NA")
          } else {
            # For multiple transcripts, perform fine-mapping
            cat("Fine-mapping", nrow(res_group), "transcripts in group", g, "\n")
            
            # Prepare data for fine-mapping
            zscores <- res_group$Z
            
            # Load models and LD matrices to calculate weighted correlation
            all_snps <- c()
            omega <- c()
            gene <- c()
            
            for (i in 1:nrow(res_group)) {
              gene_name <- res_group$Gene[i]
              tx <- res_group$Transcript[i]
              
              # Find model file
              model_file <- config$isotwas_models[grep(gene_name, config$isotwas_models)]
              if (length(model_file) == 0) {
                cat("Model not found for gene", gene_name, "\n")
                next
              }
              
              # Load model
              model <- readRDS(model_file)
              model_subset <- subset(model, Feature == tx)
              
              # Add to collections
              all_snps <- c(all_snps, as.character(model_subset$SNP))
              omega <- c(omega, as.numeric(model_subset$Weight))
              gene <- c(gene, rep(res_group$TTA[i], nrow(model_subset)))
            }
            
            # Create data frame with all SNPs and weights
            tot_df <- data.frame(
              SNP = all_snps,
              Gene = gene,
              Effect = omega
            )
            
            # Create weight matrix
            unique_snps <- unique(all_snps)
            model_df <- as.data.frame(matrix(nrow = length(unique_snps), ncol = nrow(res_group) + 1))
            colnames(model_df) <- c('SNP', res_group$TTA)
            model_df$SNP <- as.character(unique_snps)
            
            for (q in 1:nrow(res_group)) {
              cur_tot_df <- subset(tot_df, Gene == res_group$TTA[q])
              cur_tot_df$SNP <- as.character(cur_tot_df$SNP)
              
              for (i in 1:nrow(model_df)) {
                w <- which(cur_tot_df$SNP == model_df$SNP[i])
                model_df[i, q+1] <- ifelse(length(w) != 0, cur_tot_df$Effect[w], 0)
              }
            }
            
            # Get LD matrix for fine-mapping
            gene_name <- res_group$Gene[1]  # Assuming all in same gene
            ld_file <- config$ld_matrices[grep(gene_name, config$ld_matrices)]
            if (length(ld_file) == 0) {
              cat("LD matrix not found for gene", gene_name, "\n")
              next
            }
            ld_matrix <- readRDS(ld_file)
            
            # Calculate weighted correlation matrix
            Omega <- Matrix(as.matrix(model_df[, -1]))  # Remove SNP column
            
            # Extract subset of LD matrix for SNPs in the model
            snp_indices <- match(model_df$SNP, rownames(ld_matrix))
            snp_indices <- snp_indices[!is.na(snp_indices)]
            ld_subset <- ld_matrix[snp_indices, snp_indices]
            
            # Calculate weighted correlation
            wcor <- tryCatch({
              wcor <- isotwas::estimate_cor(as.matrix(Omega), as.matrix(ld_subset), intercept = TRUE)[[1]]
              diag(wcor) <- 1
              wcor[is.na(wcor)] <- 0
              wcor
            }, error = function(e) {
              cat("Error calculating weighted correlation:", e$message, "\n")
              # Fallback to correlation of 0 between transcripts
              wcor <- diag(nrow(res_group))
              wcor
            })
            
            # Calculate posterior inclusion probabilities
            pips <- tryCatch({
              # Try isotwas::calculate_pips if available
              if (exists("calculate_pips", where = "package:isotwas")) {
                isotwas::calculate_pips(zscores, wcor)$pips
              } else {
                # Implement a simple version if not available
                null_res <- nrow(res_group) * log(1 - 1e-3)
                marginal <- nrow(res_group) * log(1 - 1e-3)
                pips <- rep(0, length(zscores))
                
                # Generate all combinations up to size 3
                for (n in 1:min(3, length(zscores))) {
                  combs <- combn(1:length(zscores), n, simplify = FALSE)
                  for (subset in combs) {
                    sub_z <- zscores[subset]
                    sub_wcor <- wcor[subset, subset, drop = FALSE]
                    
                    # Calculate multivariate normal likelihood
                    log_det <- determinant(sub_wcor, logarithm = TRUE)$modulus[1]
                    inv_wcor <- solve(sub_wcor)
                    maha <- t(sub_z) %*% inv_wcor %*% sub_z
                    
                    # Compute Bayes factor
                    local <- -0.5 * log_det + 0.5 * maha
                    
                    # Update marginal likelihood
                    marginal <- log(exp(local) + exp(marginal))
                    
                    # Update PIPs
                    for (idx in subset) {
                      if (pips[idx] == 0) {
                        pips[idx] <- local
                      } else {
                        pips[idx] <- log(exp(pips[idx]) + exp(local))
                      }
                    }
                  }
                }
                
                # Normalize PIPs
                pips <- exp(pips - marginal)
                pips
              }
            }, error = function(e) {
              cat("Error calculating PIPs:", e$message, "\n")
              # Fallback to equal PIPs
              rep(1/nrow(res_group), nrow(res_group))
            })
            
            # Add PIPs to results
            res_group$pip <- pips
            
            # Sort by PIP in descending order
            res_group <- res_group[order(res_group$pip, decreasing = TRUE), ]
            
            # Determine credible set (90% confidence)
            npost <- res_group$pip / sum(res_group$pip)
            csum <- cumsum(npost)
            res_group$in_cred_set <- FALSE
            
            for (i in 1:nrow(res_group)) {
              if (i == 1 || csum[i-1] < 0.9) {
                res_group$in_cred_set[i] <- TRUE
              }
              if (csum[i] > 0.9) {
                break  # Stop once we reach 90%
              }
            }
            
            res_group$Overlap <- 'Yes'
            
            # Ensure column names match expected format
            colnames(res_group)[colnames(res_group) == "Screening Adjusted P"] <- "Screening.Adjusted.P"
            colnames(res_group)[colnames(res_group) == "Confirmation P"] <- "Confirmation.P"
            colnames(res_group)[colnames(res_group) == "Permutation P"] <- "Permutation.P"
            colnames(res_group)[colnames(res_group) == "Top GWAS SNP"] <- "Top.GWAS.SNP"
            colnames(res_group)[colnames(res_group) == "Top GWAS P"] <- "Top.GWAS.P"
            
            # Select columns for output
            cols <- c('Indication', 'Gene', 'Transcript', 'Chromosome',
                      'Start', 'End', 'Z', 'Screening.Adjusted.P',
                      'Confirmation.P', 'Permutation.P', 'Top.GWAS.SNP', 'Top.GWAS.P',
                      'pip', 'in_cred_set', 'Group', 'Overlap')
            
            if (!all(cols %in% colnames(res_group))) {
              missing_cols <- cols[!cols %in% colnames(res_group)]
              cat("Missing columns in results:", paste(missing_cols, collapse = ", "), "\n")
              
              # Add any missing columns
              for (col in missing_cols) {
                res_group[[col]] <- NA
              }
            }
            
            fwrite(res_group[, cols], output_file, append = TRUE, row.names = FALSE, 
                   quote = FALSE, sep = '\t', na = "NA")
          }
        }
      }
    }
  }
}

# Create final consolidated output
final_output <- file.path(config$output_dir, "isoTWAS_FineMap_PanCan.tsv")
if (file.exists(final_output)) {
  file.remove(final_output)
}

if (file.exists(output_file)) {
  # Read results
  results <- fread(output_file)
  
  # Select columns for final output
  cols <- c('Indication', 'Gene', 'Transcript', 'Chromosome',
            'Start', 'End', 'Z', 'Screening.Adjusted.P',
            'Confirmation.P', 'Permutation.P', 'Top.GWAS.SNP', 'Top.GWAS.P',
            'pip', 'in_cred_set')
  
  # Create timestamp for dated output
  timestamp <- format(Sys.time(), "%m%d%y")
  final_output_dated <- gsub("\\.tsv$", paste0("_", timestamp, ".tsv"), final_output)
  
  # Write final output
  fwrite(results[, cols], final_output_dated, append = FALSE, row.names = FALSE, 
         quote = FALSE, sep = '\t', na = "NA")
  
  cat("Fine-mapping completed. Final results in", final_output_dated, "\n")
} else {
  cat("No results produced. Check for errors.\n")
}