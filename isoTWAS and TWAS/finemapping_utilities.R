#########################################################################
# finemapping_utilities.R
#
# Purpose: Utility functions for statistical fine-mapping on TWAS and 
# isoTWAS results to identify causal genes/isoforms.
#########################################################################

# Load required packages
require(data.table)
require(GenomicRanges)
require(Matrix)
require(bigsnpr)
require(isotwas)

#' Read association results from a file
#'
#' @param file_path Path to the file containing association results
#' @return data.frame of association results
read_assoc_results <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste0("File not found: ", file_path))
  }
  
  results <- fread(file_path)
  return(results)
}

#' Prepare results for fine-mapping
#'
#' @param results Data frame containing association results
#' @param trait Trait name
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return Processed data frame
prepare_results <- function(results, trait, is_isotwas = TRUE) {
  # Handle missing pLI values
  results$pLI[is.na(results$pLI)] <- 0
  
  # Sort by chromosome and start position
  results <- results[order(results$Chromosome, results$Start),]
  
  # Add trait information
  results$Trait <- trait
  
  # For TWAS, set Transcript equal to Gene
  if (!is_isotwas) {
    results$Transcript <- results$Gene
  }
  
  # Create identifier columns
  results$GTA <- paste(results$Gene, results$Tissue, sep = ":")
  results$TTA <- paste(results$Transcript, results$Tissue, sep = ":")
  
  return(results)
}

#' Identify overlapping loci
#'
#' @param results Prepared results data frame
#' @return List of overlapping TTAs and GTAs
identify_overlaps <- function(results) {
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
  
  return(list(gta = keep_gta, tta = keep_tta))
}

#' Process non-overlapping associations
#'
#' @param results Prepared results data frame
#' @param overlap_ids List of overlapping IDs
#' @param output_file Output file path
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return None, results are written to file
process_nonoverlap <- function(results, overlap_ids, output_file, is_isotwas = TRUE) {
  new_res <- subset(results, !TTA %in% overlap_ids$tta)
  
  if (nrow(new_res) > 0) {
    new_res$pip <- 1
    new_res$in_cred_set <- TRUE
    new_res$Group <- paste0('ind', 1:nrow(new_res))
    new_res$Overlap <- 'No'
    
    if (is_isotwas) {
      cols <- c('Indication', 'Gene', 'HGNC', 'Transcript', 'Chromosome',
                'Start', 'End', 'Tissue', 'Z', 'Screening Adjusted P',
                'Confirmation P', 'Permutation P', 'Top GWAS SNP', 'Top GWAS P', 'pLI',
                'pip', 'in_cred_set', 'Group', 'Overlap')
    } else {
      cols <- c('Indication', 'Gene', 'HGNC', 'Chromosome',
                'Start', 'End', 'Tissue', 'Z', 'FDR', 'Permutation P',
                'Top GWAS SNP', 'Top GWAS P', 'pLI',
                'pip', 'in_cred_set', 'Group', 'Overlap')
    }
    
    fwrite(new_res[, cols], output_file, append = TRUE, row.names = FALSE, 
           quote = FALSE, sep = '\t', na = "NA")
  }
}

#' Process specific chromosome region
#'
#' @param results_base Base results data frame
#' @param chr Chromosome number
#' @param start Start position
#' @param end End position
#' @param output_file Output file path
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return Updated results_base with specific region removed
process_specific_region <- function(results_base, chr, start, end, output_file, is_isotwas = TRUE) {
  res_region <- subset(results_base, Chromosome == chr & Start >= start & End <= end)
  
  if (nrow(res_region) > 0) {
    res_region$pip <- 1
    res_region$in_cred_set <- TRUE
    res_region$Group <- paste0('ind', 1:nrow(res_region))
    res_region$Overlap <- 'No'
    
    if (is_isotwas) {
      cols <- c('Indication', 'Gene', 'HGNC', 'Transcript', 'Chromosome',
                'Start', 'End', 'Tissue', 'Z', 'Screening Adjusted P',
                'Confirmation P', 'Permutation P', 'Top GWAS SNP', 'Top GWAS P', 'pLI',
                'pip', 'in_cred_set', 'Group', 'Overlap')
    } else {
      cols <- c('Indication', 'Gene', 'HGNC', 'Chromosome',
                'Start', 'End', 'Tissue', 'Z', 'FDR', 'Permutation P',
                'Top GWAS SNP', 'Top GWAS P', 'pLI',
                'pip', 'in_cred_set', 'Group', 'Overlap')
    }
    
    fwrite(res_region[, cols], output_file, append = TRUE, row.names = FALSE, 
           quote = FALSE, sep = '\t', na = "NA")
    
    # Remove processed regions from results_base
    if (is_isotwas) {
      results_base <- subset(results_base, !(paste0(Transcript, Tissue) %in% paste0(res_region$Transcript, res_region$Tissue)))
    } else {
      results_base <- subset(results_base, !(paste0(Gene, Tissue) %in% paste0(res_region$Gene, res_region$Tissue)))
    }
  }
  
  return(results_base)
}

#' Group associations by chromosome
#'
#' @param results Results data frame
#' @return Results with group assignments
group_by_chromosome <- function(results, chr) {
  res_tot <- subset(results, Chromosome == chr)
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
  
  return(res_tot)
}

#' Process single association
#'
#' @param results Results data frame
#' @param output_file Output file path
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return None, results are written to file
process_single_assoc <- function(results, output_file, is_isotwas = TRUE) {
  results$pip <- 1
  results$in_cred_set <- TRUE
  results$Group <- paste0('ind', 1:nrow(results))
  results$Overlap <- 'No'
  
  if (is_isotwas) {
    cols <- c('Indication', 'Gene', 'HGNC', 'Transcript', 'Chromosome',
              'Start', 'End', 'Tissue', 'Z', 'Screening Adjusted P',
              'Confirmation P', 'Permutation P', 'Top GWAS SNP', 'Top GWAS P', 'pLI',
              'pip', 'in_cred_set', 'Group', 'Overlap')
  } else {
    cols <- c('Indication', 'Gene', 'HGNC', 'Chromosome',
              'Start', 'End', 'Tissue', 'Z', 'FDR', 'Permutation P',
              'Top GWAS SNP', 'Top GWAS P', 'pLI',
              'pip', 'in_cred_set', 'Group', 'Overlap')
  }
  
  fwrite(results[, cols], output_file, append = TRUE, row.names = FALSE, 
         quote = FALSE, sep = '\t', na = "NA")
}

#' Load prediction model for gene/transcript
#'
#' @param tissue Tissue name
#' @param gene Gene name
#' @param transcript Transcript name
#' @param model_dir Directory containing prediction models
#' @param is_isotwas Boolean indicating if models are isoTWAS (TRUE) or TWAS (FALSE)
#' @return Data frame with model information
load_prediction_model <- function(tissue, gene, transcript, model_dir, is_isotwas = TRUE) {
  model_type <- ifelse(is_isotwas, "isoTWAS", "TWAS")
  model_path <- file.path(model_dir, tissue, model_type, paste0(gene, "_", model_type, ".RDS"))
  
  if (!file.exists(model_path)) {
    stop(paste0("Model file not found: ", model_path))
  }
  
  model_data <- readRDS(model_path)
  model_data <- subset(model_data, Feature == transcript)
  
  model_df <- data.frame(
    SNP = model_data$SNP,
    Chromosome = model_data$Chromosome,
    Position = model_data$Position,
    Effect = model_data$Weight,
    A1 = model_data$ALT,
    A2 = model_data$REF
  )
  
  model_df <- subset(model_df, Effect != 0)
  model_df <- model_df[!duplicated(model_df$SNP), ]
  
  return(model_df)
}

#' Prepare SNP weight matrix for multiple genes/transcripts
#'
#' @param results Results data frame for a group
#' @param model_dir Directory containing prediction models
#' @param is_isotwas Boolean indicating if models are isoTWAS (TRUE) or TWAS (FALSE)
#' @return List containing model data frame and position information
prepare_snp_weights <- function(results, model_dir, is_isotwas = TRUE) {
  all_snps <- c()
  omega <- c()
  pos <- c()
  gene <- c()
  snp_chr <- c()
  
  for (i in 1:nrow(results)) {
    feature_id <- ifelse(is_isotwas, results$Transcript[i], results$Gene[i])
    model_df <- load_prediction_model(results$Tissue[i], results$Gene[i], 
                                      feature_id, model_dir, is_isotwas)
    
    all_snps <- c(all_snps, as.character(model_df$SNP))
    omega <- c(omega, as.numeric(model_df$Effect))
    gene <- c(gene, rep(results$TTA[i], nrow(model_df)))
    snp_chr <- c(snp_chr, as.numeric(model_df$Chromosome))
    pos <- c(pos, as.numeric(model_df$Position))
  }
  
  tot_df <- data.frame(
    SNP = all_snps,
    Gene = gene,
    Effect = omega,
    Chromosome = snp_chr
  )
  
  # Create weight matrix
  unique_snps <- unique(all_snps)
  model_df <- as.data.frame(matrix(nrow = length(unique_snps), ncol = nrow(results) + 1))
  colnames(model_df) <- c('SNP', results$TTA)
  model_df$SNP <- as.character(unique_snps)
  
  for (q in 1:nrow(results)) {
    cur_tot_df <- subset(tot_df, Gene == results$TTA[q])
    cur_tot_df$SNP <- as.character(cur_tot_df$SNP)
    
    for (i in 1:nrow(model_df)) {
      w <- which(cur_tot_df$SNP == model_df$SNP[i])
      model_df[i, q+1] <- ifelse(length(w) != 0, cur_tot_df$Effect[w], 0)
    }
  }
  
  # Add chromosome information
  model_df$Chromosome <- results$Chromosome[1]  # Assuming all in the same chromosome
  for (i in 1:nrow(model_df)) {
    rrr <- subset(tot_df, SNP == model_df$SNP[i])
    if (nrow(rrr) > 0) {
      model_df$Chromosome[i] <- rrr$Chromosome[1]
    }
  }
  
  min_pos <- max(c(1, min(pos) - 1e6))
  max_pos <- max(pos) + 1e6
  
  return(list(
    model_df = model_df,
    min_pos = min_pos,
    max_pos = max_pos,
    pos = pos,
    chr = results$Chromosome[1]
  ))
}

#' Generate LD matrix for SNPs
#'
#' @param weight_data Weight data from prepare_snp_weights function
#' @param reference_file Path to reference genome LD file
#' @param temp_dir Directory for temporary files
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return List containing correlation matrices and SNP information
generate_ld_matrix <- function(weight_data, reference_file, temp_dir, is_isotwas = TRUE) {
  temp_prefix <- ifelse(is_isotwas, "temp", "temptwas")
  temp_file <- file.path(temp_dir, temp_prefix)
  
  # Extract PLINK data for the region
  plink_cmd <- paste0(
    "plink --bfile ", reference_file,
    " --chr ", weight_data$chr,
    " --from-bp ", weight_data$min_pos,
    " --to-bp ", weight_data$max_pos,
    " --make-bed --out ", temp_file
  )
  
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # Read bed file and extract SNP correlation
  snps <- snp_attach(snp_readBed2(paste0(temp_file, ".bed"), backingfile = tempfile()))
  snp_set <- subset(snps$map, marker.ID %in% weight_data$model_df$SNP)
  
  # Match model SNPs with snp_set
  weight_data$model_df <- weight_data$model_df[match(snp_set$marker.ID, weight_data$model_df$SNP), ]
  
  # Extract genotype correlations
  V <- snp_cor(
    snp_attach(
      subset(
        snps,
        ind.col = which(snps$map$marker.ID %in% weight_data$model_df$SNP),
        backingfile = tempfile()
      )
    )$genotypes
  )
  
  # Create weight matrix
  Omega <- Matrix(as.matrix(weight_data$model_df[, -c(1, ncol(weight_data$model_df))]))
  
  return(list(
    Omega = Omega,
    V = V,
    model_df = weight_data$model_df
  ))
}

#' Calculate posterior inclusion probabilities
#'
#' @param zscores Vector of Z-scores
#' @param wcor Weighted correlation matrix
#' @param max_combo Maximum size of combinations to consider
#' @return List containing PIPs and other stats
calculate_pips <- function(zscores, wcor, max_combo = 3) {
  m <- length(zscores)
  null_res <- m * log(1 - 1e-3)
  marginal <- m * log(1 - 1e-3)
  
  # Generate all combinations up to max_combo
  comb_list <- list()
  for (n in 1:min(max_combo, length(zscores))) {
    comb_list <- c(comb_list, combn(1:length(zscores), n, simplify = FALSE))
  }
  
  # Initialize PIPs
  pips <- rep(0, length(zscores))
  
  # Calculate Bayes factors
  for (j in 1:length(comb_list)) {
    subset <- comb_list[[j]]
    local <- bayes_factor(zscores, idx_set = subset, wcor = wcor)
    marginal <- log(exp(local) + exp(marginal))
    
    for (idx in subset) {
      if (pips[idx] == 0) {
        pips[idx] <- local
      } else {
        pips[idx] <- log(exp(pips[idx]) + exp(local))
      }
    }
  }
  
  # Handle infinite values by reducing combination size
  iter <- max_combo - 1
  while (any(is.infinite(pips)) && iter >= 1) {
    # Reset
    null_res <- m * log(1 - 1e-3)
    marginal <- m * log(1 - 1e-3)
    
    # Generate smaller combinations
    comb_list <- list()
    for (n in 1:min(iter, length(zscores))) {
      comb_list <- c(comb_list, combn(1:length(zscores), n, simplify = FALSE))
    }
    
    # Recalculate PIPs
    pips <- rep(0, length(zscores))
    for (j in 1:length(comb_list)) {
      subset <- comb_list[[j]]
      local <- bayes_factor(zscores, idx_set = subset, wcor = wcor)
      marginal <- log(exp(local) + exp(marginal))
      
      for (idx in subset) {
        if (pips[idx] == 0) {
          pips[idx] <- local
        } else {
          pips[idx] <- log(exp(pips[idx]) + exp(local))
        }
      }
    }
    
    iter <- iter - 1
  }
  
  # Handle remaining infinite values
  if (is.infinite(marginal)) {
    marginal <- max(pips) + 1
  }
  
  # Normalize PIPs
  pips <- exp(pips - marginal)
  null_res <- exp(null_res - marginal)
  
  return(list(
    pips = pips,
    null_res = null_res,
    marginal = marginal
  ))
}

#' Determine credible set
#'
#' @param results Results data frame
#' @param pips Vector of posterior inclusion probabilities
#' @return Updated results data frame with credible set information
determine_credible_set <- function(results, pips) {
  # Add PIPs to results
  results$pip <- pips
  
  # Sort by PIP in descending order
  results <- results[order(results$pip, decreasing = TRUE), ]
  
  # Calculate normalized posterior
  npost <- results$pip / sum(results$pip)
  csum <- cumsum(npost)
  
  # Determine credible set (features with cumulative prob <= 0.9)
  results$in_cred_set <- FALSE
  
  for (i in 1:nrow(results)) {
    results$in_cred_set[i] <- TRUE
    
    if (i > 1) {
      if (csum[i] > 0.9 && csum[i-1] < 0.9) {
        results$in_cred_set[i] <- TRUE
      }
      if (csum[i] < 0.9) {
        results$in_cred_set[i] <- TRUE
      }
      if (csum[i] > 0.9 && csum[i-1] > 0.9) {
        results$in_cred_set[i] <- FALSE
      }
    }
  }
  
  results$Overlap <- 'Yes'
  return(results)
}

#' Process overlapping associations by tissue
#'
#' @param overlap_results Results data frame with overlapping associations
#' @param output_file Output file path
#' @param model_dir Directory containing prediction models
#' @param reference_file Path to reference genome LD file
#' @param temp_dir Directory for temporary files
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return None, results are written to file
process_overlaps_by_tissue <- function(overlap_results, output_file, model_dir, 
                                       reference_file, temp_dir, is_isotwas = TRUE) {
  tissues_list <- unique(overlap_results$Tissue)
  
  for (tissue in tissues_list) {
    # Get results for current tissue
    res <- subset(overlap_results, Tissue == tissue)
    chr_list <- as.numeric(unique(res$Chromosome))
    
    if (length(chr_list) >= 1) {
      for (chr in chr_list) {
        # Group by chromosome
        res_tot <- group_by_chromosome(res, chr)
        
        # Process each group
        for (group in unique(res_tot$Group)) {
          res_group <- subset(res_tot, Group == group)
          
          if (nrow(res_group) == 1) {
            # Process single association
            process_single_assoc(res_group, output_file, is_isotwas)
          } else {
            # Process multiple associations
            weight_data <- prepare_snp_weights(res_group, model_dir, is_isotwas)
            ld_data <- generate_ld_matrix(weight_data, reference_file, temp_dir, is_isotwas)
            
            # Calculate weighted correlation
            wcor <- estimate_cor(as.matrix(ld_data$Omega), as.matrix(ld_data$V), intercept = TRUE)[[1]]
            diag(wcor) <- 1
            wcor[is.na(wcor)] <- 0
            
            # Get residuals
            swld <- estimate_cor(as.matrix(ld_data$Omega), as.matrix(ld_data$V), intercept = TRUE)[[2]]
            
            # Calculate PIPs
            zscores <- res_group$Z
            zscores_resid <- get_resid(zscores, as.matrix(swld), as.matrix(wcor))[[1]]
            pip_results <- calculate_pips(zscores_resid, wcor)
            
            # Determine credible set
            res_group <- determine_credible_set(res_group, pip_results$pips)
            
            # Select output columns
            if (is_isotwas) {
              cols <- c('Indication', 'Gene', 'HGNC', 'Transcript', 'Chromosome',
                        'Start', 'End', 'Tissue', 'Z', 'Screening Adjusted P',
                        'Confirmation P', 'Permutation P', 'Top GWAS SNP', 'Top GWAS P', 'pLI',
                        'pip', 'in_cred_set', 'Group', 'Overlap')
            } else {
              cols <- c('Indication', 'Gene', 'HGNC', 'Chromosome',
                        'Start', 'End', 'Tissue', 'Z', 'FDR', 'Permutation P',
                        'Top GWAS SNP', 'Top GWAS P', 'pLI',
                        'pip', 'in_cred_set', 'Group', 'Overlap')
            }
            
            # Write to output file
            fwrite(res_group[, cols], output_file, append = TRUE, row.names = FALSE, 
                   quote = FALSE, sep = '\t', na = "NA")
          }
        }
      }
    }
  }
}

#' Consolidate results across traits
#'
#' @param traits Vector of trait names
#' @param output_dir Directory containing output files
#' @param final_output Final output file path
#' @param is_isotwas Boolean indicating if results are from isoTWAS (TRUE) or TWAS (FALSE)
#' @return Path to consolidated results file
consolidate_results <- function(traits, output_dir, final_output, is_isotwas = TRUE) {
  # Remove existing final output file if it exists
  if (file.exists(final_output)) {
    file.remove(final_output)
  }
  
  # Process each trait
  for (trait in traits) {
    if (is_isotwas) {
      trait_file <- file.path(output_dir, paste0(trait, "_isoTWAS_FOCUS_Pancan.tsv"))
      result_cols <- c('Indication', 'Gene', 'HGNC', 'Transcript', 'Chromosome',
                       'Start', 'End', 'Tissue', 'Z', 'Screening Adjusted P',
                       'Confirmation P', 'Permutation P', 'Top GWAS SNP', 'Top GWAS P', 'pLI',
                       'pip', 'in_cred_set')
    } else {
      trait_file <- file.path(output_dir, paste0(trait, "_TWAS_FOCUS_Pancan.tsv"))
      result_cols <- c('Indication', 'Gene', 'HGNC', 'Chromosome',
                       'Start', 'End', 'Tissue', 'Z', 'FDR', 'Permutation P',
                       'Top GWAS SNP', 'Top GWAS P', 'pLI',
                       'pip', 'in_cred_set')
    }
    
    if (file.exists(trait_file)) {
      trait_results <- fread(trait_file)
      trait_results <- trait_results[order(trait_results$Chromosome, trait_results$Start), ]
      
      # Write to final output
      fwrite(trait_results[, result_cols], final_output, append = TRUE,
             row.names = FALSE, quote = FALSE, sep = '\t')
    }
  }
  
  # Remove duplicates from final output
  if (file.exists(final_output)) {
    all_results <- fread(final_output)
    
    if (nrow(all_results) > 0) {
      # Remove duplicates
      all_results <- all_results[!duplicated(all_results), ]
      
      if (is_isotwas) {
        all_results$TTA <- paste(all_results$Transcript, all_results$Tissue, all_results$Indication, sep = ':')
        all_results <- all_results[!duplicated(all_results$TTA), ]
      } else {
        all_results$GTA <- paste(all_results$Gene, all_results$Tissue, all_results$Indication, sep = ':')
        all_results <- all_results[order(all_results$pip, decreasing = TRUE), ]
        all_results <- all_results[!duplicated(all_results$GTA), ]
      }
      
      # Write final deduplicated output
      timestamp <- format(Sys.time(), "%m%d%y")
      final_output_dated <- gsub("\\.tsv$", paste0("_", timestamp, ".tsv"), final_output)
      
      fwrite(all_results[, result_cols], final_output_dated, append = FALSE,
             row.names = FALSE, quote = FALSE, sep = '\t')
      
      return(final_output_dated)
    }
  }
  
  return(final_output)
}