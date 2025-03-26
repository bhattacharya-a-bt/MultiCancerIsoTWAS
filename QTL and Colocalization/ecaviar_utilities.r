#!/usr/bin/env Rscript

#########################################################################
# run_ecaviar_colocalization.R
#
# Purpose: Run eCAVIAR colocalization analysis on TWAS/isoTWAS results
# to identify genetic variants that may be causal for both GWAS traits
# and gene/isoform expression.
#
# Usage: Rscript run_ecaviar_colocalization.R [config_file] [analysis_type] [trait] [index]
#   config_file: Path to configuration YAML file
#   analysis_type: "twas" or "isotwas"
#   trait: Trait name to analyze (must match what's in the results files)
#   index: Optional index for job array (processes specific gene/isoform)
#########################################################################

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(optparse)
  library(dplyr)
})

# Set up command line options
option_list <- list(
  make_option(c("--config"), action="store", default=NULL, type="character",
              help="Path to configuration file"),
  make_option(c("--type"), action="store", default=NULL, type="character",
              help="Analysis type: 'twas' or 'isotwas'"),
  make_option(c("--trait"), action="store", default=NULL, type="character",
              help="Trait name to analyze"),
  make_option(c("--index"), action="store", default=NULL, type="integer",
              help="Index for job array (to process specific gene/isoform)")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$config) || is.null(opt$type) || is.null(opt$trait)) {
  stop("Required arguments missing. Usage: Rscript run_ecaviar_colocalization.R --config [config_file] --type [twas|isotwas] --trait [trait_name] [--index INDEX]")
}

# Read configuration file
config <- yaml::read_yaml(opt$config)
analysis_type <- opt$type
trait <- opt$trait
index <- opt$index

# Validate analysis type
if (!analysis_type %in% c("twas", "isotwas")) {
  stop("Analysis type must be either 'twas' or 'isotwas'")
}

# Set paths based on configuration
base_dir <- config$results_dir
model_dir <- config$model_dir
output_dir <- config$output_dir
temp_dir <- config$temp_dir
gwas_dir <- config$gwas_dir

# Create output directory if it doesn't exist
coloc_output_dir <- file.path(output_dir, "colocalization", analysis_type, trait)
dir.create(coloc_output_dir, recursive = TRUE, showWarnings = FALSE)

# Define input file paths
if (analysis_type == "isotwas") {
  results_file <- file.path(base_dir, config$isotwas_results_file)
  fine_map_file <- file.path(base_dir, gsub("\\.tsv$", "_FineMap_PanCan.tsv", config$isotwas_results_file))
} else {
  results_file <- file.path(base_dir, config$twas_results_file)
  fine_map_file <- file.path(base_dir, gsub("\\.tsv$", "_FineMap_PanCan.tsv", config$twas_results_file))
}

# Check if files exist
if (!file.exists(results_file)) {
  stop(paste("Results file not found:", results_file))
}
if (!file.exists(fine_map_file)) {
  stop(paste("Fine-mapping file not found:", fine_map_file))
}

# Read fine-mapping results for the specified trait
fine_map_data <- fread(fine_map_file)
if (analysis_type == "isotwas") {
  fine_map_data <- fine_map_data[fine_map_data$Indication == trait, ]
} else {
  fine_map_data <- fine_map_data[fine_map_data$Indication == trait, ]
}

# Only keep features in the 90% credible set
fine_map_data <- fine_map_data[fine_map_data$in_cred_set == TRUE, ]

# If index is provided, process only that specific gene/isoform
if (!is.null(index) && index > 0 && index <= nrow(fine_map_data)) {
  cat("Processing index", index, "of", nrow(fine_map_data), "\n")
  gene_info <- fine_map_data[index, ]
  result <- run_colocalization(gene_info, gwas_dir, model_dir, temp_dir, 
                              coloc_output_dir, analysis_type)
} else {
  # Process all genes/isoforms
  cat("Processing all", nrow(fine_map_data), "genes/isoforms\n")
  
  for (i in 1:nrow(fine_map_data)) {
    gene_info <- fine_map_data[i, ]
    result <- run_colocalization(gene_info, gwas_dir, model_dir, temp_dir, 
                                coloc_output_dir, analysis_type)
  }
  
  # Merge all results
  output_files <- list.files(coloc_output_dir, pattern = "_coloc.tsv$", full.names = TRUE)
  if (length(output_files) > 0) {
    all_results <- rbindlist(lapply(output_files, fread), fill = TRUE)
    final_output <- file.path(coloc_output_dir, paste0(trait, "_", analysis_type, "_all_coloc.tsv"))
    fwrite(all_results, final_output, sep = "\t", quote = FALSE)
    cat("All colocalization results merged to:", final_output, "\n")
  }
}

#' Get GWAS summary statistics for a region
#'
#' @param gwas_file Path to GWAS summary statistics file
#' @param chr Chromosome
#' @param start Start position
#' @param end End position
#' @param window Window size around region (bp)
#' @return Data frame with GWAS summary statistics
get_gwas_region <- function(gwas_file, chr, start, end, window = 1e6) {
  # Extend region by window size
  region_start <- max(1, start - window)
  region_end <- end + window
  
  # Read GWAS summary statistics for the region
  cmd <- paste("tabix", gwas_file, 
               paste0(chr, ":", region_start, "-", region_end))
  
  gwas_region <- try(fread(cmd = cmd))
  
  if (inherits(gwas_region, "try-error")) {
    warning(paste("No GWAS data found for region:", 
                  paste0(chr, ":", region_start, "-", region_end)))
    return(NULL)
  }
  
  # Standardize column names (assuming common GWAS format)
  if ("P" %in% colnames(gwas_region)) {
    gwas_region$pvalue <- gwas_region$P
  } else if ("PVAL" %in% colnames(gwas_region)) {
    gwas_region$pvalue <- gwas_region$PVAL
  } else if ("P_BOLT_LMM" %in% colnames(gwas_region)) {
    gwas_region$pvalue <- gwas_region$P_BOLT_LMM
  }
  
  if ("BETA" %in% colnames(gwas_region)) {
    gwas_region$beta <- gwas_region$BETA
  } else if ("ES" %in% colnames(gwas_region)) {
    gwas_region$beta <- gwas_region$ES
  }
  
  if ("SE" %in% colnames(gwas_region)) {
    gwas_region$se <- gwas_region$SE
  } else if ("SE_BOLT_LMM" %in% colnames(gwas_region)) {
    gwas_region$se <- gwas_region$SE_BOLT_LMM
  }
  
  return(gwas_region)
}

#' Get expression QTLs for a gene/isoform
#'
#' @param gene_id Gene ID
#' @param transcript_id Transcript ID (for isoTWAS)
#' @param tissue Tissue name
#' @param model_dir Directory containing prediction models
#' @param analysis_type "twas" or "isotwas"
#' @return Data frame with eQTL information
get_eqtls <- function(gene_id, transcript_id = NULL, tissue, model_dir, analysis_type) {
  # Set path to model file
  if (analysis_type == "isotwas") {
    model_path <- file.path(model_dir, tissue, "isoTWAS", 
                           paste0(gene_id, "_isoTWAS.RDS"))
  } else {
    model_path <- file.path(model_dir, tissue, "TWAS", 
                           paste0(gene_id, "_TWAS.RDS"))
  }
  
  if (!file.exists(model_path)) {
    warning(paste("Model file not found:", model_path))
    return(NULL)
  }
  
  # Load model
  model <- readRDS(model_path)
  
  # Filter for specific transcript if isoTWAS
  if (analysis_type == "isotwas" && !is.null(transcript_id)) {
    model <- subset(model, Feature == transcript_id)
  }
  
  # Extract eQTL information
  eqtls <- data.frame(
    SNP = model$SNP,
    Chromosome = model$Chromosome,
    Position = model$Position,
    Weight = model$Weight,
    REF = model$REF,
    ALT = model$ALT
  )
  
  # Remove duplicates and SNPs with zero weights
  eqtls <- subset(eqtls, Weight != 0)
  eqtls <- eqtls[!duplicated(eqtls$SNP), ]
  
  return(eqtls)
}

#' Run eCAVIAR colocalization
#'
#' @param gwas_data GWAS summary statistics
#' @param eqtl_data eQTL data
#' @param output_prefix Output file prefix
#' @param temp_dir Temporary directory
#' @param caviar_path Path to eCAVIAR executable
#' @param max_causal Maximum number of causal variants
#' @return Path to colocalization results file
run_ecaviar <- function(gwas_data, eqtl_data, output_prefix, temp_dir, 
                        caviar_path = "eCAVIAR", max_causal = 5) {
  # Create temporary files
  gwas_z_file <- file.path(temp_dir, paste0(output_prefix, "_gwas.z"))
  eqtl_z_file <- file.path(temp_dir, paste0(output_prefix, "_eqtl.z"))
  ld_file <- file.path(temp_dir, paste0(output_prefix, "_ld.txt"))
  
  # Create output paths
  output_file <- file.path(temp_dir, paste0(output_prefix, "_ecaviar"))
  final_output <- paste0(output_file, "_col.txt")
  
  # Merge GWAS and eQTL data based on variant ID
  merged_data <- merge(gwas_data, eqtl_data, by = "SNP")
  
  if (nrow(merged_data) < 10) {
    warning("Too few overlapping SNPs for colocalization")
    return(NULL)
  }
  
  # Calculate Z-scores for GWAS
  merged_data$gwas_z <- merged_data$beta / merged_data$se
  
  # Normalize eQTL weights to Z-scores (simplified approach)
  merged_data$eqtl_z <- scale(merged_data$Weight)
  
  # Write Z-scores to files
  gwas_z_df <- data.frame(SNP = merged_data$SNP, Z = merged_data$gwas_z)
  eqtl_z_df <- data.frame(SNP = merged_data$SNP, Z = merged_data$eqtl_z)
  
  fwrite(gwas_z_df, gwas_z_file, sep = " ", col.names = FALSE)
  fwrite(eqtl_z_df, eqtl_z_file, sep = " ", col.names = FALSE)
  
  # For LD, we would normally compute this from genotype data
  # But for simplicity, we'll use an identity matrix as a placeholder
  n_snps <- nrow(merged_data)
  ld_matrix <- diag(n_snps)
  ld_df <- as.data.frame(ld_matrix)
  fwrite(ld_df, ld_file, sep = " ", col.names = FALSE)
  
  # Run eCAVIAR
  cmd <- paste(caviar_path, "-z", gwas_z_file, "-z", eqtl_z_file,
               "-l", ld_file, "-o", output_file, "-c", max_causal)
  
  system(cmd)
  
  # Check if output file was created
  if (!file.exists(final_output)) {
    warning("eCAVIAR did not produce output file")
    return(NULL)
  }
  
  return(final_output)
}

#' Parse eCAVIAR results
#'
#' @param ecaviar_file Path to eCAVIAR output file
#' @param snp_data Data frame with SNP information
#' @return Data frame with colocalization results
parse_ecaviar_results <- function(ecaviar_file, snp_data) {
  if (is.null(ecaviar_file) || !file.exists(ecaviar_file)) {
    return(NULL)
  }
  
  # Read eCAVIAR output
  coloc <- fread(ecaviar_file)
  colnames(coloc) <- c("SNP", "CLPP")
  
  # Merge with SNP information
  coloc_results <- merge(coloc, snp_data, by = "SNP")
  
  # Sort by CLPP (descending)
  coloc_results <- coloc_results[order(coloc_results$CLPP, decreasing = TRUE), ]
  
  return(coloc_results)
}

# Main function to run colocalization for a gene/isoform
run_colocalization <- function(gene_info, gwas_dir, model_dir, temp_dir, 
                              coloc_output_dir, analysis_type) {
  cat("Processing", ifelse(analysis_type == "isotwas", "isoform", "gene"), 
      ":", gene_info$Gene, "\n")
  
  # Set paths
  gwas_file <- file.path(gwas_dir, paste0(gene_info$Indication, ".gwas.bgz"))
  
  if (!file.exists(gwas_file)) {
    cat("GWAS file not found:", gwas_file, "\n")
    return(NULL)
  }
  
  # Get GWAS summary stats for the region
  gwas_data <- get_gwas_region(
    gwas_file = gwas_file,
    chr = gene_info$Chromosome,
    start = gene_info$Start,
    end = gene_info$End,
    window = 1e6
  )
  
  if (is.null(gwas_data) || nrow(gwas_data) == 0) {
    cat("No GWAS data for region around", gene_info$Gene, "\n")
    return(NULL)
  }
  
  # Get eQTLs
  if (analysis_type == "isotwas") {
    eqtl_data <- get_eqtls(
      gene_id = gene_info$Gene,
      transcript_id = gene_info$Transcript,
      tissue = gene_info$Tissue,
      model_dir = model_dir,
      analysis_type = analysis_type
    )
  } else {
    eqtl_data <- get_eqtls(
      gene_id = gene_info$Gene,
      tissue = gene_info$Tissue,
      model_dir = model_dir,
      analysis_type = analysis_type
    )
  }
  
  if (is.null(eqtl_data) || nrow(eqtl_data) == 0) {
    cat("No eQTL data for", gene_info$Gene, "\n")
    return(NULL)
  }
  
  # Define output prefix
  if (analysis_type == "isotwas") {
    output_prefix <- paste(gene_info$Gene, gene_info$Transcript, 
                          gene_info$Tissue, gene_info$Indication, sep = "_")
  } else {
    output_prefix <- paste(gene_info$Gene, gene_info$Tissue, 
                          gene_info$Indication, sep = "_")
  }
  
  # Run eCAVIAR
  ecaviar_output <- run_ecaviar(
    gwas_data = gwas_data,
    eqtl_data = eqtl_data,
    output_prefix = output_prefix,
    temp_dir = temp_dir
  )
  
  # Parse results
  coloc_results <- parse_ecaviar_results(ecaviar_output, 
                                       merge(gwas_data, eqtl_data, by = "SNP"))
  
  if (is.null(coloc_results) || nrow(coloc_results) == 0) {
    cat("No colocalization results for", gene_info$Gene, "\n")
    return(NULL)
  }
  
  # Add gene information
  coloc_results$Gene <- gene_info$Gene
  coloc_results$HGNC <- gene_info$HGNC
  coloc_results$Chromosome <- gene_info$Chromosome
  coloc_results$Tissue <- gene_info$Tissue
  coloc_results$Indication <- gene_info$Indication
  
  if (analysis_type == "isotwas") {
    coloc_results$Transcript <- gene_info$Transcript
  }
  
  # Write results
  output_file <- file.path(coloc_output_dir, paste0(output_prefix, "_coloc.tsv"))
  fwrite(coloc_results, output_file, sep = "\t", quote = FALSE)
  
  cat("Colocalization results written to:", output_file, "\n")
  return(output_file)
}