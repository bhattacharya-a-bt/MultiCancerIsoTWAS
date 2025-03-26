#########################################################################
# run_finemapping.R
#
# Purpose: Run statistical fine-mapping on TWAS and isoTWAS results
# to identify causal genes/isoforms in loci with multiple associations.
#
# Usage: Rscript run_finemapping.R [config_file]
#########################################################################

# Load utility functions
source("finemapping_utilities.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a configuration file")
}

# Read configuration file
config_file <- args[1]
config <- yaml::read_yaml(config_file)

# Extract configuration values
results_dir <- config$results_dir
model_dir <- config$model_dir
reference_file <- config$reference_file
temp_dir <- config$temp_dir
output_dir <- config$output_dir
isotwas_results_file <- file.path(results_dir, config$isotwas_results_file)
twas_results_file <- file.path(results_dir, config$twas_results_file)

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#' Fine-map isoTWAS results
#'
#' @param results_file File containing isoTWAS results
#' @param traits Vector of trait names to process
#' @param output_dir Directory for output files
#' @param model_dir Directory containing prediction models
#' @param reference_file Path to reference genome LD file
#' @param temp_dir Directory for temporary files
#' @return None
finemap_isotwas <- function(results_file, traits, output_dir, model_dir, reference_file, temp_dir) {
  # Read all results
  all_results <- read_assoc_results(results_file)
  
  # Process each trait
  for (trait in traits) {
    cat("Processing isoTWAS results for trait:", trait, "\n")
    
    # Define output file
    output_file <- file.path(output_dir, paste0(trait, "_isoTWAS_FOCUS_Pancan.tsv"))
    
    # Remove existing output file if any
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    
    # Extract results for this trait
    trait_results <- subset(all_results, File == trait)
    
    # Prepare results
    trait_results <- prepare_results(trait_results, trait, is_isotwas = TRUE)
    
    # Identify overlapping loci
    overlap_ids <- identify_overlaps(trait_results)
    
    # Process non-overlapping associations
    process_nonoverlap(trait_results, overlap_ids, output_file, is_isotwas = TRUE)
    
    # Extract overlapping results
    overlap_results <- subset(trait_results, TTA %in% overlap_ids$tta)
    
    # Process overlapping associations
    if (nrow(overlap_results) > 0) {
      process_overlaps_by_tissue(overlap_results, output_file, model_dir, 
                                reference_file, temp_dir, is_isotwas = TRUE)
    }
  }
  
  # Consolidate results across traits
  final_output <- file.path(output_dir, "isoTWAS_FineMap_PanCan.tsv")
  consolidated_file <- consolidate_results(traits, output_dir, final_output, is_isotwas = TRUE)
  
  cat("isoTWAS fine-mapping complete. Final results in:", consolidated_file, "\n")
}

#' Fine-map TWAS results
#'
#' @param results_file File containing TWAS results
#' @param traits Vector of trait names to process
#' @param output_dir Directory for output files
#' @param model_dir Directory containing prediction models
#' @param reference_file Path to reference genome LD file
#' @param temp_dir Directory for temporary files
#' @return None
finemap_twas <- function(results_file, traits, output_dir, model_dir, reference_file, temp_dir) {
  # Read all results
  all_results <- read_assoc_results(results_file)
  
  # Process each trait
  for (trait in traits) {
    cat("Processing TWAS results for trait:", trait, "\n")
    
    # Define output file
    output_file <- file.path(output_dir, paste0(trait, "_TWAS_FOCUS_Pancan.tsv"))
    
    # Remove existing output file if any
    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    
    # Extract results for this trait
    trait_results <- subset(all_results, File == trait)
    
    # Prepare results
    trait_results <- prepare_results(trait_results, trait, is_isotwas = FALSE)
    
    # Identify overlapping loci
    overlap_ids <- identify_overlaps(trait_results)
    
    # Process non-overlapping associations
    process_nonoverlap(trait_results, overlap_ids, output_file, is_isotwas = FALSE)
    
    # Extract overlapping results
    overlap_results <- subset(trait_results, TTA %in% overlap_ids$tta)
    
    # Process overlapping associations
    if (nrow(overlap_results) > 0) {
      process_overlaps_by_tissue(overlap_results, output_file, model_dir, 
                                reference_file, temp_dir, is_isotwas = FALSE)
    }
  }
  
  # Consolidate results across traits
  final_output <- file.path(output_dir, "TWAS_FineMap_PanCan.tsv")
  consolidated_file <- consolidate_results(traits, output_dir, final_output, is_isotwas = FALSE)
  
  cat("TWAS fine-mapping complete. Final results in:", consolidated_file, "\n")
}

# Get list of traits to process
if (!is.null(config$selected_traits) && length(config$selected_traits) > 0) {
  selected_traits <- config$selected_traits
} else {
  # Get unique traits from results files
  isotwas_results <- read_assoc_results(isotwas_results_file)
  twas_results <- read_assoc_results(twas_results_file)
  selected_traits <- unique(c(isotwas_results$File, twas_results$File))
}

# Run fine-mapping
if (config$run_isotwas) {
  finemap_isotwas(isotwas_results_file, selected_traits, output_dir, 
                 model_dir, reference_file, temp_dir)
}

if (config$run_twas) {
  finemap_twas(twas_results_file, selected_traits, output_dir,
              model_dir, reference_file, temp_dir)
}

cat("Fine-mapping complete!\n")
