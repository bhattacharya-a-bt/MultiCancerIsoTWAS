#########################################################################
# utility_functions.R
#
# Purpose: Contains utility functions for TWAS and isoTWAS analysis
# These functions are used by both the TWAS and isoTWAS scripts to
# prepare data, load models, and perform shared operations.
#########################################################################

# Load required libraries
library(optparse)
library(data.table)
library(bigsnpr)
library(vroom)

#' Read and prepare GWAS summary statistics
#'
#' @param sumstats_file Path to the GWAS summary statistics file
#' @param bim_file Path to the LD reference bim file
#' @return Prepared GWAS summary statistics
read_gwas_sumstats <- function(sumstats_file, bim_file) {
  # Read in reference SNP information
  bim <- fread(bim_file)
  
  # Read the GWAS summary statistics
  sumstats <- fread(sumstats_file)
  
  # Ensure SNP column exists
  sumstats$SNP <- sumstats$rsid
  
  # Filter to SNPs in the reference
  sumstats <- subset(sumstats, SNP %in% bim$V2)
  
  # Standardize allele coding to uppercase
  sumstats$A1 <- toupper(sumstats$A1)
  sumstats$A2 <- toupper(sumstats$A2)
  
  return(sumstats)
}

#' Load LD matrix for a gene
#'
#' @param gene Gene name
#' @param ld_folder Path to folder containing LD matrices
#' @param tissue Tissue name
#' @param type Type of LD matrix ("isotwas" or "twas")
#' @return LD matrix for the gene
load_ld_matrix <- function(gene, ld_folder, tissue, type = "isotwas") {
  # Define path to LD file
  ld_path <- file.path(ld_folder, tissue, paste0("LD_", type), 
                       paste0(type, "_", gene, ".bed"))
  
  # Check if LD file exists
  if (file.exists(ld_path)) {
    # Load the LD matrix
    ld <- snp_attach(snp_readBed2(ld_path, backingfile = tempfile()))
    snpnames <- ld$map$marker.ID
    ld_matrix <- as.matrix(snp_cor(ld$genotypes))
    colnames(ld_matrix) <- rownames(ld_matrix) <- snpnames
    return(ld_matrix)
  } else {
    warning(paste("LD file not found for", gene, "using type", type))
    return(NULL)
  }
}

#' Prepare gene model for association testing
#'
#' @param model The gene prediction model
#' @return Prepared model
prepare_model <- function(model) {
  # Standardize column names
  model$A1 <- model$REF
  model$A2 <- model$ALT
  model$Transcript <- model$Feature
  
  # Ensure data types are correct
  if (nrow(model) > 1) {
    model <- as.data.frame(apply(model, 2, unlist))
  }
  
  model$Chromosome <- as.numeric(model$Chromosome)
  model$Weight <- as.numeric(model$Weight)
  model$R2 <- as.numeric(model$R2)
  model$Position <- as.numeric(model$Position)
  
  return(model)
}

#' Filter GWAS summary statistics to include only SNPs in models
#'
#' @param sumstats GWAS summary statistics
#' @param twas_model TWAS model
#' @param isotwas_model isoTWAS model
#' @return Filtered GWAS summary statistics
filter_sumstats <- function(sumstats, twas_model, isotwas_model) {
  # Get all SNPs from both models
  all_snps <- unique(c(twas_model$SNP, isotwas_model$SNP))
  
  # Filter summary statistics
  sumstats_filtered <- subset(sumstats, SNP %in% all_snps)
  
  # Get position information from models
  pos_info <- rbind(
    twas_model[, c('SNP', 'Chromosome', 'Position')],
    isotwas_model[, c('SNP', 'Chromosome', 'Position')]
  )
  pos_info <- pos_info[!duplicated(pos_info$SNP), ]
  
  # Merge position information with summary statistics
  sumstats_merged <- merge(sumstats_filtered, pos_info, by = 'SNP')
  
  return(sumstats_merged)
}

#' Track processed genes in a done file
#'
#' @param gene Gene name
#' @param trait Trait name
#' @param done_file Path to the done file
#' @return None, writes to file
track_gene_done <- function(gene, trait, done_file) {
  # Mark gene as done
  fwrite(
    data.frame(
      Trait = trait,
      Gene = gene,
      Index = 1
    ),
    done_file,
    sep = '\t',
    row.names = FALSE,
    quote = FALSE,
    append = TRUE
  )
}

#' Check if gene has already been processed
#'
#' @param gene Gene name
#' @param trait Trait name
#' @param done_file Path to the done file
#' @return TRUE if gene has been processed, FALSE otherwise
is_gene_done <- function(gene, trait, done_file) {
  if (!file.exists(done_file)) {
    return(FALSE)
  }
  
  # Read done file
  done_data <- vroom(done_file, show_col_types = FALSE)
  
  # Check if gene is in done file
  return(nrow(subset(done_data, Gene == gene & Trait == trait)) > 0)
}

#' Create a new done file for tracking processed genes
#'
#' @param trait Trait name
#' @param done_file Path to the done file
#' @return None, creates a file
initialize_done_file <- function(trait, done_file) {
  # Create a new done file with a test entry
  fwrite(
    data.frame(
      Trait = trait,
      Gene = 'test',
      Index = 1
    ),
    done_file,
    sep = '\t',
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
}

#' Get list of genes to process for a tissue
#'
#' @param model_folder Path to folder containing models
#' @param done_file Path to the done file
#' @param trait Trait name
#' @return Vector of gene names to process
get_genes_to_process <- function(model_folder, done_file, trait) {
  # Get all model files
  model_files <- list.files(model_folder)
  
  # Extract gene names
  all_genes <- sapply(strsplit(model_files, '_iso'), function(x) x[1])
  
  # If done file exists, filter out already processed genes
  if (file.exists(done_file)) {
    done_data <- vroom(done_file, show_col_types = FALSE)
    done_genes <- subset(done_data, Trait == trait)$Gene
    genes_to_process <- all_genes[!all_genes %in% done_genes]
    return(genes_to_process)
  }
  
  return(all_genes)
}
