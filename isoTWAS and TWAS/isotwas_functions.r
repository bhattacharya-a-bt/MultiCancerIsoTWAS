#########################################################################
# run_isotwas.R
#
# Purpose: Performs isoform-level transcriptome-wide association studies (isoTWAS)
# for specified cancer-tissue pairs using precomputed prediction models and
# GWAS summary statistics.
#
# Usage: Rscript run_isotwas.R --index [INDEX]
#########################################################################

# Load utility functions
source("utility_functions.R")

# Load required package
library(isotwas)

# Setup command line options
option_list <- list(
  make_option(c("-i", "--index"), action="store_true", default=TRUE,
              help="index value to select which cancer-tissue pair to analyze", 
              type='numeric')
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Set parameters
index <- opt$index

#' Main function to run isoTWAS analysis
#'
#' @param manifest_file Path to the manifest file
#' @param index Index of the cancer-tissue pair to analyze
#' @param base_folder Base folder for all data
#' @param gwas_folder Folder containing GWAS summary statistics
#' @param output_folder Folder for output files
run_isotwas_analysis <- function(manifest_file, index, base_folder, gwas_folder, output_folder) {
  # Read manifest file
  manifest <- as.data.frame(fread(manifest_file, header=TRUE))
  
  # Get cancer type and tissue for this index
  cancer <- manifest$Cancer[index]
  tissue <- manifest$Tissue[index]
  
  # Print what we're analyzing
  cat("Running isoTWAS for", cancer, "in", tissue, "\n")
  
  # Set paths for this analysis
  isotwas_folder <- file.path(base_folder, tissue, "TWAS")
  bim_file <- file.path(base_folder, "reference", "reference.bim")
  ld_folder <- file.path(base_folder, "LD_matrices")
  
  # Find GWAS file for this cancer
  gwas_files <- list.files(gwas_folder)
  gwas_files <- gwas_files[grepl('_MAF0.01', gwas_files)]
  trait_names <- sapply(strsplit(gwas_files, '_MAF'), function(x) x[1])
  trait <- trait_names[trait_names == cancer]
  sumstats_file <- file.path(gwas_folder, gwas_files[grep(cancer, gwas_files)])
  
  # Set output paths
  out_file <- file.path(output_folder, trait, tissue, paste0(trait, "_isoTWAS.tsv"))
  done_file <- file.path(output_folder, paste0(trait, "_", tissue, "_done.tsv"))
  
  # Create output directory if it doesn't exist
  dir.create(file.path(output_folder, trait, tissue), recursive = TRUE, showWarnings = FALSE)
  
  # Initialize done file if needed
  if (!file.exists(done_file)) {
    initialize_done_file(trait, done_file)
  }
  
  # Read GWAS summary statistics
  sumstats <- read_gwas_sumstats(sumstats_file, bim_file)
  
  # Get list of genes to process
  genes_to_process <- get_genes_to_process(isotwas_folder, done_file, trait)
  
  # Process each gene
  for (gene in genes_to_process) {
    # Check if gene has already been processed
    if (is_gene_done(gene, trait, done_file)) {
      cat("Skipping", gene, "- already processed\n")
      next
    }
    
    cat("Processing", gene, "\n")
    
    # Try to load isoTWAS model
    model_file <- file.path(isotwas_folder, paste0(gene, "_isoTWAS.RDS"))
    if (!file.exists(model_file)) {
      cat("No isoTWAS model found for", gene, "\n")
      track_gene_done(gene, trait, done_file)
      next
    }
    
    # Load and prepare isoTWAS model
    isotwas_model <- readRDS(model_file)
    isotwas_model$R2 <- unlist(isotwas_model$R2)
    isotwas_model <- prepare_model(isotwas_model)
    
    # Check if model has good prediction performance
    if (all(isotwas_model$R2 <= 0.01)) {
      cat("isoTWAS model for", gene, "has R² <= 0.01, skipping\n")
      track_gene_done(gene, trait, done_file)
      next
    }
    
    # Load LD matrix for this gene
    ld_matrix <- load_ld_matrix(gene, ld_folder, tissue, type = "isotwas")
    if (is.null(ld_matrix)) {
      cat("Could not load LD matrix for", gene, "\n")
      next
    }
    
    # Filter isoTWAS model to only include SNPs in summary statistics
    isotwas_model <- subset(isotwas_model, SNP %in% sumstats$SNP)
    
    # Keep only isoforms with good prediction performance
    isotwas_model <- subset(isotwas_model, R2 >= 0.01)
    
    if (nrow(isotwas_model) == 0) {
      cat("No isoforms with R² >= 0.01 for", gene, "\n")
      track_gene_done(gene, trait, done_file)
      next
    }
    
    # Filter summary statistics to only include SNPs in isoTWAS model
    gene_sumstats <- subset(sumstats, SNP %in% isotwas_model$SNP)
    
    # Make sure we have position information
    pos_info <- isotwas_model[, c('SNP', 'Chromosome', 'Position')]
    pos_info <- pos_info[!duplicated(pos_info$SNP), ]
    gene_sumstats <- merge(gene_sumstats, pos_info, by = 'SNP')
    
    # Process each transcript/isoform
    for (tx in unique(isotwas_model$Transcript)) {
      cat("Running isoTWAS association for", gene, "transcript", tx, "\n")
      
      # Get model for this transcript
      tx_model <- subset(isotwas_model, Transcript == tx)
      
      # Run isoTWAS association test
      tx_results <- tryCatch({
        burdenTest(
          mod = tx_model,
          ld = ld_matrix,
          gene = gene,
          sumStats = gene_sumstats,
          chr = 'Chromosome',
          pos = 'Position',
          a1 = 'A1',
          a2 = 'A2',
          Z = 'Z',
          beta = 'Beta',
          se = 'SE',
          R2cutoff = 0.01,
          alpha = 1e-3,
          nperms = 1e3,
          usePos = FALSE
        )
      }, error = function(e) {
        cat("Error in isoTWAS association for", gene, "transcript", tx, ":", e$message, "\n")
        return(NULL)
      })
      
      # Write results if successful
      if (!is.null(tx_results) && class(tx_results) == 'data.frame') {
        # Add transcript and R2 information
        tx_results$Transcript <- tx
        tx_results$R2 <- tx_model$R2[1]
        
        # Ensure column names are standardized
        colnames(tx_results) <- c('Gene', 'Transcript', 'Z', 'P', 
                                 'permute.P', 'topSNP', 'topSNP.P', 'R2')
        
        fwrite(
          tx_results,
          out_file,
          append = file.exists(out_file),
          sep = '\t',
          quote = FALSE,
          row.names = FALSE
        )
        cat("isoTWAS results for", gene, "transcript", tx, "written to", out_file, "\n")
      }
    }
    
    # Mark gene as done
    track_gene_done(gene, trait, done_file)
  }
  
  cat("isoTWAS analysis completed for", cancer, "in", tissue, "\n")
}

# Call the main function with placeholder paths
manifest_file <- "PATH_TO_MANIFEST/twas_cancer_tissue_full.tsv"
base_folder <- "PATH_TO_BASE_FOLDER"
gwas_folder <- "PATH_TO_GWAS_FOLDER"
output_folder <- "PATH_TO_OUTPUT_FOLDER"

# Run the isoTWAS analysis
run_isotwas_analysis(manifest_file, index, base_folder, gwas_folder, output_folder)
