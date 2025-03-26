#########################################################################
# run_twas.R
#
# Purpose: Performs transcriptome-wide association studies (TWAS) for
# specified cancer-tissue pairs using precomputed prediction models and
# GWAS summary statistics.
#
# Usage: Rscript run_twas.R --index [INDEX]
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

#' Main function to run TWAS analysis
#'
#' @param manifest_file Path to the manifest file
#' @param index Index of the cancer-tissue pair to analyze
#' @param base_folder Base folder for all data
#' @param gwas_folder Folder containing GWAS summary statistics
#' @param output_folder Folder for output files
run_twas_analysis <- function(manifest_file, index, base_folder, gwas_folder, output_folder) {
  # Read manifest file
  manifest <- as.data.frame(fread(manifest_file, header=TRUE))
  
  # Get cancer type and tissue for this index
  cancer <- manifest$Cancer[index]
  tissue <- manifest$Tissue[index]
  
  # Print what we're analyzing
  cat("Running TWAS for", cancer, "in", tissue, "\n")
  
  # Set paths for this analysis
  twas_folder <- file.path(base_folder, tissue, "TWAS")
  bim_file <- file.path(base_folder, "reference", "reference.bim")
  ld_folder <- file.path(base_folder, "LD_matrices")
  
  # Find GWAS file for this cancer
  gwas_files <- list.files(gwas_folder)
  gwas_files <- gwas_files[grepl('_MAF0.01', gwas_files)]
  trait_names <- sapply(strsplit(gwas_files, '_MAF'), function(x) x[1])
  trait <- trait_names[trait_names == cancer]
  sumstats_file <- file.path(gwas_folder, gwas_files[grep(cancer, gwas_files)])
  
  # Set output paths
  out_file <- file.path(output_folder, trait, tissue, paste0(trait, "_TWAS.tsv"))
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
  genes_to_process <- get_genes_to_process(twas_folder, done_file, trait)
  
  # Process each gene
  for (gene in genes_to_process) {
    # Check if gene has already been processed
    if (is_gene_done(gene, trait, done_file)) {
      cat("Skipping", gene, "- already processed\n")
      next
    }
    
    cat("Processing", gene, "\n")
    
    # Try to load TWAS model
    model_file <- file.path(twas_folder, paste0(gene, "_TWAS.RDS"))
    if (!file.exists(model_file)) {
      cat("No TWAS model found for", gene, "\n")
      track_gene_done(gene, trait, done_file)
      next
    }
    
    # Load and prepare TWAS model
    twas_model <- readRDS(model_file)
    twas_model <- prepare_model(twas_model)
    
    # Check if model has good prediction performance
    if (all(twas_model$R2 <= 0.01)) {
      cat("TWAS model for", gene, "has R² <= 0.01, skipping\n")
      track_gene_done(gene, trait, done_file)
      next
    }
    
    # Load LD matrix for this gene
    ld_matrix <- load_ld_matrix(gene, ld_folder, tissue, type = "twas")
    if (is.null(ld_matrix)) {
      cat("Could not load LD matrix for", gene, "\n")
      next
    }
    
    # Filter TWAS model to only include SNPs in summary statistics
    twas_model <- subset(twas_model, SNP %in% sumstats$SNP)
    
    # Filter summary statistics to only include SNPs in TWAS model
    gene_sumstats <- subset(sumstats, SNP %in% twas_model$SNP)
    
    # Make sure we have position information
    pos_info <- twas_model[, c('SNP', 'Chromosome', 'Position')]
    pos_info <- pos_info[!duplicated(pos_info$SNP), ]
    gene_sumstats <- merge(gene_sumstats, pos_info, by = 'SNP')
    
    # Run TWAS association test
    cat("Running TWAS association for", gene, "\n")
    gene_results <- tryCatch({
      burdenTest(
        mod = twas_model,
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
      cat("Error in TWAS association for", gene, ":", e$message, "\n")
      return(NULL)
    })
    
    # Write results if successful
    if (!is.null(gene_results) && class(gene_results) == 'data.frame') {
      gene_results$R2 <- max(twas_model$R2)
      fwrite(
        gene_results,
        out_file,
        append = file.exists(out_file),
        sep = '\t',
        quote = FALSE,
        row.names = FALSE
      )
      cat("TWAS results for", gene, "written to", out_file, "\n")
    }
    
    # Mark gene as done
    track_gene_done(gene, trait, done_file)
  }
  
  cat("TWAS analysis completed for", cancer, "in", tissue, "\n")
}

# Call the main function with placeholder paths
manifest_file <- "PATH_TO_MANIFEST/twas_cancer_tissue_full.tsv"
base_folder <- "PATH_TO_BASE_FOLDER"
gwas_folder <- "PATH_TO_GWAS_FOLDER"
output_folder <- "PATH_TO_OUTPUT_FOLDER"

# Run the TWAS analysis
run_twas_analysis(manifest_file, index, base_folder, gwas_folder, output_folder)
