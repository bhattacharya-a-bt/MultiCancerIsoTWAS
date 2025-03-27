#########################################################################
# run_sample_isotwas.R
#
# Purpose: Run isoTWAS analysis on the sample data using burdenTest function
#
# Usage: Rscript run_sample_isotwas.R
#########################################################################

# Load required packages
library(data.table)
library(isotwas)
library(boot)
library(yaml)

# Read configuration
config <- yaml::read_yaml("config_example.txt")

# Create directories if they don't exist
for (dir in c(config$results_dir, config$temp_dir, config$output_dir, config$model_dir)) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# Read GWAS summary statistics
cat("Reading GWAS summary statistics...\n")
sumstats <- fread(config$gwas_file)

# Ensure required columns exist
required_cols <- c("SNP", "A1", "A2")
if (!all(required_cols %in% colnames(sumstats))) {
  # Try to map columns if needed
  if ("rsid" %in% colnames(sumstats)) sumstats$SNP <- sumstats$rsid
  
  # Check again
  if (!all(required_cols %in% colnames(sumstats))) {
    stop("Missing required columns in GWAS summary statistics")
  }
}

# Calculate Z-score if needed
if (!("Z" %in% colnames(sumstats))) {
  if (all(c("BETA", "SE") %in% colnames(sumstats))) {
    sumstats$Z <- sumstats$BETA / sumstats$SE
  } else {
    stop("Cannot calculate Z-score: missing BETA or SE columns")
  }
}

# Function to load LD matrix for a specific gene
load_ld_matrix <- function(gene, ld_files) {
  ld_file <- ld_files[grep(gene, ld_files)]
  
  if (length(ld_file) == 0) {
    stop(paste("LD matrix not found for gene", gene))
  }
  
  cat("Loading LD matrix from", ld_file, "\n")
  ld_matrix <- readRDS(ld_file)
  return(ld_matrix)
}

# Process each isoTWAS model
results <- data.frame()

for (model_file in config$isotwas_models) {
  cat("Processing", model_file, "\n")
  
  # Extract gene name
  gene_name <- gsub("_iso.*$", "", basename(model_file))
  
  # Load model
  model <- readRDS(model_file)
  
  # Load LD matrix for this gene
  ld_matrix <- load_ld_matrix(gene_name, config$ld_matrices)
  
  # Get unique transcripts
  transcripts <- unique(model$Feature)
  
  # Process each transcript
  for (tx in transcripts) {
    cat("  Processing transcript", tx, "\n")
    
    # Extract model for this transcript
    tx_model <- subset(model, Feature == tx)
    
    # Skip if model is empty
    if (nrow(tx_model) == 0) {
      cat("  Empty model for transcript", tx, ", skipping\n")
      next
    }
    
    # Run burdenTest function
    result <- tryCatch({
      burdenTest(
        mod = tx_model,
        ld = ld_matrix,
        gene = gene_name,
        sumStats = sumstats,
        chr = 'CHR',
        pos = 'BP',
        a1 = 'A1',
        a2 = 'A2',
        a1_mod = 'ALT',
        a2_mod = 'REF',
        snpName = 'SNP',
        Z = 'Z',
        beta = 'BETA',
        se = 'SE',
        featureName = 'Feature',
        R2cutoff = 0.01,
        alpha = 0.05,  # Use 0.05 as threshold for permutation
        nperms = 1000, # 1000 permutations
        usePos = FALSE
      )
    }, error = function(e) {
      cat("  Error running burdenTest for transcript", tx, ":", e$message, "\n")
      return(NULL)
    })
    
    # Check if result is successful
    if (is.null(result) || !is.data.frame(result)) {
      cat("  Failed to get valid results for transcript", tx, "\n")
      next
    }
    
    # Get chromosome position from model
    chr <- unique(tx_model$Chromosome)
    start <- min(tx_model$Position)
    end <- max(tx_model$Position)
    r2 <- unique(tx_model$R2)
    if (length(r2) > 1) r2 <- r2[1]
    
    # Add to results
    transcript_result <- data.frame(
      File = "sample_trait",
      Gene = gene_name,
      Transcript = result$Feature,
      Z = result$Z,
      P = result$P,
      permute.P = result$permute.P,
      topSNP = result$topSNP,
      topSNP.P = result$topSNP.P,
      R2 = r2,
      Chromosome = chr,
      Start = start,
      End = end
    )
    
    results <- rbind(results, transcript_result)
    cat("  Complete for transcript", tx, "\n")
  }
}

# Write to output file
output_file <- file.path(config$results_dir, "isotwas_results.tsv")
fwrite(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("isoTWAS analysis completed. Results in", output_file, "\n")