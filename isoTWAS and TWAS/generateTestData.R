# This script performs three main tasks:
# 1. Extract SNPs from GWAS file near genes and convert hg19 to hg38
# 2. Copy isoTWAS RDS files for the genes to the specified location
# 3. Extract LD matrix for SNPs within 1Mb of the genes

# Load required libraries
library(data.table)
library(dplyr)
library(rtracklayer)  # For liftOver functionality
library(bigsnpr)      # For LD calculations

# Define file paths
gwas_file <- "/rsrch5/home/epi/bhattacharya_lab/data/munged_GWAS/BreastOverall.tsv.gz"
chain_file <- "/rsrch5/home/epi/bhattacharya_lab/data/GenomicReferences/liftover/hg19ToHg38.over.chain"
output_file <- "/rsrch5/home/epi/abhattacharya3/MultiCancerIsoTWAS/isoTWAS and TWAS/Test Data/SampleGWASData.tsv.gz"

# Define source and destination paths for isoTWAS files
isotwas_source_dir <- "/rsrch5/home/epi/bhattacharya_lab/data/TWAS_weights/GTEx_GENCODEv38/Adipose_Subcutaneous/isoTWAS"
isotwas_dest_dir <- "/rsrch5/home/epi/abhattacharya3/MultiCancerIsoTWAS/isoTWAS and TWAS/Test Data"

# Define LD reference file
ld_ref_file <- "/rsrch5/home/epi/bhattacharya_lab/data/GTEx_v8/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig.LDREF"

# Gene coordinates in hg19 with Ensembl IDs
genes_hg19 <- data.frame(
  gene = c("BABAM1", "KLF2", "BSG"),
  ensembl_id = c("ENSG00000105393", "ENSG00000127528", "ENSG00000172270"),
  chr = c(19, 19, 19),
  start = c(17378252, 16347571, 571513),  # hg19 coordinates
  end = c(17390146, 16352457, 582582)     # hg19 coordinates
)

# Gene coordinates in hg38 with Ensembl IDs
genes_hg38 <- data.frame(
  gene = c("BABAM1", "KLF2", "BSG"),
  ensembl_id = c("ENSG00000105393", "ENSG00000127528", "ENSG00000172270"),
  chr = c(19, 19, 19),
  start = c(17267376, 16327305, 572251),  # hg38 coordinates
  end = c(17281249, 16332304, 583493)     # hg38 coordinates
)

# Function to extract SNPs within range of a gene
extract_snps_near_gene <- function(gwas_data, gene_info, window_size = 1000000) {
  chr_data <- gwas_data[gwas_data$CHR == gene_info$chr, ]
  
  # Define region boundaries with 1 Mb window
  region_start <- max(1, gene_info$start - window_size)
  region_end <- gene_info$end + window_size
  
  # Filter SNPs within the region
  region_snps <- chr_data[chr_data$BP >= region_start & chr_data$BP <= region_end, ]
  
  # Add gene name and Ensembl ID as metadata
  region_snps$NearGene <- gene_info$gene
  region_snps$NearGeneEnsembl <- gene_info$ensembl_id
  
  return(region_snps)
}

# Function to lift over coordinates from hg19 to hg38
lift_over_coordinates <- function(snps_df, chain_file) {
  # Create a GRanges object from the SNPs
  gr <- GRanges(
    seqnames = paste0("chr", snps_df$CHR),
    ranges = IRanges(start = snps_df$BP, end = snps_df$BP),
    strand = "*",
    SNP = snps_df$SNP
  )
  
  # Import the chain file
  chain <- import.chain(chain_file)
  
  # Perform liftOver
  gr_hg38 <- liftOver(gr, chain)
  
  # Convert to data frame and handle potential issues
  gr_hg38_df <- data.frame()
  
  for (i in 1:length(gr_hg38)) {
    if (length(gr_hg38[[i]]) == 1) {
      # Successfully mapped to one location
      snp_entry <- data.frame(
        SNP = mcols(gr_hg38[[i]])$SNP,
        CHR_hg38 = as.character(seqnames(gr_hg38[[i]])),
        BP_hg38 = start(gr_hg38[[i]])
      )
      gr_hg38_df <- rbind(gr_hg38_df, snp_entry)
    } else if (length(gr_hg38[[i]]) > 1) {
      # Mapped to multiple locations - use the first one and add a flag
      snp_entry <- data.frame(
        SNP = mcols(gr_hg38[[i]][1])$SNP,
        CHR_hg38 = as.character(seqnames(gr_hg38[[i]][1])),
        BP_hg38 = start(gr_hg38[[i]][1]),
        Multiple_Mappings = TRUE
      )
      gr_hg38_df <- rbind(gr_hg38_df, snp_entry)
    } else {
      # Not mapped - keep the SNP but add NA values for hg38
      snp_entry <- data.frame(
        SNP = mcols(gr[[i]])$SNP,
        CHR_hg38 = NA,
        BP_hg38 = NA,
        Failed_Mapping = TRUE
      )
      gr_hg38_df <- rbind(gr_hg38_df, snp_entry)
    }
  }
  
  # Clean up chromosomes
  gr_hg38_df$CHR_hg38 <- gsub("chr", "", gr_hg38_df$CHR_hg38)
  
  return(gr_hg38_df)
}

# Function to process the GWAS file and perform liftOver
process_gwas <- function(gwas_file, genes, chain_file, output_file) {
  # Read GWAS summary statistics file
  cat("Reading GWAS summary statistics file...\n")
  gwas_data <- fread(gwas_file)
  
  cat("Extracting SNPs near genes (hg19 coordinates)...\n")
  all_region_snps <- data.frame()
  
  # Process each gene
  for (i in 1:nrow(genes)) {
    gene_info <- genes[i, ]
    cat(paste("Processing gene:", gene_info$gene, "(", gene_info$ensembl_id, ")\n"))
    
    # Extract SNPs near this gene
    region_snps <- extract_snps_near_gene(gwas_data, gene_info)
    
    # Add to combined results
    all_region_snps <- rbind(all_region_snps, region_snps)
  }
  
  # Remove duplicates (SNPs that are near multiple genes)
  all_region_snps <- all_region_snps[!duplicated(all_region_snps$SNP), ]
  
  # Sort by chromosome and position
  all_region_snps <- all_region_snps[order(all_region_snps$CHR, all_region_snps$BP), ]
  
  # Rename hg19 columns for clarity
  names(all_region_snps)[names(all_region_snps) == "CHR"] <- "CHR_hg19"
  names(all_region_snps)[names(all_region_snps) == "BP"] <- "BP_hg19"
  
  # Perform liftOver to hg38
  cat("Converting coordinates from hg19 to hg38...\n")
  hg38_coords <- lift_over_coordinates(all_region_snps, chain_file)
  
  # Merge with original data
  all_region_snps <- merge(all_region_snps, hg38_coords, by = "SNP", all.x = TRUE)
  
  # Sort by hg38 chromosome and position
  all_region_snps <- all_region_snps[order(all_region_snps$CHR_hg38, all_region_snps$BP_hg38), ]
  
  # Prepare output with only the requested columns
  output_data <- all_region_snps[, c("SNP", "CHR_hg38", "BP_hg38", "A1", "A2", "BETA", "SE", "P", "FRQ", "N", "N_CAS", "N_CON", "Z")]
  
  # Rename columns to match requested format
  names(output_data)[names(output_data) == "CHR_hg38"] <- "CHR"
  names(output_data)[names(output_data) == "BP_hg38"] <- "BP"
  
  # Add conversion success statistics
  failed_count <- sum(is.na(all_region_snps$BP_hg38))
  multiple_count <- sum(all_region_snps$Multiple_Mappings, na.rm = TRUE)
  
  cat(paste("\nConversion statistics:", 
            "\nTotal SNPs:", nrow(all_region_snps),
            "\nSuccessfully mapped to hg38:", nrow(all_region_snps) - failed_count,
            "\nFailed to map:", failed_count,
            "\nMapped to multiple locations:", multiple_count, "\n"))
  
  # Write to output file
  cat(paste("Writing", nrow(output_data), "SNPs to output file...\n"))
  fwrite(output_data, output_file, sep = "\t")
  
  cat(paste("Done! Results saved to", output_file, "\n"))
  
  return(output_data)
}

# Function to copy isoTWAS files for the specified genes
copy_isotwas_files <- function(genes, source_dir, dest_dir) {
  cat("\nCopying isoTWAS files...\n")
  
  # Ensure destination directory exists
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  
  # List all files in the source directory
  all_files <- list.files(source_dir, pattern = "_isoTWAS.RDS$", full.names = TRUE)
  
  # Filter files matching our genes' Ensembl IDs
  copied_count <- 0
  
  for (ensembl_id in genes$ensembl_id) {
    # Find files that start with this Ensembl ID
    matching_files <- grep(paste0("^", ensembl_id), basename(all_files), value = TRUE)
    
    # List files in source directory to find matches
    for (file in all_files) {
      if (grepl(paste0("^", ensembl_id), basename(file))) {
        # Get destination path
        dest_file <- file.path(dest_dir, basename(file))
        
        # Copy file
        file.copy(file, dest_file, overwrite = TRUE)
        cat(paste("Copied:", basename(file), "\n"))
        copied_count <- copied_count + 1
      }
    }
  }
  
  if (copied_count == 0) {
    cat("No matching isoTWAS files found for the specified genes.\n")
  } else {
    cat(paste("Copied", copied_count, "isoTWAS files to", dest_dir, "\n"))
  }
}

# Function to extract LD matrix for SNPs within 1Mb of genes
calculate_ld_matrix <- function(genes, ld_ref_file, window_size = 1000000, output_dir) {
  cat("\nCalculating LD matrices...\n")
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Calculate LD for each gene separately
  for (i in 1:nrow(genes)) {
    gene_info <- genes[i, ]
    gene_name <- gene_info$gene
    ensembl_id <- gene_info$ensembl_id
    
    cat(paste("Processing LD for gene:", gene_name, "(", ensembl_id, ")\n"))
    
    # Define region boundaries with 1 Mb window
    chr <- gene_info$chr
    region_start <- max(1, gene_info$start - window_size)
    region_end <- gene_info$end + window_size
    
    # Create temporary directory for the backing files
    temp_dir <- tempdir()
    
    # Extract region using PLINK
    cat(paste("  Extracting region chr", chr, ":", region_start, "-", region_end, "...\n"))
    
    # Prepare temporary file for PLINK
    temp_plink_file <- file.path(temp_dir, paste0(gene_name, "_region"))
    
    # Run PLINK to extract the region
    system(paste0(
      "plink --bfile ", ld_ref_file,
      " --chr ", chr,
      " --from-bp ", region_start,
      " --to-bp ", region_end,
      " --make-bed --out ", temp_plink_file
    ))
    
    # Read the BED file
    cat("  Reading genotype data...\n")
    snp_obj <- snp_attach(snp_readBed2(
      paste0(temp_plink_file, ".bed"),
      backingfile = tempfile(tmpdir = temp_dir)
    ))
    
    # Calculate LD matrix
    cat("  Calculating LD matrix...\n")
    ld_matrix <- snp_cor(snp_obj$genotypes)
    
    # Get SNP information
    snp_info <- snp_obj$map
    
    # Add SNP rsids as row and column names
    rsids <- snp_info$marker.ID
    rownames(ld_matrix) <- rsids
    colnames(ld_matrix) <- rsids
    
    # Save LD matrix with SNP rsids as row and column names
    ld_output_file <- file.path(output_dir, paste0(gene_name, "_LD_matrix.rds"))
    saveRDS(ld_matrix, ld_output_file)
    
    cat(paste("  Saved LD matrix with SNP rsids to", ld_output_file, "\n"))
    
    # Clean up temporary files
    unlink(paste0(temp_plink_file, "*"))
  }
  
  cat("LD matrix calculation completed.\n")
}

# Main execution block

# 1. Process GWAS file
cat("Task 1: Processing GWAS file and converting coordinates\n")
results <- process_gwas(gwas_file, genes_hg19, chain_file, output_file)

# 2. Copy isoTWAS files
cat("\nTask 2: Copying isoTWAS files\n")
copy_isotwas_files(genes_hg19, isotwas_source_dir, isotwas_dest_dir)

# 3. Calculate LD matrices (using hg38 coordinates)
cat("\nTask 3: Calculating LD matrices\n")
calculate_ld_matrix(genes_hg38, ld_ref_file, 1000000, isotwas_dest_dir)

cat("\nAll tasks completed!\n")