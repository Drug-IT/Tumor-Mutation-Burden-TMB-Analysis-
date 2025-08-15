# Base Quality Score Recalibration (BQSR) Script
# This script performs BQSR using GATK BaseRecalibrator and ApplyBQSR

# Load required libraries
library(tools)     # For file path operations
library(parallel)  # For parallel processing

# Function to perform Base Quality Score Recalibration (BQSR)
perform_bqsr <- function(input_bam, reference_genome, known_sites_vcf, 
                         output_bam = NULL, recal_table = NULL, temp_dir = tempdir()) {
  
  # Generate output filenames if not provided
  if (is.null(output_bam)) {
    base_name <- tools::file_path_sans_ext(basename(input_bam))
    output_bam <- gsub("_marked_dups", "_recalibrated", input_bam)
    if (output_bam == input_bam) {  # If no _marked_dups in name
      output_bam <- paste0(tools::file_path_sans_ext(input_bam), "_recalibrated.bam")
    }
  }
  
  if (is.null(recal_table)) {
    base_name <- tools::file_path_sans_ext(basename(input_bam))
    recal_table <- paste0(base_name, "_recal_data.table")
  }
  
  # Check if required files exist
  required_files <- c(input_bam, reference_genome, known_sites_vcf)
  for (file in required_files) {
    if (!file.exists(file)) {
      # If VCF not found, check for compressed version
      if (grepl("\\.vcf$", file) && file.exists(paste0(file, ".gz"))) {
        cat("Found compressed version:", paste0(file, ".gz"), "\n")
        if (file == known_sites_vcf) {
          known_sites_vcf <- paste0(file, ".gz")
          cat("Using compressed VCF file\n")
        }
      } else {
        stop(paste("Required file not found:", file, "\nIf compressed, use .vcf.gz extension"))
      }
    }
  }
  
  cat("=== Step 1: Creating Recalibration Table ===\n")
  
  # Step 1: BaseRecalibrator - Create recalibration table
  base_recal_cmd <- paste(
    "gatk BaseRecalibrator",
    "-I", input_bam,
    "-R", reference_genome,
    "--known-sites", known_sites_vcf,
    "-O", recal_table,
    paste("--tmp-dir", temp_dir)
  )
  
  cat("Running BaseRecalibrator command:\n")
  cat(base_recal_cmd, "\n\n")
  
  result1 <- system(base_recal_cmd)
  
  if (result1 != 0) {
    stop("BaseRecalibrator command failed")
  }
  
  cat("Successfully created recalibration table:", recal_table, "\n\n")
  
  cat("=== Step 2: Applying Base Quality Score Recalibration ===\n")
  
  # Step 2: ApplyBQSR - Apply recalibration
  apply_bqsr_cmd <- paste(
    "gatk ApplyBQSR",
    "-R", reference_genome,
    "-I", input_bam,
    "--bqsr-recal-file", recal_table,
    "-O", output_bam,
    paste("--tmp-dir", temp_dir)
  )
  
  cat("Running ApplyBQSR command:\n")
  cat(apply_bqsr_cmd, "\n\n")
  
  result2 <- system(apply_bqsr_cmd)
  
  if (result2 != 0) {
    stop("ApplyBQSR command failed")
  }
  
  cat("Successfully created recalibrated BAM:", output_bam, "\n")
  
  # Index the recalibrated BAM
  cat("Indexing recalibrated BAM...\n")
  samtools_index_cmd <- paste("samtools", "index", output_bam)
  system(samtools_index_cmd)
  
  return(list(
    recalibrated_bam = output_bam,
    recal_table = recal_table,
    status = "success"
  ))
}

# Function to perform BQSR on multiple BAM files
perform_bqsr_multiple <- function(bam_files, reference_genome, known_sites_vcf, 
                                  output_dir = "bqsr_output") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results <- list()
  
  for (bam_file in bam_files) {
    base_name <- tools::file_path_sans_ext(basename(bam_file))
    output_bam <- file.path(output_dir, paste0(base_name, "_recalibrated.bam"))
    recal_table <- file.path(output_dir, paste0(base_name, "_recal_data.table"))
    
    cat("\n=== Processing BQSR for:", bam_file, "===\n")
    
    tryCatch({
      result <- perform_bqsr(bam_file, reference_genome, known_sites_vcf,
                             output_bam, recal_table)
      results[[base_name]] <- result
    }, error = function(e) {
      cat("Error processing", bam_file, ":", e$message, "\n")
      results[[base_name]] <- list(status = "error", message = e$message)
    })
  }
  
  return(results)
}

# Function to perform BQSR in parallel
perform_bqsr_parallel <- function(bam_files, reference_genome, known_sites_vcf,
                                  output_dir = "bqsr_output", 
                                  n_cores = parallel::detectCores() - 1) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Function to process single BAM
  process_single_bam_bqsr <- function(bam_file) {
    base_name <- tools::file_path_sans_ext(basename(bam_file))
    output_bam <- file.path(output_dir, paste0(base_name, "_recalibrated.bam"))
    recal_table <- file.path(output_dir, paste0(base_name, "_recal_data.table"))
    
    tryCatch({
      result <- perform_bqsr(bam_file, reference_genome, known_sites_vcf,
                             output_bam, recal_table)
      return(result)
    }, error = function(e) {
      return(list(status = "error", message = e$message, sample = base_name))
    })
  }
  
  # Run in parallel
  cat("Processing", length(bam_files), "BAM files for BQSR using", n_cores, "cores\n")
  results <- mclapply(bam_files, process_single_bam_bqsr, mc.cores = n_cores)
  names(results) <- sapply(bam_files, function(x) tools::file_path_sans_ext(basename(x)))
  
  return(results)
}

# Function to generate BQSR quality report (optional)
generate_bqsr_report <- function(input_bam, reference_genome, known_sites_vcf, recal_table, 
                                 output_report = NULL) {
  
  if (is.null(output_report)) {
    base_name <- tools::file_path_sans_ext(basename(input_bam))
    output_report <- paste0(base_name, "_bqsr_report.pdf")
  }
  
  # Generate post-recalibration table for comparison
  post_recal_table <- gsub("\\.table$", "_post.table", recal_table)
  
  post_recal_cmd <- paste(
    "gatk BaseRecalibrator",
    "-I", input_bam,
    "-R", reference_genome,
    "--known-sites", known_sites_vcf,  # Fixed: use actual VCF file
    "-O", post_recal_table
  )
  
  # Generate comparison plots
  plot_cmd <- paste(
    "gatk AnalyzeCovariates",
    "-before", recal_table,
    "-after", post_recal_table,
    "-plots", output_report
  )
  
  cat("Generating BQSR quality report...\n")
  system(post_recal_cmd)
  system(plot_cmd)
  
  cat("BQSR report generated:", output_report, "\n")
  return(output_report)
}

# Example usage function
run_bqsr_example <- function() {
  # Example parameters - modify these for your data
  input_bam <- "SRR26456208_marked_dups.bam"
  reference_genome <- "hg38.fa"
  known_sites_vcf <- "dbsnp.vcf"
  
  cat("=== Single BAM BQSR Example ===\n")
  result <- perform_bqsr(input_bam, reference_genome, known_sites_vcf)
  print(result)
  
  # Multiple files example
  cat("\n=== Multiple BAM Files BQSR ===\n")
  bam_files <- c("sample1_marked_dups.bam", "sample2_marked_dups.bam")
  results <- perform_bqsr_multiple(bam_files, reference_genome, known_sites_vcf)
  
  # Print summary
  cat("\nBQSR Processing Summary:\n")
  for (sample in names(results)) {
    status <- results[[sample]]$status
    cat(paste("Sample:", sample, "- Status:", status, "\n"))
  }
}

# Function to check if GATK is available
check_gatk <- function() {
  result <- system("gatk --version", intern = TRUE, ignore.stderr = TRUE)
  if (length(result) == 0) {
    cat("GATK not found in PATH\n")
    cat("Install with: conda install -c bioconda gatk4\n")
    return(FALSE)
  } else {
    cat("GATK is available:\n")
    cat(paste(result, collapse = "\n"), "\n")
    return(TRUE)
  }
}

# Print usage instructions
cat("=== Base Quality Score Recalibration (BQSR) Script ===\n")
cat("Main functions available:\n")
cat("1. perform_bqsr() - Perform BQSR on single BAM file\n")
cat("2. perform_bqsr_multiple() - Process multiple BAM files sequentially\n")
cat("3. perform_bqsr_parallel() - Process multiple BAM files in parallel\n")
cat("4. generate_bqsr_report() - Generate quality assessment report\n")
cat("5. run_bqsr_example() - Run example workflow\n")
cat("6. check_gatk() - Verify GATK installation\n\n")

# Automatically check GATK availability
check_gatk()

cat("\nRequired files for BQSR:\n")
cat("- Input BAM file (preferably after duplicate marking)\n")
cat("- Reference genome (FASTA format, indexed)\n")
cat("- Known sites VCF (e.g., dbSNP, 1000 Genomes, Mills indels)\n")
cat("- GATK4 must be installed and in PATH\n")