# Duplicate Marking Script
# This script performs duplicate marking using Picard MarkDuplicates

# Load required libraries
library(tools)     # For file path operations
library(parallel)  # For parallel processing

# Function to mark duplicates using Picard MarkDuplicates
mark_duplicates <- function(input_bam, output_bam = NULL, metrics_file = NULL, 
                            remove_duplicates = FALSE, temp_dir = tempdir()) {
  
  # Generate output filenames if not provided
  if (is.null(output_bam)) {
    base_name <- tools::file_path_sans_ext(basename(input_bam))
    output_bam <- paste0(tools::file_path_sans_ext(input_bam), "_marked_dups.bam")
  }
  
  if (is.null(metrics_file)) {
    base_name <- tools::file_path_sans_ext(basename(input_bam))
    metrics_file <- paste0(base_name, "_dup_metrics.txt")
  }
  
  # Check if input file exists
  if (!file.exists(input_bam)) {
    stop(paste("Input BAM file not found:", input_bam))
  }
  
  cat("=== Marking Duplicates ===\n")
  cat("Input BAM:", input_bam, "\n")
  cat("Output BAM:", output_bam, "\n")
  cat("Metrics file:", metrics_file, "\n\n")
  
  # Construct MarkDuplicates command
  markdup_cmd <- paste(
    "picard MarkDuplicates",
    paste("INPUT=", input_bam, sep=""),
    paste("OUTPUT=", output_bam, sep=""),
    paste("METRICS_FILE=", metrics_file, sep=""),
    paste("TMP_DIR=", temp_dir, sep=""),
    "VALIDATION_STRINGENCY=LENIENT",
    "CREATE_INDEX=true"
  )
  
  if (remove_duplicates) {
    markdup_cmd <- paste(markdup_cmd, "REMOVE_DUPLICATES=true")
    cat("Note: Duplicates will be REMOVED (not just marked)\n")
  } else {
    cat("Note: Duplicates will be MARKED (not removed)\n")
  }
  
  cat("Running MarkDuplicates command:\n")
  cat(markdup_cmd, "\n\n")
  
  result <- system(markdup_cmd)
  
  if (result == 0) {
    cat("Successfully marked duplicates. Output:", output_bam, "\n")
    
    # Print duplicate statistics if metrics file exists
    if (file.exists(metrics_file)) {
      print_duplicate_stats(metrics_file)
    }
    
    return(list(
      marked_bam = output_bam,
      metrics_file = metrics_file,
      status = "success"
    ))
  } else {
    stop("MarkDuplicates command failed")
  }
}

# Function to print duplicate statistics
print_duplicate_stats <- function(metrics_file) {
  cat("\n=== Duplicate Marking Statistics ===\n")
  
  # Read metrics file and extract key statistics
  tryCatch({
    metrics_lines <- readLines(metrics_file)
    
    # Find the header line to understand column positions
    header_idx <- which(grepl("^LIBRARY", metrics_lines))
    
    if (length(header_idx) > 0) {
      # Get the data line (usually right after header)
      data_idx <- header_idx + 1
      
      if (data_idx <= length(metrics_lines)) {
        header <- strsplit(metrics_lines[header_idx], "\t")[[1]]
        values <- strsplit(metrics_lines[data_idx], "\t")[[1]]
        
        # Find relevant columns
        read_pairs_col <- which(header == "READ_PAIRS_EXAMINED")
        unpaired_col <- which(header == "UNPAIRED_READS_EXAMINED")
        dup_pairs_col <- which(header == "READ_PAIR_DUPLICATES")
        dup_unpaired_col <- which(header == "UNPAIRED_READ_DUPLICATES")
        percent_dup_col <- which(header == "PERCENT_DUPLICATION")
        
        # Calculate and display statistics
        if (length(read_pairs_col) > 0 && length(values) >= read_pairs_col) {
          read_pairs <- as.numeric(values[read_pairs_col])
          unpaired_reads <- as.numeric(values[unpaired_col])
          dup_pairs <- as.numeric(values[dup_pairs_col])
          dup_unpaired <- as.numeric(values[dup_unpaired_col])
          percent_dup <- as.numeric(values[percent_dup_col])
          
          total_reads <- (read_pairs * 2) + unpaired_reads
          total_duplicates <- (dup_pairs * 2) + dup_unpaired
          
          cat("Total read pairs examined:", format(read_pairs, big.mark = ","), "\n")
          cat("Total unpaired reads examined:", format(unpaired_reads, big.mark = ","), "\n")
          cat("Total reads examined:", format(total_reads, big.mark = ","), "\n")
          cat("Duplicate read pairs:", format(dup_pairs, big.mark = ","), "\n")
          cat("Duplicate unpaired reads:", format(dup_unpaired, big.mark = ","), "\n")
          cat("Total duplicate reads:", format(total_duplicates, big.mark = ","), "\n")
          cat("Duplication percentage:", paste0(round(percent_dup * 100, 2), "%"), "\n")
        }
      }
    }
  }, error = function(e) {
    cat("Could not parse duplicate metrics file. Raw content:\n")
    cat("Error:", e$message, "\n")
  })
  
  cat("Full metrics saved to:", metrics_file, "\n\n")
}

# Function to mark duplicates on multiple BAM files
mark_duplicates_multiple <- function(bam_files, output_dir = "marked_duplicates", 
                                     remove_duplicates = FALSE) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results <- list()
  
  for (bam_file in bam_files) {
    base_name <- tools::file_path_sans_ext(basename(bam_file))
    output_bam <- file.path(output_dir, paste0(base_name, "_marked_dups.bam"))
    metrics_file <- file.path(output_dir, paste0(base_name, "_dup_metrics.txt"))
    
    cat("\n=== Processing duplicate marking for:", bam_file, "===\n")
    
    tryCatch({
      result <- mark_duplicates(bam_file, output_bam, metrics_file, remove_duplicates)
      results[[base_name]] <- result
    }, error = function(e) {
      cat("Error processing", bam_file, ":", e$message, "\n")
      results[[base_name]] <- list(status = "error", message = e$message)
    })
  }
  
  return(results)
}

# Function to mark duplicates in parallel
mark_duplicates_parallel <- function(bam_files, output_dir = "marked_duplicates",
                                     remove_duplicates = FALSE,
                                     n_cores = parallel::detectCores() - 1) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Function to process single BAM
  process_single_bam_markdup <- function(bam_file) {
    base_name <- tools::file_path_sans_ext(basename(bam_file))
    output_bam <- file.path(output_dir, paste0(base_name, "_marked_dups.bam"))
    metrics_file <- file.path(output_dir, paste0(base_name, "_dup_metrics.txt"))
    
    tryCatch({
      result <- mark_duplicates(bam_file, output_bam, metrics_file, remove_duplicates)
      return(result)
    }, error = function(e) {
      return(list(status = "error", message = e$message, sample = base_name))
    })
  }
  
  # Run in parallel
  cat("Processing", length(bam_files), "BAM files for duplicate marking using", n_cores, "cores\n")
  results <- mclapply(bam_files, process_single_bam_markdup, mc.cores = n_cores)
  names(results) <- sapply(bam_files, function(x) tools::file_path_sans_ext(basename(x)))
  
  return(results)
}

# Function to generate duplicate marking summary report
generate_duplicate_summary <- function(results) {
  cat("\n=== Duplicate Marking Summary Report ===\n")
  
  successful <- 0
  failed <- 0
  
  for (sample in names(results)) {
    result <- results[[sample]]
    if (result$status == "success") {
      successful <- successful + 1
      cat("✓", sample, "- SUCCESS\n")
    } else {
      failed <- failed + 1
      cat("✗", sample, "- FAILED:", result$message, "\n")
    }
  }
  
  cat("\nSummary:\n")
  cat("Successful:", successful, "\n")
  cat("Failed:", failed, "\n")
  cat("Total:", length(results), "\n")
}

# Example usage function
run_markdup_example <- function() {
  # Example parameters - modify these for your data
  input_bam <- "SRR26456208_aligned_sorted.bam"
  
  cat("=== Single BAM Duplicate Marking Example ===\n")
  result <- mark_duplicates(input_bam)
  print(result)
  
  # Multiple files example
  cat("\n=== Multiple BAM Files Duplicate Marking ===\n")
  bam_files <- c("sample1_aligned_sorted.bam", "sample2_aligned_sorted.bam")
  results <- mark_duplicates_multiple(bam_files)
  
  # Print summary
  generate_duplicate_summary(results)
}

# Function to check if Picard is available
check_picard <- function() {
  result <- system("picard MarkDuplicates --version", intern = TRUE, ignore.stderr = TRUE)
  if (length(result) == 0) {
    # Try alternative command
    result2 <- system("java -jar picard.jar MarkDuplicates --version", intern = TRUE, ignore.stderr = TRUE)
    if (length(result2) == 0) {
      cat("Picard not found in PATH\n")
      cat("Install with: conda install -c bioconda picard\n")
      cat("Or download from: https://broadinstitute.github.io/picard/\n")
      return(FALSE)
    } else {
      cat("Picard is available (via java -jar):\n")
      cat(paste(result2, collapse = "\n"), "\n")
      return(TRUE)
    }
  } else {
    cat("Picard is available:\n")
    cat(paste(result, collapse = "\n"), "\n")
    return(TRUE)
  }
}

# Print usage instructions
cat("=== Duplicate Marking Script ===\n")
cat("Main functions available:\n")
cat("1. mark_duplicates() - Mark duplicates on single BAM file\n")
cat("2. mark_duplicates_multiple() - Process multiple BAM files sequentially\n")
cat("3. mark_duplicates_parallel() - Process multiple BAM files in parallel\n")
cat("4. generate_duplicate_summary() - Generate summary report\n")
cat("5. run_markdup_example() - Run example workflow\n")
cat("6. check_picard() - Verify Picard installation\n\n")

# Automatically check Picard availability
check_picard()

cat("\nRequired files for duplicate marking:\n")
cat("- Input BAM file (aligned and sorted)\n")
cat("- Picard tools must be installed and in PATH\n")
cat("- Sufficient disk space for output BAM and metrics files\n\n")

cat("Usage examples:\n")
cat('mark_duplicates("input.bam")  # Basic usage\n')
cat('mark_duplicates("input.bam", remove_duplicates = TRUE)  # Remove instead of mark\n')
cat('mark_duplicates_multiple(c("file1.bam", "file2.bam"))  # Multiple files\n')