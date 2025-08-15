# Integrated BAM Preprocessing Pipeline
# This script combines duplicate marking, BQSR, and quality assessment
# Optimized for TMB analysis with appropriate filtering strategies

# Load required libraries
library(tools)     # For file path operations
library(parallel)  # For parallel processing
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(reshape2)  # For data transformation
library(gridExtra) # For plot arrangement

# Main preprocessing function
preprocess_bam_for_tmb <- function(input_bam, reference_genome, known_sites_vcf,
                                   output_dir = "preprocessed_bams",
                                   sample_name = NULL,
                                   remove_duplicates = FALSE,
                                   generate_reports = TRUE,
                                   temp_dir = tempdir()) {
  
  # Generate sample name if not provided
  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(input_bam))
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create sample-specific subdirectory
  sample_dir <- file.path(output_dir, sample_name)
  if (!dir.exists(sample_dir)) {
    dir.create(sample_dir, recursive = TRUE)
  }
  
  cat("=== Starting BAM Preprocessing Pipeline ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Input BAM:", input_bam, "\n")
  cat("Output directory:", sample_dir, "\n\n")
  
  # Check input files
  required_files <- c(input_bam, reference_genome)
  for (file in required_files) {
    if (!file.exists(file)) {
      stop(paste("Required file not found:", file))
    }
  }
  
  # Handle compressed VCF files
  if (!file.exists(known_sites_vcf)) {
    if (file.exists(paste0(known_sites_vcf, ".gz"))) {
      known_sites_vcf <- paste0(known_sites_vcf, ".gz")
      cat("Using compressed VCF file:", known_sites_vcf, "\n")
    } else {
      stop(paste("Known sites VCF not found:", known_sites_vcf))
    }
  }
  
  # Define intermediate file names
  marked_bam <- file.path(sample_dir, paste0(sample_name, "_marked_dups.bam"))
  metrics_file <- file.path(sample_dir, paste0(sample_name, "_dup_metrics.txt"))
  recal_table <- file.path(sample_dir, paste0(sample_name, "_recal_data.table"))
  final_bam <- file.path(sample_dir, paste0(sample_name, "_preprocessed.bam"))
  
  results <- list(sample_name = sample_name)
  
  tryCatch({
    # Step 1: Mark Duplicates
    cat("=== Step 1: Marking Duplicates ===\n")
    dup_result <- mark_duplicates_internal(input_bam, marked_bam, metrics_file, 
                                           remove_duplicates, temp_dir)
    results$duplicate_marking <- dup_result
    
    # Step 2: BQSR
    cat("\n=== Step 2: Base Quality Score Recalibration ===\n")
    bqsr_result <- perform_bqsr_internal(marked_bam, reference_genome, known_sites_vcf,
                                         final_bam, recal_table, temp_dir)
    results$bqsr <- bqsr_result
    
    # Step 3: Generate Quality Reports
    if (generate_reports) {
      cat("\n=== Step 3: Generating Quality Assessment Reports ===\n")
      report_result <- generate_preprocessing_reports(input_bam, marked_bam, final_bam,
                                                      metrics_file, recal_table,
                                                      sample_dir, sample_name)
      results$reports <- report_result
    }
    
    results$status <- "success"
    results$final_bam <- final_bam
    
    cat("\n=== Preprocessing Complete ===\n")
    cat("Final preprocessed BAM:", final_bam, "\n")
    
  }, error = function(e) {
    cat("Error in preprocessing:", e$message, "\n")
    results$status <- "error"
    results$error_message <- e$message
  })
  
  return(results)
}

# Internal function for duplicate marking
mark_duplicates_internal <- function(input_bam, output_bam, metrics_file, 
                                     remove_duplicates, temp_dir) {
  
  cat("Input BAM:", input_bam, "\n")
  cat("Output BAM:", output_bam, "\n")
  cat("Metrics file:", metrics_file, "\n")
  
  # Use the detected Picard command
  if (!exists("PICARD_CMD") || is.null(PICARD_CMD)) {
    picard_cmd <- detect_picard_command()
    if (is.null(picard_cmd)) {
      stop("Picard not found. Please install Picard and ensure it's in your PATH.")
    }
    assign("PICARD_CMD", picard_cmd, envir = .GlobalEnv)
  } else {
    picard_cmd <- PICARD_CMD
  }
  
  # Construct MarkDuplicates command
  markdup_cmd <- paste(
    picard_cmd, "MarkDuplicates",
    paste("INPUT=", input_bam, sep=""),
    paste("OUTPUT=", output_bam, sep=""),
    paste("METRICS_FILE=", metrics_file, sep=""),
    paste("TMP_DIR=", temp_dir, sep=""),
    "VALIDATION_STRINGENCY=LENIENT",
    "CREATE_INDEX=true"
  )
  
  if (remove_duplicates) {
    markdup_cmd <- paste(markdup_cmd, "REMOVE_DUPLICATES=true")
    cat("Note: Duplicates will be REMOVED\n")
  } else {
    cat("Note: Duplicates will be MARKED (recommended for TMB analysis)\n")
  }
  
  cat("Running MarkDuplicates...\n")
  result <- system(markdup_cmd)
  
  if (result != 0) {
    stop("MarkDuplicates command failed")
  }
  
  # Parse and return duplicate statistics
  dup_stats <- parse_duplicate_metrics(metrics_file)
  cat("Duplicate marking completed successfully\n")
  
  return(list(
    marked_bam = output_bam,
    metrics_file = metrics_file,
    statistics = dup_stats,
    status = "success"
  ))
}

# Internal function for BQSR
perform_bqsr_internal <- function(input_bam, reference_genome, known_sites_vcf,
                                  output_bam, recal_table, temp_dir) {
  
  cat("Input BAM:", input_bam, "\n")
  cat("Reference genome:", reference_genome, "\n")
  cat("Known sites VCF:", known_sites_vcf, "\n")
  cat("Output BAM:", output_bam, "\n")
  
  # Use the detected GATK command
  if (!exists("GATK_CMD") || is.null(GATK_CMD)) {
    gatk_result <- system("gatk --version", intern = TRUE, ignore.stderr = TRUE)
    if (length(gatk_result) == 0) {
      stop("GATK not found. Please install GATK4 and ensure it's in your PATH.")
    }
    assign("GATK_CMD", "gatk", envir = .GlobalEnv)
  }
  gatk_cmd <- GATK_CMD
  
  # Step 1: BaseRecalibrator
  base_recal_cmd <- paste(
    gatk_cmd, "BaseRecalibrator",
    "-I", input_bam,
    "-R", reference_genome,
    "--known-sites", known_sites_vcf,
    "-O", recal_table,
    paste("--tmp-dir", temp_dir)
  )
  
  cat("Creating recalibration table...\n")
  result1 <- system(base_recal_cmd)
  
  if (result1 != 0) {
    stop("BaseRecalibrator command failed")
  }
  
  # Step 2: ApplyBQSR
  apply_bqsr_cmd <- paste(
    gatk_cmd, "ApplyBQSR",
    "-R", reference_genome,
    "-I", input_bam,
    "--bqsr-recal-file", recal_table,
    "-O", output_bam,
    paste("--tmp-dir", temp_dir)
  )
  
  cat("Applying base quality score recalibration...\n")
  result2 <- system(apply_bqsr_cmd)
  
  if (result2 != 0) {
    stop("ApplyBQSR command failed")
  }
  
  # Index the recalibrated BAM
  cat("Indexing recalibrated BAM...\n")
  if (!exists("SAMTOOLS_CMD") || is.null(SAMTOOLS_CMD)) {
    samtools_result <- system("samtools --version", intern = TRUE, ignore.stderr = TRUE)
    if (length(samtools_result) == 0) {
      warning("samtools not found. BAM indexing may fail.")
      assign("SAMTOOLS_CMD", "samtools", envir = .GlobalEnv)
    } else {
      assign("SAMTOOLS_CMD", "samtools", envir = .GlobalEnv)
    }
  }
  samtools_index_cmd <- paste(SAMTOOLS_CMD, "index", output_bam)
  system(samtools_index_cmd)
  
  cat("BQSR completed successfully\n")
  
  return(list(
    recalibrated_bam = output_bam,
    recal_table = recal_table,
    status = "success"
  ))
}

# Function to parse duplicate metrics
parse_duplicate_metrics <- function(metrics_file) {
  tryCatch({
    metrics_lines <- readLines(metrics_file)
    header_idx <- which(grepl("^LIBRARY", metrics_lines))
    
    if (length(header_idx) > 0) {
      data_idx <- header_idx + 1
      
      if (data_idx <= length(metrics_lines)) {
        header <- strsplit(metrics_lines[header_idx], "\t")[[1]]
        values <- strsplit(metrics_lines[data_idx], "\t")[[1]]
        
        # Extract key statistics
        read_pairs_col <- which(header == "READ_PAIRS_EXAMINED")
        unpaired_col <- which(header == "UNPAIRED_READS_EXAMINED")
        dup_pairs_col <- which(header == "READ_PAIR_DUPLICATES")
        dup_unpaired_col <- which(header == "UNPAIRED_READ_DUPLICATES")
        percent_dup_col <- which(header == "PERCENT_DUPLICATION")
        
        stats <- list(
          read_pairs_examined = as.numeric(values[read_pairs_col]),
          unpaired_reads_examined = as.numeric(values[unpaired_col]),
          duplicate_pairs = as.numeric(values[dup_pairs_col]),
          duplicate_unpaired = as.numeric(values[dup_unpaired_col]),
          percent_duplication = as.numeric(values[percent_dup_col])
        )
        
        # Calculate totals
        stats$total_reads <- (stats$read_pairs_examined * 2) + stats$unpaired_reads_examined
        stats$total_duplicates <- (stats$duplicate_pairs * 2) + stats$duplicate_unpaired
        
        return(stats)
      }
    }
    return(NULL)
  }, error = function(e) {
    cat("Warning: Could not parse duplicate metrics:", e$message, "\n")
    return(NULL)
  })
}

# Function to parse GATK recalibration table
parse_gatk_recal_table <- function(file_path) {
  lines <- readLines(file_path)
  
  # Locate sections
  recal_table1_start <- grep("^#:GATKTable:RecalTable1", lines)
  recal_table2_start <- grep("^#:GATKTable:RecalTable2", lines)
  
  result <- list()
  
  # Parse RecalTable1 (QualityScore table)
  if (length(recal_table1_start) > 0) {
    header_line <- trimws(lines[recal_table1_start + 1])
    header <- strsplit(header_line, "\\s+")[[1]]
    
    data_start <- recal_table1_start + 2
    data_end <- if (length(recal_table2_start) > 0) recal_table2_start - 1 else length(lines)
    table_lines <- lines[data_start:data_end]
    table_lines <- table_lines[!grepl("^#", table_lines) & table_lines != ""]
    
    parsed <- lapply(table_lines, function(x) strsplit(trimws(x), "\\s+")[[1]])
    parsed <- parsed[sapply(parsed, length) == length(header)]
    
    if (length(parsed) > 0) {
      table_data <- do.call(rbind, parsed)
      colnames(table_data) <- header
      df <- as.data.frame(table_data, stringsAsFactors = FALSE)
      
      numeric_cols <- c("QualityScore", "Observations", "Errors")
      for (col in numeric_cols) {
        if (col %in% colnames(df)) {
          df[[col]] <- as.numeric(df[[col]])
        }
      }
      result$QualityScore <- df
    }
  }
  
  return(result)
}

# Function to generate comprehensive preprocessing reports
generate_preprocessing_reports <- function(original_bam, marked_bam, final_bam,
                                           metrics_file, recal_table,
                                           output_dir, sample_name) {
  
  cat("Generating quality assessment reports...\n")
  
  # Generate BAM statistics comparison
  bam_stats <- compare_bam_statistics(original_bam, marked_bam, final_bam, sample_name)
  
  # Generate duplicate marking report
  dup_stats <- parse_duplicate_metrics(metrics_file)
  
  # Generate BQSR plots
  recal_data <- parse_gatk_recal_table(recal_table)
  bqsr_plots <- create_bqsr_plots(recal_data, sample_name)
  
  # Create comprehensive PDF report
  pdf_file <- file.path(output_dir, paste0(sample_name, "_preprocessing_report.pdf"))
  
  pdf(pdf_file, width = 12, height = 8)
  
  # Plot 1: BAM Statistics Overview
  if (!is.null(bam_stats)) {
    print(create_bam_stats_plot(bam_stats, sample_name))
  }
  
  # Plot 2: Duplicate Statistics
  if (!is.null(dup_stats)) {
    print(create_duplicate_stats_plot(dup_stats, sample_name))
  }
  
  # Plot 3-5: BQSR Quality Plots
  if (length(bqsr_plots) > 0) {
    for (plot in bqsr_plots) {
      print(plot)
    }
  }
  
  dev.off()
  
  # Generate text summary report
  summary_file <- file.path(output_dir, paste0(sample_name, "_preprocessing_summary.txt"))
  generate_text_summary(original_bam, final_bam, dup_stats, recal_data, 
                        summary_file, sample_name)
  
  cat("Reports generated:\n")
  cat("- PDF report:", pdf_file, "\n")
  cat("- Summary report:", summary_file, "\n")
  
  return(list(
    pdf_report = pdf_file,
    summary_report = summary_file,
    bam_statistics = bam_stats,
    duplicate_statistics = dup_stats,
    bqsr_data = recal_data
  ))
}

# Function to compare BAM statistics
compare_bam_statistics <- function(original_bam, marked_bam, final_bam, sample_name) {
  
  get_bam_stats <- function(bam_file) {
    if (!exists("SAMTOOLS_CMD")) {
      assign("SAMTOOLS_CMD", "samtools", envir = .GlobalEnv)
    }
    cmd <- paste(SAMTOOLS_CMD, "flagstat", bam_file)
    result <- system(cmd, intern = TRUE)
    
    # Parse flagstat output
    total_reads <- as.numeric(gsub(" .*", "", result[1]))
    mapped_reads <- as.numeric(gsub(" .*", "", result[grep("mapped", result)[1]]))
    
    return(list(
      total_reads = total_reads,
      mapped_reads = mapped_reads,
      mapping_rate = mapped_reads / total_reads
    ))
  }
  
  tryCatch({
    original_stats <- get_bam_stats(original_bam)
    marked_stats <- get_bam_stats(marked_bam)
    final_stats <- get_bam_stats(final_bam)
    
    return(list(
      original = original_stats,
      marked = marked_stats,
      final = final_stats
    ))
  }, error = function(e) {
    cat("Warning: Could not generate BAM statistics:", e$message, "\n")
    return(NULL)
  })
}

# Function to create BQSR plots
create_bqsr_plots <- function(recal_data, sample_name) {
  plots <- list()
  
  if ("QualityScore" %in% names(recal_data)) {
    qs_data <- recal_data$QualityScore
    
    # Calculate empirical quality
    qs_data$EmpiricalQuality <- -10 * log10(qs_data$Errors / qs_data$Observations)
    qs_data$EmpiricalQuality[is.infinite(qs_data$EmpiricalQuality)] <- NA
    
    # Plot 1: Reported vs Empirical Quality
    p1 <- ggplot(qs_data, aes(x = QualityScore, y = EmpiricalQuality)) +
      geom_point(aes(size = Observations), alpha = 0.7, color = "steelblue") +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      labs(title = paste("BQSR Quality Calibration -", sample_name),
           subtitle = "Reported vs Empirical Quality Scores",
           x = "Reported Quality Score",
           y = "Empirical Quality Score",
           size = "Observations") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plots[[1]] <- p1
    
    # Plot 2: Quality Score Distribution
    p2 <- ggplot(qs_data, aes(x = QualityScore, y = Observations)) +
      geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
      labs(title = paste("Quality Score Distribution -", sample_name),
           x = "Quality Score",
           y = "Number of Observations") +
      theme_minimal()
    
    plots[[2]] <- p2
    
    # Plot 3: Error Rate by Quality Score
    qs_data$ErrorRate <- qs_data$Errors / qs_data$Observations
    p3 <- ggplot(qs_data, aes(x = QualityScore, y = ErrorRate)) +
      geom_point(aes(size = Observations), alpha = 0.7, color = "darkred") +
      scale_y_log10() +
      labs(title = paste("Error Rate by Quality Score -", sample_name),
           x = "Quality Score",
           y = "Error Rate (log scale)",
           size = "Observations") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    plots[[3]] <- p3
  }
  
  return(plots)
}

# Function to create duplicate statistics plot
create_duplicate_stats_plot <- function(dup_stats, sample_name) {
  if (is.null(dup_stats)) return(NULL)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Category = c("Unique Reads", "Duplicate Reads"),
    Count = c(dup_stats$total_reads - dup_stats$total_duplicates, dup_stats$total_duplicates),
    Percentage = c((1 - dup_stats$percent_duplication) * 100, dup_stats$percent_duplication * 100)
  )
  
  p <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = paste0(format(Count, big.mark = ","), "\n(", 
                                 round(Percentage, 1), "%)")), 
              vjust = -0.5) +
    scale_fill_manual(values = c("Unique Reads" = "steelblue", 
                                 "Duplicate Reads" = "orange")) +
    labs(title = paste("Duplicate Read Statistics -", sample_name),
         subtitle = paste("Overall Duplication Rate:", 
                          round(dup_stats$percent_duplication * 100, 2), "%"),
         x = "Read Category",
         y = "Number of Reads") +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}

# Function to create BAM statistics plot
create_bam_stats_plot <- function(bam_stats, sample_name) {
  if (is.null(bam_stats)) return(NULL)
  
  # Create comparison data frame
  stages <- c("Original", "Marked Duplicates", "BQSR Recalibrated")
  total_reads <- c(bam_stats$original$total_reads, 
                   bam_stats$marked$total_reads, 
                   bam_stats$final$total_reads)
  mapped_reads <- c(bam_stats$original$mapped_reads,
                    bam_stats$marked$mapped_reads,
                    bam_stats$final$mapped_reads)
  
  plot_data <- data.frame(
    Stage = factor(stages, levels = stages),
    Total_Reads = total_reads,
    Mapped_Reads = mapped_reads,
    Mapping_Rate = mapped_reads / total_reads * 100
  )
  
  p <- ggplot(plot_data, aes(x = Stage)) +
    geom_bar(aes(y = Total_Reads), stat = "identity", fill = "lightblue", alpha = 0.7) +
    geom_bar(aes(y = Mapped_Reads), stat = "identity", fill = "darkblue", alpha = 0.8) +
    geom_text(aes(y = Total_Reads, label = format(Total_Reads, big.mark = ",")), 
              vjust = -0.5, size = 3) +
    labs(title = paste("BAM Processing Pipeline Statistics -", sample_name),
         subtitle = "Blue bars show mapped reads, light blue shows total reads",
         x = "Processing Stage",
         y = "Number of Reads") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Function to generate text summary report
generate_text_summary <- function(original_bam, final_bam, dup_stats, recal_data,
                                  output_file, sample_name) {
  
  sink(output_file)
  
  cat("=== BAM Preprocessing Pipeline Summary Report ===\n")
  cat("Sample:", sample_name, "\n")
  cat("Generated:", Sys.time(), "\n")
  cat("Original BAM:", original_bam, "\n")
  cat("Final BAM:", final_bam, "\n\n")
  
  # Duplicate marking summary
  if (!is.null(dup_stats)) {
    cat("=== Duplicate Marking Results ===\n")
    cat("Total reads examined:", format(dup_stats$total_reads, big.mark = ","), "\n")
    cat("Duplicate reads:", format(dup_stats$total_duplicates, big.mark = ","), "\n")
    cat("Duplication rate:", paste0(round(dup_stats$percent_duplication * 100, 2), "%"), "\n\n")
  }
  
  # BQSR summary
  if ("QualityScore" %in% names(recal_data)) {
    qs_data <- recal_data$QualityScore
    total_obs <- sum(qs_data$Observations, na.rm = TRUE)
    total_errors <- sum(qs_data$Errors, na.rm = TRUE)
    overall_error_rate <- total_errors / total_obs
    
    cat("=== Base Quality Score Recalibration Results ===\n")
    cat("Total observations:", format(total_obs, big.mark = ","), "\n")
    cat("Total errors:", format(total_errors, big.mark = ","), "\n")
    cat("Overall error rate:", sprintf("%.6f", overall_error_rate), "\n")
    cat("Overall empirical quality:", sprintf("%.2f", -10 * log10(overall_error_rate)), "\n\n")
  }
  
  cat("=== Recommendations for TMB Analysis ===\n")
  cat("1. The preprocessed BAM is ready for variant calling\n")
  cat("2. For TMB analysis, use permissive filtering to preserve somatic variants\n")
  cat("3. Focus on non-synonymous mutations in coding regions\n")
  cat("4. Consider tumor-normal paired analysis if normal sample available\n")
  cat("5. Validate TMB results with known cancer gene panels\n")
  
  sink()
}

# Function to process multiple samples in parallel
preprocess_multiple_samples <- function(bam_files, reference_genome, known_sites_vcf,
                                        output_dir = "preprocessed_bams",
                                        remove_duplicates = FALSE,
                                        generate_reports = TRUE,
                                        n_cores = parallel::detectCores() - 1) {
  
  process_single_sample <- function(bam_file) {
    sample_name <- tools::file_path_sans_ext(basename(bam_file))
    
    tryCatch({
      result <- preprocess_bam_for_tmb(bam_file, reference_genome, known_sites_vcf,
                                       output_dir, sample_name, remove_duplicates,
                                       generate_reports)
      return(result)
    }, error = function(e) {
      return(list(
        sample_name = sample_name,
        status = "error",
        error_message = e$message
      ))
    })
  }
  
  cat("Processing", length(bam_files), "samples using", n_cores, "cores\n")
  results <- mclapply(bam_files, process_single_sample, mc.cores = n_cores)
  names(results) <- sapply(bam_files, function(x) tools::file_path_sans_ext(basename(x)))
  
  return(results)
}

# Function to detect and configure Picard installation
detect_picard_command <- function() {
  # Try different Picard command variants
  picard_commands <- c(
    "picard",
    "java -jar picard.jar",
    "java -jar /usr/local/bin/picard.jar", 
    "java -jar ~/picard/picard.jar"
  )
  
  for (cmd in picard_commands) {
    test_cmd <- paste(cmd, "MarkDuplicates --version")
    
    # Test if command works by checking exit status
    exit_status <- system(test_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (exit_status == 0) {
      # Command succeeded, now get the output (which might be in stderr)
      result <- tryCatch({
        # Try to capture output - Picard often outputs version to stderr
        system(test_cmd, intern = TRUE, ignore.stderr = TRUE)
      }, error = function(e) {
        # If intern fails, try without intern and capture manually
        temp_file <- tempfile()
        system(paste(test_cmd, ">", temp_file, "2>&1"))
        output <- readLines(temp_file)
        file.remove(temp_file)
        return(output)
      })
      
      # If we still don't have output, the command worked so it's valid
      if (length(result) == 0 || all(result == "")) {
        cat("Found working Picard command:", cmd, "\n")
        return(cmd)
      }
      
      # Check if we got version info
      if (any(grepl("[0-9]+\\.[0-9]+", result))) {
        cat("Found Picard using command:", cmd, "\n")
        return(cmd)
      }
    }
  }
  
  # Try to find picard.jar in common locations
  common_paths <- c(
    "/usr/local/share/picard/picard.jar",
    "/opt/picard/picard.jar", 
    "/usr/share/picard/picard.jar",
    "~/tools/picard/picard.jar",
    "./picard.jar"
  )
  
  for (path in common_paths) {
    expanded_path <- path.expand(path)
    if (file.exists(expanded_path)) {
      cmd <- paste("java -jar", expanded_path)
      # Test this command
      exit_status <- system(paste(cmd, "MarkDuplicates --version"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      if (exit_status == 0) {
        cat("Found Picard JAR at:", expanded_path, "\n")
        return(cmd)
      }
    }
  }
  
  return(NULL)
}

# Function to check required tools with better error handling
check_required_tools <- function() {
  tools_status <- list()
  tool_commands <- list()
  
  # Check Picard with multiple methods
  cat("=== Checking Picard Installation ===\n")
  
  # First, let's do a simple test that we know works from your terminal
  cat("Testing basic picard command...\n")
  basic_test <- system("picard MarkDuplicates --version", ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  if (basic_test == 0) {
    # The command works, now get version info
    cat("✓ Picard command is working\n")
    tools_status$picard <- TRUE
    tool_commands$picard <- "picard"
    
    # Get version info by redirecting stderr to stdout
    version_cmd <- "picard MarkDuplicates --version 2>&1"
    version_result <- system(version_cmd, intern = TRUE)
    
    if (length(version_result) > 0) {
      # Extract version line (look for pattern like "2.20.4-SNAPSHOT")
      version_lines <- version_result[grepl("[0-9]+\\.[0-9]+", version_result)]
      if (length(version_lines) > 0) {
        cat("✓ Picard version:", version_lines[1], "\n")
      } else {
        cat("✓ Picard is working (version info mixed with warnings)\n")
      }
      
      # Note about locale warning if present
      if (any(grepl("warning.*locale|setlocale", version_result))) {
        cat("📝 Note: Locale warning detected (doesn't affect functionality)\n")
      }
    }
  } else {
    # Basic command failed, try alternative detection
    cat("Basic picard command failed, trying alternative methods...\n")
    picard_cmd <- detect_picard_command()
    if (!is.null(picard_cmd)) {
      tools_status$picard <- TRUE
      tool_commands$picard <- picard_cmd
      cat("✓ Found alternative Picard command:", picard_cmd, "\n")
    } else {
      tools_status$picard <- FALSE
      cat("✗ Picard not found\n")
    }
  }
  
  # Check GATK
  cat("\n=== Checking GATK Installation ===\n")
  gatk_result <- system("gatk --version", intern = TRUE, ignore.stderr = TRUE)
  if (length(gatk_result) > 0 && !grepl("command not found", paste(gatk_result, collapse = " "))) {
    tools_status$gatk <- TRUE
    tool_commands$gatk <- "gatk"
    cat("✓ GATK found:", paste(gatk_result[1:min(2, length(gatk_result))], collapse = " "), "\n")
  } else {
    tools_status$gatk <- FALSE
    cat("✗ GATK not found\n")
  }
  
  # Check samtools
  cat("\n=== Checking samtools Installation ===\n")
  samtools_result <- system("samtools --version", intern = TRUE, ignore.stderr = TRUE)
  if (length(samtools_result) > 0 && !grepl("command not found", paste(samtools_result, collapse = " "))) {
    tools_status$samtools <- TRUE
    tool_commands$samtools <- "samtools"
    cat("✓ samtools found:", paste(samtools_result[1:min(2, length(samtools_result))], collapse = " "), "\n")
  } else {
    tools_status$samtools <- FALSE
    cat("✗ samtools not found\n")
  }
  
  # Summary
  cat("\n=== Tool Status Summary ===\n")
  for (tool in names(tools_status)) {
    status <- if (tools_status[[tool]]) "✓ Available" else "✗ Not found"
    cat(paste(tool, ":", status, "\n"))
  }
  
  if (all(unlist(tools_status))) {
    cat("\n🎉 All required tools are available!\n")
    # Store the working commands globally
    assign("PICARD_CMD", tool_commands$picard, envir = .GlobalEnv)
    assign("GATK_CMD", tool_commands$gatk, envir = .GlobalEnv)
    assign("SAMTOOLS_CMD", tool_commands$samtools, envir = .GlobalEnv)
    return(TRUE)
  } else {
    cat("\n❌ Some tools are missing. Installation instructions:\n")
    
    if (!tools_status$picard) {
      cat("\n--- Picard Installation Options ---\n")
      cat("Option 1 (Conda - Recommended):\n")
      cat("  conda install -c bioconda picard\n")
      cat("Option 2 (Manual download):\n")
      cat("  wget https://github.com/broadinstitute/picard/releases/download/2.27.4/picard.jar\n")
      cat("  # Then use: java -jar picard.jar MarkDuplicates ...\n")
      cat("Option 3 (Ubuntu/Debian):\n")
      cat("  sudo apt-get install picard-tools\n")
    }
    
    if (!tools_status$gatk) {
      cat("\n--- GATK Installation Options ---\n")
      cat("Option 1 (Conda - Recommended):\n")
      cat("  conda install -c bioconda gatk4\n")
      cat("Option 2 (Manual download):\n")
      cat("  wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip\n")
      cat("  unzip gatk-4.4.0.0.zip\n")
      cat("  export PATH=$PATH:/path/to/gatk-4.4.0.0\n")
    }
    
    if (!tools_status$samtools) {
      cat("\n--- samtools Installation Options ---\n")
      cat("Option 1 (Conda - Recommended):\n")
      cat("  conda install -c bioconda samtools\n")
      cat("Option 2 (Ubuntu/Debian):\n")
      cat("  sudo apt-get install samtools\n")
    }
    
    return(FALSE)
  }
}

# Manual tool configuration (if auto-detection fails)
configure_tools_manually <- function(picard_cmd = "picard", gatk_cmd = "gatk", samtools_cmd = "samtools") {
  cat("=== Manual Tool Configuration ===\n")
  cat("Setting tool commands manually...\n")
  
  # Test each command
  cat("Testing Picard:", picard_cmd, "\n")
  picard_test <- system(paste(picard_cmd, "MarkDuplicates --version"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (picard_test == 0) {
    assign("PICARD_CMD", picard_cmd, envir = .GlobalEnv)
    cat("✓ Picard configured successfully\n")
  } else {
    cat("✗ Picard test failed\n")
  }
  
  cat("Testing GATK:", gatk_cmd, "\n")
  gatk_test <- system(paste(gatk_cmd, "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (gatk_test == 0) {
    assign("GATK_CMD", gatk_cmd, envir = .GlobalEnv)
    cat("✓ GATK configured successfully\n")
  } else {
    cat("✗ GATK test failed\n")
  }
  
  cat("Testing samtools:", samtools_cmd, "\n")
  samtools_test <- system(paste(samtools_cmd, "--version"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (samtools_test == 0) {
    assign("SAMTOOLS_CMD", samtools_cmd, envir = .GlobalEnv)
    cat("✓ samtools configured successfully\n")
  } else {
    cat("✗ samtools test failed\n")
  }
  
  cat("\nManual configuration complete. You can now run your pipeline.\n")
}

# Function to help fix locale warning (optional)
fix_locale_warning <- function() {
  cat("=== Fixing Locale Warning ===\n")
  cat("The locale warning doesn't affect Picard functionality, but you can fix it:\n\n")
  
  cat("Option 1 - Set locale in your current session:\n")
  cat('export LC_ALL="C.UTF-8"\n')
  cat('export LANG="C.UTF-8"\n\n')
  
  cat("Option 2 - Add to your ~/.bashrc (permanent fix):\n")
  cat('echo \'export LC_ALL="C.UTF-8"\' >> ~/.bashrc\n')
  cat('echo \'export LANG="C.UTF-8"\' >> ~/.bashrc\n')
  cat('source ~/.bashrc\n\n')
  
  cat("Option 3 - Install additional locales (if you want en_US.UTF-8):\n")
  cat("sudo apt-get install locales\n")
  cat("sudo locale-gen en_US.UTF-8\n")
  cat("sudo update-locale LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8\n\n")
  
  cat("The pipeline will work fine even with the warning!\n")
}

# Quick installation helper function
install_tools_conda <- function() {
  cat("=== Installing Required Tools via Conda ===\n")
  cat("This will install picard, gatk4, and samtools using conda/mamba\n")
  cat("Make sure you have conda/mamba installed first!\n\n")
  
  response <- readline(prompt = "Proceed with installation? (y/n): ")
  
  if (tolower(response) %in% c("y", "yes")) {
    cat("Installing tools...\n")
    
    # Try mamba first (faster), then conda
    if (system("which mamba", ignore.stdout = TRUE) == 0) {
      system("mamba install -c bioconda picard gatk4 samtools -y")
    } else {
      system("conda install -c bioconda picard gatk4 samtools -y")
    }
    
    cat("Installation complete! Please restart R and run check_required_tools() again.\n")
  } else {
    cat("Installation cancelled.\n")
  }
}

# Print usage instructions
cat("=== Integrated BAM Preprocessing Pipeline for TMB Analysis ===\n")
cat("This pipeline combines duplicate marking, BQSR, and quality assessment\n")
cat("Main functions:\n")
cat("1. preprocess_bam_for_tmb() - Complete preprocessing for single sample\n")
cat("2. preprocess_multiple_samples() - Process multiple samples in parallel\n")
cat("3. check_required_tools() - Verify tool installation\n")
cat("4. configure_tools_manually() - Manual tool configuration\n")
cat("5. install_tools_conda() - Quick conda installation helper\n")
cat("6. fix_locale_warning() - Help fix locale warnings\n\n")

# Check tools automatically
cat("Checking tool availability...\n")
tools_available <- check_required_tools()

if (!tools_available) {
  cat("\n💡 If auto-detection failed but tools work in terminal, try:\n")
  cat("configure_tools_manually()\n")
  cat("\n💡 Or install missing tools with: install_tools_conda()\n")
} else {
  cat("\n🚀 Ready to process BAM files! Example usage:\n")
  cat('result <- preprocess_bam_for_tmb(\n')
  cat('  input_bam = "your_sample.bam",\n')
  cat('  reference_genome = "hg38.fa",\n')
  cat('  known_sites_vcf = "dbsnp.vcf.gz"\n')
  cat(')\n')
}

cat("\nExample usage:\n")
cat('# Single sample processing\n')
cat('result <- preprocess_bam_for_tmb(\n')
cat('  input_bam = "sample.bam",\n')
cat('  reference_genome = "hg38.fa",\n')
cat('  known_sites_vcf = "dbsnp.vcf.gz",\n')
cat('  sample_name = "tumor_sample"\n')
cat(')\n\n')

cat('# Multiple samples\n')
cat('results <- preprocess_multiple_samples(\n')
cat('  bam_files = c("sample1.bam", "sample2.bam"),\n')
cat('  reference_genome = "hg38.fa",\n')
cat('  known_sites_vcf = "dbsnp.vcf.gz"\n')
cat(')\n')