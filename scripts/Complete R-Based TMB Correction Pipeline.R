# Complete R-Based TMB Correction Pipeline
# This script will extract quality metrics from your VCF and apply proper filtering

library(dplyr)
library(readr)

cat("=== R-BASED TMB CORRECTION PIPELINE ===\n\n")

# Function to parse VCF and extract quality metrics
parse_vcf_quality <- function(vcf_file) {
  cat("Step 1: Reading VCF file and extracting quality metrics...\n")
  
  # Read VCF file (skip header lines starting with #)
  vcf_lines <- readLines(vcf_file)
  
  # Find where actual data starts (after header)
  data_start <- which(!grepl("^#", vcf_lines))[1]
  vcf_data_lines <- vcf_lines[data_start:length(vcf_lines)]
  
  # Parse each line
  quality_data <- data.frame(
    Chr = character(),
    Pos = numeric(),
    Ref = character(),
    Alt = character(),
    QUAL = numeric(),
    DP = numeric(),
    AF = numeric(),
    AO = numeric(),
    RO = numeric(),
    AB = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(vcf_data_lines)) {
    if (i %% 10000 == 0) cat("  Processed", i, "variants...\n")
    
    line <- vcf_data_lines[i]
    fields <- strsplit(line, "\t")[[1]]
    
    if (length(fields) >= 8) {
      chr <- fields[1]
      pos <- as.numeric(fields[2])
      ref <- fields[4]
      alt <- fields[5]
      qual <- as.numeric(fields[6])
      info <- fields[8]
      
      # Extract metrics from INFO field
      dp <- NA
      af <- NA
      ao <- NA
      ro <- NA
      ab <- NA
      
      # Parse INFO field
      info_parts <- strsplit(info, ";")[[1]]
      for (part in info_parts) {
        if (grepl("^DP=", part)) dp <- as.numeric(sub("DP=", "", part))
        if (grepl("^AF=", part)) af <- as.numeric(sub("AF=", "", part))
        if (grepl("^AO=", part)) ao <- as.numeric(sub("AO=", "", part))
        if (grepl("^RO=", part)) ro <- as.numeric(sub("RO=", "", part))
        if (grepl("^AB=", part)) ab <- as.numeric(sub("AB=", "", part))
      }
      
      # Add to data frame
      quality_data <- rbind(quality_data, data.frame(
        Chr = chr,
        Pos = pos,
        Ref = ref,
        Alt = alt,
        QUAL = qual,
        DP = dp,
        AF = af,
        AO = ao,
        RO = ro,
        AB = ab,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat("  Extracted quality data for", nrow(quality_data), "variants\n\n")
  return(quality_data)
}

# Function to merge quality data with ANNOVAR results
merge_quality_with_annovar <- function(annovar_file, quality_data) {
  cat("Step 2: Reading ANNOVAR results and merging with quality data...\n")
  
  # Read ANNOVAR results
  annovar_data <- read.csv(annovar_file, stringsAsFactors = FALSE)
  cat("  ANNOVAR data:", nrow(annovar_data), "variants\n")
  
  # Create matching keys
  quality_data$key <- paste(quality_data$Chr, quality_data$Pos, 
                            quality_data$Ref, quality_data$Alt, sep = "_")
  annovar_data$key <- paste(annovar_data$Chr, annovar_data$Start, 
                            annovar_data$Ref, annovar_data$Alt, sep = "_")
  
  # Merge datasets
  merged_data <- merge(annovar_data, 
                       quality_data[, c("key", "QUAL", "DP", "AF", "AO", "RO", "AB")],
                       by = "key", all.x = TRUE)
  
  cat("  Merged data:", nrow(merged_data), "variants\n")
  cat("  Variants with quality data:", sum(!is.na(merged_data$DP)), "\n\n")
  
  return(merged_data)
}

# Function to apply comprehensive filtering
apply_tmb_filters <- function(merged_data, 
                              min_vaf = 0.10, 
                              max_vaf = 0.90,
                              min_coverage = 20,
                              min_alt_reads = 4,
                              min_qual = 30,
                              max_pop_freq = 0.001) {
  
  cat("Step 3: Applying TMB filtering criteria...\n")
  cat("  Filters:\n")
  cat("    VAF range:", min_vaf, "-", max_vaf, "\n")
  cat("    Min coverage:", min_coverage, "x\n")
  cat("    Min alternate reads:", min_alt_reads, "\n")
  cat("    Min QUAL score:", min_qual, "\n")
  cat("    Max population frequency:", max_pop_freq, "\n\n")
  
  # Start with all data
  filtered_data <- merged_data
  
  # Track filtering steps
  filter_stats <- data.frame(
    Step = character(),
    Variants_Remaining = numeric(),
    Variants_Removed = numeric(),
    stringsAsFactors = FALSE
  )
  
  initial_count <- nrow(filtered_data)
  filter_stats <- rbind(filter_stats, data.frame(
    Step = "Initial",
    Variants_Remaining = initial_count,
    Variants_Removed = 0
  ))
  
  # Filter 1: Functional mutations only
  functional_mutations <- c("nonsynonymous SNV", "stopgain", "stoploss", 
                            "frameshift deletion", "frameshift insertion",
                            "nonframeshift deletion", "nonframeshift insertion", 
                            "nonframeshift substitution")
  
  before_count <- nrow(filtered_data)
  filtered_data <- filtered_data[filtered_data$ExonicFunc.refGene %in% functional_mutations, ]
  after_count <- nrow(filtered_data)
  
  filter_stats <- rbind(filter_stats, data.frame(
    Step = "Functional mutations only",
    Variants_Remaining = after_count,
    Variants_Removed = before_count - after_count
  ))
  
  # Filter 2: Quality score
  if (!all(is.na(filtered_data$QUAL))) {
    before_count <- nrow(filtered_data)
    filtered_data <- filtered_data[!is.na(filtered_data$QUAL) & 
                                     filtered_data$QUAL >= min_qual, ]
    after_count <- nrow(filtered_data)
    
    filter_stats <- rbind(filter_stats, data.frame(
      Step = paste("QUAL >=", min_qual),
      Variants_Remaining = after_count,
      Variants_Removed = before_count - after_count
    ))
  }
  
  # Filter 3: Coverage
  if (!all(is.na(filtered_data$DP))) {
    before_count <- nrow(filtered_data)
    filtered_data <- filtered_data[!is.na(filtered_data$DP) & 
                                     filtered_data$DP >= min_coverage, ]
    after_count <- nrow(filtered_data)
    
    filter_stats <- rbind(filter_stats, data.frame(
      Step = paste("Coverage >=", min_coverage),
      Variants_Remaining = after_count,
      Variants_Removed = before_count - after_count
    ))
  }
  
  # Filter 4: VAF range
  if (!all(is.na(filtered_data$AF))) {
    before_count <- nrow(filtered_data)
    filtered_data <- filtered_data[!is.na(filtered_data$AF) & 
                                     filtered_data$AF >= min_vaf & 
                                     filtered_data$AF <= max_vaf, ]
    after_count <- nrow(filtered_data)
    
    filter_stats <- rbind(filter_stats, data.frame(
      Step = paste("VAF", min_vaf, "-", max_vaf),
      Variants_Remaining = after_count,
      Variants_Removed = before_count - after_count
    ))
  }
  
  # Filter 5: Minimum alternate reads
  if (!all(is.na(filtered_data$AO))) {
    before_count <- nrow(filtered_data)
    filtered_data <- filtered_data[!is.na(filtered_data$AO) & 
                                     filtered_data$AO >= min_alt_reads, ]
    after_count <- nrow(filtered_data)
    
    filter_stats <- rbind(filter_stats, data.frame(
      Step = paste("Alt reads >=", min_alt_reads),
      Variants_Remaining = after_count,
      Variants_Removed = before_count - after_count
    ))
  }
  
  # Filter 6: Population frequency (gnomAD)
  gnomad_cols <- c("gnomAD_exome_ALL", "gnomAD_genome_ALL")
  for (col in gnomad_cols) {
    if (col %in% colnames(filtered_data)) {
      before_count <- nrow(filtered_data)
      
      # Keep variants with no frequency data OR frequency <= threshold
      keep_variants <- is.na(filtered_data[[col]]) | 
        filtered_data[[col]] == "." | 
        filtered_data[[col]] == "" |
        (suppressWarnings(as.numeric(filtered_data[[col]])) <= max_pop_freq)
      
      filtered_data <- filtered_data[keep_variants, ]
      after_count <- nrow(filtered_data)
      
      filter_stats <- rbind(filter_stats, data.frame(
        Step = paste(col, "<=", max_pop_freq),
        Variants_Remaining = after_count,
        Variants_Removed = before_count - after_count
      ))
    }
  }
  
  cat("Filtering steps completed:\n")
  print(filter_stats)
  cat("\n")
  
  return(list(filtered_data = filtered_data, filter_stats = filter_stats))
}

# Function to calculate and report TMB
calculate_tmb <- function(filtered_data, exome_size_mb = 30) {
  cat("Step 4: Calculating TMB...\n")
  
  n_mutations <- nrow(filtered_data)
  tmb_score <- n_mutations / exome_size_mb
  
  cat("Final results:\n")
  cat("  Total somatic mutations:", n_mutations, "\n")
  cat("  Exome size:", exome_size_mb, "Mb\n")
  cat("  TMB Score:", round(tmb_score, 2), "mutations/Mb\n")
  
  # TMB classification
  if (tmb_score < 6) {
    classification <- "Low TMB"
  } else if (tmb_score < 20) {
    classification <- "Intermediate TMB"
  } else {
    classification <- "High TMB"
  }
  
  cat("  TMB Classification:", classification, "\n\n")
  
  # Quality summary
  if (!all(is.na(filtered_data$AF))) {
    cat("Quality summary of final variants:\n")
    cat("  VAF range:", round(min(filtered_data$AF, na.rm = TRUE), 3), "-", 
        round(max(filtered_data$AF, na.rm = TRUE), 3), "\n")
    cat("  Median VAF:", round(median(filtered_data$AF, na.rm = TRUE), 3), "\n")
  }
  
  if (!all(is.na(filtered_data$DP))) {
    cat("  Coverage range:", min(filtered_data$DP, na.rm = TRUE), "-", 
        max(filtered_data$DP, na.rm = TRUE), "x\n")
    cat("  Median coverage:", median(filtered_data$DP, na.rm = TRUE), "x\n")
  }
  
  if (!all(is.na(filtered_data$QUAL))) {
    cat("  QUAL range:", round(min(filtered_data$QUAL, na.rm = TRUE), 1), "-", 
        round(max(filtered_data$QUAL, na.rm = TRUE), 1), "\n")
    cat("  Median QUAL:", round(median(filtered_data$QUAL, na.rm = TRUE), 1), "\n")
  }
  
  return(list(
    tmb_score = tmb_score,
    classification = classification,
    n_mutations = n_mutations
  ))
}

# MAIN EXECUTION
cat("Starting TMB correction pipeline...\n\n")

# Step 1: Extract quality metrics from VCF
vcf_file <- "SRR26456208_normalized.vcf"  # Adjust path if needed
quality_data <- parse_vcf_quality(vcf_file)

# Step 2: Merge with ANNOVAR results
annovar_file <- "filtered_somatic_mutations.csv"  # Adjust path if needed
merged_data <- merge_quality_with_annovar(annovar_file, quality_data)

# Step 3: Apply comprehensive filtering
filtering_result <- apply_tmb_filters(merged_data)
filtered_data <- filtering_result$filtered_data

# Step 4: Calculate corrected TMB
tmb_result <- calculate_tmb(filtered_data)

# Step 5: Save results
output_file <- "corrected_tmb_results.csv"
write.csv(filtered_data, output_file, row.names = FALSE)
cat("Corrected TMB results saved to:", output_file, "\n\n")

# Step 6: Comparison summary
cat("=== COMPARISON SUMMARY ===\n")
cat("Original TMB: 793.63 mutations/Mb\n")
cat("Corrected TMB:", round(tmb_result$tmb_score, 2), "mutations/Mb\n")
cat("Reduction factor:", round(793.63 / tmb_result$tmb_score, 1), "x\n")
cat("Original variants: 23,809\n")
cat("Corrected variants:", tmb_result$n_mutations, "\n")
cat("Retention rate:", round(tmb_result$n_mutations / 23809 * 100, 1), "%\n")

cat("\n=== CLINICAL INTERPRETATION ===\n")
cat("TMB Classification:", tmb_result$classification, "\n")

if (tmb_result$tmb_score >= 20) {
  cat("Clinical recommendation: Excellent immunotherapy candidate\n")
} else if (tmb_result$tmb_score >= 6) {
  cat("Clinical recommendation: Good immunotherapy candidate\n")
} else {
  cat("Clinical recommendation: Limited immunotherapy benefit expected\n")
}

cat("\nPipeline completed successfully!\n")