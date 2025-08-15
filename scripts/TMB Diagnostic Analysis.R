# TMB Diagnostic Analysis Script
# This will help identify why your TMB is so high

# Load your final filtered mutations
# Replace with your actual file path
df <- read.csv("filtered_somatic_mutations.csv", stringsAsFactors = FALSE)

# Alternative: if the file has a different name, try:
# df <- read.csv("your_actual_filename.csv", stringsAsFactors = FALSE)

cat("=== TMB DIAGNOSTIC ANALYSIS ===\n\n")

# 1. Basic statistics
cat("1. BASIC STATISTICS:\n")
cat("Total mutations:", nrow(df), "\n")
cat("Reported TMB:", nrow(df) / 30, "mutations/Mb\n\n")

# 2. VAF distribution analysis
cat("2. VAF DISTRIBUTION ANALYSIS:\n")
if ("VAF" %in% colnames(df) || "AF" %in% colnames(df)) {
  vaf_col <- ifelse("VAF" %in% colnames(df), "VAF", "AF")
  vaf_values <- as.numeric(df[[vaf_col]])
  
  cat("VAF statistics:\n")
  cat("  Min VAF:", min(vaf_values, na.rm = TRUE), "\n")
  cat("  Median VAF:", median(vaf_values, na.rm = TRUE), "\n")
  cat("  Mean VAF:", mean(vaf_values, na.rm = TRUE), "\n")
  cat("  Max VAF:", max(vaf_values, na.rm = TRUE), "\n")
  
  # Count low VAF variants (likely artifacts)
  low_vaf <- sum(vaf_values < 0.1, na.rm = TRUE)
  cat("  Variants with VAF < 10%:", low_vaf, "(", round(low_vaf/nrow(df)*100, 1), "%)\n")
  
  very_low_vaf <- sum(vaf_values < 0.05, na.rm = TRUE)
  cat("  Variants with VAF < 5%:", very_low_vaf, "(", round(very_low_vaf/nrow(df)*100, 1), "%)\n")
} else {
  cat("  VAF column not found - this is a problem for quality assessment!\n")
}

# 3. Coverage analysis
cat("\n3. COVERAGE ANALYSIS:\n")
coverage_cols <- c("DP", "Coverage", "Depth", "depth", "dp")
coverage_col <- coverage_cols[coverage_cols %in% colnames(df)][1]

if (!is.na(coverage_col)) {
  coverage_values <- as.numeric(df[[coverage_col]])
  cat("Coverage statistics:\n")
  cat("  Min coverage:", min(coverage_values, na.rm = TRUE), "\n")
  cat("  Median coverage:", median(coverage_values, na.rm = TRUE), "\n")
  cat("  Mean coverage:", mean(coverage_values, na.rm = TRUE), "\n")
  
  low_cov <- sum(coverage_values < 20, na.rm = TRUE)
  cat("  Variants with coverage < 20x:", low_cov, "(", round(low_cov/nrow(df)*100, 1), "%)\n")
} else {
  cat("  Coverage column not found!\n")
}

# 4. Mutation type distribution
cat("\n4. MUTATION TYPE DISTRIBUTION:\n")
if ("ExonicFunc.refGene" %in% colnames(df)) {
  mut_types <- table(df$ExonicFunc.refGene)
  for (i in 1:length(mut_types)) {
    cat("  ", names(mut_types)[i], ":", mut_types[i], "(", round(mut_types[i]/sum(mut_types)*100, 1), "%)\n")
  }
  
  # Check for synonymous mutations (shouldn't be in TMB)
  if ("synonymous SNV" %in% names(mut_types)) {
    cat("  WARNING: Synonymous mutations found - these should NOT be in TMB calculation!\n")
  }
}

# 5. Population frequency check
cat("\n5. POPULATION FREQUENCY CHECK:\n")
gnomad_cols <- c("gnomAD_exome_ALL", "gnomAD_genome_ALL")
for (col in gnomad_cols) {
  if (col %in% colnames(df)) {
    freqs <- suppressWarnings(as.numeric(as.character(df[[col]])))
    
    # Count variants with different frequency ranges
    common_vars <- sum(freqs > 0.01, na.rm = TRUE)  # > 1%
    uncommon_vars <- sum(freqs > 0.001 & freqs <= 0.01, na.rm = TRUE)  # 0.1-1%
    rare_vars <- sum(freqs <= 0.001, na.rm = TRUE)  # <= 0.1%
    no_freq <- sum(is.na(freqs))
    
    cat("  ", col, ":\n")
    cat("    Common variants (>1%):", common_vars, "\n")
    cat("    Uncommon variants (0.1-1%):", uncommon_vars, "\n") 
    cat("    Rare variants (<=0.1%):", rare_vars, "\n")
    cat("    No frequency data:", no_freq, "\n")
    
    if (common_vars > 0) {
      cat("    WARNING: Found", common_vars, "common variants that should be filtered!\n")
    }
  }
}

# 6. Chromosome distribution (to check for potential artifacts)
cat("\n6. CHROMOSOME DISTRIBUTION:\n")
if ("Chr" %in% colnames(df)) {
  chr_dist <- table(df$Chr)
  # Show top 10 chromosomes by mutation count
  top_chrs <- head(sort(chr_dist, decreasing = TRUE), 10)
  for (i in 1:length(top_chrs)) {
    cat("  Chr", names(top_chrs)[i], ":", top_chrs[i], "mutations\n")
  }
}

cat("\n=== RECOMMENDATIONS ===\n")
cat("Based on this analysis:\n")
cat("1. If you have many low VAF (<10%) variants, increase VAF threshold\n")
cat("2. If you have common population variants, check your filtering\n")
cat("3. If you have synonymous mutations, exclude them from TMB\n")
cat("4. Consider using stricter coverage thresholds (>20x)\n")
cat("5. If possible, obtain a matched normal sample for better filtering\n")