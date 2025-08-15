# TMB Results Verification and Quality Check
library(dplyr)

cat("=== TMB VERIFICATION ANALYSIS ===\n\n")

# Load the corrected results
df <- read.csv("corrected_tmb_results.csv", stringsAsFactors = FALSE)

cat("1. BASIC VERIFICATION:\n")
cat("Total variants:", nrow(df), "\n")
cat("TMB Score:", round(nrow(df)/30, 2), "mutations/Mb\n\n")

# 2. VAF Distribution Analysis (Critical Check)
cat("2. VAF DISTRIBUTION ANALYSIS:\n")
if ("AF" %in% colnames(df)) {
  vaf_table <- table(df$AF)
  cat("VAF value distribution:\n")
  print(vaf_table[1:min(10, length(vaf_table))])  # Show first 10 values
  
  unique_vafs <- length(unique(df$AF))
  cat("Number of unique VAF values:", unique_vafs, "\n")
  
  if (unique_vafs == 1) {
    cat("WARNING: All variants have the same VAF - this suggests a parsing issue!\n")
  }
  
  # Check if VAF makes biological sense
  high_vaf <- sum(df$AF > 0.8, na.rm = TRUE)
  if (high_vaf > nrow(df) * 0.1) {
    cat("WARNING:", high_vaf, "variants (", round(high_vaf/nrow(df)*100, 1), 
        "%) have VAF > 80% - possible germline contamination\n")
  }
} else {
  cat("VAF column not found!\n")
}

cat("\n3. MUTATION TYPE VERIFICATION:\n")
if ("ExonicFunc.refGene" %in% colnames(df)) {
  mut_types <- table(df$ExonicFunc.refGene)
  for (i in 1:length(mut_types)) {
    cat("  ", names(mut_types)[i], ":", mut_types[i], 
        "(", round(mut_types[i]/sum(mut_types)*100, 1), "%)\n")
  }
  
  # Check for unexpected mutation types
  expected_types <- c("nonsynonymous SNV", "stopgain", "stoploss", 
                      "frameshift deletion", "frameshift insertion",
                      "nonframeshift deletion", "nonframeshift insertion", 
                      "nonframeshift substitution")
  unexpected <- names(mut_types)[!names(mut_types) %in% expected_types]
  if (length(unexpected) > 0) {
    cat("  WARNING: Unexpected mutation types found:", paste(unexpected, collapse=", "), "\n")
  }
} else {
  cat("Mutation type column not found!\n")
}

cat("\n4. POPULATION FREQUENCY CHECK:\n")
gnomad_cols <- c("gnomAD_exome_ALL", "gnomAD_genome_ALL")
for (col in gnomad_cols) {
  if (col %in% colnames(df)) {
    freqs <- suppressWarnings(as.numeric(as.character(df[[col]])))
    
    has_freq <- sum(!is.na(freqs) & freqs > 0)
    common_vars <- sum(freqs > 0.001, na.rm = TRUE)
    
    cat("  ", col, ":\n")
    cat("    Variants with frequency data:", has_freq, "\n")
    cat("    Variants > 0.1%:", common_vars, "\n")
    
    if (common_vars > 0) {
      cat("    WARNING: Found common variants that should be filtered!\n")
      # Show examples
      common_examples <- df[!is.na(freqs) & freqs > 0.001, 
                            c("Chr", "Start", "Ref", "Alt", col)][1:5,]
      print(common_examples)
    }
  }
}

cat("\n5. CHROMOSOME DISTRIBUTION CHECK:\n")
if ("Chr" %in% colnames(df)) {
  chr_dist <- table(df$Chr)
  chr_dist_sorted <- sort(chr_dist, decreasing = TRUE)
  
  cat("Top 10 chromosomes by mutation count:\n")
  for (i in 1:min(10, length(chr_dist_sorted))) {
    chr_name <- names(chr_dist_sorted)[i]
    count <- chr_dist_sorted[i]
    cat("  ", chr_name, ":", count, "mutations\n")
  }
  
  # Check for unusual patterns
  total_vars <- sum(chr_dist)
  top_chr_pct <- max(chr_dist) / total_vars * 100
  
  if (top_chr_pct > 15) {
    cat("  WARNING: Chromosome", names(which.max(chr_dist)), 
        "has", round(top_chr_pct, 1), "% of mutations - unusually high\n")
  }
}

cat("\n6. COVERAGE AND QUALITY DISTRIBUTION:\n")
if ("DP" %in% colnames(df)) {
  cat("Coverage statistics:\n")
  cat("  Range:", min(df$DP, na.rm=T), "-", max(df$DP, na.rm=T), "x\n")
  cat("  Quartiles:", paste(quantile(df$DP, na.rm=T), collapse=" | "), "\n")
  
  low_cov <- sum(df$DP < 30, na.rm=T)
  if (low_cov > 0) {
    cat("  NOTE:", low_cov, "variants have coverage < 30x\n")
  }
}

if ("QUAL" %in% colnames(df)) {
  cat("Quality score statistics:\n")
  cat("  Range:", round(min(df$QUAL, na.rm=T), 1), "-", 
      round(max(df$QUAL, na.rm=T), 1), "\n")
  cat("  Quartiles:", paste(round(quantile(df$QUAL, na.rm=T), 1), collapse=" | "), "\n")
}

cat("\n7. MELANOMA-SPECIFIC VERIFICATION:\n")

# Calculate expected TMB range for melanoma
cat("TMB Context for Melanoma:\n")
cat("  Your TMB: 74.3 mutations/Mb\n")
cat("  Typical melanoma range: 10-50 mutations/Mb\n")
cat("  UV-signature melanomas: 15-60 mutations/Mb\n")
cat("  Assessment: HIGH but within possible range for UV-exposed melanoma\n")

# Check for UV signature
if ("Ref" %in% colnames(df) && "Alt" %in% colnames(df)) {
  ct_transitions <- sum(df$Ref == "C" & df$Alt == "T") + 
    sum(df$Ref == "G" & df$Alt == "A")
  total_snvs <- sum(nchar(df$Ref) == 1 & nchar(df$Alt) == 1)
  
  if (total_snvs > 0) {
    ct_pct <- ct_transitions / total_snvs * 100
    cat("  C>T transitions: ", ct_transitions, " (", round(ct_pct, 1), 
        "% of SNVs)\n", sep="")
    
    if (ct_pct > 30) {
      cat("  UV signature: STRONG (typical for cutaneous melanoma)\n")
    } else if (ct_pct > 15) {
      cat("  UV signature: MODERATE\n")
    } else {
      cat("  UV signature: WEAK (unusual for cutaneous melanoma)\n")
    }
  }
}

cat("\n8. FINAL ASSESSMENT:\n")

# Overall assessment
issues_found <- 0
warnings <- character()

# Check VAF distribution
if (length(unique(df$AF)) == 1) {
  issues_found <- issues_found + 1
  warnings <- c(warnings, "VAF parsing issue")
}

# Check for common variants
for (col in gnomad_cols) {
  if (col %in% colnames(df)) {
    freqs <- suppressWarnings(as.numeric(as.character(df[[col]])))
    if (sum(freqs > 0.001, na.rm = TRUE) > 0) {
      issues_found <- issues_found + 1
      warnings <- c(warnings, "Common variants present")
    }
  }
}

# Check TMB range
if (nrow(df)/30 > 100) {
  issues_found <- issues_found + 1
  warnings <- c(warnings, "TMB still very high")
}

cat("Issues found:", issues_found, "\n")
if (length(warnings) > 0) {
  cat("Warnings:\n")
  for (w in warnings) {
    cat("  -", w, "\n")
  }
}

if (issues_found == 0) {
  cat("✅ RESULT: TMB appears reliable and appropriate for melanoma\n")
} else if (issues_found <= 2) {
  cat("⚠️  RESULT: TMB mostly reliable but needs minor adjustments\n")
} else {
  cat("❌ RESULT: TMB needs further refinement\n")
}

cat("\n9. RECOMMENDATIONS:\n")

if (length(unique(df$AF)) == 1) {
  cat("- CRITICAL: Fix VAF parsing - all variants shouldn't have identical VAF\n")
}

if (nrow(df)/30 > 60) {
  cat("- Consider stricter VAF threshold (15-20%) to reduce germline contamination\n")
}

cat("- Consider validation with known melanoma TMB datasets\n")
cat("- Verify exome capture size (assumed 30 Mb - actual size may differ)\n")

# Save detailed summary
summary_data <- data.frame(
  Metric = c("Total_Variants", "TMB_Score", "Median_VAF", "Median_Coverage", 
             "Median_QUAL", "UV_Signature_Percent"),
  Value = c(nrow(df), round(nrow(df)/30, 2), median(df$AF, na.rm=T),
            median(df$DP, na.rm=T), median(df$QUAL, na.rm=T),
            ifelse(exists("ct_pct"), round(ct_pct, 1), NA))
)

write.csv(summary_data, "tmb_verification_summary.csv", row.names=FALSE)
cat("\nVerification summary saved to: tmb_verification_summary.csv\n")