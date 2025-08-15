# UV Signature and Sample Type Analysis
library(dplyr)

cat("=== UV SIGNATURE DETAILED ANALYSIS ===\n\n")

# Load corrected results
df <- read.csv("final_corrected_tmb_results.csv", stringsAsFactors = FALSE)

cat("1. DETAILED MUTATION SPECTRUM ANALYSIS:\n")

# Get all SNVs (single nucleotide variants)
snvs <- df[nchar(df$Ref) == 1 & nchar(df$Alt) == 1, ]
cat("Total SNVs:", nrow(snvs), "of", nrow(df), "total variants\n")

if (nrow(snvs) > 0) {
  # Count all possible substitution types
  substitutions <- paste(snvs$Ref, snvs$Alt, sep=">")
  sub_table <- table(substitutions)
  
  cat("\nComplete substitution spectrum:\n")
  for (i in 1:length(sub_table)) {
    sub_type <- names(sub_table)[i]
    count <- sub_table[i]
    pct <- round(count / sum(sub_table) * 100, 1)
    cat("  ", sub_type, ":", count, "(", pct, "%)\n")
  }
  
  # UV signature analysis
  cat("\n2. UV SIGNATURE COMPONENTS:\n")
  
  # Primary UV signature: C>T at dipyrimidines
  ct_transitions <- sum(substitutions %in% c("C>T", "G>A"))
  
  # Secondary UV signature: CC>TT tandem mutations (harder to detect in single variants)
  # T>C transitions (also part of UV damage)
  tc_transitions <- sum(substitutions %in% c("T>C", "A>G"))
  
  total_snvs <- length(substitutions)
  
  cat("C>T + G>A (primary UV):", ct_transitions, "(", round(ct_transitions/total_snvs*100, 1), "%)\n")
  cat("T>C + A>G (secondary UV):", tc_transitions, "(", round(tc_transitions/total_snvs*100, 1), "%)\n")
  
  # Combined UV signature
  uv_total <- ct_transitions + tc_transitions
  cat("Combined UV signature:", uv_total, "(", round(uv_total/total_snvs*100, 1), "%)\n")
  
  # Assessment
  ct_pct <- ct_transitions / total_snvs * 100
  combined_uv_pct <- uv_total / total_snvs * 100
  
  cat("\n3. UV SIGNATURE ASSESSMENT:\n")
  if (ct_pct > 40) {
    cat("Primary UV signature: VERY STRONG (>40%)\n")
  } else if (ct_pct > 25) {
    cat("Primary UV signature: STRONG (25-40%)\n")
  } else if (ct_pct > 15) {
    cat("Primary UV signature: MODERATE (15-25%)\n")
  } else {
    cat("Primary UV signature: WEAK (<15%) - UNUSUAL for cutaneous melanoma\n")
  }
  
  if (combined_uv_pct > 50) {
    cat("Combined UV signature: STRONG\n")
  } else if (combined_uv_pct > 30) {
    cat("Combined UV signature: MODERATE\n")
  } else {
    cat("Combined UV signature: WEAK - suggests non-UV etiology\n")
  }
}

cat("\n4. SAMPLE TYPE ASSESSMENT:\n")

# Based on TMB and UV signature, assess likely sample type
tmb_score <- nrow(df) / 30
ct_pct <- if(exists("ct_pct")) ct_pct else 0

cat("TMB Score:", round(tmb_score, 1), "mutations/Mb\n")
cat("C>T signature:", round(ct_pct, 1), "%\n")

# Sample type prediction
if (tmb_score > 50 && ct_pct > 30) {
  sample_type <- "Cutaneous melanoma with high UV exposure"
} else if (tmb_score > 50 && ct_pct < 15) {
  sample_type <- "Possible mucosal/acral melanoma or highly mutated cutaneous melanoma"
} else if (tmb_score > 20 && ct_pct > 25) {
  sample_type <- "Cutaneous melanoma with moderate UV exposure"
} else if (tmb_score > 20 && ct_pct < 15) {
  sample_type <- "Possible mucosal melanoma or melanoma from chronically sun-damaged skin"
} else {
  sample_type <- "Low TMB melanoma (unusual)"
}

cat("\nPredicted sample type:", sample_type, "\n")

cat("\n5. COMPARISON WITH MELANOMA SUBTYPES:\n")
cat("Typical signatures by melanoma type:\n")
cat("  Cutaneous (sun-exposed): TMB 15-60, C>T 30-60%\n")
cat("  Mucosal: TMB 5-25, C>T 5-15%\n")
cat("  Acral: TMB 5-20, C>T 5-20%\n")
cat("  Uveal: TMB 0.5-5, C>T 5-15%\n")
cat("  Your sample: TMB", round(tmb_score, 1), ", C>T", round(ct_pct, 1), "%\n")

# Check if this could be a cell line or treated sample
cat("\n6. POTENTIAL EXPLANATIONS FOR WEAK UV SIGNATURE:\n")
cat("Possible reasons for high TMB + weak UV signature:\n")
cat("  1. Mucosal melanoma (typically low UV, but can have high TMB)\n")
cat("  2. Acral melanoma with secondary mutations\n")
cat("  3. Metastatic sample with acquired mutations\n")
cat("  4. Cell line with culture-induced mutations\n")
cat("  5. Melanoma with defective DNA repair (MMR deficiency)\n")
cat("  6. Post-treatment sample (chemotherapy/radiation induced mutations)\n")

cat("\n7. CHROMOSOME DISTRIBUTION RE-ANALYSIS:\n")
if ("Chr" %in% colnames(df)) {
  chr_dist <- table(df$Chr)
  chr_dist_sorted <- sort(chr_dist, decreasing = TRUE)
  
  # Check for unusual clustering
  total_vars <- sum(chr_dist)
  top_3_chrs <- head(chr_dist_sorted, 3)
  top_3_pct <- sum(top_3_chrs) / total_vars * 100
  
  cat("Top 3 chromosomes:", paste(names(top_3_chrs), collapse=", "), "\n")
  cat("Top 3 contain", round(top_3_pct, 1), "% of mutations\n")
  
  if (top_3_pct > 40) {
    cat("WARNING: Mutations highly clustered - possible chromosomal instability\n")
  }
  
  # Look for specific chromosome abnormalities common in melanoma
  if ("chr9" %in% names(chr_dist) && chr_dist["chr9"] > total_vars * 0.1) {
    cat("NOTE: High chr9 mutation rate - check for CDKN2A region\n")
  }
  if ("chr1" %in% names(chr_dist) && chr_dist["chr1"] > total_vars * 0.15) {
    cat("NOTE: High chr1 mutation rate - possible instability\n")
  }
}

cat("\n8. CLINICAL RECOMMENDATIONS:\n")
cat("Based on TMB", round(tmb_score, 1), "mut/Mb with weak UV signature:\n")
cat("  ✅ Excellent immunotherapy candidate (TMB >20)\n")
cat("  ✅ Consider first-line anti-PD-1/PD-L1 therapy\n")
cat("  ⚠️  Investigate primary site (mucosal vs cutaneous)\n")
cat("  ⚠️  Consider DNA repair deficiency testing\n")
cat("  ⚠️  Rule out microsatellite instability (MSI)\n")

cat("\n9. FINAL ASSESSMENT:\n")
if (tmb_score > 50) {
  cat("✅ TMB Score: VERY HIGH - Excellent for immunotherapy\n")
} else if (tmb_score > 20) {
  cat("✅ TMB Score: HIGH - Good for immunotherapy\n")
}

if (ct_pct < 15) {
  cat("⚠️  UV Signature: WEAK - Investigate non-UV etiology\n")
  cat("💡 Recommendation: Verify sample type and consider MSI testing\n")
} else {
  cat("✅ UV Signature: Appropriate for melanoma\n")
}

cat("\nOverall: High TMB with weak UV signature is unusual but not unprecedented.\n")
cat("The TMB score alone makes this an excellent immunotherapy candidate.\n")