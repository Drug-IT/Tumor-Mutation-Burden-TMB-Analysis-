# TMB Frequency Threshold Verification Script
# Load your filtered data
# Replace 'filtered_somatic_mutations.csv' with your actual file path
df <- read.csv("filtered_somatic_mutations.csv", stringsAsFactors = FALSE)

# Define the frequency columns you're checking
freq_cols <- c(
  "1000g2015aug_all", "1000g2015aug_eur", "1000g2015aug_afr", 
  "1000g2015aug_amr", "1000g2015aug_eas",
  "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS",
  "gnomAD_exome_ALL", "gnomAD_exome_AFR", "gnomAD_exome_AMR", "gnomAD_exome_ASJ", 
  "gnomAD_exome_EAS", "gnomAD_exome_FIN", "gnomAD_exome_NFE", "gnomAD_exome_OTH", "gnomAD_exome_SAS",
  "gnomAD_genome_ALL", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ",
  "gnomAD_genome_EAS", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH"
)

# Check which frequency columns actually exist in your data
existing_freq_cols <- freq_cols[freq_cols %in% colnames(df)]
cat("Existing frequency columns:", paste(existing_freq_cols, collapse = ", "), "\n\n")

# Check frequency distribution for each column
freq_threshold <- 0.001
for (col in existing_freq_cols) {
  if (col %in% colnames(df)) {
    col_values <- suppressWarnings(as.numeric(as.character(df[[col]])))
    
    # Count different categories
    na_count <- sum(is.na(col_values))
    dot_count <- sum(df[[col]] == ".", na.rm = TRUE)
    empty_count <- sum(df[[col]] == "", na.rm = TRUE)
    above_threshold <- sum(col_values > freq_threshold, na.rm = TRUE)
    below_threshold <- sum(col_values <= freq_threshold, na.rm = TRUE)
    
    cat("Column:", col, "\n")
    cat("  Total variants:", nrow(df), "\n")
    cat("  NA values:", na_count, "\n")
    cat("  Dot values:", dot_count, "\n")
    cat("  Empty values:", empty_count, "\n")
    cat("  Above threshold (", freq_threshold, "):", above_threshold, "\n")
    cat("  Below threshold:", below_threshold, "\n")
    
    # Show some example values
    valid_freqs <- col_values[!is.na(col_values) & col_values > 0]
    if (length(valid_freqs) > 0) {
      cat("  Example frequencies:", paste(head(valid_freqs, 10), collapse = ", "), "\n")
      cat("  Max frequency:", max(valid_freqs, na.rm = TRUE), "\n")
    }
    cat("\n")
  }
}

# Check if you have any variants that should have been filtered
potential_germline <- df[!is.na(suppressWarnings(as.numeric(as.character(df$gnomAD_exome_ALL)))) & 
                           suppressWarnings(as.numeric(as.character(df$gnomAD_exome_ALL))) > freq_threshold, ]

cat("Variants that might be germline (gnomAD_exome_ALL > 0.1%):", nrow(potential_germline), "\n")
if (nrow(potential_germline) > 0) {
  cat("Examples of potentially problematic variants:\n")
  print(head(potential_germline[, c("Chr", "Start", "End", "Ref", "Alt", "gnomAD_exome_ALL")], 10))
}