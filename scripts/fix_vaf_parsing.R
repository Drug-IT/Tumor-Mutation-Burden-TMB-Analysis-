# Fix VAF Parsing and Recalculate TMB
library(dplyr)

cat("=== FIXING VAF PARSING ISSUE ===\n\n")

# Function to properly parse VCF and extract real VAF values
parse_vcf_with_correct_vaf <- function(vcf_file) {
  cat("Reading VCF file and parsing VAF correctly...\n")
  
  # Read VCF file
  vcf_lines <- readLines(vcf_file)
  data_start <- which(!grepl("^#", vcf_lines))[1]
  vcf_data_lines <- vcf_lines[data_start:length(vcf_lines)]
  
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
    Real_VAF = numeric(),  # This will be the correctly calculated VAF
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(vcf_data_lines)) {
    if (i %% 10000 == 0) cat("  Processed", i, "variants...\n")
    
    line <- vcf_data_lines[i]
    fields <- strsplit(line, "\t")[[1]]
    
    if (length(fields) >= 10) {  # Make sure we have FORMAT and sample columns
      chr <- fields[1]
      pos <- as.numeric(fields[2])
      ref <- fields[4]
      alt <- fields[5]
      qual <- as.numeric(fields[6])
      info <- fields[8]
      format_field <- fields[9]
      sample_data <- fields[10]
      
      # Extract from INFO field
      dp <- NA
      af <- NA
      ao <- NA
      ro <- NA
      ab <- NA
      
      info_parts <- strsplit(info, ";")[[1]]
      for (part in info_parts) {
        if (grepl("^DP=", part)) dp <- as.numeric(sub("DP=", "", part))
        if (grepl("^AF=", part)) af <- as.numeric(sub("AF=", "", part))
        if (grepl("^AO=", part)) ao <- as.numeric(sub("AO=", "", part))
        if (grepl("^RO=", part)) ro <- as.numeric(sub("RO=", "", part))
        if (grepl("^AB=", part)) ab <- as.numeric(sub("AB=", "", part))
      }
      
      # Calculate REAL VAF from sample data (FORMAT field)
      real_vaf <- NA
      
      # Parse FORMAT field to understand structure
      format_fields <- strsplit(format_field, ":")[[1]]
      sample_values <- strsplit(sample_data, ":")[[1]]
      
      # Look for AD (Allelic Depth) field
      if ("AD" %in% format_fields && length(sample_values) == length(format_fields)) {
        ad_index <- which(format_fields == "AD")
        ad_values <- sample_values[ad_index]
        
        # Parse AD field (format: ref_depth,alt_depth)
        if (!is.na(ad_values) && grepl(",", ad_values)) {
          depths <- as.numeric(strsplit(ad_values, ",")[[1]])
          if (length(depths) >= 2 && !is.na(depths[1]) && !is.na(depths[2])) {
            ref_depth <- depths[1]
            alt_depth <- depths[2]
            total_depth <- ref_depth + alt_depth
            
            if (total_depth > 0) {
              real_vaf <- alt_depth / total_depth
            }
          }
        }
      }
      
      # If AD not available, try to use AO and RO from INFO
      if (is.na(real_vaf) && !is.na(ao) && !is.na(ro)) {
        total_reads <- ao + ro
        if (total_reads > 0) {
          real_vaf <- ao / total_reads
        }
      }
      
      # If still no VAF, use AB (allele balance) from INFO
      if (is.na(real_vaf) && !is.na(ab)) {
        real_vaf <- ab
      }
      
      # If still no VAF, use AF from INFO (but this might be the problem)
      if (is.na(real_vaf) && !is.na(af)) {
        real_vaf <- af
      }
      
      # Add to data frame
      quality_data <- rbind(quality_data, data.frame(
        Chr = chr,
        Pos = pos,
        Ref = ref,
        Alt = alt,
        QUAL = qual,
        DP = dp,
        AF = af,           # Original AF from INFO
        AO = ao,
        RO = ro,
        AB = ab,
        Real_VAF = real_vaf,  # Correctly calculated VAF
        stringsAsFactors = FALSE
      ))
    }
  }
  
  cat("  Extracted quality data for", nrow(quality_data), "variants\n")
  
  # Show VAF distribution
  cat("\nVAF Distribution Check:\n")
  cat("Original AF values (unique):", length(unique(quality_data$AF)), "\n")
  cat("Real VAF values (unique):", length(unique(quality_data$Real_VAF)), "\n")
  
  if (length(unique(quality_data$Real_VAF)) > 1) {
    cat("Real VAF range:", round(min(quality_data$Real_VAF, na.rm=T), 3), "-", 
        round(max(quality_data$Real_VAF, na.rm=T), 3), "\n")
    cat("Real VAF quartiles:", paste(round(quantile(quality_data$Real_VAF, na.rm=T), 3), collapse=" | "), "\n")
  }
  
  return(quality_data)
}

# Re-parse VCF with correct VAF calculation
vcf_file <- "SRR26456208_normalized.vcf"
cat("Step 1: Re-parsing VCF with correct VAF calculation...\n")
correct_quality_data <- parse_vcf_with_correct_vaf(vcf_file)

# Merge with ANNOVAR results using Real_VAF
cat("\nStep 2: Merging with ANNOVAR results...\n")
annovar_data <- read.csv("filtered_somatic_mutations.csv", stringsAsFactors = FALSE)

# Create matching keys
correct_quality_data$key <- paste(correct_quality_data$Chr, correct_quality_data$Pos, 
                                  correct_quality_data$Ref, correct_quality_data$Alt, sep = "_")
annovar_data$key <- paste(annovar_data$Chr, annovar_data$Start, 
                          annovar_data$Ref, annovar_data$Alt, sep = "_")

# Merge datasets
merged_corrected <- merge(annovar_data, 
                          correct_quality_data[, c("key", "QUAL", "DP", "AF", "AO", "RO", "AB", "Real_VAF")],
                          by = "key", all.x = TRUE)

cat("Merged data:", nrow(merged_corrected), "variants\n")
cat("Variants with Real_VAF data:", sum(!is.na(merged_corrected$Real_VAF)), "\n")

# Apply corrected filtering using Real_VAF
cat("\nStep 3: Applying corrected filtering with Real_VAF...\n")

# Filter with proper VAF
filtered_corrected <- merged_corrected %>%
  filter(
    # Functional mutations
    ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain", "stoploss", 
                              "frameshift deletion", "frameshift insertion",
                              "nonframeshift deletion", "nonframeshift insertion", 
                              "nonframeshift substitution"),
    
    # Quality filters using Real_VAF
    !is.na(Real_VAF) & Real_VAF >= 0.10 & Real_VAF <= 0.90,  # Use Real_VAF
    !is.na(DP) & DP >= 20,
    !is.na(AO) & AO >= 4,
    !is.na(QUAL) & QUAL >= 30,
    
    # Population frequency filters
    (is.na(gnomAD_exome_ALL) | gnomAD_exome_ALL == "." | 
       suppressWarnings(as.numeric(gnomAD_exome_ALL)) <= 0.001),
    (is.na(gnomAD_genome_ALL) | gnomAD_genome_ALL == "." | 
       suppressWarnings(as.numeric(gnomAD_genome_ALL)) <= 0.001)
  )

# Calculate corrected TMB
corrected_tmb <- nrow(filtered_corrected) / 30

cat("\n=== CORRECTED TMB RESULTS ===\n")
cat("Previous TMB (with parsing issue):", 74.3, "mutations/Mb\n")
cat("Corrected TMB (with Real_VAF):", round(corrected_tmb, 2), "mutations/Mb\n")
cat("Previous variants:", 2229, "\n")
cat("Corrected variants:", nrow(filtered_corrected), "\n")

if (nrow(filtered_corrected) > 0) {
  cat("\nCorrected quality distribution:\n")
  cat("Real_VAF range:", round(min(filtered_corrected$Real_VAF, na.rm=T), 3), "-", 
      round(max(filtered_corrected$Real_VAF, na.rm=T), 3), "\n")
  cat("Real_VAF median:", round(median(filtered_corrected$Real_VAF, na.rm=T), 3), "\n")
  cat("Coverage median:", median(filtered_corrected$DP, na.rm=T), "x\n")
  
  # Check VAF distribution
  vaf_dist <- table(cut(filtered_corrected$Real_VAF, 
                        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                        labels = c("10-20%", "20-40%", "40-60%", "60-80%", "80-90%")))
  cat("\nVAF distribution:\n")
  print(vaf_dist)
  
  # UV signature check
  snvs <- filtered_corrected[nchar(filtered_corrected$Ref) == 1 & 
                               nchar(filtered_corrected$Alt) == 1, ]
  if (nrow(snvs) > 0) {
    ct_transitions <- sum(snvs$Ref == "C" & snvs$Alt == "T") + 
      sum(snvs$Ref == "G" & snvs$Alt == "A")
    ct_pct <- ct_transitions / nrow(snvs) * 100
    cat("C>T transitions:", ct_transitions, "(", round(ct_pct, 1), "% of SNVs)\n")
    
    if (ct_pct > 30) {
      cat("UV signature: STRONG\n")
    } else if (ct_pct > 15) {
      cat("UV signature: MODERATE\n")
    } else {
      cat("UV signature: WEAK\n")
    }
  }
  
} else {
  cat("WARNING: No variants passed corrected filtering!\n")
  cat("This suggests the VAF parsing is still having issues.\n")
  
  # Diagnostic: show some example VAF values
  cat("\nDiagnostic - showing first 10 variants with Real_VAF:\n")
  sample_data <- merged_corrected[!is.na(merged_corrected$Real_VAF), 
                                  c("Chr", "Start", "Ref", "Alt", "AF", "Real_VAF", "AO", "RO", "DP")][1:10,]
  print(sample_data)
}

# Save corrected results
if (nrow(filtered_corrected) > 0) {
  write.csv(filtered_corrected, "final_corrected_tmb_results.csv", row.names = FALSE)
  cat("\nCorrected results saved to: final_corrected_tmb_results.csv\n")
}

# TMB classification
if (corrected_tmb < 6) {
  classification <- "Low TMB"
} else if (corrected_tmb < 20) {
  classification <- "Intermediate TMB"
} else {
  classification <- "High TMB"
}

cat("\nFinal Assessment:\n")
cat("TMB Classification:", classification, "\n")

if (corrected_tmb >= 20) {
  cat("Clinical interpretation: Excellent immunotherapy candidate\n")
} else if (corrected_tmb >= 6) {
  cat("Clinical interpretation: Good immunotherapy candidate\n")
} else {
  cat("Clinical interpretation: Limited immunotherapy benefit expected\n")
}

cat("\n=== SUMMARY ===\n")
cat("Issue identified: VAF parsing was using AF=0.5 for all variants\n")
cat("Solution: Used Real_VAF calculated from actual read depths\n")
cat("Result: More accurate TMB calculation\n")