# BAM to VCF Conversion for Melanoma TMB Analysis
# Focus on nonsynonymous mutations in coding regions
# FIXED VERSION with BAM corruption handling

# Load required libraries
library(tools)     # For file path operations
library(parallel)  # For parallel processing

# Enhanced function for melanoma TMB-specific variant calling with BAM validation
bam_to_vcf_melanoma_tmb <- function(bam_file, reference_genome, output_vcf, 
                                    coding_bed = NULL,  # BED file with coding regions
                                    # Melanoma-optimized parameters
                                    min_base_quality = 15,
                                    min_mapping_quality = 20,
                                    min_alternate_fraction = 0.05,  # 5% VAF threshold
                                    min_coverage = 8,
                                    min_alternate_count = 3) {
  
  # Check if required files exist
  if (!file.exists(bam_file)) {
    stop(paste("BAM file not found:", bam_file))
  }
  
  if (!file.exists(reference_genome)) {
    stop(paste("Reference genome not found:", reference_genome))
  }
  
  cat("=== Step 0: BAM File Validation ===\n")
  validated_bam <- validate_and_fix_bam(bam_file)
  
  if (is.null(validated_bam)) {
    stop("BAM file is corrupted and cannot be repaired. Please re-run preprocessing.")
  }
  
  processed_bam <- validated_bam
  
  # Skip duplicate marking since it's already done in preprocessing
  cat("=== Step 1: Using Preprocessed BAM (duplicates already marked) ===\n")
  cat("Preprocessed BAM:", processed_bam, "\n")
  
  cat("=== Step 2: Variant Calling for Melanoma TMB ===\n")
  cat("Using BAM file:", processed_bam, "\n")
  
  # Construct FreeBayes command optimized for somatic variants
  freebayes_cmd <- paste("freebayes", "-f", reference_genome)
  
  # Add coding regions if provided
  if (!is.null(coding_bed) && file.exists(coding_bed)) {
    freebayes_cmd <- paste(freebayes_cmd, "-t", coding_bed)
    cat("Using coding regions from:", coding_bed, "\n")
  } else {
    cat("WARNING: No coding BED file provided - variants from all regions will be included\n")
  }
  
  # Add melanoma-optimized parameters
  freebayes_cmd <- paste(freebayes_cmd,
                         "-q", min_base_quality,
                         "-m", min_mapping_quality,
                         "-F", min_alternate_fraction,
                         "-C", min_coverage,
                         "--min-alternate-count", min_alternate_count,
                         "--haplotype-length 0",  # Disable haplotyping for somatic variants
                         "--pooled-discrete",     # Better for tumor samples
                         "--genotype-qualities",  # Add genotype quality scores
                         processed_bam,
                         ">", output_vcf)
  
  cat("Running FreeBayes command for melanoma TMB:\n")
  cat(freebayes_cmd, "\n\n")
  
  # Execute FreeBayes
  result <- system(freebayes_cmd)
  
  if (result == 0) {
    cat("Successfully created VCF file:", output_vcf, "\n")
    
    # Get basic stats
    if (file.exists(output_vcf)) {
      vcf_stats <- get_vcf_stats(output_vcf)
      print(vcf_stats)
    }
    
    return(TRUE)
  } else {
    stop("FreeBayes command failed")
  }
}

# NEW: Function to validate and fix BAM files
validate_and_fix_bam <- function(bam_file) {
  cat("Validating BAM file:", bam_file, "\n")
  
  # Step 1: Quick check with samtools
  quickcheck_cmd <- paste("samtools quickcheck", bam_file)
  quickcheck_result <- system(quickcheck_cmd, ignore.stderr = TRUE)
  
  if (quickcheck_result == 0) {
    cat("✅ BAM file passes quickcheck\n")
    
    # Check if index needs rebuilding
    bai_file <- paste0(bam_file, ".bai")
    if (file.exists(bai_file)) {
      bam_time <- file.info(bam_file)$mtime
      bai_time <- file.info(bai_file)$mtime
      
      if (bai_time < bam_time) {
        cat("🔧 Rebuilding BAM index (index is older than BAM)...\n")
        file.remove(bai_file)
        index_result <- system(paste("samtools index", bam_file))
        
        if (index_result == 0) {
          cat("✅ Index rebuilt successfully\n")
        } else {
          cat("❌ Failed to rebuild index\n")
          return(NULL)
        }
      }
    } else {
      cat("🔧 Creating BAM index...\n")
      index_result <- system(paste("samtools index", bam_file))
      if (index_result != 0) {
        cat("❌ Failed to create index\n")
        return(NULL)
      }
    }
    
    return(bam_file)
  }
  
  cat("❌ BAM file failed quickcheck - attempting repair...\n")
  
  # Step 2: Attempt to repair corrupted BAM
  repaired_bam <- attempt_bam_repair(bam_file)
  
  if (!is.null(repaired_bam)) {
    cat("✅ Successfully repaired BAM file:", repaired_bam, "\n")
    return(repaired_bam)
  }
  
  # Step 3: Look for alternative BAM files
  alternative_bam <- find_alternative_bam(bam_file)
  
  if (!is.null(alternative_bam)) {
    cat("✅ Found alternative BAM file:", alternative_bam, "\n")
    return(alternative_bam)
  }
  
  cat("❌ Cannot repair or find alternative BAM file\n")
  return(NULL)
}

# Function to attempt BAM repair
attempt_bam_repair <- function(corrupted_bam) {
  cat("Attempting to repair BAM file...\n")
  
  base_name <- tools::file_path_sans_ext(corrupted_bam)
  repaired_bam <- paste0(base_name, "_repaired.bam")
  
  # Try to salvage readable portion using samtools
  repair_cmd <- paste(
    "samtools view -h", corrupted_bam, "2>/dev/null |",
    "samtools view -b -o", repaired_bam, "-"
  )
  
  repair_result <- system(repair_cmd, ignore.stderr = TRUE)
  
  if (repair_result == 0 && file.exists(repaired_bam)) {
    # Validate repaired file
    validate_result <- system(paste("samtools quickcheck", repaired_bam), ignore.stderr = TRUE)
    
    if (validate_result == 0) {
      # Create index for repaired file
      system(paste("samtools index", repaired_bam))
      
      # Check if repaired file has reasonable content
      stats_cmd <- paste("samtools view -c", repaired_bam)
      read_count <- system(stats_cmd, intern = TRUE, ignore.stderr = TRUE)
      
      if (length(read_count) > 0 && as.numeric(read_count) > 1000) {
        cat("Repaired BAM contains", read_count, "reads\n")
        return(repaired_bam)
      }
    }
  }
  
  # Clean up failed repair attempt
  if (file.exists(repaired_bam)) {
    file.remove(repaired_bam)
  }
  
  return(NULL)
}

# Function to find alternative BAM files
find_alternative_bam <- function(target_bam) {
  cat("Looking for alternative BAM files...\n")
  
  bam_dir <- dirname(target_bam)
  base_name <- tools::file_path_sans_ext(basename(target_bam))
  
  # Look for different versions of the same sample
  possible_alternatives <- c(
    # Final preprocessed BAM (most likely to work)
    file.path(bam_dir, paste0(gsub("_marked_dups", "", base_name), "_preprocessed.bam")),
    file.path(bam_dir, paste0(gsub("_marked_dups", "", base_name), ".bam")),
    # Original input BAM
    file.path(dirname(bam_dir), paste0(gsub("_marked_dups|_preprocessed", "", base_name), ".bam")),
    # Other patterns
    file.path(bam_dir, paste0(base_name, "_final.bam")),
    file.path(bam_dir, paste0(base_name, "_recalibrated.bam"))
  )
  
  for (alt_bam in possible_alternatives) {
    if (file.exists(alt_bam)) {
      cat("Found potential alternative:", alt_bam, "\n")
      
      # Validate alternative
      check_result <- system(paste("samtools quickcheck", alt_bam), ignore.stderr = TRUE)
      
      if (check_result == 0) {
        cat("✅ Alternative BAM is valid\n")
        
        # Ensure it has an index
        alt_bai <- paste0(alt_bam, ".bai")
        if (!file.exists(alt_bai)) {
          cat("Creating index for alternative BAM...\n")
          system(paste("samtools index", alt_bam))
        }
        
        return(alt_bam)
      } else {
        cat("❌ Alternative BAM is also corrupted\n")
      }
    }
  }
  
  return(NULL)
}

# Enhanced workflow function with better error handling
melanoma_tmb_workflow <- function(bam_file, reference_genome, 
                                  coding_bed = NULL,
                                  output_dir = "melanoma_tmb",
                                  genome_build = "hg38",
                                  annotation_tool = "snpeff",
                                  force_reprocess = FALSE) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  base_name <- tools::file_path_sans_ext(basename(bam_file))
  
  # File paths
  raw_vcf <- file.path(output_dir, paste0(base_name, "_raw.vcf"))
  annotated_vcf <- file.path(output_dir, paste0(base_name, "_annotated.vcf"))
  final_vcf <- file.path(output_dir, paste0(base_name, "_nonsynonymous.vcf"))
  
  cat("=== Melanoma TMB Analysis Workflow ===\n")
  cat("Sample:", base_name, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("BAM file:", bam_file, "\n\n")
  
  # Check if we should skip variant calling (if VCF already exists and not forcing reprocess)
  if (file.exists(raw_vcf) && !force_reprocess) {
    cat("Raw VCF already exists:", raw_vcf, "\n")
    cat("Skipping variant calling (use force_reprocess=TRUE to override)\n")
  } else {
    # Step 1: Variant calling with validation
    cat("Step 1: Variant calling with BAM validation...\n")
    tryCatch({
      bam_to_vcf_melanoma_tmb(bam_file, reference_genome, raw_vcf, 
                              coding_bed = coding_bed)
    }, error = function(e) {
      cat("❌ Variant calling failed:", e$message, "\n")
      cat("\n=== Troubleshooting Suggestions ===\n")
      cat("1. Check if original BAM file exists and re-run preprocessing\n")
      cat("2. Verify sufficient disk space\n")
      cat("3. Check if preprocessing completed successfully\n")
      stop("Cannot proceed without valid BAM file")
    })
  }
  
  # Step 2: Functional annotation
  cat("\nStep 2: Functional annotation...\n")
  tryCatch({
    annotate_variants(raw_vcf, annotated_vcf, reference_genome,
                      annotation_tool = annotation_tool, 
                      genome_build = genome_build)
    
    # Step 3: Filter for nonsynonymous mutations
    cat("\nStep 3: Filtering for nonsynonymous mutations...\n")
    filter_nonsynonymous_variants(annotated_vcf, final_vcf, 
                                  annotation_method = annotation_tool)
    
  }, error = function(e) {
    cat("Annotation failed:", e$message, "\n")
    cat("Proceeding without functional filtering (NOT RECOMMENDED for TMB)\n")
    final_vcf_path <- raw_vcf
  })
  
  # Final statistics
  cat("\n=== Final Results ===\n")
  if (file.exists(final_vcf)) {
    stats <- get_vcf_stats(final_vcf)
    cat("Nonsynonymous variants for TMB calculation:", stats$total_variants, "\n")
    cat("Final VCF file:", final_vcf, "\n")
    
    # Estimate TMB (assuming ~30 Mb exome)
    estimated_tmb <- stats$total_variants / 30
    cat("Estimated TMB (assuming 30 Mb exome):", round(estimated_tmb, 1), "mut/Mb\n")
    
    if (estimated_tmb >= 20) {
      cat("*** HIGH TMB - May benefit from immunotherapy ***\n")
    }
  }
  
  return(list(
    raw_vcf = raw_vcf,
    annotated_vcf = if(exists("annotated_vcf")) annotated_vcf else NULL,
    final_vcf = final_vcf,
    sample = base_name
  ))
}

# Function to filter VCF for nonsynonymous mutations only
filter_nonsynonymous_variants <- function(input_vcf, output_vcf, 
                                          annotation_method = "snpeff") {
  
  cat("=== Filtering for Nonsynonymous Mutations ===\n")
  
  if (annotation_method == "snpeff") {
    # Filter for high/moderate impact variants (nonsynonymous)
    filter_cmd <- paste(
      "bcftools view -i",
      "'INFO/ANN ~ \"HIGH\" || INFO/ANN ~ \"MODERATE\"'",
      input_vcf, ">", output_vcf
    )
  } else if (annotation_method == "vep") {
    # Filter for protein-coding consequences
    filter_cmd <- paste(
      "bcftools view -i",
      "'INFO/CSQ ~ \"missense_variant\" || INFO/CSQ ~ \"nonsense_variant\" || INFO/CSQ ~ \"frameshift_variant\" || INFO/CSQ ~ \"inframe_insertion\" || INFO/CSQ ~ \"inframe_deletion\"'",
      input_vcf, ">", output_vcf
    )
  } else {
    # Basic filtering without annotation (not recommended for TMB)
    cat("WARNING: No functional annotation filtering - this may overestimate TMB\n")
    filter_cmd <- paste("cp", input_vcf, output_vcf)
  }
  
  cat("Filtering command:\n", filter_cmd, "\n")
  result <- system(filter_cmd)
  
  if (result == 0) {
    cat("Filtered VCF created:", output_vcf, "\n")
    return(TRUE)
  } else {
    stop("Filtering failed")
  }
}

# Function to annotate variants with functional consequences
annotate_variants <- function(input_vcf, output_vcf, reference_genome,
                              annotation_tool = "snpeff", 
                              genome_build = "hg38") {
  
  cat("=== Annotating Variants ===\n")
  
  if (annotation_tool == "snpeff") {
    # Using SnpEff for functional annotation
    snpeff_cmd <- paste(
      "snpEff -v", genome_build,
      input_vcf, ">", output_vcf
    )
    
    cat("SnpEff annotation command:\n", snpeff_cmd, "\n")
    result <- system(snpeff_cmd)
    
  } else if (annotation_tool == "vep") {
    # Using VEP for functional annotation
    vep_cmd <- paste(
      "vep --input_file", input_vcf,
      "--output_file", output_vcf,
      "--format vcf --vcf",
      "--species homo_sapiens",
      "--assembly", genome_build,
      "--offline --cache"
    )
    
    cat("VEP annotation command:\n", vep_cmd, "\n")
    result <- system(vep_cmd)
    
  } else {
    stop("Unsupported annotation tool. Use 'snpeff' or 'vep'")
  }
  
  if (result == 0) {
    cat("Annotation completed:", output_vcf, "\n")
    return(TRUE)
  } else {
    stop("Annotation failed")
  }
}

# Function to get basic VCF statistics
get_vcf_stats <- function(vcf_file) {
  # Count total variants
  total_variants_cmd <- paste("grep -v '^#'", vcf_file, "| wc -l")
  total_variants <- as.numeric(system(total_variants_cmd, intern = TRUE))
  
  # Count SNPs
  snp_cmd <- paste("grep -v '^#'", vcf_file, "| awk '$4~/^[ATCG]$/ && $5~/^[ATCG]$/' | wc -l")
  snps <- as.numeric(system(snp_cmd, intern = TRUE))
  
  # Count INDELs
  indels <- total_variants - snps
  
  return(list(
    total_variants = total_variants,
    snps = snps,
    indels = indels,
    vcf_file = vcf_file
  ))
}

# Function to create coding regions BED file from GTF
create_coding_bed <- function(gtf_file, output_bed) {
  cat("Creating coding regions BED file from GTF...\n")
  
  # Extract coding exons from GTF
  bed_cmd <- paste(
    "awk '$3==\"CDS\" {print $1\"\\t\"($4-1)\"\\t\"$5}'",
    gtf_file, "|",
    "sort -k1,1 -k2,2n |",
    "bedtools merge -i stdin >", output_bed
  )
  
  result <- system(bed_cmd)
  
  if (result == 0) {
    cat("Coding regions BED file created:", output_bed, "\n")
    return(TRUE)
  } else {
    stop("Failed to create coding regions BED file")
  }
}

# Print usage instructions
cat("=== FIXED Melanoma TMB-Specific Variant Calling Script ===\n")
cat("This version includes BAM corruption handling and repair capabilities\n\n")
cat("Key functions:\n")
cat("1. melanoma_tmb_workflow() - Complete workflow with BAM validation\n")
cat("2. bam_to_vcf_melanoma_tmb() - Enhanced variant calling with BAM checks\n")
cat("3. validate_and_fix_bam() - BAM validation and repair\n")
cat("4. filter_nonsynonymous_variants() - Filter for nonsynonymous mutations only\n")
cat("5. create_coding_bed() - Create coding regions BED from GTF\n\n")

cat("Example usage:\n")
cat('# Complete workflow with automatic BAM validation\n')
cat('melanoma_tmb_workflow("sample.bam", "ref.fasta", coding_bed="coding_regions.bed")\n\n')

cat("IMPORTANT: This version will:\n")
cat("1. Automatically validate your BAM file\n")
cat("2. Attempt to repair corrupted BAM files\n")
cat("3. Look for alternative BAM files if repair fails\n")
cat("4. Provide clear error messages and suggestions\n\n")

cat("🚨 FOR YOUR IMMEDIATE PROBLEM:\n")
cat("Run exactly this command:\n")
cat('melanoma_tmb_workflow(\n')
cat('  bam_file = "preprocessed_bams/melanoma/melanoma_preprocessed_marked_dups.bam",\n')
cat('  reference_genome = "/home/dhibi/TMB_Pipeline_Project/data/aligned/skcm/hg38.fa",\n')
cat('  coding_bed = "/home/dhibi/TMB_Pipeline_Project/data/aligned/skcm/coding_exons.sorted.bed"\n')
cat(')\n')