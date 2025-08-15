#!/usr/bin/env Rscript
# ===============================================================================
# COMPLETE ANNOVAR + TMB ANALYSIS PIPELINE FOR CUTANEOUS MELANOMA
# ===============================================================================
# Runs ANNOVAR annotation first, then performs comprehensive TMB analysis
# Addresses germline contamination with aggressive filtering
# ===============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(maftools)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(VennDiagram)
  library(survivalROC)
  library(maxstat)
  library(GSVA)
  library(limma)
})

# ===============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS
# ===============================================================================

# File paths - MODIFY THESE ACCORDING TO YOUR SETUP
SAMPLE_ID <- "SRR26456208"  # Change this to your sample name
VCF_FILE <- "/home/dhibi/TMB_Pipeline_Project/script/melanoma_tmb/melanoma_preprocessed_marked_dups_nonsynonymous.vcf"  # Your input VCF file

# ANNOVAR paths (based on your directory structure)
ANNOVAR_DIR <- "/home/dhibi/TMB_Pipeline_Project/data/aligned/skcm/TMB_analysis/annovar"  # Directory containing ANNOVAR scripts
ANNOVAR_DB <- "/home/dhibi/TMB_Pipeline_Project/data/aligned/skcm/TMB_analysis/annovar/humandb"  # Path to humandb directory
BUILD <- "hg38"  # Genome build (hg38 or hg19)

# Working directories
WORKING_DIR <- getwd()
OUTPUT_DIR <- file.path(WORKING_DIR, "TMB_Analysis_Results")
ANNOVAR_OUTPUT_DIR <- file.path(OUTPUT_DIR, "annovar_output")

# Create output directories
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(ANNOVAR_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Logging function
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

log_message("=== COMPLETE ANNOVAR + TMB ANALYSIS PIPELINE ===")
log_message(paste("Sample ID:", SAMPLE_ID))
log_message(paste("VCF File:", VCF_FILE))
log_message(paste("ANNOVAR Directory:", ANNOVAR_DIR))
log_message(paste("Database Directory:", ANNOVAR_DB))
log_message(paste("Output Directory:", OUTPUT_DIR))

# ===============================================================================
# FUNCTION 1: RUN ANNOVAR ANNOTATION
# ===============================================================================

run_annovar_annotation <- function(vcf_file, sample_id, annovar_dir, annovar_db, 
                                   build = "hg38", output_dir) {
  log_message("STEP 1: Running ANNOVAR Annotation")
  
  # Check if VCF file exists
  if (!file.exists(vcf_file)) {
    stop("VCF file not found: ", vcf_file)
  }
  
  # Define file paths
  avinput_file <- file.path(output_dir, paste0(sample_id, ".avinput"))
  output_prefix <- file.path(output_dir, sample_id)
  
  # ANNOVAR script paths
  convert2annovar <- file.path(annovar_dir, "convert2annovar.pl")
  table_annovar <- file.path(annovar_dir, "table_annovar.pl")
  
  # Check if ANNOVAR scripts exist
  if (!file.exists(convert2annovar)) {
    stop("convert2annovar.pl not found at: ", convert2annovar)
  }
  if (!file.exists(table_annovar)) {
    stop("table_annovar.pl not found at: ", table_annovar)
  }
  
  log_message("Converting VCF to ANNOVAR input format...")
  
  # Step 1: Convert VCF to ANNOVAR input
  convert_cmd <- paste(
    "perl", convert2annovar,
    "-format vcf4",
    "-allsample",
    "-withfreq",
    vcf_file,
    ">", avinput_file
  )
  
  convert_result <- system(convert_cmd, intern = FALSE)
  if (convert_result != 0) {
    stop("Error in convert2annovar step")
  }
  
  log_message("Running comprehensive annotation...")
  
  # Step 2: Run table_annovar with comprehensive annotation
  # This includes all the databases you have in humandb
  protocols <- paste(c(
    "refGene",
    "ensGene", 
    "avsnp150",
    "1000g2015aug_all",
    "1000g2015aug_eur",
    "1000g2015aug_afr",
    "1000g2015aug_amr",
    "1000g2015aug_eas",
    "exac03",
    "gnomad_exome",
    "gnomad_genome"
  ), collapse = ",")
  
  operations <- paste(c(
    "g",  # refGene
    "g",  # ensGene
    "f",  # avsnp150
    "f",  # 1000g2015aug_all
    "f",  # 1000g2015aug_eur
    "f",  # 1000g2015aug_afr
    "f",  # 1000g2015aug_amr
    "f",  # 1000g2015aug_eas
    "f",  # exac03
    "f",  # gnomad_exome
    "f"   # gnomad_genome
  ), collapse = ",")
  
  annotate_cmd <- paste(
    "perl", table_annovar,
    avinput_file,
    annovar_db,
    "-buildver", build,
    "-out", output_prefix,
    "-remove",
    "-protocol", protocols,
    "-operation", operations,
    "-nastring .",
    "-csvout",
    "-polish"
  )
  
  annotate_result <- system(annotate_cmd, intern = FALSE)
  if (annotate_result != 0) {
    stop("Error in table_annovar step")
  }
  
  # Define output file
  annotated_file <- paste0(output_prefix, ".", build, "_multianno.csv")
  
  if (!file.exists(annotated_file)) {
    stop("ANNOVAR annotation failed - output file not found: ", annotated_file)
  }
  
  log_message(paste("ANNOVAR annotation completed successfully"))
  log_message(paste("Annotated file:", annotated_file))
  
  return(annotated_file)
}

# ===============================================================================
# FUNCTION 2: ADVANCED TMB CALCULATION WITH GERMLINE FILTERING
# ===============================================================================

calculate_tmb_advanced_filtering <- function(annovar_file, sample_id) {
  log_message("STEP 2: Advanced TMB Calculation with Germline Filtering")
  
  if (!file.exists(annovar_file)) {
    stop("ANNOVAR file not found: ", annovar_file)
  }
  
  # Read ANNOVAR output
  df <- read_csv(annovar_file, show_col_types = FALSE)
  initial_count <- nrow(df)
  log_message(paste("Initial variants loaded:", formatC(initial_count, format="d", big.mark=",")))
  
  current_df <- df
  filtering_stats <- data.frame(
    Step = character(),
    Count = integer(),
    Percentage = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Add initial count
  filtering_stats <- rbind(filtering_stats, 
                           data.frame(Step = "Initial Variants", 
                                      Count = initial_count, 
                                      Percentage = 100.0))
  
  # Step 1: Keep only exonic and splicing variants
  if ("Func.refGene" %in% colnames(current_df)) {
    current_df <- current_df %>% 
      filter(Func.refGene %in% c("exonic", "splicing"))
    count_after <- nrow(current_df)
    filtering_stats <- rbind(filtering_stats,
                             data.frame(Step = "Exonic/Splicing Only",
                                        Count = count_after,
                                        Percentage = (count_after/initial_count)*100))
    log_message(paste("After exonic/splicing filter:", formatC(count_after, format="d", big.mark=",")))
  }
  
  # Step 2: Remove synonymous mutations
  if ("ExonicFunc.refGene" %in% colnames(current_df)) {
    functional_mutations <- c(
      "nonsynonymous SNV", "stopgain", "stoploss", 
      "frameshift deletion", "frameshift insertion", 
      "nonframeshift deletion", "nonframeshift insertion", 
      "nonframeshift substitution"
    )
    current_df <- current_df %>% 
      filter(ExonicFunc.refGene %in% functional_mutations)
    count_after <- nrow(current_df)
    filtering_stats <- rbind(filtering_stats,
                             data.frame(Step = "Functional Mutations Only",
                                        Count = count_after,
                                        Percentage = (count_after/initial_count)*100))
    log_message(paste("After functional mutation filter:", formatC(count_after, format="d", big.mark=",")))
  }
  
  # Step 3: Remove common variants from population databases
  population_freq_cols <- c(
    "1000g2015aug_all", "1000g2015aug_eur", "1000g2015aug_afr", 
    "1000g2015aug_amr", "1000g2015aug_eas",
    "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH", "ExAC_SAS",
    "gnomAD_exome_ALL", "gnomAD_exome_AFR", "gnomAD_exome_AMR", "gnomAD_exome_ASJ", 
    "gnomAD_exome_EAS", "gnomAD_exome_FIN", "gnomAD_exome_NFE", "gnomAD_exome_OTH", "gnomAD_exome_SAS",
    "gnomAD_genome_ALL", "gnomAD_genome_AFR", "gnomAD_genome_AMR", "gnomAD_genome_ASJ",
    "gnomAD_genome_EAS", "gnomAD_genome_FIN", "gnomAD_genome_NFE", "gnomAD_genome_OTH"
  )
  
  # Find which frequency columns exist in the data
  existing_freq_cols <- intersect(population_freq_cols, colnames(current_df))
  log_message(paste("Found population frequency columns:", length(existing_freq_cols)))
  
  if (length(existing_freq_cols) > 0) {
    # Very stringent frequency threshold for somatic variants
    freq_threshold <- 0.001  # 0.1%
    
    for (col in existing_freq_cols) {
      col_values <- suppressWarnings(as.numeric(as.character(current_df[[col]])))
      # Keep variants with no frequency data OR frequency <= threshold
      keep_variants <- is.na(col_values) | 
        current_df[[col]] == "." | 
        current_df[[col]] == "" | 
        col_values <= freq_threshold
      current_df <- current_df[keep_variants, ]
      if (nrow(current_df) == 0) break
    }
    
    count_after <- nrow(current_df)
    filtering_stats <- rbind(filtering_stats,
                             data.frame(Step = "Population Frequency Filter",
                                        Count = count_after,
                                        Percentage = (count_after/initial_count)*100))
    log_message(paste("After population frequency filter:", formatC(count_after, format="d", big.mark=",")))
  }
  
  # Step 4: Remove known dbSNP variants
  if ("avsnp150" %in% colnames(current_df)) {
    # Keep variants NOT in dbSNP (likely somatic)
    keep_variants <- is.na(current_df$avsnp150) | 
      current_df$avsnp150 == "." | 
      current_df$avsnp150 == ""
    current_df <- current_df[keep_variants, ]
    count_after <- nrow(current_df)
    filtering_stats <- rbind(filtering_stats,
                             data.frame(Step = "dbSNP Filter",
                                        Count = count_after,
                                        Percentage = (count_after/initial_count)*100))
    log_message(paste("After dbSNP filter:", formatC(count_after, format="d", big.mark=",")))
  }
  
  # Calculate final TMB
  final_mutation_count <- nrow(current_df)
  
  # Exome size - adjust based on your sequencing approach
  exome_length_mb <- 30.0  # Standard whole exome
  # For targeted panels, you would adjust this value
  
  tmb_score <- final_mutation_count / exome_length_mb
  
  
  log_message("FILTERING SUMMARY:")
  log_message(paste("Initial variants:", formatC(initial_count, format="d", big.mark=",")))
  log_message(paste("Final mutations:", formatC(final_mutation_count, format="d", big.mark=",")))
  log_message(paste("Variants removed:", formatC(initial_count - final_mutation_count, format="d", big.mark=",")))
  log_message(paste("Retention rate:", sprintf("%.2f%%", (final_mutation_count/initial_count)*100)))
  log_message(paste("TMB Score:", sprintf("%.2f", tmb_score), "mutations/Mb"))
 
  
  return(list(
    initial_count = initial_count,
    final_mutations = current_df,
    mutation_count = final_mutation_count,
    tmb_score = tmb_score,
    exome_length_mb = exome_length_mb,
    retention_rate = (final_mutation_count/initial_count)*100,
    filtering_stats = filtering_stats
  ))
}

# ===============================================================================
# FUNCTION 3: TMB CLASSIFICATION FOR MELANOMA
# ===============================================================================

classify_tmb_melanoma <- function(tmb_score) {
  # Melanoma-specific TMB thresholds based on clinical literature
  if (tmb_score < 10) {
    classification <- "Low TMB"
    color <- "#3498db"  # Blue
    interpretation <- "Limited neoantigen burden; may have reduced immunotherapy response"
  } else if (tmb_score < 20) {
    classification <- "Intermediate TMB" 
    color <- "#f39c12"  # Orange
    interpretation <- "Moderate neoantigen burden; may benefit from immunotherapy"
  } else if (tmb_score < 50) {
    classification <- "High TMB"
    color <- "#e74c3c"  # Red
    interpretation <- "High neoantigen burden; likely to benefit from immunotherapy"
  } else {
    classification <- "Very High TMB"
    color <- "#8e44ad"  # Purple
    interpretation <- "Very high neoantigen burden; excellent immunotherapy candidate"
  }
  
  return(list(
    classification = classification, 
    color = color,
    interpretation = interpretation
  ))
}

# ===============================================================================
# FUNCTION 4: MUTATION SIGNATURE ANALYSIS
# ===============================================================================

analyze_mutation_signatures <- function(mutation_df) {
  log_message("STEP 3: Mutation Signature Analysis")
  
  if (nrow(mutation_df) == 0) {
    log_message("No mutations available for signature analysis")
    return(NULL)
  }
  
  signature_results <- list()
  
  # Analyze mutation types
  if ("ExonicFunc.refGene" %in% colnames(mutation_df)) {
    mutation_types <- table(mutation_df$ExonicFunc.refGene)
    signature_results$mutation_types <- mutation_types
    
    log_message("Mutation type distribution:")
    for (i in 1:length(mutation_types)) {
      log_message(paste("  ", names(mutation_types)[i], ":", mutation_types[i]))
    }
  }
  
  # Analyze substitution patterns if available
  if (all(c("Ref", "Alt") %in% colnames(mutation_df))) {
    snvs <- mutation_df %>% 
      filter(nchar(Ref) == 1 & nchar(Alt) == 1) %>%
      mutate(
        mutation_type = paste0(Ref, ">", Alt),
        transition = case_when(
          (Ref == "A" & Alt == "G") | (Ref == "G" & Alt == "A") |
            (Ref == "C" & Alt == "T") | (Ref == "T" & Alt == "C") ~ "Transition",
          TRUE ~ "Transversion"
        )
      )
    
    if (nrow(snvs) > 0) {
      substitution_counts <- table(snvs$mutation_type)
      transition_counts <- table(snvs$transition)
      
      signature_results$substitutions <- substitution_counts
      signature_results$transitions <- transition_counts
      
      transition_ratio <- sum(snvs$transition == "Transition") / nrow(snvs)
      log_message(paste("Transition/Transversion ratio:", sprintf("%.2f", transition_ratio)))
      
      # C>T transitions are hallmark of UV damage in melanoma
      ct_mutations <- sum(snvs$mutation_type == "C>T")
      ct_percentage <- (ct_mutations / nrow(snvs)) * 100
      log_message(paste("C>T mutations (UV signature):", ct_mutations, 
                        sprintf("(%.1f%%)", ct_percentage)))
      
      signature_results$ct_signature <- list(count = ct_mutations, percentage = ct_percentage)
      signature_results$ti_tv_ratio <- transition_ratio
    }
  }
  
  signature_results$total_mutations <- nrow(mutation_df)
  return(signature_results)
}

# ===============================================================================
# FUNCTION 5: COMPREHENSIVE PLOTTING
# ===============================================================================

create_comprehensive_plots <- function(tmb_results, mutation_analysis, sample_id, output_dir) {
  log_message("STEP 4: Creating Comprehensive Plots")
  
  tmb_class <- classify_tmb_melanoma(tmb_results$tmb_score)
  
  # Plot 1: TMB Summary with Classification
  p1 <- ggplot(data.frame(Sample = sample_id, TMB = tmb_results$tmb_score), 
               aes(x = Sample, y = TMB)) +
    geom_col(fill = tmb_class$color, alpha = 0.8, width = 0.6) +
    geom_text(aes(label = paste0(sprintf("%.1f", TMB), " mut/Mb\n", tmb_class$classification)),
              vjust = -0.5, size = 4, fontface = "bold") +
    labs(title = "Tumor Mutational Burden Analysis",
         subtitle = paste("Sample:", sample_id, "| Classification:", tmb_class$classification),
         x = "Sample", y = "TMB (mutations/Mb)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank()
    ) +
    ylim(0, max(50, tmb_results$tmb_score * 1.2))
  
  ggsave(file.path(output_dir, "TMB_summary.pdf"), p1, width = 10, height = 6)
  ggsave(file.path(output_dir, "TMB_summary.png"), p1, width = 10, height = 6, dpi = 300)
  
  # Plot 2: Filtering Pipeline Waterfall
  if (!is.null(tmb_results$filtering_stats)) {
    p2 <- ggplot(tmb_results$filtering_stats, 
                 aes(x = factor(Step, levels = rev(Step)), y = Count)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = paste0(scales::comma(Count), "\n(", sprintf("%.1f%%", Percentage), ")")), 
                hjust = -0.1, size = 3.5) +
      coord_flip() +
      labs(title = "Variant Filtering Pipeline",
           subtitle = paste("Final Retention Rate:", sprintf("%.2f%%", tmb_results$retention_rate)),
           x = "Filtering Step", y = "Number of Variants") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      scale_y_continuous(labels = scales::comma)
    
    ggsave(file.path(output_dir, "filtering_pipeline.pdf"), p2, width = 12, height = 8)
    ggsave(file.path(output_dir, "filtering_pipeline.png"), p2, width = 12, height = 8, dpi = 300)
  }
  
  # Plot 3: Mutation Type Distribution
  if (!is.null(mutation_analysis) && !is.null(mutation_analysis$mutation_types)) {
    mut_df <- data.frame(
      Type = names(mutation_analysis$mutation_types),
      Count = as.numeric(mutation_analysis$mutation_types),
      stringsAsFactors = FALSE
    )
    
    p3 <- ggplot(mut_df, aes(x = reorder(Type, Count), y = Count)) +
      geom_col(fill = "coral", alpha = 0.8) +
      geom_text(aes(label = Count), hjust = -0.1, size = 3.5) +
      coord_flip() +
      labs(title = "Mutation Type Distribution",
           subtitle = paste("Total Functional Mutations:", sum(mut_df$Count)),
           x = "Mutation Type", y = "Count") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")
      )
    
    ggsave(file.path(output_dir, "mutation_types.pdf"), p3, width = 10, height = 6)
    ggsave(file.path(output_dir, "mutation_types.png"), p3, width = 10, height = 6, dpi = 300)
  }
  
  # Plot 4: UV Signature Analysis (C>T mutations)
  if (!is.null(mutation_analysis) && !is.null(mutation_analysis$substitutions)) {
    sub_df <- data.frame(
      Substitution = names(mutation_analysis$substitutions),
      Count = as.numeric(mutation_analysis$substitutions),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        UV_signature = ifelse(Substitution == "C>T", "UV-related", "Other"),
        Color = ifelse(Substitution == "C>T", "#e74c3c", "#95a5a6")
      )
    
    p4 <- ggplot(sub_df, aes(x = reorder(Substitution, Count), y = Count, fill = Color)) +
      geom_col(alpha = 0.8) +
      scale_fill_identity() +
      geom_text(aes(label = Count), hjust = -0.1, size = 3.5) +
      coord_flip() +
      labs(title = "Substitution Pattern Analysis",
           subtitle = paste("C>T Mutations (UV signature):", 
                            ifelse(!is.null(mutation_analysis$ct_signature), 
                                   paste0(mutation_analysis$ct_signature$count, " (", 
                                          sprintf("%.1f%%", mutation_analysis$ct_signature$percentage), ")"),
                                   "N/A")),
           x = "Substitution Type", y = "Count") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold")
      )
    
    ggsave(file.path(output_dir, "substitution_patterns.pdf"), p4, width = 10, height = 6)
    ggsave(file.path(output_dir, "substitution_patterns.png"), p4, width = 10, height = 6, dpi = 300)
  }
  
  log_message("All plots saved successfully")
}

# ===============================================================================
# FUNCTION 6: GENERATE COMPREHENSIVE REPORT
# ===============================================================================

generate_comprehensive_report <- function(tmb_results, mutation_analysis, sample_id, output_dir) {
  log_message("STEP 5: Generating Comprehensive Report")
  
  tmb_class <- classify_tmb_melanoma(tmb_results$tmb_score)
  
  report_file <- file.path(output_dir, "TMB_Analysis_Comprehensive_Report.txt")
  
  report_content <- paste0(
    "=" , paste(rep("=", 80), collapse = ""), "\n",
    "TUMOR MUTATIONAL BURDEN ANALYSIS - COMPREHENSIVE REPORT\n",
    "=" , paste(rep("=", 80), collapse = ""), "\n\n",
    "Sample ID: ", sample_id, "\n",
    "Analysis Date: ", Sys.time(), "\n",
    "Pipeline Version: Complete ANNOVAR + TMB Analysis v1.0\n",
    "Cancer Type: Cutaneous Melanoma\n\n",
    
    "SUMMARY RESULTS:\n",
    "-" , paste(rep("-", 50), collapse = ""), "\n",
    "TMB Score: ", sprintf("%.2f", tmb_results$tmb_score), " mutations/Mb\n",
    "TMB Classification: ", tmb_class$classification, "\n",
    "Clinical Interpretation: ", tmb_class$interpretation, "\n",
    "Total Somatic Mutations: ", formatC(tmb_results$mutation_count, format="d", big.mark=","), "\n",
    "Exome Size Used: ", tmb_results$exome_length_mb, " Mb\n\n",
    
    "FILTERING STATISTICS:\n",
    "-" , paste(rep("-", 50), collapse = ""), "\n",
    "Initial Variants: ", formatC(tmb_results$initial_count, format="d", big.mark=","), "\n",
    "Final Somatic Mutations: ", formatC(tmb_results$mutation_count, format="d", big.mark=","), "\n",
    "Variants Filtered Out: ", formatC(tmb_results$initial_count - tmb_results$mutation_count, format="d", big.mark=","), "\n",
    "Retention Rate: ", sprintf("%.2f%%", tmb_results$retention_rate), "\n\n"
  )
  
  # Add filtering step details if available
  if (!is.null(tmb_results$filtering_stats)) {
    report_content <- paste0(report_content,
                             "DETAILED FILTERING STEPS:\n",
                             "-" , paste(rep("-", 30), collapse = ""), "\n")
    
    for (i in 1:nrow(tmb_results$filtering_stats)) {
      step <- tmb_results$filtering_stats[i, ]
      report_content <- paste0(report_content,
                               sprintf("%-25s: %8s variants (%.1f%%)\n", 
                                       step$Step, 
                                       formatC(step$Count, format="d", big.mark=","), 
                                       step$Percentage))
    }
    report_content <- paste0(report_content, "\n")
  }
  
  # Add mutation signature analysis
  if (!is.null(mutation_analysis)) {
    report_content <- paste0(report_content,
                             "MUTATION SIGNATURE ANALYSIS:\n",
                             "-" , paste(rep("-", 50), collapse = ""), "\n")
    
    if (!is.null(mutation_analysis$ct_signature)) {
      report_content <- paste0(report_content,
                               "C>T Mutations (UV signature): ", mutation_analysis$ct_signature$count, 
                               " (", sprintf("%.1f%%", mutation_analysis$ct_signature$percentage), ")\n")
    }
    
    if (!is.null(mutation_analysis$ti_tv_ratio)) {
      report_content <- paste0(report_content,
                               "Transition/Transversion Ratio: ", sprintf("%.2f", mutation_analysis$ti_tv_ratio), "\n")
    }
    
    report_content <- paste0(report_content, "\n")
  }
  
  # Add clinical context
  report_content <- paste0(report_content,
                           "CLINICAL CONTEXT FOR MELANOMA:\n",
                           "-" , paste(rep("-", 50), collapse = ""), "\n",
                           "• Melanomas typically have TMB ranges of 10-50 mut/Mb\n",
                           "• High C>T mutation rates indicate UV damage signature\n",
                           "• TMB is predictive of immunotherapy response\n",
                           "• Current sample TMB: ", tmb_class$classification, "\n\n"
  )
  
  # Add recommendations
  if (tmb_results$tmb_score < 10) {
    report_content <- paste0(report_content,
                             "CLINICAL RECOMMENDATIONS:\n",
                             "-" , paste(rep("-", 30), collapse = ""), "\n",
                             "• Consider targeted therapy options\n",
                             "• Evaluate for BRAF/NRAS mutations\n",
                             "• Immunotherapy may have limited efficacy\n",
                             "• Consider combination approaches\n\n")
  } else if (tmb_results$tmb_score < 20) {
    report_content <- paste0(report_content,
                             "CLINICAL RECOMMENDATIONS:\n",
                             "-" , paste(rep("-", 30), collapse = ""), "\n",
                             "• Good candidate for immunotherapy\n",
                             "• Consider anti-PD-1/PD-L1 therapy\n",
                             "• May benefit from combination treatments\n",
                             "• Monitor response carefully\n\n")
  } else {
    report_content <- paste0(report_content,
                             "CLINICAL RECOMMENDATIONS:\n",
                             "-" , paste(rep("-", 30), collapse = ""), "\n",
                             "• Excellent candidate for immunotherapy\n",
                             "• High likelihood of anti-PD-1/PD-L1 response\n",
                             "• Consider first-line immunotherapy\n",
                             "• Monitor for immune-related adverse events\n\n")
  }
  
  report_content <- paste0(report_content,
                           "QUALITY CONTROL MEASURES:\n",
                           "-" , paste(rep("-", 50), collapse = ""), "\n",
                           "• Comprehensive population database filtering applied\n",
                           "• Stringent frequency thresholds (0.1%) for rare variants\n",
                           "• dbSNP variants removed to eliminate germline contamination\n",
                           "• Only functional mutations counted toward TMB\n",
                           "• Multiple annotation databases cross-referenced\n\n",
                           
                           "FILES GENERATED:\n",
                           "-" , paste(rep("-", 30), collapse = ""), "\n",
                           "• TMB_summary.pdf/png - TMB score visualization\n",
                           "• filtering_pipeline.pdf/png - Detailed filtering statistics\n",
                           "• mutation_types.pdf/png - Mutation type distribution\n",
                           "• substitution_patterns.pdf/png - UV signature analysis\n",
                           "• filtered_somatic_mutations.csv - Final somatic mutation list\n",
                           "• TMB_Analysis_Comprehensive_Report.txt - This report\n",
                           "• annovar_output/ - Complete ANNOVAR annotation files\n\n",
                           
                           "REFERENCES:\n",
                           "-" , paste(rep("-", 20), collapse = ""), "\n",
                           "• Samstein et al. Tumor mutational load predicts survival after immunotherapy across multiple cancer types. Nat Genet. 2019\n",
                           "• Alexandrov et al. Signatures of mutational processes in human cancer. Nature. 2013\n",
                           "• Hugo et al. Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. Cell. 2016\n\n",
                           
                           "=" , paste(rep("=", 80), collapse = ""), "\n",
                           "END OF REPORT\n",
                           "=" , paste(rep("=", 80), collapse = ""), "\n"
  )
  
  writeLines(report_content, report_file)
  log_message(paste("Comprehensive report saved to:", report_file))
  
  # Save filtered somatic mutations
  write_csv(tmb_results$final_mutations, 
            file.path(output_dir, "filtered_somatic_mutations.csv"))
  log_message("Filtered somatic mutations saved to: filtered_somatic_mutations.csv")
  
  # Save summary statistics
  summary_stats <- data.frame(
    Metric = c("Sample_ID", "TMB_Score", "TMB_Classification", "Total_Mutations", 
               "Initial_Variants", "Retention_Rate_Percent", "Exome_Size_Mb"),
    Value = c(sample_id, sprintf("%.2f", tmb_results$tmb_score), tmb_class$classification,
              tmb_results$mutation_count, tmb_results$initial_count, 
              sprintf("%.2f", tmb_results$retention_rate), tmb_results$exome_length_mb),
    stringsAsFactors = FALSE
  )
  
  write_csv(summary_stats, file.path(output_dir, "TMB_summary_statistics.csv"))
  log_message("Summary statistics saved to: TMB_summary_statistics.csv")
}

# ===============================================================================
# MAIN EXECUTION PIPELINE
# ===============================================================================

main_complete_pipeline <- function() {
  log_message("Starting Complete ANNOVAR + TMB Analysis Pipeline")
  
  tryCatch({
    # Step 1: Run ANNOVAR annotation
    annotated_file <- run_annovar_annotation(
      vcf_file = VCF_FILE,
      sample_id = SAMPLE_ID,
      annovar_dir = ANNOVAR_DIR,
      annovar_db = ANNOVAR_DB,
      build = BUILD,
      output_dir = ANNOVAR_OUTPUT_DIR
    )
    
    # Step 2: Calculate TMB with advanced filtering
    tmb_results <- calculate_tmb_advanced_filtering(annotated_file, SAMPLE_ID)
    
    # Step 3: Analyze mutation signatures
    mutation_analysis <- analyze_mutation_signatures(tmb_results$final_mutations)
    
    # Step 4: Create comprehensive plots
    create_comprehensive_plots(tmb_results, mutation_analysis, SAMPLE_ID, OUTPUT_DIR)
    
    # Step 5: Generate comprehensive report
    generate_comprehensive_report(tmb_results, mutation_analysis, SAMPLE_ID, OUTPUT_DIR)
    
    # Final summary
    tmb_class <- classify_tmb_melanoma(tmb_results$tmb_score)
    
    log_message("COMPLETE PIPELINE FINISHED SUCCESSFULLY")
    
    log_message(paste("Sample:", SAMPLE_ID))
    log_message(paste("Initial variants:", formatC(tmb_results$initial_count, format="d", big.mark=",")))
    log_message(paste("Final somatic mutations:", formatC(tmb_results$mutation_count, format="d", big.mark=",")))
    log_message(paste("Retention rate:", sprintf("%.2f%%", tmb_results$retention_rate)))
    log_message(paste("Final TMB Score:", sprintf("%.2f", tmb_results$tmb_score), "mut/Mb"))
    log_message(paste("TMB Classification:", tmb_class$classification))
    log_message(paste("Clinical Interpretation:", tmb_class$interpretation))
   
    log_message(paste("All results saved to:", OUTPUT_DIR))
    
    
    return(list(
      tmb_results = tmb_results,
      mutation_analysis = mutation_analysis,
      output_directory = OUTPUT_DIR
    ))
    
  }, error = function(e) {
    log_message(paste("ERROR in pipeline:", e$message))
    log_message("Check your file paths and ANNOVAR installation")
    stop(e)
  })
}

# ===============================================================================
# UTILITY FUNCTIONS FOR BATCH PROCESSING
# ===============================================================================

# Function for processing multiple VCF files
process_multiple_vcf_files <- function(vcf_directory) {
  log_message("Starting batch processing of multiple VCF files")
  
  # Find all VCF files in directory
  vcf_files <- list.files(vcf_directory, pattern = "\\.vcf$", full.names = TRUE)
  
  if (length(vcf_files) == 0) {
    stop("No VCF files found in directory: ", vcf_directory)
  }
  
  log_message(paste("Found", length(vcf_files), "VCF files"))
  
  all_results <- list()
  batch_summary <- data.frame(
    Sample_ID = character(),
    TMB_Score = numeric(),
    TMB_Classification = character(),
    Total_Mutations = integer(),
    Retention_Rate = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(vcf_files)) {
    vcf_file <- vcf_files[i]
    sample_id <- tools::file_path_sans_ext(basename(vcf_file))
    
    log_message(paste("Processing sample", i, "of", length(vcf_files), ":", sample_id))
    
    # Temporarily update global variables
    original_vcf <- VCF_FILE
    original_sample <- SAMPLE_ID
    
    VCF_FILE <<- vcf_file
    SAMPLE_ID <<- sample_id
    
    tryCatch({
      # Run complete pipeline for this sample
      result <- main_complete_pipeline()
      all_results[[sample_id]] <- result
      
      # Add to batch summary
      tmb_class <- classify_tmb_melanoma(result$tmb_results$tmb_score)
      batch_summary <- rbind(batch_summary, 
                             data.frame(
                               Sample_ID = sample_id,
                               TMB_Score = result$tmb_results$tmb_score,
                               TMB_Classification = tmb_class$classification,
                               Total_Mutations = result$tmb_results$mutation_count,
                               Retention_Rate = result$tmb_results$retention_rate,
                               stringsAsFactors = FALSE
                             ))
      
    }, error = function(e) {
      log_message(paste("Error processing", sample_id, ":", e$message))
    })
    
    # Restore original values
    VCF_FILE <<- original_vcf
    SAMPLE_ID <<- original_sample
  }
  
  # Save batch summary
  batch_output_dir <- file.path(WORKING_DIR, "Batch_TMB_Analysis")
  dir.create(batch_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_csv(batch_summary, file.path(batch_output_dir, "batch_TMB_summary.csv"))
  
  # Create batch comparison plot
  if (nrow(batch_summary) > 1) {
    p_batch <- ggplot(batch_summary, aes(x = reorder(Sample_ID, TMB_Score), y = TMB_Score)) +
      geom_col(aes(fill = TMB_Classification), alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f", TMB_Score)), hjust = -0.1, size = 3) +
      coord_flip() +
      scale_fill_manual(values = c("Low TMB" = "#3498db", "Intermediate TMB" = "#f39c12", 
                                   "High TMB" = "#e74c3c", "Very High TMB" = "#8e44ad")) +
      labs(title = "TMB Comparison Across Samples",
           subtitle = paste("Total Samples:", nrow(batch_summary)),
           x = "Sample", y = "TMB (mutations/Mb)", fill = "TMB Classification") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "bottom"
      )
    
    ggsave(file.path(batch_output_dir, "batch_TMB_comparison.pdf"), p_batch, 
           width = 12, height = max(6, nrow(batch_summary) * 0.4))
    ggsave(file.path(batch_output_dir, "batch_TMB_comparison.png"), p_batch, 
           width = 12, height = max(6, nrow(batch_summary) * 0.4), dpi = 300)
  }
  
  log_message(paste("Batch processing complete. Results saved to:", batch_output_dir))
  return(list(results = all_results, summary = batch_summary))
}

# ===============================================================================
# INSTRUCTIONS AND SETUP VERIFICATION
# ===============================================================================

verify_setup <- function() {
  log_message("Verifying setup and dependencies...")
  
  errors <- c()
  
  # Check if VCF file exists
  if (!file.exists(VCF_FILE)) {
    errors <- c(errors, paste("VCF file not found:", VCF_FILE))
  }
  
  # Check ANNOVAR directory
  if (!dir.exists(ANNOVAR_DIR)) {
    errors <- c(errors, paste("ANNOVAR directory not found:", ANNOVAR_DIR))
  }
  
  # Check ANNOVAR database
  if (!dir.exists(ANNOVAR_DB)) {
    errors <- c(errors, paste("ANNOVAR database directory not found:", ANNOVAR_DB))
  }
  
  # Check ANNOVAR scripts
  convert2annovar <- file.path(ANNOVAR_DIR, "convert2annovar.pl")
  table_annovar <- file.path(ANNOVAR_DIR, "table_annovar.pl")
  
  if (!file.exists(convert2annovar)) {
    errors <- c(errors, paste("convert2annovar.pl not found:", convert2annovar))
  }
  
  if (!file.exists(table_annovar)) {
    errors <- c(errors, paste("table_annovar.pl not found:", table_annovar))
  }
  
  if (length(errors) > 0) {
    log_message("SETUP ERRORS FOUND:")
    for (error in errors) {
      log_message(paste("ERROR:", error))
    }
    log_message("\nPlease fix these issues before running the pipeline.")
    return(FALSE)
  } else {
    log_message("Setup verification passed! Ready to run pipeline.")
    return(TRUE)
  }
}

# ===============================================================================
# MAIN EXECUTION INSTRUCTIONS
# ===============================================================================


log_message("COMPLETE ANNOVAR + TMB ANALYSIS PIPELINE LOADED")

log_message("")
log_message("SETUP INSTRUCTIONS:")
log_message("1. Modify the configuration section at the top of this script:")
log_message("   - Set SAMPLE_ID to your sample name")
log_message("   - Set VCF_FILE to the path of your VCF file")
log_message("   - Set ANNOVAR_DIR to your ANNOVAR installation directory")
log_message("   - Set ANNOVAR_DB to your humandb directory")
log_message("")
log_message("2. Run the pipeline:")
log_message("   - For single sample: main_complete_pipeline()")
log_message("   - For batch processing: process_multiple_vcf_files('/path/to/vcf/directory')")
log_message("")
log_message("3. Check setup first: verify_setup()")
log_message("")
log_message("EXAMPLE CONFIGURATION:")
log_message('SAMPLE_ID <- "SKCM_001"')
log_message('VCF_FILE <- "/home/user/data/sample.vcf"')
log_message('ANNOVAR_DIR <- "/home/user/annovar"')
log_message('ANNOVAR_DB <- "/home/user/annovar/humandb"')
log_message("")
log_message("The pipeline will:")
log_message("- Run comprehensive ANNOVAR annotation")
log_message("- Apply aggressive germline filtering")
log_message("- Calculate realistic TMB for melanoma")
log_message("- Generate publication-quality plots")
log_message("- Create detailed clinical report")


# Uncomment the following line to run setup verification
# verify_setup()

# Uncomment the following line to run the complete pipeline
main_complete_pipeline()