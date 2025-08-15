# TMB Pipeline in RStudio - ANNOVAR-Based Melanoma Sample SRR26456208
# ------------------------------------------------------------------
# This script:
# 1) Loads and filters VCF
# 2) Fixes VCF header & normalizes
# 3) Quality filters
# 4) Runs ANNOVAR annotation
# 5) Filters germline (dbSNP) and coding variants
# 6) Calculates TMB
# 7) Creates visualizations
# 8) Generates summary report
# Key fixes: no circular calls, robust paths, detailed logging, error checks

# Load required packages
if (!requireNamespace("vcfR", quietly=TRUE)) stop("Please install vcfR")
if (!requireNamespace("data.table", quietly=TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
library(vcfR)
library(data.table)
library(dplyr)
library(ggplot2)

#------------------------------
# CONFIGURATION (ADJUST PATHS)
#------------------------------
vcf_file    <- "/home/dhibi/TMB_Pipeline_Project/data/aligned/skcm/SRR26456208_freebayes.vcf"
ref_genome  <- "/home/dhibi/TMB_Pipeline_Project/data/reference/hg38.fa"
humandb_dir <- "/home/dhibi/TMB_Pipeline_Project/data/aligned/skcm/TMB_analysis/annovar/humandb"
work_dir    <- "/home/dhibi/TMB_Pipeline_Project/script"
threads     <- 20

# Create output directories
dirs <- file.path(work_dir, c("filtered_vcfs", "annovar_output", "results", "plots"))
sapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# Logging helper
glog <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

#------------------------------
# STEP 1: LOAD & FILTER VCF
#------------------------------
load_and_filter_vcf <- function(vcf_path, qual=20, dp=20) {
  glog("STEP 1: Loading VCF:", vcf_path)
  if (!file.exists(vcf_path)) stop("VCF not found: ", vcf_path)

  vcf <- tryCatch(read.vcfR(vcf_path), error = function(e) stop("Failed to read VCF: ", e$message))
  fix_df <- as.data.frame(getFIX(vcf))
  fix_df$QUAL <- suppressWarnings(as.numeric(as.character(fix_df$QUAL)))

  # Extract DP from INFO field; assign NA if missing
  info_vec <- as.character(fix_df$INFO)
  dp_vals <- integer(length(info_vec))
  dp_vals[] <- NA_integer_
  has_dp <- grepl("DP=[0-9]+", info_vec)
  dp_vals[has_dp] <- as.integer(sub(".*DP=([0-9]+).*", "\\1", info_vec[has_dp]))
  fix_df$DP_INFO <- dp_vals

  # Apply filters
  filt <- fix_df %>%
    filter(!is.na(QUAL) & QUAL >= qual,
           !is.na(DP_INFO) & DP_INFO >= dp) %>%
    distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)
  glog("Variants after QC filter:", nrow(filt))

  # Store initial count
  attr(filt, "initial_count") <- nrow(fix_df)
  return(filt)
}

#------------------------------
# STEP 2: FIX HEADER & NORMALIZE
#------------------------------
preprocess_vcf <- function(input_vcf, ref, out_dir) {
  sample_id <- tools::file_path_sans_ext(basename(input_vcf))
  contig_file <- file.path(out_dir, paste0(sample_id, "_contigs.txt"))
  fai_file <- paste0(ref, ".fai")

  if (file.exists(fai_file)) {
    fai_lines <- readLines(fai_file)
    ids  <- sub("\\s.*", "", fai_lines)
    lens <- sub(".*\\s", "", fai_lines)
    head_lines <- c("##fileformat=VCFv4.2",
                    sprintf("##contig=<ID=%s,length=%s>", ids, lens))
    writeLines(head_lines, contig_file)
  } else {
    warning("FAI index not found; using basic contigs")
    basic <- sprintf("##contig=<ID=chr%s>", c(1:22, 'X','Y','M'))
    writeLines(c("##fileformat=VCFv4.2", basic), contig_file)
  }

  fixed_vcf <- file.path(out_dir, paste0(sample_id, "_fixed.vcf"))
  bash_cmd <- sprintf(
    "cat %s > %s && grep -v '^##fileformat' %s | grep '^##' >> %s && grep -v '^##' %s >> %s",
    contig_file, fixed_vcf, input_vcf, fixed_vcf, input_vcf, fixed_vcf)
  system2("bash", c("-c", bash_cmd))

  norm_vcf <- file.path(out_dir, paste0(sample_id, "_norm.vcf"))
  system2("bcftools", c("norm", "-m-both", "-f", ref, fixed_vcf, "-Ov"), stdout = norm_vcf)
  glog("Normalized VCF:", norm_vcf)
  return(norm_vcf)
}

#------------------------------
# STEP 3: QUALITY FILTER
#------------------------------
filter_quality <- function(vcf_in, out_dir) {
  sample_id <- tools::file_path_sans_ext(basename(vcf_in))
  out_vcf <- file.path(out_dir, paste0(sample_id, "_qual.vcf"))

  hdr <- system2("bcftools", c("view", "-h", vcf_in), stdout = TRUE)
  expr <- if (any(grepl("##FORMAT=<ID=DP", hdr))) {
    "QUAL>=20 && INFO/DP>=20 && FORMAT/DP[0]>=20"
  } else if (any(grepl("##INFO=<ID=DP", hdr))) {
    "QUAL>=20 && INFO/DP>=20"
  } else {
    "QUAL>=20"
  }

  system2("bcftools", c("filter", "-i", expr, vcf_in, "-Ov", "-o", out_vcf))
  gz <- paste0(out_vcf, ".gz")
  system2("bgzip", c("-c", out_vcf), stdout = gz)
  system2("tabix", c("-p", "vcf", gz))
  glog("Quality-filtered VCF saved:", gz)
  return(gz)
}
#------------------------------
# STEP 4: RUN ANNOVAR
#------------------------------
run_annovar <- function(vcf_gz, humandb, out_dir) {
  sample_id <- tools::file_path_sans_ext(basename(vcf_gz))
  avinput <- file.path(out_dir, paste0(sample_id, ".avinput"))
  system2("convert2annovar.pl", c("-format","vcf4",vcf_gz, ">",avinput))
  prot <- "refGene"; op <- "g"
  if (file.exists(file.path(humandb, "hg38_avsnp150.txt"))) {
    prot <- paste0(prot, ",avsnp150"); op <- paste0(op, ",f")
  }
  out_pref <- file.path(out_dir, sample_id)
  system2("table_annovar.pl",
    c(avinput, humandb, "-buildver","hg38",
      "-out",out_pref,
      "-remove","-protocol",prot,
      "-operation",op,"-nastring",".","-csvout"))
  annofile <- paste0(out_pref, "_hg38_multianno.csv")
  log("ANNOVAR CSV:", annofile)
  return(annofile)
}

#------------------------------
# STEP 5: TMB CALCULATION
#------------------------------
calc_tmb <- function(annovar_csv, exome_size=30) {
  log("STEP 5: Calculating TMB from:", annovar_csv)
  df <- fread(annovar_csv)
  snp_col <- grep("avsnp", names(df), ignore.case=TRUE, value=TRUE)[1]
  df2 <- if (!is.na(snp_col)) df %>% filter(.data[[snp_col]] == '.') else df
  log("After dbSNP filter:", nrow(df2))
  if ('Func.refGene' %in% names(df2)) df2 <- df2 %>% filter(Func.refGene %in% c('exonic','splicing'))
  log("Coding variants:", nrow(df2))
  if ('ExonicFunc.refGene' %in% names(df2)) df2 <- df2 %>% filter(ExonicFunc.refGene != 'synonymous SNV')
  nvar <- nrow(df2)
  score <- round(nvar / exome_size,2)
  class <- if (score<10) 'Low TMB' else if(score<20) 'Intermediate TMB' else 'High TMB'
  log("TMB Score:", score, "mut/Mb (",class,")")
  return(list(score=score, count=nvar, classification=class, df=df2))
}

#------------------------------
# STEP 6: VISUALIZATIONS
#------------------------------
plot_tmb <- function(df, out_dir) {
  if (nrow(df)==0) { warning("No variants to plot"); return() }
  p1 <- df %>% count(Func.refGene) %>%
    ggplot(aes(x=Func.refGene,y=n,fill=Func.refGene))+
    geom_col() + coord_flip() + theme_minimal() + theme(legend.position='none')
  ggsave(file.path(out_dir,"plots/func_dist.png"),p1)
}

#------------------------------
# STEP 7: REPORT
#------------------------------
write_report <- function(res, work_dir) {
  fn <- file.path(work_dir, "results", "tmb_report.txt")
  lines <- c(
    "=== TMB Analysis Report ===",
    paste("Score:",res$score),
    paste("Count:",res$count),
    paste("Class:",res$classification),
    paste("Date:",Sys.time())
  )
  writeLines(lines, fn)
  log("Report written to", fn)
}

#------------------------------
# MAIN EXECUTION
#------------------------------
log("Starting pipeline...")
filtered <- load_and_filter_vcf(vcf_file)
norm_vcf <- preprocess_vcf(vcf_file, ref_genome, file.path(work_dir,"filtered_vcfs"))
qc_vcf <- filter_quality(norm_vcf, file.path(work_dir,"filtered_vcfs"))
anno_csv <- run_annovar(qc_vcf, humandb_dir, file.path(work_dir,"annovar_output"))
res <- calc_tmb(anno_csv)
plot_tmb(res$df, work_dir)
write_report(res, work_dir)
log("Pipeline completed successfully")
