#!/bin/bash

# TMB Pipeline Organization Script
# Run this from your project root directory

echo "=== TMB Pipeline Organization ==="

# Set the source directory
SOURCE_DIR="/home/dhibi/TMB_Pipeline_Project/script/Final/DrugIT"

# Check if source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory $SOURCE_DIR not found!"
    echo "Current directory: $(pwd)"
    echo "Please run this script from your project root or update SOURCE_DIR"
    exit 1
fi

echo "Source directory found: $SOURCE_DIR"

# 1. Create the target directory structure
echo "Creating directory structure..."
mkdir -p scripts data/raw data/processed data/reference results/figures results/tables results/reports

# 2. Move R scripts to scripts/ directory (check if files exist first)
echo "Moving R scripts..."
R_SCRIPTS=("bwa_alignment_gui.R" "preprocess_bam_for_tmb.R" "tmb_pipeline.R" "tmb_accuracy_metrics.R" "uv_signature_analysis.R" "tmb_verification.R" "fix_vaf_parsing.R" "mark_duplicates.R" "bqsr.R")

for script in "${R_SCRIPTS[@]}"; do
    if [ -f "$SOURCE_DIR/$script" ]; then
        mv "$SOURCE_DIR/$script" scripts/
        echo "  ✓ Moved $script"
    else
        echo "  ⚠ Not found: $script"
    fi
done

# 3. Move result files to appropriate locations
echo "Moving result files..."
RESULT_FILES=("final_corrected_tmb_results.csv" "corrected_tmb_results.csv" "tmb_verification_summary.csv" "filtered_somatic_mutations.csv")

for file in "${RESULT_FILES[@]}"; do
    if [ -f "$SOURCE_DIR/$file" ]; then
        mv "$SOURCE_DIR/$file" results/tables/
        echo "  ✓ Moved $file to results/tables/"
    else
        echo "  ⚠ Not found: $file"
    fi
done

# Move reports
REPORT_FILES=("snpEff_summary.html" "snpEff_genes.txt")
for file in "${REPORT_FILES[@]}"; do
    if [ -f "$SOURCE_DIR/$file" ]; then
        mv "$SOURCE_DIR/$file" results/reports/
        echo "  ✓ Moved $file to results/reports/"
    else
        echo "  ⚠ Not found: $file"
    fi
done

# 4. Move FastQC reports to results
echo "Moving QC reports..."
QC_FILES=("SRR26456208_1_fastqc.html" "SRR26456208_1_fastqc.zip")
for file in "${QC_FILES[@]}"; do
    if [ -f "$SOURCE_DIR/$file" ]; then
        mv "$SOURCE_DIR/$file" results/reports/
        echo "  ✓ Moved $file to results/reports/"
    else
        echo "  ⚠ Not found: $file"
    fi
done

# 5. Move processed data
echo "Moving processed data..."
if [ -f "$SOURCE_DIR/SRR26456208_normalized.vcf" ]; then
    mv "$SOURCE_DIR/SRR26456208_normalized.vcf" data/processed/
    echo "  ✓ Moved SRR26456208_normalized.vcf to data/processed/"
else
    echo "  ⚠ Not found: SRR26456208_normalized.vcf"
fi

# 6. Move large files to excluded directory
echo "Moving large files to excluded directory..."
mkdir -p .gitignore_excluded

# Move BAM files
find "$SOURCE_DIR" -name "*.bam*" -exec mv {} .gitignore_excluded/ \; 2>/dev/null
find "$SOURCE_DIR" -name "*.bai" -exec mv {} .gitignore_excluded/ \; 2>/dev/null

# 7. Clean up RStudio and temp files
echo "Cleaning up temporary files..."
rm -rf "$SOURCE_DIR/.Rproj.user/" 2>/dev/null || true
rm -f "$SOURCE_DIR/.RData" 2>/dev/null || true
rm -f "$SOURCE_DIR/.Rhistory" 2>/dev/null || true
rm -f "$SOURCE_DIR/*.Rproj" 2>/dev/null || true

echo "=== Organization Complete ==="
echo "Files organized into:"
echo "  scripts/     - All R scripts (flat structure)"
echo "  data/        - Raw and processed data"  
echo "  results/     - Tables, figures, and reports"
echo "  .gitignore_excluded/ - Large files excluded from git"
echo ""
echo "Remaining files in source directory:"
ls -la "$SOURCE_DIR" 2>/dev/null || echo "Source directory is empty or removed"
