# Raw Data Download Instructions

This directory contains instructions for downloading the raw sequencing data used in the TMB analysis pipeline.

## 📊 Dataset Information

**Study**: WES of patient with hereditary cancer (melanoma)  
**Database**: SRA - NCBI  
**Primary Sample**: `SRR26456208` (cutaneous melanoma)  
**Secondary Sample**: `SRR26456207` (used for alignment artifacts)  
**Data Type**: Whole Exome Sequencing (WES)  
**Platform**: Illumina  

## 🔗 Data Source

**NCBI SRA Link**: [WES of patient with hereditary cancer (melanoma) - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra)

## 💾 Download Raw FASTQ Files

### Prerequisites
```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Or using apt (Ubuntu/Debian)
sudo apt-get install sra-toolkit
```

### Download Commands

```bash
# Navigate to raw data directory
cd data/raw/

# Download primary sample (SRR26456208)
fastq-dump --split-files --gzip SRR26456208

# Download secondary sample (SRR26456207) - optional
fastq-dump --split-files --gzip SRR26456207

# Verify downloads
ls -lh *.fastq.gz
```

### Expected Output Files
```
SRR26456208_1.fastq.gz  # Forward reads
SRR26456208_2.fastq.gz  # Reverse reads
SRR26456207_1.fastq.gz  # Forward reads (optional)
SRR26456207_2.fastq.gz  # Reverse reads (optional)
```

## 🧬 Alternative Download Method (faster)

```bash
# Using prefetch + fasterq-dump for better performance
prefetch SRR26456208
fasterq-dump --split-files SRR26456208
gzip *.fastq

# Clean up temporary files
rm -rf SRR26456208/
```

## 📏 File Size Information

**Approximate sizes**:
- `SRR26456208_1.fastq.gz`: ~2-4 GB
- `SRR26456208_2.fastq.gz`: ~2-4 GB
- **Total**: ~4-8 GB per sample

⚠️ **Storage Warning**: Ensure you have sufficient disk space (≥20 GB recommended) before downloading.

## ✅ Quality Control

After download, verify file integrity:
```bash
# Check file sizes
ls -lh *.fastq.gz

# Quick read count verification
zcat SRR26456208_1.fastq.gz | echo $((`wc -l`/4))
zcat SRR26456208_2.fastq.gz | echo $((`wc -l`/4))
```

## 📋 Next Steps

After downloading:
1. Run quality control: `fastqc *.fastq.gz`
2. Proceed with alignment using `scripts/bwa_alignment_gui.R`
3. Continue with TMB analysis pipeline

## 🔗 Related Data

- **Reference Genome**: See `data/reference/README.md` for hg38 download instructions
- **Processed Data**: Check `data/processed/` for intermediate files
