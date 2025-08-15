#!/bin/bash

# Create README files for data directories
# create_data_readmes.sh

echo "=== Creating Data README Files ==="

# Create data/raw/README.md
echo "Creating raw data README..."
cat << 'EOF' > data/raw/README.md
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
EOF

# Create data/reference/README.md
echo "Creating reference data README..."
cat << 'EOF' > data/reference/README.md
# Reference Genome Installation Instructions

This directory contains instructions for downloading and setting up the human reference genome used in the TMB analysis pipeline.

## 🧬 Reference Genome Information

**Build**: hg38 (GRCh38)  
**Source**: UCSC Genome Browser  
**Assembly**: GRCh38.p14 (latest patch)  
**Size**: ~3.2 GB (compressed), ~3.2 GB (uncompressed)  

## 📥 Download hg38 Reference Genome

### Method 1: Direct Download (Recommended)

```bash
# Navigate to reference directory
cd data/reference/

# Download hg38 reference genome (compressed)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Decompress
gunzip hg38.fa.gz

# Verify download
ls -lh hg38.fa
```

### Method 2: Using rsync (Alternative)

```bash
# Navigate to reference directory
cd data/reference/

# Download using rsync (more reliable for large files)
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ./

# Decompress
gunzip hg38.fa.gz
```

## 🔧 Create BWA Index

The TMB pipeline requires BWA indices for alignment:

```bash
# Create BWA index (takes 30-60 minutes)
bwa index hg38.fa

# This creates the following index files:
# hg38.fa.amb
# hg38.fa.ann  
# hg38.fa.bwt
# hg38.fa.pac
# hg38.fa.sa
```

## 🛠️ Create Additional Indices

### Samtools Index
```bash
# Create .fai index for samtools
samtools faidx hg38.fa
```

### GATK Dictionary
```bash
# Create sequence dictionary for GATK tools
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict

# Or using Picard (if GATK not available)
java -jar picard.jar CreateSequenceDictionary R=hg38.fa O=hg38.dict
```

## 📊 Expected File Structure

After complete setup:
```
data/reference/
├── hg38.fa           # Main reference genome (3.2 GB)
├── hg38.fa.amb       # BWA index file
├── hg38.fa.ann       # BWA index file  
├── hg38.fa.bwt       # BWA index file (~3.2 GB)
├── hg38.fa.fai       # Samtools index
├── hg38.fa.pac       # BWA index file (~800 MB)
├── hg38.fa.sa        # BWA index file (~1.6 GB)
├── hg38.dict         # GATK sequence dictionary
└── README.md         # This file
```

**Total Size**: ~9-10 GB

## ✅ Verification

### Check File Integrity
```bash
# Verify reference genome
head -5 hg38.fa
tail -5 hg38.fa

# Check BWA index
bwa mem hg38.fa 2>&1 | head -10

# Verify samtools index
samtools faidx hg38.fa chr1:1-1000
```

## ⚠️ Important Notes

1. **Disk Space**: Ensure ≥15 GB free space before starting
2. **Memory**: BWA indexing requires ~6 GB RAM
3. **Time**: Complete setup takes 1-2 hours depending on hardware
4. **Consistency**: Always use the same reference version throughout your pipeline

## 🔗 Pipeline Integration

This reference genome integrates with:
- `scripts/bwa_alignment_gui.R` - BWA-MEM alignment
- `scripts/preprocess_bam_for_tmb.R` - BAM preprocessing  
- `scripts/tmb_pipeline.R` - Variant calling and annotation

## 📚 Additional Resources

- [UCSC hg38 Information](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38)
- [GRCh38 Release Notes](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)
- [BWA Manual](http://bio-bwa.sourceforge.net/bwa.shtml)
EOF

echo "✓ Created data/raw/README.md"
echo "✓ Created data/reference/README.md"
echo ""
echo "README files created successfully!"
echo "These provide detailed instructions for:"
echo "  - Downloading SRR26456208 melanoma WES data"
echo "  - Installing hg38 reference genome and indices"
echo "  - Quality control and verification steps"
