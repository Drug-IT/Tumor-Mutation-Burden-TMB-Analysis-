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
