# Tumor-Mutation-Burden-TMB-Analysis-
 Tumor Mutation Burden (TMB) measures how many mutations exist in a tumor’s DNA. Higher TMB can indicate that a patient might respond better to immunotherapy. This project analyze mutation data, annotate which mutations affect proteins, and calculate TMB by counting important mutations in protein-coding regions all using R. 
# Tumor Mutation Burden (TMB) Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive end-to-end pipeline for calculating Tumor Mutation Burden (TMB) from Whole Exome Sequencing (WES) data, with specialized optimization for cutaneous melanoma analysis.

## 🔬 What This Pipeline Does

- **Alignment**: BWA-MEM alignment with interactive Shiny GUI
- **Variant Calling**: FreeBayes-based variant detection
- **Annotation**: ANNOVAR-powered functional annotation
- **Filtering**: Multi-tier somatic variant filtering with VAF correction
- **TMB Calculation**: Accurate TMB scoring with melanoma-specific thresholds
- **Quality Control**: Comprehensive QC metrics and UV signature analysis
- **Reporting**: Automated plots and detailed analysis reports

## 🚀 Quick Start

### Prerequisites

```bash
# Core tools (install via conda/mamba)
conda install -c bioconda bwa samtools bcftools freebayes

# R packages
Rscript -e "install.packages(c('shiny', 'dplyr', 'readr', 'ggplot2', 'maftools'))"

# ANNOVAR (requires license - download separately)
# Download from: https://annovar.openbioinformatics.org/
Installation
bashgit clone git@github.com:Drug-IT/Tumor-Mutation-Burden-TMB-Analysis-.git
cd Tumor-Mutation-Burden-TMB-Analysis-
chmod +x scripts/install_requirements.sh
./scripts/install_requirements.sh
Run Pipeline
bash# 1. Alignment (interactive GUI)
Rscript scripts/bwa_alignment_gui.R
```

## Full TMB analysis from VCF
```
Rscript scripts/tmb_pipeline.R
📁 Project Structure
├── scripts/                    # All executable scripts
│   ├── bwa_alignment_gui.R    # Interactive alignment GUI
│   ├── tmb_pipeline.R         # Main TMB analysis pipeline
│   ├── tmb_accuracy_metrics.R # QC and validation metrics
│   ├── uv_signature_analysis.R# UV signature detection
│   └── ...                    # Additional analysis scripts
├── data/
│   ├── raw/                   # Input FASTQ files
│   └── processed/             # Normalized VCFs, filtered data
├── results/
│   ├── tables/                # CSV results and statistics
│   ├── figures/               # Generated plots and visualizations
│   └── reports/               # HTML reports and summaries
├── requirements.txt           # Python dependencies
└── README.md                  # This file

```
## Sample Analysis
**Validation Sample:** SRR26456208 (cutaneous melanoma WES)

**Reference Genome:** hg38

**TMB Result:** Calculated with VAF correction and somatic filtering

## Key Features
#### Advanced VAF Correction

Corrects unreliable VCF AF values using FORMAT/AD depths

Ensures biologically sound VAF filtering (10-90%)

Melanoma-Specific Analysis

UV signature detection (C>T transitions)

Melanoma-appropriate TMB thresholds

Specialized somatic variant filtering

## Quality Control

Comprehensive filtering metrics

Coverage and quality distribution analysis

Population frequency validation

## Configuration
```
Edit the configuration block in scripts/tmb_pipeline.R:
rSAMPLE_ID <- "your_sample_id"
VCF_FILE  <- "path/to/normalized.vcf"
ANNOVAR_DIR <- "path/to/annovar"
BUILD <- "hg38"
```

## Results
The pipeline generates:

**TMB Score:** Mutations per megabase with confidence intervals

**Filtered Variants:** High-quality somatic mutations table

**QC Plots:** VAF distributions, filtering waterfalls, mutation signatures

**Comprehensive Report:** Detailed analysis summary

## Validation
Validated on melanoma WES data with:

✅ VAF parsing accuracy verification

✅ Population frequency filtering validation

✅ UV signature confirmation for melanoma samples

✅ TMB concordance with expected ranges

## Citation
If you use this pipeline, please cite:
TMB Analysis Pipeline for Cutaneous Melanoma
Dorra Dhibi, Drug-IT Startup
https://github.com/Drug-IT/Tumor-Mutation-Burden-TMB-Analysis-

## Author
Dorra Dhibi
Drug-IT Startup
📧 Contact: GitHub Issues
📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
🤝 Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

For questions or issues, please open a GitHub issue or contact Drug-IT startup.
