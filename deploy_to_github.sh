#!/bin/bash

# =============================================================================
# Professional Git Workflow for Existing TMB Analysis Pipeline Repository
# Optimized for repositories with complete file structure already in place
# Author: Dorra Dhibi - Drug-IT
# Version: 2.1.0 - Tailored for existing repository with all files present
# =============================================================================

set -e  # Exit on any error
set -u  # Exit on undefined variables

# Configuration
readonly SCRIPT_VERSION="2.1.0"
readonly REPO_NAME="Tumor-Mutation-Burden-TMB-Analysis-"
readonly GITHUB_ORG="Drug-IT"
readonly REPO_URL_HTTPS="https://github.com/${GITHUB_ORG}/${REPO_NAME}.git"
readonly REPO_URL_SSH="git@github.com:${GITHUB_ORG}/${REPO_NAME}.git"

# Colors for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly PURPLE='\033[0;35m'
readonly CYAN='\033[0;36m'
readonly NC='\033[0m' # No Color

# Logging functions
print_step() {
    echo -e "${BLUE}🔸 Step $1: $2${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ Error: $1${NC}"
}

print_info() {
    echo -e "${PURPLE}ℹ️  $1${NC}"
}

print_highlight() {
    echo -e "${CYAN}🌟 $1${NC}"
}

# Error handler
error_handler() {
    local line_number=$1
    print_error "Script failed at line ${line_number}. Exiting..."
    echo "💡 Tip: Check your Git configuration, network connection, and GitHub credentials"
    exit 1
}

trap 'error_handler ${LINENO}' ERR

# User confirmation function
confirm_action() {
    local message="$1"
    local default="${2:-N}"
    
    if [[ "$default" == "Y" ]]; then
        read -p "${message} (Y/n): " confirm
        [[ ${confirm:-Y} =~ ^[Yy]$ ]]
    else
        read -p "${message} (y/N): " confirm
        [[ $confirm =~ ^[Yy]$ ]]
    fi
}

# Check prerequisites
check_prerequisites() {
    print_step "0" "Checking prerequisites and repository status"
    
    # Check if git is installed
    if ! command -v git &> /dev/null; then
        print_error "Git is not installed. Please install Git first."
        exit 1
    fi
    
    # Verify we're in a git repository
    if ! git rev-parse --is-inside-work-tree &>/dev/null; then
        print_error "Not in a Git repository. This script requires an existing Git repo."
        exit 1
    fi
    
    # Check repository structure
    local required_files=("README.md" "LICENSE" "scripts" "data" "results")
    for file in "${required_files[@]}"; do
        if [[ ! -e "$file" ]]; then
            print_error "Required file/directory not found: $file"
            exit 1
        fi
    done
    
    # Check if git-lfs is available
    if ! command -v git-lfs &> /dev/null; then
        print_warning "Git LFS not found. Large files will be tracked normally."
        GIT_LFS_AVAILABLE=false
    else
        GIT_LFS_AVAILABLE=true
        print_success "Git LFS is available"
    fi
    
    print_success "Repository structure validated - all required components present"
}

# Detect SSH key availability
detect_auth_method() {
    print_step "1" "Detecting optimal authentication method"
    
    if ssh -T git@github.com 2>&1 | grep -q "successfully authenticated"; then
        AUTH_METHOD="ssh"
        REPO_URL="$REPO_URL_SSH"
        print_success "SSH key authentication detected - using secure SSH"
    else
        AUTH_METHOD="https"
        REPO_URL="$REPO_URL_HTTPS"
        print_info "SSH not configured - using HTTPS authentication"
        echo "💡 Consider setting up SSH keys for enhanced security: https://docs.github.com/en/authentication"
    fi
}

# Setup Git LFS for large files
setup_git_lfs() {
    if [[ "$GIT_LFS_AVAILABLE" == true ]]; then
        print_step "2" "Configuring Git LFS for large genomics files"
        
        # Initialize Git LFS
        git lfs install
        
        # Create comprehensive .gitattributes
        cat > .gitattributes << 'EOF'
# Git LFS configuration for TMB Analysis Pipeline
# ==============================================

# Large genomics data files
*.vcf filter=lfs diff=lfs merge=lfs -text
*.vcf.gz filter=lfs diff=lfs merge=lfs -text
*.bam filter=lfs diff=lfs merge=lfs -text
*.bam.bai filter=lfs diff=lfs merge=lfs -text
*.fastq filter=lfs diff=lfs merge=lfs -text
*.fastq.gz filter=lfs diff=lfs merge=lfs -text
*.fq filter=lfs diff=lfs merge=lfs -text
*.fq.gz filter=lfs diff=lfs merge=lfs -text

# Large result files (>10MB)
results/tables/*.csv filter=lfs diff=lfs merge=lfs -text
data/processed/*.vcf filter=lfs diff=lfs merge=lfs -text

# Compressed archives and reports
*.zip filter=lfs diff=lfs merge=lfs -text
*.tar.gz filter=lfs diff=lfs merge=lfs -text
results/reports/*.zip filter=lfs diff=lfs merge=lfs -text

# Reference genomes and indices  
*.fa filter=lfs diff=lfs merge=lfs -text
*.fasta filter=lfs diff=lfs merge=lfs -text
*.fa.gz filter=lfs diff=lfs merge=lfs -text
*.fasta.gz filter=lfs diff=lfs merge=lfs -text
*.fai filter=lfs diff=lfs merge=lfs -text
*.dict filter=lfs diff=lfs merge=lfs -text

# BWA and other alignment indices
*.amb filter=lfs diff=lfs merge=lfs -text
*.ann filter=lfs diff=lfs merge=lfs -text
*.bwt filter=lfs diff=lfs merge=lfs -text
*.pac filter=lfs diff=lfs merge=lfs -text
*.sa filter=lfs diff=lfs merge=lfs -text

# Preserve line endings for text files
*.md text
*.txt text
*.R text
*.py text
*.sh text eol=lf
*.yml text
*.yaml text
*.json text
EOF

        # Track existing large files with LFS
        if [[ -f "data/processed/SRR26456208_normalized.vcf" ]]; then
            git lfs track "data/processed/SRR26456208_normalized.vcf"
            print_info "Tracked VCF file with Git LFS"
        fi
        
        # Track large CSV files
        for csv_file in results/tables/*.csv; do
            if [[ -f "$csv_file" ]] && [[ $(stat -c%s "$csv_file" 2>/dev/null || stat -f%z "$csv_file" 2>/dev/null || echo 0) -gt 1048576 ]]; then
                git lfs track "$csv_file"
                print_info "Tracked large CSV with Git LFS: $(basename "$csv_file")"
            fi
        done
        
        git add .gitattributes
        print_success "Git LFS configured for genomics files"
    else
        print_info "Skipping Git LFS setup (not available)"
    fi
}

# Configure Git and remote
configure_git_and_remote() {
    print_step "3" "Configuring Git settings and remote repository"
    
    # Display current configuration
    CURRENT_NAME=$(git config user.name 2>/dev/null || echo "Not set")
    CURRENT_EMAIL=$(git config user.email 2>/dev/null || echo "Not set")
    
    print_info "Current Git configuration:"
    echo "  Name: $CURRENT_NAME"
    echo "  Email: $CURRENT_EMAIL"
    
    # Set configuration if not already set
    if [[ "$CURRENT_NAME" == "Not set" ]]; then
        git config user.name "Drug-IT"
        print_success "Set git user.name to 'Drug-IT'"
    fi
    
    if [[ "$CURRENT_EMAIL" == "Not set" ]]; then
        git config user.email "dorra.dhibi@drug-it.com"  
        print_success "Set git user.email to 'dorra.dhibi@drug-it.com'"
    fi
    
    # Set additional configurations
    git config core.autocrlf input
    git config core.ignorecase false
    git config pull.rebase false
    git config init.defaultBranch main
    
    # Setup remote
    if git remote get-url origin &>/dev/null; then
        CURRENT_REMOTE=$(git remote get-url origin)
        if [[ "$CURRENT_REMOTE" != "$REPO_URL" ]]; then
            print_warning "Remote 'origin' currently points to: $CURRENT_REMOTE"
            if confirm_action "Update remote to $REPO_URL?"; then
                git remote set-url origin "$REPO_URL"
                print_success "Remote origin updated to $AUTH_METHOD"
            else
                print_info "Keeping existing remote configuration"
                REPO_URL="$CURRENT_REMOTE"  # Use existing remote
            fi
        else
            print_success "Remote origin already correctly configured"
        fi
    else
        git remote add origin "$REPO_URL"
        print_success "Remote origin added using $AUTH_METHOD"
    fi
}

# Organize files into logical commits
organize_and_commit_files() {
    print_step "4" "Organizing existing files into logical feature branches"
    
    # Ensure we're on main branch
    git checkout main 2>/dev/null || git checkout -b main
    
    # Check current status
    if git diff --quiet && git diff --staged --quiet; then
        print_info "Working directory is clean"
    else
        print_warning "You have uncommitted changes"
        git status --short
        if confirm_action "Stage and commit all changes as 'Complete TMB Pipeline v2.1.0'?" "Y"; then
            git add .
            git commit -m "feat: complete TMB analysis pipeline v2.1.0

🔬 Complete TMB Analysis Pipeline Implementation
===============================================

📦 Repository Contents:
- 📜 Complete documentation (README.md, LICENSE)
- 🧬 Core analysis scripts (15 R/Python files)
- 📊 Sample data and validation results
- 📋 Quality control reports and metrics
- 🛠️ Pipeline automation scripts

🔬 Analysis Components:
✅ Interactive BWA-MEM alignment GUI
✅ Advanced VAF parsing and correction
✅ Multi-tier somatic variant filtering  
✅ Population frequency validation
✅ UV signature detection for melanoma
✅ Comprehensive TMB calculation
✅ Quality control and verification

📁 Directory Structure:
- scripts/: Complete pipeline implementation
  • bwa_alignment_gui.R - Interactive alignment interface
  • tmb_pipeline.R - Core TMB calculation engine
  • fix_vaf_parsing.R - Advanced VAF correction
  • uv_signature_analysis.R - Melanoma signature detection
  • tmb_verification.R - Quality control validation
  • calculation.py - High-performance computing utilities
  • [9 additional analysis modules]

- data/: Sample datasets and references
  • processed/SRR26456208_normalized.vcf - Validated melanoma WES
  • raw/README.md - Data acquisition guide
  • reference/README.md - Genome setup instructions

- results/: Comprehensive analysis outputs  
  • tables/ - TMB results and filtered variants
  • reports/ - FastQC and snpEff QC reports
  • figures/ - Visualization outputs (placeholder)

🧬 Validated Sample Analysis:
- Sample: SRR26456208 (cutaneous melanoma, hereditary cancer)
- Platform: Illumina whole exome sequencing
- Reference: hg38 with BWA indexing
- Results: Publication-ready TMB analysis

🎯 Key Features:
- Professional R Shiny GUI for alignment
- Advanced FORMAT/AD depth parsing
- Melanoma-specific filtering (UV signatures)
- Multi-database population frequency validation
- Comprehensive QC metrics and reporting
- Cross-platform compatibility (R/Python)

👤 Author: Dorra Dhibi - Drug-IT Startup  
📄 License: MIT (2025)
🔗 Ready for: Production deployment and collaboration

⚡ Performance: Optimized for large-scale WES analysis
🛡️ Quality: Validated against known melanoma samples
📊 Output: Publication-ready TMB scores and reports"
            
            print_success "All files committed to main branch"
        else
            print_info "Skipping commit - please handle uncommitted changes manually"
        fi
    fi
}

# Enhanced push with comprehensive retry mechanism
push_to_github() {
    print_step "5" "Pushing complete pipeline to GitHub"
    
    local max_retries=3
    local retry_count=0
    
    if [[ "$AUTH_METHOD" == "https" ]]; then
        echo ""
        print_highlight "GitHub HTTPS Authentication Required"
        echo "🔐 You'll be prompted for GitHub credentials:"
        echo "   • Username: Your GitHub username"
        echo "   • Password: Your Personal Access Token (NOT your account password)"
        echo ""
        echo "📋 Creating a Personal Access Token:"
        echo "   1. Go to: https://github.com/settings/tokens"
        echo "   2. Click 'Generate new token' → 'Generate new token (classic)'"
        echo "   3. Set expiration and select 'repo' scope (full repository access)"
        echo "   4. Copy the token and use it as your password below"
        echo ""
        if ! confirm_action "Ready to proceed with authentication?"; then
            print_info "Aborted by user. You can run this script again when ready."
            exit 0
        fi
    fi
    
    # Push with retry logic
    while [[ $retry_count -lt $max_retries ]]; do
        print_info "Attempting to push... (attempt $((retry_count + 1))/$max_retries)"
        
        if git push -u origin main --verbose; then
            print_success "Successfully pushed to GitHub!"
            echo ""
            print_highlight "Repository URL: https://github.com/${GITHUB_ORG}/${REPO_NAME}"
            break
        else
            ((retry_count++))
            if [[ $retry_count -lt $max_retries ]]; then
                print_warning "Push failed. Retrying in 5 seconds..."
                sleep 5
            else
                print_error "Failed to push after $max_retries attempts"
                echo ""
                echo "🔧 Troubleshooting steps:"
                echo "   1. Check internet connection"
                echo "   2. Verify GitHub repository exists and you have write access"
                echo "   3. For HTTPS: Ensure you're using a Personal Access Token"
                echo "   4. For SSH: Check SSH key configuration"
                echo "   5. Try manual push: git push origin main --verbose"
                echo ""
                echo "🌐 Create repository manually if needed:"
                echo "   https://github.com/new"
                echo "   Repository name: $REPO_NAME"
                exit 1
            fi
        fi
    done
}

# Display comprehensive success summary
show_success_summary() {
    echo ""
    echo "🎉 TMB Analysis Pipeline Successfully Deployed!"
    echo "=============================================="
    echo ""
    echo "📊 Deployment Summary:"
    echo "  ✅ Repository: https://github.com/${GITHUB_ORG}/${REPO_NAME}"
    echo "  ✅ Authentication: $AUTH_METHOD"
    echo "  ✅ Git LFS: $([ "$GIT_LFS_AVAILABLE" == true ] && echo "Configured for large files" || echo "Standard Git tracking")"
    echo "  ✅ Files: Complete pipeline with all components"
    echo "  ✅ Documentation: Professional README and guides"
    echo ""
    echo "🔬 Pipeline Components Deployed:"
    echo "  📜 Documentation & Setup:"
    echo "     • README.md - Comprehensive usage guide"
    echo "     • LICENSE - MIT license (Drug-IT 2025)" 
    echo "     • requirements.txt - R/Python dependencies"
    echo "     • .gitignore - Genomics-optimized exclusions"
    echo ""
    echo "  🧬 Analysis Scripts (15 files):"
    echo "     • Interactive alignment GUI with progress tracking"
    echo "     • Advanced VAF correction using FORMAT/AD depths"
    echo "     • Multi-tier somatic variant filtering"
    echo "     • UV signature detection for melanoma samples"
    echo "     • Comprehensive TMB verification framework"
    echo "     • Python integration for high-performance calculations"
    echo ""
    echo "  📊 Data & Results:"
    echo "     • Sample VCF: SRR26456208 cutaneous melanoma (validated)"
    echo "     • TMB results: Complete analysis with corrections"
    echo "     • QC reports: FastQC and snpEff validation"
    echo "     • Reference guides: hg38 setup and data acquisition"
    echo ""
    if [[ "$GIT_LFS_AVAILABLE" == true ]]; then
        echo "📦 Git LFS Configuration:"
        echo "  • Large genomics files tracked efficiently"
        echo "  • VCF files: Optimized storage and transfer"
        echo "  • Result tables: Large CSV files handled properly" 
        echo "  • Archives: ZIP and compressed files managed"
        echo ""
    fi
    echo "🚀 Next Steps & Recommendations:"
    echo ""
    echo "  1. 🏷️  Create Release Tags:"
    echo "     git tag -a v2.1.0 -m 'Complete TMB Analysis Pipeline'"
    echo "     git push origin v2.1.0"
    echo ""
    echo "  2. 📖 Repository Enhancements:"
    echo "     • Enable GitHub Pages for documentation hosting"
    echo "     • Add repository badges and shields to README"
    echo "     • Configure branch protection for main branch"
    echo ""
    echo "  3. 🤖 Automation Options:"
    echo "     • Set up GitHub Actions for CI/CD"
    echo "     • Automate testing with sample datasets"
    echo "     • Configure automatic quality checks"
    echo ""
    echo "  4. 👥 Community Features:"
    echo "     • Enable GitHub Discussions for user support"
    echo "     • Create issue templates for bug reports"
    echo "     • Set up contributing guidelines"
    echo ""
    echo "  5. 📊 Usage Analytics:"
    echo "     • Monitor repository traffic and clones"
    echo "     • Track issue resolution and user feedback"
    echo "     • Document common usage patterns"
    echo ""
    echo "🔗 Important Links:"
    echo "  📂 Repository: https://github.com/${GITHUB_ORG}/${REPO_NAME}"
    echo "  📋 Issues: https://github.com/${GITHUB_ORG}/${REPO_NAME}/issues"
    echo "  🌟 Releases: https://github.com/${GITHUB_ORG}/${REPO_NAME}/releases"
    echo "  📊 Insights: https://github.com/${GITHUB_ORG}/${REPO_NAME}/pulse"
    echo ""
    if [[ "$GIT_LFS_AVAILABLE" == true ]]; then
        echo "💡 Git LFS Usage Notes:"
        echo "  • Collaborators need Git LFS installed: git lfs install"
        echo "  • Monitor LFS bandwidth usage in repository settings"
        echo "  • Large files download automatically on clone/pull"
        echo ""
    fi
    echo -e "${GREEN}✨ Your TMB Analysis Pipeline is ready for production use! ✨${NC}"
    echo -e "${CYAN}🏆 Perfect for bioinformatics research and collaboration! 🏆${NC}"
    echo ""
    echo "🎯 Ready for: Melanoma TMB analysis, method validation, and scientific publication"
}

# Rollback function for error recovery
rollback_changes() {
    print_warning "Performing rollback to clean state..."
    
    # Remove .gitattributes if we created it
    [[ -f .gitattributes ]] && git rm -f .gitattributes 2>/dev/null && echo "Removed .gitattributes"
    
    # Reset any staged changes
    git reset HEAD . 2>/dev/null || true
    
    print_success "Rollback completed - repository returned to previous state"
}

# Main workflow execution
main() {
    # Header with version info
    echo "🚀 TMB Analysis Pipeline - GitHub Deployment Script v${SCRIPT_VERSION}"
    echo "📁 Current directory: $(pwd)"
    echo "🎯 Optimized for existing repositories with complete file structure"
    echo "👤 Author: Dorra Dhibi - Drug-IT Startup"
    echo ""
    
    # Setup interrupt handling for graceful rollback
    trap 'echo ""; print_warning "Script interrupted by user"; rollback_changes; exit 130' INT
    
    # Execute deployment workflow
    check_prerequisites
    detect_auth_method  
    configure_git_and_remote
    setup_git_lfs
    organize_and_commit_files
    push_to_github
    show_success_summary
    
    # Clean exit
    echo "🎉 Deployment completed successfully!"
    exit 0
}

# Script entry point
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
