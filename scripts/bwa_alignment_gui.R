# BWA-MEM Alignment Interactive GUI (Complete Version)
# Author: Your Name
# Date: 2025
# 
# HOW TO RUN THIS SCRIPT:
# 1. Save this file as "bwa_alignment_gui.R"
# 2. In R/RStudio, run: source("bwa_alignment_gui.R")
# 3. Or from command line: Rscript bwa_alignment_gui.R
# 4. The GUI will open in your web browser

# =============================================================================
# PACKAGE INSTALLATION AND LOADING
# =============================================================================

# Function to safely install and load packages
install_and_load <- function(packages) {
  success_status <- list()
  
  for (pkg in packages) {
    tryCatch({
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("Installing", pkg, "...\n")
        install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org/")
        success <- require(pkg, character.only = TRUE, quietly = TRUE)
        success_status[[pkg]] <- success
        if (success) {
          cat("✅", pkg, "installed and loaded successfully\n")
        } else {
          cat("❌", pkg, "installation failed\n")
        }
      } else {
        success_status[[pkg]] <- TRUE
        cat("✅", pkg, "already available\n")
      }
    }, error = function(e) {
      cat("❌ Error with package", pkg, ":", e$message, "\n")
      success_status[[pkg]] <- FALSE
    })
  }
  
  return(success_status)
}

# Install required packages
cat("=== PACKAGE INSTALLATION ===\n")
required_packages <- c("shiny")
optional_packages <- c("shinyFiles")

# Install essential packages
package_status <- install_and_load(required_packages)

# Check if shiny loaded successfully
if (!package_status[["shiny"]]) {
  stop("❌ Failed to install/load shiny package. Please install manually with:\n   install.packages('shiny')")
}

# Try to load optional packages
cat("\n=== OPTIONAL PACKAGES ===\n")
use_shinyfiles <- FALSE
tryCatch({
  if (require("shinyFiles", quietly = TRUE)) {
    use_shinyfiles <- TRUE
    cat("✅ shinyFiles loaded - Enhanced file browser available\n")
  } else {
    cat("⚠️  Attempting to install shinyFiles...\n")
    install.packages("shinyFiles", dependencies = TRUE, repos = "https://cloud.r-project.org/")
    if (require("shinyFiles", quietly = TRUE)) {
      use_shinyfiles <- TRUE
      cat("✅ shinyFiles installed and loaded successfully\n")
    } else {
      cat("❌ shinyFiles installation failed\n")
    }
  }
}, error = function(e) {
  cat("❌ Error loading shinyFiles:", e$message, "\n")
})

if (!use_shinyfiles) {
  cat("ℹ️  Using basic file input mode (manual path entry)\n")
  cat("   To enable file browser, install shinyFiles:\n")
  cat("   install.packages('shinyFiles')\n")
}

# Load packages
library(shiny)
if (use_shinyfiles) {
  library(shinyFiles)
}

# =============================================================================
# CORE ALIGNMENT FUNCTIONS
# =============================================================================

check_alignment_tools <- function() {
  tools <- c("bwa", "samtools")
  tool_status <- list()
  
  for (tool in tools) {
    check <- try(system(paste("which", tool), intern = TRUE, ignore.stderr = TRUE), silent = TRUE)
    available <- !inherits(check, "try-error") && length(check) > 0
    tool_status[[tool]] <- available
  }
  
  missing_tools <- names(tool_status)[!unlist(tool_status)]
  
  if (length(missing_tools) > 0) {
    return(list(success = FALSE, missing = missing_tools, all_status = tool_status))
  }
  
  return(list(success = TRUE, missing = character(0), all_status = tool_status))
}

index_reference_genome <- function(reference_path) {
  if (!file.exists(reference_path)) {
    return(list(success = FALSE, message = paste("Reference genome not found:", reference_path)))
  }
  
  # Check if already indexed
  index_files <- paste0(reference_path, c(".amb", ".ann", ".bwt", ".pac", ".sa"))
  
  if (all(file.exists(index_files))) {
    return(list(success = TRUE, message = "Reference genome already indexed"))
  }
  
  index_cmd <- paste("bwa index", shQuote(reference_path))
  result <- system(index_cmd)
  
  if (result == 0) {
    return(list(success = TRUE, message = "BWA index created successfully"))
  } else {
    return(list(success = FALSE, message = "Failed to create BWA index"))
  }
}

run_bwa_mem <- function(sample_id, fastq1, fastq2, reference_path, 
                        output_dir, threads = 4, platform = "ILLUMINA") {
  
  # Validate inputs
  if (!file.exists(fastq1)) {
    return(list(success = FALSE, message = paste("Forward reads not found:", fastq1)))
  }
  
  if (!file.exists(fastq2)) {
    return(list(success = FALSE, message = paste("Reverse reads not found:", fastq2)))
  }
  
  if (!file.exists(reference_path)) {
    return(list(success = FALSE, message = paste("Reference genome not found:", reference_path)))
  }
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Output file path
  output_sam <- file.path(output_dir, paste0(sample_id, "_aligned.sam"))
  
  # Create read group information
  read_group <- sprintf("@RG\\tID:%s\\tSM:%s\\tPL:%s\\tLB:%s\\tPU:%s", 
                        sample_id, sample_id, platform, sample_id, sample_id)
  
  # Build BWA-MEM command
  bwa_cmd <- paste(
    "bwa mem",
    "-t", threads,
    "-M",
    "-R", shQuote(read_group),
    shQuote(reference_path),
    shQuote(fastq1),
    shQuote(fastq2),
    ">", shQuote(output_sam)
  )
  
  # Execute alignment
  start_time <- Sys.time()
  result_code <- system(bwa_cmd)
  end_time <- Sys.time()
  
  if (result_code == 0 && file.exists(output_sam)) {
    elapsed_time <- round(difftime(end_time, start_time, units = "mins"), 2)
    file_size <- file.info(output_sam)$size / 1024^2
    message <- sprintf("BWA-MEM completed in %.2f minutes (%.1f MB)", elapsed_time, file_size)
    return(list(success = TRUE, message = message, output_file = output_sam))
  }
  
  return(list(success = FALSE, message = "BWA-MEM alignment failed"))
}

sam_to_sorted_bam <- function(sam_file, output_dir, sample_id, 
                              threads = 4, keep_sam = FALSE) {
  
  if (!file.exists(sam_file)) {
    return(list(success = FALSE, message = paste("SAM file not found:", sam_file)))
  }
  
  # Output file paths
  temp_bam <- file.path(output_dir, paste0(sample_id, "_temp.bam"))
  sorted_bam <- file.path(output_dir, paste0(sample_id, "_aligned_sorted.bam"))
  
  # Convert SAM to BAM
  sam_to_bam_cmd <- sprintf("samtools view -@ %d -bS %s > %s", 
                            threads, shQuote(sam_file), shQuote(temp_bam))
  
  result1 <- system(sam_to_bam_cmd)
  if (result1 != 0) {
    return(list(success = FALSE, message = "Failed to convert SAM to BAM"))
  }
  
  # Sort BAM file
  sort_cmd <- sprintf("samtools sort -@ %d -o %s %s", 
                      threads, shQuote(sorted_bam), shQuote(temp_bam))
  
  result2 <- system(sort_cmd)
  if (result2 != 0) {
    return(list(success = FALSE, message = "Failed to sort BAM file"))
  }
  
  # Index sorted BAM
  index_cmd <- paste("samtools index", shQuote(sorted_bam))
  result3 <- system(index_cmd)
  
  # Clean up temporary files
  if (file.exists(temp_bam)) file.remove(temp_bam)
  if (!keep_sam && file.exists(sam_file)) file.remove(sam_file)
  
  if (file.exists(sorted_bam)) {
    file_size <- file.info(sorted_bam)$size / 1024^2
    message <- sprintf("Final BAM created successfully (%.1f MB)", file_size)
    return(list(success = TRUE, message = message, output_file = sorted_bam))
  }
  
  return(list(success = FALSE, message = "Failed to create sorted BAM"))
}

get_alignment_stats <- function(bam_file) {
  if (!file.exists(bam_file)) {
    return(list(success = FALSE, message = paste("BAM file not found:", bam_file)))
  }
  
  # Run samtools flagstat and capture output
  flagstat_cmd <- paste("samtools flagstat", shQuote(bam_file))
  stats_output <- system(flagstat_cmd, intern = TRUE)
  
  return(list(success = TRUE, stats = paste(stats_output, collapse = "\n")))
}

# Complete alignment pipeline
align_sample <- function(sample_id, fastq1, fastq2, reference_path, 
                         output_dir, threads = 4, platform = "ILLUMINA") {
  
  # Step 1: Check tools
  tool_check <- check_alignment_tools()
  if (!tool_check$success) {
    return(list(success = FALSE, message = paste("Missing tools:", paste(tool_check$missing, collapse = ", "))))
  }
  
  # Step 2: Index reference if needed
  index_result <- index_reference_genome(reference_path)
  if (!index_result$success) {
    return(index_result)
  }
  
  # Step 3: Run BWA-MEM
  bwa_result <- run_bwa_mem(sample_id, fastq1, fastq2, reference_path, 
                            output_dir, threads, platform)
  
  if (!bwa_result$success) {
    return(bwa_result)
  }
  
  # Step 4: Convert to sorted BAM
  bam_result <- sam_to_sorted_bam(bwa_result$output_file, output_dir, sample_id, 
                                  threads, FALSE)
  
  if (!bam_result$success) {
    return(bam_result)
  }
  
  # Step 5: Get statistics
  stats_result <- get_alignment_stats(bam_result$output_file)
  
  return(list(
    success = TRUE, 
    message = "Alignment completed successfully",
    output_file = bam_result$output_file,
    stats = if(stats_result$success) stats_result$stats else "Stats unavailable"
  ))
}

# =============================================================================
# SHINY GUI APPLICATION
# =============================================================================

# Define UI
ui <- fluidPage(
  titlePanel(paste("BWA-MEM Alignment Pipeline", 
                   if(use_shinyfiles) "(Enhanced File Browser)" else "(Basic Mode)")),
  
  tags$head(
    tags$style(HTML("
            .well { background-color: #f8f9fa; }
            .btn-primary { background-color: #007bff; border-color: #007bff; }
            .btn-success { background-color: #28a745; border-color: #28a745; }
            .btn-info { background-color: #17a2b8; border-color: #17a2b8; }
            .result-box { background-color: #f8f9fa; padding: 15px; border: 1px solid #dee2e6; margin: 10px 0; }
            pre { background-color: #f4f4f4; padding: 10px; border: 1px solid #ddd; }
            .file-path { background-color: #e9ecef; padding: 5px; margin: 5px 0; border-radius: 3px; font-family: monospace; }
            .mode-info { background-color: #d4edda; padding: 10px; border: 1px solid #c3e6cb; border-radius: 5px; margin-bottom: 15px; }
        "))
  ),
  
  # Mode information
  div(class = "mode-info",
      if(use_shinyfiles) {
        p("✅ Enhanced Mode: File browser enabled for easy file selection")
      } else {
        p("ℹ️  Basic Mode: Manual file path entry. Install 'shinyFiles' package for file browser functionality.")
      }
  ),
  
  sidebarLayout(
    sidebarPanel(width = 4,
                 h3("Configuration"),
                 
                 wellPanel(
                   h4("Sample Information"),
                   textInput("sample_id", "Sample ID:", value = "SRR26456208"),
                   
                   # Adaptive file input
                   h5("Forward Reads (R1):"),
                   if (use_shinyfiles) {
                     tagList(
                       shinyFilesButton("fastq1_browse", "Browse for R1 File", 
                                        "Select forward reads file", multiple = FALSE,
                                        style = "width: 100%; margin-bottom: 5px;"),
                       div(id = "fastq1_path", class = "file-path", "No file selected")
                     )
                   } else {
                     textInput("fastq1_manual", "Forward Reads Path:", 
                               value = "", 
                               placeholder = "/path/to/your/sample_R1.fastq.gz")
                   },
                   
                   h5("Reverse Reads (R2):"),
                   if (use_shinyfiles) {
                     tagList(
                       shinyFilesButton("fastq2_browse", "Browse for R2 File", 
                                        "Select reverse reads file", multiple = FALSE,
                                        style = "width: 100%; margin-bottom: 5px;"),
                       div(id = "fastq2_path", class = "file-path", "No file selected")
                     )
                   } else {
                     textInput("fastq2_manual", "Reverse Reads Path:", 
                               value = "", 
                               placeholder = "/path/to/your/sample_R2.fastq.gz")
                   }
                 ),
                 
                 wellPanel(
                   h4("Reference & Output"),
                   
                   h5("Reference Genome:"),
                   if (use_shinyfiles) {
                     tagList(
                       shinyFilesButton("reference_browse", "Browse for Reference", 
                                        "Select reference genome file", multiple = FALSE,
                                        style = "width: 100%; margin-bottom: 5px;"),
                       div(id = "reference_path", class = "file-path", "No file selected")
                     )
                   } else {
                     textInput("reference_manual", "Reference Genome Path:", 
                               value = "/home/dhibi/TMB_Pipeline_Project/reference/hg38.fa",
                               placeholder = "/path/to/reference.fa")
                   },
                   
                   h5("Output Directory:"),
                   if (use_shinyfiles) {
                     tagList(
                       shinyDirButton("output_browse", "Browse for Output Directory", 
                                      "Select output directory",
                                      style = "width: 100%; margin-bottom: 5px;"),
                       div(id = "output_path", class = "file-path", "/home/dhibi/TMB_Pipeline_Project/data/aligned")
                     )
                   } else {
                     textInput("output_manual", "Output Directory:", 
                               value = "/home/dhibi/TMB_Pipeline_Project/data/aligned")
                   },
                   
                   numericInput("threads", "Number of Threads:", value = 4, min = 1, max = 32),
                   selectInput("platform", "Sequencing Platform:", 
                               choices = c("ILLUMINA", "PACBIO", "IONTORRENT"), selected = "ILLUMINA")
                 ),
                 
                 wellPanel(
                   h4("Actions"),
                   actionButton("check_system", "Check System", class = "btn-info", style = "width: 100%; margin-bottom: 10px;"),
                   actionButton("check_files", "Check Input Files", class = "btn-info", style = "width: 100%; margin-bottom: 10px;"),
                   actionButton("run_alignment", "Run Alignment", class = "btn-success", style = "width: 100%; margin-bottom: 10px;")
                 )
    ),
    
    mainPanel(width = 8,
              tabsetPanel(
                tabPanel("System Check",
                         br(),
                         div(class = "result-box",
                             h4("System Status"),
                             verbatimTextOutput("system_check_result")
                         ),
                         
                         div(class = "result-box",
                             h4("Package Status"),
                             verbatimTextOutput("package_status")
                         ),
                         
                         div(class = "result-box",
                             h4("Installation Instructions"),
                             p("Required Tools: BWA (Burrows-Wheeler Aligner) and SAMtools"),
                             tags$pre("# Ubuntu/Debian:\nsudo apt update\nsudo apt install bwa samtools\n\n# CentOS/RHEL:\nsudo yum install bwa samtools\n\n# macOS (with Homebrew):\nbrew install bwa samtools\n\n# For enhanced file browser:\ninstall.packages('shinyFiles')")
                         )
                ),
                
                tabPanel("File Check",
                         br(),
                         div(class = "result-box",
                             h4("Selected Files"),
                             verbatimTextOutput("selected_files")
                         ),
                         
                         div(class = "result-box",
                             h4("Input File Validation"),
                             verbatimTextOutput("file_check_result")
                         )
                ),
                
                tabPanel("Alignment Results",
                         br(),
                         div(class = "result-box",
                             h4("Alignment Status"),
                             verbatimTextOutput("alignment_result")
                         ),
                         
                         div(class = "result-box",
                             h4("Alignment Statistics"),
                             verbatimTextOutput("alignment_stats")
                         )
                ),
                
                tabPanel("Batch Processing",
                         br(),
                         div(class = "result-box",
                             h4("Batch Configuration"),
                             p("For batch processing, create a CSV file with columns: sample_id, fastq1, fastq2"),
                             fileInput("sample_sheet", "Upload Sample Sheet (CSV):", accept = c(".csv")),
                             
                             if (use_shinyfiles) {
                               tagList(
                                 h5("Batch Reference Genome:"),
                                 shinyFilesButton("batch_reference_browse", "Browse for Reference", 
                                                  "Select reference genome for batch", multiple = FALSE),
                                 div(id = "batch_reference_path", class = "file-path", "No file selected"),
                                 
                                 h5("Batch Output Directory:"),
                                 shinyDirButton("batch_output_browse", "Browse for Output Directory", 
                                                "Select output directory for batch"),
                                 div(id = "batch_output_path", class = "file-path", "No directory selected")
                               )
                             } else {
                               tagList(
                                 textInput("batch_reference", "Reference Genome:", 
                                           value = "/home/dhibi/TMB_Pipeline_Project/reference/hg38.fa"),
                                 textInput("batch_output", "Output Directory:", 
                                           value = "/home/dhibi/TMB_Pipeline_Project/data/aligned")
                               )
                             },
                             
                             actionButton("run_batch", "Run Batch Alignment", class = "btn-success")
                         ),
                         
                         div(class = "result-box",
                             h4("Sample Preview"),
                             tableOutput("sample_preview")
                         ),
                         
                         div(class = "result-box",
                             h4("Batch Results"),
                             verbatimTextOutput("batch_results")
                         )
                )
              )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Initialize reactive values
  values <- reactiveValues(
    fastq1_path = "",
    fastq2_path = "",
    reference_path = "",
    output_dir = "/home/dhibi/TMB_Pipeline_Project/data/aligned",
    batch_reference_path = "",
    batch_output_dir = ""
  )
  
  # Set up file system for shinyFiles (if available)
  if (use_shinyfiles) {
    roots <- c(home = "~", root = "/")
    
    # File browsers
    shinyFileChoose(input, "fastq1_browse", roots = roots, session = session,
                    filetypes = c("fastq", "fq", "gz"))
    shinyFileChoose(input, "fastq2_browse", roots = roots, session = session,
                    filetypes = c("fastq", "fq", "gz"))
    shinyFileChoose(input, "reference_browse", roots = roots, session = session,
                    filetypes = c("fa", "fasta", "fna"))
    shinyFileChoose(input, "batch_reference_browse", roots = roots, session = session,
                    filetypes = c("fa", "fasta", "fna"))
    
    # Directory browsers
    shinyDirChoose(input, "output_browse", roots = roots, session = session)
    shinyDirChoose(input, "batch_output_browse", roots = roots, session = session)
    
    # Update file paths when files are selected
    observeEvent(input$fastq1_browse, {
      if (!is.null(input$fastq1_browse) && !is.integer(input$fastq1_browse)) {
        values$fastq1_path <- as.character(parseFilePaths(roots, input$fastq1_browse)$datapath)
        runjs(sprintf("document.getElementById('fastq1_path').innerHTML = '%s';", values$fastq1_path))
      }
    })
    
    observeEvent(input$fastq2_browse, {
      if (!is.null(input$fastq2_browse) && !is.integer(input$fastq2_browse)) {
        values$fastq2_path <- as.character(parseFilePaths(roots, input$fastq2_browse)$datapath)
        runjs(sprintf("document.getElementById('fastq2_path').innerHTML = '%s';", values$fastq2_path))
      }
    })
    
    observeEvent(input$reference_browse, {
      if (!is.null(input$reference_browse) && !is.integer(input$reference_browse)) {
        values$reference_path <- as.character(parseFilePaths(roots, input$reference_browse)$datapath)
        runjs(sprintf("document.getElementById('reference_path').innerHTML = '%s';", values$reference_path))
      }
    })
    
    observeEvent(input$output_browse, {
      if (!is.null(input$output_browse) && !is.integer(input$output_browse)) {
        values$output_dir <- as.character(parseDirPath(roots, input$output_browse))
        runjs(sprintf("document.getElementById('output_path').innerHTML = '%s';", values$output_dir))
      }
    })
    
    observeEvent(input$batch_reference_browse, {
      if (!is.null(input$batch_reference_browse) && !is.integer(input$batch_reference_browse)) {
        values$batch_reference_path <- as.character(parseFilePaths(roots, input$batch_reference_browse)$datapath)
        runjs(sprintf("document.getElementById('batch_reference_path').innerHTML = '%s';", values$batch_reference_path))
      }
    })
    
    observeEvent(input$batch_output_browse, {
      if (!is.null(input$batch_output_browse) && !is.integer(input$batch_output_browse)) {
        values$batch_output_dir <- as.character(parseDirPath(roots, input$batch_output_browse))
        runjs(sprintf("document.getElementById('batch_output_path').innerHTML = '%s';", values$batch_output_dir))
      }
    })
  } else {
    # Handle manual input when shinyFiles is not available
    observeEvent(input$fastq1_manual, {
      values$fastq1_path <- input$fastq1_manual
    })
    
    observeEvent(input$fastq2_manual, {
      values$fastq2_path <- input$fastq2_manual
    })
    
    observeEvent(input$reference_manual, {
      values$reference_path <- input$reference_manual
    })
    
    observeEvent(input$output_manual, {
      values$output_dir <- input$output_manual
    })
    
    observeEvent(input$batch_reference, {
      values$batch_reference_path <- input$batch_reference
    })
    
    observeEvent(input$batch_output, {
      values$batch_output_dir <- input$batch_output
    })
  }
  
  # Display selected files
  output$selected_files <- renderText({
    paste(
      "Forward reads (R1):", ifelse(values$fastq1_path != "", values$fastq1_path, "Not selected"),
      "Reverse reads (R2):", ifelse(values$fastq2_path != "", values$fastq2_path, "Not selected"),
      "Reference genome:", ifelse(values$reference_path != "", values$reference_path, "Not selected"),
      "Output directory:", values$output_dir,
      sep = "\n"
    )
  })
  
  # Package status output
  output$package_status <- renderText({
    paste(
      "=== PACKAGE STATUS ===",
      paste("✅ shiny:", packageVersion("shiny")),
      if(use_shinyfiles) paste("✅ shinyFiles:", packageVersion("shinyFiles")) else "❌ shinyFiles: Not available",
      "",
      "Mode:", if(use_shinyfiles) "Enhanced (File Browser)" else "Basic (Manual Input)",
      sep = "\n"
    )
  })
  
  # System check
  observeEvent(input$check_system, {
    tool_check <- check_alignment_tools()
    
    output$system_check_result <- renderText({
      if (tool_check$success) {
        "✅ All required tools are available:\n• BWA: Available\n• SAMtools: Available\n\nSystem is ready for alignment!"
      } else {
        paste("❌ Missing tools:", paste(tool_check$missing, collapse = ", "), 
              "\n\nPlease install missing tools using the commands in the Installation Instructions section.")
      }
    })
  })
  
  # Check input files
  observeEvent(input$check_files, {
    files_to_check <- c(
      "Forward reads" = values$fastq1_path,
      "Reverse reads" = values$fastq2_path,
      "Reference genome" = values$reference_path
    )
    
    output$file_check_result <- renderText({
      results <- sapply(files_to_check, function(x) x != "" && file.exists(x))
      result_text <- c()
      
      for (i in seq_along(files_to_check)) {
        if (files_to_check[i] == "") {
          status <- "❌ Not selected"
        } else if (results[i]) {
          size_mb <- round(file.info(files_to_check[i])$size / 1024^2, 1)
          status <- sprintf("✅ Found (%.1f MB)", size_mb)
        } else {
          status <- "❌ File not found"
        }
        result_text <- c(result_text, paste(names(files_to_check)[i], ":", status))
      }
      
      all_found <- all(results) && all(files_to_check != "")
      summary <- ifelse(all_found, 
                        "\n✅ All files found - ready to proceed!",
                        "\n❌ Some files are missing or not selected")
      
      paste(c(result_text, summary), collapse = "\n")
    })
  })
  
  # Run single alignment
  observeEvent(input$run_alignment, {
    # Validate inputs
    if (values$fastq1_path == "" || values$fastq2_path == "" || values$reference_path == "") {
      output$alignment_result <- renderText("❌ Please select all required files first")
      return()
    }
    
    # Clear previous results
    output$alignment_result <- renderText("🔄 Running alignment... This may take several minutes.")
    output$alignment_stats <- renderText("")
    
    # Run alignment
    result <- align_sample(
      sample_id = input$sample_id,
      fastq1 = values$fastq1_path,
      fastq2 = values$fastq2_path,
      reference_path = values$reference_path,
      output_dir = values$output_dir,
      threads = input$threads,
      platform = input$platform
    )
    
    # Update results
    if (result$success) {
      output$alignment_result <- renderText({
        paste("✅ ALIGNMENT SUCCESSFUL!\n\n", 
              "Output file:", result$output_file, "\n",
              "Status:", result$message, "\n\n",
              "You can now proceed to variant calling or other downstream analysis.")
      })
      
      output$alignment_stats <- renderText({
        result$stats
      })
    } else {
      output$alignment_result <- renderText({
        paste("❌ ALIGNMENT FAILED\n\n", 
              "Error:", result$message, "\n\n",
              "Please check:\n",
              
              "• File paths are correct\n",
              "• BWA and samtools are installed\n",
              "• Sufficient disk space available\n",
              "• Input files are not corrupted")
      })
    }
  })
  
  # Sample sheet preview
  output$sample_preview <- renderTable({
    if (!is.null(input$sample_sheet)) {
      df <- read.csv(input$sample_sheet$datapath)
      head(df, 10)  # Show first 10 rows
    }
  })
  
  # Batch processing
  observeEvent(input$run_batch, {
    if (is.null(input$sample_sheet)) {
      output$batch_results <- renderText("❌ Please upload a sample sheet first")
      return()
    }
    
    batch_ref <- if (use_shinyfiles) values$batch_reference_path else input$batch_reference
    batch_out <- if (use_shinyfiles) values$batch_output_dir else input$batch_output
    
    if (batch_ref == "" || batch_out == "") {
      output$batch_results <- renderText("❌ Please select reference genome and output directory")
      return()
    }
    
    output$batch_results <- renderText("🔄 Running batch alignment... This will take a while.")
    
    # Read sample sheet
    samples_df <- read.csv(input$sample_sheet$datapath)
    
    if (!all(c("sample_id", "fastq1", "fastq2") %in% colnames(samples_df))) {
      output$batch_results <- renderText("❌ CSV must have columns: sample_id, fastq1, fastq2")
      return()
    }
    
    # Process each sample
    results <- list()
    progress_text <- c()
    
    for (i in 1:nrow(samples_df)) {
      sample_row <- samples_df[i, ]
      progress_text <- c(progress_text, paste("Processing", sample_row$sample_id, "..."))
      
      result <- align_sample(
        sample_id = sample_row$sample_id,
        fastq1 = sample_row$fastq1,
        fastq2 = sample_row$fastq2,
        reference_path = batch_ref,
        output_dir = batch_out,
        threads = input$threads
      )
      
      results[[sample_row$sample_id]] <- result$success
      status <- ifelse(result$success, "✅ Success", "❌ Failed")
      progress_text <- c(progress_text, paste(sample_row$sample_id, ":", status))
    }
    
    # Summary
    successful <- sum(unlist(results))
    total <- length(results)
    summary_text <- paste("\n=== BATCH SUMMARY ===",
                          paste("Completed:", successful, "/", total, "samples"),
                          sep = "\n")
    
    output$batch_results <- renderText({
      paste(c(progress_text, summary_text), collapse = "\n")
    })
  })
}

# =============================================================================
# LAUNCH THE APPLICATION
# =============================================================================

# Print startup message
cat("=== BWA ALIGNMENT GUI WITH FILE BROWSER STARTING ===\n")
cat("Loading Shiny interface...\n")

if (use_shinyfiles) {
  cat("✅ File browser functionality enabled (shinyFiles package loaded)\n")
} else {
  cat("⚠️  Basic file upload mode (shinyFiles package not available)\n")
  cat("   You can upload files directly or install shinyFiles for better browsing:\n")
  cat("   install.packages('shinyFiles')\n")
}

cat("\nTo run this script again in the future:\n")
cat("1. Save this file as 'bwa_alignment_gui.R'\n")
cat("2. In R/RStudio: source('bwa_alignment_gui.R')\n")
cat("3. Or from terminal: Rscript bwa_alignment_gui.R\n\n")
launch_app <- function() {
  shinyApp(ui = ui, server = server, 
           options = list(launch.browser = TRUE, 
                          host = "0.0.0.0", 
                          port = 8080))
}

# Then call it

launch_app()
# Launch the Shiny app
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))

cat("=== GUI LAUNCHED ===\n")
cat("The interface should open in your web browser.\n")
cat("If it doesn't open automatically, check the R console for the local URL.\n")
cat("\nFeatures:\n")
cat("• File browser for selecting FASTQ files and reference genome\n")
cat("• Directory browser for output location\n")
cat("• Real-time file validation\n")
cat("• Complete BWA-MEM alignment pipeline\n")
cat("• Batch processing support\n")
# Add this function to your file
