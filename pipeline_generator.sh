#!/bin/bash
# =============================================================================
#
# RNA-seq Pipeline Generator Script (Final Fix)
#
# This script creates a complete, flexible RNA-seq analysis pipeline structure
# for submission to an HTCondor cluster.
#
# Key Features:
# - Correct R Wrappers: Shell wrappers now correctly source the config file.
# - Robust R Scripts: R scripts now fail correctly if inputs are missing.
# - Output Validation: Explicitly checks if STAR/Salmon created their output.
# - Simplified Config: Uses a single '3_samples.tsv' for sample paths and design.
# - Conda-Aware: Automatically adds 'conda activate' to each script.
# - Interactive: Asks the user before overwriting existing config files.
#
# =============================================================================

set -e

echo "🚀 Starting RNA-seq pipeline generation..."

# --- Create Core Directories ---
echo "-> Creating directory structure..."
mkdir -p scripts
mkdir -p condor_submit_files
mkdir -p results/{fastqc,trimmed_reads,star_alignment,salmon_quant,featurecounts,deseq2_star,deseq2_salmon,fgsea_star,fgsea_salmon,multiqc}
mkdir -p logs/{fastqc,trim,star,salmon,featurecounts,merge,deseq2,fgsea,multiqc}

# =============================================================================
#
# GENERATE TEMPLATE CONFIGURATION FILES (with overwrite prompt)
#
# =============================================================================

# --- Generate 2_config.sh ---
if [ -f "2_config.sh" ]; then
    read -p "-> File '2_config.sh' already exists. Overwrite? (y/n) " -n 1 -r
    echo # Move to new line
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "   -> Overwriting '2_config.sh'..."
        cat > 2_config.sh << 'EOF'
#!/bin/bash
# =============================================================================
# Configuration file for the RNA-seq pipeline
# EDIT ALL PATHS AND VARIABLES TO MATCH YOUR SYSTEM AND PROJECT
# =============================================================================

# --- Aligner Selection ---
# Choose between "STAR" or "SALMON". This is the master switch for the pipeline.
export ALIGNER="STAR"

# --- Project Directories ---
# Use absolute paths to avoid issues with Condor's working directories.
export PROJECT_DIR="$(pwd)" # This assumes you run from the project root. Best to set an absolute path.
export SAMPLE_SHEET="${PROJECT_DIR}/3_samples.tsv"

# --- Reference and Annotation Files ---
export GENOME_FASTA="/path/to/your/genome.fa"
export GTF_FILE="/path/to/your/annotation.gtf"
export STAR_INDEX_DIR="/path/to/your/star_index"
export SALMON_INDEX_DIR="/path/to/your/salmon_index"
export GSEA_GENESETS_GMT="/path/to/your/c2.cp.kegg.v7.5.1.symbols.gmt" # Example GMT file

# --- Tool-specific Parameters ---
# Find this path with: find $(conda info --base)/envs/qc_env -name "TruSeq3-PE-2.fa"
export TRIMMOMATIC_ADAPTERS="/path/to/your/trimmomatic/adapters/TruSeq3-PE-2.fa"
# For featureCounts: 0=unstranded, 1=stranded, 2=reversely stranded.
# Modern Illumina kits (like NEBNext Ultra II Directional) are reversely stranded.
export STRANDEDNESS="2"
EOF
    else
        echo "   -> Skipping '2_config.sh'."
    fi
else
    echo "-> Generating template '2_config.sh'..."
    cat > 2_config.sh << 'EOF'
#!/bin/bash
# =============================================================================
# Configuration file for the RNA-seq pipeline
# EDIT ALL PATHS AND VARIABLES TO MATCH YOUR SYSTEM AND PROJECT
# =============================================================================

# --- Aligner Selection ---
# Choose between "STAR" or "SALMON". This is the master switch for the pipeline.
export ALIGNER="STAR"

# --- Project Directories ---
# Use absolute paths to avoid issues with Condor's working directories.
export PROJECT_DIR="$(pwd)" # This assumes you run from the project root. Best to set an absolute path.
export RESULTS_DIR="${PROJECT_DIR}/results"
export SAMPLE_SHEET="${PROJECT_DIR}/3_samples.tsv"

# --- Reference and Annotation Files ---
export GENOME_FASTA="/path/to/your/genome.fa"
export GTF_FILE="/path/to/your/annotation.gtf"
export STAR_INDEX_DIR="/path/to/your/star_index"
export SALMON_INDEX_DIR="/path/to/your/salmon_index"
export GSEA_GENESETS_GMT="/path/to/your/c2.cp.kegg.v7.5.1.symbols.gmt" # Example GMT file

# --- Tool-specific Parameters ---
# Find this path with: find $(conda info --base)/envs/qc_env -name "TruSeq3-PE-2.fa"
export TRIMMOMATIC_ADAPTERS="/path/to/your/trimmomatic/adapters/TruSeq3-PE-2.fa"
# For featureCounts: 0=unstranded, 1=stranded, 2=reversely stranded.
# Modern Illumina kits (like NEBNext Ultra II Directional) are reversely stranded.
export STRANDEDNESS="2"
EOF
fi

# --- Generate 3_samples.tsv ---
if [ -f "3_samples.tsv" ]; then
    read -p "-> File '3_samples.tsv' already exists. Overwrite? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "   -> Overwriting '3_samples.tsv'..."
        cat > 3_samples.tsv << 'EOF'
# Define the DESeq2 model and contrast at the top of this file.
# Use column names from the header row below.
#DESEQ_DESIGN=~condition
#DESEQ_CONTRAST=condition,treatment,control
sample_id	fastq_1	fastq_2	condition
sample_A1	/path/to/fastqs/A1_R1.fastq.gz	/path/to/fastqs/A1_R2.fastq.gz	control
sample_A2	/path/to/fastqs/A2_R1.fastq.gz	/path/to/fastqs/A2_R2.fastq.gz	control
sample_B1	/path/to/fastqs/B1_R1.fastq.gz	/path/to/fastqs/B1_R2.fastq.gz	treatment
sample_B2	/path/to/fastqs/B2_R1.fastq.gz	/path/to/fastqs/B2_R2.fastq.gz	treatment
EOF
    else
        echo "   -> Skipping '3_samples.tsv'."
    fi
else
    echo "-> Generating template '3_samples.tsv'..."
    cat > 3_samples.tsv << 'EOF'
# Define the DESeq2 model and contrast at the top of this file.
# Use column names from the header row below.
#DESEQ_DESIGN=~condition
#DESEQ_CONTRAST=condition,treatment,control
sample_id	fastq_1	fastq_2	condition
sample_A1	/path/to/fastqs/A1_R1.fastq.gz	/path/to/fastqs/A1_R2.fastq.gz	control
sample_A2	/path/to/fastqs/A2_R1.fastq.gz	/path/to/fastqs/A2_R2.fastq.gz	control
sample_B1	/path/to/fastqs/B1_R1.fastq.gz	/path/to/fastqs/B1_R2.fastq.gz	treatment
sample_B2	/path/to/fastqs/B2_R1.fastq.gz	/path/to/fastqs/B2_R2.fastq.gz	treatment
EOF
fi


# =============================================================================
#
# GENERATE EXECUTABLE SCRIPTS (scripts/) - These are always overwritten
#
# =============================================================================
echo "-> Generating/updating core pipeline scripts (with Conda activation)..."

# --- Conda activation snippet ---
CONDA_ACTIVATION='
# --- Conda Environment Activation ---
source $(conda info --base)/etc/profile.d/conda.sh
if [ $? -ne 0 ]; then
    echo "Error: Failed to source conda.sh. Make sure Conda is installed correctly." >&2
    exit 1
fi
conda activate %s
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment: %s" >&2
    exit 1
fi
# --- End Activation ---
'

# --- 01_fastqc.sh ---
cat > scripts/01_fastqc.sh << EOF
#!/bin/bash
set -e
SAMPLE_ID=\$1
echo "--- FASTQC SCRIPT ---"
echo "--- Received SAMPLE_ID: '\${SAMPLE_ID}' ---"
$(printf "$CONDA_ACTIVATION" "qc_env" "qc_env")
source 2_config.sh
FASTQ1=\$(awk -F'\\t' -v id="\$SAMPLE_ID" '\$1 == id {print \$2}' \${SAMPLE_SHEET})
FASTQ2=\$(awk -F'\\t' -v id="\$SAMPLE_ID" '\$1 == id {print \$3}' \${SAMPLE_SHEET})
echo "Running FastQC on \${SAMPLE_ID}"
fastqc --threads 2 -o results/fastqc \${FASTQ1} \${FASTQ2}
EOF

# --- 02_trim.sh ---
cat > scripts/02_trim.sh << EOF
#!/bin/bash
set -e
SAMPLE_ID=\$1
echo "--- TRIMMOMATIC SCRIPT ---"
echo "--- Received SAMPLE_ID: '\${SAMPLE_ID}' ---"
$(printf "$CONDA_ACTIVATION" "qc_env" "qc_env")
source 2_config.sh
FASTQ1=\$(awk -F'\\t' -v id="\$SAMPLE_ID" '\$1 == id {print \$2}' \${SAMPLE_SHEET})
FASTQ2=\$(awk -F'\\t' -v id="\$SAMPLE_ID" '\$1 == id {print \$3}' \${SAMPLE_SHEET})
OUT_PAIRED_1=results/trimmed_reads/\${SAMPLE_ID}_1_paired.fq.gz
OUT_UNPAIRED_1=results/trimmed_reads/\${SAMPLE_ID}_1_unpaired.fq.gz
OUT_PAIRED_2=results/trimmed_reads/\${SAMPLE_ID}_2_paired.fq.gz
OUT_UNPAIRED_2=results/trimmed_reads/\${SAMPLE_ID}_2_unpaired.fq.gz
echo "Running Trimmomatic on \${SAMPLE_ID}"
trimmomatic PE -threads 4 \\
    \${FASTQ1} \${FASTQ2} \\
    \${OUT_PAIRED_1} \${OUT_UNPAIRED_1} \\
    \${OUT_PAIRED_2} \${OUT_UNPAIRED_2} \\
    ILLUMINACLIP:\${TRIMMOMATIC_ADAPTERS}:2:30:10 \\
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
EOF

# --- 03a_star_align.sh ---
cat > scripts/03a_star_align.sh << EOF
#!/bin/bash
set -e
SAMPLE_ID=\$1
echo "--- STAR ALIGN SCRIPT ---"
echo "--- Received SAMPLE_ID: '\${SAMPLE_ID}' ---"
$(printf "$CONDA_ACTIVATION" "align_env" "align_env")
source 2_config.sh
READ1=results/trimmed_reads/\${SAMPLE_ID}_1_paired.fq.gz
READ2=results/trimmed_reads/\${SAMPLE_ID}_2_paired.fq.gz
OUT_PREFIX=results/star_alignment/\${SAMPLE_ID}_
BAM_FILE=\${OUT_PREFIX}Aligned.sortedByCoord.out.bam

echo "Running STAR alignment on \${SAMPLE_ID}"
STAR --runThreadN 8 \\
     --genomeDir \${STAR_INDEX_DIR} \\
     --readFilesIn \${READ1} \${READ2} \\
     --readFilesCommand zcat \\
     --outFileNamePrefix \${OUT_PREFIX} \\
     --outSAMtype BAM SortedByCoordinate \\
     --outSAMattributes Standard

# --- CRITICAL VALIDATION STEP ---
echo "Validating STAR output..."
if [ ! -f "\${BAM_FILE}" ]; then
    echo "ERROR: STAR did not produce the output BAM file: \${BAM_FILE}" >&2
    echo "Check the STAR log file for errors: \${OUT_PREFIX}Log.final.out" >&2
    exit 1
fi
echo "STAR output BAM file found. Alignment successful."
EOF

# --- 03b_salmon_quant.sh ---
cat > scripts/03b_salmon_quant.sh << EOF
#!/bin/bash
set -e
SAMPLE_ID=\$1
echo "--- SALMON QUANT SCRIPT ---"
echo "--- Received SAMPLE_ID: '\${SAMPLE_ID}' ---"
$(printf "$CONDA_ACTIVATION" "align_env" "align_env")
source 2_config.sh
READ1=results/trimmed_reads/\${SAMPLE_ID}_1_paired.fq.gz
READ2=results/trimmed_reads/\${SAMPLE_ID}_2_paired.fq.gz
OUT_DIR=results/salmon_quant/\${SAMPLE_ID}
JSON_INFO=\${OUT_DIR}/aux_info/meta_info.json

echo "Running Salmon quantification on \${SAMPLE_ID}"
salmon quant -i \${SALMON_INDEX_DIR} \\
             -l A \\
             -1 \${READ1} \\
             -2 \${READ2} \\
             --validateMappings \\
             --gcBias \\
             --seqBias \\
             -p 8 \\
             -o \${OUT_DIR}

# --- CRITICAL VALIDATION STEP ---
echo "Validating Salmon output..."
if [ ! -f "\${JSON_INFO}" ]; then
    echo "ERROR: Salmon did not produce the output meta_info.json file in \${OUT_DIR}" >&2
    exit 1
fi
echo "Salmon output found. Quantification successful."
EOF

# --- 04a_featurecounts.sh ---
cat > scripts/04a_featurecounts.sh << EOF
#!/bin/bash
set -e
SAMPLE_ID=\$1
echo "--- FEATURECOUNTS SCRIPT ---"
echo "--- Received SAMPLE_ID: '\${SAMPLE_ID}' ---"
$(printf "$CONDA_ACTIVATION" "align_env" "align_env")
source 2_config.sh
BAM_FILE=results/star_alignment/\${SAMPLE_ID}_Aligned.sortedByCoord.out.bam
OUT_FILE=results/featurecounts/\${SAMPLE_ID}_counts.txt

echo "Running featureCounts on: \${SAMPLE_ID}"
echo "Input BAM: \${BAM_FILE}"
echo "Annotation GTF: \${GTF_FILE}"
echo "Strandedness: \${STRANDEDNESS}"

featureCounts -p -s \${STRANDEDNESS} \\
              -a \${GTF_FILE} \\
              -o \${OUT_FILE} \\
              -T 4 \\
              -g gene_id \\
              \${BAM_FILE}
EOF

# --- R SCRIPT WRAPPERS AND PURE R SCRIPTS ---

# --- 05a_merge_counts_star.sh (Wrapper) ---
cat > scripts/05a_merge_counts_star.sh << EOF
#!/bin/bash
set -e
echo "--- MERGE COUNTS WRAPPER SCRIPT ---"
$(printf "$CONDA_ACTIVATION" "r_env" "r_env")
source 2_config.sh
echo "Executing R script..."
Rscript scripts/05a_merge_counts_star.R
EOF

# --- 05a_merge_counts_star.R (Pure R) ---
cat > scripts/05a_merge_counts_star.R << 'EOF'
suppressPackageStartupMessages(library(data.table))

# Define paths relative to the project root
featurecounts_dir <- "results/featurecounts"
output_file <- file.path(featurecounts_dir, "merged_gene_counts.csv")
files_to_merge <- list.files(featurecounts_dir, pattern = "*_counts.txt$", full.names = TRUE)

# --- ROBUSTNESS CHECK ---
if (length(files_to_merge) == 0) {
    write("ERROR: No featureCounts files (*_counts.txt) found to merge.", stderr())
    q(status = 1)
}

dt_list <- lapply(files_to_merge, function(f) {
    dt <- fread(f, skip = 1, sep = "\t")
    dt <- dt[, c(1, 7)]
    sample_name <- basename(names(dt)[2])
    sample_name <- sub("_Aligned.sortedByCoord.out.bam", "", sample_name)
    setnames(dt, c("Geneid", sample_name))
    dt
})

merged_dt <- Reduce(function(...) merge(..., by = "Geneid", all = TRUE), dt_list)
fwrite(merged_dt, file = output_file)
cat("Successfully merged", length(files_to_merge), "files into", output_file, "\n")
EOF

# --- 06a_deseq2_star.sh (Wrapper) ---
cat > scripts/06a_deseq2_star.sh << EOF
#!/bin/bash
set -e
echo "--- DESEQ2 STAR WRAPPER SCRIPT ---"
$(printf "$CONDA_ACTIVATION" "r_env" "r_env")
source 2_config.sh
echo "Executing R script..."
Rscript scripts/06a_deseq2_star.R
EOF

# --- 06a_deseq2_star.R (Pure R) ---
cat > scripts/06a_deseq2_star.R << 'EOF'
suppressPackageStartupMessages({
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
})

# Read configuration from environment variables for specific files
SAMPLE_SHEET <- Sys.getenv("SAMPLE_SHEET")

# Define paths relative to project root
COUNT_MATRIX_PATH <- "results/featurecounts/merged_gene_counts.csv"
RESULTS_DIR_DESEQ <- "results/deseq2_star"

# --- ROBUSTNESS CHECK ---
if (!file.exists(COUNT_MATRIX_PATH)) {
    write(paste("ERROR: Input count matrix not found at:", COUNT_MATRIX_PATH), stderr())
    q(status = 1)
}

dir.create(RESULTS_DIR_DESEQ, showWarnings = FALSE, recursive = TRUE)

design_info <- readLines(SAMPLE_SHEET)
design_formula_line <- grep("^#DESEQ_DESIGN=", design_info, value = TRUE)
contrast_line <- grep("^#DESEQ_CONTRAST=", design_info, value = TRUE)
design_formula <- as.formula(gsub("#DESEQ_DESIGN=", "", design_formula_line))
contrast_params <- strsplit(gsub("#DESEQ_CONTRAST=", "", contrast_line), ",")[[1]]

count_data <- read_csv(COUNT_MATRIX_PATH, show_col_types = FALSE)
sample_info <- read_tsv(SAMPLE_SHEET, comment = "#", show_col_types = FALSE) %>%
    column_to_rownames("sample_id")
count_matrix <- count_data %>%
    column_to_rownames("Geneid") %>%
    select(all_of(rownames(sample_info)))

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = design_formula)
dds <- DESeq(dds)
res <- results(dds, contrast = contrast_params)
res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% arrange(padj)
write_csv(res_df, file.path(RESULTS_DIR_DESEQ, "deseq2_results.csv"))

vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup=contrast_params[1], returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=!!sym(contrast_params[1]))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + ggtitle("PCA Plot (STAR)")
ggsave(file.path(RESULTS_DIR_DESEQ, "PCA_plot_star.pdf"), pca_plot)
print("DESeq2 (STAR) analysis complete.")
EOF

# --- 06b_deseq2_salmon.sh (Wrapper) ---
cat > scripts/06b_deseq2_salmon.sh << EOF
#!/bin/bash
set -e
echo "--- DESEQ2 SALMON WRAPPER SCRIPT ---"
$(printf "$CONDA_ACTIVATION" "r_env" "r_env")
source 2_config.sh
echo "Executing R script..."
Rscript scripts/06b_deseq2_salmon.R
EOF

# --- 06b_deseq2_salmon.R (Pure R) ---
cat > scripts/06b_deseq2_salmon.R << 'EOF'
suppressPackageStartupMessages({
    library(DESeq2)
    library(tximport)
    library(readr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(rtracklayer)
})

# Read configuration from environment variables
SAMPLE_SHEET <- Sys.getenv("SAMPLE_SHEET")
GTF_FILE <- Sys.getenv("GTF_FILE")

# Define paths relative to project root
RESULTS_DIR_SALMON <- "results/salmon_quant"
RESULTS_DIR_DESEQ <- "results/deseq2_salmon"
dir.create(RESULTS_DIR_DESEQ, showWarnings = FALSE, recursive = TRUE)

design_info <- readLines(SAMPLE_SHEET)
design_formula_line <- grep("^#DESEQ_DESIGN=", design_info, value = TRUE)
contrast_line <- grep("^#DESEQ_CONTRAST=", design_info, value = TRUE)
design_formula <- as.formula(gsub("#DESEQ_DESIGN=", "", design_formula_line))
contrast_params <- strsplit(gsub("#DESEQ_CONTRAST=", "", contrast_line), ",")[[1]]

sample_info <- read_tsv(SAMPLE_SHEET, comment = "#", show_col_types = FALSE)
files <- file.path(RESULTS_DIR_SALMON, sample_info$sample_id, "quant.sf")
names(files) <- sample_info$sample_id

# --- ROBUSTNESS CHECK ---
if (any(!file.exists(files))) {
    write("ERROR: Not all Salmon quant.sf files were found.", stderr())
    q(status = 1)
}

gtf <- rtracklayer::import(GTF_FILE)
tx2gene <- as.data.frame(mcols(gtf)[, c("transcript_id", "gene_name")]) %>%
    filter(!is.na(transcript_id) & !is.na(gene_name)) %>% distinct()

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
sample_table <- as.data.frame(sample_info)
rownames(sample_table) <- sample_table$sample_id
dds <- DESeqDataSetFromTximport(txi, colData=sample_table, design=design_formula)
dds <- DESeq(dds)
res <- results(dds, contrast = contrast_params)
res_df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>% arrange(padj)
write_csv(res_df, file.path(RESULTS_DIR_DESEQ, "deseq2_results.csv"))

vsd <- vst(dds, blind=FALSE)
pca_data <- plotPCA(vsd, intgroup=contrast_params[1], returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color=!!sym(contrast_params[1]))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + ggtitle("PCA Plot (Salmon)")
ggsave(file.path(RESULTS_DIR_DESEQ, "PCA_plot_salmon.pdf"), pca_plot)
print("DESeq2 (Salmon) analysis complete.")
EOF

# --- 07a_fgsea_star.sh (Wrapper) ---
cat > scripts/07a_fgsea_star.sh << EOF
#!/bin/bash
set -e
echo "--- FGSEA STAR WRAPPER SCRIPT ---"
$(printf "$CONDA_ACTIVATION" "r_env" "r_env")
source 2_config.sh
echo "Executing R script..."
Rscript scripts/07a_fgsea_star.R
EOF

# --- 07a_fgsea_star.R (Pure R) ---
cat > scripts/07a_fgsea_star.R << 'EOF'
suppressPackageStartupMessages({
    library(fgsea)
    library(readr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
})

# Read configuration from environment variables
GSEA_GENESETS_GMT <- Sys.getenv("GSEA_GENESETS_GMT")

# Define paths relative to project root
DESEQ_RESULTS_PATH <- "results/deseq2_star/deseq2_results.csv"
RESULTS_DIR_FGSEA <- "results/fgsea_star"

# --- ROBUSTNESS CHECK ---
if (!file.exists(DESEQ_RESULTS_PATH)) {
    write(paste("ERROR: Input DESeq2 results not found at:", DESEQ_RESULTS_PATH), stderr())
    q(status = 1)
}

dir.create(RESULTS_DIR_FGSEA, showWarnings = FALSE, recursive = TRUE)

deseq_res <- read_csv(DESEQ_RESULTS_PATH, show_col_types = FALSE)
ranked_stats <- deseq_res %>%
    filter(!is.na(stat)) %>%
    distinct(gene_id, .keep_all=TRUE) %>%
    pull(stat, name=gene_id)

pathways <- gmtPathways(GSEA_GENESETS_GMT)
fgsea_res <- fgsea(pathways, ranked_stats, minSize=15, maxSize=500)
fgsea_res_tbl <- fgsea_res %>% as_tibble() %>% arrange(desc(NES))
write_tsv(fgsea_res_tbl, file.path(RESULTS_DIR_FGSEA, "fgsea_results.tsv"))

top_pathways_up <- fgsea_res[ES > 0][head(order(pval), n=10), pathway]
top_pathways_down <- fgsea_res[ES < 0][head(order(pval), n=10), pathway]
top_pathways <- c(top_pathways_up, rev(top_pathways_down))
pdf(file.path(RESULTS_DIR_FGSEA, "fgsea_top_pathways.pdf"))
plotGseaTable(pathways[top_pathways], ranked_stats, fgsea_res, gseaParam=0.5)
dev.off()
print("fGSEA (STAR) analysis complete.")
EOF

# --- 07b_fgsea_salmon.sh (Wrapper) ---
cat > scripts/07b_fgsea_salmon.sh << EOF
#!/bin/bash
set -e
echo "--- FGSEA SALMON WRAPPER SCRIPT ---"
$(printf "$CONDA_ACTIVATION" "r_env" "r_env")
source 2_config.sh
echo "Executing R script..."
Rscript scripts/07b_fgsea_salmon.R
EOF

# --- 07b_fgsea_salmon.R (Pure R) ---
cat > scripts/07b_fgsea_salmon.R << 'EOF'
suppressPackageStartupMessages({
    library(fgsea)
    library(readr)
    library(dplyr)
    library(tibble)
    library(ggplot2)
})

# Read configuration from environment variables
GSEA_GENESETS_GMT <- Sys.getenv("GSEA_GENESETS_GMT")

# Define paths relative to project root
DESEQ_RESULTS_PATH <- "results/deseq2_salmon/deseq2_results.csv"
RESULTS_DIR_FGSEA <- "results/fgsea_salmon"

# --- ROBUSTNESS CHECK ---
if (!file.exists(DESEQ_RESULTS_PATH)) {
    write(paste("ERROR: Input DESeq2 results not found at:", DESEQ_RESULTS_PATH), stderr())
    q(status = 1)
}

dir.create(RESULTS_DIR_FGSEA, showWarnings = FALSE, recursive = TRUE)

deseq_res <- read_csv(DESEQ_RESULTS_PATH, show_col_types = FALSE)
ranked_stats <- deseq_res %>%
    filter(!is.na(stat)) %>%
    distinct(gene_id, .keep_all=TRUE) %>%
    pull(stat, name=gene_id)

pathways <- gmtPathways(GSEA_GENESETS_GMT)
fgsea_res <- fgsea(pathways, ranked_stats, minSize=15, maxSize=500)
fgsea_res_tbl <- fgsea_res %>% as_tibble() %>% arrange(desc(NES))
write_tsv(fgsea_res_tbl, file.path(RESULTS_DIR_FGSEA, "fgsea_results.tsv"))

top_pathways_up <- fgsea_res[ES > 0][head(order(pval), n=10), pathway]
top_pathways_down <- fgsea_res[ES < 0][head(order(pval), n=10), pathway]
top_pathways <- c(top_pathways_up, rev(top_pathways_down))
pdf(file.path(RESULTS_DIR_FGSEA, "fgsea_top_pathways.pdf"))
plotGseaTable(pathways[top_pathways], ranked_stats, fgsea_res, gseaParam=0.5)
dev.off()
print("fGSEA (Salmon) analysis complete.")
EOF

# --- 08_multiqc.sh ---
cat > scripts/08_multiqc.sh << EOF
#!/bin/bash
set -e
$(printf "$CONDA_ACTIVATION" "multiqc" "multiqc")
source 2_config.sh
echo "Running MultiQC..."
multiqc results --outdir results/multiqc
EOF

# --- Make all scripts executable ---
chmod +x scripts/*.sh

# =============================================================================
#
# GENERATE CONDOR SUBMIT FILES (condor_submit_files/) - Always overwritten
#
# =============================================================================
echo "-> Generating/updating Condor submit files..."

# --- Base submit file template ---
SUB_TEMPLATE="
universe       = vanilla
executable     = scripts/%s
arguments      = \$(sample_id)
request_cpus   = %s
request_memory = %s
output         = logs/%s/\$(sample_id).out
error          = logs/%s/\$(sample_id).err
log            = logs/condor.log
getenv         = True
# All jobs run from the main project directory
initialdir     = ${PWD}
queue
"
# --- No-args submit file template ---
SUB_TEMPLATE_NOARGS="
universe       = vanilla
executable     = scripts/%s
request_cpus   = %s
request_memory = %s
output         = logs/%s/%s.out
error          = logs/%s/%s.err
log            = logs/condor.log
getenv         = True
# All jobs run from the main project directory
initialdir     = ${PWD}
queue
"

printf "$SUB_TEMPLATE" "01_fastqc.sh" 2 4GB "fastqc" "fastqc" > condor_submit_files/fastqc.sub
printf "$SUB_TEMPLATE" "02_trim.sh" 4 8GB "trim" "trim" > condor_submit_files/trim.sub
printf "$SUB_TEMPLATE" "03a_star_align.sh" 8 32GB "star" "star" > condor_submit_files/star.sub
printf "$SUB_TEMPLATE" "03b_salmon_quant.sh" 8 16GB "salmon" "salmon" > condor_submit_files/salmon.sub
printf "$SUB_TEMPLATE" "04a_featurecounts.sh" 4 8GB "featurecounts" "featurecounts" > condor_submit_files/featurecounts.sub
# Point to the new .sh wrappers for R scripts
printf "$SUB_TEMPLATE_NOARGS" "05a_merge_counts_star.sh" 1 8GB "merge" "merge_star" "merge" "merge_star" > condor_submit_files/merge_counts.sub
printf "$SUB_TEMPLATE_NOARGS" "06a_deseq2_star.sh" 1 16GB "deseq2" "deseq2_star" "deseq2" "deseq2_star" > condor_submit_files/deseq2_star.sub
printf "$SUB_TEMPLATE_NOARGS" "06b_deseq2_salmon.sh" 1 16GB "deseq2" "deseq2_salmon" "deseq2" "deseq2_salmon" > condor_submit_files/deseq2_salmon.sub
printf "$SUB_TEMPLATE_NOARGS" "07a_fgsea_star.sh" 1 8GB "fgsea" "fgsea_star" "fgsea" "fgsea_star" > condor_submit_files/fgsea_star.sub
printf "$SUB_TEMPLATE_NOARGS" "07b_fgsea_salmon.sh" 1 8GB "fgsea" "fgsea_salmon" "fgsea" "fgsea_salmon" > condor_submit_files/fgsea_salmon.sub
printf "$SUB_TEMPLATE_NOARGS" "08_multiqc.sh" 1 4GB "multiqc" "multiqc" "multiqc" "multiqc" > condor_submit_files/multiqc.sub

# =============================================================================
#
# GENERATE WORKFLOW MANAGEMENT SCRIPTS - Always overwritten
#
# =============================================================================
echo "-> Generating/updating workflow management scripts..."

# --- generate_dag.sh ---
cat > generate_dag.sh << 'EOF'
#!/bin/bash
set -e
source 2_config.sh
if [[ "${ALIGNER}" != "STAR" && "${ALIGNER}" != "SALMON" ]]; then
    echo "Error: ALIGNER variable in 2_config.sh must be 'STAR' or 'SALMON'" >&2
    exit 1
fi
echo "Generating DAG for the ${ALIGNER} pipeline..."
DAG_FILE="rnaseq.dag"
rm -f ${DAG_FILE} ${DAG_FILE}.*
echo "# RNA-seq Analysis DAG for ${ALIGNER}" > ${DAG_FILE}
# Read sample IDs, explicitly skipping comments and the header line
SAMPLE_LIST=($(awk -F'\t' '!/^#/ && $1 != "sample_id" {print $1}' ${SAMPLE_SHEET}))
# --- Per-sample jobs ---
for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
    echo "JOB fastqc_${SAMPLE_ID} condor_submit_files/fastqc.sub" >> ${DAG_FILE}
    echo "VARS fastqc_${SAMPLE_ID} sample_id=\"${SAMPLE_ID}\"" >> ${DAG_FILE}
    echo "JOB trim_${SAMPLE_ID} condor_submit_files/trim.sub" >> ${DAG_FILE}
    echo "VARS trim_${SAMPLE_ID} sample_id=\"${SAMPLE_ID}\"" >> ${DAG_FILE}
done
# --- Alignment-specific jobs and dependencies ---
if [ "${ALIGNER}" == "STAR" ]; then
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo "JOB star_${SAMPLE_ID} condor_submit_files/star.sub" >> ${DAG_FILE}
        echo "VARS star_${SAMPLE_ID} sample_id=\"${SAMPLE_ID}\"" >> ${DAG_FILE}
        echo "JOB fcounts_${SAMPLE_ID} condor_submit_files/featurecounts.sub" >> ${DAG_FILE}
        echo "VARS fcounts_${SAMPLE_ID} sample_id=\"${SAMPLE_ID}\"" >> ${DAG_FILE}
    done
    echo "JOB merge_counts condor_submit_files/merge_counts.sub" >> ${DAG_FILE}
    echo "JOB deseq2_star condor_submit_files/deseq2_star.sub" >> ${DAG_FILE}
    echo "JOB fgsea_star condor_submit_files/fgsea_star.sub" >> ${DAG_FILE}
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo "PARENT fastqc_${SAMPLE_ID} CHILD trim_${SAMPLE_ID}" >> ${DAG_FILE}
        echo "PARENT trim_${SAMPLE_ID} CHILD star_${SAMPLE_ID}" >> ${DAG_FILE}
        echo "PARENT star_${SAMPLE_ID} CHILD fcounts_${SAMPLE_ID}" >> ${DAG_FILE}
    done
    FCOUNTS_DEPS=$(printf "fcounts_%s " "${SAMPLE_LIST[@]}")
    echo "PARENT ${FCOUNTS_DEPS} CHILD merge_counts" >> ${DAG_FILE}
    echo "PARENT merge_counts CHILD deseq2_star" >> ${DAG_FILE}
    echo "PARENT deseq2_star CHILD fgsea_star" >> ${DAG_FILE}
elif [ "${ALIGNER}" == "SALMON" ]; then
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo "JOB salmon_${SAMPLE_ID} condor_submit_files/salmon.sub" >> ${DAG_FILE}
        echo "VARS salmon_${SAMPLE_ID} sample_id=\"${SAMPLE_ID}\"" >> ${DAG_FILE}
    done
    echo "JOB deseq2_salmon condor_submit_files/deseq2_salmon.sub" >> ${DAG_FILE}
    echo "JOB fgsea_salmon condor_submit_files/fgsea_salmon.sub" >> ${DAG_FILE}
    for SAMPLE_ID in "${SAMPLE_LIST[@]}"; do
        echo "PARENT fastqc_${SAMPLE_ID} CHILD trim_${SAMPLE_ID}" >> ${DAG_FILE}
        echo "PARENT trim_${SAMPLE_ID} CHILD salmon_${SAMPLE_ID}" >> ${DAG_FILE}
    done
    SALMON_DEPS=$(printf "salmon_%s " "${SAMPLE_LIST[@]}")
    echo "PARENT ${SALMON_DEPS} CHILD deseq2_salmon" >> ${DAG_FILE}
    echo "PARENT deseq2_salmon CHILD fgsea_salmon" >> ${DAG_FILE}
fi
# --- Final MultiQC job ---
LAST_JOB_STAR="fgsea_star"
LAST_JOB_SALMON="fgsea_salmon"
ALL_SAMPLES_TRIMMED=$(printf "trim_%s " "${SAMPLE_LIST[@]}")
echo "JOB multiqc condor_submit_files/multiqc.sub" >> ${DAG_FILE}
if [ "${ALIGNER}" == "STAR" ]; then
    echo "PARENT ${ALL_SAMPLES_TRIMMED} ${LAST_JOB_STAR} CHILD multiqc" >> ${DAG_FILE}
else
    echo "PARENT ${ALL_SAMPLES_TRIMMED} ${LAST_JOB_SALMON} CHILD multiqc" >> ${DAG_FILE}
fi
echo "DAG generation complete: ${DAG_FILE}"
EOF

# --- run_dag_submission.sh ---
cat > run_dag_submission.sh << 'EOF'
#!/bin/bash
set -e
echo "--------------------------------------"
echo "--- RNA-seq Pipeline Submitter ---"
echo "--------------------------------------"
echo "[1/2] Generating DAG file..."
bash generate_dag.sh
echo ""
echo "[2/2] Submitting DAG to HTCondor..."
condor_submit_dag rnaseq.dag
echo ""
echo "Pipeline submitted successfully! Monitor progress with 'condor_q'."
echo "--------------------------------------"
EOF

chmod +x generate_dag.sh run_dag_submission.sh

echo ""
echo "✅ Pipeline generation complete!"
echo "➡️ Next Steps:"
echo "   1. Run this script: bash ./pipeline_generator.sh"
echo "   2. Edit '2_config.sh' with your system-specific paths."
echo "   3. Prepare your sample sheet and design in '3_samples.tsv'."
echo "   4. Run 'bash run_dag_submission.sh' to start the analysis."
