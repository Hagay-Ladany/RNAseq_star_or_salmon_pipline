# **RNA-seq HTCondor Pipeline**

This document provides instructions for setting up and running a flexible RNA-seq analysis pipeline using HTCondor. The pipeline has been simplified to use a single file for sample and design configuration.
*note, the R scripts of deseq2 and fgsea are only serve as a template, interative R run over the merged_gene_counts.csv is still needed.

## **Table of Contents**

1. [Environment Setup](#1-environment-setup-already-done)  
2. [Pipeline Overview](#2-pipeline-overview)  
3. [Initial Setup: Running the Generator](#3-initial-setup-running-the-generator)  
4. [Configuration: The 3_samples.tsv File](#4-configuration-the-3_samplestsv-file)  
5. [Running the Pipeline](#5-running-the-pipeline)  
6. [Output Structure](#6-output-structure)  
7. [Troubleshooting & FAQ](#7-troubleshooting--faq)

## **1. Environment Setup (already done!)**

###### *to avoid conflicts*: the eviroments packages where installed one by one without specifieng the version.

To ensure maximum stability and avoid conflicts, we will use four separate, specialized Conda environments.

Environment 1: Quality Control (qc_env)  
This environment contains tools for QC and trimming.  
mamba create -n qc_env -c bioconda -c conda-forge \  
  fastqc \  
  trimmomatic

Environment 2: Alignment & Counting (align_env)  
This environment contains the STAR aligner, featureCounts, and the Salmon aligner.  
mamba create -n align_env -c bioconda -c conda-forge \  
  star \  
  subread \  
  salmon

Environment 3: R Post-Processing (r_env)  
This environment contains R and the specific packages for statistical analysis.  
mamba create -n r_env -c conda-forge -c bioconda \  
  r-base \  
  r-data.table \  
  bioconductor-deseq2 \  
  r-readr \  
  r-dplyr \  
  r-tibble \  
  r-ggplot2 \  
  r-ggrepel \  
  bioconductor-fgsea

Environment 4: Reporting (multiqc) [MultiQC  did not like the qc_env for some reason]

  conda create --name multiqc python=3.5  
  conda activate multiqc  
  conda install multiqc

The pipeline scripts will automatically activate the correct environment for each step.

to get the packages/apps version run the '*get_env_data.sh*' which result in a **version_report.txt**

## **2. Pipeline Overview**

The workflow is managed by a DAG (Directed Acyclic Graph) submitted to an HTCondor cluster. The choice between the STAR or Salmon pathway is controlled by a single variable in the 2_config.sh file.

### **Workflow for STAR**

FastQC -> Trimmomatic -> STAR -> featureCounts -> Merge Counts -> DESeq2 -> fGSEA -> MultiQC

### **Workflow for Salmon**

FastQC -> Trimmomatic -> Salmon -> tximport & DESeq2 -> fGSEA -> MultiQC

## **3. Initial Setup: Running the Generator**

A generator script is provided to create the entire directory structure and all the necessary scripts for the pipeline.

**To run it, execute the following command in your terminal:**

bash ./pipeline_generator.sh

This will create the scripts, condor_submit_files, logs, and results directories. It will also create template configuration files (2_config.sh and 3_samples.tsv). The script will ask for confirmation before overwriting these two files if they already exist.

## **4. Configuration: The 3_samples.tsv File**

This is now the **only file** you need to edit to describe your experiment. It is a tab-separated file that combines sample paths, metadata, and the DESeq2 design formula.

**Structure:**

1. **DESeq2 Directives (Top of File):**  
   * #DESEQ_DESIGN=: Specifies the design formula (e.g., ~condition or ~batch + condition). Use the column names from the header row.  
   * #DESEQ_CONTRAST=: Specifies the contrast for DESeq2. The format is factor,level_to_compare,base_level (e.g., condition,treatment,control).  
2. **Header Row:**  
   * Must start with sample_id, fastq_1, and fastq_2.  
   * Add any additional columns needed for your experimental design (e.g., condition, batch, timepoint).  
3. **Data Rows:**  
   * Each row corresponds to one sample, with values for each column defined in the header.

**Example 3_samples.tsv:**

#DESEQ_DESIGN=~batch + condition  
#DESEQ_CONTRAST=condition,treatment,control  
sample_id    fastq_1    fastq_2    condition    batch  
sample_A1    /path/to/A1_R1.fastq.gz    /path/to/A1_R2.fastq.gz    control    B1  
sample_A2    /path/to/A2_R1.fastq.gz    /path/to/A2_R2.fastq.gz    control    B2  
sample_B1    /path/to/B1_R1.fastq.gz    /path/to/B1_R2.fastq.gz    treatment    B1  
sample_B2    /path/to/B2_R1.fastq.gz    /path/to/B2_R2.fastq.gz    treatment    B2

This single file now completely describes the experiment for the entire pipeline.

## **5. Running the Pipeline**

After editing 2_config.sh and 3_samples.tsv, submit the pipeline to the HTCondor cluster with:

bash ./run_dag_submission.sh

You can monitor the progress of your jobs using condor_q.

## **6. Output Structure**

All results are saved in the results/ directory, with subdirectories for each step.

## **7. Troubleshooting & FAQ**

**Q: A job failed and the pipeline stopped. How do I resume it without starting over?**

A: This is a key feature of the DAG workflow. You do not need to restart from scratch.

1. **Diagnose and Fix:** First, check the .err file for the failed job in the logs/ subdirectories to understand the error. For example, if fcounts_sample_A1 failed, check logs/featurecounts/sample_A1.err. Fix the underlying issue (e.g., update a script via the generator or correct a path in 2_config.sh).  

2. **Find the Rescue File:** In your main project directory, Condor will have created a rescue file named something like rnaseq.dag.rescue.001. Use the one with the highest number if there are multiple.  

3. **Submit the Rescue DAG:** Instead of the normal submission script, run condor_submit_dag directly on the rescue file:  
   condor_submit_dag rnaseq.dag.rescue.001
   
   The DAGMan will read this file, see which jobs already completed successfully, and intelligently resume the pipeline from the point of failure.

**Q: The featureCounts step failed with an error about "gene identifier attribute". What's wrong?**

A: This is a very common issue and means there is a mismatch between the GTF annotation file and the script.

* **The Cause:** The featureCounts script was told to group reads by a certain attribute (e.g., -g gene_name), but your GTF file uses a different attribute as its primary identifier (e.g., gene_id).  
* **The Fix:** The generator script has been updated to use -g gene_id, which is the standard for Ensembl annotations. If you use a different annotation source in the future, you may need to check your GTF file and adjust the -g parameter in scripts/04a_featurecounts.sh accordingly.
