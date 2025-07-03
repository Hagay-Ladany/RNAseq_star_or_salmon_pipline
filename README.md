# **RNA-seq HTCondor Pipeline**

This document provides instructions for setting up and running a flexible RNA-seq analysis pipeline using HTCondor. The pipeline has been simplified to use a single file for sample and design configuration.
*note, the R scripts of deseq2 and fgsea are only serve as a template, interative R run over the merged_gene_counts.csv is still needed.

## **Table of Contents**

1. [Environment Setup](## **1\. Environment Setup (already done!)**)  
2. [Pipeline Overview](https://www.google.com/search?q=%23pipeline-overview)  
3. [Initial Setup: Running the Generator](https://www.google.com/search?q=%23initial-setup-running-the-generator)  
4. [Configuration: The 3\_samples.tsv File](https://www.google.com/search?q=%23configuration-the-3_samplestsv-file)  
5. [Running the Pipeline](https://www.google.com/search?q=%23running-the-pipeline)  
6. [Output Structure](https://www.google.com/search?q=%23output-structure)  
7. [Troubleshooting & FAQ](https://www.google.com/search?q=%23troubleshooting--faq)

## **1\. Environment Setup (already done!)**

###### *to avoid conflicts*: the eviroments packages where installed one by one without specifieng the version.

To ensure maximum stability and avoid conflicts, we will use four separate, specialized Conda environments.

Environment 1: Quality Control (qc\_env)  
This environment contains tools for QC and trimming.  
mamba create \-n qc\_env \-c bioconda \-c conda-forge \\  
  fastqc \\  
  trimmomatic

Environment 2: Alignment & Counting (align\_env)  
This environment contains the STAR aligner, featureCounts, and the Salmon aligner.  
mamba create \-n align\_env \-c bioconda \-c conda-forge \\  
  star \\  
  subread \\  
  salmon

Environment 3: R Post-Processing (r\_env)  
This environment contains R and the specific packages for statistical analysis.  
mamba create \-n r\_env \-c conda-forge \-c bioconda \\  
  r-base \\  
  r-data.table \\  
  bioconductor-deseq2 \\  
  r-readr \\  
  r-dplyr \\  
  r-tibble \\  
  r-ggplot2 \\  
  r-ggrepel \\  
  bioconductor-fgsea

Environment 4: Reporting (multiqc) \[MultiQC  did not like the qc\_env for some reason\]

  conda create \--name multiqc python=3.5  
  conda activate multiqc  
  conda install multiqc

The pipeline scripts will automatically activate the correct environment for each step.

to get the packages/apps version run the '*<b>et_env_data.sh</b>*' which result in a **<b>version_report.txt</b>*

## **2\. Pipeline Overview**

The workflow is managed by a DAG (Directed Acyclic Graph) submitted to an HTCondor cluster. The choice between the STAR or Salmon pathway is controlled by a single variable in the 2\_config.sh file.

### **Workflow for STAR**

FastQC \-\> Trimmomatic \-\> STAR \-\> featureCounts \-\> Merge Counts \-\> DESeq2 \-\> fGSEA \-\> MultiQC

### **Workflow for Salmon**

FastQC \-\> Trimmomatic \-\> Salmon \-\> tximport & DESeq2 \-\> fGSEA \-\> MultiQC

## **3\. Initial Setup: Running the Generator**

A generator script is provided to create the entire directory structure and all the necessary scripts for the pipeline.

**To run it, execute the following command in your terminal:**

bash ./pipeline\_generator.sh

This will create the scripts, condor\_submit\_files, logs, and results directories. It will also create template configuration files (2\_config.sh and 3\_samples.tsv). The script will ask for confirmation before overwriting these two files if they already exist.

## **4\. Configuration: The 3\_samples.tsv File**

This is now the **only file** you need to edit to describe your experiment. It is a tab-separated file that combines sample paths, metadata, and the DESeq2 design formula.

**Structure:**

1. **DESeq2 Directives (Top of File):**  
   * \#DESEQ\_DESIGN=: Specifies the design formula (e.g., \~condition or \~batch \+ condition). Use the column names from the header row.  
   * \#DESEQ\_CONTRAST=: Specifies the contrast for DESeq2. The format is factor,level\_to\_compare,base\_level (e.g., condition,treatment,control).  
2. **Header Row:**  
   * Must start with sample\_id, fastq\_1, and fastq\_2.  
   * Add any additional columns needed for your experimental design (e.g., condition, batch, timepoint).  
3. **Data Rows:**  
   * Each row corresponds to one sample, with values for each column defined in the header.

**Example 3\_samples.tsv:**

\#DESEQ\_DESIGN=\~batch \+ condition  
\#DESEQ\_CONTRAST=condition,treatment,control  
sample\_id    fastq\_1    fastq\_2    condition    batch  
sample\_A1    /path/to/A1\_R1.fastq.gz    /path/to/A1\_R2.fastq.gz    control    B1  
sample\_A2    /path/to/A2\_R1.fastq.gz    /path/to/A2\_R2.fastq.gz    control    B2  
sample\_B1    /path/to/B1\_R1.fastq.gz    /path/to/B1\_R2.fastq.gz    treatment    B1  
sample\_B2    /path/to/B2\_R1.fastq.gz    /path/to/B2\_R2.fastq.gz    treatment    B2

This single file now completely describes the experiment for the entire pipeline.

## **5\. Running the Pipeline**

After editing 2\_config.sh and 3\_samples.tsv, submit the pipeline to the HTCondor cluster with:

bash ./run\_dag\_submission.sh

You can monitor the progress of your jobs using condor\_q.

## **6\. Output Structure**

All results are saved in the results/ directory, with subdirectories for each step.

## **7\. Troubleshooting & FAQ**

**Q: A job failed and the pipeline stopped. How do I resume it without starting over?**

A: This is a key feature of the DAG workflow. You do not need to restart from scratch.

1. **Diagnose and Fix:** First, check the .err file for the failed job in the logs/ subdirectories to understand the error. For example, if fcounts\_sample\_A1 failed, check logs/featurecounts/sample\_A1.err. Fix the underlying issue (e.g., update a script via the generator or correct a path in 2\_config.sh).  

2. **Find the Rescue File:** In your main project directory, Condor will have created a rescue file named something like rnaseq.dag.rescue.001. Use the one with the highest number if there are multiple.  

3. **Submit the Rescue DAG:** Instead of the normal submission script, run condor\_submit\_dag directly on the rescue file:  
   condor\_submit\_dag rnaseq.dag.rescue.001
   
   The DAGMan will read this file, see which jobs already completed successfully, and intelligently resume the pipeline from the point of failure.

**Q: The featureCounts step failed with an error about "gene identifier attribute". What's wrong?**

A: This is a very common issue and means there is a mismatch between the GTF annotation file and the script.

* **The Cause:** The featureCounts script was told to group reads by a certain attribute (e.g., \-g gene\_name), but your GTF file uses a different attribute as its primary identifier (e.g., gene\_id).  
* **The Fix:** The generator script has been updated to use \-g gene\_id, which is the standard for Ensembl annotations. If you use a different annotation source in the future, you may need to check your GTF file and adjust the \-g parameter in scripts/04a\_featurecounts.sh accordingly.
