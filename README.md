# RNA-seq Analysis Pipeline

This repository contains an end-to-end RNA-seq analysis pipeline that processes raw sequencing data (FASTQ) through quality control, alignment, read counting, and differential expression analysis using R and Bioconductor tools.

---

## Project Summary

- **Data Source**: Public RNA-seq dataset from NCBI-SRA (SRR390728)
- **Platform**: Ubuntu via WSL (Windows Subsystem for Linux)
- **Environment**: Conda-based bioinformatics setup
- **Languages Used**: R and Bash
- **Focus**: Reproducible analysis, modular structure, clear visualization

---

## Tools and Packages Used

### Command-line tools (via Conda)

- **FastQC**: Quality control of FASTQ files
- **Bowtie2**: Alignment to reference genome
- **Samtools**: Format conversion, sorting, indexing
- **SRA-tools**: Downloading sequencing data

### R / Bioconductor packages

- `Rsamtools`, `GenomicAlignments`, `Gviz`: BAM file processing and read visualization
- `Rsubread`: Feature counting using GTF
- `DESeq2`: Differential expression analysis
- `BiocManager`: To install all required Bioconductor packages

---

## Project Structure
rnaseq_pipeline/

├── data/ # FASTQ input files

├── ref/ # Reference genome and annotation (excluded from Git)

├── results/

│ ├── qc/ # FASTQC results

│ └── alignments/ # SAM, BAM, BAI files

├── RNAseq_analysis.R # BAM read coverage visualization

├── gene_counting.R # Gene-level count matrix

├── deseq2_analysis.R # DESeq2-based differential expression

├── rnaseq_pipeline.Rproj # RStudio project file

└── .gitignore # Git ignore rules


## RNA-seq Pipeline – Step-by-Step Instructions

### 1. Create and Activate Conda Environment

    ```bash
    
    conda create -n rnaseq_pipeline -c bioconda -c conda-forge fastqc bowtie2 samtools sra-tools -y
    conda activate rnaseq_pipeline

### 2. Download Raw FASTQ Files

    ```bash
    
    cd ~/projects/rnaseq_pipeline/data

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/008/SRR390728/SRR390728_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/008/SRR390728/SRR390728_2.fastq.gz

### 3. Quality Check with FastQC

    ```bash
    
    mkdir -p results/qc
    fastqc data/SRR390728_1.fastq.gz data/SRR390728_2.fastq.gz -o results/qc/

### 4. Reference Genome Setup

**Place the following files in ref/:**
- Bowtie2 index: GRCh38_noalt_as.*.bt2
- GTF: Homo_sapiens.GRCh38.111.gtf

### 5. Align Reads Using Bowtie2

     ```bash
     
     mkdir -p results/alignments
     bowtie2 -x ref/GRCh38_noalt_as \
     -1 data/SRR390728_1.fastq.gz \
     -2 data/SRR390728_2.fastq.gz \
     -S results/alignments/SRR390728.sam
 
### 6. Convert SAM to Sorted BAM and Index

      ```bash

      samtools view -bS results/alignments/SRR390728.sam > results/alignments/SRR390728.bam
      samtools sort results/alignments/SRR390728.bam -o results/alignments/SRR390728_sorted.bam
      samtools index results/alignments/SRR390728_sorted.bam

### 7. BAM Read Coverage Visualization in R

       ```r
       
      # Install packages
        install.packages("BiocManager")
        BiocManager::install(c("Rsamtools", "GenomicAlignments", "Gviz"))

      # Load and plot coverage
        library(Rsamtools)
        library(GenomicAlignments)
        library(Gviz)

        bamfile <- "results/alignments/SRR390728_sorted.bam"
        reads <- readGAlignments(bamfile)
        cov <- coverage(reads)

      # Plot (example for chr1)
        chr1_cov <- cov[["chr1"]]
 

### 8. Generate Gene Counts with featureCounts

       ```r
       
       BiocManager::install("Rsubread")
       library(Rsubread)

       fc <- featureCounts(
       files = "results/alignments/SRR390728_sorted.bam",
       annot.ext = "ref/Homo_sapiens.GRCh38.111.gtf",
       isGTFAnnotationFile = TRUE,
       GTF.featureType = "exon",
       GTF.attrType = "gene_id",
       useMetaFeatures = TRUE
       )

       head(fc$counts)

### 9. Differential Expression Using DESeq2

       ```r
       
       BiocManager::install("DESeq2")
       library(DESeq2)

    # Example simulated count matrix
      counts <- matrix(rnbinom(1000, mu=10, size=1), ncol=6)
      colnames(counts) <- c("A1", "A2", "A3", "B1", "B2", "B3")

      coldata <- data.frame(
      row.names = colnames(counts),
      condition = factor(c("control", "control", "control", "treated", "treated", "treated"))
      )

      dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
      dds <- DESeq(dds)
      res <- results(dds)

    # MA Plot
      plotMA(res)

    # Volcano Plot
      plot(res$log2FoldChange, -log10(res$pvalue),
      pch = 20, col = ifelse(res$padj < 0.05, "red", "black"),
      xlab = "log2 Fold Change", ylab = "-log10 p-value",
      main = "Volcano Plot")











