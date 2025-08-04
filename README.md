# RNA-seq Analysis Pipeline

This project presents an end-to-end RNA-seq data analysis pipeline, built on Ubuntu using WSL (Windows Subsystem for Linux), Conda, command-line bioinformatics tools, and R with Bioconductor packages. The workflow starts from raw FASTQ files and progresses through alignment, quality control, read counting, visualization, and differential expression analysis.

---

## Objective

To build a complete, reproducible pipeline for RNA-seq analysis using publicly available data (SRR390728), and analyze differential gene expression using DESeq2 in R.

---

## Tools and Technologies Used

### Environment
- **Operating System**: Ubuntu via WSL (Windows Subsystem for Linux)
- **Shell**: Bash
- **Conda**: For managing bioinformatics tools and environments

### Command-Line Tools (Installed via Conda)
- **FastQC**: Quality control for FASTQ files
- **Bowtie2**: Read alignment to reference genome
- **Samtools**: BAM file processing
- **SRA-tools (fasterq-dump)**: For downloading FASTQ files (alternative: wget)
  
### R and Bioconductor Packages
- `Rsamtools`, `GenomicAlignments`, `Gviz`: Alignment processing and visualization
- `DESeq2`: Differential gene expression analysis

---

## Directory Structure

rnaseq_pipeline/
├── data/ # Raw input FASTQ files (not included in repo)
├── ref/ # Reference genome and index files (excluded from Git)
├── results/
│ ├── qc/ # FastQC outputs
│ └── alignments/ # BAM, BAI, and SAM files
├── scripts/ # (Optional: Additional scripts)
├── RNAseq_analysis.R # Reads BAM, computes coverage, plots read distribution
├── gene_counting.R # Counts reads overlapping genes using GTF
├── deseq2_analysis.R # Performs DE analysis and visualizations
├── rnaseq_pipeline.Rproj # RStudio project file
└── .gitignore # Excludes large binary/reference files


---

## Data Acquisition

Raw data was downloaded from the SRA archive:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/008/SRR390728/SRR390728_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/008/SRR390728/SRR390728_2.fastq.gz

Pipeline Workflow
Step 1: Set up Environment
```bash
conda create -n rnaseq_pipeline -c bioconda -c conda-forge fastqc bowtie2 samtools sra-tools -y
conda activate rnaseq_pipeline

Step 2: Quality Control
```bash
fastqc data/SRR390728_1.fastq data/SRR390728_2.fastq -o results/qc/

Step 3: Read Alignment
```bash
bowtie2 -x ref/GRCh38_noalt_as \
  -1 data/SRR390728_1.fastq \
  -2 data/SRR390728_2.fastq \
  -S results/alignments/SRR390728.sam

Step 4: Convert and Sort BAM
```bash
samtools view -bS results/alignments/SRR390728.sam | samtools sort -o results/alignments/SRR390728_sorted.bam
samtools index results/alignments/SRR390728_sorted.bam

R-Based Analysis
Step 5: Launch RStudio via WSL
Open rnaseq_pipeline.Rproj inside WSL path and run the scripts in order:

1. RNAseq_analysis.R
   Reads the BAM file
   Computes coverage for chr1
   Plots genome-wide histogram using Gviz

2. gene_counting.R
   Reads aligned reads using GenomicAlignments
   Loads annotation GTF (GRCh38)
   Counts reads by gene

3. deseq2_analysis.R
   Creates a synthetic count matrix (or loads from counts)
   Runs DESeq2
   Produces MA plot and Volcano plot



