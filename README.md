# WES Analysis Pipeline (Somatic + Germline + SV)

## Overview

This is a bioinformatics pipeline designed to process Whole Exome Sequencing (WES) data.

I built this to simulate a full cancer genomics workflow, taking raw FASTQ files all the way to analysis-ready data. It runs three parallel tracks:
1.  **Somatic Mutation Calling** (Tumor vs. Normal)
2.  **Germline Haplotype Calling**
3.  **Structural Variant Calling**

### Key Features
* **Reproducible:** Fully containerized dependencies via Conda/Mamba.
* **Robust:** Includes a self-validating setup script that fetches real 1000 Genomes data (Father/Daughter pair) for end-to-end testing.
* **Standardized:** Implements GATK Best Practices.
* **Efficient:** Features dynamic resource profiling, automated intermediate file cleanup, and isolated environments for legacy tool integration.

## Pipeline Steps

1. **Preprocessing**
   - FastP: Adapter trimming and quality filtering.
   - BWA-MEM: Alignment to the GRCh38 reference genome.
   - GATK: Duplicate marking and Base Quality Score Recalibration (BQSR).

2. **Variant Calling**
   - Mutect2: Calls somatic SNVs and Indels by comparing Tumor samples against matched Normals.
   - HaplotypeCaller: Calls germline variants on the normal samples.
   - Manta: Detects structural variants (translocations, large deletions).

3. **Analysis**
   - vcf2maf: Converts the hierarchical VCF output into a tabular MAF (Mutation Annotation Format) file, making it easier to analyze in R or Python.

## How to Run

### 1. Setup environment
This pipeline uses mamba (or conda) to manage software versions.

mamba env create -f environment.yaml
conda activate wes_pipeline

### 2. Download data
setup.sh fetches the Reference Genome (GRCh38) and a specific Tumor/Normal pair (Father/Daughter) from the 1000 Genomes Project. It includes checksum validation to ensure file integrity.

chmod +x setup.sh
./setup.sh

### 3. Run it!
Run the pipeline locally or on a cluster.

# Run on a local machine with 4 cores
snakemake --cores 4 --use-conda --conda-frontend conda

## Outputs
After a successful run, check these files:

| File | Description |
| :--- | :--- |
| `results/qc/multiqc_report.html` | **Start Here.** Summary of alignment rates, duplication, and quality. |
| `results/maf/Tumor_1000G.maf` | **Analysis Ready.** Tabular file of somatic mutations (ready for Oncoplots). |
| `results/vcf/somatic/` | Raw/Filtered VCFs from Mutect2. |
| `results/sv/manta/` | Structural Variant VCFs. |

## Project Structure

- **Snakefile**: The workflow logic. Defines the rules and dependencies.
- **config/config.yaml**: Configuration for file paths and run parameters.
- **config/resources.yaml**: Hardware specifications (RAM/CPU) decoupled from the pipeline logic.
- **envs/**: Isolated environment files for specific tools.
