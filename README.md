# CountsToFPKM

This repository provides two R scripts for converting raw gene expression counts to FPKM values using gene lengths retrieved from Ensembl. A straightforward script for quick processing or a modular script for more complex workflows can be choosen.

## Background

FPKM (Fragments Per Kilobase of transcript per Million mapped reads) is a common method used to normalize RNA-seq count data, accounting for both sequencing depth and gene length. These scripts automate the conversion of raw counts into FPKM values by retrieving gene lengths from Ensembl.

## Repository Contents

### 1. CountsToFPKM_Simple.R

This script provides a quick solution for converting raw counts to FPKM values without complex setup.  
Appropriate for processing a single dataset.

- **Input**: A tab-separated file where:
  - Rows represent genes (Ensembl gene IDs as row names).
  - Columns represent samples (sample names as column headers).
- **Output**: A tab-separated file where:
  - Rows represent genes.
  - Columns represent FPKM values for each sample.

### 2. CountsToFPKM_Modular.R

This is a more flexible version with separate functions for each step of the pipeline.  
It’s designed for reuse, modification, or extended workflows.

- **Input**: Same format as above (tab-separated file with gene IDs as rows and samples as columns).
- **Output**: Same format as above (tab-separated file with FPKM values).
- **Functions**:
  - `read_raw_counts`: Reads the input raw counts file.
  - `get_gene_lengths`: Queries Ensembl for gene lengths using Ensembl gene IDs.
  - `convert_to_fpkm`: Calculates FPKM values using raw counts and gene lengths.
  - `run_pipeline`: Main function that runs the entire process.

## Example Input

Example input file (`ROSMAP_DLPFC_Counts.tsv`):

| Gene ID        | Sample_1 | Sample_2 | Sample_3 |
|----------------|----------|----------|----------|
| ENSG00000123415| 1500     | 1800     | 1750     |
| ENSG00000234567| 2300     | 2100     | 2200     |
| ENSG00000345678| 1100     | 900      | 950      |

- The first column contains Ensembl gene IDs.
- The remaining columns contain raw counts for each sample.

## Example Output

Example output file (`FPKM_output.tsv`):

| Gene ID        | Sample_1 | Sample_2 | Sample_3 |
|----------------|----------|----------|----------|
| ENSG00000123415| 5.23     | 6.14     | 5.90     |
| ENSG00000234567| 8.45     | 7.85     | 8.11     |
| ENSG00000345678| 3.25     | 2.89     | 3.01     |

- The output file contains the same genes and samples as the input file.
- The values represent the calculated FPKM for each gene in each sample.

## Requirements

- **R** (version 4.0 or later)
- Required libraries: `biomaRt`, `edgeR`

## Usage Instructions

### Running the Simple Script

1. Open `CountsToFPKM_Simple.R`.
2. Update the following lines with the paths to your input and output files:
   ```r
   counts_file <- "path/to/your/ROSMAP_DLPFC_Counts.tsv"
   output_file <- "path/to/your/FPKM_output.tsv"

