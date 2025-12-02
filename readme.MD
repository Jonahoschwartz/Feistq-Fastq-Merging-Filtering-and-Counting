# FASTQ Processing Pipeline

A pipeline for processing paired-end sequencing data through merging, filtering, trimming, and counting matched reads.

## Overview

This pipeline processes FASTQ files through a series of steps: unzipping, merging paired reads, quality filtering, trimming based on target sequences, and counting. Originally written by Phil Huss and heavily modified/maintained by JOS, it's designed to handle multiple pairs of sequencing files separately and works especially well with files from the Biotech Core's sequencing service.

**Note:** This code uses paths configured for the Scarcity server and includes a Condor-formatted submit file. While it can be adapted for other Condor servers, you'll need to modify the fastp path accordingly.

## Pipeline Steps

1. Unzip FASTQ files
2. Filter and merge R1 and R2 reads using Fastp
3. Remove reads with summed error >1 (using precompiled C code from Phil)
4. Trim reads based on exact 5'/3' sequence matches
5. Sort and count reads
6. Separate sequences into two files based on count threshold

## Setup Instructions

### 1. Clone and Prepare Repository
```bash
# Clone repo into your Scarcity folder
git clone [repo-url] [your-folder-name]
cd [your-folder-name]
```

### 2. Add Data Files
Place your FASTQ files in the `data/` folder. The pipeline accepts both direct uploads and symlinks from Data Cat. Keep the compiled `fastq_filter_reporter` executable in this folder.

### 3. Configure Submit File
Fill out the `process.submit` submit file with the following parameters:

- **FIVEPRIME_SEQUENCE**: The 15-base sequence immediately 5' of the variable region  
  Example: `GAAAGGAGGAACTAT`

- **BARCODE_CONSTANT**: The constant sequence between the variable region and barcodes  
  For most applications, use the same as THREEPRIME_SEQUENCE

- **THREEPRIME_SEQUENCE**: The 15-base sequence immediately 3' of the random barcode  
  Example: `GGTATGCAAACGCAA`

- **READ_THRESHOLD**: Minimum read count for high-confidence sequences  
  Example: `5`

**Important:** You may need to change the `../data` path to the absolute path of your data folder.

### 4. Prepare Headers File
Add your sample headers to `header.txt`, one per line. Use the library name column from the core submission file.

**Critical:** No header should be a substring of another. For example:
- ❌ Bad: `sample_1` and `sample_12`
- ✅ Good: `sample_1_` and `sample_12_`

### 5. Submit Job
```bash
condor_submit process.submit
```

## Output Structure

Each read pair generates its own directory containing:

| File | Description |
|------|-------------|
| `*.fastq.gz` | Original zipped FASTQ files |
| `*.html`, `*.json` | Fastp quality reports |
| `*_failed_reads.fastq` | Reads that failed Fastp filtering |
| `merged_reads.fastq` | Successfully merged reads from Fastp |
| `read_distribution.txt` | Frequency distribution of read counts |
| `trimmed_PCR_sequences.txt` | Reads after trimming primer sequences |
| `counted_PCR_sequences.txt` | Trimmed sequences with count data |
| `*_reads_above_threshold.txt` | High-confidence sequences (≥ threshold) |
| `*_reads_below_threshold.txt` | Low-confidence sequences (< threshold) |

Condor report files will be placed outside these directories.

## Troubleshooting

If `reads_below_threshold.txt` is empty, check the intermediate output files to identify where reads were filtered out during the pipeline.

---

**Credits:** Original pipeline desing by Phil Huss, modified and maintained by Jonah O'Mara Schwartz