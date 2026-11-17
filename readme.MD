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
| `*_fastp.html`, `*_fastp.json` | Fastp quality reports for that sample |
| `*_failed_reads.fastq` | Reads that failed Fastp filtering |
| `merged_reads.fastq` | Successfully merged reads from Fastp |
| `read_distribution.txt` | Frequency distribution of read counts |
| `trimmed_PCR_sequences.txt` | Reads after trimming primer sequences |
| `counted_PCR_sequences.txt` | Trimmed sequences with count data |
| `*_reads_above_threshold.txt` | High-confidence sequences (more reads than the threshold) |
| `*_reads_below_threshold.txt` | Low-confidence sequences (at or below the threshold) |

Condor report files (`*_processing.err`, `.log`, `.out`) are moved into each sample's directory along with its results.

## Packaging Results for Download

After all jobs finish (check with `condor_q`), run from the repo folder:

```bash
./prep_download.sh
```

This creates `download_package/` and `download_package.zip` containing:

| File | Description |
|------|-------------|
| `reads_above_threshold/` | Every sample's `*_reads_above_threshold.txt` |
| `fastp_reports/` | Every sample's fastp `.html` and `.json` report |
| `QC_summary.tsv` | One row per sample with read counts at each pipeline step (opens in Excel) |
| `README.txt` | Description of each QC column |

Samples that failed or never ran are still listed in `QC_summary.tsv`, with a note explaining what is missing.

Options: `-o NAME` for a different package name, `-d DIR` for a different data folder, `-i FILE` for a different header file, `-x` to leave out the fastp reports, `-f` to overwrite a previous package. Run `./prep_download.sh -h` for details.

Download the zip from your own computer with:

```bash
scp [user]@scarcity-ap-1.glbrc.org:/full/path/to/repo/download_package.zip .
```

## Resetting the Pipeline

To rerun the pipeline on the same data, run from the repo folder:

```bash
./reset_pipeline.sh -n    # dry run: lists every move and delete, changes nothing
./reset_pipeline.sh       # asks you to type "yes" before changing anything
condor_submit process.submit
```

This moves the original `*.fastq.gz` files out of each sample folder back into `data/`, then deletes the sample folders and any pipeline leftovers in `data/`. The `.fastq.gz` files, `fastq_filter_reporter`, and anything it does not recognise as pipeline output are never deleted. **All results are permanently deleted**, so run `./prep_download.sh` first if you want to keep them, and make sure no jobs are still running.

Options: `-n` dry run, `-y` skip the confirmation, `-d DIR` for a different data folder.

## Troubleshooting

If `reads_below_threshold.txt` is empty, check the intermediate output files to identify where reads were filtered out during the pipeline.

---

**Credits:** Original pipeline desing by Phil Huss, modified and maintained by Jonah O'Mara Schwartz