#!/bin/bash

# Function to display usage and script parameters
usage() {
    echo "Usage: $0 -f FIVEPRIME_SEQUENCE -b BARCODE_CONSTANT -t THREEPRIME_SEQUENCE -r READ_THRESHOLD -i HEADER"
    echo
    echo "Options:"
    echo "  -f FIVEPRIME_SEQUENCE   The five-prime sequence before the variable region (e.g., 'GAAAGGAGGAACTAT')"
    echo "                          Typically 15 bases immediately 5' of the random region."
    echo "  -b BARCODE_CONSTANT     The barcode constant between the variable region and the barcodes"
    echo "                          (e.g., 'ATAACGAGTCTTCGTACCGGAAGACTCCTTG'). Includes the BbsI sites for Rec3 constructs."
    echo "  -t THREEPRIME_SEQUENCE  The three-prime sequence after the barcodes (e.g., 'GGTATGCAAACGCAA')."
    echo "                          Typically 15 bases immediately 3' of the random barcode."
    echo "  -r READ_THRESHOLD       The read threshold (e.g., 5). Reads below this threshold will be removed"
    echo "                          prior to downstream processing."
    echo "  -i HEADER               The file header for this experiment."
    echo "  -h                      Display this help message."
    exit 1
}

# Parse named parameters
while getopts "f:b:t:r:i:h" opt; do
  case $opt in
    f) FIVEPRIME_SEQUENCE="$OPTARG" ;;
    b) BARCODE_CONSTANT="$OPTARG" ;;
    t) THREEPRIME_SEQUENCE="$OPTARG" ;;
    r) READ_THRESHOLD="$OPTARG" ;;
    i) HEADER="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Check if all required variables are set
if [ -z "$FIVEPRIME_SEQUENCE" ] || \
   [ -z "$BARCODE_CONSTANT" ] || \
   [ -z "$THREEPRIME_SEQUENCE" ] || \
   [ -z "$READ_THRESHOLD" ] || \
   [ -z "$HEADER" ]; then
    usage
fi

# Print the variables to verify
echo "FIVEPRIME_SEQUENCE=${FIVEPRIME_SEQUENCE}"
echo "BARCODE_CONSTANT=${BARCODE_CONSTANT}"
echo "THREEPRIME_SEQUENCE=${THREEPRIME_SEQUENCE}"
echo "READ_THRESHOLD=${READ_THRESHOLD}"
echo "HEADER=${HEADER}"

# Capture the start time
start_time=$(date)

# Unzip fastq files
gunzip -fk "${HEADER}"*.gz

echo "Merging Started using fastp"

/opt/bifxapps/fastp-0.23.2/fastp \
    -i "${HEADER}"*R1_001.fastq \
    -I "${HEADER}"*R2_001.fastq \
    -m --merged_out "${HEADER}_merged_reads.fastq" \
    --failed_out "${HEADER}_failed_reads.fastq" \
    -c --overlap_len_require 4 \
    -j "${HEADER}_fastp.json" \
    -h "${HEADER}_fastp.html" \
    --thread 8

echo "Completed Merging"

# Take merge data from fastp
# Reports are named per sample so jobs running at the same time in the shared
# data folder cannot overwrite or pick up each other's fastp.json/fastp.html
JSON_FILE="${HEADER}_fastp.json"

# Extract total_reads from the after_filtering section under summary
total_reads=$(jq '.summary.after_filtering.total_reads // "No data found"' "$JSON_FILE")
echo "Total Reads: $total_reads"

# Extract and display filtering results
echo "Filtering Result:"
jq '.filtering_result // "No filtering results found"' "$JSON_FILE"

echo "Started Quality Filter"

# Quality filter to remove reads with summed error over 1
./fastq_filter_reporter "${HEADER}_merged_reads.fastq" "${HEADER}_quality_filtered_reads.fastq"

echo "Completed Quality Filter"

# Constants, input and output files
INPUT_FILE="${HEADER}_quality_filtered_reads.fastq"
TRIMMED_FILE="${HEADER}_trimmed_PCR_sequences.txt"
COUNTED_FILE="${HEADER}_counted_PCR_sequences.txt"
DISTRIBUTION_FILE="${HEADER}_read_distribution.txt"
BELOW_THRESHOLD_FILE="${HEADER}_reads_below_threshold.txt"
ABOVE_THRESHOLD_FILE="${HEADER}_reads_above_threshold.txt"

# Count total lines before processing (optional debug step)
total_lines_before=$(wc -l < "$INPUT_FILE")
echo "Total lines before processing: $total_lines_before"

# 🔁 Extract sequence between the 5' and 3' primers
awk -v five="$FIVEPRIME_SEQUENCE" -v three="$THREEPRIME_SEQUENCE" '
{
    pattern = five "[ACGTN]*" three
    if (match($0, pattern, m)) {
        print m[0]
    }
}' "$INPUT_FILE" > "$TRIMMED_FILE"

# Sort and count unique sequences
sort "$TRIMMED_FILE" | uniq -c | sort -nr > "$COUNTED_FILE"

# Extract numbers, sort them, and count occurrences
awk '{print $1}' "$COUNTED_FILE" | sort -n | uniq -c > "$DISTRIBUTION_FILE"

# Remove lines based on threshold
awk -v threshold="$READ_THRESHOLD" '
  {
    if ($1 <= threshold)
      print $0 > "'"$BELOW_THRESHOLD_FILE"'"
    else
      print $0
  }
' "$COUNTED_FILE" > "$ABOVE_THRESHOLD_FILE"

# Count the lines in each file
below_threshold_count=$(wc -l < "$BELOW_THRESHOLD_FILE")
above_threshold_count=$(wc -l < "$ABOVE_THRESHOLD_FILE")

# Print the results
echo "Sequences with reads at or below $READ_THRESHOLD: $below_threshold_count"
echo "Sequences with reads above $READ_THRESHOLD: $above_threshold_count"

# Create output directory and move results
mkdir -p "${HEADER}"

# Clean up raw fastq files
rm "${HEADER}"*_001.fastq 

# Move outputs into the results folder
mv "${HEADER}_failed_reads.fastq" "${HEADER}/"
mv "${HEADER}"* "${HEADER}/"

# Report script end time
end_time=$(date)
echo "Script started at: $start_time"
echo "Script ended at:   $end_time"

