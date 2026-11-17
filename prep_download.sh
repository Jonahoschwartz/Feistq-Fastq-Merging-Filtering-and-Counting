#!/bin/bash
#
# prep_download.sh
#
# Gathers the results of the FASTQ merging/filtering/counting pipeline into one
# small zipped folder that is easy to download from the cluster.
#
# The package contains:
#   reads_above_threshold/   every <sample>_reads_above_threshold.txt
#   fastp_reports/           every <sample>_fastp.html and _fastp.json
#   QC_summary.tsv           one row per sample: read counts at each pipeline step
#   README.txt               what each QC column means
#
# Run it from the repo folder (the one containing data/) after all jobs finish:
#   ./prep_download.sh
#
# Options:
#   -d DATA_DIR   folder holding the per-sample result folders   (default: data)
#   -o NAME       name of the package folder and archive         (default: download_package)
#   -i HEADERS    header file used for the run, to flag samples that
#                 produced no output   (default: header.txt beside DATA_DIR)
#   -x            leave the fastp reports out (smaller archive)
#   -f            overwrite a package previously made by this script
#   -h            show this help

usage() {
    sed -n '3,24p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-1}"
}

DATA_DIR="data"
OUT_NAME="download_package"
HEADER_FILE=""
FORCE=0
REPORTS=1

while getopts "d:o:i:xfh" opt; do
  case $opt in
    d) DATA_DIR="$OPTARG" ;;
    o) OUT_NAME="$OPTARG" ;;
    i) HEADER_FILE="$OPTARG" ;;
    x) REPORTS=0 ;;
    f) FORCE=1 ;;
    h) usage 0 ;;
    *) usage 1 ;;
  esac
done

DATA_DIR="${DATA_DIR%/}"
OUT_NAME="${OUT_NAME%/}"
SUFFIX="_reads_above_threshold.txt"
[ -z "$HEADER_FILE" ] && HEADER_FILE="$(dirname "$DATA_DIR")/header.txt"

if [ -z "$DATA_DIR" ] || [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: data folder '$DATA_DIR' not found. Run this from the repo folder or use -d." >&2
    exit 1
fi
if [ -z "$OUT_NAME" ]; then
    echo "ERROR: -o needs a name." >&2
    exit 1
fi

# Refuse to clobber anything this script did not make
for existing in "$OUT_NAME" "$OUT_NAME.zip" "$OUT_NAME.tar.gz"; do
    [ -e "$existing" ] || continue
    if [ "$FORCE" -ne 1 ]; then
        echo "ERROR: '$existing' already exists. Use -f to overwrite or -o for a new name." >&2
        exit 1
    fi
    if [ -d "$existing" ] && [ ! -f "$existing/QC_summary.tsv" ]; then
        echo "ERROR: '$existing' was not made by this script; not deleting it. Pick another name with -o." >&2
        exit 1
    fi
done
rm -rf "$OUT_NAME" "$OUT_NAME.zip" "$OUT_NAME.tar.gz"

mkdir -p "$OUT_NAME/reads_above_threshold" || exit 1
if [ "$REPORTS" -eq 1 ]; then mkdir -p "$OUT_NAME/fastp_reports" || exit 1; fi
QC="$OUT_NAME/QC_summary.tsv"

# Pull one number out of the Condor .err/.out text (last match wins).
#   $1 = section header to wait for ("" = anywhere), $2 = line label, $3... = files
grab() {
    local section="$1" label="$2"
    shift 2
    if [ $# -eq 0 ]; then echo "NA"; return; fi
    awk -v section="$section" -v label="$label" '
        BEGIN { insec = (section == "") ; val = "NA" }
        { sub(/\r$/, "") }
        section != "" && index($0, section) == 1 { insec = 1; next }
        insec && index($0, label) == 1 {
            v = substr($0, length(label) + 1)
            gsub(/[^0-9.]/, "", v)
            if (v != "") val = v
            if (section != "") insec = 0
        }
        END { print val }
    ' "$@" 2>/dev/null
}

# "<unique sequences> <total reads>" for a "count sequence" file; "NA NA" if missing
tally() {
    if [ -f "$1" ]; then
        awk 'NF { n++; s += $1 } END { printf "%.0f %.0f\n", n, s }' "$1"
    else
        echo "NA NA"
    fi
}

pct() {  # percent of $1 out of $2, or NA
    awk -v a="$1" -v b="$2" 'BEGIN {
        if (a == "NA" || b == "NA" || b + 0 == 0) print "NA"
        else printf "%.2f\n", 100 * a / b
    }'
}

# Copy a sample's fastp reports into the package and say what happened:
#   yes = copied, no = none found, "wrong sample" = report is from another sample.
# fastp records its command line in both reports, so each one is checked against
# the sample name. Runs made before the pipeline named reports per sample could
# end up holding another sample's report; those are left out.
copy_fastp_reports() {
    local sample="$1" dir="$2" ext f copied=0 wrong=0
    for ext in html json; do
        f="$dir/${sample}_fastp.$ext"
        [ -f "$f" ] || continue
        if grep -Fq -- "--merged_out ${sample}_merged_reads.fastq " "$f"; then
            if [ "$REPORTS" -eq 1 ]; then cp -p "$f" "$OUT_NAME/fastp_reports/" || exit 1; fi
            copied=$((copied + 1))
        else
            wrong=$((wrong + 1))
        fi
    done
    if   [ "$wrong" -gt 0 ];  then echo "wrong sample"
    elif [ "$copied" -gt 0 ] && [ "$REPORTS" -eq 1 ]; then echo "yes"
    elif [ "$copied" -gt 0 ]; then echo "left out (-x)"
    else echo "no"
    fi
}

n_ok=0
n_problem=0
n_wrong=0
seen=""

# Write one QC row.  $1 = sample, $2 = folder with its files, $3 = above-threshold file ("" if none), $4 = note
qc_row() {
    local sample="$1" dir="$2" above="$3" notes="$4"
    local logs=() ext lines
    local raw_pairs passed lowq manyN short merged threshold qfilt dup insert report
    local uniq_all reads_matched uniq_above reads_above uniq_below reads_below

    # Condor logs: fastp writes its summary to stderr, the pipeline echoes to stdout.
    # They normally sit in the sample folder but can be left in the data folder.
    for ext in err out; do
        if   [ -f "$dir/${sample}_processing.$ext" ];      then logs+=( "$dir/${sample}_processing.$ext" )
        elif [ -f "$DATA_DIR/${sample}_processing.$ext" ]; then logs+=( "$DATA_DIR/${sample}_processing.$ext" )
        else notes="${notes}missing .$ext; "
        fi
    done

    raw_pairs=$(grab "Read1 before filtering:" "total reads:" "${logs[@]}")
    passed=$(grab "" "reads passed filter:" "${logs[@]}")
    lowq=$(grab "" "reads failed due to low quality:" "${logs[@]}")
    manyN=$(grab "" "reads failed due to too many N:" "${logs[@]}")
    short=$(grab "" "reads failed due to too short:" "${logs[@]}")
    merged=$(grab "" "Read pairs merged:" "${logs[@]}")
    dup=$(grab "" "Duplication rate:" "${logs[@]}")
    insert=$(grab "" "Insert size peak (evaluated by paired-end reads):" "${logs[@]}")
    threshold=$(grab "" "READ_THRESHOLD=" "${logs[@]}")

    report=$(copy_fastp_reports "$sample" "$dir")
    if [ "$report" = "wrong sample" ]; then
        notes="${notes}fastp report in folder is from another sample (old pipeline run), left out; "
        n_wrong=$((n_wrong + 1))
    fi

    # Reads surviving the summed-error filter: prefer the filter's own report,
    # fall back to the line count the pipeline echoes (4 lines per read).
    qfilt=$(grab "" "Total reads passed filter:" "${logs[@]}")
    if [ "$qfilt" = "NA" ]; then
        lines=$(grab "" "Total lines before processing:" "${logs[@]}")
        [ "$lines" != "NA" ] && qfilt=$(( ${lines%%.*} / 4 ))
    fi

    read -r uniq_all reads_matched < <(tally "$dir/${sample}_counted_PCR_sequences.txt")
    if [ -n "$above" ]; then
        read -r uniq_above reads_above < <(tally "$above")
        if [ -f "$dir/${sample}_reads_below_threshold.txt" ]; then
            read -r uniq_below reads_below < <(tally "$dir/${sample}_reads_below_threshold.txt")
        else
            # the pipeline only creates this file if something falls below threshold
            uniq_below=0; reads_below=0
        fi
        [ "$uniq_above" = "0" ] && notes="${notes}no sequences above threshold; "
    else
        uniq_above=NA; reads_above=NA; uniq_below=NA; reads_below=NA
    fi
    [ "$raw_pairs" = "NA" ] && [ ${#logs[@]} -gt 0 ] && notes="${notes}fastp summary not found in logs; "

    printf '%s\t' "$sample" "$raw_pairs" "$passed" "$lowq" "$manyN" "$short" \
        "$merged" "$(pct "$merged" "$raw_pairs")" "$dup" "$insert" "$qfilt" \
        "$reads_matched" "$uniq_all" "$threshold" \
        "$uniq_above" "$reads_above" "$uniq_below" "$reads_below" \
        "$(pct "$reads_above" "$raw_pairs")" "$report" >> "$QC"
    printf '%s\n' "${notes%; }" >> "$QC"

    seen="${seen}${sample}
"
}

printf '%s\t' sample raw_read_pairs fastp_reads_passed fastp_low_quality fastp_too_many_N \
    fastp_too_short merged_reads pct_pairs_merged duplication_pct insert_size_peak \
    quality_filtered_reads \
    reads_with_primer_match unique_sequences read_threshold \
    unique_above_threshold reads_above_threshold \
    unique_below_threshold reads_below_threshold pct_raw_pairs_in_final \
    fastp_report > "$QC"
printf 'notes\n' >> "$QC"

# Every sample folder the pipeline made
shopt -s nullglob
for dir in "$DATA_DIR"/*/; do
    dir="${dir%/}"
    name=$(basename "$dir")
    above_files=( "$dir"/*"$SUFFIX" )

    if [ ${#above_files[@]} -eq 0 ]; then
        # Only report it if it looks like a pipeline folder (holds <name>_* files)
        leftovers=( "$dir/${name}_"* )
        if [ ${#leftovers[@]} -gt 0 ]; then
            echo "WARNING: $name has no *$SUFFIX file (job failed or still running?)" >&2
            qc_row "$name" "$dir" "" "NO reads_above_threshold file; "
            n_problem=$((n_problem + 1))
        fi
        continue
    fi

    for above in "${above_files[@]}"; do
        base=$(basename "$above")
        sample="${base%"$SUFFIX"}"
        cp -p "$above" "$OUT_NAME/reads_above_threshold/" || exit 1
        qc_row "$sample" "$dir" "$above" ""
        n_ok=$((n_ok + 1))
    done
done
shopt -u nullglob

# Samples that were submitted but never produced a folder
if [ -f "$HEADER_FILE" ]; then
    while IFS= read -r header || [ -n "$header" ]; do
        header="${header%$'\r'}"
        header="${header#"${header%%[![:space:]]*}"}"
        header="${header%"${header##*[![:space:]]}"}"
        [ -z "$header" ] && continue
        if ! printf '%s' "$seen" | grep -Fxq -- "$header"; then
            echo "WARNING: $header is in $HEADER_FILE but has no output folder" >&2
            qc_row "$header" "$DATA_DIR" "" "in header file but NO output folder; "
            n_problem=$((n_problem + 1))
        fi
    done < "$HEADER_FILE"
fi

if [ "$n_ok" -eq 0 ]; then
    echo "ERROR: no *$SUFFIX files found under '$DATA_DIR'/*/ . Nothing to package." >&2
    rm -rf "$OUT_NAME"
    exit 1
fi

cat > "$OUT_NAME/README.txt" <<EOF
Download package made by prep_download.sh on $(date)
Source: $(cd "$DATA_DIR" && pwd)
Samples packaged: $n_ok
Samples with no result file: $n_problem

reads_above_threshold/   one <sample>_reads_above_threshold.txt per sample
                         (columns: read count, sequence)
fastp_reports/           fastp's own report for each sample; open the .html in a browser
QC_summary.tsv           tab-separated, one row per sample (opens in Excel)

QC_summary.tsv columns, in pipeline order
  raw_read_pairs           read pairs in the original FASTQs (fastp "Read1 before filtering")
  fastp_reads_passed       individual reads (R1 + R2) passing fastp filters
  fastp_low_quality        individual reads failed by fastp for low quality
  fastp_too_many_N         individual reads failed by fastp for too many N
  fastp_too_short          individual reads failed by fastp for being too short
  merged_reads             read pairs fastp merged into one read
  pct_pairs_merged         merged_reads / raw_read_pairs
  duplication_pct          fastp's estimated duplication rate (%)
  insert_size_peak         most common insert size fastp saw (bases)
  quality_filtered_reads   merged reads kept by the summed-error (>1) filter
  reads_with_primer_match  reads with both the 5' and 3' sequences (kept after trimming)
  unique_sequences         distinct trimmed sequences
  read_threshold           READ_THRESHOLD used; sequences need MORE reads than this to pass
  unique_above_threshold   distinct sequences above threshold (rows in the delivered file)
  reads_above_threshold    reads those sequences account for
  unique_below_threshold   distinct sequences at or below threshold
  reads_below_threshold    reads those sequences account for
  pct_raw_pairs_in_final   reads_above_threshold / raw_read_pairs
  fastp_report             yes = in fastp_reports/; no = none found; "wrong sample" = the
                           report in the folder came from a different sample, so it was left out
  notes                    anything missing or odd for that sample

fastp numbers are read from each sample's own Condor _processing.err/.out files;
the counts from quality filtering onward are computed from the result files.
NA means the value could not be found.
EOF

# Zip it (fall back to tar.gz if zip is not installed)
parent=$(dirname "$OUT_NAME")
name=$(basename "$OUT_NAME")
if command -v zip >/dev/null 2>&1; then
    archive="$OUT_NAME.zip"
    ( cd "$parent" && zip -rq "$name.zip" "$name" ) || { echo "ERROR: zip failed" >&2; exit 1; }
else
    archive="$OUT_NAME.tar.gz"
    tar -czf "$archive" -C "$parent" "$name" || { echo "ERROR: tar failed" >&2; exit 1; }
fi

echo "Packaged $n_ok samples"
[ "$n_problem" -gt 0 ] && echo "  $n_problem sample(s) had no result file; see the notes column in QC_summary.tsv"
[ "$n_wrong" -gt 0 ] && echo "  $n_wrong sample(s) had a fastp report from another sample (old pipeline run); left out"
echo "  folder : $OUT_NAME/"
echo "  archive: $archive ($(du -h "$archive" | cut -f1))"
echo
echo "Download from your own computer with something like:"
echo "  scp $(whoami)@$(hostname -f 2>/dev/null || hostname):$(cd "$parent" && pwd)/$(basename "$archive") ."
