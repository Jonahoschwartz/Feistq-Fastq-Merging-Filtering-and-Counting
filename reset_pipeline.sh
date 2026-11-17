#!/bin/bash
#
# reset_pipeline.sh
#
# Puts the data folder back to how it was before the pipeline ran, so the
# jobs can be resubmitted:
#   1. moves every original *.fastq.gz out of its sample folder back into data/
#   2. deletes the sample folders and everything the pipeline wrote in them
#   3. deletes pipeline leftovers sitting loose in data/ (failed or partial jobs)
#
# The *.fastq.gz files, fastq_filter_reporter and anything it does not
# recognise as pipeline output are never deleted.
#
# Run it from the repo folder (the one containing data/):
#   ./reset_pipeline.sh
# Save anything you want first (e.g. run ./prep_download.sh), and make sure no
# pipeline jobs are still running (condor_q).
#
# Options:
#   -d DATA_DIR   data folder to reset            (default: data)
#   -n            dry run: list what would happen, change nothing
#   -y            do not ask for confirmation
#   -h            show this help

usage() {
    sed -n '3,24p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-1}"
}

DATA_DIR="data"
DRY=0
YES=0

while getopts "d:nyh" opt; do
  case $opt in
    d) DATA_DIR="$OPTARG" ;;
    n) DRY=1 ;;
    y) YES=1 ;;
    h) usage 0 ;;
    *) usage 1 ;;
  esac
done

DATA_DIR="${DATA_DIR%/}"
if [ -z "$DATA_DIR" ] || [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: data folder '$DATA_DIR' not found. Run this from the repo folder or use -d." >&2
    exit 1
fi

# Files the pipeline writes, by ending. Only these are ever deleted from data/ itself.
OUTPUT_ENDINGS=(
    _merged_reads.fastq _failed_reads.fastq _quality_filtered_reads.fastq
    _trimmed_PCR_sequences.txt _counted_PCR_sequences.txt _read_distribution.txt
    _reads_above_threshold.txt _reads_below_threshold.txt
    _processing.err _processing.log _processing.out
    _fastp.json _fastp.html
    _001.fastq
)

shopt -s nullglob

move_src=()      # fastq.gz files to move back
remove_dirs=()   # sample folders to delete
remove_files=()  # loose leftovers in data/ to delete
skipped=0

for dir in "$DATA_DIR"/*/; do
    dir="${dir%/}"
    name=$(basename "$dir")

    # Only touch folders the pipeline made: they hold files named <folder>_*
    owned=( "$dir/${name}_"* )
    [ ${#owned[@]} -eq 0 ] && continue

    gz=( "$dir"/*.fastq.gz )
    if [ ${#gz[@]} -eq 0 ]; then
        echo "SKIPPING $name: no .fastq.gz inside, leaving the folder alone" >&2
        skipped=$((skipped + 1))
        continue
    fi

    clash=0
    for f in "${gz[@]}"; do
        target="$DATA_DIR/$(basename "$f")"
        if [ -e "$target" ] || [ -L "$target" ]; then
            echo "SKIPPING $name: $(basename "$f") already exists in $DATA_DIR/" >&2
            clash=1
        fi
    done
    if [ "$clash" -eq 1 ]; then
        skipped=$((skipped + 1))
        continue
    fi

    move_src+=( "${gz[@]}" )
    remove_dirs+=( "$dir" )
done

for ending in "${OUTPUT_ENDINGS[@]}"; do
    for f in "$DATA_DIR"/*"$ending"; do
        [ -f "$f" ] || [ -L "$f" ] || continue
        remove_files+=( "$f" )
    done
done
for f in "$DATA_DIR/fastp.json" "$DATA_DIR/fastp.html"; do
    [ -f "$f" ] && remove_files+=( "$f" )
done

shopt -u nullglob

if [ ${#move_src[@]} -eq 0 ] && [ ${#remove_files[@]} -eq 0 ]; then
    echo "Nothing to reset in '$DATA_DIR'."
    [ "$skipped" -gt 0 ] && echo "$skipped folder(s) were skipped, see messages above."
    exit 0
fi

echo "Reset plan for $(cd "$DATA_DIR" && pwd):"
echo "  move back into data folder : ${#move_src[@]} .fastq.gz file(s)"
echo "  delete sample folders      : ${#remove_dirs[@]}"
echo "  delete loose leftovers     : ${#remove_files[@]} file(s)"
[ "$skipped" -gt 0 ] && echo "  skipped folders            : $skipped (see messages above)"

if [ "$DRY" -eq 1 ]; then
    echo
    for f in "${move_src[@]}";     do echo "would move    $f"; done
    for d in "${remove_dirs[@]}";  do echo "would delete  $d/"; done
    for f in "${remove_files[@]}"; do echo "would delete  $f"; done
    echo
    echo "Dry run only, nothing was changed."
    exit 0
fi

if [ "$YES" -ne 1 ]; then
    echo
    echo "All pipeline results in those folders will be permanently deleted."
    read -r -p "Type 'yes' to continue: " answer
    if [ "$answer" != "yes" ]; then
        echo "Cancelled, nothing was changed."
        exit 1
    fi
fi

failed=0

# Move the originals back, then delete a folder only once it holds no .fastq.gz
for dir in "${remove_dirs[@]}"; do
    ok=1
    for f in "$dir"/*.fastq.gz; do
        mv -n -- "$f" "$DATA_DIR/" || ok=0
        if [ -e "$f" ] || [ -L "$f" ]; then ok=0; fi
    done
    if [ "$ok" -eq 1 ]; then
        rm -rf -- "$dir" || failed=$((failed + 1))
    else
        echo "ERROR: could not move every .fastq.gz out of $dir, folder kept" >&2
        failed=$((failed + 1))
    fi
done

for f in "${remove_files[@]}"; do
    rm -f -- "$f" || failed=$((failed + 1))
done

if [ "$failed" -gt 0 ]; then
    echo "Finished with $failed problem(s), see messages above." >&2
    exit 1
fi

echo "Done. ${#move_src[@]} .fastq.gz file(s) are back in $DATA_DIR/ and ${#remove_dirs[@]} sample folder(s) were removed."
echo "Resubmit with: condor_submit process.submit"
