#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# split_by_ref.sh
# -----------------------------------------------------------------------------
# Split one or more BAM files in the current directory into per-reference BAMs.
#
# Output structure:
#   OUTROOT/
#     <reference_1>/
#       REFERENCE_NAME.txt
#       sampleA.bam
#       sampleA.bam.bai
#       sampleB.bam
#       ...
#     <reference_2>/
#       ...
#
# This is useful when a BAM contains alignments to many references/contigs and
# you want each target separated for manual inspection, downstream plotting, or
# per-reference statistics.
#
# Design notes
# ------------
# - Only references with mapped reads (>0 in samtools idxstats) are considered.
# - If two different references collapse to the same sanitised folder name,
#   a short MD5 suffix is appended so they remain unique.
# - Existing output BAMs are skipped rather than overwritten.
# - Input BAMs are assumed to be coordinate-sorted, which is the common case for
#   mapping outputs; extracting a single reference preserves sort order.
#
# Configuration is done through environment variables rather than command-line
# flags to keep the script simple for batch use:
#   THREADS=16 OUTROOT=by_ref DRYRUN=1 ./split_by_ref.sh
#
# Requirements
# ------------
# - samtools
# - awk
# - sort
# - md5sum (coreutils)
# -----------------------------------------------------------------------------
set -euo pipefail

# Split each BAM into per-reference BAMs and place them under one folder per reference.
# Requirements: samtools, awk, sort, md5sum (coreutils)

THREADS="${THREADS:-8}"                 # override: THREADS=16 ./split_by_ref.sh
OUTROOT="${OUTROOT:-split_by_reference}" # override: OUTROOT=by_ref ./split_by_ref.sh
DRYRUN="${DRYRUN:-0}"                   # override: DRYRUN=1 ./split_by_ref.sh

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing required command: $1" >&2; exit 1; }
}

safe_name() {
  # Make reference name filesystem-safe
  # Keep alnum . _ - ; convert everything else to _
  echo "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

run() {
  if [[ "$DRYRUN" == "1" ]]; then
    echo "[DRYRUN] $*"
  else
    eval "$@"
  fi
}

need_cmd samtools
need_cmd awk
need_cmd sort
need_cmd md5sum

shopt -s nullglob
bams=( *.bam )
if (( ${#bams[@]} == 0 )); then
  echo "ERROR: No .bam files found in $(pwd)" >&2
  exit 1
fi

mkdir -p "$OUTROOT"
LOG="$OUTROOT/split.log"
echo "Logging to: $LOG"
echo "Started: $(date)" | tee "$LOG"
echo "THREADS=$THREADS OUTROOT=$OUTROOT DRYRUN=$DRYRUN" | tee -a "$LOG"
echo "Found ${#bams[@]} BAM(s)." | tee -a "$LOG"

# Ensure indexes exist, and collect a global reference list (only refs with mapped reads > 0)
tmp_refs="$(mktemp)"
trap 'rm -f "$tmp_refs"' EXIT

echo "Indexing BAMs if needed and collecting references..." | tee -a "$LOG"
for bam in "${bams[@]}"; do
  # Check for an index alongside BAM (samtools may create .bam.bai or .bai depending on version)
  has_index=0
  [[ -f "${bam}.bai" ]] && has_index=1
  [[ -f "${bam%.bam}.bai" ]] && has_index=1
  [[ -f "${bam}.csi" ]] && has_index=1
  [[ -f "${bam%.bam}.csi" ]] && has_index=1

  if [[ "$has_index" -eq 0 ]]; then
    echo "  Indexing: $bam" | tee -a "$LOG"
    run "samtools index -@ \"$THREADS\" \"$bam\""
  fi

  # Collect refs with mapped reads > 0 (exclude the '*' line)
  samtools idxstats "$bam" \
    | awk '$1!="*" && $3>0 {print $1}' >> "$tmp_refs"
done

# Unique reference names
refs_file="$OUTROOT/references.txt"
sort -u "$tmp_refs" > "$refs_file"
ref_count=$(wc -l < "$refs_file" | awk '{print $1}')
echo "Unique references with mapped reads (across all BAMs): $ref_count" | tee -a "$LOG"
echo "Reference list written to: $refs_file" | tee -a "$LOG"

# Create per-reference directories with collision-safe names
declare -A REF_DIR_MAP  # original_ref -> dir_name
declare -A DIR_ORIG_MAP # dir_name -> original_ref

echo "Creating reference folders..." | tee -a "$LOG"
while IFS= read -r ref; do
  [[ -z "$ref" ]] && continue
  dir="$(safe_name "$ref")"

  # Handle rare collisions (different refs sanitizing to same folder name)
  if [[ -n "${DIR_ORIG_MAP[$dir]+x}" && "${DIR_ORIG_MAP[$dir]}" != "$ref" ]]; then
    h="$(printf "%s" "$ref" | md5sum | awk '{print substr($1,1,8)}')"
    dir="${dir}_${h}"
  fi

  REF_DIR_MAP["$ref"]="$dir"
  DIR_ORIG_MAP["$dir"]="$ref"

  run "mkdir -p \"$OUTROOT/$dir\""
  # Keep a record of the exact reference name
  if [[ "$DRYRUN" == "0" ]]; then
    printf "%s\n" "$ref" > "$OUTROOT/$dir/REFERENCE_NAME.txt"
  else
    echo "[DRYRUN] write $OUTROOT/$dir/REFERENCE_NAME.txt"
  fi
done < "$refs_file"

# Split per BAM
echo "Splitting BAMs..." | tee -a "$LOG"
for bam in "${bams[@]}"; do
  sample="${bam%.bam}"

  # Only refs that actually have reads in THIS bam
  mapfile -t bam_refs < <(samtools idxstats "$bam" | awk '$1!="*" && $3>0 {print $1}')

  echo "  $bam -> ${#bam_refs[@]} reference(s) with reads" | tee -a "$LOG"

  for ref in "${bam_refs[@]}"; do
    dir="${REF_DIR_MAP[$ref]}"
    outbam="$OUTROOT/$dir/${sample}.bam"

    # Skip if already exists
    if [[ -f "$outbam" ]]; then
      echo "    [skip] exists: $outbam" | tee -a "$LOG"
      continue
    fi

    echo "    extracting ref='$ref' -> $outbam" | tee -a "$LOG"

    # Assumes input BAM is coordinate-sorted (typical for mapping outputs).
    # Filtering to one reference preserves coordinate order, so indexing works without re-sorting.
    run "samtools view -@ \"$THREADS\" -bh \"$bam\" \"$ref\" -o \"$outbam\""
    run "samtools index -@ \"$THREADS\" \"$outbam\""
  done
done

echo "Done: $(date)" | tee -a "$LOG"
echo "Outputs are in: $OUTROOT" | tee -a "$LOG"
echo "Tip: set DRYRUN=1 first to preview actions." | tee -a "$LOG"

