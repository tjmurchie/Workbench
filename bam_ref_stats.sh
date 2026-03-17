#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# bam_ref_stats.sh
# -----------------------------------------------------------------------------
# Summarise mapping support for one or more BAM files, either across all
# references in each BAM or for a specified set of contigs/references.
#
# This script is designed as a lightweight command-line QC/reporting utility for
# workflows where you want a quick table of:
#   - mapped read counts above a chosen MAPQ threshold
#   - coverage breadth (% covered bases)
#   - mean depth
#   - total covered bases
#   - reference length
#
# Typical use cases
# -----------------
# 1. Compare how strongly several samples map to a candidate reference genome.
# 2. Summarise one BAM against many specific contigs in a panel or bait set.
# 3. Generate a simple sortable TSV-like table for downstream plotting or QC.
#
# Behaviour notes
# ---------------
# - Only PRIMARY alignments are counted.
# - Secondary and supplementary alignments are excluded.
# - If --template is used, only read1 records are counted, which is often
#   useful when you want approximate template counts from paired-end libraries.
# - Coverage and depth are computed with samtools coverage.
# - If no references are provided, the script reports one aggregate line per BAM.
#
# Requirements
# ------------
# - bash
# - samtools
# - awk, find, sort
#
# This file is intentionally kept as a standalone shell script rather than a
# packaged tool so it can be dropped directly into existing command-line
# workflows.
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bam_ref_stats.sh [options] -i <bam_dir|bam_glob|bam_file> [-i ...] [--ref <contig>]... [--ref-fasta <ref.fa>]...

Options:
  -i, --input       Input BAMs: a directory, a glob, or a BAM file (repeatable)
  -q, --mapq        MAPQ threshold (default: 30)
  -t, --threads     Threads for samtools view/index (default: 8)
  --ref             Reference contig name to summarize (repeatable)
  --ref-fasta       Reference FASTA file; contig names taken from FASTA headers (repeatable)
  --template        Count templates for paired-end by counting only read1 (R1) records
  --no-index        Do not attempt to create .bai indexes
  -h, --help        Show this help

Output columns (sorted by mapped_qMAPQ descending):
  mapped_qXX  cov_pct  mean_depth  covered_bases  ref_len  sample  ref  bam

Notes:
  - Counts are PRIMARY alignments only (excludes secondary 0x100 and supplementary 0x800).
  - If you don't pass --ref/--ref-fasta, the script aggregates across all references in the BAM.
EOF
}

MAPQ=30
THREADS=8
NO_INDEX=0
TEMPLATE=0

inputs=()
refs=()

# Parse args (supports repeated -i/--input, --ref, --ref-fasta)
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) inputs+=("$2"); shift 2 ;;
    -q|--mapq) MAPQ="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    --ref) refs+=("$2"); shift 2 ;;
    --ref-fasta)
      fa="$2"
      if [[ ! -f "$fa" ]]; then
        echo "ERROR: ref fasta not found: $fa" >&2
        exit 2
      fi
      # add ALL contig names from this fasta (first token after '>')
      while IFS= read -r name; do refs+=("$name"); done < <(awk '/^>/{sub(/^>/,"",$1); print $1}' "$fa")
      shift 2
      ;;
    --template) TEMPLATE=1; shift ;;
    --no-index) NO_INDEX=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *)
      # allow bare inputs without -i
      inputs+=("$1"); shift ;;
  esac
done

if [[ ${#inputs[@]} -eq 0 ]]; then
  usage >&2
  exit 2
fi

# Expand inputs into a bam list
bam_list=()
for p in "${inputs[@]}"; do
  if [[ -d "$p" ]]; then
    while IFS= read -r f; do bam_list+=("$f"); done < <(find "$p" -type f -name "*.bam" ! -name "*.bai" | sort)
  elif [[ -f "$p" ]]; then
    bam_list+=("$p")
  else
    # treat as glob/pattern
    while IFS= read -r f; do bam_list+=("$f"); done < <(compgen -G "$p" | sort || true)
  fi
done

if [[ ${#bam_list[@]} -eq 0 ]]; then
  echo "ERROR: no BAMs found from inputs." >&2
  exit 2
fi

# Flags:
#  -F 0x900 excludes secondary(0x100) + supplementary(0x800) => primary only
#  -F 0x904 also excludes unmapped(0x4) => primary & mapped
# template mode adds -f 0x40 (read1 only)
VIEW_F_PRIMARY="0x900"
VIEW_F_MAPPED="0x904"

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

# Header
echo -e "mapped_q${MAPQ}\tcov_pct\tmean_depth\tcovered_bases\tref_len\tsample\tref\tbam" > "$tmp"

for bam in "${bam_list[@]}"; do
  # index if needed (for samtools coverage)
  if [[ "$NO_INDEX" -eq 0 ]]; then
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
      samtools index -@ "$THREADS" "$bam" >/dev/null
    fi
  fi

  sample="$(basename "$bam" .bam)"

  if [[ ${#refs[@]} -gt 0 ]]; then
    # one line per requested reference contig
    for ref in "${refs[@]}"; do
      if [[ "$TEMPLATE" -eq 1 ]]; then
        mapped_q=$(samtools view -@ "$THREADS" -c -f 0x40 -F "$VIEW_F_MAPPED" -q "$MAPQ" "$bam" "$ref" 2>/dev/null || echo 0)
      else
        mapped_q=$(samtools view -@ "$THREADS" -c -F "$VIEW_F_MAPPED" -q "$MAPQ" "$bam" "$ref" 2>/dev/null || echo 0)
      fi

      # coverage/depth for that ref; if ref absent, emit zeros
      read -r cov_pct mean_dp cov_bases ref_len < <(
        samtools coverage -q "$MAPQ" -r "$ref" "$bam" 2>/dev/null |
        awk 'NR==2{
               len=$3-$2+1;
               printf "%.2f %.6f %d %d\n", $6, $7, $5, len
             }
             END{
               if (NR<2) printf "0.00 0.000000 0 0\n"
             }'
      )

      echo -e "${mapped_q}\t${cov_pct}\t${mean_dp}\t${cov_bases}\t${ref_len}\t${sample}\t${ref}\t${bam}" >> "$tmp"
    done
  else
    # aggregate across all contigs in BAM
    if [[ "$TEMPLATE" -eq 1 ]]; then
      mapped_q=$(samtools view -@ "$THREADS" -c -f 0x40 -F "$VIEW_F_MAPPED" -q "$MAPQ" "$bam" 2>/dev/null || echo 0)
    else
      mapped_q=$(samtools view -@ "$THREADS" -c -F "$VIEW_F_MAPPED" -q "$MAPQ" "$bam" 2>/dev/null || echo 0)
    fi

    read -r cov_pct mean_dp cov_bases ref_len < <(
      samtools coverage -q "$MAPQ" "$bam" 2>/dev/null |
      awk 'NR>1{
             len=$3-$2+1;
             totlen+=len;
             covbases+=$5;
             depthsum+=$7*len
           }
           END{
             covpct=(totlen?covbases/totlen*100:0);
             meandp=(totlen?depthsum/totlen:0);
             printf "%.2f %.6f %d %d\n", covpct, meandp, covbases, totlen
           }'
    )

    echo -e "${mapped_q}\t${cov_pct}\t${mean_dp}\t${cov_bases}\t${ref_len}\t${sample}\tALL\t${bam}" >> "$tmp"
  fi
done

# Sort by mapped reads desc (skip header)
{ head -n 1 "$tmp"; tail -n +2 "$tmp" | sort -t$'\t' -k1,1nr; }
