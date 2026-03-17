#!/usr/bin/env python3
"""
eager_read_summary.py
=====================
Extract per-sample read count summaries from nf-core/eager ancient DNA pipeline
output, suitable for manuscript supplementary materials (SOM) tables.

Overview
--------
This script parses the MultiQC data files that nf-core/eager produces at the end
of a run.  It walks the standard eager results/ directory layout, automatically
locates the MultiQC data directory, and pulls read counts from each processing
step into a single summary table.

No external dependencies are required beyond the Python 3.6+ standard library.

Pipeline steps tracked
----------------------
The columns in the output table follow the eager processing order:

  1. Raw_Reads          Total demultiplexed reads (R1 count x 2 for paired-end).
                        Source: multiqc_fastqc.txt (pre-trimming FastQC).

  2. After_PolyG_Trim   Reads surviving fastp poly-G tail removal.
                        Only present when fastp was run (e.g. NextSeq/NovaSeq data).
                        Source: results/FastP/*_fastp.json, or multiqc_general_stats.txt.

  3. After_Clipping     Reads retained after AdapterRemoval adapter trimming,
                        quality trimming, and paired-end read merging/collapsing.
                        Source: multiqc_adapter_removal.txt ("retained" column).

  4. %_Collapsed        Percentage of input reads that were successfully merged
                        into single collapsed reads by AdapterRemoval (typical for
                        short aDNA fragments where R1/R2 overlap).
                        Source: multiqc_adapter_removal.txt ("percent_collapsed").

  5. Mapped_Reads       Reads that aligned to the reference genome(s) before any
                        quality filtering.
                        Source: multiqc_samtools_flagstat.txt (pre-filter).

  6. Mapping_%          Percentage of input reads that mapped (mapped / total).
                        Source: multiqc_samtools_flagstat.txt ("mapped_passed_pct").

  7. After_Qual_Filter  Reads remaining after samtools quality/length filtering.
                        Only present when a post-filter flagstat file exists.
                        Source: multiqc_samtools_flagstat_1.txt (post-filter).

  8. Endogenous_%       Endogenous DNA content as reported by endorSpy.
                        Source: results/endorspy/*_endogenous_dna_mqc.json,
                        or multiqc_general_stats.txt.

  9. After_Dedup        Unique (non-duplicate) reads after PCR duplicate removal.
                        Calculated as: examined reads - duplicate reads.
                        Source: multiqc_picard_dups.txt (Picard MarkDuplicates)
                        or multiqc_dedup.txt (DeDup).

  10. Duplication_%     Percentage of reads that were PCR duplicates (0-100 scale).
                        Note: Picard stores this as a 0-1 fraction internally;
                        this script converts it to percent for consistency.
                        Source: multiqc_picard_dups.txt ("PERCENT_DUPLICATION").

  11. Mean_ReadLen      Mean read length of mapped reads after deduplication.
                        Source: multiqc_damageprofiler_metrics.txt.

  12. Median_ReadLen    Median read length of mapped reads after deduplication.
                        Source: multiqc_damageprofiler_metrics.txt.

Columns are automatically omitted if the corresponding tool was not run in a
given eager execution (e.g. no fastp column if poly-G trimming was skipped).

Expected directory layout
-------------------------
The script expects a standard nf-core/eager results/ directory:

    results/
    ├── multiqc/
    │   └── <run_name>_multiqc_report_data/   (or multiqc_data/)
    │       ├── multiqc_general_stats.txt
    │       ├── multiqc_fastqc.txt
    │       ├── multiqc_adapter_removal.txt
    │       ├── multiqc_samtools_flagstat.txt
    │       ├── multiqc_samtools_flagstat_1.txt
    │       ├── multiqc_picard_dups.txt
    │       ├── multiqc_damageprofiler_metrics.txt
    │       └── multiqc_sources.txt
    ├── FastP/          (or fastp/)
    │   └── *_fastp.json
    └── endorspy/
        └── *_endogenous_dna_mqc.json

You can point the script at the results/ directory itself, or at the parent
directory that contains results/.

Usage
-----
  # Print a formatted table to the terminal:
  eager_read_summary.py /path/to/eager/results

  # Save as a tab-separated file (for Excel, R, etc.):
  eager_read_summary.py /path/to/eager/results -o summary.tsv

  # Save as CSV:
  eager_read_summary.py /path/to/eager/results --csv -o summary.csv

  # Write TSV only, suppress the terminal table:
  eager_read_summary.py /path/to/eager/results --no-table -o summary.tsv

  # Show help:
  eager_read_summary.py -h

How pre- vs post-filter flagstat files are distinguished
--------------------------------------------------------
MultiQC numbers duplicate module outputs (e.g. multiqc_samtools_flagstat.txt and
multiqc_samtools_flagstat_1.txt).  This script determines which is pre-filter
and which is post-filter by:

  1. Checking multiqc_sources.txt for "pre-samtools filter" / "post-samtools
     filter" module labels (the method eager uses internally).
  2. Falling back to a heuristic: the file with larger total read counts across
     samples is treated as the pre-filter file.

If only one flagstat file exists, it is treated as pre-filter and the
After_Qual_Filter column is omitted.

Sample name normalisation
-------------------------
MultiQC appends suffixes like _R1_001, _R2_001, and _L0_polyg to sample names
depending on which tool generated the data.  This script strips those suffixes
to merge data across tools into a single row per sample.
"""

import argparse
import csv
import glob
import json
import os
import sys
from collections import OrderedDict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def safe_int(val):
    """Convert a value to int, returning None on failure."""
    if val is None or val == '' or val == 'nan':
        return None
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def safe_float(val):
    """Convert a value to float, returning None on failure."""
    if val is None or val == '' or val == 'nan':
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def fmt_int(val):
    """Format an integer with comma separators, or 'NA'."""
    if val is None:
        return 'NA'
    return f'{val:,}'


def fmt_pct(val, decimals=2):
    """Format a float as a percentage string, or 'NA'."""
    if val is None:
        return 'NA'
    return f'{val:.{decimals}f}'


def read_tsv(filepath):
    """Read a tab-delimited file and return a list of row dicts."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return list(reader)


def normalize_sample_name(name):
    """Strip common eager/MultiQC suffixes to recover the base sample name.

    Suffixes removed: _R1_001, _R2_001, _R1, _R2, _L0_polyg
    """
    for suffix in ('_R1_001', '_R2_001', '_R1', '_R2', '_L0_polyg'):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


# ---------------------------------------------------------------------------
# Auto-detect the MultiQC data directory
# ---------------------------------------------------------------------------

def find_multiqc_data_dir(results_dir):
    """Locate the MultiQC data directory inside an eager results folder.

    Tries several common naming patterns:
      results/multiqc/*_multiqc_report_data/
      results/multiqc/multiqc_data/
    Also handles being pointed directly at the data dir or at the parent of
    results/.
    """
    candidates = glob.glob(os.path.join(results_dir, 'multiqc', '*multiqc*data*'))
    if candidates:
        return candidates[0]

    candidates = glob.glob(os.path.join(results_dir, 'multiqc', 'multiqc_data'))
    if candidates:
        return candidates[0]

    # User pointed directly at the data dir
    if os.path.isfile(os.path.join(results_dir, 'multiqc_general_stats.txt')):
        return results_dir

    # User pointed at the parent containing results/
    parent = os.path.join(results_dir, '..')
    candidates = glob.glob(os.path.join(parent, 'multiqc', '*multiqc*data*'))
    if candidates:
        return candidates[0]

    return None


# ---------------------------------------------------------------------------
# Identify which flagstat file is pre- vs post-filter
# ---------------------------------------------------------------------------

def identify_flagstat_files(data_dir):
    """Return (pre_filter_path, post_filter_path) for samtools flagstat.

    See module docstring for the detection strategy.
    """
    base = os.path.join(data_dir, 'multiqc_samtools_flagstat')
    file_a = base + '.txt'
    file_b = base + '_1.txt'

    if not os.path.exists(file_a):
        return None, None
    if not os.path.exists(file_b):
        return file_a, None

    # Strategy 1: check multiqc_sources.txt for module labels
    sources_path = os.path.join(data_dir, 'multiqc_sources.txt')
    if os.path.exists(sources_path):
        pre_module = False
        post_module = False
        with open(sources_path) as fh:
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                module = parts[0]
                if 'pre-samtools filter' in module or 'pre_samtools_filter' in module:
                    pre_module = True
                if 'post-samtools filter' in module or 'post_samtools_filter' in module:
                    post_module = True
        if pre_module and post_module:
            # eager writes pre-filter first -> file_a, post-filter second -> file_b
            return file_a, file_b

    # Strategy 2: the file with larger total_passed sums is pre-filter
    rows_a = read_tsv(file_a)
    rows_b = read_tsv(file_b)
    sum_a = sum(safe_int(r.get('total_passed', 0)) or 0 for r in rows_a)
    sum_b = sum(safe_int(r.get('total_passed', 0)) or 0 for r in rows_b)

    if sum_a >= sum_b:
        return file_a, file_b
    else:
        return file_b, file_a


# ---------------------------------------------------------------------------
# Parsers for each MultiQC module file
# ---------------------------------------------------------------------------

def parse_fastqc_raw(data_dir):
    """Get raw read counts from the pre-trimming FastQC MultiQC file.

    Returns {sample: raw_total_reads} where raw_total_reads = R1 count * 2
    (paired-end).  Falls back to multiqc_fastqc_1.txt if needed.
    """
    results = {}
    for fname in ('multiqc_fastqc.txt', 'multiqc_fastqc_1.txt'):
        filepath = os.path.join(data_dir, fname)
        rows = read_tsv(filepath)
        for row in rows:
            name = row.get('Sample', '')
            if '_R1' not in name:
                continue
            sample = normalize_sample_name(name)
            total = safe_int(row.get('Total Sequences'))
            if total is not None and sample not in results:
                results[sample] = total * 2
        if results:
            break
    return results


def parse_fastp(results_dir, data_dir):
    """Get fastp poly-G trimming stats.

    Strategy 1: Parse individual fastp JSON files from results/FastP/ (or fastp/).
    Strategy 2: Fall back to multiqc_general_stats.txt (_L0_polyg rows).

    Returns {sample: {'before': N, 'after': N, 'pct_surviving': F}}
    """
    results = {}

    # Strategy 1: individual JSON files
    fastp_dir = os.path.join(results_dir, 'FastP')
    if not os.path.isdir(fastp_dir):
        fastp_dir = os.path.join(results_dir, 'fastp')
    if os.path.isdir(fastp_dir):
        for jf in glob.glob(os.path.join(fastp_dir, '*_fastp.json')):
            try:
                with open(jf) as fh:
                    data = json.load(fh)
                basename = os.path.basename(jf)
                sample = basename.replace('_fastp.json', '')
                sample = normalize_sample_name(sample)

                before = data.get('summary', {}).get('before_filtering', {}).get('total_reads')
                after = data.get('summary', {}).get('after_filtering', {}).get('total_reads')
                pct = None
                if before and after and before > 0:
                    pct = (after / before) * 100.0
                results[sample] = {
                    'before': safe_int(before),
                    'after': safe_int(after),
                    'pct_surviving': pct,
                }
            except (json.JSONDecodeError, KeyError):
                continue
        if results:
            return results

    # Strategy 2: multiqc_general_stats.txt
    gs_path = os.path.join(data_dir, 'multiqc_general_stats.txt')
    rows = read_tsv(gs_path)
    for row in rows:
        name = row.get('Sample', '')
        if '_L0_polyg' not in name and '_polyg' not in name:
            continue
        sample = normalize_sample_name(name)
        passed = None
        pct = None
        for col, val in row.items():
            if 'filtering_result_passed_filter_reads' in col:
                passed = safe_int(val)
            if 'pct_surviving' in col:
                pct = safe_float(val)
        if passed is not None:
            results[sample] = {
                'before': None,
                'after': passed,
                'pct_surviving': pct,
            }

    return results


def parse_adapter_removal(data_dir):
    """Parse AdapterRemoval trimming/merging statistics.

    Returns {sample: {reads_in, retained, pct_collapsed, pct_discarded, pct_aligned}}
    """
    results = {}
    filepath = os.path.join(data_dir, 'multiqc_adapter_removal.txt')
    for row in read_tsv(filepath):
        sample = normalize_sample_name(row.get('Sample', ''))
        results[sample] = {
            'reads_in': safe_int(row.get('reads_total')),
            'retained': safe_int(row.get('retained')),
            'pct_collapsed': safe_float(row.get('percent_collapsed')),
            'pct_discarded': safe_float(row.get('percent_discarded')),
            'pct_aligned': safe_float(row.get('percent_aligned')),
        }
    return results


def parse_flagstat(filepath):
    """Parse a samtools flagstat MultiQC file.

    Returns {sample: {total, mapped, mapped_pct}}
    """
    results = {}
    for row in read_tsv(filepath):
        sample = normalize_sample_name(row.get('Sample', ''))
        results[sample] = {
            'total': safe_int(row.get('total_passed')) or safe_int(row.get('flagstat_total')),
            'mapped': safe_int(row.get('mapped_passed')),
            'mapped_pct': safe_float(row.get('mapped_passed_pct')),
        }
    return results


def parse_dedup(data_dir):
    """Parse deduplication stats from Picard MarkDuplicates or DeDup.

    Tries Picard first (multiqc_picard_dups.txt), then DeDup (multiqc_dedup.txt).
    Picard's PERCENT_DUPLICATION is a 0-1 fraction and is converted to 0-100%.

    Returns {sample: {examined, duplicates, unique, pct_dup}}
    """
    results = {}

    # Try Picard MarkDuplicates
    filepath = os.path.join(data_dir, 'multiqc_picard_dups.txt')
    if os.path.exists(filepath):
        for row in read_tsv(filepath):
            sample = normalize_sample_name(row.get('Sample', ''))
            unpaired = safe_int(row.get('UNPAIRED_READS_EXAMINED')) or 0
            pairs = safe_int(row.get('READ_PAIRS_EXAMINED')) or 0
            examined = unpaired + (pairs * 2)

            unpaired_dups = safe_int(row.get('UNPAIRED_READ_DUPLICATES')) or 0
            pair_dups = safe_int(row.get('READ_PAIR_DUPLICATES')) or 0
            duplicates = unpaired_dups + (pair_dups * 2)

            unique = examined - duplicates
            pct = safe_float(row.get('PERCENT_DUPLICATION'))
            if pct is not None:
                pct = pct * 100.0  # Picard stores as 0-1 fraction
            results[sample] = {
                'examined': examined,
                'duplicates': duplicates,
                'unique': unique,
                'pct_dup': pct,
            }
        if results:
            return results

    # Try DeDup
    filepath = os.path.join(data_dir, 'multiqc_dedup.txt')
    if os.path.exists(filepath):
        for row in read_tsv(filepath):
            sample = normalize_sample_name(row.get('Sample', ''))
            total = safe_int(row.get('total_reads'))
            dups = safe_int(row.get('reverse_removed')) or 0
            unique = safe_int(row.get('forward_reads')) or (total - dups if total else None)
            pct = safe_float(row.get('dup_rate'))
            if pct is not None:
                pct = pct * 100.0  # DeDup also stores as fraction
            results[sample] = {
                'examined': total,
                'duplicates': dups,
                'unique': unique,
                'pct_dup': pct,
            }

    return results


def parse_endorspy(results_dir, data_dir):
    """Parse endogenous DNA percentage from endorSpy.

    Strategy 1: Parse individual endorspy JSON files from results/endorspy/.
    Strategy 2: Fall back to multiqc_general_stats.txt (base sample rows).

    Returns {sample: endogenous_pct}
    """
    results = {}

    # Strategy 1: individual JSON files
    endo_dir = os.path.join(results_dir, 'endorspy')
    if os.path.isdir(endo_dir):
        for jf in glob.glob(os.path.join(endo_dir, '*_endogenous_dna_mqc.json')):
            try:
                with open(jf) as fh:
                    data = json.load(fh)
                for sample, vals in data.get('data', {}).items():
                    sample = normalize_sample_name(sample)
                    results[sample] = safe_float(vals.get('endogenous_dna'))
            except (json.JSONDecodeError, KeyError):
                continue
        if results:
            return results

    # Strategy 2: general_stats
    gs_path = os.path.join(data_dir, 'multiqc_general_stats.txt')
    for row in read_tsv(gs_path):
        name = row.get('Sample', '')
        if '_R1' in name or '_R2' in name or '_polyg' in name:
            continue
        sample = normalize_sample_name(name)
        for col, val in row.items():
            if 'endorspy' in col and 'endogenous_dna' in col and 'post' not in col:
                v = safe_float(val)
                if v is not None:
                    results[sample] = v
                    break

    return results


def parse_damage_profiler(data_dir):
    """Parse DamageProfiler mean and median read lengths.

    Returns {sample: {mean_readlen, median_readlen}}
    """
    results = {}
    filepath = os.path.join(data_dir, 'multiqc_damageprofiler_metrics.txt')
    for row in read_tsv(filepath):
        sample = normalize_sample_name(row.get('Sample', ''))
        results[sample] = {
            'mean_readlen': safe_float(row.get('mean_readlength')),
            'median_readlen': safe_float(row.get('median')),
        }
    return results


# ---------------------------------------------------------------------------
# Main assembly
# ---------------------------------------------------------------------------

def build_summary(results_dir):
    """Assemble the full per-sample summary table.

    Reads all available MultiQC module files, merges them by normalised sample
    name, and returns (column_headers, rows_as_ordered_dicts).

    Columns for pipeline steps that were not run are automatically omitted.
    """
    data_dir = find_multiqc_data_dir(results_dir)
    if data_dir is None:
        print(f'Error: Could not find MultiQC data directory in {results_dir}',
              file=sys.stderr)
        print('Expected to find results/multiqc/*_multiqc_report_data/ or '
              'results/multiqc/multiqc_data/', file=sys.stderr)
        sys.exit(1)

    print(f'Reading MultiQC data from: {data_dir}', file=sys.stderr)

    # Parse all sources
    raw_reads = parse_fastqc_raw(data_dir)
    fastp_data = parse_fastp(results_dir, data_dir)
    ar_data = parse_adapter_removal(data_dir)
    pre_file, post_file = identify_flagstat_files(data_dir)
    pre_flag = parse_flagstat(pre_file) if pre_file else {}
    post_flag = parse_flagstat(post_file) if post_file else {}
    dedup_data = parse_dedup(data_dir)
    endorspy_data = parse_endorspy(results_dir, data_dir)
    damage_data = parse_damage_profiler(data_dir)

    # Collect all sample names across every source
    all_samples = set()
    for d in (raw_reads, fastp_data, ar_data, pre_flag, post_flag,
              dedup_data, endorspy_data, damage_data):
        all_samples.update(d.keys())
    samples = sorted(all_samples)

    if not samples:
        print('Error: No samples found. Check that the results directory is correct.',
              file=sys.stderr)
        sys.exit(1)

    print(f'Found {len(samples)} samples', file=sys.stderr)

    # Determine which optional steps were run
    has_fastp = bool(fastp_data)
    has_post_filter = bool(post_flag)
    has_dedup = bool(dedup_data)
    has_endorspy = bool(endorspy_data)
    has_damage = bool(damage_data)

    # Build column list (only include columns for steps that were executed)
    columns = [('sample', 'Sample')]
    columns.append(('raw_reads', 'Raw_Reads'))
    if has_fastp:
        columns.append(('after_fastp', 'After_PolyG_Trim'))
    columns.append(('after_clipping', 'After_Clipping'))
    columns.append(('pct_collapsed', '%_Collapsed'))
    columns.append(('mapped', 'Mapped_Reads'))
    columns.append(('mapping_pct', 'Mapping_%'))
    if has_post_filter:
        columns.append(('after_filter', 'After_Qual_Filter'))
    if has_endorspy:
        columns.append(('endogenous_pct', 'Endogenous_%'))
    if has_dedup:
        columns.append(('after_dedup', 'After_Dedup'))
        columns.append(('dup_pct', 'Duplication_%'))
    if has_damage:
        columns.append(('mean_readlen', 'Mean_ReadLen'))
        columns.append(('median_readlen', 'Median_ReadLen'))

    headers = [c[1] for c in columns]

    # Build one row per sample
    rows = []
    for sample in samples:
        row = OrderedDict()
        row['Sample'] = sample

        row['Raw_Reads'] = raw_reads.get(sample)

        if has_fastp:
            fd = fastp_data.get(sample, {})
            row['After_PolyG_Trim'] = fd.get('after')

        ar = ar_data.get(sample, {})
        row['After_Clipping'] = ar.get('retained')
        row['%_Collapsed'] = ar.get('pct_collapsed')

        pf = pre_flag.get(sample, {})
        row['Mapped_Reads'] = pf.get('mapped')
        row['Mapping_%'] = pf.get('mapped_pct')

        if has_post_filter:
            pof = post_flag.get(sample, {})
            row['After_Qual_Filter'] = pof.get('total')

        if has_endorspy:
            row['Endogenous_%'] = endorspy_data.get(sample)

        if has_dedup:
            dd = dedup_data.get(sample, {})
            row['After_Dedup'] = dd.get('unique')
            row['Duplication_%'] = dd.get('pct_dup')

        if has_damage:
            dm = damage_data.get(sample, {})
            row['Mean_ReadLen'] = dm.get('mean_readlen')
            row['Median_ReadLen'] = dm.get('median_readlen')

        rows.append(row)

    return headers, rows


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def print_table(headers, rows, file=sys.stdout):
    """Print a human-readable aligned table to the console."""
    int_cols = {'Raw_Reads', 'After_PolyG_Trim', 'After_Clipping',
                'Mapped_Reads', 'After_Qual_Filter', 'After_Dedup'}
    pct_cols = {'%_Collapsed', 'Mapping_%', 'Endogenous_%', 'Duplication_%'}
    float_cols = {'Mean_ReadLen', 'Median_ReadLen'}

    str_rows = []
    for row in rows:
        sr = OrderedDict()
        for h in headers:
            val = row.get(h)
            if h == 'Sample':
                sr[h] = str(val)
            elif h in int_cols:
                sr[h] = fmt_int(val)
            elif h in pct_cols:
                sr[h] = fmt_pct(val)
            elif h in float_cols:
                sr[h] = fmt_pct(val, 1) if val is not None else 'NA'
            else:
                sr[h] = str(val) if val is not None else 'NA'
        str_rows.append(sr)

    widths = {}
    for h in headers:
        widths[h] = len(h)
        for sr in str_rows:
            widths[h] = max(widths[h], len(sr.get(h, '')))

    hdr_line = '  '.join(h.rjust(widths[h]) if h != 'Sample' else h.ljust(widths[h])
                         for h in headers)
    print(hdr_line, file=file)
    print('-' * len(hdr_line), file=file)

    for sr in str_rows:
        line = '  '.join(sr[h].rjust(widths[h]) if h != 'Sample' else sr[h].ljust(widths[h])
                         for h in headers)
        print(line, file=file)


def write_delimited(headers, rows, filepath, delimiter='\t'):
    """Write the summary table as a TSV or CSV file."""
    with open(filepath, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter=delimiter)
        writer.writerow(headers)
        for row in rows:
            out = []
            for h in headers:
                val = row.get(h)
                if val is None:
                    out.append('NA')
                elif isinstance(val, float):
                    out.append(f'{val:.4f}')
                else:
                    out.append(val)
            writer.writerow(out)
    print(f'Written to {filepath}', file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Summarise per-sample read counts from nf-core/eager output.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s /path/to/eager/results
  %(prog)s /path/to/eager/results -o summary.tsv
  %(prog)s /path/to/eager/results --csv -o summary.csv
  %(prog)s /path/to/eager/results --no-table -o summary.tsv

The script auto-detects the MultiQC data directory inside the eager results
folder and parses FastQC, fastp, AdapterRemoval, samtools flagstat, Picard/DeDup,
endorSpy, and DamageProfiler outputs.  Columns for tools that were not run are
automatically omitted.
""")

    parser.add_argument('results_dir',
                        help='Path to the eager results/ directory')
    parser.add_argument('-o', '--output',
                        help='Output file path (TSV by default, CSV with --csv)')
    parser.add_argument('--csv', action='store_true',
                        help='Write output as CSV instead of TSV')
    parser.add_argument('--no-table', action='store_true',
                        help='Suppress the formatted console table')

    args = parser.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    if not os.path.isdir(results_dir):
        print(f'Error: {results_dir} is not a directory', file=sys.stderr)
        sys.exit(1)

    headers, rows = build_summary(results_dir)

    if not args.no_table:
        print(file=sys.stderr)
        print_table(headers, rows)
        print(file=sys.stderr)

    if args.output:
        delim = ',' if args.csv else '\t'
        write_delimited(headers, rows, args.output, delimiter=delim)


if __name__ == '__main__':
    main()
