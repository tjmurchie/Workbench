<p align="center">
  <img src="assets/Workbench.png" alt="Workbench header" width="700">
</p>

# Workbench
Workbench is a collection of lightweight utility scripts for bioinformatics, plotting, and file summarization that don't warrant standalone repositories.

The tools here are intentionally practical and task-oriented: each script is meant to solve a common analysis or file-handling problem quickly, without requiring a large package or framework.

## Quick start

Examples:

```bash
# eager summary table
python3 eager_read_summary.py /path/to/eager/results -o eager_summary.tsv

# MEGAN .rma6 summary
python3 MEGAN-rma6_summary.py *.rma6 -o rma6_summary.tsv

# edit-distance histograms
Rscript plot_edit_distance_from_bam.R --outdir EditDistance *.bam

# mapping summaries against all references in each BAM
bash bam_ref_stats.sh -i '*.bam' > bam_ref_stats.tsv

# mapping summaries against selected references
bash bam_ref_stats.sh -i '*.bam' --ref NC_012920.1 --ref chrM > mt_stats.tsv

# split BAMs by reference
THREADS=8 OUTROOT=split_by_reference bash split_by_ref.sh

# recreate legacy mapDamage terminal-damage files
python3 recreate_mapDamage_5p3p_freqs.py misincorporation.txt --max-pos 25
```

---

## Included tools

### `eager_read_summary.py`
Summarises read-count progression across an `nf-core/eager` run by parsing the final MultiQC data files.

Useful for:
- generating manuscript/SOM-ready read-count tables
- checking how many reads survive each major pipeline step
- comparing sample performance across a run

Typical outputs include columns such as raw reads, post-clipping reads, mapped reads, duplication, endogenous content, and read length metrics.

### `MEGAN-rma6_summary.py`
Extracts metadata, LCA settings, read counts, match counts, and broad taxonomy counts from one or more MEGAN `.rma6` files.

Useful for:
- confirming that all `.rma6` files were generated with the same MEGAN/LCA settings
- capturing project-wide summaries without opening the MEGAN GUI
- logging assignment statistics for QC and reporting

### `plot_edit_distance_from_bam.R`
Reads one or more BAM files, extracts `NM:i` edit-distance tags, writes the raw edit distances to text, and plots edit-distance histograms in both PNG and SVG format.

Useful for:
- quickly visualising edit-distance distributions from mapped BAMs
- comparing mapping quality profiles across files
- generating publication-ready vector output alongside raster output

### `bam_ref_stats.sh`
Computes per-reference or whole-BAM mapping summaries using `samtools`, including mapped-read counts above a MAPQ threshold, coverage breadth, mean depth, and covered bases.

Useful for:
- quickly comparing how strongly samples map to one or more references
- summarising candidate references in competitive-mapping workflows
- reporting breadth/depth-style mapping metrics in a simple tabular format

### `split_by_ref.sh`
Splits BAM files into per-reference BAMs, writing one directory per reference and one BAM per sample within each directory.

Useful for:
- separating multi-reference BAMs into per-contig/per-target BAMs
- organising downstream analyses by reference sequence
- preparing input files for manual inspection or per-reference plotting

### `recreate_mapDamage_5p3p_freqs.py`
Recreates legacy `mapDamage` output files (`5pCtoT_freq.txt` and `3pGtoA_freq.txt`) from newer `misincorporation.txt` files.

Useful for:
- restoring output tables expected by older downstream workflows
- comparing modern mapDamage output with legacy projects
- generating simple terminal-damage frequency files for plotting or archiving

---

## Requirements

Workbench is not a single packaged software environment. Each script has its own lightweight requirements.

### General
- Linux/macOS shell environment recommended
- Python 3 for Python scripts
- R for the R plotting script

### External tools used by specific scripts
- `samtools` for BAM-processing scripts
- `rma2info` for `MEGAN-rma6_summary.py`
- `ggplot2` for `plot_edit_distance_from_bam.R`
- `pandas` for `MEGAN-rma6_summary.py`

Refer to the header documentation inside each script for exact usage and assumptions.

---

## Suggested repository layout

This repository is intentionally flat so the scripts are easy to browse and run directly:

```text
Workbench/
├── README.md
├── .gitignore
├── LICENSE
├── assets/
│   ├── Workbench.png
│   ├── Workbench_editable.svg
│   └── Workbench_editable_preview.png
├── eager_read_summary.py
├── MEGAN-rma6_summary.py
├── plot_edit_distance_from_bam.R
├── bam_ref_stats.sh
├── split_by_ref.sh
└── recreate_mapDamage_5p3p_freqs.py
```

---

## Usage philosophy

These scripts are meant to be:
- easy to inspect
- easy to edit for one-off project needs
- usable directly from the command line
- understandable months later with clear in-file documentation

They are not currently packaged as a formal Python/R library, and that is intentional.

---

## Notes


- `assets/Workbench.png` is the current GitHub header banner.
- `assets/Workbench_editable.svg` is the editable vector version from the earlier custom pixel-art pass.
- `assets/Workbench_editable_preview.png` is the PNG preview that matches that editable SVG.
- Most scripts are designed for local command-line use rather than installation through `pip` or `conda`.
- Output files and intermediate products are usually better kept out of version control; a starter `.gitignore` is included for that purpose.
- If a script grows substantially in scope, it can later be promoted into its own dedicated repository.

---

## Author

Tyler Murchie


## Licensing

Workbench is released under the MIT License.

See [`LICENSE`](LICENSE) for the full text.
