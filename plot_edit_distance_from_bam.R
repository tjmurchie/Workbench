#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_edit_distance_from_bam.R
# -----------------------------------------------------------------------------
# Extract NM edit-distance tags from one or more BAM files and generate:
#   1. a plain-text list of per-read edit distances
#   2. a PNG histogram
#   3. an SVG histogram
#
# The NM tag stores the edit distance between a read and the reference sequence
# for that alignment. Lower values generally indicate closer sequence identity,
# although interpretation always depends on mapping context, read length,
# filtering, and the expected biology of the sample.
#
# This script is intentionally simple and file-based:
# - it does not modify BAMs
# - it does not require indexing
# - it streams `samtools view` output and pulls NM tags directly from records
# - it writes both raster and vector plots for convenience
#
# Typical use cases
# -----------------
# - quick QC checks of mapped BAMs
# - comparing candidate mapping targets
# - generating a lightweight histogram for reports or figure drafting
#
# Requirements
# ------------
# - R
# - ggplot2
# - samtools
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This script requires ggplot2. Install it with install.packages('ggplot2')")
  }
})

library(ggplot2)

print_usage <- function() {
  cat(
"Usage:
  Rscript plot_edit_distance_from_bam.R [--samtools /path/to/samtools] [--outdir DIR] file1.bam [file2.bam ...]

What it does:
  - runs samtools view on each BAM
  - extracts NM:i:<value> edit distance tags
  - writes a text file of edit distances
  - writes histogram plots as PNG and SVG

Arguments:
  --samtools PATH   Optional path to samtools executable.
                    If omitted, the script tries:
                    1) samtools from PATH
                    2) /mnt/Genomics/Working/tyler.murchie/legacy-dependencies/MakefileScripts/samtools-patched/sam

  --outdir DIR      Output directory. Default: EditDistance

Input:
  One or more BAM files with NM:i tags present in the SAM records.

Outputs per BAM:
  <basename>.edit-distance.txt
  <basename>.edit-distance.png
  <basename>.edit-distance.svg

Examples:
  Rscript plot_edit_distance_from_bam.R sample.bam

  Rscript plot_edit_distance_from_bam.R \
    --samtools /mnt/Genomics/Working/tyler.murchie/legacy-dependencies/MakefileScripts/samtools-patched/sam \
    --outdir EditDistance_2 \
    MoccasinMerged_Human-Mapped.min24MQ30.bam

  Rscript plot_edit_distance_from_bam.R \
    --samtools /mnt/Genomics/Working/tyler.murchie/legacy-dependencies/MakefileScripts/samtools-patched/sam \
    --outdir EditDistance_2 \
    *.bam
", sep = "")
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  print_usage()
  quit(save = "no", status = 0)
}

samtools <- NULL
outdir <- "EditDistance"
bam_files <- character(0)

i <- 1
while (i <= length(args)) {
  arg <- args[i]

  if (arg == "--samtools") {
    if (i == length(args)) stop("Missing value after --samtools")
    samtools <- args[i + 1]
    i <- i + 2
  } else if (arg == "--outdir") {
    if (i == length(args)) stop("Missing value after --outdir")
    outdir <- args[i + 1]
    i <- i + 2
  } else if (startsWith(arg, "--")) {
    stop(paste("Unknown option:", arg))
  } else {
    bam_files <- c(bam_files, arg)
    i <- i + 1
  }
}

if (length(bam_files) == 0) {
  stop("No BAM files provided. Run with --help for usage.")
}

# Resolve samtools in a predictable order so the script works both on
# systems with a standard PATH install and in older local environments where
# a patched binary was used.
resolve_samtools <- function(user_path = NULL) {
  candidates <- character(0)

  if (!is.null(user_path)) {
    candidates <- c(candidates, user_path)
  }

  path_samtools <- Sys.which("samtools")
  if (nzchar(path_samtools)) {
    candidates <- c(candidates, unname(path_samtools))
  }

  patched <- "/mnt/Genomics/Working/tyler.murchie/legacy-dependencies/MakefileScripts/samtools-patched/sam"
  candidates <- c(candidates, patched)

  candidates <- unique(candidates)

  for (cand in candidates) {
    if (file.exists(cand) || identical(basename(cand), "samtools")) {
      test <- tryCatch(
        system2(cand, "--version", stdout = TRUE, stderr = TRUE),
        error = function(e) NULL
      )
      if (!is.null(test)) return(cand)
    }
  }

  stop(
    "Could not find a working samtools executable. ",
    "Provide one with --samtools /path/to/samtools"
  )
}

samtools <- resolve_samtools(samtools)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Stream SAM records and pull out NM:i:<int> tags without creating
# intermediate SAM files on disk. This keeps the script simple and avoids
# temporary-file overhead for typical BAM sizes used in quick QC.
extract_nm_values <- function(bam, samtools_path) {
  cmd <- paste(shQuote(samtools_path), "view", shQuote(bam))
  con <- pipe(cmd, open = "r")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  values <- integer(0)

  repeat {
    lines <- readLines(con, n = 100000, warn = FALSE)
    if (length(lines) == 0) break

    hits <- regmatches(lines, regexpr("NM:i:[0-9]+", lines, perl = TRUE))
    keep <- nzchar(hits)

    if (any(keep)) {
      nm <- suppressWarnings(as.integer(sub("NM:i:", "", hits[keep], fixed = TRUE)))
      nm <- nm[!is.na(nm)]
      if (length(nm) > 0) values <- c(values, nm)
    }
  }

  values
}

# Build a one-bin-per-integer histogram so edit distances 0, 1, 2, ... are
# shown as discrete categories rather than smoothed bins.
make_plot <- function(values, sample_name) {
  max_ed <- max(values, na.rm = TRUE)

  ggplot(data.frame(EditDistance = values), aes(x = EditDistance)) +
    geom_histogram(
      binwidth = 1,
      boundary = -0.5,
      color = "black",
      fill = "white"
    ) +
    scale_x_continuous(breaks = 0:max_ed) +
    labs(
      title = paste("Edit distance histogram:", sample_name),
      subtitle = "Low edit distance suggests a better mapping target",
      x = "Edit distance (NM tag)",
      y = "Read count"
    ) +
    theme_bw(base_size = 12)
}

for (bam in bam_files) {
  if (!file.exists(bam)) {
    warning(sprintf("Skipping missing file: %s", bam))
    next
  }

  sample_name <- sub("\\.[Bb][Aa][Mm]$", "", basename(bam))
  message("Processing: ", bam)

  nm_values <- extract_nm_values(bam, samtools)

  if (length(nm_values) == 0) {
    warning(sprintf("No NM tags found in %s. No output written.", bam))
    next
  }

  txt_file <- file.path(outdir, paste0(sample_name, ".edit-distance.txt"))
  png_file <- file.path(outdir, paste0(sample_name, ".edit-distance.png"))
  svg_file <- file.path(outdir, paste0(sample_name, ".edit-distance.svg"))

  write.table(
    nm_values,
    file = txt_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  p <- make_plot(nm_values, sample_name)

  grDevices::png(filename = png_file, width = 1800, height = 1200, res = 200)
  print(p)
  grDevices::dev.off()

  grDevices::svg(filename = svg_file, width = 10, height = 7)
  print(p)
  grDevices::dev.off()

  message("  Wrote: ", txt_file)
  message("  Wrote: ", png_file)
  message("  Wrote: ", svg_file)
}
