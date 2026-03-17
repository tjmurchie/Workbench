#!/usr/bin/env python3
"""
MEGAN-rma6_summary.py — Extract metadata and counts from MEGAN .rma6 files
=====================================================================

Batch-extracts LCA settings, read/match counts, and major-node taxonomy
summaries from .rma6 files produced by MEGAN's blast2rma or daa2rma tools.
Outputs a single TSV with one row per file, making it easy to verify that
all files in a project share identical LCA parameters and to get at-a-glance
count summaries without opening the MEGAN GUI.

Data sources
------------
Each .rma6 file is queried three ways:

1. **Direct binary read** — The embedded metadata block (starting at
   ``@Creator``) is extracted from the raw bytes. This is the only reliable
   way to recover the full ``@Parameters`` string (``rma2info -m`` leaves it
   blank) and other header fields like ``@Creator``, ``@CreationDate``,
   ``@BlastMode``, and ``@Contaminants``. The metadata is stored as Java
   modified UTF-8 inside the .rma6 container, so the standard ``strings``
   utility does not find it.

2. **``rma2info -l``** — Fast listing that provides ``Number of reads``
   and ``Number of matches``.

3. **``rma2info -c2c Taxonomy``** — Two calls are made:
   - Summarized (``-s``) with names (``-n``): provides counts for major
     taxonomy nodes (Bacteria, Eukaryota, etc.) where each node's count
     includes all descendant reads.
   - Non-summarized (no ``-s``), IDs only: provides the special node
     counts for Not Assigned (taxid -2), No Hits (taxid -1), and
     Contaminants (taxid -6, present only in ContamRem/ConRemv files).
     This call is only made for files that have ``@Contaminants`` in
     their binary metadata, to avoid unnecessary overhead.

Output columns (one row per .rma6 file)
---------------------------------------
**Identification:**
    File, FilePath

**Metadata:**
    Creator, CreationDate, BlastMode

**Counts:**
    NumReads, NumMatches

**LCA settings:**
    minScore, maxExpected, minPercentIdentity, topPercent,
    minSupportPercent, minSupport, lcaAlgorithm, minPercentReadToCover,
    minPercentReferenceToCover, minReadLength, longReads, pairedReads,
    identityFilter, contaminantFilter, readAssignmentMode

**Taxonomy summary:**
    NotAssigned, NoHits, Contaminants, Assigned, cellular_organisms,
    Bacteria, Archaea, Eukaryota, Viruses, Viroids, Viridiplantae,
    Metazoa, Fungi

The ``Assigned`` column matches MEGAN's status-bar "Assigned" count:
``Assigned = NumReads - NotAssigned - NoHits - Contaminants``.
``Contaminants`` is the number of reads placed on disabled (contaminant)
taxa (MEGAN internal taxid -6); it is 0 for non-ContamRem files and
reflects the ``@Contaminants`` taxid list for ContamRem/ConRemv files.

Dependencies
------------
- Python 3.6+ with ``pandas``
- ``rma2info`` on PATH (ships with MEGAN, typically at
  ``~/software/megan7/tools/``)

Usage
-----
::

    # All files in current directory
    rma6_summary.py *.rma6 -o rma6_summary.tsv

    # Files from multiple locations
    rma6_summary.py /path/A/*.rma6 /path/B/*.rma6 -o combined.tsv

    # Verbose progress
    rma6_summary.py *.rma6 -o summary.tsv -v

Author
------
Tyler Murchie, with assistance from Claude (Anthropic).
"""

import argparse
import os
import re
import subprocess
import sys
from collections import OrderedDict

import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# LCA setting keys, in the column order for the output TSV.
LCA_KEYS = [
    "minScore",
    "maxExpected",
    "minPercentIdentity",
    "topPercent",
    "minSupportPercent",
    "minSupport",
    "lcaAlgorithm",
    "minPercentReadToCover",
    "minPercentReferenceToCover",
    "minReadLength",
    "longReads",
    "pairedReads",
    "identityFilter",
    "contaminantFilter",
    "readAssignmentMode",
]

# Major taxonomy nodes to report in the output.  Display names use
# underscores; MEGAN names use spaces (e.g. "cellular organisms").
TAXONOMY_NODES = [
    "cellular_organisms",
    "Bacteria",
    "Archaea",
    "Eukaryota",
    "Viruses",
    "Viroids",
    "Viridiplantae",
    "Metazoa",
    "Fungi",
]

# Metadata keys we care about from the binary @Key\tValue header block.
METADATA_KEYS = {
    "Creator",
    "CreationDate",
    "BlastMode",
    "Parameters",
    "Contaminants",
}

# Lines emitted by rma2info that are not data (version banner, etc.).
_RMA2INFO_NOISE_PREFIXES = (
    "Version",
    "Author",
    "Copyright",
    "Java",
    "Loading",
    "Total time",
    "Peak memory",
)


# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------

def parse_args():
    """Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        ``rma6_files`` : list of str — input .rma6 file paths.
        ``output``     : str — output TSV path.
        ``verbose``    : bool — print per-file progress to stderr.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract metadata and counts from MEGAN .rma6 files into a TSV.  "
            "Produces one row per file with LCA settings, read/match counts, "
            "and summarised taxonomy counts for major nodes."
        ),
    )
    parser.add_argument(
        "rma6_files",
        nargs="+",
        help="One or more .rma6 files (shell globs are expanded by the shell)",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file path",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print per-file progress to stderr",
    )
    return parser.parse_args()


def extract_binary_metadata(path):
    """Read the .rma6 binary to extract the ``@Key\\tValue`` metadata block.

    MEGAN stores header metadata as Java modified UTF-8 with length prefixes
    inside the .rma6 container.  The standard ``strings`` utility cannot
    reliably extract it, so we read the raw bytes and locate the block
    starting from ``@Creator``.

    Parameters
    ----------
    path : str
        Path to the .rma6 file.

    Returns
    -------
    dict
        ``{key: value}`` for each key in :data:`METADATA_KEYS` found in the
        header.  Missing keys are simply absent from the dict.
    """
    result = {}
    try:
        with open(path, "rb") as f:
            data = f.read()

        # The metadata block starts with @Creator and runs for a few hundred
        # bytes.  We grab a generous 4 KiB window to be safe.
        idx = data.find(b"@Creator")
        if idx < 0:
            return result

        chunk = data[idx : idx + 4096].decode("ascii", errors="replace")

        for line in chunk.split("\n"):
            line = line.strip()
            if not line.startswith("@"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            key = parts[0][1:]          # strip leading '@'
            value = parts[1].strip()
            if key in METADATA_KEYS:
                result[key] = value
    except Exception as e:
        print(
            f"  WARNING: binary metadata extraction failed for {path}: {e}",
            file=sys.stderr,
        )
    return result


def parse_parameters(param_str):
    """Parse a MEGAN ``@Parameters`` value string into individual settings.

    The parameter string looks like::

        minScore=50.0 maxExpected='1.0E-7' minPercentIdentity='96.0' ...
        fNames= { Taxonomy }

    Three value formats are handled:
    - **Single-quoted**: ``maxExpected='1.0E-7'`` → ``1.0E-7``
    - **Brace-enclosed**: ``fNames={ Taxonomy }`` → ``Taxonomy``
    - **Plain token**: ``minScore=50.0`` → ``50.0``

    Parameters
    ----------
    param_str : str
        The raw ``@Parameters`` value.

    Returns
    -------
    dict
        ``{setting_name: value_string}`` for every ``key=value`` pair found.
    """
    result = {}
    if not param_str:
        return result
    for m in re.finditer(r"(\w+)=('([^']*)'|\{([^}]*)\}|(\S+))", param_str):
        key = m.group(1)
        # Pick the matching capture group: quoted → braced → plain.
        if m.group(3) is not None:
            value = m.group(3)
        elif m.group(4) is not None:
            value = m.group(4).strip()
        else:
            value = m.group(5)
        result[key] = value
    return result


def _parse_rma2info_taxon_lines(output):
    """Parse tab-separated ``name\\tcount`` lines from rma2info output.

    Filters out the version/loading banner lines that rma2info prints to
    stdout/stderr.

    Parameters
    ----------
    output : str
        Combined stdout + stderr from an rma2info invocation.

    Returns
    -------
    dict
        ``{taxon_name_or_id: int_count}``
    """
    counts = {}
    for line in output.splitlines():
        line = line.strip()
        if not line or line.startswith(_RMA2INFO_NOISE_PREFIXES):
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            continue
        try:
            counts[parts[0].strip()] = int(float(parts[1]))
        except ValueError:
            continue
    return counts


def extract_rma2info_stats(path):
    """Run ``rma2info -l`` and return read/match counts.

    Parameters
    ----------
    path : str
        Path to the .rma6 file.

    Returns
    -------
    dict
        Keys ``NumReads`` and ``NumMatches`` (int), when available.
    """
    result = {}
    try:
        proc = subprocess.run(
            ["rma2info", "-i", path, "-l"],
            capture_output=True, text=True, timeout=120,
        )
        output = proc.stdout + proc.stderr
        for line in output.splitlines():
            if "Number of reads:" in line:
                val = line.split(":")[-1].strip().replace(",", "")
                result["NumReads"] = int(val)
            elif "Number of matches:" in line:
                val = line.split(":")[-1].strip().replace(",", "")
                result["NumMatches"] = int(val)
    except Exception as e:
        print(
            f"  WARNING: rma2info -l failed for {path}: {e}",
            file=sys.stderr,
        )
    return result


def extract_taxonomy_summarized(path):
    """Run ``rma2info -c2c Taxonomy -n -s -u false`` (summarised, named).

    Returns counts for the special nodes (Not Assigned, No Hits) and for
    each major taxonomy node listed in :data:`TAXONOMY_NODES`.

    Parameters
    ----------
    path : str
        Path to the .rma6 file.

    Returns
    -------
    dict
        Keys: ``NotAssigned``, ``NoHits``, and each entry from
        :data:`TAXONOMY_NODES` (int, defaulting to 0 when absent).
    """
    result = {}
    try:
        proc = subprocess.run(
            ["rma2info", "-i", path, "-c2c", "Taxonomy",
             "-n", "-s", "-u", "false"],
            capture_output=True, text=True, timeout=300,
        )
        taxon_counts = _parse_rma2info_taxon_lines(proc.stdout + proc.stderr)

        # MEGAN special taxids:  -1 = No Hits,  -2 = Not Assigned
        result["NotAssigned"] = taxon_counts.get("-2", 0)
        result["NoHits"] = taxon_counts.get("-1", 0)

        for node in TAXONOMY_NODES:
            # Column names use underscores; MEGAN uses spaces.
            lookup = node.replace("_", " ")
            result[node] = taxon_counts.get(lookup, taxon_counts.get(node, 0))

    except Exception as e:
        print(
            f"  WARNING: rma2info taxonomy (summarised) failed for {path}: {e}",
            file=sys.stderr,
        )
    return result


def extract_contaminant_count(path):
    """Run ``rma2info -c2c Taxonomy -u false`` (non-summarised, IDs only).

    Looks specifically for taxid **-6**, which is MEGAN's internal node for
    reads assigned to disabled (contaminant) taxa.  This node only appears
    in ContamRem / ConRemv files and is only visible in non-summarised mode.

    Parameters
    ----------
    path : str
        Path to the .rma6 file.

    Returns
    -------
    int
        Number of reads assigned to contaminant taxa, or 0 if the node is
        absent.
    """
    try:
        proc = subprocess.run(
            ["rma2info", "-i", path, "-c2c", "Taxonomy", "-u", "false"],
            capture_output=True, text=True, timeout=300,
        )
        taxon_counts = _parse_rma2info_taxon_lines(proc.stdout + proc.stderr)
        return taxon_counts.get("-6", 0)
    except Exception as e:
        print(
            f"  WARNING: rma2info contaminant count failed for {path}: {e}",
            file=sys.stderr,
        )
        return 0


# ---------------------------------------------------------------------------
# Per-file orchestration
# ---------------------------------------------------------------------------

def process_one_file(path, verbose=False):
    """Extract all data for one .rma6 file.

    Runs three (or four, for ContamRem files) fast extractions and merges
    the results into a single ordered row.  Each extraction is independently
    try/excepted so that a failure in one section (e.g. rma2info timeout)
    fills that section with ``NA`` while still capturing data from the
    other sources.

    Parameters
    ----------
    path : str
        Path to the .rma6 file.
    verbose : bool
        If True, print the filename to stderr before processing.

    Returns
    -------
    OrderedDict
        One row of data, keyed by column name.
    """
    row = OrderedDict()
    basename = os.path.basename(path)
    abspath = os.path.abspath(path)

    row["File"] = basename
    row["FilePath"] = abspath

    if verbose:
        print(f"  Processing: {basename}", file=sys.stderr)

    # --- 1. Binary metadata (Creator, CreationDate, BlastMode, Parameters) --
    meta = extract_binary_metadata(path)
    row["Creator"] = meta.get("Creator", "NA")
    row["CreationDate"] = meta.get("CreationDate", "NA")
    row["BlastMode"] = meta.get("BlastMode", "NA")

    # --- 2. rma2info -l  (read / match counts) ---
    stats = extract_rma2info_stats(path)
    row["NumReads"] = stats.get("NumReads", "NA")
    row["NumMatches"] = stats.get("NumMatches", "NA")

    # --- 3. LCA parameters from the @Parameters string ---
    params = parse_parameters(meta.get("Parameters", ""))
    for key in LCA_KEYS:
        row[key] = params.get(key, "NA")

    # --- 4. Taxonomy counts (summarised) ---
    tax = extract_taxonomy_summarized(path)
    row["NotAssigned"] = tax.get("NotAssigned", "NA")
    row["NoHits"] = tax.get("NoHits", "NA")

    # --- 5. Contaminant count (non-summarised, only for ContamRem files) ---
    #
    # MEGAN stores reads assigned to disabled (contaminant) taxa under the
    # internal taxid -6.  This node only exists in ContamRem/ConRemv files
    # and is only visible in non-summarised mode.  We detect ContamRem files
    # by the presence of @Contaminants in the binary metadata and make the
    # extra rma2info call only when needed.
    has_contaminants = "Contaminants" in meta
    if has_contaminants:
        contam_count = extract_contaminant_count(path)
    else:
        contam_count = 0
    row["Contaminants"] = contam_count

    # Assigned = NumReads - NotAssigned - NoHits - Contaminants
    # This matches MEGAN's status-bar "Assigned" count.
    try:
        row["Assigned"] = (
            int(row["NumReads"])
            - int(row["NotAssigned"])
            - int(row["NoHits"])
            - int(row["Contaminants"])
        )
    except (ValueError, TypeError):
        row["Assigned"] = "NA"

    # Major taxonomy nodes (summarised counts).
    for node in TAXONOMY_NODES:
        row[node] = tax.get(node, 0)

    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Entry point: parse args, process files, write TSV."""
    args = parse_args()

    # Validate that all input paths exist.
    valid_files = []
    for path in args.rma6_files:
        if not os.path.isfile(path):
            print(f"WARNING: File not found, skipping: {path}", file=sys.stderr)
            continue
        valid_files.append(path)

    if not valid_files:
        print("ERROR: No valid .rma6 files provided.", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(valid_files)} .rma6 file(s)...", file=sys.stderr)

    rows = []
    for path in sorted(valid_files):
        row = process_one_file(path, verbose=args.verbose)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
