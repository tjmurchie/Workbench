#!/usr/bin/env python3
"""
recreate_mapDamage_5p3p_freqs.py
================================
Recreate legacy mapDamage output files from a modern ``misincorporation.txt``
file.

This script writes:
  - ``5pCtoT_freq.txt``  (position, 5' C>T frequency)
  - ``3pGtoA_freq.txt``  (position, 3' G>A frequency)

Why this exists
---------------
Older versions of mapDamage produced these two summary tables directly. Newer
versions provide the richer ``misincorporation.txt`` table instead. Some legacy
analysis or plotting workflows still expect the older files, so this script
reconstructs them in a way that matches the legacy calculation.

Typical usage
-------------
::

    python3 recreate_mapDamage_5p3p_freqs.py misincorporation.txt --max-pos 25

Computation (matches legacy behaviour)
--------------------------------------
  5pC>T(pos) = sum_{rows End=5p, Pos=pos} (C>T) / sum_{rows End=5p, Pos=pos} (C)
  3pG>A(pos) = sum_{rows End=3p, Pos=pos} (G>A) / sum_{rows End=3p, Pos=pos} (G)

Compatibility
-------------
Works with both the older ``misincorporation.txt`` format
(``Chr/End/Std/Pos/...``) and newer versions that include additional columns
such as ``Sample`` and ``Library``.
"""

from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Tuple, Iterable

def read_rows(misincorp_path: Path) -> Iterable[dict]:
    """Yield dict rows from misincorporation.txt, skipping comment lines."""
    with misincorp_path.open("r", encoding="utf-8", errors="replace") as fh:
        # Skip leading comment lines
        header = None
        # Find header line
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            header = line.split("\t")
            break
        if header is None:
            raise SystemExit(f"ERROR: no header/data found in {misincorp_path}")

        reader = csv.DictReader(fh, fieldnames=header, delimiter="\t")
        for row in reader:
            # DictReader will include the raw strings; skip empty rows
            if not row or all(v is None or v == "" for v in row.values()):
                continue
            yield row

def as_int(x: str) -> int:
    return int(float(x))

def as_float(x: str) -> float:
    return float(x)

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Recreate legacy 5pCtoT_freq.txt and 3pGtoA_freq.txt from misincorporation.txt"
    )
    ap.add_argument("misincorporation", type=Path, help="Path to misincorporation.txt")
    ap.add_argument("-o", "--outdir", type=Path, default=None,
                    help="Output directory (default: misincorporation.txt parent directory)")
    ap.add_argument("--max-pos", type=int, default=25,
                    help="Maximum position to write (default: 25; legacy default)")
    ap.add_argument("--all-positions", action="store_true",
                    help="Write all positions present in the file (ignores --max-pos)")
    args = ap.parse_args()

    mis = args.misincorporation
    if not mis.exists():
        raise SystemExit(f"ERROR: file not found: {mis}")

    outdir = args.outdir if args.outdir is not None else mis.parent
    outdir.mkdir(parents=True, exist_ok=True)

    # Accumulators: (End, Pos) -> sums
    sum_C: Dict[Tuple[str,int], float] = {}
    sum_CtoT: Dict[Tuple[str,int], float] = {}
    sum_G: Dict[Tuple[str,int], float] = {}
    sum_GtoA: Dict[Tuple[str,int], float] = {}

    positions_seen = set()

    for row in read_rows(mis):
        end = row.get("End")
        pos_s = row.get("Pos")
        if end is None or pos_s is None:
            continue

        pos = as_int(pos_s)
        positions_seen.add(pos)

        # Columns are named exactly like "C>T" and "G>A" in mapDamage output.
        # Some mapDamage builds may output as "C>T" etc; if missing, treat as 0.
        C = as_float(row.get("C", "0") or "0")
        Ct = as_float(row.get("C>T", "0") or "0")
        G = as_float(row.get("G", "0") or "0")
        Ga = as_float(row.get("G>A", "0") or "0")

        key = (end, pos)
        sum_C[key] = sum_C.get(key, 0.0) + C
        sum_CtoT[key] = sum_CtoT.get(key, 0.0) + Ct
        sum_G[key] = sum_G.get(key, 0.0) + G
        sum_GtoA[key] = sum_GtoA.get(key, 0.0) + Ga

    if args.all_positions:
        pos_list = sorted(positions_seen)
    else:
        pos_list = list(range(1, args.max_pos + 1))

    # Write 5pCtoT_freq.txt
    p1 = outdir / "5pCtoT_freq.txt"
    with p1.open("w", encoding="utf-8") as fh:
        fh.write("pos\t5pC>T\n")
        for pos in pos_list:
            key = ("5p", pos)
            den = sum_C.get(key, 0.0)
            num = sum_CtoT.get(key, 0.0)
            freq = (num / den) if den > 0 else float("nan")
            fh.write(f"{pos}\t{freq:.15g}\n")

    # Write 3pGtoA_freq.txt
    p2 = outdir / "3pGtoA_freq.txt"
    with p2.open("w", encoding="utf-8") as fh:
        fh.write("pos\t3pG>A\n")
        for pos in pos_list:
            key = ("3p", pos)
            den = sum_G.get(key, 0.0)
            num = sum_GtoA.get(key, 0.0)
            freq = (num / den) if den > 0 else float("nan")
            fh.write(f"{pos}\t{freq:.15g}\n")

    print(f"Wrote: {p1}")
    print(f"Wrote: {p2}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
