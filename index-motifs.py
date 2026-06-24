#!/usr/bin/env python3

"""
Build a raw motif index from a FASTA proteome.

A motif pattern is a compressed description of related protein windows.

A letter means that position is fixed.
A dot means that position is open.

So:

    A.....G

means:

    any 7-amino-acid window that starts with A and ends with G

The Rust code starts with broad endpoint patterns. Then it asks whether one
internal position should be filled to create more specific child patterns:

    A.....G
        A.P...G
        A.K...G
        A...L.G

The main control is:

    --specificity

Specificity means exactly two things in this version.

First: length-scaled depth.

Instead of a hard-coded max_fixed value, the maximum number of internal fixed
positions is derived from motif length:

    internal_positions = k - 2
    internal_budget = round(internal_positions * specificity)

Endpoints do not count because they define the starting parent pattern.

Second: parent-level split coverage.

A parent should not split just because one tiny child exists. The children from
one proposed split position must collectively explain enough of the parent:

    split_coverage =
        unique proteins covered by all children from that split
        /
        parent proteins

The required split coverage is derived from specificity:

    required_split_coverage = 1 - specificity

So at specificity 0.35, children must cover at least 65% of the parent protein
set for the split to happen.

Lower specificity keeps motifs broader.
Higher specificity allows deeper and smaller substructure.

The splitter avoids an artificial left-side bias. If several split positions
are equally good, it does not simply pick the earliest internal position. It
scores split positions by protein coverage and child structure first, then uses
centrality as a final tie-break.

This script writes only the raw atlas:

    motifs.tsv
    postings.bin
    proteins.tsv
    build_info.tsv

No family hierarchy.
No cleaned motif table.
No cross-length merge.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build a raw skeleton motif index.")

    p.add_argument("fasta", type=Path)
    p.add_argument("outdir", type=Path)

    p.add_argument("--min", dest="min_k", type=int, default=5)
    p.add_argument("--max", dest="max_k", type=int, default=12)
    p.add_argument("--specificity", type=float, default=0.35)

    p.add_argument("--min-proteins", type=int, default=2)
    p.add_argument("--include-trembl", action="store_true")
    p.add_argument("--include-other", action="store_true")
    p.add_argument("--exclude-aa", default="*BJOXZU")
    p.add_argument("--no-expand-duplicates", action="store_true")

    p.add_argument("--threads", type=int)
    p.add_argument("--chunk-size", type=int, default=64)
    p.add_argument(
        "--branch-min-hits",
        type=int,
        default=8192,
        help="Performance-only threshold for parallelizing recursive branches. Does not change output.",
    )
    p.add_argument("--overwrite", action="store_true")

    return p.parse_args()


def main() -> int:
    args = parse_args()

    here = Path(__file__).resolve().parent
    rust_dir = here / "motif-indexing"

    cargo = shutil.which("cargo")
    if cargo is None:
        print("ERROR: cargo not found on PATH", file=sys.stderr)
        return 1

    subprocess.run([cargo, "build", "--release"], cwd=rust_dir, check=True)

    binary = rust_dir / "target" / "release" / "motif_index"
    if sys.platform.startswith("win"):
        binary = binary.with_suffix(".exe")

    cmd = [
        str(binary),
        str(args.fasta),
        str(args.outdir),
        "--min",
        str(args.min_k),
        "--max",
        str(args.max_k),
        "--specificity",
        str(args.specificity),
        "--min-proteins",
        str(args.min_proteins),
        "--exclude-aa",
        args.exclude_aa,
        "--chunk-size",
        str(args.chunk_size),
        "--branch-min-hits",
        str(args.branch_min_hits),
    ]

    if args.include_trembl:
        cmd.append("--include-trembl")
    if args.include_other:
        cmd.append("--include-other")
    if args.no_expand_duplicates:
        cmd.append("--no-expand-duplicates")
    if args.threads is not None:
        cmd.extend(["--threads", str(args.threads)])
    if args.overwrite:
        cmd.append("--overwrite")

    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
