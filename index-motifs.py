#!/usr/bin/env python3

"""
Build a raw motif index from a FASTA proteome.

This is only the launcher. The real indexing logic is in:

    motif-indexing/src/main.rs

The backend uses:

    --specificity

as one dynamic control for:

    internal motif depth
    split coverage acceptance
    child survival threshold

There is no center-position preference.
There is no min_proteins floor.
There is no amino-acid frequency/enrichment scoring.
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

    # Backward-compatible ignored argument.
    # Old commands can still contain --min-proteins, but it is not passed to Rust.
    p.add_argument("--min-proteins", type=int, default=None, help=argparse.SUPPRESS)

    p.add_argument("--include-trembl", action="store_true")
    p.add_argument("--include-other", action="store_true")
    p.add_argument("--exclude-aa", default="*BJOXZU")
    p.add_argument("--no-expand-duplicates", action="store_true")

    p.add_argument("--threads", type=int)

    p.add_argument(
        "--chunk-size",
        type=int,
        default=16,
        help="Number of endpoint buckets processed before writing. Lower uses less RAM; higher may improve throughput.",
    )

    p.add_argument(
        "--branch-min-hits",
        type=int,
        default=2048,
        help="Performance-only threshold for parallelizing recursive branches. Does not change intended output.",
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
