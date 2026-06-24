#!/usr/bin/env python3

"""
Build a skeleton motif atlas.

Output:
  proteins.tsv
  motifs.tsv
  postings.bin
  build_info.tsv

Specificity logic:
  internal_budget = round((k - 2) * specificity)

  At each next internal fixed position:
    step = next_internal_fixed / internal_budget
    child_min = max(min_proteins, ceil(parent_support * specificity * step))

So --specificity controls both:
  1. how many internal positions a skeleton may fix
  2. how much parent support a new child skeleton needs


Build a skeleton motif atlas from a FASTA proteome.

A skeleton motif is a partially fixed peptide pattern:

    A...G
    A.P.G
    A.P.LG

The endpoints are fixed first. Internal positions are filled only when doing so
creates a pattern that is supported by enough proteins.

The main tuning flag is:

    --specificity 0.35

Specificity controls how detailed motifs are allowed to become.

It dynamically controls two linked things:

1. How many internal positions may be fixed

   For a motif length k:

       internal_positions = k - 2
       internal_budget = round(internal_positions * specificity)

   Endpoints do not count toward this budget.

   Example with --specificity 0.35:

       k=6   -> 4 internal slots  -> about 1 slot may be fixed
       k=12  -> 10 internal slots -> about 4 slots may be fixed
       k=25  -> 23 internal slots -> about 8 slots may be fixed

2. How much evidence is needed to create a more-specific child motif

   As a skeleton becomes more specific, the support required for another child
   also increases.

       progress = next_internal_fixed_position / internal_budget
       required_fraction = specificity * progress
       child_min = max(min_proteins, ceil(parent_support * required_fraction))

   This prevents broad parent motifs from spawning many tiny children while also
   preventing already-specific motifs from fragmenting endlessly.

In plain terms:

    --specificity controls the allowed detail level of the atlas.

Lower values produce broader, more general skeletons.
Higher values allow more detailed skeletons, but require stronger support as
the motif becomes more specific.

The default, --specificity 0.35, means the atlas may fix roughly 35% of the
internal positions, with child-support requirements tightening as that budget
is used.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build skeleton motif index.")

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
