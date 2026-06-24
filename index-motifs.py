#!/usr/bin/env python3

"""
Skeleton-first motif atlas wrapper
==================================

This Python file is the user-facing entry point for the motif atlas builder.

It does not perform motif enumeration itself. It finds, builds, and runs the Rust
program in:

    motif-indexing/src/main.rs

The Rust program performs the actual skeleton-first motif indexing.

Why this design exists
----------------------

Earlier exact-peptide-first designs tried to build tables like:

    exact motif -> protein
    scaffold motif -> protein
    scaffold motif -> exact motif
    occurrence rows
    DuckDB aggregation tables

That made the build extremely large because every exact peptide/protein
relationship was materialized. Even with Rust doing the enumeration, the output
became dominated by disk I/O and database aggregation.

The current design changes the game:

    The primary object is not an exact peptide.
    The primary object is a skeleton motif.

For example:

    A...G
    A.P.G
    A.PQG

Each skeleton motif stores a compact protein posting list. Exact peptide examples
and exact peptide support rows are deliberately not stored.

High-level motif collection process
-----------------------------------

1. Load the FASTA.

   By default, only SwissProt entries are used. These are recognized by FASTA
   accessions beginning with:

       sp|

   TrEMBL entries can be included with:

       --include-trembl

   Other accessions can be included with:

       --include-other

2. Deduplicate identical protein sequences.

   If multiple FASTA records have the exact same amino-acid sequence, the Rust
   builder enumerates the sequence once.

   By default, the resulting motif hits are assigned back to all duplicate
   accessions. This keeps protein identity while avoiding repeated work.

   To index only representative sequences, use:

       --no-expand-duplicates

3. Process peptide lengths.

   The default range is:

       --min-k 5
       --max-k 12

   Each k is processed as its own stage. This gives clear progress messages:

       [k=5] collecting endpoint skeleton groups...
       [k=5] refining skeletons...
       [k=6] collecting endpoint skeleton groups...
       ...

4. Collect endpoint skeleton groups.

   For each valid k-length window, the first and last amino acids define the
   initial broad skeleton group.

   Example for k=5:

       A---G becomes A...G

   This produces at most 20 x 20 endpoint groups for canonical amino acids.

5. Refine skeletons inside each endpoint group.

   Inside an endpoint group, the Rust builder recursively tests internal fixed
   positions.

   Example:

       A...G
       A.P.G
       A.PQG

   Internal positions are added only when they create useful protein-set
   structure. This is the controlled replacement for the older exact-peptide
   explosion.

6. Store only motifs that appear in enough proteins.

   The default is:

       --min-proteins 2

   This means singleton motifs are not kept. This is not a broad frequency
   filter; it is only a minimum visibility rule.

7. Write compact motif postings.

   The final index stores:

       motif_text -> compressed list of protein_ids

   The compressed protein lists are written into:

       postings.bin

   The offset and byte length for each motif's posting list are stored in:

       motifs.tsv

Output files
------------

The output directory contains:

    proteins.tsv
        protein_id
        accession
        source
        representative_protein_id
        is_representative

    motifs.tsv
        motif_id
        motif_text
        posting_offset
        posting_bytes

    postings.bin
        Concatenated delta-varint encoded sorted protein_id lists.

    build_info.tsv
        Build settings and high-level counts.

Not stored
----------

The current design deliberately does not store:

    exact peptide examples
    exact peptide support rows
    occurrence counts
    motif length
    protein count
    total count
    protein -> motif reverse index
    DuckDB database
    Parquet exports
    raw build TSV event logs

Reasons:

    motif length is inherent in motif_text

    protein count is recoverable by decoding the posting list

    total occurrence count is not needed for presence-based motif grouping

    exact peptide examples are debugging material, not part of the atlas

    protein -> motif can be derived later if needed

    raw event logs recreate the disk-space problem this design avoids

Default usage
-------------

    python ~/motif-quantification/index_motifs.py \\
      ~/data/proteomics/fastas/proteomes/Human_Homo_sapien.fasta \\
      human-skeleton-motifs \\
      --threads 8 \\
      --overwrite

Useful knobs
------------

Make a smaller motif atlas:

    --min-proteins 5

Allow more specific skeletons:

    --max-fixed 5

Allow broader skeletons only:

    --max-fixed 3

Include TrEMBL:

    --include-trembl

Index representative sequences only:

    --no-expand-duplicates

Current conceptual model
------------------------

This builder is intended to create a compact reusable skeleton motif atlas.

It favors:

    coverage through broad skeletons
    specificity through recursive internal-position refinement
    storage efficiency through compressed protein postings

It avoids:

    exhaustive exact peptide indexing
    massive motif/protein edge tables
    post-hoc database aggregation as the core workflow
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
RUST_PROJECT_DIR = SCRIPT_DIR / "motif-indexing"
RUST_MANIFEST = RUST_PROJECT_DIR / "Cargo.toml"
RUST_BINARY_RELEASE = RUST_PROJECT_DIR / "target" / "release" / "motif_index"
RUST_BINARY_DEBUG = RUST_PROJECT_DIR / "target" / "debug" / "motif_index"


def fail(message: str, exit_code: int = 1) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(exit_code)


def run_command(command: list[str], cwd: Path | None = None) -> None:
    print("+ " + " ".join(command), file=sys.stderr)

    completed = subprocess.run(command, cwd=cwd)

    if completed.returncode != 0:
        fail(
            f"command failed with exit code {completed.returncode}: {' '.join(command)}",
            completed.returncode,
        )


def newest_mtime(paths: list[Path]) -> float:
    mtimes = []

    for path in paths:
        if path.is_file():
            mtimes.append(path.stat().st_mtime)
        elif path.is_dir():
            for child in path.rglob("*"):
                if child.is_file():
                    mtimes.append(child.stat().st_mtime)

    return max(mtimes) if mtimes else 0.0


def rust_binary_is_stale(binary: Path) -> bool:
    if not binary.exists():
        return True

    source_mtime = newest_mtime(
        [
            RUST_PROJECT_DIR / "src",
            RUST_PROJECT_DIR / "Cargo.toml",
            RUST_PROJECT_DIR / "Cargo.lock",
        ]
    )

    return source_mtime > binary.stat().st_mtime


def build_rust_binary(debug: bool, clean_build: bool, no_build: bool) -> Path:
    if not RUST_MANIFEST.exists():
        fail(f"Rust project not found: {RUST_MANIFEST}")

    binary = RUST_BINARY_DEBUG if debug else RUST_BINARY_RELEASE

    if no_build:
        if not binary.exists():
            fail(f"Rust binary does not exist and --no-build was used: {binary}")
        return binary

    if clean_build:
        run_command(["cargo", "clean", "--manifest-path", str(RUST_MANIFEST)])

    if rust_binary_is_stale(binary):
        if shutil.which("cargo") is None:
            fail("cargo was not found on PATH")

        command = ["cargo", "build", "--manifest-path", str(RUST_MANIFEST)]

        if not debug:
            command.append("--release")

        run_command(command)

    if not binary.exists():
        fail(f"Rust build finished, but expected binary was not found: {binary}")

    return binary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compact skeleton-first motif atlas.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("fasta", help="Input protein FASTA.")
    parser.add_argument("outdir", help="Output folder.")

    parser.add_argument(
        "--min",
        dest="min_k",
        type=int,
        default=5,
        help="Minimum peptide/window length.",
    )

    parser.add_argument(
        "--max",
        dest="max_k",
        type=int,
        default=12,
        help="Maximum peptide/window length.",
    )

    parser.add_argument(
        "--include-trembl",
        action="store_true",
        help="Include tr| TrEMBL entries. Default is SwissProt-only.",
    )

    parser.add_argument(
        "--include-other",
        action="store_true",
        help="Include FASTA entries that are neither sp| nor tr|.",
    )

    parser.add_argument(
        "--exclude-aa",
        default="*BJOXZU",
        help="Amino acid symbols that invalidate a window.",
    )

    parser.add_argument(
        "--min-proteins",
        type=int,
        default=2,
        help="Minimum number of proteins required to keep a motif.",
    )

    parser.add_argument(
        "--max-fixed",
        type=int,
        default=4,
        help="Maximum fixed positions in a skeleton, including both endpoints.",
    )

    parser.add_argument(
        "--no-expand-duplicates",
        action="store_true",
        help="Only index representative sequences, not duplicate accessions.",
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=64,
        help="Representative protein groups per Rayon task.",
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Rust/Rayon worker threads.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output folder if it already exists.",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Build/run the debug Rust binary instead of release.",
    )

    parser.add_argument(
        "--clean-build",
        action="store_true",
        help="Run cargo clean before building.",
    )

    parser.add_argument(
        "--no-build",
        action="store_true",
        help="Do not build; use the existing Rust binary.",
    )

    parser.add_argument(
        "--print-command",
        action="store_true",
        help="Print the final Rust command and exit without running it.",
    )

    return parser.parse_args()


def build_rust_args(args: argparse.Namespace) -> list[str]:
    rust_args = [
        str(Path(args.fasta).expanduser()),
        str(Path(args.outdir).expanduser()),
        "--min-k",
        str(args.min_k),
        "--max-k",
        str(args.max_k),
        "--exclude-aa",
        args.exclude_aa,
        "--min-proteins",
        str(args.min_proteins),
        "--max-fixed",
        str(args.max_fixed),
        "--chunk-size",
        str(args.chunk_size),
    ]

    if args.include_trembl:
        rust_args.append("--include-trembl")

    if args.include_other:
        rust_args.append("--include-other")

    if args.no_expand_duplicates:
        rust_args.append("--no-expand-duplicates")

    if args.threads is not None:
        rust_args.extend(["--threads", str(args.threads)])

    if args.overwrite:
        rust_args.append("--overwrite")

    return rust_args


def main() -> None:
    args = parse_args()

    binary = build_rust_binary(
        debug=args.debug,
        clean_build=args.clean_build,
        no_build=args.no_build,
    )

    command = [str(binary), *build_rust_args(args)]

    if args.print_command:
        print(" ".join(command))
        return

    run_command(command)


if __name__ == "__main__":
    main()
