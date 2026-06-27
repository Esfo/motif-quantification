#!/usr/bin/env python3
"""Run MS1 distribution determination over a project's centroid mzML files.

Wraps distributions/index_ms1.py (the per-file detector that writes the sqlite
the GUI reads). One sqlite is produced per mzML in <project>/distributions/, and
each file is processed in its own subprocess so memory is released between files
(the detector is memory-heavy on large runs).

Usage:
    python index-distributions.py --project /path/to/PXDxxxxx
    python index-distributions.py --mzml-dir <dir> --out-dir <dir> [--overwrite]

Extra args after `--` are forwarded to index_ms1.py, e.g.:
    python index-distributions.py --project P -- --max-charge 6 --workers 8
"""

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INDEXER = HERE / "distributions" / "index_ms1.py"
RUST_CRATE = HERE / "detector-rs"
RUST_BIN = RUST_CRATE / "target" / "release" / "ms1-detector"


def ensure_rust_binary():
    """Build the Rust detector if needed and return the binary path. The Rust
    engine is a validated drop-in for index_ms1.py (same sqlite output) and is
    much faster; the line + edge stages dominate runtime and are native there."""
    if not RUST_BIN.exists():
        print("building rust detector (cargo build --release)...", file=sys.stderr)
        result = subprocess.run(["cargo", "build", "--release"], cwd=str(RUST_CRATE))
        if result.returncode != 0:
            raise SystemExit("cargo build failed")
    return RUST_BIN


def find_centroid_mzmls(mzml_dir):
    files = sorted(mzml_dir.glob("*.centroid.mzML")) + sorted(mzml_dir.glob("*.centroid.mzml"))
    return [f for f in files if ".centroid." in f.name]


def stem_of(path):
    name = path.name
    for suffix in (".centroid.mzML", ".centroid.mzml"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--project", type=Path,
                        help="project dir; uses <project>/mzML and writes <project>/distributions")
    parser.add_argument("--mzml-dir", type=Path, help="override: directory of centroid mzML files")
    parser.add_argument("--out-dir", type=Path, help="override: directory to write the sqlite files")
    parser.add_argument("--overwrite", action="store_true", help="rebuild existing sqlite files")
    parser.add_argument("--engine", choices=("python", "rust"), default="python",
                        help="detector engine: python (index_ms1.py) or rust "
                             "(detector-rs/ms1-detector, faster, same sqlite)")
    parser.add_argument("forward", nargs="*",
                        help="args forwarded to the detector (place after --)")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.mzml_dir:
        mzml_dir = args.mzml_dir
    elif args.project:
        mzml_dir = args.project / "mzML"
    else:
        raise SystemExit("provide --project or --mzml-dir")

    if args.out_dir:
        out_dir = args.out_dir
    elif args.project:
        out_dir = args.project / "distributions"
    else:
        raise SystemExit("provide --project or --out-dir")

    rust_bin = ensure_rust_binary() if args.engine == "rust" else None
    if args.engine == "python" and not INDEXER.exists():
        raise SystemExit(f"missing indexer: {INDEXER}")
    if not mzml_dir.is_dir():
        raise SystemExit(f"missing mzML dir: {mzml_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    mzmls = find_centroid_mzmls(mzml_dir)
    if not mzmls:
        raise SystemExit(f"no *.centroid.mzML found in {mzml_dir}")

    print(f"{len(mzmls)} files -> {out_dir}")
    failures = []

    for i, mzml in enumerate(mzmls, 1):
        out = out_dir / f"{stem_of(mzml)}.distributions.sqlite"
        if out.exists() and not args.overwrite:
            print(f"[{i}/{len(mzmls)}] skip (exists): {out.name}")
            continue

        if args.engine == "rust":
            cmd = [str(rust_bin), str(mzml), "--out", str(out)]
        else:
            cmd = [sys.executable, str(INDEXER), str(mzml), "--out", str(out)]
        if args.overwrite:
            cmd.append("--overwrite")
        cmd.extend(args.forward)

        print(f"[{i}/{len(mzmls)}] {mzml.name}")
        # Separate process per file -> memory is freed between files.
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"  FAILED ({result.returncode}): {mzml.name}", file=sys.stderr)
            failures.append(mzml.name)

    if failures:
        raise SystemExit(f"{len(failures)} file(s) failed: {', '.join(failures)}")
    print("done")


if __name__ == "__main__":
    main()
