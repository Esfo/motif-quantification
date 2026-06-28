# ms1-detector (Rust)

A native Rust port of `distributions/index_ms1.py` — the MS1 isotope-distribution
detector that writes the `*.distributions.sqlite` the GUI reads. It is a
**validated drop-in**: same algorithm, same SQLite schema, byte-equivalent output
(parity-tested to floating-point noise), and much faster (the line and edge
stages, which dominate runtime, are native + rayon-parallel).

## Build

```
cd detector-rs
cargo build --release          # binary at target/release/ms1-detector
```

## Run

```
ms1-detector <file.centroid.mzML> --out <file.distributions.sqlite> --overwrite
# options: --min-intensity <f>  --threads <n>  --progress <n>
```

Or through the project wrapper (builds the binary on first use):

```
python index-distributions.py --project /path/to/PXDxxxxx --engine rust -- --overwrite
```

## Parity

The Python pipeline stays the reference implementation. Each stage was validated
against it:

- `peaks.rs`  vs `peaks.py` — 300 random arrays, 0 mismatches.
- line model  vs `LineModel` — synthetic scans, identical line/feature values.
- edges+dists vs `build_isotope_edges`/`build_distributions` — identical counts,
  charge histogram, and per-distribution rows.
- end-to-end  — real mzML through both; identical table counts + charge
  histogram; distribution/feature rows match to < 1e-11.

If you change the Python algorithm, re-port the changed stage and re-validate
(the `*_probe` binaries under `src/bin/` and the harness in git history show how).

## Layout

- `config.rs`        — Config (mirrors the Python dataclass) + constants.
- `peaks.rs`         — port of `peaks.py` (axis_peaks).
- `linemodel.rs`     — tracker, history-weighted death, fragment merge, split.
- `distributions.rs` — edges (data-derived charge), distributions, competition, analytes.
- `store.rs`         — SQLite writer (replicates `store.py`).
- `main.rs`          — mzdata mzML reader + orchestration.
