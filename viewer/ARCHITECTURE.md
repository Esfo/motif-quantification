# Viewer architecture & roadmap

This document is the map for the `/viewer` GUI: what each module does, the data
contracts it depends on, what is wired today, and what is deliberately staged.
It also records the proposed repo reorganization so the rebuild stays coherent.

## Running it

```bash
pip install PySide6 pyqtgraph numpy pyteomics lxml psims PyOpenGL distinctipy
python viewer/app.py                  # then double-click / Ctrl+O to open a folder
python viewer/app.py --reorganized /path/to/searches/reorganized
```

Open accepts **either** a project directory (containing `searches/reorganized`,
`distributions/`, `experimental-setup`) **or** a `reorganized` directory
directly. When given a project dir the viewer auto-detects the distributions
sqlite and the experimental-setup file.

## Module map (`viewer/`)

| module | responsibility |
|---|---|
| `app.py` | CLI entry, `QApplication`, SIGINT→quit |
| `main_window.py` | window chrome, 4-tab shell, folder open/reload, theme, dock-layout persistence, distributions/experimental auto-detect |
| `ms_viewer_tab.py` | **Tab 1** dock workspace (lists, panels 1–3, table 1) |
| `proteins_tab.py` | **Tab 2** protein sequence viewer: file selector + FDR + protein list, panel 1 (horizontal sequence, peptide rectangles coloured by q-value, All button) over panel 2 (same protein verticalized per-file, click a column to load it in panel 1) |
| `quant_tab.py` | **Tab 3** Quantitative Comparisons: role-assignment design panel + group/titration visualization + DE volcano over the top half; sortable feature table with a Peptides⇄Proteins switch across the bottom |
| `quant_model.py` | builds/caches the feature×file quantity matrix (peptide or protein level; protein roll-up sum/median/unique) from the reorganized quant tables |
| `de_stats.py` | pure-Python differential-expression stats: log2, Welch's & paired t-tests (Student-t p via incomplete beta, no scipy), Benjamini-Hochberg FDR |
| `session.py` | reads the `reorganized/` tables (files, PSMs, peptides, proteins, quant) + the project FASTA (protein sequences) and in-silico tryptic digestion (Tab 2) |
| `mzml_store.py` | indexed mzML reader: metadata (no array decode), random-access scans, `extract_xics`, `extract_region` |
| `region_view.py` | standalone bounded 2D heatmap + 3D region (`RegionWorker` thread is reused by Tab 1) |
| `plots.py` | pyqtgraph helpers (spectrum sticks, points, traces, bars, short labels) |
| `theming.py` | light/dark palettes applied live to plots + GL |
| `chemistry.py` | elemental masses/abundances, AA & fragment compositions, derived isotope tables |
| `isotopes.py` | theoretical peptide isotope distribution (ported `descending_partial_products`) |
| `distributions_db.py` | read-only reader for the MS1 distributions sqlite |
| `experimental.py` | reader for `experimental-setup` (design groups/contrasts) |
| `motifs.py` | reader for the motif skeleton index (motifs.tsv + varint postings.bin + proteins.tsv) |

## Data contracts

### `searches/reorganized/` (from `reorganize-results.py`)
- `files.tsv` — filename, run_dir, mzml_path, fraction, replicate, counts
- `peptides.tsv`, `proteins.tsv` — global tables (proteins.tsv packs all peptides
  into one field → csv field limit is raised in `session.py`)
- `by_file/<run>/{psms,peptides,proteins,peptide_quant,scan_lookup}.tsv`
- profile mzML always sits next to its centroid (`X.mzML` beside `X.centroid.mzML`)

### `distributions/*.sqlite` (from `distributions/store.py`)
- `distributions(distribution_id, charge, neutral_mass, mono_mz, rt_*, ms1_*, n_members, score, quality)`
- `distribution_members(distribution_id, feature_id, isotope_index, member_score)`
- `features(feature_id, line_id, mz_*, rt_*, ms1_*, height, area, n_points, quality)`
- `lines`, `scans`, `analytes`, `analyte_members`
- **A "line" in the spec = a `feature` row** (one isotope trace's peak); a
  "distribution" = its `distribution_members`. Raw per-point arrays are *not*
  stored — the 3D points come from `mzml_store.extract_region` on the window.

### `experimental-setup` (csv)
`filename,condition,fraction,replicate,pair_id` — `filename` matches the mzML
stem. `experimental.py` exposes `group_by(*cols)` / `filenames_for(**filters)`
so tabs 3/4 can define arbitrary groups and contrasts.

### `…/motifs/<index>/`
`build_info.tsv`, `motifs.tsv` (motif_id, motif_text like `A.....G`,
protein_count, posting_offset, posting_bytes), `postings.bin` (LEB128 varint,
count-prefixed, delta-encoded protein ids), `proteins.tsv`.

## Tab 1 — wired vs staged

**Wired now:** dock layout (drag/float/resize, 3-column default, versioned
persistence + reset); single-file selector; three cross-linked lists with All
buttons (All preserves selection+scroll; LFQ-only peptides labelled); evidence
read off the UI thread (`EvidenceWorker`, latest-wins). Panels are **window-driven**
(`self.window = [mz_min, mz_max, rt_start, rt_end]` is the source of truth):
- Panel 1 2D shows **every datapoint** in the window (m/z vs intensity,
  `extract_points`); profile draws per-scan curves, centroid draws dots; only the
  m/z axis is interactive, wheel over the y-axis strip scrolls intensity.
- Panel 1 3D renders the raw points (coloured by intensity) + interpolated surface
  (mapped to actual rt/mz so points and surface align), with m/z/time labels.
- Panel 2 is a **connect-the-dots** view: raw points (m/z x, RT y) as dots + thin
  connecting lines, coloured per sqlite **distribution** (stable colour per
  distribution_id); X-linked to panel 1, drag/zoom reloads the window.
- A thin clickable **MS2 strip** sits left of panel 2 (horizontal lines at each
  MS2 RT, shared y); clicking loads that MS2 spectrum into panel 3.
- Panel 3: **charge-comparison grid** (columns = analyte charge states, rows =
  retention time / peak area / charge distances / cross-charge / intensity sum %
  / adjacency / ppm-to-mean / ppm-error) when the match maps to a distribution;
  else the isotope overlay; MS2 spectrum on an MS2 click. Cross-charge rows use
  base-mass (`mz*z - proton*z`) nearest alignment; the RT row uses the window
  raw points filtered to each feature.
- Table 1: line metrics from the sqlite distribution members.
- Dock layout autosaves (4s) and restores (version-gated).

### Panel 3 MS1 — charge-comparison grid (staged, next major build)

Replaces the single isotope overlay with the multi-distribution charge grid from
the charge-state-determination code (columns = charge states of one analyte,
rows = metrics). Data source: `DistributionsDB.charge_group(distribution_id)` →
`{charge: {distribution, features}}` via `analyte_members`. Row → data mapping:
- **retention time** — raw points per feature (needs `extract_points` per line)
- **peak area** — `features.area` (log y)
- **charge distances** — `diff(feature mz) * charge` (≈ 1.0)
- **cross-charge** — intensity ratios of aligned isotopes across charges
- **intensity sum %** — feature height / summed aligned intensities
- **adjacency** — adjacent isotope intensity ratios (symlog)
- **ppm to mean / ppm error** — base-mass alignment errors across charges

Old terms → new schema: `chargedistgroups[analyte]` = `charge_group`;
`distributionmasses/intensities` = member `features` (`mz_mean`/`height`);
`linesofdistributions` = `distribution_members`; `trackedgroups`/`regions` raw
points = `extract_points` over each feature's window. Alignment across charges
reuses the base-mass (`mz*z - proton*z`) nearest-match logic from the original.

**Staged (each has a home in the structure):**
1. **Distribution/line selection + colouring** — click a distribution in panel
   1/2, then a line within it; colour it (and its panel-2 trace) a user-chosen
   colour, others a second colour. Needs the sqlite overlay drawing + a colour
   settings dropdown on the top bar.
2. **Charge search** — left/right arrows step the m/z by ±1 charge at the same
   RT, lock on, light up a real matching distribution if one exists, allow a
   manual charge override, and keep back/forward navigation history.
3. **3D as a true surface** — colour the *area between* datapoints as a
   continuous surface (datapoints keep their own colour), grid lines at every
   scan/mass, spheres at measured points. Colour gradient with user-set
   min/max colour pickers and a normal/log toggle on the top bar.
4. **Profile↔centroid peak linking** — when Panel 1 switches to profile, run the
   `centroid-mzml.py` peak-finding so every profile datapoint is attributed to
   the centroid point it came from, and colour/aggregate accordingly. This is
   the heaviest staged piece and drives Table 1's profile-mode metrics.
5. **Panel 3 MS2** — clicking a sampled MS2 point shows fragment-ion isotope
   compositions (labelled by ion type+number) over the MS2 spectrum, with
   sequence coverage beneath; Table 2 lists other candidate PSMs for that MS1
   distribution with their coverage. Fragment isotopes reuse the
   `fragmentation_compositions` + fragment descending-products logic (to be
   ported into `isotopes.py` next to the MS1 path).

## Tabs 2–4 (scaffolded)

- **Proteins** — whole sequences with peptide coverage coloured by q-value on the
  shared q-value colour scale; single-file or verticalized side-by-side; block
  chunks proportional to peptide length when zoomed out.
- **Quantitative Comparisons** (wired) — quantitative peptide/protein comparison
  across files with differential expression. Reads `experimental-setup` and lets
  **any** column take a *role* (Group / Series-axis / Replicate / Pair / Ignore),
  so the contrast and the titration/time ordering are user-defined rather than
  hard-coded to condition/fraction/replicate. Top half: design/role panel (left),
  per-group or along-series visualization of the selected feature (center), DE
  volcano (right, log2FC vs -log10 FDR, click a point to select). Bottom half: a
  sortable feature table with a Peptides⇄Proteins switch (protein roll-up
  sum/median/unique) that re-analyzes everything above. Stats: Welch's or paired
  t-test on log2 quantities + BH-FDR (`de_stats.py`, `quant_model.py`). Staged
  next: worker-thread DE for large peptide sets, saved contrasts, MA/heatmap views.
- **Motifs** — DE at the motif level (proteins represented by shared skeleton
  motif); include/exclude refinement saved to a new sibling folder (e.g.
  `motif-sets/`) alongside `distributions/` and `searches/`.

## Isotope math provenance

`isotopes.py` is a faithful port of the provided `individual_element_binomial_walk`
/ `descending_partial_products` exploration onto `chemistry.py`. Validated:
angiotensin II `DRVYIHPF` → neutral monoisotopic **1045.5345**, base peak at M+0,
z=2 mono m/z **523.77**, ~1.0029 Da group spacing. Changes from the base script:
extracted globals into `chemistry.py`, removed the file/pickle I/O and matplotlib,
and exposed `peptide_isotope_distribution` / `peptide_isotope_mzs`. The
fragment-ion isotope path (for Panel 3 MS2) is not ported yet.

## Proposed repo reorganization

Current redundancy: elemental constants and the isotope walk exist in
`distributions/elementalcomponents.py` and the provided scripts; the viewer now
has its own `chemistry.py`/`isotopes.py`. Proposed end state (not yet executed to
avoid breaking the working pipeline):

```
chem/                     # shared, importable by pipeline + viewer
  elements.py             # = viewer/chemistry.py (single source)
  isotopes.py             # MS1 + fragment isotope distributions
pipeline/                 # the existing scripts, unchanged behavior
  centroid_mzml.py, index_motifs.py, reorganize_results.py, fasta_updater.py, ...
  distributions/          # index_ms1.py, store.py, ... importing chem/
viewer/                   # GUI only, importing chem/
```

Migration is incremental: keep `distributions/elementalcomponents.py` as a thin
re-export of `chem/elements.py` so nothing breaks, then delete duplicates once
imports are switched. Until then `viewer/chemistry.py` is the source of truth for
the GUI and is kept numerically identical to the pipeline constants.
