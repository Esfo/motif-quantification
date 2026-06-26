# Phased plan & requirement coverage

This tracks every requirement from the spec against status, plus the distribution-
pipeline faithfulness work and the Rust question. Status: ✅ done · 🟡 partial ·
⬜ todo.

## Tab 1 — MS viewing

| Requirement | Status | Notes |
|---|---|---|
| Movable/resizable dock panes; Ctrl+drag move; edge/corner resize | ✅ | QDockWidget |
| Auto-save layout + remember personalization; reset-to-default | ✅ | geometry+state, autosave 4s, version-gated; Reset layout action |
| Panel 1: 2D⇄3D and centroid⇄profile fixed-size toggles | ✅ | labels flip, size fixed |
| Panel 1 2D: m/z × intensity, every datapoint; profile=curves | ✅ | wheel-over-y-axis scales intensity, baseline pinned at 0 |
| Panel 1 3D: m/z × time × intensity, points + surface, axis labels | ✅ | surface aligned to actual rt/mz |
| Panel 1 lines/distributions colored differently | 🟡 | Panel 2 colors by distribution; Panel 1 2D/3D color-by-distribution ⬜ |
| Panel 2: m/z × time map = panel 1 window; zoom/move reloads panel 1 | ✅ | window-driven, X-linked |
| Panel 2: dots + thin connecting lines colored by distribution | ✅ | from sqlite distributions |
| Panel 1↔2 m/z sync; 3D move only rotates (no axis move) | ✅ | |
| 3D color gradient: min/max color pickers + normal/log toggle (top bar) | ⬜ | color settings dropdown |
| Select a distribution, then a line within it (toggle) | ⬜ | needs click-pick on panels 1/2 |
| Selected distribution/line highlighted (user-chosen colors) across 2D/3D/panel-2 | ⬜ | color settings dropdown |
| 3D area-fill surface from the data points (true peak surface) | 🟡 | interpolated grid surface; per-distribution area-fill ⬜ |
| Charge search: ◀/▶ ±1 charge at same RT, lock on | ✅ | recenters m/z on same neutral mass |
| Charge search lights up a real same-charge distribution if it exists | ⬜ | needs DB match on stepped m/z |
| Manual charge override (type + navigate) | 🟡 | `set_charge` exists; toolbar field ⬜ |
| Back/forward navigation history | ✅ | ⟲/⟳ |
| Table 1: per-line metrics (AUC, max/sum I, n pts, min/max time, min/max/mean mass, RT) | 🟡 | from sqlite features; sum-intensity + peptide/protein/q/files-dropdown ⬜ |
| Table 1 adjusts when panel 1 → profile (peak-finding link) | ⬜ | profile↔centroid peak linking |
| Lists: proteins / peptides / PSMs, no file column, "All" buttons | ✅ | cross-link both directions; All preserves selection+scroll |
| Panel 3 MS1: charge-comparison grid (8 rows × charges) | 🟡 | rendered; colors/readability + cross-charge correctness need tuning |
| Panel 3 MS1 shared-mass plots synchronized; double-click resets | ✅ | per-column x-link + dbl-click reset |
| Panel 3 MS1: theoretical isotope dist (descending_partial_products) twin-plotted, height-matched | 🟡 | isotope overlay exists separately; integrate into grid ⬜ |
| MS2 points visible/clickable in panel 1 (2D+3D) and panel 2 | 🟡 | clickable MS2 strip left of panel 2; markers in panels 1/2 ⬜ |
| Panel 3 MS2: matched fragment ions w/ isotope comps, labeled type+number | ⬜ | needs fragment isotope port |
| Panel 3 MS2: sequence coverage below | ⬜ | coverage_print port |
| Table 2: other candidate PSMs + their coverage (MS2 only) | 🟡 | dock added, lists candidate PSMs; coverage staged |
| profile→centroid peak-finding links every profile point to its centroid line | ⬜ | the "SUPER important" one; runs centroid-mzml.py logic |

## Tab 2 — Protein viewing
Whole protein sequences; peptide coverage colored by q (shared q color scale);
single-file or verticalized side-by-side; block chunks when zoomed out. — ⬜ scaffolded tab only.

## Tab 3 — File-by-file comparison
Quant comparison across files; time series + DE; reads `experimental-setup`;
arbitrary hierarchical treatment grouping. — ⬜ scaffold; `experimental.py` reader done.

## Tab 4 — Motif quantification
DE at the motif level (proteins represented by shared skeleton motif); peptides
organized within skeletons; include/exclude refinement saved to a new sibling
folder. — ⬜ scaffold; `motifs.py` reader done.

## Distribution detection faithfulness (the "SOLID process")

A read-only audit compared `distributions/{linemodel,distributionassembly,
chargehandling}.py` to the reference. The structure matches; three divergences
plausibly degrade quality (need validation on a known-good file before applying):

1. **Intensity-step gating inverted** (`distributionassembly.py` ~240–278): defaults
   to `steplimit=0.5`, only tightens to `newinclimit=0.1` when reversing decrease→
   increase. Reference: strict `newinclimit` when increasing, lenient `steplimit`
   when decreasing. → loosens acceptance of rising isotopes.
2. **Asymmetric acdiff tolerance** (`distributionassembly.py` ~229–231): lower bound
   `-(diffcut*chargetolerance + masswidthlimit)` vs upper `diffcut + masswidthlimit`
   → preferentially rejects low-m/z isotopomers.
3. **Charge tolerance globally averaged** (`chargehandling.py` ~113–114):
   `ctol = mean(|proton − massdiffs|)` per charge instead of per-distribution → flattens
   local variation.

Medium: early-scan `roundcutoff` starts at 0 (inflated early), and early-line
acceptance hardened to `roundcutoff*2`.

## Rust question
The line model + distribution assembly + charge handling are the memory/time
hot spots on big files. A Rust port (the motif indexer is already Rust) would cut
memory and runtime, but is a multi-week effort and must reproduce the reference
exactly. Recommended: first restore faithfulness in Python (small, reviewable
diffs, validate on a reference file), then port the validated algorithm to Rust.

## Proposed phase order
- **P1 (now):** GUI bug-fixes ✅ (this batch).
- **P2:** distribution faithfulness — apply audited fixes behind validation on a reference file.
- **P3:** Panel-1/2/3 distribution **selection + coloring** + color settings dropdown (ties panels together).
- **P4:** profile↔centroid peak linking (Table 1 profile mode, point attribution).
- **P5:** Panel-3 MS2 (fragment isotopes + sequence coverage) + Table 2 coverage.
- **P6:** Tabs 2–4 (protein coverage, file DE, motif DE) on the experimental-setup + motif readers.
- **P7:** Rust port of the validated distribution pipeline.
