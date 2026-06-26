# Master roadmap — motif-quantification viewer & pipeline

Single source of truth for everything discussed. Supersedes `viewer/PLAN.md`.
Legend: ✅ done · 🟡 partial · ⬜ todo.

## Reference files (user-provided in `/examples`)

| Reference | Drives |
|---|---|
| `examples/panel-3-plot.py` | Panel 3 MS1 charge-comparison grid |
| `examples/linemodel.py`, `examples/distributionassembly.py`, `examples/chargehandling.py` | Distribution generation (line model → assembly → charge handling) |
| `examples/peptidefragmentscoring.py` | MS2 fragment isotopic distributions |
| `examples/libraryadditions.py` | MS1 (peptide) isotopic distributions |

(The faithful reference algorithm currently also lives in `distributions/linemodel.py`
etc., which are not run; `/examples` copies are the authoritative reference.)

---

## ✅ Already done

- Viewer foundation: PySide6 app, open project/reorganized folder, last-dir memory,
  light/dark theme, Ctrl+C quit, scratch crash fixes.
- Shared chemistry + MS1 peptide isotope math (`viewer/chemistry.py`, `viewer/isotopes.py`),
  validated (angiotensin II → 1045.5345, z=2 523.77).
- Readers: distributions sqlite, `experimental-setup`, motif index (`viewer/*`).
- Tab 1 dock workspace: movable/resizable panes, layout+geometry persistence (autosave,
  version-gated), reset layout.
- Single-file lists (protein/peptide/PSM) with All buttons (preserve selection+scroll),
  cross-linking, LFQ-only labelling.
- Window-driven panels: every-datapoint Panel 1 (2D), threaded evidence (no UI freeze),
  Panel 1↔2 m/z X-link + x/y-axis alignment, profile vs centroid.
- Panel 2 connect-the-dots coloured by sqlite distribution.
- Panel 1 intensity wheel-scroll pinned to baseline 0.
- Charge search arrows + back/forward history.
- Panel 3 charge-comparison grid (8 rows × charges), per-column mass-x sync + dbl-click reset.
- MS2 strip left of Panel 2 (clickable) + Table 2 dock.
- Pipeline: top-level `index-distributions.py` driver + `execution.xsh` wiring (per-file sqlite).
- Pipeline: ported faithful acdiff + charge-tolerance + intensity-step gating into `index_ms1.py` edge acceptance.
- Docs: `viewer/ARCHITECTURE.md`, `distributions/PIPELINE_FAITHFULNESS.md`.

---

## Stage 1 — Distribution faithfulness (the "SOLID process")  ⬜ priority

Symptoms: distributions look like artifacts; **more distributions (387,883) than lines
(321,237)** → features reused across many spurious distributions; data "disappears" on
zoom. Root cause: `index_ms1.py` (active, writes sqlite) is a reimplementation missing
whole reference stages. Goal: make it faithful to `examples/linemodel.py` +
`distributionassembly.py` + `chargehandling.py`, keep the sqlite schema. Validate each
step on a real file (user runs it).

1. ✅ acdiff acceptance (asymmetric, proton-spacing, charge_tolerance).
2. ✅ intensity-step gating (new_inc_limit / step_limit).
3. ⬜ **Adaptive `roundcutoff`** per scan (moving avg of the knee of sorted match
   distances); feeds `masswidthlimit = roundcutoff*2` (replace the `mass_width_limit`
   proxy currently used by acdiff).
4. ⬜ **Moving-avg diff (`madiff`) per trace** + the **3-tier line acceptance**
   (in-range / madiff-stable / roundcutoff tiers) replacing the flat ppm/abs gate.
5. ⬜ **Dead-signal counter** (close after `>deadsignal`, halve on match) replacing the
   hard `max_gap_scans`; add NN **tie-break** (in-range, then intensity).
6. ⬜ **Line-correction merge** (sub-`subisomax` 0.013377 lines, non-redundant timepoints,
   intersection_merge) — recovers fragmented isotope envelopes.
7. ⬜ **3-tier RT overlap** (encompassed `>newinclimit`, partial `>0.5`) + **`overlap_counts`**
   path RT-geometry scoring + **3-tier pair ranking** in distribution assembly.
8. ⬜ **Charge handling**: base-mass alignment + **intensity-rank-order gating**
   (`Σ|rankdiff| ≤ size-1`), **RT-overlap majority** gate, **adjacent-charge-only** search,
   active nodist up/down matching → fixes the distribution over-generation.
9. ⬜ Sanity asserts/log: distributions ≤ feature pairs; warn if distributions > lines.
10. ⬜ Re-validate `peaks.py` split-trace peak detection vs reference expectations.

---

## Stage 2 — Tab 1 GUI bugs (reported)  ⬜

11. ⬜ **Panel 3 3D plot**:
    - remove the central GL axis line / xyz gnomon; keep only m/z + time labels;
    - labels **theme-adaptive** (visible in light & dark);
    - labels **attached** to the axis ends (not free-floating);
    - **start x-aligned** with Panel 2 orientation, and add an **"align/reset 3D" button**;
    - fix **upside-down / stuck orientation** (clamp/reset camera; the reset button fixes it).
12. ⬜ **Profile 3D zoom lag** — decimate points / cap surface resolution on zoom; reuse
    cached region; avoid full re-render per wheel tick.
13. ⬜ **MS2 RT lines (Panel 2 strip)**: must get **WIDER, not thinner** when zooming RT,
    and not disappear. Render as RT-data-space rectangles/`InfiniteLine`s (perspective-correct),
    not fixed-pixel strokes; keep clickable.
14. ⬜ **Reach the Panel 3 MS2 view** — clicking an MS2 point isn't switching the stack to
    the MS2 spectrum reliably; ensure MS2 click always shows Panel 3 MS2 + Table 2.
15. ⬜ Remove the leftover **"Panel 3" / "loading Panel 3"** title text.
16. ⬜ **"loading… <context>" everywhere**: every plot shows a loading label while its
    worker runs (Panel 1 2D/3D, Panel 2, Panel 3, MS2), cleared on draw.

---

## Stage 3 — Tab 1 distribution interaction + coloring  ⬜

17. ⬜ **Select a distribution** by clicking on Panel 1/2; then select a **line** within it
    (toggle back to distribution on re-click).
18. ⬜ **Highlight** the selection across Panel 1 (2D+3D), Panel 2, and Panel 3 in a
    user-chosen color; all other distributions a second user-chosen color.
19. ⬜ **Color settings dropdown** on the top bar: selected/other colors, 3D gradient min/max
    full color pickers, normal/log color scale toggle (same fixed-size toggle style).
20. ⬜ 3D **area-fill surface** colored per distribution (continuous surface between the
    distribution's points), datapoints keep their point color.
21. ⬜ Charge search: light up a **real matching distribution** of the stepped charge if one
    exists in the sqlite (auto-select it); manual charge text field on the top bar.
22. ⬜ Table 1: add **sum intensity**, and when a line/distribution maps to a Sage ID:
    **peptide / protein(s) / q-value / files** (files = click-to-expand vertical list).

---

## Stage 4 — Profile ↔ centroid peak linking  ⬜ ("SUPER important")

23. ⬜ When Panel 1 switches to **profile**, run the `centroid-mzml.py` peak-finding so every
    profile datapoint is attributed to the centroid point/line it came from.
24. ⬜ Color profile points by their assigned centroid distribution; Table 1 recomputes in
    profile mode from the attributed points.

---

## Stage 5 — Panel 3 MS1 (full, per `examples/panel-3-plot.py`)  🟡

25. 🟡 8-row grid exists; ⬜ make it faithful to `panel-3-plot.py` (row scales: peak-area log,
    charge-distances ylim, cross-charge log, intensity-sum% log, adjacency symlog;
    white-on-gray styling adapted to theme; spine hiding; per-charge column titles `z(distid)`).
26. ⬜ Integrate the **theoretical MS1 isotope distribution** (`libraryadditions.py` /
    `isotopes.py`) as a twin-plotted overlay, height-matched to the experimental peak.

---

## Stage 6 — Panel 3 MS2 + sequence coverage  ⬜

27. ⬜ MS2 point click → Panel 3 MS2 spectrum with **matched fragment ions**, each labeled by
    **ion type + number**, with **fragment isotopic distributions** (`peptidefragmentscoring.py`).
28. ⬜ Link the sampled MS1 distribution to its MS2 scan's identified peptide via the search
    info (if not already linked).
29. ⬜ **Sequence coverage** display below the MS2 plot; **Table 2** = other candidate PSMs for
    that MS1 distribution with their sequence coverage.

---

## Stage 7 — Tabs 2–4  ⬜

30. ⬜ **Tab 2 Protein viewing**: whole protein sequences; peptide coverage colored by q-value
    (shared q color scale = the 3D gradient); single-file or verticalized side-by-side across
    files; proportional block chunks when zoomed out, zoomable.
31. ⬜ **Tab 3 File-by-file comparison**: quant across files; time series + differential
    expression; reads `experimental-setup`; arbitrary hierarchical treatment grouping/contrasts.
32. ⬜ **Tab 4 Motif quantification**: DE at the motif level (proteins represented by shared
    skeleton motif from the motif index); organize peptides within skeletons; include/exclude
    sequences to narrow motif sets; **save narrowed motif sets** to a new sibling folder
    (e.g. `motif-sets/`). Tabs 3 & 4 both consume `experimental-setup`.

---

## Stage 8 — Rust port  ⬜

33. ⬜ Once the Python distribution pipeline is faithful **and validated**, port the hot path
    (line model → assembly → charge handling, and the isotopic-distribution generation) to
    Rust (the motif indexer is already Rust). Reproduce the validated Python exactly; same
    sqlite output; bound memory; one process per file. Language: Rust unless a specific stage
    is better served otherwise.

---

## Suggested execution order
Stage 1 (3–10) → Stage 2 bugs → Stage 3 selection/coloring → Stage 5 panel-3 polish →
Stage 6 MS2 → Stage 4 profile linking → Stage 7 tabs 2–4 → Stage 8 Rust. Stage 2 GUI bugs
can interleave with Stage 1 since they're independent.
