# Master specification & roadmap — motif-quantification viewer & pipeline

The complete, detailed spec for everything discussed, preserved at full fidelity,
with status on each item. Single source of truth (supersedes `viewer/PLAN.md`).

Status legend: ✅ done · 🟡 partial · ⬜ todo.

---

## 0. Reference files (`/examples`, user-provided)

| Reference file | Used for |
|---|---|
| `examples/panel-3-plot.py` | Panel 3 **MS1** charge-comparison grid (the 8-row × N-charge plot from the charge-state-determination code) |
| `examples/linemodel.py` | Distribution generation — **stage 1**, the line model |
| `examples/distributionassembly.py` | Distribution generation — **stage 2**, isotope distribution assembly |
| `examples/chargehandling.py` | Distribution generation — **stage 3**, charge-state grouping into analytes |
| `examples/peptidefragmentscoring.py` | **MS2** fragment isotopic distributions (b/y ion compositions, fragment descending-products, ion scoring) |
| `examples/libraryadditions.py` | **MS1** peptide isotopic distributions (`descending_partial_products`) |

The faithful reference algorithm is the authority for the pipeline; the active
`distributions/index_ms1.py` must be made to reproduce it (keeping the sqlite output).

---

# TAB 1 — MS VIEWING

## 1.1 Dock layout & panes
- ⬜ A **default pane setup**; provide a **"reset to default" pane** option that returns to this default.
- ✅ Panes are **completely movable and resizable**: Ctrl+click+drag to move, drag edges/corners to resize.
- ✅ Layout **auto-saves** and **remembers the user's prior personalization** for next launch (geometry + dock state).
- ✅ Window geometry persists; opening a file no longer resets the layout.
- Default arrangement: left = the 3 lists; middle column = Panel 1 (top) / Panel 2 (mid) / Table 1 (below Panel 2); right column = Panel 3 (with Table 2 below it when MS2). ✅ structure, ⬜ "reset to default" exactly matching this.

## 1.2 Left lists (single-file view)
- ✅ Three lists: **list 1 = proteins**, **list 2 = peptides**, **list 3 = PSMs**. **No file column / no file info** — this is a single-file view.
- ✅ A **file selector** picks the file first (must be chosen; defaults to first file) — required so the lists aren't empty.
- ✅ Each list title has an **"All" button**.
- ✅ Cross-linking: click a **protein** → peptide list shows that protein's peptides; click a **peptide** → protein list shows its proteins **and** PSM list shows its PSMs; same relationship PSMs↔peptides.
- ✅ **All** button restores the full unfiltered list (stops showing only associated entries) and **preserves the list's current selection + scroll position**.
- ✅ Peptides with no PSM in this file (LFQ-only / quantified-but-not-identified) are labeled — this is expected, not a bug.

## 1.3 Panel 1 — spectrum (2D⇄3D, centroid⇄profile)
- ✅ Starts as **non-profile (centroid)** and **non-3D (2D)**.
- ✅ Two **toggle buttons**: one shows **"3D"/"2D"**, one shows **"profile"/"centroid"** — each label shows **what it will switch to** when pressed. ✅ Buttons are **fixed size/location** so they don't move when the text changes.
- ✅ So Panel 1 is either 2D or 3D, and either profile or centroid.
- **2D view**: ✅ m/z on **x-axis**, intensity on **y-axis**; ⬜ different lines/distributions **colored differently** (from the sqlite distributions).
  - ✅ Shows **every datapoint** (not averages). Profile = continuous per-scan curves; centroid = dots.
  - ⬜ Profile dots/curves should **scale dynamically on zoom** so they stay visible (currently can get sparse/laggy).
  - ✅ Only the **m/z (x) axis** is interactive when 2D — **no vertical drag**; y auto-scales.
  - ✅ Wheel **over the y-axis label strip** (left of the axis) scrolls **intensity (y)**; ✅ pinned so the **baseline stays at 0** (data is always > 0; the baseline never lifts).
- **3D view**: m/z, time, intensity.
  - ⬜ **Grid lines on the m/z and time axes** marking where each scan was; ⬜ measured datapoints shown as **spherical points** along the grid at their intensity height.
  - ⬜ Points **colored by intensity** via a **gradient**; ⬜ the **3D peaks are a continuous surface** built from the **area between the 3D datapoints** (the surface legitimately represents each time×mass datapoint), with the **datapoints keeping their own point color** on top.
  - 🟡 surface renders; ⬜ make it the true point-to-point area-fill, and ⬜ profile points must **align with the surface** (was misaligned).
  - **3D bugs (reported):** ⬜ remove the **central gnomon/axis line** (the giant line aiming up) — only m/z and time need labeling; ⬜ axis **labels must be theme-adaptive** (currently only visible in dark mode); ⬜ labels must **attach to the axis ends** (they float freely now); ⬜ it must **start x-aligned with Panel 2** orientation; ⬜ add an **"align/reset 3D" button** to snap it back; ⬜ fix the **upside-down / stuck orientation** (can't recover currently); ⬜ reduce **zoom lag** (decimate/cap on zoom, reuse cached region).
  - ✅ In 3D, moving/rotating **does not move the m/z or time axes** — it only changes the 3D perspective. ✅ When the perspective is rearranged and Panel 2 moves, the perspective stays the same even though the data/axes change.

## 1.4 Panel 2 — m/z × time map
- ✅ The **mass and time window of Panel 1 = the exact region shown in Panel 2.** Panel 2 is zoomable/movable and **shifts everything in Panel 1**.
- ✅ Axes: m/z on **x** (swapped to align with Panel 1's x), time on **y**. (Original ask was m/z on y / time on x; swapped per later request so Panel 1↔2 x-axes align — ✅ now aligned.)
- ✅ Moving Panel 1 (in 2D) pans horizontally on the **mass axis** and moves Panel 2's mass axis too (synchronized). ✅ Moving Panel 2 realigns Panel 1's (2D) mass axis.
- ✅ Zoom/drag on Panel 2 (both axes + scroll-to-zoom) **reloads the data** for the viewed region.
- ✅ Render as **connect-the-dots**: individual datapoints as **dots** for time vs mass, **plus thin connecting lines**; the **line is thinner than the dots and the same color**; **each distribution a separate color** (colors from the sqlite). Points not in any distribution = faint gray. (Like the matplotlib time-vs-mass "connect the dots" reference.)
- ⬜ Data should **not disappear on zoom** — investigate (re-extract on zoom + the distribution over-generation may be making artifacts look like data).

## 1.5 MS2 strip (left of Panel 2)
- ✅ A **tall thin plot just to the left of Panel 2** (fits in the space Panel 1's wider y-axis labels create). It shares Panel 2's **time (y) axis**.
- ✅ MS2 scans shown as **horizontal lines that align with the time axis**, clickable.
- ⬜ **Bug:** the MS2 RT lines **disappear when zooming** and look like arrows. They must get **WIDER, not thinner, as you zoom** (perspective-correct: render in RT data-space, not fixed-pixel strokes), and stay clickable.

## 1.6 Distribution & line selection + coloring
- ⬜ From **Panel 1 or Panel 2**, the user can **select a distribution**.
- ⬜ After selecting a distribution, the user can **select a single line trace** within it (clarifies the selection, **replaces** the distribution selection); **clicking again brings it back to the distribution**.
- ⬜ A selected distribution/line is **highlighted in a color**; **all other distributions** are a different color. Both colors **user-selectable** via a **color settings drop-down**.
- ⬜ The selection is colored the same in **the 2D distribution, the 3D distribution, AND the Panel 2 line-style distribution** simultaneously.
- ⬜ In 3D, the **selected distribution's area** is colored as a **continuous surface** (not the individual datapoints); the datapoints keep their selected point color.

## 1.7 Charge search
- 🟡 With a distribution or line selected, **"charge search"**: charge taken from the distribution, or **assumed 1 by default** if none.
- ✅ **Left/right arrows** step to a m/z to look for a mass distribution at the **same time point**, **one charge higher/lower**.
- 🟡 It **"locks on"** to that m/z and **returns the user to the original** when they navigate back.
- ⬜ If a **legitimately-marked distribution of the same charge exists in the sqlite**, it **lights up the same color** and is **selected by default**.
- 🟡 The **assumed charge is user-editable** (type it, then navigate). ⬜ Add the manual charge field on the top bar.
- ⬜ Navigation **"locks"** in that m/z navigation **unless** Panel 2 is moved, the charge is overwritten, or either panel is moved.
- ✅ **Back / Forward history buttons** track the user's navigation so they can re-trace actions even after moving a pane or changing charge.

## 1.8 Table 1 (below Panel 2) — line metrics
- 🟡 Columns, **one row per line of a distribution**: **trapezoidal AUC**, **max intensity**, **sum intensity**, **number of data points**, **min time**, **max time**, **min mass**, **max mass**, **mean mass**, **retention time (timepoint of the highest-intensity point)**.
  - 🟡 AUC/max-I/n-pts/min-max time/min-max-mean mass/RT come from sqlite; ⬜ **sum intensity** still to add.
- ⬜ **If** the line/distribution is assigned to a peptide/protein from the Sage search, also show: **peptide**, **protein(s)**, **q value**, **file(s)** — where **file(s) is a click-to-expand dropdown listing all files VERTICALLY** for the files of the group being visualized.
- ⬜ When Panel 1 switches to **profile**, Table 1 **adjusts accordingly** (profile-mode metrics, see 1.9).

## 1.9 Profile ↔ centroid peak linking ("SUPER important")
- ⬜ When Panel 1 switches to **profile**, all the **peak-finding processes from `centroid-mzml.py` run in the background** and **link the output centroid point** of the matched distributions/lines **to EVERY profile datapoint that that specific centroided point came from**.
- ⬜ So **all profile datapoints are accounted for** this way, in Table 1 and the visualizations, and their **colors are assigned based on the centroid data** that each peak-finding process assigns the profile data to.

## 1.10 Panel 3 — MS1 view (charge-comparison grid)
- 🟡 When an **MS1 distribution** is clicked, Panel 3 shows the plot from `examples/panel-3-plot.py`: the **charge-state-determination grid** that **links multiple distributions (charge states of one analyte) into one plot of many comparisons** — columns = charge states, **8 rows**: retention time / peak area / charge distances / cross-charge / intensity sum % / adjacency / ppm-to-mean / ppm-error.
  - 🟡 grid renders; ⬜ make it **faithful to `panel-3-plot.py`** (per-row scales — peak-area log, charge-distance ylim, cross-charge log, intensity-sum% log, adjacency symlog; the white-on-gray styling adapted to theme; spine hiding; per-column titles `z(distid)`); ⬜ fix colors/readability.
- ⬜ Use **`descending_partial_products` (`libraryadditions.py` / `isotopes.py`)** to compute the **expected isotopic distribution of the peptide** if a peptide is being viewed and found via the search, and **twin-plot** it with the experimental on **different x-axes**, scaled so the **theoretically-most-abundant isotope and the max experimentally-measured datapoint are at the same height**.
- ✅ All MS1 Panel 3 plots that **share the mass x-axis are synchronized** when moved; ✅ **double-click resets** them all.
- ⬜ Remove the leftover **"Panel 3"** title text.

## 1.11 Panel 3 — MS2 view
- ⬜ **Bug:** currently **can't reach the Panel 3 MS2 view** — clicking an MS2 point must reliably switch Panel 3 to the MS2 spectrum.
- ⬜ Panel 3 MS2 is triggered by the user **clicking a sampled MS2 point**.
- ⬜ **MS2 points must be visible in both the 2D and 3D Panel 1, and in Panel 2**, standing out as a **'start' point / a distinct color**.
- ⬜ When an **identified peptide is assigned to that MS2 spectrum OR to the distribution sampled during that MS2 scan** (link the two via the search info if not already linked), **visualize the theoretical distribution of that specific MS1 progenitor**.
- ⬜ Label the **MS1 progenitor isotopic distribution** and its **fragment isotopic compositions**, labeling the **ions by both type and number** (use `examples/peptidefragmentscoring.py`).
- ⬜ Below the MS2 plot: **sequence coverage** of the peptide (from the sequence-coverage logic).
- ⬜ **Table 2** (only appears below Panel 3 when MS2): **other candidate PSMs** for that MS1 distribution **with their relevant sequence coverage** (optional panel). 🟡 dock + candidate listing exists; ⬜ coverage column.

## 1.12 Top-bar controls
- ✅ ± m/z and ± RT window controls (live).
- ✅ Reset zoom; ✅ charge ◀/▶ + history ⟲/⟳; ✅ theme toggle.
- ⬜ **Color settings drop-down** on the top bar: selected-distribution color, other-distributions color, **3D gradient min/max via a full color selector for both values**, and a **normal/log color-scale switch** (same fixed-size switch style as the 2D/3D and profile/centroid toggles).
- ⬜ Manual **charge** entry field; ⬜ **"align/reset 3D"** button.
- ⬜ **"loading… <context>"** label shown above **every** plot while its data worker runs, cleared when drawn — must happen **everywhere** (Panel 1 2D/3D, Panel 2, Panel 3, MS2), not just some places.

---

# TAB 2 — PROTEIN VIEWING ⬜
- ⬜ Show **entire protein sequences**; **peptide coverage** of individual proteins shown as **colored-in regions** of the protein sequence.
- ⬜ Region color **represents q-value**, using the **same color scale used for the 3D points** in the profile plots (shared q-value gradient).
- ⬜ View for a **single file**, OR the sequence **verticalized and displayed side-by-side against all other files**.
- ⬜ If the per-residue sequence is too hard to see, display **block chunks proportionate to the peptide length** they represent; peptides can be **zoomed in on** to clarify.

---

# TAB 3 — FILE-BY-FILE COMPARISON ⬜
- ⬜ Quantitative data **comparing peptides and proteins across files** — this is where the **quant work** is done.
- ⬜ **Time series** and **differential expression** analysis.
- ⬜ Reads the `experimental-setup` file (`filename,condition,fraction,replicate,pair_id`) to compare files; **modulate columns across files** as any treatment vs any other; **group multiple treatments hierarchically** to compare them any way.

---

# TAB 4 — MOTIF QUANTIFICATION ⬜
- ⬜ Quantify the **motifs** found via `index-motifs.py`. **Time series + DE at the MOTIF level**, where proteins with specific motifs are **represented by that motif**.
- ⬜ Look for **changes in expression of groups of proteins that all share a specific motif**.
- ⬜ A functional database links proteins to **skeleton motifs**; **organize the different peptides within these skeletons**.
- ⬜ **Include/exclude specific sequences** to narrow the protein lists; **save that narrowed motif set** within the database — a **new folder at the same level as `/distributions` and `/searches`** (e.g. `motif-sets/`).
- ⬜ Tabs 3 & 4 both **read `experimental-setup`** for grouping/contrasts.

---

# DISTRIBUTION DETECTION (the "SOLID process")

Symptoms (observed on `Tanya_Skin_NaCl_18`): 321,237 lines → 647,536 features →
**387,883 distributions** (more distributions than lines) → 338,377 analytes. The
distributions look like **artifacts**; data **disappears on zoom**. The user wants the
detection **reworked to a faithful version of the original** reference.

Root cause: the **active** pipeline `distributions/index_ms1.py` (writes the sqlite) is a
from-scratch reimplementation **missing whole reference stages**. The faithful reference is
`examples/linemodel.py` + `distributionassembly.py` + `chargehandling.py`. Keep the sqlite
schema (`distributions/store.py`); the GUI needs it. Validate each step on a real file
(user runs it — I can't run the pipeline here).

### Stage-1 line model (`linemodel.py`)
- ✅ acdiff acceptance (asymmetric, proton-spacing, `charge_tolerance`) — ported.
- ✅ intensity-step gating (`new_inc_limit`/`step_limit`) — ported.
- ⬜ **Adaptive `roundcutoff`** per scan (moving avg of the knee of sorted match distances; init 0); feeds `masswidthlimit = roundcutoff*2` (replace the current `mass_width_limit` proxy).
- ⬜ **Per-trace moving-avg diff (`madiff`)** + the **3-tier line acceptance** (in-range OR dist<range/2; if n≥`minmovinginds` accept if `nmadiff≤madiff`; n>1 → `d≤roundcutoff+range`, n==1 → `d≤roundcutoff*2`), replacing the flat ppm/abs gate.
- ⬜ **Dead-signal counter** (close after `>deadsignal`, halve on match) replacing the hard `max_gap_scans`.
- ⬜ NN **tie-break** for equidistant masses (prefer in-range, else closer in intensity).
- ⬜ **Line-correction merge**: after closing, merge sub-`subisomax` (0.01337851739·(1+chargetolerance)) lines with **non-redundant timepoints** via `intersection_merge` + directional-graph checks — recovers fragmented isotope envelopes.

### Stage-2 distribution assembly (`distributionassembly.py`)
- ⬜ **3-tier RT overlap** (encompassed → `>newinclimit`, partial → `>0.5`) replacing flat overlap/union.
- ⬜ **`overlap_counts`** path RT-geometry scoring + the **3-tier pair ranking** replacing single-best-edge-per-left. **This is the main fix for the over-generation** (distributions > lines).
- ⬜ Explicit **charge ±1 spread refinement**; **masswidthlimit clamping** of feature m/z ranges.

### Stage-3 charge handling (`chargehandling.py`)
- ⬜ Base-mass (`mz*z − proton*z`) alignment + **intensity-rank-order gating** (`Σ|rankdiff| ≤ size−1`).
- ⬜ **RT-overlap majority** gate (`overpass > matchables/2`); **adjacent-charge-only** search; active **nodist up/down matching**.

### Pipeline plumbing
- ✅ `index-distributions.py` driver + `execution.xsh` step (per-file sqlite in `<project>/distributions/`, subprocess per file for memory).
- ⬜ **Sanity checks/log**: warn when distributions > lines; assert members ≥ `min_distribution_members`; report per-stage counts.
- ⬜ Re-check `peaks.py` split-trace peak detection against reference behavior.

---

# RUST PORT ⬜
- ⬜ Once the Python distribution pipeline **and** the isotopic-distribution generation are
  **faithful and validated**, port the hot path (line model → assembly → charge handling,
  and `descending_partial_products` MS1/MS2 isotope generation) to **Rust** (the motif
  indexer is already Rust) — or whatever language is appropriate per stage. Must reproduce
  the validated Python exactly, write the same sqlite, bound memory, one process per file.

---

# Suggested execution order
1. **Distribution faithfulness** (line model roundcutoff + line-correction merge → assembly overlap_counts/ranking → charge rank/overlap gating). Top priority — fixes the artifacts.
2. **Tab-1 GUI bugs** (3D labels/axis/orientation/reset, MS2 RT lines wider-not-thinner, reach Panel 3 MS2, remove "Panel 3" text, loading-everywhere, profile zoom lag/upside-down) — can interleave with (1).
3. **Distribution selection + coloring** + color-settings dropdown.
4. **Panel 3 MS1** faithful to `panel-3-plot.py` + theoretical isotope overlay.
5. **Panel 3 MS2** (fragment isotopes type+number, sequence coverage, Table 2).
6. **Profile↔centroid peak linking** + Table 1 profile mode + peptide/protein/q/files columns.
7. **Tabs 2–4** (protein coverage, file DE, motif DE).
8. **Rust port**.
