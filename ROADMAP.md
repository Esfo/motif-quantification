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
| `examples/sequencecoverageconcept.py` | **Panel 3 MS2 sequence coverage** + **Table 2** coverage concept (which residues each matched ion covers, the coverage/divider-string logic). **NOT Tab 2's protein coverage — that is a different concept.** ✅ *file provided; ported to `viewer/coverage.py`* |

The faithful reference algorithm is the authority for the pipeline; the active
`distributions/index_ms1.py` must be made to reproduce it (keeping the sqlite output).

---

## 0.1 Original vision (the foundational goals)
The project searches shotgun proteomics files via **Sage** (run through `execution.xsh`)
and does **differential expression analysis via motifs**. The GUI must visualize all the
data the pipeline produces. The four founding goals:
1. ⬜ Visualize the **original MS1 profile data**, and put the **supposed (theoretical) MS1
   distribution against the experimental one** (theoretical isotope envelope vs measured).
2. ⬜ View the **MS2 fragments** and provide a way to **quantify**.
3. ⬜ **Group multiple proteins together via motif** and do **DE**.
4. ✅ Visualize the things the pipeline produces by file — the MS1 distributions found and the
   Sage search results.
(Isotopic-distribution *calculation* was deferred at the start; it is now in scope —
`viewer/isotopes.py` + `examples/libraryadditions.py`.)

## 0.2 Cross-cutting requirements (apply everywhere)
- ⬜ **Seamless & crisp** — overall polish is an explicit, standing requirement.
- ⬜ **Shared color gradient**: there is one user-defined **min/max color gradient** (the color-
  settings dropdown, 1.12). It colors the **3D points by intensity** AND is the **same q-value
  color scale** used for **Tab 2** protein coverage. Changing it changes both.
- ⬜ **Document everything built, and every change/deviation from the provided reference/base
  code** (`panel-3-plot.py`, `linemodel/distributionassembly/chargehandling.py`,
  `peptidefragmentscoring.py`, `libraryadditions.py`) — the user asked for this explicitly.
- 🟡 **Repo organization**: no redundant code; reorganize so it runs sensibly without losing
  functionality; isotopic/chemistry functions live in importable shared modules (done:
  `viewer/chemistry.py`, `viewer/isotopes.py`, removed dead `views.py`/`workspace.py`); ⬜ the
  broader pipeline reorg (shared `chem/` package; thin re-exports) per `viewer/ARCHITECTURE.md`.
- ✅ **PySide6 desktop** app (not web).
- Data locations (not hard-coded; auto-detected from the project): the **distributions sqlite**
  lives in `<project>/distributions/` (sibling of `/searches` and its `/reorganized`); an
  `*.inspect.json` may sit alongside it; the **motif index** has the structure
  `human-proteome-skeletons/{build_info.tsv,motifs.tsv,postings.bin,proteins.tsv}` under a
  motifs dir (e.g. `~/data/proteomics/motifs/`); `experimental-setup` sits in the project root.
- "**single-file view so far**" (Tab 1) — a multi-file Tab-1 mode may come later.

## 0.3 Miscellaneous small requests (mostly done — listed so nothing is lost)
- ✅ **Open via a folder dialog, not flags**; accept the **project folder** directly (auto-finds
  `searches/reorganized`, `distributions/`, `experimental-setup`); **remember the last location**.
- ✅ **Empty GUI by default** — start filled with empty widgets (not a blank placeholder);
  **double-click the empty area** to open the folder dialog (in addition to File ▸ Open).
- ✅ **Ctrl+C from the CLI quits** cleanly (no force-kill needed).
- ✅ **Instant light/dark theme** toggle (button press, applies to all plots + GL). Dark-mode axis/label/title text is **pure white** (was a dim grey the user found illegible). On toggle, the **data-bearing panels are re-rendered from cache** so theme-dependent colours follow the theme everywhere: Panel 1 & 3 plot **titles**, the **MS2 spectrum data** (was stuck grey in light mode → now theme fg), and the **charge-grid axes/labels/values** all adapt now (they were baked at draw time before).
- ✅ Removed the **"Panel 1"/"Panel 2"** dock title text (and now the **"Panel 3"** text too — see 1.10).
- ✅ **Bounded region** — the profile/region view is the ID's ± m/z / ± RT window, **never the
  entire spectrum** (which was unreadable).
- ✅ Removed the **orange 3D background** (was the height-color shader).
- ✅ **No UI freeze** — selection reads run on a worker thread; rapid changes are latest-wins. **Table 1's line-metric query** (distribution lookup + members on the sqlite) also runs on its **own background thread** (`Table1Worker`, latest-wins by token, opens its own read-only connection); panel 3's charge grid redraws if the distribution id resolves after it first drew. The per-scan reads inside `extract_region`/`extract_points` are now **read in parallel** (a thread pool of independent mzML readers — the base64+zlib decode releases the GIL, so this genuinely parallelises), speeding up the initial load and every grow/zoom/move. (Threads rather than processes: avoids per-process file re-indexing and array IPC, while still parallelising the decode.)
- ✅ Sync model decision: **"forget the old generic sync"** — synchronization is now per-shared-
  axis (Panel 1↔2 m/z; Panel 3 columns on mass), not a single global lock.
- ✅ Pipeline driver runs from a **top-level `.py`** (`index-distributions.py`) via `execution.xsh`,
  writing to the project's `/distributions/`.
- ⬜ A short **"how to run"** note / README for launching the viewer and the pipeline.
- ⬜ **General "seamless & crisp"** polish pass once features land (see 0.2).

---

# TAB 1 — MS VIEWING

## 1.1 Dock layout & panes
- ⬜ A **default pane setup**; provide a **"reset to default" pane** option that returns to this default.
- ✅ Panes are **resize-only** (per the user): drag the splitter edges to resize, but they can **no longer be moved, floated, re-docked, or closed** (all docks set to `NoDockWidgetFeatures`), so the layout can't be accidentally torn apart.
- ✅ Layout **auto-saves** and **remembers the user's prior personalization** for next launch (geometry + dock state).
- ✅ Window geometry persists; opening a file no longer resets the layout. **Bug fixed:** `app.py` was calling `window.resize(1600×950)` *after* the saved geometry was restored, clobbering it every launch — it now only applies the default size when nothing was restored, and the dock layout is re-applied once on first `showEvent` (a restore before the first show was being discarded as the nested tab settled).
- Default arrangement: left = the 3 lists; middle column = Panel 1 (top) / Panel 2 (mid) / Table 1 (below Panel 2); right column = Panel 3 (with Table 2 below it when MS2). ✅ structure, ⬜ "reset to default" exactly matching this.

## 1.2 Left lists (single-file view)
- ✅ Three lists: **list 1 = proteins**, **list 2 = peptides**, **list 3 = PSMs**. **No file column / no file info** — this is a single-file view.
- ✅ A **file selector** picks the file first (must be chosen; defaults to first file) — required so the lists aren't empty.
- ✅ Each list title has an **"All" button**.
- ✅ Cross-linking: click a **protein** → peptide list shows that protein's peptides; click a **peptide** → protein list shows its proteins **and** PSM list shows its PSMs; same relationship PSMs↔peptides.
- ✅ **Auto-load when a selection resolves to one item**: selecting a peptide with a single PSM auto-loads that PSM (no extra click); a protein that resolves to a single peptide auto-selects it (which can then auto-load its single PSM).
- ✅ **All** button restores the full unfiltered list (stops showing only associated entries) and **preserves the list's current selection + scroll position**.
- ✅ The lists show **only entries identified (with a PSM) in the selected file**. Peptides quantified-but-not-identified in this file (LFQ-only / match-between-runs transfers) and proteins with no file-identified peptide are **excluded** (per the user; supersedes the earlier "label LFQ-only" behaviour). Identified set = plain sequences from this file's PSMs (`identified_peptides` / `identified_proteins`).

## 1.3 Panel 1 — spectrum (2D⇄3D, centroid⇄profile)
- ✅ Starts as **non-profile (centroid)** and **non-3D (2D)**.
- ✅ Two **toggle buttons**: one shows **"3D"/"2D"**, one shows **"profile"/"centroid"** — each label shows **what it will switch to** when pressed. ✅ Buttons are **fixed size/location** so they don't move when the text changes.
- ✅ So Panel 1 is either 2D or 3D, and either profile or centroid.
- **2D view**: ✅ m/z on **x-axis**, intensity on **y-axis**; ⬜ different lines/distributions **colored differently** (from the sqlite distributions).
  - ✅ Shows **every datapoint** (not averages). Both centroid and **profile draw as dots** (m/z vs intensity) — profile is just denser. (Earlier profile-as-per-scan-curves was a regression the user rejected; reverted to dots.)
  - ✅ **Points no longer disappear when zooming into the 2D plot** — `clipToView` + 'peak' auto-downsampling were culling the scatter on zoom-in; both disabled (the window is bounded so the point count stays manageable).
  - ✅ **2D datapoints are coloured by their panel-2 distribution** (one scatter per distribution, grey for points in no distribution shown only when noise is on) — matching panel 2.
  - ✅ **Panel 1 is now filtered to panel 2's visible window** (m/z AND RT). Panel 1 collapses RT, so when panel 2 was zoomed to a narrow RT band, panel 1 still showed points from the whole loaded RT range (lots of colours) while panel 2 showed one — that mismatch is fixed; panel 1 (and the 3D) re-filter on every zoom/pan (debounced). Panel 1 dots **inverse-scale on zoom** too (rescaled after every redraw).
  - ⬜ Profile dots/curves should **scale dynamically on zoom** so they stay visible (currently can get sparse/laggy).
  - ✅ Only the **m/z (x) axis** is interactive when 2D — **no vertical drag**; y auto-scales.
  - ✅ Wheel **over the y-axis label strip** (left of the axis) scrolls **intensity (y)**; ✅ pinned so the **baseline stays at 0** (data is always > 0; the baseline never lifts).
- **3D view**: m/z, time, intensity.
  - ⬜ **Grid lines on the m/z and time axes** marking where each scan was; ⬜ measured datapoints shown as **spherical points** along the grid at their intensity height.
  - ⬜ Points **colored by intensity** via a **gradient**; ⬜ the **3D peaks are a continuous surface** built from the **area between the 3D datapoints** (the surface legitimately represents each time×mass datapoint), with the **datapoints keeping their own point color** on top.
  - 🟡 surface renders; ⬜ make it the true point-to-point area-fill, and ⬜ profile points must **align with the surface** (was misaligned).
  - **3D placement (front-on, 2D-like):** ✅ pyqtgraph's `fov` is the **horizontal** FOV (this was the alignment bug — the data went off-screen). Now m/z → x in **[−1, 1]** fills the pane width with `distance = 1/tan(fov/2)`, and intensity is scaled to the pane's height/width so it fills vertically with the **0 baseline at the bottom**. A **left spacer = the MS2 strip width** makes the 3D content start at the same screen x as panel 2's plot and panel 1's 2D plot. ✅ The **m/z/time text labels were removed** (per the user — they looked bad and weren't helping). ✅ surface removed; just datapoints. ✅ **Spawn and "align 3D" both load a FRONT-ON, near-orthographic view that looks exactly like the 2D plot**: m/z on the horizontal axis (aligned with panel 2), intensity vertical with the **0-intensity baseline pinned to the BOTTOM** of the pane (intensity → GL-z mapped to [−1,+1]); time → GL-y (depth). Camera elevation 0, azimuth −90, small FOV (≈4°) for near-orthographic so m/z is linear and lines up with panel 2; from here you orbit up to reveal the time dimension, and **align 3D returns to exactly this**. ✅ **m/z and time labels are GL text pinned to the ends of their data axes** (they move/rotate WITH the scene, so they always mark the correct axis); intensity is deliberately unlabelled. (The earlier static side-QLabels were removed.) ✅ **Noise toggle now applies to the 3D too** (drops unassigned points when off). ✅ Normalised to panel 2's **visible** window and re-rendered as panel 2 moves; footprint aspect-scaled to fill the pane. ✅ Datapoints coloured by the panel-2 distribution colour, white at the peak tips, log/linear via the colour toggle. 🟡 lag reduced (5k cap, surface gone).
  - ✅ In 3D, moving/rotating **does not move the m/z or time axes** — it only changes the 3D perspective. ✅ When the perspective is rearranged and Panel 2 moves, the perspective stays the same even though the data/axes change.

## 1.4 Panel 2 — m/z × time map
- ✅ The **mass and time window of Panel 1 = the exact region shown in Panel 2.** Panel 2 is zoomable/movable and **shifts everything in Panel 1**.
- ✅ Axes: m/z on **x** (swapped to align with Panel 1's x), time on **y**. (Original ask was m/z on y / time on x; swapped per later request so Panel 1↔2 x-axes align — ✅ now aligned.)
- ✅ Moving Panel 1 (in 2D) pans horizontally on the **mass axis** and moves Panel 2's mass axis too (synchronized). ✅ Moving Panel 2 realigns Panel 1's (2D) mass axis.
- ✅ Zoom/drag on Panel 2 (both axes + scroll-to-zoom) reloads the data only when the view **leaves the cached region**.
- ✅ Render as **connect-the-dots**: individual datapoints as **dots** for time vs mass, **plus thin connecting lines**; the **line is thinner than the dots and the same color**; **each distribution a separate color** (colors from the sqlite). Points not in any distribution = faint gray. **Bug fixed:** the connecting line now follows **each individual line (feature) along its own RT-sorted trace** — it no longer jumps across the different lines of a distribution (NaN-separated polyline per distribution). Consolidated to one curve + one scatter **per distribution** (was per feature) for speed.
- ✅ **Data no longer disappears on zoom** (was the BIGGEST bug). Root cause: every zoom re-extracted only the *visible* window, and an RT view narrower than the MS1 scan spacing returned **zero scans** → blank. Fix: extract a **padded region** (`_padded`, m/z ×2, RT ×2.2) and **cache its extent** (`_loaded_window`); zoom/pan **within** the cache is now a pure view operation (no re-extraction, never empty), reloading only when the view leaves the cached region.
- ✅ Zooming in **inverse-scales the datapoint size** (`_rescale_points`, √ of the cached-vs-view span, clamped 1–4×) so dots stay visible as you zoom in, in both Panel 1 (2D) and Panel 2.
- ✅ Clicking a **distribution's dots in Panel 2 selects it** and brings up the matching **MS1 Panel 3** (charge grid / isotope overlay) for that distribution.
- ✅ **Noise on/off toggle** (fixed-size switch on the Panel 1 bar, like 2D/3D & profile/centroid): "noise off" drops all points **not in any distribution** from Panel 1 and Panel 2 (redrawn from cache, no re-extraction). **Default is noise OFF.**
- ✅ **Distribution colours use `distinctipy`** (a pool of 48 visually-distinct colours, black/white excluded since white is the 3D peak-tip colour), shared across Panel 1, Panel 2 and the Panel 3 MS1 grid; falls back to the fixed palette if `distinctipy` isn't installed (`pip install distinctipy`).
- ✅ **Dot/line sizing** follows the user's matplotlib reference (tiny dots `s≈0.02`, line width `≈0.2`): Panel 2 distribution dots base **3 px** / connecting line **0.5 px**, gray noise dots **1.5 px**, Panel 1 dots **3 px** — all still inverse-scaled on zoom-in.
- ✅ Removed pyqtgraph's in-plot **auto-range "A" button** from every panel (Panel 1/2/3, MS2 strip, charge-grid cells) — fit-to-data made no sense for the window-driven panels and the buttons overlapped the data.

## 1.5 MS2 strip (left of Panel 2)
- ✅ A **tall thin plot just to the left of Panel 2**. It shares Panel 2's **time (y) axis** and is **the sole RT ruler** for the row: Panel 2's own left axis is value-less, so the two RT axes **can never overlap** and this left strip is always visible.
- ✅ MS2 scans shown as **horizontal lines that align with the time axis**, clickable. **These are the only MS2 trigger markers — they live ONLY on this left strip, never inside the Panel 2 plot.** Lines are **green** when the scan has a PSM passing the FDR acceptance criteria and **red** otherwise; the FDR % is an editable "acceptance criteria" field on the panel-1 bar (default 0.1%).
- ✅ The strip **only shows MS2 scans visible within Panel 2's current view** — RT *and* precursor-m/z both inside the view (`_refresh_ms2_visible`, updated on every zoom/pan). This fixes the "tons of MS2 lines" (it was showing every scan in the padded RT range regardless of precursor m/z).
- ✅ **Hovering** an MS2 line draws a **yellow/orange isolation band on Panel 2** (replacing the earlier star): a 3 px horizontal band at the scan's RT (same width as the left strip's RT bands) spanning the **exact isolation m/z range** the MS2 selected for (read from the mzML), in the strip orange at ~50% opacity.
- ✅ The MS2 lines' **thickness inverse-scales with RT zoom** (thicker as you zoom in, fixed minimum) so they never fade to nothing.
- ✅ **Bug fixed:** the MS2 RT lines are now **solid horizontal lines** (one NaN-separated `PlotCurveItem`, fixed **3 px** width) spanning the strip at each RT. Fixed-pixel width keeps them **always visible and a consistent size at any zoom** (they never collapse to dots or vanish); zooming in just **spreads them apart** so individual scans become distinguishable. (The earlier data-space `LinearRegionItem` bands resized inconsistently across reloads — replaced.) Selection = **click anywhere on the strip → nearest line by RT**.

## 1.6 Distribution & line selection + coloring
- Distributions and their member **lines** come from the **sqlite in `<project>/distributions/`**
  (the per-file `*.distributions.sqlite`); a "line" = a feature/isotope trace, a "distribution"
  = its grouped members.
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
- 🟡 When Panel 1 is in **profile**, the **peak-finding (`axis_peaks`, the same centroiding peak-detection from `peaks.py`) runs per scan**; each peak's **apex is the centroid**, matched to a sqlite feature, and **every profile datapoint under that peak (left..right) is assigned to that feature's distribution** — so profile points are linked back to the centroid they'd reduce to and **coloured to match the lines/distributions** (in Panel 1 2D, Panel 2, and the 3D), instead of all reading as noise (`_assignment_profile`).
- 🟡 Profile points under a matched peak are accounted for and coloured by the centroid's distribution; ⬜ still to: surface this in **Table 1** (profile-mode metrics), run it **off the UI thread** for big windows, and tie the peak boundaries to the actual `centroid-mzml.py` output rather than re-detecting.

## 1.10 Panel 3 — MS1 view (charge-comparison grid)
- 🟡 When an **MS1 distribution** is clicked, Panel 3 shows the plot from `examples/panel-3-plot.py`: the **charge-state-determination grid** that **links multiple distributions (charge states of one analyte) into one plot of many comparisons** — columns = charge states, **8 rows**: retention time / peak area / charge distances / cross-charge / intensity sum % / adjacency / ppm-to-mean / ppm-error.
  - ✅ Each **row shares one y-axis across all charge columns** (y-linked + aligned to the left-most cell, which is the only one showing values), so the single left y-axis represents the whole row; columns stay x-linked per charge.
  - ✅ The **retention-time row now plots for every charge** (incl. the left-most): each charge's raw points are read from the store over that charge's own m/z×RT span (`_features_points`, cached) instead of relying on the panel-1 window, which only covered the selected charge.
  - ✅ Each row's y-axis **fits the UNION of all its columns' data** and each column's x-axis fits its own m/z, both set **synchronously on build** (no deferred pass) — so the grid **opens directly in the right state with no flicker** through a wrong auto-fit, and double-click reset is deterministic (same fit). The ppm rows now show their full +/- range.
  - ✅ **Y axes show exactly three ticks** (bottom / middle / top of the live range, via `ThreeTickAxis`) so the rows aren't cluttered; they update dynamically on zoom.
  - Faithful check (`panel-3-plot.py`): `sharey='row'` + `sharex='col'` confirmed; the **cross-charge row's offset/overlapping bars are faithful** — the reference plots a bar per other-charge ratio offset by `diffgen*nc`, so the columns intentionally sit side-by-side. (Inner-column y values are hidden to reduce the clutter.)
  - ✅ Colours now match the reference: single-series rows (peak area, charge distances, intensity-sum%, adjacency, ppm-to-mean) use **that charge's own colour**; the **cross-charge and ppm-error rows colour each bar by the OTHER charge it compares to** (ref `cols[nc]`), with the bars sub-divided side-by-side so each comparison is distinct. **All bars/dots are fully opaque** and use a **uniform bar width** per row.
  - ✅ Peak area / cross-charge / intensity-sum% are **linear bars glued to a 0 baseline**, with the **bottom locked at 0** (`setLimits(yMin=0)`) so dragging can't separate it. Scrolling **over the y-axis** zooms intensity with the baseline pinned (it "goes down" rather than zooming symmetrically); scrolling **over the plot** zooms the m/z (x) axis only (cells are `setMouseEnabled(x=True, y=False)`, wheel-over-axis handled in `eventFilter`).
  - ✅ **Every cell now shows exactly 3 y ticks** (bottom / middle / top), set explicitly and updated live on range change (the previous custom axis wasn't rendering them reliably).
  - ✅ **All charge columns are equal width** (equal column stretch + fixed axis width), so no column's plots are wider than another's.
  - 🟡 grid renders; ⬜ make it **faithful to `panel-3-plot.py`** (per-row scales — peak-area log, charge-distance ylim, cross-charge log, intensity-sum% log, adjacency symlog; the white-on-gray styling adapted to theme; spine hiding; per-column titles `z(distid)`); ⬜ fix colors/readability.
- ⬜ Use **`descending_partial_products` (`libraryadditions.py` / `isotopes.py`)** to compute the **expected isotopic distribution of the peptide** if a peptide is being viewed and found via the search, and **twin-plot** it with the experimental on **different x-axes**, scaled so the **theoretically-most-abundant isotope and the max experimentally-measured datapoint are at the same height**.
- ✅ All MS1 Panel 3 plots that **share the mass x-axis are synchronized** when moved; ✅ **double-click resets** them all.
- ✅ **Panel 3 (single-plot) zoom** matches Panel 1: dragging/scrolling inside the plot zooms the **m/z (x) axis only**; wheel over the **y-axis strip zooms intensity** with the **baseline pinned at 0**; **double-click resets** the zoom.
- ✅ The MS1 Panel 3 now takes the **full Panel 3 + Table 2 space**: Table 2 is **hidden unless Panel 3 is in MS2 mode** (see 1.11 / Table 2).
- ✅ Charge-grid **y-axis labelling cleaned up**: only the **left-most column** shows the row-name label + y tick values (inner columns hide their y values); axis text is **white**; **SI/scientific-notation prefixes disabled**; left axis widened so numbers don't overlap the row label.
- ✅ Removed the leftover **"Panel 3"** title text (both the in-panel caption and the dock title).

## 1.11 Panel 3 — MS2 view
- ✅ **Bug fixed:** clicking an MS2 point now reliably switches Panel 3 to the MS2 spectrum and a `_panel3_mode` flag **keeps it on MS2** across background reloads (was snapping back to MS1).
- ✅ Panel 3 MS2 is triggered by the user **clicking a sampled MS2 point** on the **left MS2 RT strip** (1.5). The MS2 spectrum is **grounded at y=0** (baseline pinned to the bottom of the axis, no gap).
- ⬜ **MS2 trigger markers live on the left MS2 RT strip only** — per the user, they must **not** sit inside the Panel 2 plot (an earlier Panel-2 red-triangle overlay was wrong and has been removed). Overlaying clickable MS2 points on the 2D/3D Panel 1 is still open but, if added, must follow this same "distinct start-point" rule without cluttering the data plots.
- ⬜ When an **identified peptide is assigned to that MS2 spectrum OR to the distribution sampled during that MS2 scan** (link the two via the search info if not already linked), **visualize the theoretical distribution of that specific MS1 progenitor**.
- 🟡 Label the **MS1 progenitor isotopic distribution** and its **fragment isotopic compositions**, labeling the **ions by both type and number** (use `examples/peptidefragmentscoring.py`). ✅ Selecting a candidate in **Table 2** annotates the **actual MS2 spectrum**: ported `fragment_element_binomial_walk` / `fragment_descending_partial_products` + `nearest_neighbors_ppm_tolerance` into `viewer/fragments.py` (dividingthreshold 0.1, subisotopomericdepth 0.5). Only the **precursor isotopes inside the MS2 isolation window** (from the mzML) seed the fragment-isotope walk; the theoretical b/y ions are matched to the **real peaks** via nearest-neighbour at the **search fragment ppm (20)**. Matched **real** peaks plot **green** (labelled ion + the MS1 `M+k` isotope it came from, theme fg); the remaining real peaks plot **red**; unmatched theoretical ions are **not** drawn. The **Table 2 coverage column uses the SAME generate-and-match** (fragments → annotate → matchcounts), so the value and the green ions always agree. Fragment generation always includes the **monoisotopic** ions (the isolation window only ADDS co-isolated M+1/M+2), so a candidate is never spuriously zeroed. Selecting a candidate also **auto-selects its MS1 distribution** in panel 2 (dotted border) without leaving the MS2 view, and the **isolation band stays** on panel 2 (persistent, not hover-only). ✅ Selecting a candidate also overlays that peptide's **theoretical MS1 isotope distribution on panel 1 (2D only)** — a 50%-opacity bar chart at the theoretical m/z, height-normalized so the tallest theoretical bar matches the tallest experimental peak, with a **raw ⇄ summed-M+N** switch (`viewer/isotopes.peptide_isotope_bars`, learned from `libraryadditions.py`). ⬜ the MS1 progenitor envelope twin-plot for the *charge-grid* panel-3 path is still pending (1.10).
- ⬜ Below the MS2 plot: **sequence coverage** of the peptide (the coverage/divider-string logic —
  `coverage_print`-style), per **`examples/sequencecoverageconcept.py`** (reference to be added by the
  user). This is the **MS2 fragment coverage**, distinct from Tab 2's protein coverage.
- ⬜ **Table 2** (only appears below Panel 3 when MS2): the **other peptides this MS1 distribution
  could have matched to** (the other candidate PSMs) **with their relevant sequence coverage**, so
  the user can judge whether one peptide is distinguishable from another. Optional panel. 🟡 dock +
  candidate listing exists, and **Table 2 now only appears when Panel 3 is in MS2 mode** (hidden for
  MS1 so the MS1 view takes the full space); ⬜ "could-have-matched" candidate logic; ✅ the **coverage
  column + score** per **`examples/sequencecoverageconcept.py`** — ported to `viewer/coverage.py`: each
  candidate peptide's theoretical b/y fragment ions are matched (±ppm, charges 1–2) against the MS2
  spectrum, rendered as the divider string (`PEP|T|I|DER`) with matched-ion count + residue coverage,
  and scored by the reference's **`secondfinalmetrics` = `matchcounts`** (unique matched ions + Σ
  segment-length·cover-count) — higher = more/longer/more-redundant coverage = more confident. Table 2
  columns are **peptide / q / coverage**, where coverage IS the matchcounts value, ranked highest-first
  (header-click re-sorts). The MS2 fragment coverage concept, **distinct from Tab 2's protein
  coverage**. ⬜ modified residues are matched on their plain (unmodified) mass for now — a modified
  identified peptide can therefore read 0 until mod mass shifts are applied.

## 1.12 Top-bar controls
- ✅ ± m/z and ± RT window controls (live).
- ✅ charge ◀/▶ + history ⟲/⟳; ✅ theme toggle. (Removed the top-bar **Reset zoom** and the **duplicate Align-3D** button — the 3D align button lives above the 3D plot.)
- 🟡 **Color settings**: ✅ a **log/linear colour-scale switch** for the 3D intensity colouring (fixed-size toggle on the Panel 1 bar, same style as 2D/3D & profile/centroid); ✅ a **raw ⇄ summed-M+N switch** for the theoretical MS1 overlay (same FlipButton style); ⬜ full color-settings drop-down still to do (selected/other-distribution colours, 3D gradient min/max pickers).
- ⬜ Manual **charge** entry field; ✅ **"align/reset 3D"** button (top bar + Panel 1 toolbar).
- 🟡 **"loading… <context>"** label (rendered in **black**, not the old amber) now shown above **Panel 1, Panel 2, and Panel 3** while the evidence worker runs, cleared when drawn (`_set_loading`); ⬜ Panel 3 MS2 / charge-grid sub-loads not yet separately labelled.

---

# TAB 2 — PROTEIN VIEWING 🟡
- ✅ Show **entire protein sequences** (from the project FASTA — reorganized tables carry none); each **tryptic peptide the search attempted** is an outlined rectangle (protease cut segments), with the letters written inside. Segments outside the search length bounds carry no rectangle (plain letters = "not searched").
- ✅ Rectangle **background represents q-value** on a green(best)→red(worst) gradient (`proteins_tab.q_color`); attempted-but-unidentified peptides stay uncoloured. *(Not yet wired to the shared 3D-points color-settings dropdown — see 0.2; standalone gradient for now.)*
- ✅ Same as MS Data: **file selector** + a scrollable **FDR spin box** (near the left, default = search `q_max`) governing the protein list, and a single **protein list**.
- ✅ **Panel 1** (horizontal, wrapping, click-drag to scroll) over **panel 2** (same protein **verticalized side-by-side across every file**, N→C top-to-bottom, centred, file names on a 45° slant; click-drag scrolls, click a column to load that file). Panel 1's bar holds the **All** button, the **FDR** spin box, and a **Light/Dark** toggle (right of the FDR); the **colour-bar legend** (q-value/FDR → colour) floats inside panel 1's top-right.
- ✅ **Double-click a peptide in panel 1** → jump to the MS Data tab reconstructed on that identification: panel 3 MS2 spectrum, panel 1 theoretical-MS1 overlay, panel 2 distribution selected + MS2 isolation band. In All mode it targets the best file for that peptide.
- ✅ **Every** tryptic peptide is boxed in both panels (searched → solid outline, outside the search length bounds → dotted); identified peptides are filled by q-value, the rest left unfilled.
- ✅ Panel 2 **zooms** (Ctrl-wheel) all the way down to fit the whole protein — residues collapse to filled colour blocks (FDR colour preserved) while **file names stay full size**.
- ✅ **Lin/Log colour** toggle (like MS Data's) rescales the FDR colour gradient; the colour bar's percentage ticks move with it.
- ✅ **Light/dark** adaptive colours for the panels only (the left file list keeps the app palette); theme stays in sync across tabs.
- ✅ Panel sizes are **remembered** (h/v splitters persisted by the main window, like the MS Data dock layout); opening/reloading a folder keeps the current tab.
- ✅ Protein **table** (sorted value | protein name) with a "Sort By" row (metric dropdown + asc/desc toggle): % Coverage, Protein Length, Total Identified Peptides, FDR, Spectral Count (PSMs). Changing the sort jumps the view to the top.

---

# TAB 3 — QUANTITATIVE COMPARISONS 🟡
- ✅ Quantitative data **comparing peptides and proteins across files** — the **quant work**. `quant_model.py` builds a feature×file quantity matrix at the **peptide** level (charge/mod variants summed on `peptide_plain`, unique flag per peptide) or **protein** level (**unique quantification** — unique peptides only). A **Peptides⇄Proteins** switch + **unique-only** filter over the bottom table.
- ✅ **Fully generic over the design — NO hard-coded column names.** The first `experimental-setup` column is the filename; **every** other column is a generic category. The only special designation is optional and **user-made**: mark one column as the **replicate** column (its runs are averaged; every other column is **compared, never averaged**). Works for any setup file regardless of which columns exist.
- ✅ **Reactive — no run button.** Changing a category, the replicate column, the contrast, the level, or the unique filter recomputes immediately.
- ✅ **All-features fold-change scatter** (top-right): x = **log2 fold change** (a plain log difference between two chosen category values, **no statistical test**), y = mean log2 abundance; every peptide/protein is a point; click one to select it.
- ✅ **Faceted per-feature view** (top-left) driven by a **growable organizer pseudo-table** (no fixed number of levels): each row picks a category + a mode (**split → columns / split ↓ rows / x-axis**); add layers in any order to nest the panels by depth (split by condition → time-series x-axis inside each). Replicate runs at the same x get a mean line.
- ✅ **Long-form table**: columns are **the experimental-setup categories** (every design column **except the filename** — each file is described by its category values) plus the feature id, a **unique/non-unique** flag (filterable to uniques only), and the quantity — one row per (feature, file), sortable.
- ✅ **Robust filename join**: search-table filenames (e.g. `…​.centroid.mzML`) are matched to the design's `filename` column by stripped stem, so grouping/fold-change work even when the extensions differ (this was silently collapsing every file into one empty group).
- ✅ **Optional normalization**: a **median-center** mode shifts each file so its median log2 quantity matches the grand median (standard label-free correction for per-run loading/intensity), applied consistently to the table, facets, and fold-change; **none** by default.
- ✅ **Persisted state** (like the panel layouts, via QSettings): the organizer layers, replicate/normalize/log-Y, compare A/B, level, and unique filter are restored on reopen. The organizer **never auto-fills** — it starts empty (or however the user left it).
- ✅ **Adaptive light/dark** (matches the other tabs): **all marks white on dark** (incl. the replicate mean line), dark-outlined on light; `(x0.0001)` SI axis prefixes disabled.
- ✅ **Slanted x labels**: panel-1 leaf plots draw their category/filename x labels at 45° (`RotatedAxisItem`) when they'd otherwise overlap (heuristic on label count/width), so long filenames stay legible; horizontal when few/short.
- ✅ **Split sub-plots are self-describing + share axes**: each faceted leaf is titled with its full split path (`condition = keloid / fraction = NaCl`), all leaves use one **global x ordering** (same category at the same position), and their **X and Y views are linked to a common range** so heights/positions are directly comparable across the whole grid.
- ✅ **Table default layout**: columns auto-expand to their content (sampled, so 100k+ rows don't stall), cells + headers centred, and spare viewport width is spread across columns (re-applied on resize) so the table fills the width instead of clumping left.
- ⬜ Worker-thread recompute for large peptide sets; **save named layouts / narrowed feature sets**; clustered-heatmap view; per-feature link into the MS Data tab.

---

# TAB 4 — MOTIF QUANTIFICATION ⬜
- Motif index location is **auto-detected, not hard-coded**; expected structure is a
  `human-proteome-skeletons/{build_info.tsv, motifs.tsv, postings.bin, proteins.tsv}` dir under a
  motifs folder (e.g. `~/data/proteomics/motifs/`). Reader exists: `viewer/motifs.py`.
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

**Intent of the rework — derive constants, don't hardcode them.** The *point* of
rewriting the distribution functions is to get rid of the magic numbers
(`step_limit`, `new_inc_limit`, `charge_tolerance`, `mass_width_limit`,
`roundcutoff`, etc.) by making them **data-adaptive / derived from the data** rather
than fixed. The reference already shows the mechanism: `roundcutoff` is a per-scan
moving average of the knee of sorted match distances (init 0), `masswidthlimit =
roundcutoff*2`, and `madiff` is a per-trace moving average — i.e. the thresholds
should fall out of the signal's own spacing/intensity statistics. Prefer that for
every constant we can. **If a given value genuinely can't be derived, keeping it as a
tunable constant (with a sensible default + CLI flag) is an acceptable fallback** —
this is not a hard requirement, just the guiding intent. Note in code which constants
ended up derived vs. left fixed and why.

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
