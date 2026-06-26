# Restoring faithful MS1 distribution detection

## Root cause (corrected)

The **active** pipeline is `index_ms1.py` (imports only `peaks.py` + `store.py`,
writes the sqlite). It is a from-scratch reimplementation. The **reference**
(known-good) algorithm lives in `linemodel.py` (`line_model`),
`distributionassembly.py` (`distribution_assembly`, `overlap_counts`), and
`chargehandling.py` (`charge_handling`) — these are faithful to the user's pasted
reference but are **not imported / not run**.

So the bad distributions are not a few changed constants in the reference files;
they are **missing algorithm stages** in `index_ms1.py`. Decision (user): keep the
sqlite output (the GUI needs it), make the algorithm faithful to the reference.

## Faithful changes to port INTO index_ms1.py (keep sqlite)

Ordered by impact (from the code-to-code comparison). Each must reproduce the
reference exactly; validation is on the user's machine against a known-good file.

### Stage 1 — line model (LineModel / process_scan)
1. **Adaptive `roundcutoff`** (ref linemodel.py:201–209): per-scan moving average
   of the knee cutoff of sorted match distances; init 0. Currently absent.
2. **Moving-avg diff tracking** (`groupdifftoma`): add running `madiff` per trace.
3. **3-tier acceptance** (ref:224–266) replacing the single ppm/abs tolerance:
   (a) within line mass range OR dist < range/2; (b) if n ≥ `minmovinginds`,
   accept if `nmadiff ≤ madiff`; (c) n>1: `d ≤ roundcutoff + range`, n==1:
   `d ≤ roundcutoff*2`.
4. **Dead-signal counter** (ref:302–303) replacing the hard `max_gap_scans`:
   per-line counter, close after `> deadsignal`, halve on match.
5. **Tie-break** (ref:121–156): equidistant masses → prefer the one inside the
   line's range, else closer in intensity.
6. **Line-correction merge** (ref:310–500): after closing, merge sub-`subisomax`
   (0.01337851739·(1+chargetolerance)) lines with non-redundant timepoints via
   intersection_merge + directional-graph time checks. Currently absent.

### Stage 2 — distribution assembly (edge building / charge)
7. **acdiff acceptance** (ref:226–231): `expdiff=proton/charge`, `diffcut=expdiff*
   chargetolerance`, accept if `-(diffcut*chargetolerance+masswidthlimit) < acdiff
   ≤ diffcut+masswidthlimit`; `masswidthlimit=roundcutoff*2`. Currently a flat
   ppm/abs gate.
8. **Intensity-step gating** (ref:240–278): `intensitypercdiff=|ΔI|/(ΣI)/2`,
   gate on `steplimit`=0.5 / `newinclimit`=0.1 (tighten when reversing decrease→
   increase). Currently absent.
9. **3-tier RT overlap** (ref:191–209): encompassed → `>newinclimit`, partial →
   `>0.5`. Currently flat overlap/union.
10. **overlap_counts scoring** (ref:11–117, used at 344) for path RT geometry, and
    the **3-tier pair ranking** (ref:416–437). Currently single best-edge-per-left.
11. **charge ±1 spread refinement** (ref:215–222) — make explicit.

### Stage 3 — charge handling (analyte grouping)
12. **Intensity-rank order gating** (ref:185–188): align base masses
    (`mz*z - proton*z`), require `Σ|rankdiff| ≤ size-1`. Currently absent.
13. **RT-overlap majority gate** (ref): `overpass > matchables/2`. Currently soft score.
14. **Adjacent-charge-only search** (ref:98–99) + active nodist up/down matching
    (ref:200–355). Currently tests all charges.

## Constants to introduce (reference values)
`subisomax=0.01337851739` (·(1+chargetolerance)), `newinclimit=0.1`,
`steplimit=0.5`, `masswidthlimit=roundcutoff*2`, plus params `minpoints`,
`minmovinginds`, `deadsignal`, `chargetolerance` (currently `index_ms1.py` uses
`line_mz_ppm=8`, `line_mz_abs=0.002`, `isotope_mz_ppm=10`, `isotope_mz_abs=0.004`,
`max_gap_scans=2`, `min_trace_points=4`, `min_edge_score=0.30`, `charge_mass_ppm=12`).

## Approach
Port stage by stage, smallest reviewable diffs, preserving the sqlite schema
(`store.py`). After each stage the user runs the pipeline on a real file and
compares distributions to the prior known-good output. Once Python is faithful
and validated, the validated algorithm can be ported to Rust (separate phase).
