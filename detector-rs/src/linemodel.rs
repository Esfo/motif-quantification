//! Port of LineModel in index_ml1.py: streaming tracker -> history-weighted
//! death -> fragment merge -> whole-line emission with 2-sigma peak split.

use std::collections::HashMap;

use crate::config::Config;
use crate::peaks::{axis_peaks, moving_average};

#[derive(Clone)]
pub struct LineRow {
    pub line_id: i64,
    pub mz_mean: f64,
    pub mz_min: f64,
    pub mz_max: f64,
    pub rt_start: f64,
    pub rt_end: f64,
    pub ms1_start: i64,
    pub ms1_end: i64,
    pub n_points: i64,
}

#[derive(Clone)]
pub struct FeatureRow {
    pub feature_id: i64,
    pub line_id: i64,
    pub mz_mean: f64,
    pub mz_min: f64,
    pub mz_max: f64,
    pub rt_start: f64,
    pub rt_apex: f64,
    pub rt_end: f64,
    pub ms1_start: i64,
    pub ms1_apex: i64,
    pub ms1_end: i64,
    pub height: f64,
    pub area: f64,
    pub n_points: i64,
    pub quality: f64,
}

struct Trace {
    scans: Vec<i64>,
    rts: Vec<f64>,
    mzs: Vec<f64>,
    intensities: Vec<f64>,
    mean_mz: f64,
    last_scan: i64,
}

impl Trace {
    fn new(ms1_index: i64, rt: f64, mz: f64, intensity: f64) -> Self {
        Trace {
            scans: vec![ms1_index],
            rts: vec![rt],
            mzs: vec![mz],
            intensities: vec![intensity],
            mean_mz: mz,
            last_scan: ms1_index,
        }
    }
    fn append(&mut self, ms1_index: i64, rt: f64, mz: f64, intensity: f64) {
        let n = self.mzs.len() as f64;
        self.mean_mz = (self.mean_mz * n + mz) / (n + 1.0);
        self.last_scan = ms1_index;
        self.scans.push(ms1_index);
        self.rts.push(rt);
        self.mzs.push(mz);
        self.intensities.push(intensity);
    }
}

struct Fragment {
    scans: Vec<i64>,
    rts: Vec<f64>,
    mzs: Vec<f64>,
    intensities: Vec<f64>,
    mean_mz: f64,
    start: i64,
    end: i64,
}

struct Chain {
    mean_mz: f64,
    end: i64,
    n: usize,
    frags: Vec<usize>, // indices into fragments
}

pub struct LineModel {
    cfg: Config,
    active: HashMap<i64, Trace>,
    misses: HashMap<i64, i64>,
    closed: Vec<Trace>,
    next_line_id: i64,
    next_feature_id: i64,
    pub lines: Vec<LineRow>,
    pub features: Vec<FeatureRow>,
}

impl LineModel {
    pub fn new(cfg: Config) -> Self {
        LineModel {
            cfg,
            active: HashMap::new(),
            misses: HashMap::new(),
            closed: Vec::new(),
            next_line_id: 0,
            next_feature_id: 0,
            lines: Vec::new(),
            features: Vec::new(),
        }
    }

    fn mz_tol(&self, mz: f64) -> f64 {
        self.cfg.line_mz_abs.max(mz * self.cfg.line_mz_ppm / 1_000_000.0)
    }

    pub fn process_scan(&mut self, ms1_index: i64, rt: f64, mzs: &[f64], intensities: &[f64]) {
        if mzs.is_empty() {
            let ids: Vec<i64> = self.active.keys().copied().collect();
            for id in ids {
                *self.misses.entry(id).or_insert(0) += 1;
            }
            self.close_dead_lines();
            return;
        }

        // sort centroids by m/z ascending
        let mut order: Vec<usize> = (0..mzs.len()).collect();
        order.sort_by(|&a, &b| mzs[a].partial_cmp(&mzs[b]).unwrap());

        // active sorted by mean_mz (stable by id for ties)
        let mut active_sorted: Vec<(f64, i64)> =
            self.active.iter().map(|(&id, t)| (t.mean_mz, id)).collect();
        active_sorted.sort_by(|a, b| {
            a.0.partial_cmp(&b.0).unwrap().then(a.1.cmp(&b.1))
        });
        let active_mz: Vec<f64> = active_sorted.iter().map(|x| x.0).collect();

        let active_at_start: Vec<i64> = active_sorted.iter().map(|x| x.1).collect();
        let mut matched_existing: std::collections::HashSet<i64> = std::collections::HashSet::new();
        let mut used_lines: std::collections::HashSet<i64> = std::collections::HashSet::new();

        for &oi in &order {
            let mz = mzs[oi];
            let intensity = intensities[oi];
            let tol = self.mz_tol(mz);

            // searchsorted [mz-tol, mz+tol)
            let left = lower_bound(&active_mz, mz - tol);
            let right = upper_bound(&active_mz, mz + tol);

            let mut best_id: Option<i64> = None;
            let mut best_score = f64::INFINITY;
            for loc in left..right {
                let id = active_sorted[loc].1;
                if used_lines.contains(&id) {
                    continue;
                }
                let trace = &self.active[&id];
                let gap = ms1_index - trace.last_scan;
                if gap < 1 || gap > self.cfg.deadsignal + 1 {
                    continue;
                }
                let score = (trace.mean_mz - mz).abs() / tol + gap as f64 * 0.03;
                if score < best_score {
                    best_score = score;
                    best_id = Some(id);
                }
            }

            match best_id {
                None => {
                    let id = self.next_line_id;
                    self.next_line_id += 1;
                    self.active.insert(id, Trace::new(ms1_index, rt, mz, intensity));
                    self.misses.insert(id, 0);
                    used_lines.insert(id);
                }
                Some(id) => {
                    self.active.get_mut(&id).unwrap().append(ms1_index, rt, mz, intensity);
                    matched_existing.insert(id);
                    used_lines.insert(id);
                }
            }
        }

        for id in active_at_start {
            if matched_existing.contains(&id) {
                let m = self.misses.entry(id).or_insert(0);
                *m /= 2;
            } else {
                *self.misses.entry(id).or_insert(0) += 1;
            }
        }

        self.close_dead_lines();
    }

    fn close_dead_lines(&mut self) {
        let dead: Vec<i64> = self
            .active
            .iter()
            .filter(|(id, _)| *self.misses.get(id).unwrap_or(&0) > self.cfg.deadsignal)
            .map(|(&id, _)| id)
            .collect();
        for id in dead {
            let t = self.active.remove(&id).unwrap();
            self.misses.remove(&id);
            self.closed.push(t);
        }
    }

    pub fn finalize(&mut self) {
        // flush active (sorted by id, like Python's sorted(self.active))
        let mut ids: Vec<i64> = self.active.keys().copied().collect();
        ids.sort_unstable();
        for id in ids {
            let t = self.active.remove(&id).unwrap();
            self.closed.push(t);
        }
        self.misses.clear();
        self.merge_and_emit();
    }

    fn merge_and_emit(&mut self) {
        let mut fragments: Vec<Fragment> = Vec::with_capacity(self.closed.len());
        for t in self.closed.drain(..) {
            if t.scans.is_empty() {
                continue;
            }
            // sort by scan
            let mut ord: Vec<usize> = (0..t.scans.len()).collect();
            ord.sort_by(|&a, &b| t.scans[a].cmp(&t.scans[b]));
            let scans: Vec<i64> = ord.iter().map(|&i| t.scans[i]).collect();
            let rts: Vec<f64> = ord.iter().map(|&i| t.rts[i]).collect();
            let mzs: Vec<f64> = ord.iter().map(|&i| t.mzs[i]).collect();
            let intensities: Vec<f64> = ord.iter().map(|&i| t.intensities[i]).collect();
            let start = scans[0];
            let end = *scans.last().unwrap();
            fragments.push(Fragment {
                mean_mz: t.mean_mz,
                start,
                end,
                scans,
                rts,
                mzs,
                intensities,
            });
        }

        // sort fragments by (start, mean_mz)
        let mut frag_order: Vec<usize> = (0..fragments.len()).collect();
        frag_order.sort_by(|&a, &b| {
            fragments[a]
                .start
                .cmp(&fragments[b].start)
                .then(fragments[a].mean_mz.partial_cmp(&fragments[b].mean_mz).unwrap())
        });

        let chains = self.chain_fragments(&fragments, &frag_order);

        // materialise lines (concat each chain's fragments, sort by scan)
        let mut materialised: Vec<(Vec<i64>, Vec<f64>, Vec<f64>, Vec<f64>)> = Vec::new();
        for chain in &chains {
            let mut scans = vec![];
            let mut rts = vec![];
            let mut mzs = vec![];
            let mut ints = vec![];
            for &fi in &chain.frags {
                let f = &fragments[fi];
                scans.extend_from_slice(&f.scans);
                rts.extend_from_slice(&f.rts);
                mzs.extend_from_slice(&f.mzs);
                ints.extend_from_slice(&f.intensities);
            }
            let mut ord: Vec<usize> = (0..scans.len()).collect();
            ord.sort_by(|&a, &b| scans[a].cmp(&scans[b]));
            let scans: Vec<i64> = ord.iter().map(|&i| scans[i]).collect();
            let rts: Vec<f64> = ord.iter().map(|&i| rts[i]).collect();
            let mzs: Vec<f64> = ord.iter().map(|&i| mzs[i]).collect();
            let ints: Vec<f64> = ord.iter().map(|&i| ints[i]).collect();
            if scans.len() < self.cfg.min_trace_points {
                continue;
            }
            materialised.push((scans, rts, mzs, ints));
        }

        // 2-sigma split threshold over line RT spans
        let split_threshold = if materialised.is_empty() {
            f64::INFINITY
        } else {
            let spans: Vec<f64> = materialised
                .iter()
                .map(|(_, rts, _, _)| rts[rts.len() - 1] - rts[0])
                .collect();
            let mean = spans.iter().sum::<f64>() / spans.len() as f64;
            let var = spans.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / spans.len() as f64;
            mean + self.cfg.line_split_sigma * var.sqrt()
        };

        for (line_id, (scans, rts, mzs, ints)) in materialised.into_iter().enumerate() {
            let span = rts[rts.len() - 1] - rts[0];
            self.emit_line(line_id as i64, &scans, &rts, &mzs, &ints, span > split_threshold);
        }
    }

    fn chain_fragments(&self, fragments: &[Fragment], frag_order: &[usize]) -> Vec<Chain> {
        let gap = self.cfg.line_merge_gap_scans;
        let bin_width = 0.05f64;
        let mut active_by_bin: HashMap<i64, Vec<usize>> = HashMap::new(); // bin -> chain idx
        let mut chains: Vec<Chain> = Vec::new();
        let bin_of = |mz: f64| (mz / bin_width).floor() as i64;

        for &fi in frag_order {
            let f_start = fragments[fi].start;
            let f_mz = fragments[fi].mean_mz;
            let f_end = fragments[fi].end;
            let f_n = fragments[fi].scans.len();
            let tol = self
                .cfg
                .line_merge_mz_abs
                .max(f_mz * self.cfg.line_merge_mz_ppm / 1_000_000.0);
            let home = bin_of(f_mz);

            let mut best: Option<usize> = None;
            let mut best_d = f64::INFINITY;

            for b in [home - 1, home, home + 1] {
                if let Some(bucket) = active_by_bin.get_mut(&b) {
                    let mut kept = Vec::with_capacity(bucket.len());
                    for &ci in bucket.iter() {
                        if f_start - chains[ci].end > gap {
                            continue; // stale, drop
                        }
                        kept.push(ci);
                        if f_start <= chains[ci].end {
                            continue; // overlap -> not a continuation
                        }
                        let d = (chains[ci].mean_mz - f_mz).abs();
                        if d <= tol && d < best_d {
                            best_d = d;
                            best = Some(ci);
                        }
                    }
                    *bucket = kept;
                }
            }

            match best {
                None => {
                    let ci = chains.len();
                    chains.push(Chain {
                        mean_mz: f_mz,
                        end: f_end,
                        n: f_n,
                        frags: vec![fi],
                    });
                    active_by_bin.entry(home).or_default().push(ci);
                }
                Some(ci) => {
                    let old_bin = bin_of(chains[ci].mean_mz);
                    let n0 = chains[ci].n;
                    chains[ci].mean_mz = (chains[ci].mean_mz * n0 as f64 + f_mz * f_n as f64)
                        / (n0 + f_n) as f64;
                    chains[ci].n = n0 + f_n;
                    chains[ci].end = chains[ci].end.max(f_end);
                    chains[ci].frags.push(fi);
                    let new_bin = bin_of(chains[ci].mean_mz);
                    if new_bin != old_bin {
                        if let Some(bucket) = active_by_bin.get_mut(&old_bin) {
                            if let Some(pos) = bucket.iter().position(|&x| x == ci) {
                                bucket.remove(pos);
                            }
                        }
                        active_by_bin.entry(new_bin).or_default().push(ci);
                    }
                }
            }
        }

        chains
    }

    fn emit_line(
        &mut self,
        line_id: i64,
        scans: &[i64],
        rts: &[f64],
        mzs: &[f64],
        ints: &[f64],
        do_split: bool,
    ) {
        let (mz_min, mz_max) = min_max(mzs);
        let (rt_min, rt_max) = min_max(rts);
        self.lines.push(LineRow {
            line_id,
            mz_mean: mean(mzs),
            mz_min,
            mz_max,
            rt_start: rt_min,
            rt_end: rt_max,
            ms1_start: *scans.iter().min().unwrap(),
            ms1_end: *scans.iter().max().unwrap(),
            n_points: scans.len() as i64,
        });

        if do_split {
            self.split_trace(line_id, scans, rts, mzs, ints);
        } else {
            self.emit_whole_feature(line_id, scans, rts, mzs, ints);
        }
    }

    fn emit_whole_feature(&mut self, line_id: i64, scans: &[i64], rts: &[f64], mzs: &[f64], ints: &[f64]) {
        let area = if rts.len() > 1 { trapezoid(ints, rts) } else { ints[0] };
        let height = max(ints);
        let local_apex = argmax(ints);
        let total: f64 = ints.iter().sum();
        let mz_mean = if total > 0.0 {
            mzs.iter().zip(ints).map(|(m, i)| m * i).sum::<f64>() / total
        } else {
            mean(mzs)
        };
        let quality = (area.max(0.0)).ln_1p() * (scans.len() as f64).sqrt();
        let (mz_min, mz_max) = min_max(mzs);
        self.features.push(FeatureRow {
            feature_id: self.next_feature_id,
            line_id,
            mz_mean,
            mz_min,
            mz_max,
            rt_start: min(rts),
            rt_apex: rts[local_apex],
            rt_end: max(rts),
            ms1_start: *scans.iter().min().unwrap(),
            ms1_apex: scans[local_apex],
            ms1_end: *scans.iter().max().unwrap(),
            height,
            area,
            n_points: scans.len() as i64,
            quality,
        });
        self.next_feature_id += 1;
    }

    fn split_trace(&mut self, line_id: i64, scans: &[i64], rts: &[f64], mzs: &[f64], ints: &[f64]) {
        if ints.len() < self.cfg.min_trace_points {
            return;
        }
        let smoothed = moving_average(ints, self.cfg.smooth_points);
        let peaks = axis_peaks(&smoothed, self.cfg.peak_mindist);
        if peaks.is_empty() {
            return;
        }
        let trace_height = max(ints);
        // candidate peaks: (apex_height, left, apex, right)
        let mut candidates: Vec<(f64, usize, usize, usize)> = Vec::new();
        for (left, apex, right) in peaks {
            if right <= left {
                continue;
            }
            let sub_count = right - left;
            if sub_count < self.cfg.min_peak_points {
                continue;
            }
            let apex_height = ints[apex];
            let edge_height = ints[left].max(ints[right - 1]);
            let prominence = apex_height - edge_height;
            if apex_height < self.cfg.min_peak_height {
                continue;
            }
            if self.cfg.min_peak_prominence_fraction > 0.0
                && trace_height > 0.0
                && prominence < trace_height * self.cfg.min_peak_prominence_fraction
            {
                continue;
            }
            candidates.push((apex_height, left, apex, right));
        }
        if candidates.is_empty() {
            return;
        }
        // sort by apex_height desc
        candidates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        if self.cfg.max_trace_peaks > 0 && candidates.len() > self.cfg.max_trace_peaks {
            candidates.truncate(self.cfg.max_trace_peaks);
        }
        // sort by left
        candidates.sort_by(|a, b| a.1.cmp(&b.1));

        // (valley merge skipped: min_split_valley_fraction default 0.0)

        for (_, left, _apex, right) in candidates {
            let sub_scans = &scans[left..right];
            let sub_rts = &rts[left..right];
            let sub_mzs = &mzs[left..right];
            let sub_ints = &ints[left..right];
            if sub_scans.len() < self.cfg.min_peak_points {
                continue;
            }
            let width = max(sub_rts) - min(sub_rts);
            if width < self.cfg.min_peak_width {
                continue;
            }
            if self.cfg.max_peak_width > 0.0 && width > self.cfg.max_peak_width {
                continue;
            }
            let area = if sub_rts.len() > 1 { trapezoid(sub_ints, sub_rts) } else { sub_ints[0] };
            if area < self.cfg.min_peak_area {
                continue;
            }
            let height = max(sub_ints);
            if height < self.cfg.min_peak_height {
                continue;
            }
            let local_apex = argmax(sub_ints);
            let total: f64 = sub_ints.iter().sum();
            let mz_mean = if total > 0.0 {
                sub_mzs.iter().zip(sub_ints).map(|(m, i)| m * i).sum::<f64>() / total
            } else {
                mean(sub_mzs)
            };
            let quality = (area.max(0.0)).ln_1p() * (sub_scans.len() as f64).sqrt();
            let (mz_min, mz_max) = min_max(sub_mzs);
            self.features.push(FeatureRow {
                feature_id: self.next_feature_id,
                line_id,
                mz_mean,
                mz_min,
                mz_max,
                rt_start: min(sub_rts),
                rt_apex: sub_rts[local_apex],
                rt_end: max(sub_rts),
                ms1_start: *sub_scans.iter().min().unwrap(),
                ms1_apex: sub_scans[local_apex],
                ms1_end: *sub_scans.iter().max().unwrap(),
                height,
                area,
                n_points: sub_scans.len() as i64,
                quality,
            });
            self.next_feature_id += 1;
        }
    }
}

fn lower_bound(a: &[f64], x: f64) -> usize {
    // first index with a[i] >= x  (searchsorted side="left")
    let mut lo = 0;
    let mut hi = a.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

fn upper_bound(a: &[f64], x: f64) -> usize {
    // first index with a[i] > x  (searchsorted side="right")
    let mut lo = 0;
    let mut hi = a.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] <= x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

fn trapezoid(y: &[f64], x: &[f64]) -> f64 {
    let mut s = 0.0;
    for i in 1..y.len() {
        s += (x[i] - x[i - 1]) * (y[i] + y[i - 1]) / 2.0;
    }
    s
}

fn mean(a: &[f64]) -> f64 {
    a.iter().sum::<f64>() / a.len() as f64
}
fn max(a: &[f64]) -> f64 {
    a.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}
fn min(a: &[f64]) -> f64 {
    a.iter().copied().fold(f64::INFINITY, f64::min)
}
fn min_max(a: &[f64]) -> (f64, f64) {
    (min(a), max(a))
}
fn argmax(a: &[f64]) -> usize {
    let mut bi = 0;
    for i in 1..a.len() {
        if a[i] > a[bi] {
            bi = i;
        }
    }
    bi
}
