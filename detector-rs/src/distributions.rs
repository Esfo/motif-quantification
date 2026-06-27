//! Port of the edge + distribution stages from index_ml1.py.

#![allow(dead_code)]
use rayon::prelude::*;

use crate::config::{Config, C13_DELTA, PROTON};
use crate::linemodel::{FeatureRow, FeatureTrace};

const AVERAGINE_LAMBDA_PER_DA: f64 = 0.000594;
const MAX_ISOTOPE_GAP: i64 = 2;
const MAX_OFFSET_SEARCH: i64 = 3;
const MAX_Q: i64 = 2;

pub struct Features {
    pub mz: Vec<f64>,
    pub rt_start: Vec<f64>,
    pub rt_apex: Vec<f64>,
    pub rt_end: Vec<f64>,
    pub ms1_start: Vec<i64>,
    pub ms1_apex: Vec<i64>,
    pub ms1_end: Vec<i64>,
    pub height: Vec<f64>,
    pub traces: Vec<FeatureTrace>,
    // m/z-sorted view
    pub order: Vec<usize>,
    pub sorted_mz: Vec<f64>,
}

impl Features {
    pub fn build(rows: &[FeatureRow], traces: &[FeatureTrace]) -> Self {
        let mz: Vec<f64> = rows.iter().map(|f| f.mz_mean).collect();
        let mut order: Vec<usize> = (0..rows.len()).collect();
        order.sort_by(|&a, &b| mz[a].partial_cmp(&mz[b]).unwrap());
        let sorted_mz: Vec<f64> = order.iter().map(|&i| mz[i]).collect();
        Features {
            rt_start: rows.iter().map(|f| f.rt_start).collect(),
            rt_apex: rows.iter().map(|f| f.rt_apex).collect(),
            rt_end: rows.iter().map(|f| f.rt_end).collect(),
            ms1_start: rows.iter().map(|f| f.ms1_start).collect(),
            ms1_apex: rows.iter().map(|f| f.ms1_apex).collect(),
            ms1_end: rows.iter().map(|f| f.ms1_end).collect(),
            height: rows.iter().map(|f| f.height).collect(),
            traces: traces.to_vec(),
            mz,
            order,
            sorted_mz,
        }
    }

    pub fn from_rows(rows: &[FeatureRow]) -> Self {
        Features::build(rows, &[])
    }

    fn mz_tol(&self, cfg: &Config, mz: f64) -> f64 {
        cfg.isotope_mz_abs.max(mz * cfg.isotope_mz_ppm / 1_000_000.0)
    }

    /// Feature ids whose m/z is within tolerance of target.
    fn features_near(&self, cfg: &Config, target: f64) -> Vec<usize> {
        let tol = self.mz_tol(cfg, target);
        let lo = lower_bound(&self.sorted_mz, target - tol);
        let hi = upper_bound_le(&self.sorted_mz, target + tol);
        (lo..hi).map(|k| self.order[k]).collect()
    }

    fn features_in_range(&self, lo_mz: f64, hi_mz: f64) -> Vec<usize> {
        let lo = lower_bound(&self.sorted_mz, lo_mz);
        let hi = upper_bound_le(&self.sorted_mz, hi_mz);
        (lo..hi).map(|k| self.order[k]).collect()
    }

    fn neighbours_above(&self, seed: usize, max_dmz: f64) -> Vec<usize> {
        let seed_mz = self.mz[seed];
        let lo = upper_bound_le(&self.sorted_mz, seed_mz);
        let hi = upper_bound_le(&self.sorted_mz, seed_mz + max_dmz);
        (lo..hi).map(|k| self.order[k]).filter(|&j| j != seed).collect()
    }

    fn coelutes(&self, a: usize, b: usize) -> bool {
        self.rt_end[a].min(self.rt_end[b]) - self.rt_start[a].max(self.rt_start[b]) > 0.0
    }

    /// Cosine of two elution traces aligned on shared scans (0..1).
    fn trace_similarity(&self, a: usize, b: usize) -> f64 {
        if a == b {
            return 1.0;
        }
        let ta = &self.traces[a];
        let tb = &self.traces[b];
        if ta.scans.is_empty() || tb.scans.is_empty() {
            return 0.0;
        }
        let lo = *ta.scans.iter().min().unwrap().min(tb.scans.iter().min().unwrap());
        let hi = *ta.scans.iter().max().unwrap().max(tb.scans.iter().max().unwrap());
        let n = (hi - lo + 1) as usize;
        if n == 0 {
            return 0.0;
        }
        let mut va = vec![0.0f64; n];
        let mut vb = vec![0.0f64; n];
        for (s, i) in ta.scans.iter().zip(&ta.intensities) {
            va[(*s - lo) as usize] = *i;
        }
        for (s, i) in tb.scans.iter().zip(&tb.intensities) {
            vb[(*s - lo) as usize] = *i;
        }
        let na = va.iter().map(|x| x * x).sum::<f64>().sqrt();
        let nb = vb.iter().map(|x| x * x).sum::<f64>().sqrt();
        if na <= 0.0 || nb <= 0.0 {
            return 0.0;
        }
        let dot: f64 = va.iter().zip(&vb).map(|(x, y)| x * y).sum();
        (dot / (na * nb)).clamp(0.0, 1.0)
    }

    /// Consensus elution trace; returns each member's cosine to it + the mean.
    fn consensus_trace_scores(&self, fids: &[usize]) -> (Vec<f64>, f64) {
        let traces: Vec<&FeatureTrace> = fids.iter().map(|&f| &self.traces[f]).collect();
        let lo = traces.iter().map(|t| *t.scans.iter().min().unwrap()).min().unwrap();
        let hi = traces.iter().map(|t| *t.scans.iter().max().unwrap()).max().unwrap();
        let n = (hi - lo + 1) as usize;
        if n == 0 {
            return (vec![1.0; fids.len()], 1.0);
        }
        let mut vecs: Vec<Vec<f64>> = Vec::with_capacity(traces.len());
        let mut consensus = vec![0.0f64; n];
        for t in &traces {
            let mut v = vec![0.0f64; n];
            for (s, i) in t.scans.iter().zip(&t.intensities) {
                v[(*s - lo) as usize] = *i;
            }
            for k in 0..n {
                consensus[k] += v[k];
            }
            vecs.push(v);
        }
        let cn = consensus.iter().map(|x| x * x).sum::<f64>().sqrt();
        if cn <= 0.0 {
            return (vec![1.0; fids.len()], 1.0);
        }
        let mut scores = Vec::with_capacity(vecs.len());
        for v in &vecs {
            let vn = v.iter().map(|x| x * x).sum::<f64>().sqrt();
            let s = if vn > 0.0 {
                let dot: f64 = v.iter().zip(&consensus).map(|(x, y)| x * y).sum();
                (dot / (vn * cn)).clamp(0.0, 1.0)
            } else {
                0.0
            };
            scores.push(s);
        }
        let mean = scores.iter().sum::<f64>() / scores.len() as f64;
        (scores, mean)
    }
}

/// Normalised (max=1) averagine isotope intensities for indices 0..n_peaks-1.
fn averagine_envelope(neutral_mass: f64, n_peaks: usize) -> Vec<f64> {
    let lam = neutral_mass.max(1.0) * AVERAGINE_LAMBDA_PER_DA;
    let mut p = vec![0.0f64; n_peaks];
    let mut max = 0.0f64;
    for k in 0..n_peaks {
        let logp = -lam + k as f64 * lam.ln() - ln_factorial(k);
        let v = logp.exp();
        p[k] = v;
        if v > max {
            max = v;
        }
    }
    if max > 0.0 {
        for v in p.iter_mut() {
            *v /= max;
        }
    }
    p
}

fn ln_factorial(k: usize) -> f64 {
    // lgamma(k+1)
    let mut s = 0.0;
    for i in 2..=k {
        s += (i as f64).ln();
    }
    s
}

#[derive(Clone)]
pub struct Edge {
    pub left: usize,
    pub right: usize,
    pub charge: i64,
    pub score: f64,
}

fn rt_overlap(f: &Features, l: usize, r: usize) -> f64 {
    let overlap = f.rt_end[l].min(f.rt_end[r]) - f.rt_start[l].max(f.rt_start[r]);
    let union = f.rt_end[l].max(f.rt_end[r]) - f.rt_start[l].min(f.rt_start[r]);
    if union <= 0.0 {
        0.0
    } else {
        (overlap / union).max(0.0)
    }
}

fn intensity_score(f: &Features, l: usize, r: usize) -> f64 {
    let lh = f.height[l];
    let rh = f.height[r];
    if lh <= 0.0 || rh <= 0.0 {
        return 0.0;
    }
    let ratio = rh / lh;
    let mut score = 1.0 - (ratio.log2().abs() / 4.0).min(1.0);
    if rh > lh * 4.0 {
        score *= 0.75;
    }
    score.clamp(0.0, 1.0)
}

fn score_edge(f: &Features, cfg: &Config, l: usize, r: usize, charge: i64) -> Option<Edge> {
    let expdiff = PROTON / charge as f64;
    let observed = f.mz[r] - f.mz[l];
    let acdiff = expdiff - observed;
    let diffcut = expdiff * cfg.charge_tolerance;
    let masswidthlimit = cfg.mass_width_limit;
    let lower = -(diffcut * cfg.charge_tolerance + masswidthlimit);
    let upper = diffcut + masswidthlimit;
    if !(acdiff > lower && acdiff <= upper) {
        return None;
    }
    let mz_error = -acdiff;
    let tolerance = if upper > 0.0 { upper } else { cfg.isotope_mz_abs.max(1e-9) };

    let neutral_mass = (f.mz[l] - PROTON) * charge as f64;
    if neutral_mass <= 0.0 || neutral_mass > cfg.max_neutral_mass {
        return None;
    }

    let lh = f.height[l];
    let rh = f.height[r];
    if lh > 0.0 && rh > 0.0 {
        let ratio = lh.max(rh) / lh.min(rh);
        if ratio > cfg.max_adjacent_intensity_ratio {
            return None;
        }
    }

    let left_width = f.rt_end[l] - f.rt_start[l];
    let right_width = f.rt_end[r] - f.rt_start[r];
    let width = left_width.max(right_width).max(1e-9);
    let max_shift = cfg.max_apex_shift.max(cfg.max_apex_shift_width_fraction * width);
    let rt_shift = f.rt_apex[r] - f.rt_apex[l];
    if rt_shift.abs() > max_shift {
        return None;
    }

    let mz_score = 1.0 - (mz_error.abs() / tolerance).min(1.0);
    let shift_score = 1.0 - (rt_shift.abs() / max_shift).min(1.0);
    let overlap_score = rt_overlap(f, l, r);
    let i_score = intensity_score(f, l, r);
    let score = mz_score * 0.50 + shift_score * 0.20 + overlap_score * 0.15 + i_score * 0.15;
    if score < cfg.min_edge_score {
        return None;
    }
    Some(Edge { left: l, right: r, charge, score })
}

fn lower_bound(a: &[f64], x: f64) -> usize {
    let (mut lo, mut hi) = (0, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < x { lo = mid + 1; } else { hi = mid; }
    }
    lo
}
fn upper_bound_le(a: &[f64], x: f64) -> usize {
    let (mut lo, mut hi) = (0, a.len());
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] <= x { lo = mid + 1; } else { hi = mid; }
    }
    lo
}

pub fn build_edges(f: &Features, cfg: &Config) -> Vec<Edge> {
    let n = f.mz.len();
    (0..n)
        .into_par_iter()
        .flat_map_iter(|sorted_left| {
            let left_i = f.order[sorted_left];
            let left_mz = f.mz[left_i];
            let tol = cfg
                .isotope_mz_abs
                .max((left_mz + C13_DELTA) * cfg.isotope_mz_ppm / 1_000_000.0);
            // lo = searchsorted(sorted_mz, left_mz, "right"); hi = ... left_mz+C13+tol "right"
            let lo = upper_bound_le(&f.sorted_mz, left_mz);
            let hi = upper_bound_le(&f.sorted_mz, left_mz + C13_DELTA + tol);
            // best per derived charge
            let mut best_by_charge: std::collections::HashMap<i64, Edge> =
                std::collections::HashMap::new();
            for sr in lo..hi {
                let right_i = f.order[sr];
                if right_i == left_i {
                    continue;
                }
                let spacing = f.mz[right_i] - left_mz;
                if spacing <= 0.0 {
                    continue;
                }
                let charge = (C13_DELTA / spacing).round() as i64;
                if charge < 1 {
                    continue;
                }
                if let Some(edge) = score_edge(f, cfg, left_i, right_i, charge) {
                    match best_by_charge.get(&charge) {
                        Some(e) if e.score >= edge.score => {}
                        _ => {
                            best_by_charge.insert(charge, edge);
                        }
                    }
                }
            }
            best_by_charge.into_values().collect::<Vec<_>>()
        })
        .collect()
}

#[derive(Clone)]
pub struct MemberRow {
    pub feature_id: usize,
    pub isotope_index: i64,
    pub member_score: f64,
    pub mz_residual: f64,
    pub intensity_observed: f64,
    pub intensity_expected: f64,
    pub trace_score: f64,
}

#[derive(Clone)]
pub struct DistRow {
    pub distribution_id: i64,
    pub charge: i64,
    pub neutral_mass: f64,
    pub mono_mz: f64,
    pub rt_start: f64,
    pub rt_apex: f64,
    pub rt_end: f64,
    pub ms1_start: i64,
    pub ms1_apex: i64,
    pub ms1_end: i64,
    pub n_members: i64,
    pub score: f64,
    pub quality: f64,
    pub mz_score: f64,
    pub iso_score: f64,
    pub trace_score: f64,
    pub missing_score: f64,
    pub interloper_score: f64,
    pub mono_offset: i64,
    pub n_missing_interior: i64,
    pub n_interlopers: i64,
    pub ambiguity_score: f64,
    pub status: String,
    pub members: Vec<MemberRow>,
}

fn min_members_for_charge(charge: i64, cfg: &Config) -> usize {
    if charge == 1 {
        cfg.min_members_charge_one
    } else {
        cfg.min_distribution_members
    }
}

/// Plausible charges for an envelope seeded at `seed`, read off the spacing to
/// coeluting neighbours for isotope-index gaps q in 1..=MAX_Q.
fn derive_charges(f: &Features, cfg: &Config, seed: usize) -> Vec<i64> {
    let mut charges: Vec<i64> = Vec::new();
    let max_dmz = MAX_Q as f64 * C13_DELTA + f.mz_tol(cfg, f.mz[seed] + C13_DELTA);
    for j in f.neighbours_above(seed, max_dmz) {
        let dmz = f.mz[j] - f.mz[seed];
        if dmz <= 0.0 || !f.coelutes(seed, j) {
            continue;
        }
        for q in 1..=MAX_Q {
            let z = (q as f64 * C13_DELTA / dmz).round() as i64;
            if z < 1 {
                continue;
            }
            let expected = q as f64 * C13_DELTA / z as f64;
            if (dmz - expected).abs() <= f.mz_tol(cfg, f.mz[j]) && !charges.contains(&z) {
                charges.push(z);
            }
        }
    }
    charges
}

/// Occupy observed features on the z-lattice anchored at the seed (rel index 0).
fn grow_lattice(f: &Features, cfg: &Config, seed: usize, z: i64) -> std::collections::BTreeMap<i64, usize> {
    use std::collections::BTreeMap;
    let seed_mz = f.mz[seed];
    let spacing = C13_DELTA / z as f64;
    let min_sim = cfg.min_trace_similarity;
    let mut occupied: BTreeMap<i64, usize> = BTreeMap::new();
    occupied.insert(0, seed);
    let mut used: std::collections::HashSet<usize> = std::collections::HashSet::new();
    used.insert(seed);

    for &direction in &[1i64, -1i64] {
        let mut misses = 0;
        let mut k = direction;
        while misses <= MAX_ISOTOPE_GAP {
            let target = seed_mz + k as f64 * spacing;
            if target <= 0.0 {
                break;
            }
            let mut best_j: Option<usize> = None;
            let mut best_sim = min_sim;
            for j in f.features_near(cfg, target) {
                if used.contains(&j) || !f.coelutes(seed, j) {
                    continue;
                }
                let sim = f.trace_similarity(seed, j);
                if sim >= best_sim {
                    best_sim = sim;
                    best_j = Some(j);
                }
            }
            match best_j {
                Some(j) => {
                    occupied.insert(k, j);
                    used.insert(j);
                    misses = 0;
                }
                None => misses += 1,
            }
            k += direction;
        }
    }
    occupied
}

struct Scored {
    charge: i64,
    total: f64,
    mz_score: f64,
    iso_score: f64,
    trace_score: f64,
    missing_score: f64,
    interloper_score: f64,
    missing_interior: i64,
    interlopers: i64,
    neutral_mass: f64,
    mono_mz: f64,
    mono_offset: i64,
    members: Vec<MemberRow>,
    fids: Vec<usize>,
}

/// Fit monoisotope offset and score the whole envelope. None if below the floor.
fn score_envelope(f: &Features, cfg: &Config, occupied: &std::collections::BTreeMap<i64, usize>, z: i64) -> Option<Scored> {
    let rel: Vec<i64> = occupied.keys().copied().collect();
    let fids: Vec<usize> = rel.iter().map(|r| occupied[r]).collect();
    let n_obs = fids.len();
    if n_obs < min_members_for_charge(z, cfg) {
        return None;
    }
    let mz: Vec<f64> = fids.iter().map(|&i| f.mz[i]).collect();
    let heights: Vec<f64> = fids.iter().map(|&i| f.height[i]).collect();
    let hmax = heights.iter().copied().fold(0.0f64, f64::max);
    let obs_norm: Vec<f64> = if hmax > 0.0 {
        heights.iter().map(|h| h / hmax).collect()
    } else {
        heights.clone()
    };
    let spacing = C13_DELTA / z as f64;
    let base_rel = rel[0];
    let (member_trace_scores, trace_score) = f.consensus_trace_scores(&fids);

    let mut best: Option<Scored> = None;
    for offset in 0..=MAX_OFFSET_SEARCH {
        let iso_idx: Vec<i64> = rel.iter().map(|r| r - base_rel + offset).collect();
        let mono_mz = mz[0] - iso_idx[0] as f64 * spacing;
        let neutral_mass = (mono_mz - PROTON) * z as f64;
        if neutral_mass <= 0.0 || neutral_mass > cfg.max_neutral_mass {
            continue;
        }
        let span = (*iso_idx.iter().max().unwrap() + 1).max(n_obs as i64) as usize;
        let expected = averagine_envelope(neutral_mass, span);

        // m/z lattice residual
        let mut resid = vec![0.0f64; n_obs];
        let mut tol = vec![0.0f64; n_obs];
        let mut mz_acc = 0.0;
        for i in 0..n_obs {
            let lattice = mono_mz + iso_idx[i] as f64 * spacing;
            resid[i] = (mz[i] - lattice).abs();
            tol[i] = f.mz_tol(cfg, mz[i]);
            mz_acc += (1.0 - resid[i] / tol[i]).clamp(0.0, 1.0);
        }
        let mz_score = mz_acc / n_obs as f64;

        // averagine intensity agreement
        let exp_obs: Vec<f64> = iso_idx.iter().map(|&k| expected[k as usize]).collect();
        let denom: f64 = exp_obs.iter().map(|e| e * e).sum();
        let scale = if denom > 0.0 {
            obs_norm.iter().zip(&exp_obs).map(|(o, e)| o * e).sum::<f64>() / denom
        } else {
            0.0
        };
        let pred: Vec<f64> = exp_obs.iter().map(|e| scale * e).collect();
        let num: f64 = obs_norm.iter().zip(&pred).map(|(o, p)| o * p).sum();
        let na = obs_norm.iter().map(|x| x * x).sum::<f64>().sqrt();
        let nb = pred.iter().map(|x| x * x).sum::<f64>().sqrt();
        let iso_score = if na > 0.0 && nb > 0.0 { (num / (na * nb)).max(0.0) } else { 0.0 };

        // missing-expected-peak penalty (detection-aware)
        let present: std::collections::HashSet<i64> = iso_idx.iter().copied().collect();
        let obs_floor = obs_norm.iter().copied().fold(f64::INFINITY, f64::min);
        let imin = *iso_idx.iter().min().unwrap();
        let imax = *iso_idx.iter().max().unwrap();
        let mut missing_interior = 0i64;
        for k in imin..=imax {
            if !present.contains(&k) && scale * expected[k as usize] > obs_floor {
                missing_interior += 1;
            }
        }
        if !present.contains(&0) && scale * expected[0] > obs_floor.max(0.4) {
            missing_interior += 1;
        }
        let missing_score = 1.0 / (1.0 + missing_interior as f64);

        // interloper penalty (kills decimated aliases)
        let lattice_mzs: Vec<f64> = iso_idx.iter().map(|&k| mono_mz + k as f64 * spacing).collect();
        let fid_set: std::collections::HashSet<usize> = fids.iter().copied().collect();
        let mz_lo = mz.iter().copied().fold(f64::INFINITY, f64::min);
        let mz_hi = mz.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let margin = spacing * 0.5;
        let mut interlopers = 0i64;
        for j in f.features_in_range(mz_lo - 1e-6, mz_hi + 1e-6) {
            if fid_set.contains(&j) || !f.coelutes(fids[0], j) {
                continue;
            }
            if f.trace_similarity(fids[0], j) < cfg.min_trace_similarity {
                continue;
            }
            let jmz = f.mz[j];
            let mut min_dist = f64::INFINITY;
            let mut on_lattice = false;
            for &lm in &lattice_mzs {
                let d = (lm - jmz).abs();
                if d < min_dist {
                    min_dist = d;
                }
                if d <= f.mz_tol(cfg, jmz) {
                    on_lattice = true;
                }
            }
            if !on_lattice && min_dist < margin {
                interlopers += 1;
            }
        }
        let interloper_score = 1.0 / (1.0 + interlopers as f64);

        let total = 0.22 * mz_score
            + 0.34 * iso_score
            + 0.18 * trace_score
            + 0.13 * missing_score
            + 0.13 * interloper_score;

        let better = match &best {
            None => true,
            Some(b) => total > b.total,
        };
        if better {
            let members: Vec<MemberRow> = (0..n_obs)
                .map(|i| MemberRow {
                    feature_id: fids[i],
                    isotope_index: iso_idx[i],
                    member_score: (1.0 - resid[i] / tol[i]).clamp(0.0, 1.0),
                    mz_residual: resid[i],
                    intensity_observed: obs_norm[i],
                    intensity_expected: exp_obs[i],
                    trace_score: member_trace_scores[i],
                })
                .collect();
            best = Some(Scored {
                charge: z,
                total,
                mz_score,
                iso_score,
                trace_score,
                missing_score,
                interloper_score,
                missing_interior,
                interlopers,
                neutral_mass,
                mono_mz,
                mono_offset: offset,
                members,
                fids: fids.clone(),
            });
        }
    }
    best
}

pub fn build_distributions(f: &Features, cfg: &Config) -> Vec<DistRow> {
    let n = f.mz.len();
    let min_total = cfg.min_envelope_score;

    // Generate candidate envelopes across seeds/charges (parallel over seeds).
    let raw: Vec<Scored> = (0..n)
        .into_par_iter()
        .flat_map_iter(|seed| {
            let mut out: Vec<Scored> = Vec::new();
            for z in derive_charges(f, cfg, seed) {
                let occupied = grow_lattice(f, cfg, seed, z);
                if occupied.len() < min_members_for_charge(z, cfg) {
                    continue;
                }
                if let Some(s) = score_envelope(f, cfg, &occupied, z) {
                    if s.total >= min_total {
                        out.push(s);
                    }
                }
            }
            out
        })
        .collect();

    // Dedup by (charge, sorted feature set).
    let mut seen: std::collections::HashSet<(i64, Vec<usize>)> = std::collections::HashSet::new();
    let mut candidates: Vec<(Scored, i64)> = Vec::new();
    for s in raw {
        let zc = s.charge;
        let mut key_fids = s.fids.clone();
        key_fids.sort_unstable();
        if seen.insert((zc, key_fids)) {
            candidates.push((s, zc));
        }
    }

    compete(f, cfg, candidates)
}

fn compete(f: &Features, cfg: &Config, mut candidates: Vec<(Scored, i64)>) -> Vec<DistRow> {
    candidates.sort_by(|a, b| b.0.total.partial_cmp(&a.0.total).unwrap());
    let margin = cfg.ambiguity_margin;
    let mut claimed: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut winner_of: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut runner_up: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    let mut kept: Vec<(Scored, i64)> = Vec::new();

    for (s, z) in candidates.into_iter() {
        let fids = s.fids.clone();
        let conflict: Vec<usize> = fids.iter().copied().filter(|x| claimed.contains(x)).collect();
        if !conflict.is_empty() {
            for cfid in conflict {
                if let Some(&wi) = winner_of.get(&cfid) {
                    if kept[wi].1 != z {
                        let e = runner_up.entry(wi).or_insert(0.0);
                        if s.total > *e {
                            *e = s.total;
                        }
                    }
                }
            }
            continue;
        }
        let idx = kept.len();
        for &fid in &fids {
            claimed.insert(fid);
            winner_of.insert(fid, idx);
        }
        kept.push((s, z));
    }

    let mut rows: Vec<DistRow> = Vec::new();
    for (idx, (s, z)) in kept.iter().enumerate() {
        let fids: Vec<usize> = s.fids.clone();
        let runner = *runner_up.get(&idx).unwrap_or(&0.0);
        let ambiguity_score = if s.total > 0.0 { runner / s.total } else { 0.0 };
        let status = if runner > 0.0 && (s.total - runner) < margin * s.total {
            "ambiguous".to_string()
        } else {
            "validated".to_string()
        };
        let apex_local = {
            let mut bi = 0;
            for i in 1..fids.len() {
                if f.height[fids[i]] > f.height[fids[bi]] {
                    bi = i;
                }
            }
            fids[bi]
        };
        rows.push(DistRow {
            distribution_id: -1,
            charge: *z,
            neutral_mass: s.neutral_mass,
            mono_mz: s.mono_mz,
            rt_start: fids.iter().map(|&p| f.rt_start[p]).fold(f64::INFINITY, f64::min),
            rt_apex: f.rt_apex[apex_local],
            rt_end: fids.iter().map(|&p| f.rt_end[p]).fold(f64::NEG_INFINITY, f64::max),
            ms1_start: fids.iter().map(|&p| f.ms1_start[p]).min().unwrap(),
            ms1_apex: f.ms1_apex[apex_local],
            ms1_end: fids.iter().map(|&p| f.ms1_end[p]).max().unwrap(),
            n_members: fids.len() as i64,
            score: s.total,
            quality: s.total,
            mz_score: s.mz_score,
            iso_score: s.iso_score,
            trace_score: s.trace_score,
            missing_score: s.missing_score,
            interloper_score: s.interloper_score,
            mono_offset: s.mono_offset,
            n_missing_interior: s.missing_interior,
            n_interlopers: s.interlopers,
            ambiguity_score,
            status,
            members: s.members.clone(),
        });
    }

    rows.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.partial_cmp(&b.neutral_mass).unwrap())
            .then(a.rt_apex.partial_cmp(&b.rt_apex).unwrap())
            .then(a.mono_mz.partial_cmp(&b.mono_mz).unwrap())
    });
    for (i, d) in rows.iter_mut().enumerate() {
        d.distribution_id = i as i64;
    }
    rows
}


#[derive(Clone)]
pub struct AnalyteRow {
    pub analyte_id: i64,
    pub neutral_mass: f64,
    pub rt_start: f64,
    pub rt_apex: f64,
    pub rt_end: f64,
    pub ms1_start: i64,
    pub ms1_apex: i64,
    pub ms1_end: i64,
    pub charge_min: i64,
    pub charge_max: i64,
    pub n_distributions: i64,
    pub score: f64,
    pub members: Vec<(i64, i64)>, // (distribution_id, charge)
}

fn distribution_rt_score(a: &DistRow, b: &DistRow) -> f64 {
    let overlap = a.rt_end.min(b.rt_end) - a.rt_start.max(b.rt_start);
    let union = a.rt_end.max(b.rt_end) - a.rt_start.min(b.rt_start);
    let raw = if union <= 0.0 { 0.0 } else { (overlap / union).max(0.0) };
    let apex_gap = (a.rt_apex - b.rt_apex).abs();
    let width = (a.rt_end - a.rt_start).max(b.rt_end - b.rt_start).max(1e-9);
    let apex_score = (1.0 - apex_gap / width).max(0.0);
    (raw * 0.4 + apex_score * 0.6).clamp(0.0, 1.0)
}

pub fn build_analytes(cfg: &Config, dists: &[DistRow]) -> Vec<AnalyteRow> {
    if dists.is_empty() {
        return vec![];
    }
    // Seed-based grouping (replaces transitive UnionFind, which chained distinct
    // RT peaks into giant analytes). Best-quality distribution first; each
    // unclaimed seed pulls in the single best-coeluting distribution at each
    // OTHER charge within neutral-mass tolerance -> at most one per charge, tight
    // in RT.
    let mut order: Vec<usize> = (0..dists.len()).collect();
    order.sort_by(|&a, &b| dists[a].neutral_mass.partial_cmp(&dists[b].neutral_mass).unwrap());
    let sorted: Vec<&DistRow> = order.iter().map(|&i| &dists[i]).collect();
    let masses: Vec<f64> = sorted.iter().map(|d| d.neutral_mass).collect();

    let mut by_quality: Vec<usize> = (0..sorted.len()).collect();
    by_quality.sort_by(|&a, &b| sorted[b].quality.partial_cmp(&sorted[a].quality).unwrap());

    let mut claimed = vec![false; sorted.len()];
    let mut groups: Vec<Vec<usize>> = Vec::new();

    for &idx in &by_quality {
        if claimed[idx] {
            continue;
        }
        claimed[idx] = true;
        let seed = sorted[idx];
        let mut members = vec![idx];
        let mut charges_used: std::collections::HashSet<i64> = std::collections::HashSet::new();
        charges_used.insert(seed.charge);

        let tol = 0.002_f64.max(seed.neutral_mass * cfg.charge_mass_ppm / 1_000_000.0);
        let start = lower_bound(&masses, seed.neutral_mass - tol);
        let end = upper_bound_le(&masses, seed.neutral_mass + tol);

        let mut candidates: Vec<(f64, usize)> = Vec::new();
        for j in start..end {
            if claimed[j] || j == idx || charges_used.contains(&sorted[j].charge) {
                continue;
            }
            let rt_score = distribution_rt_score(seed, sorted[j]);
            if rt_score >= cfg.min_charge_group_rt_score {
                candidates.push((rt_score, j));
            }
        }
        candidates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
        for (_, j) in candidates {
            if claimed[j] || charges_used.contains(&sorted[j].charge) {
                continue;
            }
            claimed[j] = true;
            charges_used.insert(sorted[j].charge);
            members.push(j);
        }
        groups.push(members);
    }

    let mut analytes = vec![];
    for members in groups.iter() {
        let mut wsum = 0.0;
        let mut wmass = 0.0;
        let mut score_sum = 0.0;
        let mut charge_min = i64::MAX;
        let mut charge_max = i64::MIN;
        let mut rt_start = f64::INFINITY;
        let mut rt_end = f64::NEG_INFINITY;
        let mut ms1_start = i64::MAX;
        let mut ms1_end = i64::MIN;
        let mut apex_idx = members[0];
        let mut best_quality = f64::NEG_INFINITY;
        for &mi in members {
            let d = sorted[mi];
            let w = d.score.max(1e-6);
            wsum += w;
            wmass += d.neutral_mass * w;
            score_sum += d.score;
            charge_min = charge_min.min(d.charge);
            charge_max = charge_max.max(d.charge);
            rt_start = rt_start.min(d.rt_start);
            rt_end = rt_end.max(d.rt_end);
            ms1_start = ms1_start.min(d.ms1_start);
            ms1_end = ms1_end.max(d.ms1_end);
            if d.quality > best_quality {
                best_quality = d.quality;
                apex_idx = mi;
            }
        }
        let apex = sorted[apex_idx];
        analytes.push(AnalyteRow {
            analyte_id: analytes.len() as i64,
            neutral_mass: wmass / wsum,
            rt_start,
            rt_apex: apex.rt_apex,
            rt_end,
            ms1_start,
            ms1_apex: apex.ms1_apex,
            ms1_end,
            charge_min,
            charge_max,
            n_distributions: members.len() as i64,
            score: score_sum / members.len() as f64,
            members: members.iter().map(|&mi| (sorted[mi].distribution_id, sorted[mi].charge)).collect(),
        });
    }
    analytes
}
