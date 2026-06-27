//! Port of the edge + distribution stages from index_ml1.py.

use rayon::prelude::*;

use crate::config::{Config, C13_DELTA, PROTON};
use crate::linemodel::FeatureRow;

pub struct Features {
    pub mz: Vec<f64>,
    pub rt_start: Vec<f64>,
    pub rt_apex: Vec<f64>,
    pub rt_end: Vec<f64>,
    pub ms1_start: Vec<i64>,
    pub ms1_apex: Vec<i64>,
    pub ms1_end: Vec<i64>,
    pub height: Vec<f64>,
    // m/z-sorted view
    pub order: Vec<usize>,
    pub sorted_mz: Vec<f64>,
}

impl Features {
    pub fn from_rows(rows: &[FeatureRow]) -> Self {
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
            mz,
            order,
            sorted_mz,
        }
    }
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

fn envelope_shape_score(heights: &[f64]) -> f64 {
    if heights.len() < 3 {
        return 1.0;
    }
    let mut valleys = 0;
    for i in 1..heights.len() - 1 {
        let lower = heights[i - 1].min(heights[i + 1]);
        if heights[i] < heights[i - 1] && heights[i] < heights[i + 1] && heights[i] < 0.7 * lower {
            valleys += 1;
        }
    }
    1.0 / (1.0 + valleys as f64)
}

fn coherent_rt_path(f: &Features, cfg: &Config, path: &[usize]) -> bool {
    if path.len() < 3 {
        return true;
    }
    let m = path.len();
    let ys: Vec<f64> = path.iter().map(|&p| f.rt_apex[p]).collect();
    let mut widths: Vec<f64> = path
        .iter()
        .map(|&p| (f.rt_end[p] - f.rt_start[p]).max(1e-9))
        .collect();
    // linear least squares y = slope*x + intercept, x = 0..m-1
    let xbar = (m as f64 - 1.0) / 2.0;
    let ybar = ys.iter().sum::<f64>() / m as f64;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for i in 0..m {
        let dx = i as f64 - xbar;
        sxy += dx * (ys[i] - ybar);
        sxx += dx * dx;
    }
    let slope = if sxx != 0.0 { sxy / sxx } else { 0.0 };
    let intercept = ybar - slope * xbar;
    let mut residual = 0.0f64;
    for i in 0..m {
        let pred = slope * i as f64 + intercept;
        residual = residual.max((ys[i] - pred).abs());
    }
    widths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_w = if m % 2 == 1 {
        widths[m / 2]
    } else {
        (widths[m / 2 - 1] + widths[m / 2]) / 2.0
    };
    let allowed = cfg.max_apex_shift.max(cfg.max_apex_shift_width_fraction * median_w);
    residual <= allowed
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
    pub members: Vec<(usize, i64, f64)>, // (feature_id, isotope_index, member_score)
}

fn dist_worker_for_charge(f: &Features, cfg: &Config, charge: i64, edges: &[&Edge]) -> Vec<DistRow> {
    use std::collections::{HashMap, HashSet};
    // best_out: left -> edge (first seen in score-desc order)
    let mut idx: Vec<usize> = (0..edges.len()).collect();
    idx.sort_by(|&a, &b| edges[b].score.partial_cmp(&edges[a].score).unwrap());
    let mut best_out: HashMap<usize, &Edge> = HashMap::new();
    let mut incoming: HashSet<usize> = HashSet::new();
    for &i in &idx {
        let e = edges[i];
        if !best_out.contains_key(&e.left) {
            best_out.insert(e.left, e);
            incoming.insert(e.right);
        }
    }
    let mut starts: Vec<usize> = best_out
        .values()
        .map(|e| e.left)
        .filter(|l| !incoming.contains(l))
        .collect();
    starts.sort_unstable();

    let min_members = min_members_for_charge(charge, cfg);

    let mut rows = vec![];
    let mut used_paths: HashSet<Vec<usize>> = HashSet::new();

    for start in starts {
        let mut path = vec![start];
        let mut path_edges: Vec<&Edge> = vec![];
        let mut current = start;
        loop {
            let edge = match best_out.get(&current) {
                Some(e) => *e,
                None => break,
            };
            let next = edge.right;
            if path.contains(&next) {
                break;
            }
            let mut tentative = path.clone();
            tentative.push(next);
            if !coherent_rt_path(f, cfg, &tentative) {
                break;
            }
            path_edges.push(edge);
            path = tentative;
            current = next;
        }
        if path.len() < min_members {
            continue;
        }
        if used_paths.contains(&path) {
            continue;
        }
        used_paths.insert(path.clone());

        let heights: Vec<f64> = path.iter().map(|&p| f.height[p]).collect();
        let apex_local = {
            let mut bi = 0;
            for i in 1..heights.len() {
                if heights[i] > heights[bi] {
                    bi = i;
                }
            }
            bi
        };
        let apex_feature = path[apex_local];
        let mono = path[0];
        let score = if path_edges.is_empty() {
            0.0
        } else {
            path_edges.iter().map(|e| e.score).sum::<f64>() / path_edges.len() as f64
        };
        let shape = envelope_shape_score(&heights);
        let quality = score * (path.len() as f64).sqrt() * shape;

        let members: Vec<(usize, i64, f64)> = path
            .iter()
            .enumerate()
            .map(|(iso, &p)| {
                let ms = if iso == 0 { 1.0 } else { path_edges[iso - 1].score };
                (p, iso as i64, ms)
            })
            .collect();

        rows.push(DistRow {
            distribution_id: -1,
            charge,
            neutral_mass: (f.mz[mono] - PROTON) * charge as f64,
            mono_mz: f.mz[mono],
            rt_start: path.iter().map(|&p| f.rt_start[p]).fold(f64::INFINITY, f64::min),
            rt_apex: f.rt_apex[apex_feature],
            rt_end: path.iter().map(|&p| f.rt_end[p]).fold(f64::NEG_INFINITY, f64::max),
            ms1_start: path.iter().map(|&p| f.ms1_start[p]).min().unwrap(),
            ms1_apex: f.ms1_apex[apex_feature],
            ms1_end: path.iter().map(|&p| f.ms1_end[p]).max().unwrap(),
            n_members: path.len() as i64,
            score,
            quality,
            members,
        });
    }
    rows
}

fn min_members_for_charge(charge: i64, cfg: &Config) -> usize {
    if charge == 1 {
        cfg.min_members_charge_one
    } else if charge <= cfg.high_charge_threshold {
        cfg.min_distribution_members
    } else {
        let extra = ((charge - cfg.high_charge_threshold - 1) / 5).max(0) as usize;
        cfg.high_charge_min_members + extra
    }
}

fn charge_prior(charge: i64, cfg: &Config) -> f64 {
    if charge == 1 {
        cfg.charge_one_score_penalty
    } else {
        1.0 / (1.0 + cfg.charge_prior_strength * (charge - 2) as f64)
    }
}

fn resolve_competition(cfg: &Config, rows: Vec<DistRow>) -> Vec<DistRow> {
    let rank = |r: &DistRow| -> f64 { r.quality * charge_prior(r.charge, cfg) };
    let mut idx: Vec<usize> = (0..rows.len()).collect();
    idx.sort_by(|&a, &b| rank(&rows[b]).partial_cmp(&rank(&rows[a])).unwrap());
    let mut claimed: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut kept = vec![];
    for i in idx {
        let r = &rows[i];
        if r.members.iter().any(|(fid, _, _)| claimed.contains(fid)) {
            continue;
        }
        for (fid, _, _) in &r.members {
            claimed.insert(*fid);
        }
        kept.push(rows[i].clone());
    }
    kept
}

pub fn build_distributions(f: &Features, cfg: &Config, edges: &[Edge]) -> Vec<DistRow> {
    use std::collections::BTreeMap;
    let mut by_charge: BTreeMap<i64, Vec<&Edge>> = BTreeMap::new();
    for e in edges {
        by_charge.entry(e.charge).or_default().push(e);
    }
    let mut all_rows = vec![];
    for (charge, ce) in &by_charge {
        all_rows.extend(dist_worker_for_charge(f, cfg, *charge, ce));
    }
    let mut kept = resolve_competition(cfg, all_rows);
    // assign distribution_id order: sort by (charge, neutral_mass, rt_apex, mono_mz)
    kept.sort_by(|a, b| {
        a.charge
            .cmp(&b.charge)
            .then(a.neutral_mass.partial_cmp(&b.neutral_mass).unwrap())
            .then(a.rt_apex.partial_cmp(&b.rt_apex).unwrap())
            .then(a.mono_mz.partial_cmp(&b.mono_mz).unwrap())
    });
    for (i, d) in kept.iter_mut().enumerate() {
        d.distribution_id = i as i64;
    }
    kept
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
