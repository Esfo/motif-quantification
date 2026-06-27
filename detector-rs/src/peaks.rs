//! Faithful port of distributions/peaks.py (axis_peaks + helpers).

pub fn moving_average(array: &[f64], width: usize) -> Vec<f64> {
    let n = array.len();
    if width <= 1 || n < width {
        return array.to_vec();
    }
    // np.convolve(array, ones(width)/width, mode="same")
    let full_len = n + width - 1;
    let mut full = vec![0.0f64; full_len];
    let inv = 1.0 / width as f64;
    for (k, &a) in array.iter().enumerate() {
        for j in 0..width {
            full[k + j] += a * inv;
        }
    }
    let start = (width - 1) / 2;
    full[start..start + n].to_vec()
}

fn local_extrema(narray: &[f64]) -> (Vec<usize>, Vec<usize>) {
    // Returns (mins, maxes) interior-or-edge per peaks.py logic.
    let n = narray.len();
    if n == 0 {
        return (vec![], vec![]);
    }
    let mut fmax = vec![false; n];
    let mut bmax = vec![false; n];
    let mut fmin = vec![false; n];
    let mut bmin = vec![false; n];
    for i in 0..n - 1 {
        fmax[i] = narray[i] > narray[i + 1];
        fmin[i] = narray[i] < narray[i + 1];
    }
    // backward checks shifted by one with the first element duplicated
    bmax[0] = fmax[0];
    bmin[0] = fmin[0];
    for i in 1..n {
        bmax[i] = narray[i] > narray[i - 1];
        bmin[i] = narray[i] < narray[i - 1];
    }
    fmax[n - 1] = bmax[n - 1];
    fmin[n - 1] = bmin[n - 1];
    let mut mins = vec![];
    let mut maxes = vec![];
    for i in 0..n {
        if fmin[i] && bmin[i] {
            mins.push(i);
        }
        if fmax[i] && bmax[i] {
            maxes.push(i);
        }
    }
    (mins, maxes)
}

pub fn minpoint_reduction(barray: &[f64], mindist: usize) -> Vec<usize> {
    let total = barray.len();
    let mut mask = vec![false; total];
    let mut extramaxes: Vec<usize> = vec![];
    let mut last_maxes: Vec<usize> = vec![];
    let mut last_narray_argmax = 0usize;

    loop {
        // narray = barray[~mask]; keep mapping kept-index -> original-index
        let mut orig_index = vec![];
        let mut narray = vec![];
        for i in 0..total {
            if !mask[i] {
                orig_index.push(i);
                narray.push(barray[i]);
            }
        }
        if narray.is_empty() {
            return vec![];
        }
        // track argmax of narray for the fallback
        let mut am = 0usize;
        for i in 1..narray.len() {
            if narray[i] > narray[am] {
                am = i;
            }
        }
        last_narray_argmax = orig_index[am];

        let (mut mins, maxes) = local_extrema(&narray);
        last_maxes = maxes.iter().map(|&m| orig_index[m]).collect();

        let mut extremas: Vec<usize> = mins.iter().chain(maxes.iter()).copied().collect();
        extremas.sort_unstable();
        extremas.dedup();
        if extremas.is_empty() {
            break;
        }

        // adjacency: extremadistances[i][j] = |e_i - e_j| < mindist, diagonal false
        let m = extremas.len();
        let mut any_adjacent = false;
        let mut isolated = vec![true; m]; // ~extremadistances.any(axis=0)
        for i in 0..m {
            for j in 0..m {
                if i == j {
                    continue;
                }
                let d = (extremas[i] as isize - extremas[j] as isize).unsigned_abs();
                if d < mindist {
                    isolated[i] = false;
                    any_adjacent = true;
                }
            }
        }

        let maxes_set: std::collections::HashSet<usize> = maxes.iter().copied().collect();
        let mins_set: std::collections::HashSet<usize> = mins.iter().copied().collect();

        // separatedextremas = extremas where isolated
        let mut newmask_kept: std::collections::HashSet<usize> =
            mins.iter().copied().collect(); // indices into narray that are mins
        let mut maintained_mins: Vec<usize> = vec![];
        for (idx, &e) in extremas.iter().enumerate() {
            if !isolated[idx] {
                continue;
            }
            if maxes_set.contains(&e) {
                // maxestomaintain: map narray index e -> original index
                extramaxes.push(orig_index[e]);
            } else if mins_set.contains(&e) {
                // keep this min out of the mask
                newmask_kept.remove(&e);
                maintained_mins.push(e);
            }
        }
        // mins after removing maintained ones
        if !maintained_mins.is_empty() {
            let maintained: std::collections::HashSet<usize> =
                maintained_mins.iter().copied().collect();
            mins.retain(|x| !maintained.contains(x));
        }

        if any_adjacent && !mins.is_empty() {
            // mask the kept-min original indices that are still in newmask_kept
            for &k in &newmask_kept {
                mask[orig_index[k]] = true;
            }
        } else {
            break;
        }
    }

    let mut maxes_orig = last_maxes;
    if maxes_orig.is_empty() {
        maxes_orig = vec![last_narray_argmax];
    }
    let mut fmaxes: Vec<usize> = maxes_orig.into_iter().chain(extramaxes.into_iter()).collect();
    fmaxes.sort_unstable();
    fmaxes.dedup();
    fmaxes
}

pub fn boundary_finding(fmaxes: &[usize], array: &[f64]) -> Vec<(usize, usize, usize)> {
    if fmaxes.is_empty() {
        return vec![];
    }
    let n = array.len();
    let mut anchors = vec![0usize];
    anchors.extend_from_slice(fmaxes);
    anchors.push(n - 1);

    let mut peakbounds: Vec<(Option<usize>, Option<usize>)> = vec![];

    for idx in 0..anchors.len() - 1 {
        let left_anchor = anchors[idx];
        let right_anchor = anchors[idx + 1] + 1;

        if idx > 0 {
            // right bound of previous peak
            let series: Vec<f64> = array[left_anchor..right_anchor.min(n)].to_vec();
            let mut acc = f64::INFINITY;
            let mut trim = vec![];
            for &v in &series {
                acc = acc.min(v);
                trim.push(v <= acc);
            }
            // np.trim_zeros(trim, 'b').size : length after trimming trailing false
            let mut end = trim.len();
            while end > 0 && !trim[end - 1] {
                end -= 1;
            }
            let nr = left_anchor + end;
            let rseries: Vec<f64> = array[left_anchor..nr.min(n)].to_vec();
            let rightbound = if !rseries.is_empty() {
                let mut minv = rseries[0];
                let mut mini = 0usize;
                for (k, &v) in rseries.iter().enumerate() {
                    if v < minv {
                        minv = v;
                        mini = k;
                    }
                }
                left_anchor + mini + 1
            } else {
                left_anchor + 1
            };
            if let Some(last) = peakbounds.last_mut() {
                last.1 = Some(rightbound);
            }
        }

        if idx < anchors.len() - 1 - 1 {
            // left bound of next peak
            let series: Vec<f64> = array[left_anchor..right_anchor.min(n)].to_vec();
            // flip, min-accumulate, flip
            let mut acc = f64::INFINITY;
            let mut trim = vec![false; series.len()];
            for k in (0..series.len()).rev() {
                acc = acc.min(series[k]);
                trim[k] = series[k] <= acc;
            }
            // np.trim_zeros(trim, 'f').size : length after trimming leading false
            let mut startf = 0usize;
            while startf < trim.len() && !trim[startf] {
                startf += 1;
            }
            let leftestimate = trim.len() - startf;
            let nl = right_anchor.saturating_sub(leftestimate);
            let lseries: Vec<f64> = array[nl.min(n)..right_anchor.min(n)].to_vec();
            let leftbound = if !lseries.is_empty() {
                let mut minv = lseries[0];
                let mut mini = 0usize;
                for (k, &v) in lseries.iter().enumerate() {
                    if v <= minv {
                        minv = v;
                        mini = k;
                    }
                }
                nl + mini
            } else {
                left_anchor
            };
            peakbounds.push((Some(leftbound), None));
        }
    }

    let mut out = vec![];
    for (i, &apex) in fmaxes.iter().enumerate() {
        if i < peakbounds.len() {
            let (l, r) = peakbounds[i];
            let left = l.unwrap_or(0);
            let right = r.unwrap_or(n - 1);
            out.push((left, apex, right));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

pub fn axis_peaks(array: &[f64], mindist: usize) -> Vec<(usize, usize, usize)> {
    if array.is_empty() {
        return vec![];
    }
    let maxes = minpoint_reduction(array, mindist);
    boundary_finding(&maxes, array)
}
