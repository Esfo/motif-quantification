use std::io::{self, BufRead, Write};
use std::collections::BTreeMap;
use ms1_detector::config::Config;
use ms1_detector::linemodel::LineModel;
use ms1_detector::distributions::{Features, build_edges, build_distributions, build_analytes};

fn main() {
    let cfg = Config::default();
    let mut model = LineModel::new(cfg.clone());
    for (idx, line) in io::stdin().lock().lines().enumerate() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let rt: f64 = match it.next() { Some(t) => t.parse().unwrap(), None => continue };
        let mut mzs = vec![]; let mut ints = vec![];
        for tok in it {
            let mut p = tok.split(':');
            mzs.push(p.next().unwrap().parse().unwrap());
            ints.push(p.next().unwrap().parse().unwrap());
        }
        model.process_scan(idx as i64, rt, &mzs, &ints);
    }
    model.finalize();
    let f = Features::from_rows(&model.features);
    let edges = build_edges(&f, &cfg);
    let dists = build_distributions(&f, &cfg, &edges);
    let analytes = build_analytes(&cfg, &dists);

    let stdout = io::stdout(); let mut out = stdout.lock();
    let mut hist: BTreeMap<i64, i64> = BTreeMap::new();
    for d in &dists { *hist.entry(d.charge).or_insert(0) += 1; }
    writeln!(out, "FEATURES {} EDGES {} DISTS {} ANALYTES {}", model.features.len(), edges.len(), dists.len(), analytes.len()).unwrap();
    for (c, n) in &hist { writeln!(out, "z={c} {n}").unwrap(); }
    let mut rows: Vec<String> = dists.iter().map(|d| {
        format!("{} {:.5} {:.5} {:.4} {} {:.5} {:.5}", d.charge, d.neutral_mass, d.mono_mz, d.rt_apex, d.n_members, d.score, d.quality)
    }).collect();
    rows.sort();
    for r in rows { writeln!(out, "{r}").unwrap(); }
    let mut arows: Vec<String> = analytes.iter().map(|a| {
        format!("A {:.5} {:.4} {} {} {} {:.5}", a.neutral_mass, a.rt_apex, a.charge_min, a.charge_max, a.n_distributions, a.score)
    }).collect();
    arows.sort();
    for r in arows { writeln!(out, "{r}").unwrap(); }
}
