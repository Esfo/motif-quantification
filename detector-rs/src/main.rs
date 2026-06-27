use anyhow::Result;
use clap::Parser;
use std::time::Instant;

use mzdata::io::MzMLReader;
use mzdata::prelude::*;

mod config;
mod peaks;
mod linemodel;
mod distributions;
mod store;

use config::Config;
use distributions::{build_analytes, build_distributions, build_edges, Features};
use linemodel::LineModel;
use store::{compress_f32, write_db, ScanPoints, ScanRow};

#[derive(Parser, Debug)]
#[command(about = "MS1 distribution detector (Rust port of distributions/index_ms1.py)")]
struct Args {
    mzml: String,
    #[arg(long)]
    out: String,
    #[arg(long, default_value_t = false)]
    overwrite: bool,
    #[arg(long, default_value_t = 0.0)]
    min_intensity: f64,
    #[arg(long, default_value_t = 0)]
    threads: usize,
    #[arg(long, default_value_t = 500)]
    progress: i64,
    /// Store per-scan centroids in the sqlite (scan_points) so the GUI reads raw
    /// points from the db instead of re-decoding the mzML. --store-points=false
    /// keeps the sqlite small.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    store_points: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .ok();
    }

    let started = Instant::now();
    let cfg = Config::default();
    let mut model = LineModel::new(cfg.clone());
    let mut scans: Vec<ScanRow> = Vec::new();
    let mut points: Vec<ScanPoints> = Vec::new();

    let t_line = Instant::now();
    let reader = MzMLReader::open_path(&args.mzml)?;
    let mut ms1_index: i64 = 0;
    for spec in reader {
        if spec.ms_level() != 1 {
            continue;
        }
        let rt = spec.start_time();
        let arrays = match spec.arrays.as_ref() {
            Some(a) => a,
            None => continue,
        };
        let mz_arr = arrays.mzs()?;
        let int_arr = arrays.intensities()?;
        let mut mzs: Vec<f64> = Vec::with_capacity(mz_arr.len());
        let mut ints: Vec<f64> = Vec::with_capacity(int_arr.len());
        for (m, i) in mz_arr.iter().zip(int_arr.iter()) {
            let iv = *i as f64;
            if args.min_intensity > 0.0 && iv < args.min_intensity {
                continue;
            }
            mzs.push(*m);
            ints.push(iv);
        }
        let tic: f64 = ints.iter().sum();
        scans.push(ScanRow {
            ms1_index,
            spectrum_index: spec.index() as i64,
            scan_id: spec.id().to_string(),
            rt,
            tic,
            n_points: mzs.len() as i64,
        });
        if args.store_points {
            let mz_f32: Vec<f32> = mzs.iter().map(|v| *v as f32).collect();
            let int_f32: Vec<f32> = ints.iter().map(|v| *v as f32).collect();
            points.push(ScanPoints {
                ms1_index,
                n: mzs.len() as i64,
                mz_blob: compress_f32(&mz_f32),
                intensity_blob: compress_f32(&int_f32),
            });
        }
        model.process_scan(ms1_index, rt, &mzs, &ints);
        if args.progress > 0 && ms1_index > 0 && ms1_index % args.progress == 0 {
            eprintln!("scans={ms1_index} active+closed (streaming)");
        }
        ms1_index += 1;
    }
    model.finalize();
    let line_secs = t_line.elapsed().as_secs_f64();

    let t_edge = Instant::now();
    let f = Features::from_rows(&model.features);
    let edges = build_edges(&f, &cfg);
    let edge_secs = t_edge.elapsed().as_secs_f64();

    let t_dist = Instant::now();
    let dists = build_distributions(&f, &cfg, &edges);
    let dist_secs = t_dist.elapsed().as_secs_f64();

    let t_charge = Instant::now();
    let analytes = build_analytes(&cfg, &dists);
    let charge_secs = t_charge.elapsed().as_secs_f64();

    let t_write = Instant::now();
    let params = vec![
        ("script".to_string(), "\"detector-rs/ms1-detector\"".to_string()),
        ("script_version".to_string(), "\"0.1.0\"".to_string()),
        ("mzml".to_string(), serde_json_string(&args.mzml)),
        (
            "counts".to_string(),
            format!(
                "{{\"scans\":{},\"lines\":{},\"features\":{},\"edges\":{},\"distributions\":{},\"analytes\":{}}}",
                scans.len(), model.lines.len(), model.features.len(), edges.len(), dists.len(), analytes.len()
            ),
        ),
    ];
    write_db(&args.out, args.overwrite, &params, &scans, &model.lines, &model.features, &dists, &analytes, &points)?;
    let write_secs = t_write.elapsed().as_secs_f64();

    let total = started.elapsed().as_secs_f64();
    println!(
        "{{\n  \"out\": \"{}\",\n  \"seconds\": {:.3},\n  \"stage_seconds\": {{\"line\": {:.3}, \"edges\": {:.3}, \"distributions\": {:.3}, \"charge\": {:.3}, \"write\": {:.3}}},\n  \"scans\": {},\n  \"lines\": {},\n  \"features\": {},\n  \"best_edges\": {},\n  \"distributions\": {},\n  \"analytes\": {}\n}}",
        args.out, total, line_secs, edge_secs, dist_secs, charge_secs, write_secs,
        scans.len(), model.lines.len(), model.features.len(), edges.len(), dists.len(), analytes.len()
    );

    // charge histogram + query
    let mut hist: std::collections::BTreeMap<i64, i64> = std::collections::BTreeMap::new();
    for d in &dists {
        *hist.entry(d.charge).or_insert(0) += 1;
    }
    let total_d = dists.len().max(1) as f64;
    eprintln!("distributions by charge:");
    for (c, n) in &hist {
        eprintln!("  z={:<3} {:>8}  ({:.1}%)", c, n, 100.0 * *n as f64 / total_d);
    }
    eprintln!(
        "sqlite3 {} \"SELECT charge, COUNT(*) AS n FROM distributions GROUP BY charge ORDER BY charge;\"",
        args.out
    );
    Ok(())
}

fn serde_json_string(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            _ => out.push(c),
        }
    }
    out.push('"');
    out
}
