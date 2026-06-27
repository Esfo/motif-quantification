//! SQLite writer replicating distributions/store.py schema exactly.

use anyhow::Result;
use rusqlite::{params, Connection};

use crate::distributions::{AnalyteRow, DistRow};
use crate::linemodel::{FeatureRow, LineRow};

pub struct ScanRow {
    pub ms1_index: i64,
    pub spectrum_index: i64,
    pub scan_id: String,
    pub rt: f64,
    pub tic: f64,
    pub n_points: i64,
}

/// One scan's centroids, zlib-compressed little-endian f32 arrays, so the GUI
/// can pull a window's raw points straight from sqlite instead of re-decoding
/// the mzML (Option B).
pub struct ScanPoints {
    pub ms1_index: i64,
    pub n: i64,
    pub mz_blob: Vec<u8>,
    pub intensity_blob: Vec<u8>,
}

/// One feature's chromatographic trace, zlib-compressed little-endian arrays
/// (scans i32, rts/mzs/intensities f32) — mirrors distributions/store.py.
pub struct FeatureTraceRow {
    pub feature_id: i64,
    pub n: i64,
    pub scans_blob: Vec<u8>,
    pub rts_blob: Vec<u8>,
    pub mzs_blob: Vec<u8>,
    pub intensities_blob: Vec<u8>,
}

fn zlib_compress(bytes: &[u8]) -> Vec<u8> {
    use flate2::write::ZlibEncoder;
    use flate2::Compression;
    use std::io::Write;
    let mut enc = ZlibEncoder::new(Vec::new(), Compression::fast());
    enc.write_all(bytes).expect("zlib write");
    enc.finish().expect("zlib finish")
}

pub fn compress_f32(values: &[f32]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    zlib_compress(&bytes)
}

pub fn compress_i32(values: &[i32]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(values.len() * 4);
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    zlib_compress(&bytes)
}

const SCHEMA: &str = r#"
CREATE TABLE IF NOT EXISTS parameters (key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS scans (ms1_index INTEGER PRIMARY KEY, spectrum_index INTEGER NOT NULL, scan_id TEXT NOT NULL, rt REAL NOT NULL, tic REAL NOT NULL, n_points INTEGER NOT NULL);
CREATE TABLE IF NOT EXISTS scan_points (ms1_index INTEGER PRIMARY KEY, n INTEGER NOT NULL, mz BLOB NOT NULL, intensity BLOB NOT NULL);
CREATE TABLE IF NOT EXISTS lines (line_id INTEGER PRIMARY KEY, mz_mean REAL NOT NULL, mz_min REAL NOT NULL, mz_max REAL NOT NULL, rt_start REAL NOT NULL, rt_end REAL NOT NULL, ms1_start INTEGER NOT NULL, ms1_end INTEGER NOT NULL, n_points INTEGER NOT NULL);
CREATE TABLE IF NOT EXISTS features (feature_id INTEGER PRIMARY KEY, line_id INTEGER NOT NULL, mz_mean REAL NOT NULL, mz_min REAL NOT NULL, mz_max REAL NOT NULL, rt_start REAL NOT NULL, rt_apex REAL NOT NULL, rt_end REAL NOT NULL, ms1_start INTEGER NOT NULL, ms1_apex INTEGER NOT NULL, ms1_end INTEGER NOT NULL, height REAL NOT NULL, area REAL NOT NULL, n_points INTEGER NOT NULL, quality REAL NOT NULL);
CREATE TABLE IF NOT EXISTS feature_traces (feature_id INTEGER PRIMARY KEY, n INTEGER NOT NULL, scans BLOB NOT NULL, rts BLOB NOT NULL, mzs BLOB NOT NULL, intensities BLOB NOT NULL);
CREATE TABLE IF NOT EXISTS distributions (distribution_id INTEGER PRIMARY KEY, charge INTEGER NOT NULL, neutral_mass REAL NOT NULL, mono_mz REAL NOT NULL, rt_start REAL NOT NULL, rt_apex REAL NOT NULL, rt_end REAL NOT NULL, ms1_start INTEGER NOT NULL, ms1_apex INTEGER NOT NULL, ms1_end INTEGER NOT NULL, n_members INTEGER NOT NULL, score REAL NOT NULL, quality REAL NOT NULL, mz_score REAL NOT NULL DEFAULT 0, iso_score REAL NOT NULL DEFAULT 0, trace_score REAL NOT NULL DEFAULT 0, missing_score REAL NOT NULL DEFAULT 0, interloper_score REAL NOT NULL DEFAULT 0, mono_offset INTEGER NOT NULL DEFAULT 0, n_missing_interior INTEGER NOT NULL DEFAULT 0, n_interlopers INTEGER NOT NULL DEFAULT 0, ambiguity_score REAL NOT NULL DEFAULT 0, status TEXT NOT NULL DEFAULT 'validated');
CREATE TABLE IF NOT EXISTS distribution_members (distribution_id INTEGER NOT NULL, feature_id INTEGER NOT NULL, isotope_index INTEGER NOT NULL, member_score REAL NOT NULL, mz_residual REAL NOT NULL DEFAULT 0, intensity_observed REAL NOT NULL DEFAULT 0, intensity_expected REAL NOT NULL DEFAULT 0, trace_score REAL NOT NULL DEFAULT 0, PRIMARY KEY (distribution_id, feature_id));
CREATE TABLE IF NOT EXISTS analytes (analyte_id INTEGER PRIMARY KEY, neutral_mass REAL NOT NULL, rt_start REAL NOT NULL, rt_apex REAL NOT NULL, rt_end REAL NOT NULL, ms1_start INTEGER NOT NULL, ms1_apex INTEGER NOT NULL, ms1_end INTEGER NOT NULL, charge_min INTEGER NOT NULL, charge_max INTEGER NOT NULL, n_distributions INTEGER NOT NULL, score REAL NOT NULL);
CREATE TABLE IF NOT EXISTS analyte_members (analyte_id INTEGER NOT NULL, distribution_id INTEGER NOT NULL, charge INTEGER NOT NULL, PRIMARY KEY (analyte_id, distribution_id));
"#;

const INDEXES: &str = r#"
CREATE INDEX IF NOT EXISTS idx_scans_rt ON scans(rt);
CREATE INDEX IF NOT EXISTS idx_lines_mz ON lines(mz_mean);
CREATE INDEX IF NOT EXISTS idx_lines_rt ON lines(rt_start, rt_end);
CREATE INDEX IF NOT EXISTS idx_lines_ms1 ON lines(ms1_start, ms1_end);
CREATE INDEX IF NOT EXISTS idx_features_mz ON features(mz_mean);
CREATE INDEX IF NOT EXISTS idx_features_rt ON features(rt_apex);
CREATE INDEX IF NOT EXISTS idx_features_window ON features(rt_start, rt_end);
CREATE INDEX IF NOT EXISTS idx_features_ms1 ON features(ms1_start, ms1_end);
CREATE INDEX IF NOT EXISTS idx_distributions_mass ON distributions(neutral_mass);
CREATE INDEX IF NOT EXISTS idx_distributions_rt ON distributions(rt_apex);
CREATE INDEX IF NOT EXISTS idx_distributions_window ON distributions(rt_start, rt_end);
CREATE INDEX IF NOT EXISTS idx_distributions_ms1 ON distributions(ms1_start, ms1_end);
CREATE INDEX IF NOT EXISTS idx_distributions_charge ON distributions(charge);
CREATE INDEX IF NOT EXISTS idx_distribution_members_dist ON distribution_members(distribution_id);
CREATE INDEX IF NOT EXISTS idx_distribution_members_feature ON distribution_members(feature_id);
CREATE INDEX IF NOT EXISTS idx_analytes_mass ON analytes(neutral_mass);
CREATE INDEX IF NOT EXISTS idx_analytes_rt ON analytes(rt_apex);
CREATE INDEX IF NOT EXISTS idx_analytes_window ON analytes(rt_start, rt_end);
CREATE INDEX IF NOT EXISTS idx_analyte_members_analyte ON analyte_members(analyte_id);
CREATE INDEX IF NOT EXISTS idx_analyte_members_distribution ON analyte_members(distribution_id);
"#;

#[allow(clippy::too_many_arguments)]
pub fn write_db(
    path: &str,
    overwrite: bool,
    params_json: &[(String, String)],
    scans: &[ScanRow],
    lines: &[LineRow],
    features: &[FeatureRow],
    feature_traces: &[FeatureTraceRow],
    dists: &[DistRow],
    analytes: &[AnalyteRow],
    points: &[ScanPoints],
) -> Result<()> {
    if overwrite && std::path::Path::new(path).exists() {
        std::fs::remove_file(path)?;
    }
    let mut conn = Connection::open(path)?;
    conn.execute_batch(
        "PRAGMA journal_mode=WAL; PRAGMA synchronous=NORMAL; PRAGMA temp_store=MEMORY; PRAGMA foreign_keys=ON;",
    )?;
    conn.execute_batch(SCHEMA)?;

    let tx = conn.transaction()?;
    {
        let mut s = tx.prepare("INSERT OR REPLACE INTO parameters(key,value) VALUES (?,?)")?;
        for (k, v) in params_json {
            s.execute(params![k, v])?;
        }
        let mut s = tx.prepare("INSERT INTO scans(ms1_index,spectrum_index,scan_id,rt,tic,n_points) VALUES (?,?,?,?,?,?)")?;
        for r in scans {
            s.execute(params![r.ms1_index, r.spectrum_index, r.scan_id, r.rt, r.tic, r.n_points])?;
        }
        let mut s = tx.prepare("INSERT INTO scan_points(ms1_index,n,mz,intensity) VALUES (?,?,?,?)")?;
        for p in points {
            s.execute(params![p.ms1_index, p.n, p.mz_blob, p.intensity_blob])?;
        }
        let mut s = tx.prepare("INSERT INTO lines(line_id,mz_mean,mz_min,mz_max,rt_start,rt_end,ms1_start,ms1_end,n_points) VALUES (?,?,?,?,?,?,?,?,?)")?;
        for r in lines {
            s.execute(params![r.line_id, r.mz_mean, r.mz_min, r.mz_max, r.rt_start, r.rt_end, r.ms1_start, r.ms1_end, r.n_points])?;
        }
        let mut s = tx.prepare("INSERT INTO features(feature_id,line_id,mz_mean,mz_min,mz_max,rt_start,rt_apex,rt_end,ms1_start,ms1_apex,ms1_end,height,area,n_points,quality) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)")?;
        for r in features {
            s.execute(params![r.feature_id, r.line_id, r.mz_mean, r.mz_min, r.mz_max, r.rt_start, r.rt_apex, r.rt_end, r.ms1_start, r.ms1_apex, r.ms1_end, r.height, r.area, r.n_points, r.quality])?;
        }
        let mut s = tx.prepare("INSERT INTO feature_traces(feature_id,n,scans,rts,mzs,intensities) VALUES (?,?,?,?,?,?)")?;
        for t in feature_traces {
            s.execute(params![t.feature_id, t.n, t.scans_blob, t.rts_blob, t.mzs_blob, t.intensities_blob])?;
        }
        let mut s = tx.prepare("INSERT INTO distributions(distribution_id,charge,neutral_mass,mono_mz,rt_start,rt_apex,rt_end,ms1_start,ms1_apex,ms1_end,n_members,score,quality,mz_score,iso_score,trace_score,missing_score,interloper_score,mono_offset,n_missing_interior,n_interlopers,ambiguity_score,status) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)")?;
        for r in dists {
            s.execute(params![r.distribution_id, r.charge, r.neutral_mass, r.mono_mz, r.rt_start, r.rt_apex, r.rt_end, r.ms1_start, r.ms1_apex, r.ms1_end, r.n_members, r.score, r.quality, r.mz_score, r.iso_score, r.trace_score, r.missing_score, r.interloper_score, r.mono_offset, r.n_missing_interior, r.n_interlopers, r.ambiguity_score, r.status])?;
        }
        let mut s = tx.prepare("INSERT INTO distribution_members(distribution_id,feature_id,isotope_index,member_score,mz_residual,intensity_observed,intensity_expected,trace_score) VALUES (?,?,?,?,?,?,?,?)")?;
        for r in dists {
            for m in &r.members {
                s.execute(params![r.distribution_id, m.feature_id as i64, m.isotope_index, m.member_score, m.mz_residual, m.intensity_observed, m.intensity_expected, m.trace_score])?;
            }
        }
        let mut s = tx.prepare("INSERT INTO analytes(analyte_id,neutral_mass,rt_start,rt_apex,rt_end,ms1_start,ms1_apex,ms1_end,charge_min,charge_max,n_distributions,score) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)")?;
        for a in analytes {
            s.execute(params![a.analyte_id, a.neutral_mass, a.rt_start, a.rt_apex, a.rt_end, a.ms1_start, a.ms1_apex, a.ms1_end, a.charge_min, a.charge_max, a.n_distributions, a.score])?;
        }
        let mut s = tx.prepare("INSERT INTO analyte_members(analyte_id,distribution_id,charge) VALUES (?,?,?)")?;
        for a in analytes {
            for (did, charge) in &a.members {
                s.execute(params![a.analyte_id, did, charge])?;
            }
        }
    }
    tx.commit()?;
    conn.execute_batch(INDEXES)?;
    Ok(())
}
