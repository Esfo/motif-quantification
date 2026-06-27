//! Detector configuration — mirrors the Config dataclass in index_ml1.py.

pub const PROTON: f64 = 1.007276554940804;
pub const C13_DELTA: f64 = 1.00335483507;

#[derive(Clone, Debug)]
// Some fields mirror the Python Config for parity but are not read by the Rust
// path (valley merge disabled by default; step/new_inc replaced by the ratio gate).
#[allow(dead_code)]
pub struct Config {
    // line model
    pub line_mz_ppm: f64,
    pub line_mz_abs: f64,
    pub deadsignal: i64,
    pub line_merge_mz_ppm: f64,
    pub line_merge_mz_abs: f64,
    pub line_merge_gap_scans: i64,
    pub line_split_sigma: f64,
    pub min_split_valley_fraction: f64,
    pub min_trace_points: usize,
    pub peak_mindist: usize,
    pub smooth_points: usize,
    pub min_peak_points: usize,
    pub min_peak_height: f64,
    pub min_peak_area: f64,
    pub min_peak_width: f64,
    pub max_peak_width: f64,
    pub min_peak_prominence_fraction: f64,
    pub max_trace_peaks: usize,
    // isotope / distributions
    pub isotope_mz_ppm: f64,
    pub isotope_mz_abs: f64,
    pub max_neutral_mass: f64,
    pub max_apex_shift: f64,
    pub max_apex_shift_width_fraction: f64,
    pub min_edge_score: f64,
    pub min_distribution_members: usize,
    pub min_members_charge_one: usize,
    pub charge_one_score_penalty: f64,
    pub max_adjacent_intensity_ratio: f64,
    pub charge_tolerance: f64,
    pub mass_width_limit: f64,
    pub step_limit: f64,
    pub new_inc_limit: f64,
    // charge grouping
    pub charge_mass_ppm: f64,
    pub min_charge_group_rt_score: f64,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            line_mz_ppm: 8.0,
            line_mz_abs: 0.002,
            deadsignal: 3,
            line_merge_mz_ppm: 10.0,
            line_merge_mz_abs: 0.004,
            line_merge_gap_scans: 6,
            line_split_sigma: 2.0,
            min_split_valley_fraction: 0.0,
            min_trace_points: 4,
            peak_mindist: 2,
            smooth_points: 3,
            min_peak_points: 4,
            min_peak_height: 0.0,
            min_peak_area: 0.0,
            min_peak_width: 0.0,
            max_peak_width: 6.0,
            min_peak_prominence_fraction: 0.02,
            max_trace_peaks: 0,
            isotope_mz_ppm: 10.0,
            isotope_mz_abs: 0.004,
            max_neutral_mass: 8000.0,
            max_apex_shift: 0.15,
            max_apex_shift_width_fraction: 0.50,
            min_edge_score: 0.30,
            min_distribution_members: 2,
            min_members_charge_one: 3,
            charge_one_score_penalty: 0.85,
            max_adjacent_intensity_ratio: 10.0,
            charge_tolerance: 0.1,
            mass_width_limit: 0.002,
            step_limit: 0.5,
            new_inc_limit: 0.1,
            charge_mass_ppm: 12.0,
            min_charge_group_rt_score: 0.10,
        }
    }
}
