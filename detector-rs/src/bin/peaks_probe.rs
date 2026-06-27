use std::io::{self, BufRead, Write};
use ms1_detector::peaks::axis_peaks;

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let arr: Vec<f64> = line.split_whitespace().filter_map(|t| t.parse().ok()).collect();
        let peaks = axis_peaks(&arr, 2);
        let s: Vec<String> = peaks.iter().map(|(l,a,r)| format!("{l},{a},{r}")).collect();
        writeln!(out, "{}", s.join(";")).unwrap();
    }
}
