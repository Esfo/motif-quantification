use std::io::{self, BufRead, Write};
use ms1_detector::config::Config;
use ms1_detector::linemodel::LineModel;

fn main() {
    // each stdin line: "RT mz:int mz:int ..."  (ms1_index = line number)
    let mut model = LineModel::new(Config::default());
    let stdin = io::stdin();
    for (idx, line) in stdin.lock().lines().enumerate() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let rt: f64 = match it.next() { Some(t) => t.parse().unwrap(), None => continue };
        let mut mzs = vec![];
        let mut ints = vec![];
        for tok in it {
            let mut p = tok.split(':');
            let mz: f64 = p.next().unwrap().parse().unwrap();
            let int: f64 = p.next().unwrap().parse().unwrap();
            mzs.push(mz); ints.push(int);
        }
        model.process_scan(idx as i64, rt, &mzs, &ints);
    }
    model.finalize();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    let mut rows: Vec<String> = model.features.iter().map(|f| {
        format!("{:.5} {:.5} {:.3} {:.3} {}", f.mz_mean, f.rt_apex, f.area, f.height, f.n_points)
    }).collect();
    rows.sort();
    writeln!(out, "LINES {} FEATURES {}", model.lines.len(), model.features.len()).unwrap();
    for r in rows { writeln!(out, "{r}").unwrap(); }
}
