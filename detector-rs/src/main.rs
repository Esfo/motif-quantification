use anyhow::Result;
use mzdata::prelude::*;
use mzdata::io::MzMLReader;

fn main() -> Result<()> {
    let path = std::env::args().nth(1).expect("mzml path");
    let reader = MzMLReader::open_path(&path)?;
    let mut ms1 = 0usize;
    let mut first_points = 0usize;
    for spec in reader {
        if spec.ms_level() != 1 { continue; }
        let _rt: f64 = spec.start_time();
        let arrays = match spec.arrays.as_ref() { Some(a) => a, None => continue };
        let mzs = arrays.mzs()?;
        let ints = arrays.intensities()?;
        if ms1 == 0 { first_points = mzs.len(); }
        let _ = (&mzs, &ints);
        ms1 += 1;
    }
    println!("ms1_scans={ms1} first_points={first_points}");
    Ok(())
}
