#!/usr/bin/env python3
"""Generate a small synthetic MS1 mzML with known isotope envelopes.

Used as a deterministic fixture for the envelope-first distribution builder:
we plant peptides of known monoisotopic mass, charge and elution apex, render
their averagine isotope envelopes across a Gaussian chromatographic peak, add a
little noise, and write centroided MS1 scans. Because the ground truth is known
(see the printed manifest), it lets us check that the detector recovers the
right charge — in particular that 2+ dominates when the planted set is 2+-heavy.

No randomness primitives that the harness forbids are used at import time; a
fixed seed is taken from the CLI so runs are reproducible.
"""
import argparse
import math

import numpy as np
from psims.mzml import MzMLWriter

PROTON = 1.007276554940804
C13 = 1.00335483507
AVERAGINE_LAMBDA_PER_DA = 0.000594  # Poisson mean of the averagine isotope model


def averagine_intensities(neutral_mass, n_peaks):
    lam = neutral_mass * AVERAGINE_LAMBDA_PER_DA
    ks = np.arange(n_peaks)
    # Poisson pmf without scipy: exp(-lam) lam^k / k!
    logp = -lam + ks * math.log(lam) - np.array([math.lgamma(k + 1) for k in ks])
    p = np.exp(logp)
    return p / p.max()


def gaussian(x, mu, sigma):
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def build(seed, n_scans):
    rng = np.random.default_rng(seed)
    rts = np.linspace(0.0, n_scans * 0.01, n_scans)  # ~0.6s spacing

    # (neutral_mass, charge, apex_scan, peak_intensity, n_isotopes, rt_sigma_scans)
    # 2+-heavy planted set so the correct histogram is 2+-dominant.
    plant = []
    masses_2 = [1200, 1450, 1600, 1850, 2100, 2400, 980, 1320, 1750, 2250, 1100, 1550]
    for i, m in enumerate(masses_2):
        plant.append((float(m), 2, 30 + i * 18, 1.0e6, 6, 2.2))
    for i, m in enumerate([2600, 3000, 3400, 2800]):
        plant.append((float(m), 3, 60 + i * 40, 7.0e5, 7, 2.0))
    for i, m in enumerate([4200, 4800]):
        plant.append((float(m), 4, 90 + i * 70, 5.0e5, 8, 1.8))
    plant.append((900.0, 1, 120, 4.0e5, 3, 2.0))

    # Faint, short envelopes near the detection floor: real distributions whose
    # member lines are weak/short, meant to be missed by the strict primary pass
    # and picked up by the relaxed recovery pass.
    for i, m in enumerate([1380, 1700, 2050]):
        plant.append((float(m), 2, 70 + i * 55, 2.2e4, 4, 1.1))

    manifest = []
    # scan_index -> list of (mz, intensity)
    peaks_per_scan = [[] for _ in range(n_scans)]

    for (mass, z, apex, amp, n_iso, sigma) in plant:
        mono_mz = (mass + z * PROTON) / z
        iso = averagine_intensities(mass, n_iso)
        manifest.append((mass, z, apex, mono_mz))
        for scan in range(n_scans):
            shape = gaussian(scan, apex, sigma)
            if shape < 0.02:
                continue
            for k in range(n_iso):
                mz = mono_mz + k * C13 / z
                inten = amp * shape * iso[k]
                inten *= 1.0 + 0.03 * rng.standard_normal()
                # small m/z jitter (centroiding precision)
                mz += 0.0008 * rng.standard_normal()
                if inten > 1.0e3:
                    peaks_per_scan[scan].append((mz, inten))

    # background noise peaks
    for scan in range(n_scans):
        for _ in range(40):
            mz = 400 + rng.random() * 1600
            inten = rng.random() * 5.0e3
            peaks_per_scan[scan].append((float(mz), float(inten)))

    return rts, peaks_per_scan, manifest


def write_mzml(path, rts, peaks_per_scan):
    with MzMLWriter(open(path, "wb")) as writer:
        writer.controlled_vocabularies()
        with writer.run(id="synth"):
            n = len(peaks_per_scan)
            with writer.spectrum_list(count=n):
                for i, (rt, peaks) in enumerate(zip(rts, peaks_per_scan)):
                    peaks = sorted(peaks)
                    mzs = np.array([p[0] for p in peaks], dtype=np.float64)
                    its = np.array([p[1] for p in peaks], dtype=np.float64)
                    writer.write_spectrum(
                        mzs, its, id=f"scan={i+1}", centroided=True,
                        scan_start_time=rt, params=[{"ms level": 1}],
                    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--scans", type=int, default=320)
    args = ap.parse_args()
    rts, peaks, manifest = build(args.seed, args.scans)
    write_mzml(args.out, rts, peaks)
    print(f"wrote {args.out} scans={args.scans} planted={len(manifest)}")
    from collections import Counter
    c = Counter(z for (_, z, _, _) in manifest)
    print("planted charge histogram:", dict(sorted(c.items())))
    for mass, z, apex, mono in manifest:
        print(f"  mass={mass:.1f} z={z} apex_scan={apex} mono_mz={mono:.4f}")


if __name__ == "__main__":
    main()
