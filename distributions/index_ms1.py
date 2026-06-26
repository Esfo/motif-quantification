#!/usr/bin/env python3

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from time import time

import numpy as np
from pyteomics import mzml

try:
    from .peaks import axis_peaks, moving_average
    from .store import connect_db, create_indexes, init_schema, write_parameters, write_rows
except ImportError:
    from peaks import axis_peaks, moving_average
    from store import connect_db, create_indexes, init_schema, write_parameters, write_rows

try:
    from .elementalcomponents import proton
except ImportError:
    try:
        from elementalcomponents import proton
    except ImportError:
        proton = 1.007276554940804


C13_DELTA = 1.00335483507
SCRIPT_VERSION = "0.3.0"

_EDGE = {}


@dataclass(slots=True)
class Config:
    line_mz_ppm: float = 8.0
    line_mz_abs: float = 0.002
    max_gap_scans: int = 2

    min_trace_points: int = 4
    peak_mindist: int = 2
    smooth_points: int = 3

    min_peak_points: int = 4
    min_peak_height: float = 0.0
    min_peak_area: float = 0.0
    min_peak_width: float = 0.0
    max_peak_width: float = 6.0
    min_peak_prominence_fraction: float = 0.02
    max_trace_peaks: int = 0

    max_charge: int = 6
    isotope_mz_ppm: float = 10.0
    isotope_mz_abs: float = 0.004
    max_neutral_mass: float = 8000.0

    max_apex_shift: float = 0.15
    max_apex_shift_width_fraction: float = 0.50
    min_edge_score: float = 0.30
    min_distribution_members: int = 2

    charge_mass_ppm: float = 12.0
    min_charge_group_rt_score: float = 0.10

    # Reference (distributionassembly.py) isotope-edge acceptance, ported faithfully:
    # asymmetric acdiff tolerance around proton-spacing, plus intensity-step gating.
    # mass_width_limit stands in for the reference masswidthlimit (=roundcutoff*2),
    # which this pipeline does not track per-scan; tune if needed.
    charge_tolerance: float = 0.1
    mass_width_limit: float = 0.002
    step_limit: float = 0.5
    new_inc_limit: float = 0.1


@dataclass(slots=True)
class Trace:
    line_id: int
    scans: list
    rts: list
    mzs: list
    intensities: list
    mean_mz: float
    min_mz: float
    max_mz: float
    last_scan: int

    @classmethod
    def create(cls, line_id, ms1_index, rt, mz_value, intensity):
        mz_value = float(mz_value)
        intensity = float(intensity)

        return cls(
            line_id=int(line_id),
            scans=[int(ms1_index)],
            rts=[float(rt)],
            mzs=[mz_value],
            intensities=[intensity],
            mean_mz=mz_value,
            min_mz=mz_value,
            max_mz=mz_value,
            last_scan=int(ms1_index),
        )

    def append(self, ms1_index, rt, mz_value, intensity):
        mz_value = float(mz_value)
        intensity = float(intensity)
        n = len(self.mzs)

        self.mean_mz = (self.mean_mz * n + mz_value) / (n + 1)
        self.min_mz = min(self.min_mz, mz_value)
        self.max_mz = max(self.max_mz, mz_value)
        self.last_scan = int(ms1_index)

        self.scans.append(int(ms1_index))
        self.rts.append(float(rt))
        self.mzs.append(mz_value)
        self.intensities.append(intensity)


@dataclass(slots=True)
class Feature:
    feature_id: int
    line_id: int
    mz_mean: float
    mz_min: float
    mz_max: float
    rt_start: float
    rt_apex: float
    rt_end: float
    ms1_start: int
    ms1_apex: int
    ms1_end: int
    height: float
    area: float
    n_points: int
    quality: float


@dataclass(slots=True)
class Distribution:
    distribution_id: int
    charge: int
    neutral_mass: float
    mono_mz: float
    rt_start: float
    rt_apex: float
    rt_end: float
    ms1_start: int
    ms1_apex: int
    ms1_end: int
    n_members: int
    score: float
    quality: float
    members: list


@dataclass(slots=True)
class Analyte:
    analyte_id: int
    neutral_mass: float
    rt_start: float
    rt_apex: float
    rt_end: float
    ms1_start: int
    ms1_apex: int
    ms1_end: int
    charge_min: int
    charge_max: int
    n_distributions: int
    score: float
    members: list


class LineModel:
    def __init__(self, config):
        self.config = config
        self.active = {}
        self.next_line_id = 0
        self.next_feature_id = 0
        self.lines = []
        self.features = []

    def mz_tolerance(self, mz_value):
        return max(
            self.config.line_mz_abs,
            mz_value * self.config.line_mz_ppm / 1_000_000.0,
        )

    def process_scan(self, ms1_index, rt, mzs, intensities):
        mzs = np.asarray(mzs, dtype=np.float64).reshape(-1)
        intensities = np.asarray(intensities, dtype=np.float64).reshape(-1)

        if mzs.size == 0:
            self.close_dead_lines(ms1_index)
            return

        order = np.argsort(mzs)
        mzs = mzs[order]
        intensities = intensities[order]

        active_ids = np.fromiter(self.active.keys(), dtype=np.int64, count=len(self.active))

        if active_ids.size:
            active_mz = np.array(
                [self.active[int(line_id)].mean_mz for line_id in active_ids],
                dtype=np.float64,
            )
            active_order = np.argsort(active_mz)
            active_ids = active_ids[active_order]
            active_mz = active_mz[active_order]
        else:
            active_mz = np.empty(0, dtype=np.float64)

        used_lines = set()

        for mz_value, intensity in zip(mzs, intensities):
            mz_value = float(mz_value)
            intensity = float(intensity)
            tolerance = self.mz_tolerance(mz_value)

            best_line_id = None
            best_score = math.inf

            if active_mz.size:
                left = np.searchsorted(active_mz, mz_value - tolerance, side="left")
                right = np.searchsorted(active_mz, mz_value + tolerance, side="right")

                for loc in range(left, right):
                    line_id = int(active_ids[loc])

                    if line_id in used_lines:
                        continue

                    trace = self.active[line_id]
                    gap = ms1_index - trace.last_scan

                    if gap < 1 or gap > self.config.max_gap_scans + 1:
                        continue

                    score = abs(trace.mean_mz - mz_value) / tolerance + gap * 0.03

                    if score < best_score:
                        best_score = score
                        best_line_id = line_id

            if best_line_id is None:
                line_id = self.next_line_id
                self.next_line_id += 1
                self.active[line_id] = Trace.create(line_id, ms1_index, rt, mz_value, intensity)
                used_lines.add(line_id)
            else:
                self.active[best_line_id].append(ms1_index, rt, mz_value, intensity)
                used_lines.add(best_line_id)

        self.close_dead_lines(ms1_index)

    def close_dead_lines(self, ms1_index):
        dead = []

        for line_id, trace in self.active.items():
            if ms1_index - trace.last_scan > self.config.max_gap_scans:
                dead.append(line_id)

        for line_id in dead:
            trace = self.active.pop(line_id)
            self.finalize_trace(trace)

    def close_all(self):
        for line_id in sorted(self.active):
            self.finalize_trace(self.active[line_id])

        self.active.clear()

    def finalize_trace(self, trace):
        if len(trace.scans) < self.config.min_trace_points:
            return

        scans = np.asarray(trace.scans, dtype=np.int32)
        rts = np.asarray(trace.rts, dtype=np.float64)
        mzs = np.asarray(trace.mzs, dtype=np.float64)
        intensities = np.asarray(trace.intensities, dtype=np.float64)

        self.lines.append(
            {
                "line_id": trace.line_id,
                "mz_mean": float(mzs.mean()),
                "mz_min": float(mzs.min()),
                "mz_max": float(mzs.max()),
                "rt_start": float(rts.min()),
                "rt_end": float(rts.max()),
                "ms1_start": int(scans.min()),
                "ms1_end": int(scans.max()),
                "n_points": int(scans.size),
            }
        )

        self.split_trace(trace.line_id, scans, rts, mzs, intensities)

    def split_trace(self, line_id, scans, rts, mzs, intensities):
        if intensities.size < self.config.min_trace_points:
            return

        smoothed = moving_average(intensities, self.config.smooth_points)
        peaks = axis_peaks(smoothed, mindist=self.config.peak_mindist)

        if not peaks:
            return

        trace_height = float(intensities.max())
        candidate_peaks = []

        for left, apex, right in peaks:
            left = int(left)
            apex = int(apex)
            right = int(right)

            if right <= left:
                continue

            sub_count = right - left

            if sub_count < self.config.min_peak_points:
                continue

            apex_height = float(intensities[apex])
            edge_height = float(max(intensities[left], intensities[right - 1]))
            prominence = apex_height - edge_height

            if apex_height < self.config.min_peak_height:
                continue

            if (
                self.config.min_peak_prominence_fraction > 0
                and trace_height > 0
                and prominence < trace_height * self.config.min_peak_prominence_fraction
            ):
                continue

            candidate_peaks.append((apex_height, left, apex, right))

        if not candidate_peaks:
            return

        candidate_peaks.sort(reverse=True)

        if self.config.max_trace_peaks > 0:
            candidate_peaks = candidate_peaks[: self.config.max_trace_peaks]

        candidate_peaks.sort(key=lambda x: x[1])

        for _, left, apex, right in candidate_peaks:
            sub_scans = scans[left:right]
            sub_rts = rts[left:right]
            sub_mzs = mzs[left:right]
            sub_intensities = intensities[left:right]

            if sub_scans.size < self.config.min_peak_points:
                continue

            width = float(sub_rts.max() - sub_rts.min())

            if width < self.config.min_peak_width:
                continue

            if self.config.max_peak_width > 0 and width > self.config.max_peak_width:
                continue

            if sub_rts.size > 1:
                area = float(np.trapezoid(sub_intensities, sub_rts))
            else:
                area = float(sub_intensities[0])

            if area < self.config.min_peak_area:
                continue

            height = float(sub_intensities.max())

            if height < self.config.min_peak_height:
                continue

            local_apex = int(sub_intensities.argmax())
            total_intensity = float(sub_intensities.sum())

            if total_intensity > 0:
                mz_mean = float((sub_mzs * sub_intensities).sum() / total_intensity)
            else:
                mz_mean = float(sub_mzs.mean())

            quality = float(np.log1p(max(area, 0.0)) * math.sqrt(sub_scans.size))

            self.features.append(
                Feature(
                    feature_id=self.next_feature_id,
                    line_id=int(line_id),
                    mz_mean=mz_mean,
                    mz_min=float(sub_mzs.min()),
                    mz_max=float(sub_mzs.max()),
                    rt_start=float(sub_rts.min()),
                    rt_apex=float(sub_rts[local_apex]),
                    rt_end=float(sub_rts.max()),
                    ms1_start=int(sub_scans.min()),
                    ms1_apex=int(sub_scans[local_apex]),
                    ms1_end=int(sub_scans.max()),
                    height=height,
                    area=area,
                    n_points=int(sub_scans.size),
                    quality=quality,
                )
            )

            self.next_feature_id += 1


def scan_id(scan):
    return scan.get("id") or f"index={scan.get('index', 0)}"


def scan_rt(scan):
    try:
        return float(np.real(scan["scanList"]["scan"][0]["scan start time"]))
    except Exception:
        return None


def scan_ms_level(scan):
    return int(scan.get("ms level", 0) or 0)


def stream_ms1(mzml_path, min_intensity):
    ms1_index = 0

    with mzml.MzML(str(mzml_path), dtype=np.float64) as reader:
        for spectrum_index, scan in enumerate(reader):
            if scan_ms_level(scan) != 1:
                continue

            rt = scan_rt(scan)

            if rt is None:
                continue

            mzs = np.asarray(scan.get("m/z array", []), dtype=np.float64)
            intensities = np.asarray(scan.get("intensity array", []), dtype=np.float64)

            if mzs.size != intensities.size:
                raise ValueError(
                    f"{scan_id(scan)}: m/z and intensity arrays differ: {mzs.size} != {intensities.size}"
                )

            if min_intensity > 0:
                keep = intensities >= min_intensity
                mzs = mzs[keep]
                intensities = intensities[keep]

            yield {
                "ms1_index": ms1_index,
                "spectrum_index": int(spectrum_index),
                "scan_id": scan_id(scan),
                "rt": float(rt),
                "mzs": mzs,
                "intensities": intensities,
                "tic": float(intensities.sum()),
                "n_points": int(mzs.size),
            }

            ms1_index += 1


def feature_rows(features):
    for feature in features:
        yield {
            "feature_id": feature.feature_id,
            "line_id": feature.line_id,
            "mz_mean": feature.mz_mean,
            "mz_min": feature.mz_min,
            "mz_max": feature.mz_max,
            "rt_start": feature.rt_start,
            "rt_apex": feature.rt_apex,
            "rt_end": feature.rt_end,
            "ms1_start": feature.ms1_start,
            "ms1_apex": feature.ms1_apex,
            "ms1_end": feature.ms1_end,
            "height": feature.height,
            "area": feature.area,
            "n_points": feature.n_points,
            "quality": feature.quality,
        }


def _set_edge_context(features, config):
    global _EDGE

    expected_ids = list(range(len(features)))

    actual_ids = [feature.feature_id for feature in features]

    if actual_ids != expected_ids:
        raise RuntimeError("feature ids must be dense zero-based integers for this version")

    _EDGE = {
        "feature_id": np.asarray(actual_ids, dtype=np.int64),
        "mz_mean": np.asarray([f.mz_mean for f in features], dtype=np.float64),
        "rt_start": np.asarray([f.rt_start for f in features], dtype=np.float64),
        "rt_apex": np.asarray([f.rt_apex for f in features], dtype=np.float64),
        "rt_end": np.asarray([f.rt_end for f in features], dtype=np.float64),
        "ms1_start": np.asarray([f.ms1_start for f in features], dtype=np.int64),
        "ms1_apex": np.asarray([f.ms1_apex for f in features], dtype=np.int64),
        "ms1_end": np.asarray([f.ms1_end for f in features], dtype=np.int64),
        "height": np.asarray([f.height for f in features], dtype=np.float64),
        "area": np.asarray([f.area for f in features], dtype=np.float64),
        "n_points": np.asarray([f.n_points for f in features], dtype=np.int64),
        "quality": np.asarray([f.quality for f in features], dtype=np.float64),
        "config": asdict(config),
    }

    order = np.argsort(_EDGE["mz_mean"])
    _EDGE["order"] = order
    _EDGE["sorted_mz"] = _EDGE["mz_mean"][order]


def _rt_overlap(left_i, right_i):
    rt_start = _EDGE["rt_start"]
    rt_end = _EDGE["rt_end"]

    overlap = min(rt_end[left_i], rt_end[right_i]) - max(rt_start[left_i], rt_start[right_i])
    union = max(rt_end[left_i], rt_end[right_i]) - min(rt_start[left_i], rt_start[right_i])

    if union <= 0:
        return 0.0

    return max(0.0, float(overlap / union))


def _intensity_score(left_i, right_i):
    left_height = _EDGE["height"][left_i]
    right_height = _EDGE["height"][right_i]

    if left_height <= 0 or right_height <= 0:
        return 0.0

    ratio = right_height / left_height
    score = 1.0 - min(1.0, abs(math.log2(ratio)) / 4.0)

    if right_height > left_height * 4.0:
        score *= 0.75

    return max(0.0, min(1.0, float(score)))


def _score_edge(left_i, right_i, charge):
    cfg = _EDGE["config"]
    mz = _EDGE["mz_mean"]
    rt_start = _EDGE["rt_start"]
    rt_apex = _EDGE["rt_apex"]
    rt_end = _EDGE["rt_end"]

    # Faithful to the reference (distributionassembly.py:226-231): acdiff is the
    # signed deviation of the observed isotope spacing from proton/charge, with an
    # ASYMMETRIC tolerance scaled by charge_tolerance (+ mass_width_limit).
    expdiff = proton / charge
    observed = mz[right_i] - mz[left_i]
    acdiff = expdiff - observed
    diffcut = expdiff * cfg["charge_tolerance"]
    masswidthlimit = cfg["mass_width_limit"]
    lower = -(diffcut * cfg["charge_tolerance"] + masswidthlimit)
    upper = diffcut + masswidthlimit

    if not (acdiff > lower and acdiff <= upper):
        return None

    mz_error = -acdiff
    tolerance = upper if upper > 0 else max(cfg["isotope_mz_abs"], 1e-9)

    neutral_mass = (mz[left_i] - proton) * charge

    if neutral_mass <= 0 or neutral_mass > cfg["max_neutral_mass"]:
        return None

    # Reference intensity-step gating (distributionassembly.py:240-278): the
    # fractional intensity difference between adjacent isotopes is gated more
    # strictly when intensity is increasing (new_inc_limit) than decreasing
    # (step_limit).
    left_h = _EDGE["height"][left_i]
    right_h = _EDGE["height"][right_i]
    denom = left_h + right_h
    if denom > 0:
        intensitypercdiff = abs(right_h - left_h) / denom / 2.0
        ratiocheck = cfg["new_inc_limit"] if right_h > left_h else cfg["step_limit"]
        if intensitypercdiff > ratiocheck:
            return None

    left_width = rt_end[left_i] - rt_start[left_i]
    right_width = rt_end[right_i] - rt_start[right_i]
    width = max(left_width, right_width, 1e-9)

    max_shift = max(
        cfg["max_apex_shift"],
        cfg["max_apex_shift_width_fraction"] * width,
    )

    rt_shift = rt_apex[right_i] - rt_apex[left_i]

    if abs(rt_shift) > max_shift:
        return None

    mz_score = 1.0 - min(1.0, abs(mz_error) / tolerance)
    shift_score = 1.0 - min(1.0, abs(rt_shift) / max_shift)
    overlap_score = _rt_overlap(left_i, right_i)
    i_score = _intensity_score(left_i, right_i)

    score = (
        mz_score * 0.50
        + shift_score * 0.20
        + overlap_score * 0.15
        + i_score * 0.15
    )

    if score < cfg["min_edge_score"]:
        return None

    return {
        "left_feature_id": int(left_i),
        "right_feature_id": int(right_i),
        "charge": int(charge),
        "isotope_step": 1,
        "mz_error": float(mz_error),
        "mz_error_ppm": float(mz_error / mz[right_i] * 1_000_000.0),
        "rt_shift": float(rt_shift),
        "rt_overlap": float(overlap_score),
        "intensity_score": float(i_score),
        "score": float(score),
    }


def _edge_worker(start, end, store_edges):
    cfg = _EDGE["config"]
    order = _EDGE["order"]
    sorted_mz = _EDGE["sorted_mz"]
    mz = _EDGE["mz_mean"]

    best_edges = []
    stored_edges = []

    for sorted_left in range(start, end):
        left_i = int(order[sorted_left])
        left_mz = mz[left_i]

        for charge in range(1, cfg["max_charge"] + 1):
            target = left_mz + C13_DELTA / charge

            tolerance = max(
                cfg["isotope_mz_abs"],
                target * cfg["isotope_mz_ppm"] / 1_000_000.0,
            )

            left = np.searchsorted(sorted_mz, target - tolerance, side="left")
            right = np.searchsorted(sorted_mz, target + tolerance, side="right")

            best = None

            for sorted_right in range(left, right):
                right_i = int(order[sorted_right])

                if right_i == left_i:
                    continue

                if mz[right_i] <= left_mz:
                    continue

                edge = _score_edge(left_i, right_i, charge)

                if edge is None:
                    continue

                if store_edges == "all":
                    stored_edges.append(edge)

                if best is None or edge["score"] > best["score"]:
                    best = edge

            if best is not None:
                best_edges.append(best)

                if store_edges == "best":
                    stored_edges.append(best)

    return best_edges, stored_edges


def build_isotope_edges(features, config, workers, chunk_size, store_edges, progress):
    if not features:
        return [], []

    _set_edge_context(features, config)

    n_features = len(features)
    chunks = [
        (start, min(start + chunk_size, n_features))
        for start in range(0, n_features, chunk_size)
    ]

    all_best = []
    all_stored = []

    if workers <= 1:
        for chunk_n, (start, end) in enumerate(chunks, start=1):
            best, stored = _edge_worker(start, end, store_edges)
            all_best.extend(best)
            all_stored.extend(stored)

            if progress:
                print(
                    f"edge_chunks={chunk_n}/{len(chunks)} best_edges={len(all_best)} stored_edges={len(all_stored)}",
                    file=sys.stderr,
                    flush=True,
                )

        return all_best, all_stored

    try:
        context = mp.get_context("fork")
    except ValueError:
        context = None

    if context is None:
        return build_isotope_edges(
            features=features,
            config=config,
            workers=1,
            chunk_size=chunk_size,
            store_edges=store_edges,
            progress=progress,
        )

    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
        futures = [
            executor.submit(_edge_worker, start, end, store_edges)
            for start, end in chunks
        ]

        for chunk_n, future in enumerate(as_completed(futures), start=1):
            best, stored = future.result()
            all_best.extend(best)
            all_stored.extend(stored)

            if progress:
                print(
                    f"edge_chunks={chunk_n}/{len(chunks)} best_edges={len(all_best)} stored_edges={len(all_stored)}",
                    file=sys.stderr,
                    flush=True,
                )

    return all_best, all_stored


def _coherent_rt_path_ids(path, config_dict):
    if len(path) < 3:
        return True

    rt_apex = _EDGE["rt_apex"]
    rt_start = _EDGE["rt_start"]
    rt_end = _EDGE["rt_end"]

    path_array = np.asarray(path, dtype=np.int64)
    xs = np.arange(path_array.size, dtype=np.float64)
    ys = rt_apex[path_array]

    widths = np.maximum(
        rt_end[path_array] - rt_start[path_array],
        1e-9,
    )

    slope, intercept = np.polyfit(xs, ys, 1)
    predicted = slope * xs + intercept
    residual = float(np.max(np.abs(ys - predicted)))

    allowed = max(
        config_dict["max_apex_shift"],
        config_dict["max_apex_shift_width_fraction"] * float(np.median(widths)),
    )

    return residual <= allowed


def _distribution_worker_for_charge(charge, edges, config_dict):
    mz_mean = _EDGE["mz_mean"]
    rt_start = _EDGE["rt_start"]
    rt_apex = _EDGE["rt_apex"]
    rt_end = _EDGE["rt_end"]
    ms1_start = _EDGE["ms1_start"]
    ms1_apex = _EDGE["ms1_apex"]
    ms1_end = _EDGE["ms1_end"]
    height = _EDGE["height"]

    best_out = {}
    incoming = set()

    for edge in sorted(edges, key=lambda x: x["score"], reverse=True):
        left_feature_id = edge["left_feature_id"]
        right_feature_id = edge["right_feature_id"]

        key = (left_feature_id, charge)

        if key not in best_out:
            best_out[key] = edge
            incoming.add((right_feature_id, charge))

    starts = []

    for edge in best_out.values():
        if (edge["left_feature_id"], charge) not in incoming:
            starts.append(edge["left_feature_id"])

    rows = []
    used_paths = set()

    for start_feature_id in sorted(starts):
        path = [int(start_feature_id)]
        path_edges = []
        current = int(start_feature_id)

        while True:
            edge = best_out.get((current, charge))

            if edge is None:
                break

            next_feature_id = int(edge["right_feature_id"])

            if next_feature_id in path:
                break

            tentative = path + [next_feature_id]

            if not _coherent_rt_path_ids(tentative, config_dict):
                break

            path_edges.append(edge)
            path = tentative
            current = next_feature_id

        if len(path) < config_dict["min_distribution_members"]:
            continue

        key = (charge, tuple(path))

        if key in used_paths:
            continue

        used_paths.add(key)

        path_array = np.asarray(path, dtype=np.int64)
        apex_feature_id = int(path_array[np.argmax(height[path_array])])
        mono_feature_id = int(path_array[0])

        score = float(np.mean([edge["score"] for edge in path_edges])) if path_edges else 0.0
        quality = float(score * math.sqrt(len(path)))

        rows.append(
            {
                "charge": int(charge),
                "neutral_mass": float((mz_mean[mono_feature_id] - proton) * charge),
                "mono_mz": float(mz_mean[mono_feature_id]),
                "rt_start": float(np.min(rt_start[path_array])),
                "rt_apex": float(rt_apex[apex_feature_id]),
                "rt_end": float(np.max(rt_end[path_array])),
                "ms1_start": int(np.min(ms1_start[path_array])),
                "ms1_apex": int(ms1_apex[apex_feature_id]),
                "ms1_end": int(np.max(ms1_end[path_array])),
                "n_members": int(len(path)),
                "score": score,
                "quality": quality,
                "members": [
                    {
                        "feature_id": int(feature_id),
                        "isotope_index": int(isotope_index),
                        "member_score": 1.0 if isotope_index == 0 else float(path_edges[isotope_index - 1]["score"]),
                    }
                    for isotope_index, feature_id in enumerate(path)
                ],
            }
        )

    return rows


def build_distributions(features, best_edges, config, workers=1, progress=False):
    if not best_edges:
        return []

    config_dict = asdict(config)
    edges_by_charge = {}

    for edge in best_edges:
        edges_by_charge.setdefault(edge["charge"], []).append(edge)

    jobs = sorted(edges_by_charge.items())
    all_rows = []

    if workers <= 1 or len(jobs) <= 1:
        for job_n, (charge, edges) in enumerate(jobs, start=1):
            rows = _distribution_worker_for_charge(charge, edges, config_dict)
            all_rows.extend(rows)

            if progress:
                print(
                    f"distribution_charge={job_n}/{len(jobs)} charge={charge} rows={len(rows)} total={len(all_rows)}",
                    file=sys.stderr,
                    flush=True,
                )
    else:
        try:
            context = mp.get_context("fork")
        except ValueError:
            context = None

        if context is None:
            return build_distributions(
                features=features,
                best_edges=best_edges,
                config=config,
                workers=1,
                progress=progress,
            )

        max_workers = min(workers, len(jobs))

        with ProcessPoolExecutor(max_workers=max_workers, mp_context=context) as executor:
            futures = {
                executor.submit(_distribution_worker_for_charge, charge, edges, config_dict): charge
                for charge, edges in jobs
            }

            for done_n, future in enumerate(as_completed(futures), start=1):
                charge = futures[future]
                rows = future.result()
                all_rows.extend(rows)

                if progress:
                    print(
                        f"distribution_charge={done_n}/{len(jobs)} charge={charge} rows={len(rows)} total={len(all_rows)}",
                        file=sys.stderr,
                        flush=True,
                    )

    all_rows.sort(
        key=lambda row: (
            row["charge"],
            row["neutral_mass"],
            row["rt_apex"],
            row["mono_mz"],
        )
    )

    distributions = []

    for distribution_id, row in enumerate(all_rows):
        distributions.append(
            Distribution(
                distribution_id=distribution_id,
                charge=row["charge"],
                neutral_mass=row["neutral_mass"],
                mono_mz=row["mono_mz"],
                rt_start=row["rt_start"],
                rt_apex=row["rt_apex"],
                rt_end=row["rt_end"],
                ms1_start=row["ms1_start"],
                ms1_apex=row["ms1_apex"],
                ms1_end=row["ms1_end"],
                n_members=row["n_members"],
                score=row["score"],
                quality=row["quality"],
                members=row["members"],
            )
        )

    return distributions


def distribution_rt_score(a, b):
    overlap = min(a.rt_end, b.rt_end) - max(a.rt_start, b.rt_start)
    union = max(a.rt_end, b.rt_end) - min(a.rt_start, b.rt_start)

    if union <= 0:
        raw = 0.0
    else:
        raw = max(0.0, overlap / union)

    apex_gap = abs(a.rt_apex - b.rt_apex)
    width = max(a.rt_end - a.rt_start, b.rt_end - b.rt_start, 1e-9)
    apex_score = max(0.0, 1.0 - apex_gap / width)

    return max(0.0, min(1.0, raw * 0.4 + apex_score * 0.6))


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, item):
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]

        return item

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)

        if left_root != right_root:
            self.parent[right_root] = left_root


def build_analytes(distributions, config):
    if not distributions:
        return []

    sorted_dists = sorted(distributions, key=lambda x: x.neutral_mass)
    masses = np.asarray([dist.neutral_mass for dist in sorted_dists], dtype=np.float64)
    uf = UnionFind(len(sorted_dists))

    for i, dist in enumerate(sorted_dists):
        tolerance = max(0.002, dist.neutral_mass * config.charge_mass_ppm / 1_000_000.0)
        start = np.searchsorted(masses, dist.neutral_mass - tolerance, side="left")
        end = np.searchsorted(masses, dist.neutral_mass + tolerance, side="right")

        for j in range(start, end):
            if i == j:
                continue

            other = sorted_dists[j]

            if other.charge == dist.charge:
                continue

            if distribution_rt_score(dist, other) >= config.min_charge_group_rt_score:
                uf.union(i, j)

    groups = {}

    for i, dist in enumerate(sorted_dists):
        groups.setdefault(uf.find(i), []).append(dist)

    analytes = []

    for members in groups.values():
        weights = np.asarray([max(member.score, 1e-6) for member in members], dtype=np.float64)
        masses = np.asarray([member.neutral_mass for member in members], dtype=np.float64)
        charges = [member.charge for member in members]
        apex_member = max(members, key=lambda x: x.quality)

        analytes.append(
            Analyte(
                analyte_id=len(analytes),
                neutral_mass=float((masses * weights).sum() / weights.sum()),
                rt_start=float(min(member.rt_start for member in members)),
                rt_apex=float(apex_member.rt_apex),
                rt_end=float(max(member.rt_end for member in members)),
                ms1_start=int(min(member.ms1_start for member in members)),
                ms1_apex=int(apex_member.ms1_apex),
                ms1_end=int(max(member.ms1_end for member in members)),
                charge_min=int(min(charges)),
                charge_max=int(max(charges)),
                n_distributions=len(members),
                score=float(np.mean([member.score for member in members])),
                members=[
                    {
                        "distribution_id": member.distribution_id,
                        "charge": member.charge,
                    }
                    for member in members
                ],
            )
        )

    return analytes


def edge_rows(edges):
    for edge_id, edge in enumerate(edges):
        yield {
            "edge_id": edge_id,
            "left_feature_id": edge["left_feature_id"],
            "right_feature_id": edge["right_feature_id"],
            "charge": edge["charge"],
            "isotope_step": edge["isotope_step"],
            "mz_error": edge["mz_error"],
            "mz_error_ppm": edge["mz_error_ppm"],
            "rt_shift": edge["rt_shift"],
            "rt_overlap": edge["rt_overlap"],
            "intensity_score": edge["intensity_score"],
            "score": edge["score"],
        }


def distribution_rows(distributions):
    for distribution in distributions:
        yield {
            "distribution_id": distribution.distribution_id,
            "charge": distribution.charge,
            "neutral_mass": distribution.neutral_mass,
            "mono_mz": distribution.mono_mz,
            "rt_start": distribution.rt_start,
            "rt_apex": distribution.rt_apex,
            "rt_end": distribution.rt_end,
            "ms1_start": distribution.ms1_start,
            "ms1_apex": distribution.ms1_apex,
            "ms1_end": distribution.ms1_end,
            "n_members": distribution.n_members,
            "score": distribution.score,
            "quality": distribution.quality,
        }


def distribution_member_rows(distributions):
    for distribution in distributions:
        for member in distribution.members:
            yield {
                "distribution_id": distribution.distribution_id,
                "feature_id": member["feature_id"],
                "isotope_index": member["isotope_index"],
                "member_score": member["member_score"],
            }


def analyte_rows(analytes):
    for analyte in analytes:
        yield {
            "analyte_id": analyte.analyte_id,
            "neutral_mass": analyte.neutral_mass,
            "rt_start": analyte.rt_start,
            "rt_apex": analyte.rt_apex,
            "rt_end": analyte.rt_end,
            "ms1_start": analyte.ms1_start,
            "ms1_apex": analyte.ms1_apex,
            "ms1_end": analyte.ms1_end,
            "charge_min": analyte.charge_min,
            "charge_max": analyte.charge_max,
            "n_distributions": analyte.n_distributions,
            "score": analyte.score,
        }


def analyte_member_rows(analytes):
    for analyte in analytes:
        for member in analyte.members:
            yield {
                "analyte_id": analyte.analyte_id,
                "distribution_id": member["distribution_id"],
                "charge": member["charge"],
            }


def make_config(args):
    return Config(
        line_mz_ppm=args.line_mz_ppm,
        line_mz_abs=args.line_mz_abs,
        max_gap_scans=args.max_gap_scans,
        min_trace_points=args.min_trace_points,
        peak_mindist=args.peak_mindist,
        smooth_points=args.smooth_points,
        min_peak_points=args.min_peak_points,
        min_peak_height=args.min_peak_height,
        min_peak_area=args.min_peak_area,
        min_peak_width=args.min_peak_width,
        max_peak_width=args.max_peak_width,
        min_peak_prominence_fraction=args.min_peak_prominence_fraction,
        max_trace_peaks=args.max_trace_peaks,
        max_charge=args.max_charge,
        isotope_mz_ppm=args.isotope_mz_ppm,
        isotope_mz_abs=args.isotope_mz_abs,
        max_neutral_mass=args.max_neutral_mass,
        max_apex_shift=args.max_apex_shift,
        max_apex_shift_width_fraction=args.max_apex_shift_width_fraction,
        min_edge_score=args.min_edge_score,
        min_distribution_members=args.min_distribution_members,
        charge_mass_ppm=args.charge_mass_ppm,
        min_charge_group_rt_score=args.min_charge_group_rt_score,
        charge_tolerance=args.charge_tolerance,
        mass_width_limit=args.mass_width_limit,
        step_limit=args.step_limit,
        new_inc_limit=args.new_inc_limit,
    )


def resolve_workers(value):
    if value > 0:
        return value

    cpu_count = os.cpu_count() or 2
    return max(1, min(8, cpu_count - 1))


def run(args):
    started = time()
    stage_started = time()

    config = make_config(args)
    workers = resolve_workers(args.workers)

    model = LineModel(config)
    scans = []

    for scan in stream_ms1(args.mzml, args.min_intensity):
        scans.append(
            {
                "ms1_index": scan["ms1_index"],
                "spectrum_index": scan["spectrum_index"],
                "scan_id": scan["scan_id"],
                "rt": scan["rt"],
                "tic": scan["tic"],
                "n_points": scan["n_points"],
            }
        )

        model.process_scan(
            ms1_index=scan["ms1_index"],
            rt=scan["rt"],
            mzs=scan["mzs"],
            intensities=scan["intensities"],
        )

        if args.progress and scan["ms1_index"] > 0 and scan["ms1_index"] % args.progress == 0:
            print(
                f"scans={scan['ms1_index']} active_lines={len(model.active)} "
                f"lines={len(model.lines)} features={len(model.features)}",
                file=sys.stderr,
                flush=True,
            )

    model.close_all()

    line_seconds = time() - stage_started
    stage_started = time()

    print(
        f"line_stage scans={len(scans)} lines={len(model.lines)} features={len(model.features)} seconds={line_seconds:.3f}",
        file=sys.stderr,
        flush=True,
    )

    best_edges, stored_edges = build_isotope_edges(
        features=model.features,
        config=config,
        workers=workers,
        chunk_size=args.edge_chunk_size,
        store_edges=args.store_edges,
        progress=bool(args.progress),
    )

    edge_seconds = time() - stage_started
    stage_started = time()

    print(
        f"edge_stage best_edges={len(best_edges)} stored_edges={len(stored_edges)} workers={workers} seconds={edge_seconds:.3f}",
        file=sys.stderr,
        flush=True,
    )

    distributions = build_distributions(
        features=model.features,
        best_edges=best_edges,
        config=config,
        workers=workers,
        progress=bool(args.progress),
    )

    distribution_seconds = time() - stage_started
    stage_started = time()

    print(
        f"distribution_stage distributions={len(distributions)} seconds={distribution_seconds:.3f}",
        file=sys.stderr,
        flush=True,
    )

    analytes = build_analytes(distributions, config)

    charge_seconds = time() - stage_started
    stage_started = time()

    print(
        f"charge_stage analytes={len(analytes)} seconds={charge_seconds:.3f}",
        file=sys.stderr,
        flush=True,
    )

    conn = connect_db(args.out, overwrite=args.overwrite)
    store_edges = args.store_edges != "none"

    try:
        init_schema(conn, store_edges=store_edges)

        write_parameters(
            conn,
            {
                "script": "distributions/index_ms1.py",
                "script_version": SCRIPT_VERSION,
                "mzml": str(Path(args.mzml).resolve()),
                "store_edges": args.store_edges,
                "workers": workers,
                "config": asdict(config),
                "counts": {
                    "scans": len(scans),
                    "lines": len(model.lines),
                    "features": len(model.features),
                    "best_edges": len(best_edges),
                    "stored_edges": len(stored_edges),
                    "distributions": len(distributions),
                    "analytes": len(analytes),
                },
            },
        )

        write_rows(conn, "scans", scans)
        write_rows(conn, "lines", model.lines)
        write_rows(conn, "features", feature_rows(model.features))

        if store_edges:
            write_rows(conn, "isotope_edges", edge_rows(stored_edges))

        write_rows(conn, "distributions", distribution_rows(distributions))
        write_rows(conn, "distribution_members", distribution_member_rows(distributions))
        write_rows(conn, "analytes", analyte_rows(analytes))
        write_rows(conn, "analyte_members", analyte_member_rows(analytes))

        create_indexes(conn, store_edges=store_edges)
        conn.commit()
    finally:
        conn.close()

    write_seconds = time() - stage_started
    total_seconds = time() - started

    print(
        json.dumps(
            {
                "out": str(args.out),
                "seconds": round(total_seconds, 3),
                "stage_seconds": {
                    "line": round(line_seconds, 3),
                    "edges": round(edge_seconds, 3),
                    "distributions": round(distribution_seconds, 3),
                    "charge": round(charge_seconds, 3),
                    "write": round(write_seconds, 3),
                },
                "scans": len(scans),
                "lines": len(model.lines),
                "features": len(model.features),
                "best_edges": len(best_edges),
                "stored_edges": len(stored_edges),
                "distributions": len(distributions),
                "analytes": len(analytes),
                "workers": workers,
                "store_edges": args.store_edges,
            },
            indent=2,
        )
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a broad library-free MS1 distribution deconvolution SQLite index from centroid mzML."
    )

    parser.add_argument("mzml", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--progress", type=int, default=500)

    parser.add_argument("--workers", type=int, default=0)
    parser.add_argument("--edge-chunk-size", type=int, default=20000)
    parser.add_argument("--store-edges", choices=("none", "best", "all"), default="none")

    parser.add_argument("--min-intensity", type=float, default=0.0)

    parser.add_argument("--line-mz-ppm", type=float, default=8.0)
    parser.add_argument("--line-mz-abs", type=float, default=0.002)
    parser.add_argument("--max-gap-scans", type=int, default=2)

    parser.add_argument("--min-trace-points", type=int, default=4)
    parser.add_argument("--peak-mindist", type=int, default=2)
    parser.add_argument("--smooth-points", type=int, default=3)

    parser.add_argument("--min-peak-points", type=int, default=4)
    parser.add_argument("--min-peak-height", type=float, default=0.0)
    parser.add_argument("--min-peak-area", type=float, default=0.0)
    parser.add_argument("--min-peak-width", type=float, default=0.0)
    parser.add_argument("--max-peak-width", type=float, default=6.0)
    parser.add_argument("--min-peak-prominence-fraction", type=float, default=0.02)
    parser.add_argument("--max-trace-peaks", type=int, default=0)

    parser.add_argument("--max-charge", type=int, default=6)
    parser.add_argument("--isotope-mz-ppm", type=float, default=10.0)
    parser.add_argument("--isotope-mz-abs", type=float, default=0.004)
    parser.add_argument("--max-neutral-mass", type=float, default=8000.0)

    parser.add_argument("--max-apex-shift", type=float, default=0.15)
    parser.add_argument("--max-apex-shift-width-fraction", type=float, default=0.50)
    parser.add_argument("--min-edge-score", type=float, default=0.30)
    parser.add_argument("--min-distribution-members", type=int, default=2)

    parser.add_argument("--charge-mass-ppm", type=float, default=12.0)
    parser.add_argument("--min-charge-group-rt-score", type=float, default=0.10)

    # reference isotope-edge acceptance (acdiff + intensity-step gating)
    parser.add_argument("--charge-tolerance", type=float, default=0.1)
    parser.add_argument("--mass-width-limit", type=float, default=0.002)
    parser.add_argument("--step-limit", type=float, default=0.5)
    parser.add_argument("--new-inc-limit", type=float, default=0.1)

    return parser.parse_args()


def main():
    args = parse_args()

    if not args.mzml.exists():
        raise SystemExit(f"missing mzML: {args.mzml}")

    if args.out.exists() and not args.overwrite:
        raise SystemExit(f"output exists; use --overwrite: {args.out}")

    run(args)


if __name__ == "__main__":
    main()
