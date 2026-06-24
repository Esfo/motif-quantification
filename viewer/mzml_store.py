import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pyteomics import mzml


def scan_id(scan):
    return scan.get("id") or f"index={scan.get('index', 0)}"


def scan_number_from_id(value):
    if value is None:
        return ""

    value = str(value)
    match = re.search(r"scan=(\d+)", value)

    if match:
        return match.group(1)

    return value


def scan_ms_level(scan):
    return int(scan.get("ms level", 0) or 0)


def scan_rt(scan):
    try:
        return float(np.real(scan["scanList"]["scan"][0]["scan start time"]))
    except Exception:
        return None


def precursor_mz(scan):
    try:
        precursor = scan["precursorList"]["precursor"][0]
        selected = precursor["selectedIonList"]["selectedIon"][0]
        value = selected.get("selected ion m/z")

        if value is not None:
            return float(value)
    except Exception:
        pass

    return None


def scan_arrays(scan):
    mz_value = scan.get("m/z array")
    intensity_value = scan.get("intensity array")

    if mz_value is None:
        mza = np.array([], dtype=np.float64)
    else:
        mza = np.asarray(mz_value, dtype=np.float64).reshape(-1)

    if intensity_value is None:
        inten = np.array([], dtype=np.float64)
    else:
        inten = np.asarray(intensity_value, dtype=np.float64).reshape(-1)

    if mza.size != inten.size:
        sid = scan_id(scan)
        raise ValueError(f"{sid}: m/z and intensity arrays differ: {mza.size} != {inten.size}")

    return mza, inten


@dataclass(slots=True)
class ScanSummary:
    number: str
    spectrum_id: str
    level: int
    rt: float | None
    precursor_mz: float | None


class MzmlStore:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.loaded = False

        self.summaries = []
        self.by_number = {}
        self.by_id = {}
        self.ms1 = []

    def load_metadata(self):
        if self.loaded:
            return

        if not self.path.exists():
            raise FileNotFoundError(self.path)

        summaries = []

        with mzml.MzML(str(self.path), dtype=np.float64) as reader:
            for scan in reader:
                sid = scan_id(scan)
                number = scan_number_from_id(sid)
                level = scan_ms_level(scan)
                rt = scan_rt(scan)
                pmz = precursor_mz(scan)

                summary = ScanSummary(
                    number=number,
                    spectrum_id=sid,
                    level=level,
                    rt=rt,
                    precursor_mz=pmz,
                )

                summaries.append(summary)

        self.summaries = summaries
        self.by_number = {summary.number: summary for summary in summaries if summary.number}
        self.by_id = {summary.spectrum_id: summary for summary in summaries if summary.spectrum_id}
        self.ms1 = [summary for summary in summaries if summary.level == 1 and summary.rt is not None]
        self.loaded = True

    def get_summary_by_number(self, number):
        self.load_metadata()
        return self.by_number.get(str(number))

    def get_scan_by_number(self, number):
        self.load_metadata()
        summary = self.get_summary_by_number(number)

        if summary is None:
            return None

        return self.get_scan_by_id(summary.spectrum_id)

    def get_scan_by_id(self, spectrum_id):
        spectrum_id = str(spectrum_id)

        with mzml.MzML(str(self.path), dtype=np.float64) as reader:
            if hasattr(reader, "get_by_id"):
                try:
                    return reader.get_by_id(spectrum_id)
                except Exception:
                    pass

            for scan in reader:
                if scan_id(scan) == spectrum_id:
                    return scan

        return None

    def nearest_ms1_by_rt(self, rt):
        self.load_metadata()

        if rt is None or not self.ms1:
            return None

        return min(self.ms1, key=lambda x: abs((x.rt or 0.0) - rt))

    def preceding_ms1_for_scan(self, number):
        self.load_metadata()

        try:
            selected = int(number)
        except Exception:
            return None

        best = None
        best_number = None

        for summary in self.ms1:
            try:
                current = int(summary.number)
            except Exception:
                continue

            if current <= selected and (best_number is None or current > best_number):
                best = summary
                best_number = current

        return best

    def ms1_summaries_in_rt(self, rt_start, rt_end):
        self.load_metadata()

        return [
            summary
            for summary in self.ms1
            if summary.rt is not None and rt_start <= summary.rt <= rt_end
        ]

    def extract_xics(self, targets, rt_start, rt_end, ppm=10.0, abs_tol=0.01):
        targets = [float(x) for x in targets]
        rts = []
        traces = [[] for _ in targets]

        with mzml.MzML(str(self.path), dtype=np.float64) as reader:
            for scan in reader:
                if scan_ms_level(scan) != 1:
                    continue

                rt = scan_rt(scan)

                if rt is None or rt < rt_start or rt > rt_end:
                    continue

                mza, inten = scan_arrays(scan)

                if mza.size == 0:
                    continue

                order = np.argsort(mza)
                mza = mza[order]
                inten = inten[order]

                rts.append(rt)

                for target_i, target in enumerate(targets):
                    tolerance = max(abs_tol, target * ppm / 1_000_000.0)
                    left = np.searchsorted(mza, target - tolerance, side="left")
                    right = np.searchsorted(mza, target + tolerance, side="right")

                    if right > left:
                        value = float(np.sum(inten[left:right]))
                    else:
                        value = 0.0

                    traces[target_i].append(value)

        return {
            "rts": np.asarray(rts, dtype=np.float64),
            "targets": np.asarray(targets, dtype=np.float64),
            "traces": [np.asarray(trace, dtype=np.float64) for trace in traces],
        }

    def scan_window_by_number(self, number, mz_min, mz_max):
        scan = self.get_scan_by_number(number)

        if scan is None:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)

        return self.scan_window(scan, mz_min, mz_max)

    def scan_window(self, scan, mz_min, mz_max):
        mza, inten = scan_arrays(scan)

        if mza.size == 0:
            return mza, inten

        keep = (mza >= mz_min) & (mza <= mz_max)
        return mza[keep], inten[keep]
