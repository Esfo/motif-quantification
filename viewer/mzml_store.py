import queue
import re
from concurrent.futures import ThreadPoolExecutor
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
        self.ms2 = []

        self._data_reader = None

    def data_reader(self):
        """Persistent indexed reader for random-access scan reads.

        Reused across calls so the (potentially multi-GB) file is indexed once
        rather than reopened per scan.
        """
        if self._data_reader is None:
            self._data_reader = mzml.MzML(str(self.path), dtype=np.float64, use_index=True)

        return self._data_reader

    def load_metadata(self):
        if self.loaded:
            return

        if not self.path.exists():
            raise FileNotFoundError(self.path)

        summaries = []

        # decode_binary=False skips base64/zlib decoding of the m/z and
        # intensity arrays: metadata only needs id, ms level, and RT, so we
        # avoid decoding gigabytes of peak data just to enumerate scans.
        with mzml.MzML(str(self.path), use_index=True, decode_binary=False) as reader:
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
        self.ms2 = [summary for summary in summaries if summary.level == 2 and summary.rt is not None]
        self.loaded = True

    def ms2_in_rt(self, rt_start, rt_end):
        """MS2 scan summaries (rt, precursor m/z) inside an RT window."""
        self.load_metadata()
        return [s for s in getattr(self, "ms2", [])
                if s.rt is not None and rt_start <= s.rt <= rt_end]

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
        reader = self.data_reader()

        try:
            return reader.get_by_id(spectrum_id)
        except Exception:
            return None

    def ms1_in_rt(self, rt_start, rt_end):
        self.load_metadata()

        return [
            summary
            for summary in self.ms1
            if summary.rt is not None and rt_start <= summary.rt <= rt_end
        ]

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

        reader = self.data_reader()

        for summary in self.ms1_in_rt(rt_start, rt_end):
            scan = reader.get_by_id(summary.spectrum_id)

            if scan is None:
                continue

            mza, inten = scan_arrays(scan)

            if mza.size == 0:
                continue

            order = np.argsort(mza)
            mza = mza[order]
            inten = inten[order]

            rts.append(summary.rt)

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

    # ---- parallel scan reading -------------------------------------------

    def _read_scans(self, summaries, workers=4):
        """Read (rt, mz_array, intensity_array) for each summary CONCURRENTLY.

        The per-scan cost is dominated by base64+zlib decoding, which releases
        the GIL, so a thread pool with independent readers (one open file handle
        each) genuinely parallelises it. Returns a list aligned to ``summaries``;
        entries are (rt, None, None) for scans that couldn't be read. Falls back
        to a single serial reader if the pool can't be built.
        """
        summaries = list(summaries)
        if not summaries:
            return []
        n = min(workers, len(summaries))

        # Lazy pool of independent readers (each its own file handle), reused
        # across calls so the index is only built once per reader.
        if not hasattr(self, "_reader_pool"):
            self._reader_pool = queue.Queue()
            self._reader_pool_size = 0
            self._read_executor = None
        while self._reader_pool_size < n:
            try:
                self._reader_pool.put(mzml.MzML(str(self.path), dtype=np.float64, use_index=True))
                self._reader_pool_size += 1
            except Exception:
                break

        if self._reader_pool_size == 0:
            # Fall back to the single shared reader, serially.
            reader = self.data_reader()
            out = []
            for s in summaries:
                try:
                    scan = reader.get_by_id(s.spectrum_id)
                    if scan is None:
                        out.append((s.rt, None, None))
                    else:
                        mza, inten = scan_arrays(scan)
                        out.append((s.rt, mza, inten))
                except Exception:
                    out.append((s.rt, None, None))
            return out

        def work(summary):
            reader = self._reader_pool.get()
            try:
                scan = reader.get_by_id(summary.spectrum_id)
                if scan is None:
                    return (summary.rt, None, None)
                mza, inten = scan_arrays(scan)
                return (summary.rt, mza, inten)
            except Exception:
                return (summary.rt, None, None)
            finally:
                self._reader_pool.put(reader)

        if self._read_executor is None:
            self._read_executor = ThreadPoolExecutor(max_workers=max(1, self._reader_pool_size))
        return list(self._read_executor.map(work, summaries))

    def extract_region(
        self,
        mz_min,
        mz_max,
        rt_start,
        rt_end,
        mz_bins=600,
        mode="profile",
    ):
        """Rasterize MS1 data in an RT x m/z window onto a regular grid.

        Returns rts (one per MS1 scan in range), a shared mz grid, and a
        Z matrix of shape (n_scans, mz_bins). Profile data is linearly
        interpolated onto the grid (it is continuous); centroid data is
        rasterized by assigning each peak to its nearest bin and keeping the
        max (interpolation would blur or zero out sparse spikes).

        Nothing here models or fits anything; it only bins measured points.
        """
        mz_min = float(mz_min)
        mz_max = float(mz_max)
        mz_bins = max(2, int(mz_bins))

        mz_grid = np.linspace(mz_min, mz_max, mz_bins)
        edges = np.linspace(mz_min, mz_max, mz_bins + 1)

        rts = []
        rows = []

        for rt, mza, inten in self._read_scans(self.ms1_in_rt(rt_start, rt_end)):
            if mza is None:
                continue

            if mza.size:
                keep = (mza >= mz_min) & (mza <= mz_max)
                mza = mza[keep]
                inten = inten[keep]

            row = np.zeros(mz_bins, dtype=np.float64)

            if mza.size:
                order = np.argsort(mza)
                mza = mza[order]
                inten = inten[order]

                if mode == "profile":
                    row = np.interp(mz_grid, mza, inten, left=0.0, right=0.0)
                else:
                    idx = np.clip(np.searchsorted(edges, mza, side="right") - 1, 0, mz_bins - 1)
                    np.maximum.at(row, idx, inten)

            rts.append(rt)
            rows.append(row)

        if rows:
            z = np.vstack(rows)
        else:
            z = np.zeros((0, mz_bins), dtype=np.float64)

        return {
            "rts": np.asarray(rts, dtype=np.float64),
            "mz_grid": mz_grid,
            "z": z,
        }

    def extract_points(self, mz_min, mz_max, rt_start, rt_end):
        """Every measured point in an RT x m/z window, across all MS1 scans.

        Returns parallel arrays (mz, rt, intensity) of the raw datapoints -- not
        binned or averaged -- so panel 1 can show every point of a line in 2D
        (m/z vs intensity) and 3D (m/z, rt, intensity).
        """
        mzs, rts, intens = [], [], []

        for rt, mza, inten in self._read_scans(self.ms1_in_rt(rt_start, rt_end)):
            if mza is None or mza.size == 0:
                continue
            keep = (mza >= mz_min) & (mza <= mz_max)
            if not keep.any():
                continue
            kept_mz = mza[keep]
            mzs.append(kept_mz)
            intens.append(inten[keep])
            rts.append(np.full(kept_mz.size, rt, dtype=np.float64))

        if mzs:
            return {
                "mz": np.concatenate(mzs),
                "rt": np.concatenate(rts),
                "intensity": np.concatenate(intens),
            }
        return {
            "mz": np.array([], dtype=np.float64),
            "rt": np.array([], dtype=np.float64),
            "intensity": np.array([], dtype=np.float64),
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
