"""Option B reader: serve a window's raw centroids from the distributions sqlite
(scan_points table) instead of re-decoding the mzML on every zoom.

PointsStore subclasses MzmlStore and overrides only the data source, so the
tested extract_points / extract_region / scan_window binning logic is reused.
"""

import sqlite3
import zlib

import numpy as np

try:
    from .mzml_store import MzmlStore, ScanSummary
except ImportError:
    from mzml_store import MzmlStore, ScanSummary


def _decode(blob):
    return np.frombuffer(zlib.decompress(blob), dtype="<f4").astype(np.float64)


class PointsStore(MzmlStore):
    _is_points = True

    def __init__(self, sqlite_path):
        super().__init__(sqlite_path)  # sets self.path / empty caches
        self._conn = None

    @staticmethod
    def has_points(sqlite_path):
        try:
            conn = sqlite3.connect(f"file:{sqlite_path}?mode=ro", uri=True)
            try:
                row = conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' AND name='scan_points'"
                ).fetchone()
                if row is None:
                    return False
                return conn.execute("SELECT 1 FROM scan_points LIMIT 1").fetchone() is not None
            finally:
                conn.close()
        except Exception:
            return False

    def _connect(self):
        if self._conn is None:
            # check_same_thread=False: the extraction worker runs on a QThread.
            self._conn = sqlite3.connect(
                f"file:{self.path}?mode=ro", uri=True, check_same_thread=False
            )
        return self._conn

    def load_metadata(self):
        if self.loaded:
            return
        rows = self._connect().execute(
            "SELECT ms1_index, rt FROM scans ORDER BY ms1_index"
        ).fetchall()
        # number == spectrum_id == str(ms1_index): the scan_points key.
        self.ms1 = [
            ScanSummary(number=str(i), spectrum_id=str(i), level=1, rt=float(rt),
                        precursor_mz=None)
            for (i, rt) in rows
        ]
        self.ms2 = []
        self.summaries = list(self.ms1)
        self.by_number = {s.number: s for s in self.ms1}
        self.by_id = {s.spectrum_id: s for s in self.ms1}
        self.loaded = True

    def _read_scans(self, summaries, workers=4):
        conn = self._connect()
        out = []
        for s in summaries:
            row = conn.execute(
                "SELECT mz, intensity FROM scan_points WHERE ms1_index = ?",
                (int(s.number),),
            ).fetchone()
            if row is None:
                out.append((s.rt, None, None))
            else:
                out.append((s.rt, _decode(row[0]), _decode(row[1])))
        return out

    def scan_window_by_number(self, number, mz_min, mz_max):
        row = self._connect().execute(
            "SELECT mz, intensity FROM scan_points WHERE ms1_index = ?",
            (int(number),),
        ).fetchone()
        if row is None:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)
        mz = _decode(row[0])
        inten = _decode(row[1])
        keep = (mz >= mz_min) & (mz <= mz_max)
        return mz[keep], inten[keep]
