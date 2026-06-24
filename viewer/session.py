import csv
import json
import re
import sqlite3
from pathlib import Path


PROTON = 1.007276554940804
C13_DELTA = 1.00335483507


def read_tsv(path):
    path = Path(path)

    if not path.exists():
        return []

    with path.open(newline="", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_json(path):
    path = Path(path)

    if not path.exists():
        return {}

    return json.loads(path.read_text(errors="replace"))


def safe_float(value, default=None):
    if value is None or value == "":
        return default

    try:
        return float(value)
    except Exception:
        return default


def safe_int(value, default=None):
    if value is None or value == "":
        return default

    try:
        return int(float(value))
    except Exception:
        return default


def safe_dir_name(filename):
    name = filename

    for suffix in [".centroid.mzML", ".centroid.mzml", ".mzML", ".mzml"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def strip_ms_suffix(path):
    name = Path(path).name
    lowered = name.lower()

    for suffix in (
        ".centroid.mzml.gz",
        ".centroid.mzml",
        ".mzml.gz",
        ".mzml",
        ".raw",
        ".wiff.scan",
        ".wiff2",
        ".wiff",
        ".d",
        ".lcd",
        ".baf",
        ".tdf_bin",
        ".tdf",
    ):
        if lowered.endswith(suffix):
            return name[: -len(suffix)]

    return Path(name).stem


def peptide_mass(row):
    calc_mass = safe_float(row.get("calc_mass"))
    exp_mass = safe_float(row.get("exp_mass"))

    if calc_mass is not None and calc_mass > 0:
        return calc_mass

    if exp_mass is not None and exp_mass > 0:
        return exp_mass

    return None


def peptide_charge(row):
    charge = safe_int(row.get("charge"))

    if charge is None or charge <= 0:
        return None

    return charge


def peptide_rt(row):
    rt = safe_float(row.get("rt"))

    if rt is not None:
        return rt

    rt = safe_float(row.get("aligned_rt"))

    if rt is not None:
        return rt

    return None


def isotope_mzs(neutral_mass, charge, n=6):
    mono_mz = neutral_mass / charge + PROTON
    step = C13_DELTA / charge

    return [mono_mz + i * step for i in range(n)]


class ViewerSession:
    def __init__(
        self,
        reorganized,
        distribution_db=None,
        centroid_dir=None,
        profile_dir=None,
    ):
        self.reorganized = Path(reorganized).resolve()
        self.manifest = read_json(self.reorganized / "manifest.json")

        manifest_mzml_dir = self.manifest.get("mzml_dir")
        manifest_by_file_dir = self.manifest.get("by_file_dir")

        self.mzml_dir = Path(centroid_dir).resolve() if centroid_dir else None

        if self.mzml_dir is None and manifest_mzml_dir:
            self.mzml_dir = Path(manifest_mzml_dir).resolve()

        self.by_file_dir = (
            Path(manifest_by_file_dir).resolve()
            if manifest_by_file_dir
            else self.reorganized / "by_file"
        )

        self.profile_dir = Path(profile_dir).resolve() if profile_dir else None
        self.distribution_db = Path(distribution_db).resolve() if distribution_db else None

        self.file_rows = read_tsv(self.reorganized / "files.tsv")
        self.file_by_name = {row.get("filename", ""): row for row in self.file_rows}

        self._psm_cache = {}

    def files(self):
        return self.file_rows

    def load_psms(self, filename):
        if filename in self._psm_cache:
            return self._psm_cache[filename]

        file_row = self.file_by_name.get(filename, {})
        run_dir = file_row.get("run_dir") or safe_dir_name(filename)

        candidates = [
            self.by_file_dir / run_dir / "scan_lookup.tsv",
            self.by_file_dir / run_dir / "psms.tsv",
        ]

        rows = []

        for path in candidates:
            rows = read_tsv(path)

            if rows:
                break

        self._psm_cache[filename] = rows
        return rows

    def centroid_path(self, filename):
        file_row = self.file_by_name.get(filename, {})
        mzml_path = file_row.get("mzml_path")

        if mzml_path:
            path = Path(mzml_path)

            if path.exists():
                return path.resolve()

        if self.mzml_dir is not None:
            path = self.mzml_dir / filename

            if path.exists():
                return path.resolve()

        return None

    def profile_path(self, filename):
        centroid = self.centroid_path(filename)

        if centroid is None:
            return None

        stem = strip_ms_suffix(centroid)

        candidates = []

        if self.profile_dir is not None:
            candidates.extend(
                [
                    self.profile_dir / f"{stem}.mzML",
                    self.profile_dir / f"{stem}.mzml",
                ]
            )

        candidates.extend(
            [
                centroid.with_name(f"{stem}.mzML"),
                centroid.with_name(f"{stem}.mzml"),
            ]
        )

        for path in candidates:
            if path.exists() and path.resolve() != centroid.resolve():
                return path.resolve()

        return None

    def distribution_candidates(
        self,
        neutral_mass,
        charge,
        rt,
        ppm=20.0,
        rt_window=1.0,
        limit=50,
    ):
        if self.distribution_db is None or not self.distribution_db.exists():
            return []

        if neutral_mass is None or charge is None or rt is None:
            return []

        tolerance = max(0.01, neutral_mass * ppm / 1_000_000.0)

        conn = sqlite3.connect(self.distribution_db)
        conn.row_factory = sqlite3.Row

        try:
            rows = conn.execute(
                """
                SELECT
                    distribution_id,
                    charge,
                    neutral_mass,
                    mono_mz,
                    rt_start,
                    rt_apex,
                    rt_end,
                    ms1_start,
                    ms1_apex,
                    ms1_end,
                    n_members,
                    score,
                    quality,
                    ABS(neutral_mass - ?) AS mass_error,
                    ABS(rt_apex - ?) AS rt_error
                FROM distributions
                WHERE charge = ?
                  AND neutral_mass BETWEEN ? AND ?
                  AND rt_apex BETWEEN ? AND ?
                ORDER BY rt_error ASC, mass_error ASC, quality DESC
                LIMIT ?
                """,
                (
                    neutral_mass,
                    rt,
                    charge,
                    neutral_mass - tolerance,
                    neutral_mass + tolerance,
                    rt - rt_window,
                    rt + rt_window,
                    limit,
                ),
            ).fetchall()

            return [dict(row) for row in rows]
        finally:
            conn.close()
