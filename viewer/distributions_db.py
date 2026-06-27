"""Reader for the MS1 distributions sqlite (the ``distributions/`` output).

Schema mirrors ``distributions/store.py``: ``scans``, ``lines``, ``features``,
``distributions``, ``distribution_members``, ``analytes``, ``analyte_members``.
A "line" the spec talks about is a ``feature`` here (one isotope trace's peak);
a "distribution" is the group of features linked through ``distribution_members``
with an ``isotope_index``.

This reader is read-only and stateless beyond a cached connection, so it is safe
to query from the UI thread for the small windowed lookups the viewer needs.
Raw per-point arrays are *not* in the DB; the viewer re-extracts those from the
mzML window (see ``mzml_store.extract_region``).
"""

import sqlite3
from pathlib import Path


class DistributionsDB:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self._conn = None

    def connect(self):
        if self._conn is None:
            self._conn = sqlite3.connect(f"file:{self.path}?mode=ro", uri=True)
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def close(self):
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def parameters(self):
        rows = self.connect().execute("SELECT key, value FROM parameters").fetchall()
        return {r["key"]: r["value"] for r in rows}

    def distributions_in_window(self, mz_min=None, mz_max=None, rt_start=None,
                                rt_end=None, charge=None, limit=2000):
        """Distributions whose mono m/z and RT apex fall in the window."""
        clauses = []
        params = []

        if mz_min is not None:
            clauses.append("mono_mz >= ?")
            params.append(mz_min)
        if mz_max is not None:
            clauses.append("mono_mz <= ?")
            params.append(mz_max)
        if rt_start is not None:
            clauses.append("rt_end >= ?")
            params.append(rt_start)
        if rt_end is not None:
            clauses.append("rt_start <= ?")
            params.append(rt_end)
        if charge is not None:
            clauses.append("charge = ?")
            params.append(int(charge))

        where = ("WHERE " + " AND ".join(clauses)) if clauses else ""
        params.append(limit)

        rows = self.connect().execute(
            f"SELECT * FROM distributions {where} ORDER BY quality DESC LIMIT ?", params
        ).fetchall()
        return [dict(r) for r in rows]

    def distribution(self, distribution_id):
        row = self.connect().execute(
            "SELECT * FROM distributions WHERE distribution_id = ?", (distribution_id,)
        ).fetchone()
        return dict(row) if row else None

    def distribution_members(self, distribution_id):
        """Member features (the isotope lines) of a distribution, ordered by isotope."""
        rows = self.connect().execute(
            """
            SELECT f.*, m.isotope_index, m.member_score
            FROM distribution_members m
            JOIN features f ON f.feature_id = m.feature_id
            WHERE m.distribution_id = ?
            ORDER BY m.isotope_index
            """,
            (distribution_id,),
        ).fetchall()
        return [dict(r) for r in rows]

    def feature(self, feature_id):
        row = self.connect().execute(
            "SELECT * FROM features WHERE feature_id = ?", (feature_id,)
        ).fetchone()
        return dict(row) if row else None

    def all_lines(self):
        """Every feature ("line") with the charge of the distribution it belongs
        to (0 when it is in none). Used by the table-1 'lines' tab."""
        rows = self.connect().execute(
            """
            SELECT f.feature_id, f.mz_mean, f.mz_min, f.mz_max, f.rt_apex,
                   f.rt_start, f.rt_end, f.n_points, f.area, f.height,
                   COALESCE(d.charge, 0) AS charge
            FROM features f
            LEFT JOIN distribution_members m ON m.feature_id = f.feature_id
            LEFT JOIN distributions d ON d.distribution_id = m.distribution_id
            ORDER BY f.feature_id
            """
        ).fetchall()
        return [dict(r) for r in rows]

    def all_distributions(self):
        # AUC = sum of member-feature areas = area under the combined isotope peak
        # (integral of the summed per-timepoint signal across all lines).
        rows = self.connect().execute(
            """
            SELECT d.*, IFNULL(a.auc, 0.0) AS auc
            FROM distributions d
            LEFT JOIN (
                SELECT m.distribution_id, SUM(f.area) AS auc
                FROM distribution_members m
                JOIN features f ON f.feature_id = m.feature_id
                GROUP BY m.distribution_id
            ) a ON a.distribution_id = d.distribution_id
            """
        ).fetchall()
        return [dict(r) for r in rows]

    def all_analytes(self):
        rows = self.connect().execute("SELECT * FROM analytes").fetchall()
        return [dict(r) for r in rows]

    def all_analytes_multicharge(self):
        """Analytes that span more than one charge state, each with one
        representative member distribution_id (for a one-click panel-3 load) and
        an AUC = total area under the combined peak of every line of every member
        distribution."""
        rows = self.connect().execute(
            """
            SELECT a.*,
                   (SELECT am.distribution_id FROM analyte_members am
                    WHERE am.analyte_id = a.analyte_id LIMIT 1) AS rep_distribution_id,
                   IFNULL(s.auc, 0.0) AS auc
            FROM analytes a
            LEFT JOIN (
                SELECT am.analyte_id, SUM(f.area) AS auc
                FROM analyte_members am
                JOIN distribution_members m ON m.distribution_id = am.distribution_id
                JOIN features f ON f.feature_id = m.feature_id
                GROUP BY am.analyte_id
            ) s ON s.analyte_id = a.analyte_id
            WHERE a.charge_max > a.charge_min
            ORDER BY a.analyte_id
            """
        ).fetchall()
        return [dict(r) for r in rows]

    def analyte_distributions(self, analyte_id):
        """The member distributions of an analyte (charge region)."""
        rows = self.connect().execute(
            """
            SELECT d.*, IFNULL(a.auc, 0.0) AS auc
            FROM analyte_members am
            JOIN distributions d ON d.distribution_id = am.distribution_id
            LEFT JOIN (
                SELECT m.distribution_id, SUM(f.area) AS auc
                FROM distribution_members m
                JOIN features f ON f.feature_id = m.feature_id
                GROUP BY m.distribution_id
            ) a ON a.distribution_id = d.distribution_id
            WHERE am.analyte_id = ?
            ORDER BY d.charge
            """,
            (analyte_id,),
        ).fetchall()
        return [dict(r) for r in rows]

    def line(self, line_id):
        row = self.connect().execute(
            "SELECT * FROM lines WHERE line_id = ?", (line_id,)
        ).fetchone()
        return dict(row) if row else None

    def scan_rt(self, ms1_index):
        row = self.connect().execute(
            "SELECT rt FROM scans WHERE ms1_index = ?", (ms1_index,)
        ).fetchone()
        return row["rt"] if row else None

    def analyte_for_distribution(self, distribution_id):
        row = self.connect().execute(
            "SELECT analyte_id FROM analyte_members WHERE distribution_id = ?",
            (distribution_id,),
        ).fetchone()
        return row["analyte_id"] if row else None

    def charge_group(self, distribution_id):
        """The analyte's distributions across charge states, with their member
        features (isotope lines). This is the data source for the panel-3
        charge-comparison grid: {charge: {distribution: {...}, features: [...]}}.

        Falls back to just the given distribution when it has no analyte.
        """
        analyte_id = self.analyte_for_distribution(distribution_id)
        if analyte_id is None:
            dist = self.distribution(distribution_id)
            if dist is None:
                return {}
            return {dist["charge"]: {"distribution": dist,
                                     "features": self.distribution_members(distribution_id)}}

        rows = self.connect().execute(
            "SELECT distribution_id, charge FROM analyte_members WHERE analyte_id = ? ORDER BY charge",
            (analyte_id,),
        ).fetchall()

        # An analyte can hold more than one distribution at the same charge. The
        # grid shows one per charge, so make sure the explicitly requested
        # distribution wins its own charge slot (otherwise panel 3 would render a
        # same-charge sibling -- a different colour/size than the one clicked).
        group = {}
        for row in rows:
            did = row["distribution_id"]
            charge = row["charge"]
            existing = group.get(charge)
            if existing is not None:
                if existing["distribution"]["distribution_id"] == distribution_id:
                    continue  # never displace the clicked distribution
                if did != distribution_id:
                    continue  # keep first-seen unless this row IS the clicked one
            group[charge] = {
                "distribution": self.distribution(did),
                "features": self.distribution_members(did),
            }
        return group
