"""Reader for the project ``experimental-setup`` design file.

Format (csv, header row):
    filename,condition,fraction,replicate,pair_id

The ``filename`` column matches the mzML stem (without ``.centroid``/``.mzML``).
Tabs 3 (file-by-file comparison) and 4 (motif quantification) use this to group
runs into conditions/fractions/pairs for time-series and differential
expression. Grouping here is deliberately generic: callers pick any column(s)
to define groups and contrasts, so the design can be sliced arbitrarily.
"""

import csv
from collections import defaultdict
from pathlib import Path

COLUMNS = ["filename", "condition", "fraction", "replicate", "pair_id"]


class ExperimentalSetup:
    def __init__(self, rows):
        self.rows = rows
        self.by_filename = {r.get("filename", ""): r for r in rows}

    @classmethod
    def load(cls, path):
        path = Path(path)
        if not path.exists():
            return cls([])

        with path.open(newline="", errors="replace") as f:
            # tolerate both comma and tab just in case
            sample = f.read(2048)
            f.seek(0)
            delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
            rows = list(csv.DictReader(f, delimiter=delimiter))

        # normalize whitespace
        cleaned = []
        for row in rows:
            cleaned.append({(k or "").strip(): (v or "").strip() for k, v in row.items()})
        return cls(cleaned)

    def is_empty(self):
        return not self.rows

    def columns(self):
        """Column names actually present, in file order, ``filename`` first.

        The design file is not restricted to the canonical five columns — any
        extra column (a batch, a timepoint, a dose…) is exposed here so the
        Quantitative Comparisons tab can assign it a role."""
        seen = []
        for row in self.rows:
            for key in row.keys():
                if key and key not in seen:
                    seen.append(key)
        # keep filename first if present, otherwise preserve discovery order
        if "filename" in seen:
            seen = ["filename"] + [c for c in seen if c != "filename"]
        return seen

    def values(self, column):
        seen = []
        for row in self.rows:
            v = row.get(column, "")
            if v and v not in seen:
                seen.append(v)
        return seen

    def group_by(self, *columns):
        """Map each combination of column values to the filenames in that group."""
        groups = defaultdict(list)
        for row in self.rows:
            key = tuple(row.get(c, "") for c in columns)
            groups[key].append(row.get("filename", ""))
        return dict(groups)

    def filenames_for(self, **filters):
        """Filenames whose row matches every column=value filter."""
        out = []
        for row in self.rows:
            if all(row.get(k, "") == v for k, v in filters.items()):
                out.append(row.get("filename", ""))
        return out

    def condition_of(self, filename):
        return self.by_filename.get(filename, {}).get("condition", "")
