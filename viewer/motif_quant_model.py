"""Motif-group quantity matrix for the Motifs tab.

Reads the ``motif-sets/`` folder produced by the off-GUI ``quantify.py`` stage
(``motifs.tsv`` + ``motif_quant.tsv``) and exposes the same ``{feature:
{filename: quantity}}`` shape as ``quant_model.QuantModel`` so the Motifs tab can
reuse the Quantitative-Comparisons machinery wholesale.

A "feature" here is a skeleton motif (e.g. ``A.....G``). Its per-file quantity is
the SUM of the AUC-based quantities of every observed protein the motif groups
(computed upstream in ``quantify.py``). Only motifs with more than one observed
protein are written, so the tab only ever sees genuine multi-protein groups.

The feature key is the motif text (unique per motif; on the rare collision the
motif_id is appended) so it doubles as the display label everywhere the tab shows
a feature.
"""

import csv
from pathlib import Path


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class MotifQuantModel:
    def __init__(self, session):
        self.session = session
        self._matrix = None
        self._observed = {}       # feature -> observed protein count
        self._accessions = {}     # feature -> [accession, ...]
        self._motif_id = {}       # feature -> motif_id
        self.directory = session.motif_sets_dir if session else None

    # ---- availability ----------------------------------------------------

    def is_available(self):
        d = self.directory
        return bool(d and (d / "motif_quant.tsv").exists()
                    and (d / "motifs.tsv").exists())

    def filenames(self):
        return [r.get("filename", "") for r in (self.session.files() or [])]

    # ---- loading ---------------------------------------------------------

    def _read_tsv(self, name):
        path = self.directory / name
        if not path.exists():
            return []
        with path.open(newline="", errors="replace") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def _ensure(self):
        if self._matrix is not None:
            return
        matrix = {}
        # motif_id -> feature label, from motifs.tsv (the group metadata).
        label_by_id = {}
        for row in self._read_tsv("motifs.tsv"):
            try:
                motif_id = int(row.get("motif_id", ""))
            except (TypeError, ValueError):
                continue
            text = row.get("motif_text", "") or str(motif_id)
            label = text
            # Disambiguate the rare case of two motifs sharing a skeleton text.
            if label in self._observed and self._motif_id.get(label) != motif_id:
                label = f"{text}#{motif_id}"
            label_by_id[motif_id] = label
            self._motif_id[label] = motif_id
            self._observed[label] = int(row.get("observed_count", 0) or 0)
            self._accessions[label] = [
                a for a in (row.get("observed_accessions", "") or "").split(";") if a]

        for row in self._read_tsv("motif_quant.tsv"):
            try:
                motif_id = int(row.get("motif_id", ""))
            except (TypeError, ValueError):
                continue
            label = label_by_id.get(motif_id)
            if label is None:
                continue
            q = _to_float(row.get("quantity"))
            if q is None:
                continue
            fname = row.get("filename", "")
            matrix.setdefault(label, {})[fname] = q
        self._matrix = matrix

    # ---- QuantModel-compatible API --------------------------------------

    def matrix(self, level=None):
        """The motif-group quantity matrix. ``level`` is accepted for API parity
        with QuantModel but ignored (motifs have a single level)."""
        self._ensure()
        return self._matrix

    def observed_count(self, feature):
        self._ensure()
        return self._observed.get(feature, 0)

    def accessions(self, feature):
        self._ensure()
        return self._accessions.get(feature, [])

    def max_observed(self):
        self._ensure()
        return max(self._observed.values(), default=0)
