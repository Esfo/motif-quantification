"""Feature-quantity matrix for the Quantitative Comparisons tab.

Turns the reorganized per-file quant tables into ``{feature: {filename:
quantity}}`` matrices, at the **peptide** or **protein** level. Nothing here
knows anything about the experimental design — it only maps features to their
per-file quantities; the tab layers the (fully user-defined) grouping on top.

**Quantity source.** When the off-GUI ``quantify.py`` stage has written a
``quant/`` folder next to ``searches/``, the AUC-based tables there are used:
``protein_quant.tsv`` (a protein = the charge-distribution AUC of its single most
abundant unique peptide, the same number the MS Data tab shows) and
``peptide_auc.tsv`` (per-peptide AUC). Otherwise it falls back to the Sage-LFQ
``quantity`` column of the reorganized ``peptide_quant.tsv``.

Peptide level keys on ``peptide_plain`` (charge/mod variants summed). Each
peptide also carries a ``unique`` flag (maps to exactly one protein).

Protein level is **unique quantification**: a protein's per-file quantity is the
sum of its *unique* peptides only (peptides that map to that protein alone), so
shared peptides never double-count across proteins.

Matrices are cached so switching Peptides⇄Proteins is instant.
"""

import csv
import re


def _split_proteins(value):
    if not value:
        return []
    return [p for p in re.split(r"[;,\s]+", value.strip()) if p]


def _to_float(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class QuantModel:
    def __init__(self, session):
        self.session = session
        self._peptide = None
        self._protein = None
        self._pep_unique = {}
        self._pep_proteins = {}

    def filenames(self):
        return [r.get("filename", "") for r in (self.session.files() or [])]

    # ---- AUC tables (quantify.py) ----------------------------------------

    def _quant_dir(self):
        return getattr(self.session, "quant_dir", None)

    def _read_quant_tsv(self, name):
        d = self._quant_dir()
        if d is None:
            return None
        path = d / name
        if not path.exists():
            return None
        with path.open(newline="", errors="replace") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    # ---- peptide-level ---------------------------------------------------

    def _ensure_peptide(self):
        if self._peptide is not None:
            return
        auc_rows = self._read_quant_tsv("peptide_auc.tsv")
        if auc_rows is not None:
            self._build_peptide_from_auc(auc_rows)
            return

        matrix = {}
        pep_proteins = {}
        for row in self.session.all_quant_rows():
            pep = row.get("peptide_plain", "")
            if not pep:
                continue
            q = _to_float(row.get("quantity"))
            if q is None:
                continue
            fname = row.get("filename", "")
            per_file = matrix.setdefault(pep, {})
            per_file[fname] = per_file.get(fname, 0.0) + q
            if pep not in pep_proteins:
                pep_proteins[pep] = _split_proteins(row.get("proteins", ""))
        self._peptide = matrix
        self._pep_proteins = pep_proteins
        self._pep_unique = {p: (len(prots) == 1) for p, prots in pep_proteins.items()}

    def _build_peptide_from_auc(self, rows):
        """Peptide matrix from ``quant/peptide_auc.tsv`` (AUC per peptide/file)."""
        matrix = {}
        pep_proteins = {}
        for row in rows:
            pep = row.get("peptide_plain", "")
            if not pep:
                continue
            q = _to_float(row.get("auc"))
            if q is None:
                continue
            fname = row.get("filename", "")
            matrix.setdefault(pep, {})[fname] = q
            if pep not in pep_proteins:
                pep_proteins[pep] = _split_proteins(row.get("proteins", ""))
        self._peptide = matrix
        self._pep_proteins = pep_proteins
        self._pep_unique = {p: (len(prots) == 1) for p, prots in pep_proteins.items()}

    def peptide_matrix(self):
        self._ensure_peptide()
        return self._peptide

    def peptide_is_unique(self, peptide_plain):
        self._ensure_peptide()
        return self._pep_unique.get(peptide_plain, False)

    def proteins_for_peptide(self, peptide_plain):
        self._ensure_peptide()
        return self._pep_proteins.get(peptide_plain, [])

    # ---- protein-level (unique quant) -----------------------------------

    def _ensure_protein(self):
        if self._protein is not None:
            return
        prot_rows = self._read_quant_tsv("protein_quant.tsv")
        if prot_rows is not None:
            matrix = {}
            for row in prot_rows:
                prot = row.get("protein", "")
                q = _to_float(row.get("quantity"))
                if not prot or q is None:
                    continue
                matrix.setdefault(prot, {})[row.get("filename", "")] = q
            self._protein = matrix
            return

        self._ensure_peptide()
        matrix = {}
        for pep, per_file in self._peptide.items():
            prots = self._pep_proteins.get(pep, [])
            if len(prots) != 1:
                continue  # unique quantification: skip shared peptides
            prot = prots[0]
            dest = matrix.setdefault(prot, {})
            for fname, q in per_file.items():
                dest[fname] = dest.get(fname, 0.0) + q
        self._protein = matrix

    def protein_matrix(self):
        self._ensure_protein()
        return self._protein

    def matrix(self, level):
        return self.peptide_matrix() if level == "peptide" else self.protein_matrix()
