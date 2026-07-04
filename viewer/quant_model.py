"""Feature-quantity matrix for the Quantitative Comparisons tab.

Turns the reorganized per-file quant tables into a ``{feature: {filename:
quantity}}`` matrix at either the **peptide** or **protein** level, so the tab
can compare any feature across any grouping of files.

Peptide level keys on ``peptide_plain`` (mod/charge variants collapse together,
summed). Protein level rolls peptide quantities up to each protein id the
peptide maps to (the quant row's ``proteins`` field), with three roll-up modes:

    sum     — sum every mapped peptide's quantity (default, standard label-free)
    median  — median of the mapped peptides' quantities (robust to one big peptide)
    unique  — sum only peptides that map to a single protein (razor/unique)

The matrix is cached per (level, rollup) so switching Peptides⇄Proteins and back
is instant.
"""

import re
import statistics


def _split_proteins(value):
    """Protein-id list from a quant row's ``proteins`` field.

    reorganize-results.py joins protein ids; tolerate ``;`` / ``,`` / whitespace
    separators and drop blanks."""
    if not value:
        return []
    return [p for p in re.split(r"[;,\s]+", value.strip()) if p]


def _to_float(value):
    if value is None or value == "":
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f


class QuantModel:
    """Builds and caches feature×file quantity matrices from a ViewerSession."""

    ROLLUPS = ("sum", "median", "unique")

    def __init__(self, session):
        self.session = session
        self._cache = {}

    def filenames(self):
        """Filenames that have any quant rows, in files.tsv order."""
        return [r.get("filename", "") for r in (self.session.files() or [])]

    def _peptide_matrix(self):
        """{peptide_plain: {filename: summed quantity}}."""
        matrix = {}
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
        return matrix

    def _protein_matrix(self, rollup):
        """{protein_id: {filename: rolled-up quantity}} for the given roll-up.

        Works from the raw quant rows so per-file peptide sets are respected. For
        ``median`` we collect each peptide's per-file summed quantity first, then
        take the median across peptides per (protein, file)."""
        # protein -> file -> list of per-peptide quantities in that file
        collected = {}
        # First collapse charge/mod variants to one quantity per (peptide, file).
        pep_file = {}
        pep_proteins = {}
        for row in self.session.all_quant_rows():
            pep = row.get("peptide_plain", "")
            if not pep:
                continue
            q = _to_float(row.get("quantity"))
            if q is None:
                continue
            fname = row.get("filename", "")
            pep_file.setdefault((pep, fname), 0.0)
            pep_file[(pep, fname)] += q
            if pep not in pep_proteins:
                pep_proteins[pep] = _split_proteins(row.get("proteins", ""))

        for (pep, fname), q in pep_file.items():
            proteins = pep_proteins.get(pep, [])
            if rollup == "unique" and len(proteins) != 1:
                continue
            for prot in proteins:
                collected.setdefault(prot, {}).setdefault(fname, []).append(q)

        matrix = {}
        for prot, files in collected.items():
            per_file = {}
            for fname, quantities in files.items():
                if rollup == "median":
                    per_file[fname] = statistics.median(quantities)
                else:  # sum, unique
                    per_file[fname] = sum(quantities)
            matrix[prot] = per_file
        return matrix

    def matrix(self, level="peptide", rollup="sum"):
        """Feature matrix for ``level`` (``peptide``/``protein``). ``rollup`` is
        ignored for the peptide level. Cached."""
        key = ("peptide", None) if level == "peptide" else ("protein", rollup)
        if key not in self._cache:
            if level == "peptide":
                self._cache[key] = self._peptide_matrix()
            else:
                self._cache[key] = self._protein_matrix(rollup)
        return self._cache[key]

    def proteins_for_peptide(self, peptide_plain):
        for row in self.session.all_quant_rows():
            if row.get("peptide_plain", "") == peptide_plain:
                return _split_proteins(row.get("proteins", ""))
        return []
