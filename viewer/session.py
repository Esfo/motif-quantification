import csv
import json
import re
import sqlite3
import sys
from pathlib import Path


PROTON = 1.007276554940804
C13_DELTA = 1.00335483507


# Reorganized proteins.tsv packs every peptide of a protein into one field,
# which can blow past csv's default 128 KB field cap. Raise it as high as the
# platform allows.
def _raise_csv_limit():
    limit = sys.maxsize

    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = int(limit // 10)


_raise_csv_limit()


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


def read_fasta(path):
    """Yield (header, sequence) pairs from a FASTA file (whitespace-collapsed)."""
    path = Path(path)
    if not path.exists():
        return []

    entries = []
    header = None
    seq = []

    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None and seq:
                entries.append((header, "".join(seq)))
            header = line[1:]
            seq = []
        else:
            seq.append(line)

    if header is not None and seq:
        entries.append((header, "".join(seq)))

    return entries


def protein_keys(header):
    """Identifier aliases for a FASTA header, matching execution.xsh so the
    reorganized proteins.tsv ``protein_id`` (a Sage/Percolator accession) maps
    back to its sequence. e.g. ``sp|P12345|NAME_HUMAN foo`` → the full header,
    the first token, and the ``|``-split parts (incl. the bare accession)."""
    first = header.split()[0]
    keys = {header, first}

    if "|" in first:
        parts = first.split("|")
        keys.update(parts)

    return {k for k in keys if k}


# In-silico digestion parameters mirror the Sage ``enzyme`` block in
# execution.xsh (trypsin: cleave after K/R but not before P, up to 2 missed
# cleavages, peptides 7–50 residues). These define which peptides the search
# *attempted* to identify — the outlined rectangles in the Proteins tab.
DIGEST_CLEAVE_AT = "KR"
DIGEST_RESTRICT = "P"
DIGEST_MISSED_CLEAVAGES = 2
DIGEST_MIN_LEN = 7
DIGEST_MAX_LEN = 50


def cleavage_sites(seq, cleave_at=DIGEST_CLEAVE_AT, restrict=DIGEST_RESTRICT):
    sites = [0]
    for i, aa in enumerate(seq[:-1]):
        if aa in cleave_at and (not restrict or seq[i + 1] not in restrict):
            sites.append(i + 1)
    sites.append(len(seq))
    return sites


def digest_peptides(seq, missed_cleavages=DIGEST_MISSED_CLEAVAGES,
                    min_len=DIGEST_MIN_LEN, max_len=DIGEST_MAX_LEN):
    """All peptides the tryptic search attempted, as ``(start, end, peptide)``.

    Ported from ``make_flank_index`` in execution.xsh: for each cleavage
    window allow up to ``missed_cleavages`` skipped sites, keep peptides within
    the search length bounds. ``start``/``end`` are 0-based residue indices into
    ``seq`` (half-open), so the caller can position them on the sequence."""
    if not seq:
        return []
    sites = cleavage_sites(seq)
    out = []
    seen = set()
    for i in range(len(sites) - 1):
        last = min(i + missed_cleavages + 1, len(sites) - 1)
        for j in range(i, last):
            start = sites[i]
            end = sites[j + 1]
            if end - start < min_len or end - start > max_len:
                continue
            key = (start, end)
            if key in seen:
                continue
            seen.add(key)
            out.append((start, end, seq[start:end]))
    return out


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
        self.reorganized = Path(reorganized).resolve() if reorganized else None
        self.is_empty = self.reorganized is None

        self.manifest = read_json(self.reorganized / "manifest.json") if self.reorganized else {}

        manifest_mzml_dir = self.manifest.get("mzml_dir")
        manifest_by_file_dir = self.manifest.get("by_file_dir")

        self.mzml_dir = Path(centroid_dir).resolve() if centroid_dir else None

        if self.mzml_dir is None and manifest_mzml_dir:
            self.mzml_dir = Path(manifest_mzml_dir).resolve()

        if manifest_by_file_dir:
            self.by_file_dir = Path(manifest_by_file_dir).resolve()
        elif self.reorganized:
            self.by_file_dir = self.reorganized / "by_file"
        else:
            self.by_file_dir = None

        self.profile_dir = Path(profile_dir).resolve() if profile_dir else None
        self.distribution_db = Path(distribution_db).resolve() if distribution_db else None

        self.file_rows = read_tsv(self.reorganized / "files.tsv") if self.reorganized else []
        self.file_by_name = {row.get("filename", ""): row for row in self.file_rows}

        self._psm_cache = {}
        self._file_table_cache = {}
        self._global_peptides = None
        self._global_proteins = None
        self._quant_cache = None
        self._fasta_index = None
        self._searches_manifest = None
        self._peptide_q_cache = {}

    def files(self):
        return self.file_rows

    @property
    def distributions_dir(self):
        # Mirrors viewer.main_window.find_distributions_db: <project>/distributions,
        # where <project> is the parent of searches/reorganized.
        if self.reorganized is None:
            return None
        return self.reorganized.parent.parent / "distributions"

    @property
    def quant_dir(self):
        # <project>/quant, written by quantify.py (AUC peptide/protein quant).
        if self.reorganized is None:
            return None
        return self.reorganized.parent.parent / "quant"

    @property
    def motif_sets_dir(self):
        # <project>/motif-sets, written by quantify.py (motif grouping + quant).
        if self.reorganized is None:
            return None
        return self.reorganized.parent.parent / "motif-sets"

    def distributions_db_for(self, filename):
        """Path to the distributions sqlite produced for THIS file, or None.

        index-distributions.py names each sqlite ``<stem>.distributions.sqlite``
        where ``stem`` strips a ``.centroid.mzML`` suffix; the viewer's file names
        may or may not carry that suffix, so try the exact name first and then a
        stem glob. Returns None when no sqlite matches (the file has no
        distributions of its own — we must not fall back to another file's).
        """
        dist_dir = self.distributions_dir
        if not filename or dist_dir is None or not dist_dir.is_dir():
            return None

        name = Path(filename).name
        for suffix in (".centroid.mzML", ".centroid.mzml", ".mzML", ".mzml"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        else:
            name = Path(name).stem

        exact = dist_dir / f"{name}.distributions.sqlite"
        if exact.exists():
            return exact

        matches = sorted(dist_dir.glob(f"{name}*.distributions.sqlite"))
        return matches[0] if matches else None

    def summary(self):
        return self.manifest

    # ---- protein sequences (FASTA) --------------------------------------

    @property
    def searches_dir(self):
        # reorganized = <project>/searches/reorganized, so searches = its parent.
        return self.reorganized.parent if self.reorganized else None

    def _searches_manifest_data(self):
        if self._searches_manifest is None:
            path = self.searches_dir / "manifest.json" if self.searches_dir else None
            self._searches_manifest = read_json(path) if path else {}
        return self._searches_manifest

    def _find_fasta(self):
        """Locate a protein FASTA for this project, most-portable first.

        The reorganized tables carry no sequences, so the Proteins tab needs the
        original FASTA. Prefer the in-project decoy FASTA (``searches/.decoy.fasta``,
        written by execution.xsh — its forward entries are the real proteins),
        then the paths recorded in ``searches/manifest.json``, then any ``*.fasta``
        beside the searches dir. Returns ``(path, decoy_tag)`` or ``(None, None)``.
        """
        manifest = self._searches_manifest_data()
        decoy_tag = manifest.get("decoy_tag") or "rev_"

        candidates = []
        if self.searches_dir:
            candidates.append(self.searches_dir / ".decoy.fasta")
        for key in ("decoy_fasta", "proteome_fasta"):
            value = manifest.get(key)
            if value:
                candidates.append(Path(value))
        if self.searches_dir and self.searches_dir.is_dir():
            candidates.extend(sorted(self.searches_dir.glob("*.fasta")))
            candidates.extend(sorted(self.searches_dir.glob("*.fa")))

        for path in candidates:
            try:
                if path and Path(path).exists():
                    return Path(path), decoy_tag
            except Exception:
                continue
        return None, None

    def fasta_index(self):
        """Map every protein-id alias → its amino-acid sequence.

        Decoy entries (headers starting with the decoy tag) are skipped so only
        real proteins resolve. Built once and cached; empty when no FASTA found."""
        if self._fasta_index is not None:
            return self._fasta_index

        index = {}
        path, decoy_tag = self._find_fasta()
        if path is not None:
            for header, seq in read_fasta(path):
                if decoy_tag and header.startswith(decoy_tag):
                    continue
                for key in protein_keys(header):
                    index.setdefault(key, seq)

        self._fasta_index = index
        return index

    @property
    def has_sequences(self):
        return bool(self.fasta_index())

    def protein_sequence(self, protein_id):
        """Amino-acid sequence for a reorganized ``protein_id``, or None.

        Tries the id as-is, then its ``|``-split parts (the reorganized id may be
        a bare accession while the FASTA header is ``sp|ACC|NAME``, or vice-versa)."""
        if not protein_id:
            return None
        index = self.fasta_index()
        seq = index.get(protein_id)
        if seq is not None:
            return seq
        for part in str(protein_id).split("|"):
            part = part.strip()
            if part and part in index:
                return index[part]
        return None

    def peptide_q_for_file(self, filename):
        """Map plain peptide sequence → best percolator q-value in this file.

        Drives the Proteins-tab colouring: an in-silico peptide that appears here
        is coloured by its q-value, one that doesn't stays uncoloured."""
        if filename in self._peptide_q_cache:
            return self._peptide_q_cache[filename]

        out = {}
        for row in self.file_peptides(filename or ""):
            plain = row.get("peptide_plain") or ""
            if not plain:
                continue
            q = safe_float(row.get("best_percolator_q"))
            if plain not in out:
                out[plain] = q
            elif q is not None and (out[plain] is None or q < out[plain]):
                out[plain] = q
        self._peptide_q_cache[filename] = out
        return out

    def global_peptides(self):
        if self._global_peptides is None:
            self._global_peptides = (
                read_tsv(self.reorganized / "peptides.tsv") if self.reorganized else []
            )

        return self._global_peptides

    def global_proteins(self):
        if self._global_proteins is None:
            self._global_proteins = (
                read_tsv(self.reorganized / "proteins.tsv") if self.reorganized else []
            )

        return self._global_proteins

    def _file_table(self, filename, name):
        if self.by_file_dir is None:
            return []

        key = (filename, name)

        if key in self._file_table_cache:
            return self._file_table_cache[key]

        file_row = self.file_by_name.get(filename, {})
        run_dir = file_row.get("run_dir") or safe_dir_name(filename)
        rows = read_tsv(self.by_file_dir / run_dir / name)

        self._file_table_cache[key] = rows
        return rows

    def file_peptides(self, filename):
        return self._file_table(filename, "peptides.tsv")

    def file_proteins(self, filename):
        return self._file_table(filename, "proteins.tsv")

    def file_quant(self, filename):
        return self._file_table(filename, "peptide_quant.tsv")

    def all_quant_rows(self):
        if self._quant_cache is not None:
            return self._quant_cache

        rows = []

        for file_row in self.file_rows:
            filename = file_row.get("filename", "")

            for row in self.file_quant(filename):
                merged = dict(row)
                merged["filename"] = filename
                rows.append(merged)

        self._quant_cache = rows
        return rows

    def quant_for_peptide(self, peptide_plain):
        """Quantities for a peptide across every file, keyed by filename.

        Returns {filename: total_quantity} summed over charge states, using the
        modification-stripped sequence so charge/mod variants collapse together.
        """
        totals = {}

        for row in self.all_quant_rows():
            if row.get("peptide_plain", "") != peptide_plain:
                continue

            quantity = safe_float(row.get("quantity"), 0.0) or 0.0
            filename = row.get("filename", "")
            totals[filename] = totals.get(filename, 0.0) + quantity

        return totals

    def load_psms(self, filename):
        if self.by_file_dir is None:
            return []

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
