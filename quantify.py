#!/usr/bin/env python3
"""AUC-based peptide/protein quantification + motif grouping (off-GUI stage).

This is the pipeline counterpart of the viewer's quantitative comparisons. It
replaces the Sage-LFQ ``quantity`` column with the **charge-distribution AUC**
the GUI shows in the MS Data tab, rolls that up to proteins via their single
most abundant unique peptide, and groups proteins by shared skeleton motif.

Inputs (a reorganized project; see ``reorganize-results.py`` / ``index-motifs.py``):
    <project>/searches/reorganized/   files.tsv + by_file/<run>/psms.tsv
    <project>/distributions/          one <stem>.distributions.sqlite per file
    <motif index>/                    motifs.tsv, postings.bin, proteins.tsv

Outputs:
    <project>/quant/peptide_auc.tsv   filename, peptide_plain, proteins, unique, auc
    <project>/quant/protein_quant.tsv filename, protein, quantity, best_peptide
    <project>/quant/manifest.json
    <motif-sets>/motifs.tsv           motif_id, motif_text, protein_count,
                                      observed_count, observed_accessions
    <motif-sets>/motif_quant.tsv      filename, motif_id, quantity
    <motif-sets>/manifest.json

Quantification, matching the GUI (``ms_viewer_tab._select_distribution_for_candidate``):
    mono m/z   = neutral_mass / z + PROTON          (neutral from search ``calc_mass``)
    m/z window = mono_mz +/- max(mono_mz * ppm/1e6, mz_floor)   ppm=10, floor=0.02 Th
    rt window  = distribution RT band overlaps [rt - rt_half, rt + rt_half]  rt_half=0.8
    charge     = exact
    best pick  = min(|dist.mono_mz - theo|, |dist.rt_apex - rt|)
    AUC        = SUM(features.area) over the whole analyte (all charge states of the
                 matched distribution); falls back to the single distribution's
                 members when the distribution has no analyte.

A peptide's per-file AUC sums the AUCs of the distinct analytes/distributions its
PSMs matched (so multiple charges of the same species are not double counted). A
protein's per-file quantity is the per-file AUC of its most abundant unique
peptide (the unique peptide with the greatest total AUC across files), so the
same peptide is tracked across samples. A motif's per-file quantity is the SUM of
its observed member proteins' quantities.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import sys
from pathlib import Path

# Proton mass; identical to viewer/session.py and chemistry.PROTON.
PROTON = 1.007276554940804

# GUI defaults (ms_viewer_tab.py): precursor_ppm=10, 0.02 Th floor, rt_half=0.8.
DEFAULT_PPM = 10.0
DEFAULT_MZ_FLOOR = 0.02
DEFAULT_RT_HALF = 0.8
DEFAULT_Q_MAX = 0.01


# --------------------------------------------------------------------------
# small IO helpers
# --------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fields})


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(obj, f, indent=2)


def _to_float(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def split_proteins(value: str) -> list[str]:
    if not value:
        return []
    return [p for p in re.split(r"[;,\s]+", value.strip()) if p]


# --------------------------------------------------------------------------
# per-file distribution AUC lookup
# --------------------------------------------------------------------------

class DistributionAUCIndex:
    """In-memory view of one file's distributions sqlite, keyed for matching.

    Loads distribution geometry (mono_mz, rt band, charge) plus the AUC of each
    distribution and each analyte, so peptide matching is a pure-Python scan with
    no per-candidate SQL.
    """

    def __init__(self, conn: sqlite3.Connection):
        self.by_charge: dict[int, list[dict]] = {}
        self.dist_auc: dict[int, float] = {}
        self.analyte_auc: dict[int, float] = {}
        self.analyte_of: dict[int, int] = {}
        self._load(conn)

    @classmethod
    def open(cls, path: Path):
        conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        try:
            return cls(conn)
        finally:
            conn.close()

    def _load(self, conn: sqlite3.Connection) -> None:
        # Per-distribution AUC = sum of member-feature areas.
        for row in conn.execute(
            """
            SELECT m.distribution_id AS did, SUM(f.area) AS auc
            FROM distribution_members m
            JOIN features f ON f.feature_id = m.feature_id
            GROUP BY m.distribution_id
            """
        ):
            self.dist_auc[row["did"]] = row["auc"] or 0.0

        # distribution -> analyte, and per-analyte AUC (all charge states).
        if _has_table(conn, "analyte_members"):
            for row in conn.execute(
                "SELECT analyte_id, distribution_id FROM analyte_members"
            ):
                self.analyte_of[row["distribution_id"]] = row["analyte_id"]
            for row in conn.execute(
                """
                SELECT am.analyte_id AS aid, SUM(f.area) AS auc
                FROM analyte_members am
                JOIN distribution_members m ON m.distribution_id = am.distribution_id
                JOIN features f ON f.feature_id = m.feature_id
                GROUP BY am.analyte_id
                """
            ):
                self.analyte_auc[row["aid"]] = row["auc"] or 0.0

        for row in conn.execute(
            "SELECT distribution_id, charge, mono_mz, rt_start, rt_end, rt_apex, "
            "quality FROM distributions"
        ):
            self.by_charge.setdefault(int(row["charge"] or 0), []).append({
                "distribution_id": row["distribution_id"],
                "mono_mz": row["mono_mz"],
                "rt_start": row["rt_start"],
                "rt_end": row["rt_end"],
                "rt_apex": row["rt_apex"],
            })

    def match(self, mono_mz: float, rt: float, charge: int,
              ppm: float, mz_floor: float, rt_half: float):
        """Return (key, auc) for the best matching species, or None.

        ``key`` is ``("analyte", id)`` or ``("dist", id)`` so a caller can
        de-duplicate species that several charges/PSMs map onto.
        """
        tol = max(mono_mz * ppm / 1e6, mz_floor)
        lo, hi = mono_mz - tol, mono_mz + tol
        rt_lo, rt_hi = rt - rt_half, rt + rt_half

        best = None
        best_key = None
        for d in self.by_charge.get(int(charge), ()):
            d_mz = d["mono_mz"]
            if d_mz is None or d_mz < lo or d_mz > hi:
                continue
            # RT band overlap, mirroring distributions_in_window.
            if d["rt_end"] is not None and d["rt_end"] < rt_lo:
                continue
            if d["rt_start"] is not None and d["rt_start"] > rt_hi:
                continue
            apex = d["rt_apex"] if d["rt_apex"] is not None else rt
            key = (abs(d_mz - mono_mz), abs(apex - rt))
            if best_key is None or key < best_key:
                best_key = key
                best = d

        if best is None:
            return None

        did = best["distribution_id"]
        aid = self.analyte_of.get(did)
        if aid is not None:
            return ("analyte", aid), self.analyte_auc.get(aid, 0.0)
        return ("dist", did), self.dist_auc.get(did, 0.0)


def _has_table(conn: sqlite3.Connection, name: str) -> bool:
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?", (name,)
    ).fetchone()
    return row is not None


# --------------------------------------------------------------------------
# peptide + protein quantification
# --------------------------------------------------------------------------

def neutral_mass_of(psm: dict):
    """Neutral monoisotopic mass for a PSM: prefer search ``calc_mass`` (includes
    modifications), fall back to ``exp_mass``."""
    return _to_float(psm.get("calc_mass")) or _to_float(psm.get("exp_mass"))


def is_target(psm: dict, q_max: float) -> bool:
    label = (psm.get("label") or "").strip()
    if label and label.lstrip("-").isdigit() and int(float(label)) < 0:
        return False  # decoy
    q = _to_float(psm.get("percolator_q"))
    if q is None:
        q = _to_float(psm.get("sage_peptide_q"))
    if q is not None and q > q_max:
        return False
    return True


def quantify_file(psms: list[dict], auc_index: DistributionAUCIndex,
                  ppm: float, mz_floor: float, rt_half: float, q_max: float):
    """Peptide -> AUC for one file.

    Returns ``{peptide_plain: {"auc": float, "proteins": [..], "unique": bool}}``.
    A peptide's AUC is the sum of the distinct matched species' AUCs across all
    its PSMs/charges (each analyte or distribution counted once).
    """
    # peptide_plain -> {"species": {key: auc}, "proteins": [...]}
    acc: dict[str, dict] = {}
    for psm in psms:
        if not is_target(psm, q_max):
            continue
        pep = psm.get("peptide_plain", "")
        if not pep:
            continue
        neutral = neutral_mass_of(psm)
        rt = _to_float(psm.get("rt"))
        charge = psm.get("charge")
        try:
            z = int(float(charge))
        except (TypeError, ValueError):
            z = 0
        if neutral is None or rt is None or z <= 0:
            continue

        entry = acc.setdefault(pep, {"species": {}, "proteins": None})
        if entry["proteins"] is None:
            entry["proteins"] = split_proteins(psm.get("proteins", ""))

        mono_mz = neutral / z + PROTON
        hit = auc_index.match(mono_mz, rt, z, ppm, mz_floor, rt_half)
        if hit is None:
            continue
        key, auc = hit
        # Keep max seen for a species key (identical across PSMs anyway).
        entry["species"][key] = auc

    out = {}
    for pep, entry in acc.items():
        proteins = entry["proteins"] or []
        total = sum(entry["species"].values())
        out[pep] = {
            "auc": total,
            "proteins": proteins,
            "unique": len(proteins) == 1,
        }
    return out


def roll_up_proteins(per_file_peptides: dict[str, dict]):
    """Choose one most-abundant unique peptide per protein, then read its per-file
    AUC as the protein's per-file quantity.

    ``per_file_peptides`` is ``{filename: {peptide_plain: {auc, proteins, unique}}}``.
    Returns ``(protein_matrix, best_peptide)`` where ``protein_matrix`` is
    ``{protein: {filename: quantity}}`` and ``best_peptide`` is
    ``{protein: peptide_plain}``.
    """
    # protein -> unique peptide -> total AUC across files, and per-file AUC.
    totals: dict[str, dict[str, float]] = {}
    per_file: dict[str, dict[str, dict[str, float]]] = {}
    for fname, peptides in per_file_peptides.items():
        for pep, info in peptides.items():
            if not info["unique"]:
                continue
            prot = info["proteins"][0]
            prot_totals = totals.setdefault(prot, {})
            prot_totals[pep] = prot_totals.get(pep, 0.0) + info["auc"]
            per_file.setdefault(prot, {}).setdefault(pep, {})[fname] = info["auc"]

    protein_matrix: dict[str, dict[str, float]] = {}
    best_peptide: dict[str, str] = {}
    for prot, pep_totals in totals.items():
        # Most abundant unique peptide; ties broken by sequence for determinism.
        best = max(pep_totals.items(), key=lambda kv: (kv[1], kv[0]))[0]
        best_peptide[prot] = best
        protein_matrix[prot] = dict(per_file[prot].get(best, {}))
    return protein_matrix, best_peptide


# --------------------------------------------------------------------------
# motif grouping + quantification
# --------------------------------------------------------------------------

def _decode_varints(buf):
    value = 0
    shift = 0
    for byte in buf:
        value |= (byte & 0x7F) << shift
        if byte & 0x80:
            shift += 7
        else:
            yield value
            value = 0
            shift = 0


def decode_posting_list(buf) -> list[int]:
    """Decode one motif's posting block (count-prefixed, delta-encoded).

    Mirrors ``viewer/motifs.decode_posting_list`` exactly.
    """
    it = _decode_varints(buf)
    try:
        count = next(it)
    except StopIteration:
        return []
    ids = []
    prev = 0
    for i, delta in enumerate(it):
        if i >= count:
            break
        prev = delta if i == 0 else prev + delta
        ids.append(prev)
    return ids


def load_motif_accessions(motif_dir: Path):
    """Yield ``(motif_id, motif_text, [accession, ...])`` for every motif."""
    proteins = read_tsv(motif_dir / "proteins.tsv")
    accession_by_id = {}
    for row in proteins:
        try:
            pid = int(row.get("protein_id", ""))
        except (TypeError, ValueError):
            continue
        accession_by_id[pid] = row.get("accession", str(pid))

    postings_path = motif_dir / "postings.bin"
    with postings_path.open("rb") as f:
        for row in read_tsv(motif_dir / "motifs.tsv"):
            try:
                motif_id = int(row.get("motif_id", ""))
                offset = int(row.get("posting_offset", "") or 0)
                nbytes = int(row.get("posting_bytes", "") or 0)
            except (TypeError, ValueError):
                continue
            if nbytes <= 0:
                yield motif_id, row.get("motif_text", ""), []
                continue
            f.seek(offset)
            pids = decode_posting_list(f.read(nbytes))
            accs = [accession_by_id.get(pid, str(pid)) for pid in pids]
            yield motif_id, row.get("motif_text", ""), accs


def quantify_motifs(motif_dir: Path, protein_matrix: dict[str, dict[str, float]],
                    filenames: list[str], min_observed: int = 2):
    """Group quantified proteins by motif and SUM their per-file quantities.

    Only motifs whose observed (quantified) protein set has more than one member
    are kept (``min_observed`` defaults to 2 -- "shows up more than once").

    Returns ``(motif_rows, motif_quant_rows)``.
    """
    quantified = set(protein_matrix)
    motif_rows = []
    motif_quant_rows = []
    for motif_id, motif_text, accessions in load_motif_accessions(motif_dir):
        observed = [a for a in accessions if a in quantified]
        if len(observed) < min_observed:
            continue
        motif_rows.append({
            "motif_id": motif_id,
            "motif_text": motif_text,
            "protein_count": len(accessions),
            "observed_count": len(observed),
            "observed_accessions": ";".join(observed),
        })
        for fname in filenames:
            total = sum(protein_matrix[a].get(fname, 0.0) for a in observed)
            if total > 0:
                motif_quant_rows.append({
                    "filename": fname,
                    "motif_id": motif_id,
                    "quantity": total,
                })
    return motif_rows, motif_quant_rows


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def distributions_sqlite_for(distributions_dir: Path, filename: str):
    """Resolve a file's distributions sqlite (mirrors session.distributions_db_for)."""
    stem = filename
    for suffix in (".centroid.mzML", ".centroid.mzml", ".mzML", ".mzml",
                   ".raw", ".RAW"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    exact = distributions_dir / f"{stem}.distributions.sqlite"
    if exact.exists():
        return exact
    matches = sorted(distributions_dir.glob(f"{stem}*.distributions.sqlite"))
    return matches[0] if matches else None


def run(reorganized: Path, distributions_dir: Path, out_dir: Path,
        motif_dir: Path = None, motif_out: Path = None,
        ppm: float = DEFAULT_PPM, mz_floor: float = DEFAULT_MZ_FLOOR,
        rt_half: float = DEFAULT_RT_HALF, q_max: float = DEFAULT_Q_MAX,
        min_observed: int = 2):
    files = read_tsv(reorganized / "files.tsv")
    if not files:
        raise SystemExit(f"no files.tsv under {reorganized}")

    by_file_dir = reorganized / "by_file"
    filenames = [r.get("filename", "") for r in files if r.get("filename")]

    per_file_peptides: dict[str, dict] = {}
    matched_counts = {}
    for file_row in files:
        filename = file_row.get("filename", "")
        run_dir = file_row.get("run_dir") or filename
        psms = read_tsv(by_file_dir / run_dir / "psms.tsv")

        sqlite_path = distributions_sqlite_for(distributions_dir, filename)
        if sqlite_path is None:
            print(f"  WARN no distributions sqlite for {filename}", file=sys.stderr)
            per_file_peptides[filename] = {}
            continue
        auc_index = DistributionAUCIndex.open(sqlite_path)
        peptides = quantify_file(psms, auc_index, ppm, mz_floor, rt_half, q_max)
        per_file_peptides[filename] = peptides
        matched_counts[filename] = sum(1 for p in peptides.values() if p["auc"] > 0)
        print(f"  {filename}: {matched_counts[filename]}/{len(peptides)} peptides matched")

    # ---- peptide table -------------------------------------------------
    peptide_rows = []
    for filename, peptides in per_file_peptides.items():
        for pep, info in sorted(peptides.items()):
            peptide_rows.append({
                "filename": filename,
                "peptide_plain": pep,
                "proteins": ";".join(info["proteins"]),
                "unique": "1" if info["unique"] else "0",
                "auc": repr(info["auc"]),
            })
    write_tsv(out_dir / "peptide_auc.tsv", peptide_rows,
              ["filename", "peptide_plain", "proteins", "unique", "auc"])

    # ---- protein table -------------------------------------------------
    protein_matrix, best_peptide = roll_up_proteins(per_file_peptides)
    protein_rows = []
    for prot in sorted(protein_matrix):
        for filename in filenames:
            q = protein_matrix[prot].get(filename)
            if q is None:
                continue
            protein_rows.append({
                "filename": filename,
                "protein": prot,
                "quantity": repr(q),
                "best_peptide": best_peptide.get(prot, ""),
            })
    write_tsv(out_dir / "protein_quant.tsv", protein_rows,
              ["filename", "protein", "quantity", "best_peptide"])

    write_json(out_dir / "manifest.json", {
        "reorganized": str(reorganized),
        "distributions_dir": str(distributions_dir),
        "n_files": len(filenames),
        "n_proteins": len(protein_matrix),
        "params": {"ppm": ppm, "mz_floor": mz_floor, "rt_half": rt_half,
                   "q_max": q_max, "auc": "analyte"},
    })
    print(f"quant: {len(protein_matrix)} proteins -> {out_dir}")

    # ---- motifs --------------------------------------------------------
    if motif_dir is not None:
        motif_out = motif_out or (out_dir.parent / "motif-sets")
        motif_rows, motif_quant_rows = quantify_motifs(
            motif_dir, protein_matrix, filenames, min_observed=min_observed)
        write_tsv(motif_out / "motifs.tsv", motif_rows,
                  ["motif_id", "motif_text", "protein_count",
                   "observed_count", "observed_accessions"])
        write_tsv(motif_out / "motif_quant.tsv", motif_quant_rows,
                  ["filename", "motif_id", "quantity"])
        write_json(motif_out / "manifest.json", {
            "motif_index": str(motif_dir),
            "n_motifs": len(motif_rows),
            "min_observed": min_observed,
        })
        print(f"motifs: {len(motif_rows)} groups (>{min_observed - 1} observed) "
              f"-> {motif_out}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project", type=Path,
                   help="project dir with searches/reorganized, distributions/")
    p.add_argument("--reorganized", type=Path, help="override reorganized dir")
    p.add_argument("--distributions", type=Path, help="override distributions dir")
    p.add_argument("--out-dir", type=Path, help="override quant output dir")
    p.add_argument("--motif-index", type=Path,
                   help="motif skeleton index dir (enables motif grouping)")
    p.add_argument("--motif-out", type=Path, help="override motif-sets output dir")
    p.add_argument("--ppm", type=float, default=DEFAULT_PPM)
    p.add_argument("--mz-floor", type=float, default=DEFAULT_MZ_FLOOR)
    p.add_argument("--rt-half", type=float, default=DEFAULT_RT_HALF)
    p.add_argument("--q-max", type=float, default=DEFAULT_Q_MAX)
    p.add_argument("--min-observed", type=int, default=2,
                   help="minimum observed proteins for a motif to be kept")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if args.reorganized:
        reorganized = args.reorganized
    elif args.project:
        reorganized = args.project / "searches" / "reorganized"
    else:
        raise SystemExit("provide --project or --reorganized")

    if args.distributions:
        distributions_dir = args.distributions
    elif args.project:
        distributions_dir = args.project / "distributions"
    else:
        raise SystemExit("provide --project or --distributions")

    if args.out_dir:
        out_dir = args.out_dir
    elif args.project:
        out_dir = args.project / "quant"
    else:
        out_dir = reorganized.parent.parent / "quant"

    run(reorganized=reorganized, distributions_dir=distributions_dir,
        out_dir=out_dir, motif_dir=args.motif_index, motif_out=args.motif_out,
        ppm=args.ppm, mz_floor=args.mz_floor, rt_half=args.rt_half,
        q_max=args.q_max, min_observed=args.min_observed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
