#!/usr/bin/env python3

import argparse
import csv
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Paths:
    mzml_dir: Path
    sage_dir: Path
    percolator_dir: Path
    out_dir: Path
    by_file_dir: Path
    sage_psms: Path
    sage_lfq: Path
    percolator_psms: Path
    percolator_peptides: Path
    percolator_proteins: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reorganize Sage and Percolator outputs into compact file-first scan, peptide, protein, and quantity tables."
    )

    parser.add_argument("--mzml-dir", required=True, type=Path)
    parser.add_argument("--sage-dir", required=True, type=Path)
    parser.add_argument("--percolator-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--q-max", type=float, default=0.01)
    parser.add_argument("--include-zero-quant", action="store_true")

    return parser.parse_args()


def make_paths(args: argparse.Namespace) -> Paths:
    mzml_dir = args.mzml_dir.resolve()
    sage_dir = args.sage_dir.resolve()
    percolator_dir = args.percolator_dir.resolve()
    out_dir = args.out_dir.resolve()

    return Paths(
        mzml_dir=mzml_dir,
        sage_dir=sage_dir,
        percolator_dir=percolator_dir,
        out_dir=out_dir,
        by_file_dir=out_dir / "by_file",
        sage_psms=sage_dir / "results.sage.tsv",
        sage_lfq=sage_dir / "lfq.tsv",
        percolator_psms=percolator_dir / "psms.tsv",
        percolator_peptides=percolator_dir / "peptides.tsv",
        percolator_proteins=percolator_dir / "proteins.tsv",
    )


def require_files(paths: Paths) -> None:
    inputs = [
        paths.sage_psms,
        paths.sage_lfq,
        paths.percolator_psms,
        paths.percolator_peptides,
        paths.percolator_proteins,
    ]

    missing = [str(path) for path in inputs if not path.exists()]

    if missing:
        raise SystemExit("missing input files:\n" + "\n".join(missing))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def f_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None

    try:
        return float(value)
    except ValueError:
        return None


def keep_q(row: dict[str, str], field: str, q_max: float) -> bool:
    value = f_float(row.get(field))

    if value is None:
        return False

    return value <= q_max


def split_items(value: str | None) -> list[str]:
    if not value:
        return []

    items = []

    for item in re.split(r"[;, ]+", value):
        item = item.strip()

        if item:
            items.append(item)

    return items


def is_decoy(value: str | None) -> bool:
    return any(item.startswith("decoy_") for item in split_items(value))


def strip_flanks(peptide: str | None) -> str:
    if not peptide:
        return ""

    if len(peptide) >= 5 and peptide[1] == "." and peptide[-2] == ".":
        return peptide[2:-2]

    return peptide


def strip_mods(peptide: str | None) -> str:
    value = strip_flanks(peptide)
    value = re.sub(r"\[[^\]]*\]", "", value)
    value = re.sub(r"\([^\)]*\)", "", value)
    value = re.sub(r"\{[^\}]*\}", "", value)
    value = re.sub(r"[^A-Z]", "", value)

    return value


def left_flank(peptide: str | None) -> str:
    if peptide and len(peptide) >= 5 and peptide[1] == "." and peptide[-2] == ".":
        return peptide[0]

    return ""


def right_flank(peptide: str | None) -> str:
    if peptide and len(peptide) >= 5 and peptide[1] == "." and peptide[-2] == ".":
        return peptide[-1]

    return ""


def scan_number(scan_native: str | None) -> str:
    if not scan_native:
        return ""

    match = re.search(r"scan=(\d+)", scan_native)

    if match:
        return match.group(1)

    return scan_native


def scan_key(filename: str, scan: str) -> str:
    return f"{filename}:scan={scan}"


def safe_dir_name(filename: str) -> str:
    name = filename

    for suffix in [".centroid.mzML", ".centroid.mzml", ".mzML", ".mzml"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break

    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def mzml_path(paths: Paths, filename: str) -> str:
    path = paths.mzml_dir / filename

    if path.exists():
        return str(path)

    return ""


def fraction_name(filename: str) -> str:
    if "_GuHCl_" in filename:
        return "GuHCl"

    if "_NaCl_" in filename:
        return "NaCl"

    return ""


def replicate_name(filename: str) -> str:
    match = re.search(r"_(\d+)(?:_|\.centroid)", filename)

    if match:
        return match.group(1)

    return ""


def sort_key(value: str) -> tuple[int, int | str]:
    if value.isdigit():
        return (0, int(value))

    return (1, value)


def join_unique(values: list[str]) -> str:
    items = sorted(set(value for value in values if value), key=sort_key)
    return ";".join(items)


def protein_aliases(protein_id: str) -> set[str]:
    aliases = {protein_id}

    for item in split_items(protein_id):
        aliases.add(item)

    return aliases


def best_float(rows: list[dict[str, Any]], field: str, default: str = "") -> str:
    values = []

    for row in rows:
        value = f_float(str(row.get(field, "")))

        if value is not None:
            values.append(value)

    if not values:
        return default

    return str(min(values))


def best_score(rows: list[dict[str, Any]], field: str, default: str = "") -> str:
    values = []

    for row in rows:
        value = f_float(str(row.get(field, "")))

        if value is not None:
            values.append(value)

    if not values:
        return default

    return str(max(values))


def build_file_rows(paths: Paths, filenames: list[str]) -> list[dict[str, Any]]:
    rows = []

    for filename in filenames:
        rows.append({
            "filename": filename,
            "run_dir": safe_dir_name(filename),
            "mzml_path": mzml_path(paths, filename),
            "fraction": fraction_name(filename),
            "replicate": replicate_name(filename),
        })

    return rows


def build_psms(
    paths: Paths,
    sage_psms: list[dict[str, str]],
    percolator_psms: list[dict[str, str]],
    q_max: float,
) -> list[dict[str, Any]]:
    sage_by_id = {row.get("psm_id", ""): row for row in sage_psms}
    rows = []

    for row in percolator_psms:
        if not keep_q(row, "q-value", q_max):
            continue

        if is_decoy(row.get("proteinIds", "")):
            continue

        psm_id = row.get("PSMId", "")
        sage = sage_by_id.get(psm_id, {})

        filename = row.get("filename", "") or sage.get("filename", "")
        run_dir = safe_dir_name(filename)
        scan_native = sage.get("scannr", "")
        scan = scan_number(scan_native)

        peptide_flanked = row.get("peptide", "")
        peptide = strip_flanks(peptide_flanked)
        peptide_plain = strip_mods(peptide_flanked)

        rows.append({
            "filename": filename,
            "run_dir": run_dir,
            "mzml_path": mzml_path(paths, filename),
            "scan": scan,
            "scan_native": scan_native,
            "scan_key": scan_key(filename, scan),
            "psm_id": psm_id,

            "peptide": peptide,
            "peptide_plain": peptide_plain,
            "peptide_flanked": peptide_flanked,
            "flank_left": left_flank(peptide_flanked),
            "flank_right": right_flank(peptide_flanked),
            "proteins": row.get("proteinIds", "") or sage.get("proteins", ""),

            "percolator_score": row.get("score", ""),
            "percolator_q": row.get("q-value", ""),
            "percolator_pep": row.get("posterior_error_prob", ""),

            "sage_spectrum_q": sage.get("spectrum_q", ""),
            "sage_peptide_q": sage.get("peptide_q", ""),
            "sage_protein_q": sage.get("protein_q", ""),
            "sage_discriminant_score": sage.get("sage_discriminant_score", ""),
            "sage_posterior_error": sage.get("posterior_error", ""),

            "rank": sage.get("rank", ""),
            "label": sage.get("label", ""),
            "charge": sage.get("charge", ""),
            "exp_mass": sage.get("expmass", ""),
            "calc_mass": sage.get("calcmass", ""),
            "precursor_ppm": sage.get("precursor_ppm", ""),
            "fragment_ppm": sage.get("fragment_ppm", ""),
            "rt": sage.get("rt", ""),
            "aligned_rt": sage.get("aligned_rt", ""),
            "predicted_rt": sage.get("predicted_rt", ""),
            "delta_rt_model": sage.get("delta_rt_model", ""),
            "hyperscore": sage.get("hyperscore", ""),
            "delta_next": sage.get("delta_next", ""),
            "delta_best": sage.get("delta_best", ""),
            "matched_peaks": sage.get("matched_peaks", ""),
            "matched_intensity_pct": sage.get("matched_intensity_pct", ""),
            "ms2_intensity": sage.get("ms2_intensity", ""),
        })

    rows.sort(key=lambda x: (x["filename"], int(x["scan"]) if str(x["scan"]).isdigit() else 0, x["psm_id"]))
    return rows


def build_quant(
    paths: Paths,
    sage_lfq: list[dict[str, str]],
    q_max: float,
    include_zero_quant: bool,
) -> list[dict[str, Any]]:
    if not sage_lfq:
        return []

    base = {
        "peptide",
        "charge",
        "proteins",
        "q_value",
        "score",
        "spectral_angle",
    }

    filenames = [field for field in sage_lfq[0].keys() if field not in base]
    rows = []

    for row in sage_lfq:
        if not keep_q(row, "q_value", q_max):
            continue

        if is_decoy(row.get("proteins", "")):
            continue

        peptide = row.get("peptide", "")
        peptide_plain = strip_mods(peptide)

        for filename in filenames:
            quantity = row.get(filename, "")

            if not include_zero_quant and f_float(quantity) == 0.0:
                continue

            rows.append({
                "filename": filename,
                "run_dir": safe_dir_name(filename),
                "mzml_path": mzml_path(paths, filename),
                "peptide": peptide,
                "peptide_plain": peptide_plain,
                "charge": row.get("charge", ""),
                "proteins": row.get("proteins", ""),
                "quantity": quantity,
                "sage_lfq_q": row.get("q_value", ""),
                "sage_lfq_score": row.get("score", ""),
                "spectral_angle": row.get("spectral_angle", ""),
            })

    rows.sort(key=lambda x: (x["filename"], x["peptide"], x["charge"]))
    return rows


def build_global_peptides(
    percolator_peptides: list[dict[str, str]],
    psms: list[dict[str, Any]],
    quant: list[dict[str, Any]],
    q_max: float,
) -> list[dict[str, Any]]:
    psms_by_peptide: dict[str, list[dict[str, Any]]] = defaultdict(list)
    quant_by_peptide: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in psms:
        psms_by_peptide[row["peptide"]].append(row)

    for row in quant:
        quant_by_peptide[row["peptide"]].append(row)

    rows = []

    for row in percolator_peptides:
        if not keep_q(row, "q-value", q_max):
            continue

        if is_decoy(row.get("proteinIds", "")):
            continue

        peptide_flanked = row.get("peptide", "")
        peptide = strip_flanks(peptide_flanked)
        peptide_plain = strip_mods(peptide_flanked)

        linked_psms = psms_by_peptide.get(peptide, [])
        linked_quant = quant_by_peptide.get(peptide, [])

        rows.append({
            "peptide": peptide,
            "peptide_plain": peptide_plain,
            "example_flanked_peptide": peptide_flanked,
            "proteins": row.get("proteinIds", ""),
            "best_psm_id": row.get("PSMId", ""),
            "best_filename": row.get("filename", ""),
            "percolator_score": row.get("score", ""),
            "percolator_q": row.get("q-value", ""),
            "percolator_pep": row.get("posterior_error_prob", ""),
            "n_psms": len(linked_psms),
            "n_files": len(set(x["filename"] for x in linked_psms if x["filename"])),
            "files": join_unique([x["filename"] for x in linked_psms]),
            "scan_keys": join_unique([x["scan_key"] for x in linked_psms]),
            "psm_ids": join_unique([x["psm_id"] for x in linked_psms]),
            "has_lfq": "1" if linked_quant else "0",
        })

    rows.sort(key=lambda x: (f_float(str(x["percolator_q"])) or 1.0, x["peptide"]))
    return rows


def build_global_proteins(
    percolator_proteins: list[dict[str, str]],
    psms: list[dict[str, Any]],
    q_max: float,
) -> list[dict[str, Any]]:
    psms_by_protein: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in psms:
        for protein in split_items(row["proteins"]):
            psms_by_protein[protein].append(row)

    rows = []

    for row in percolator_proteins:
        if not keep_q(row, "q-value", q_max):
            continue

        protein_id = row.get("ProteinId", "")

        if is_decoy(protein_id):
            continue

        linked_psms = []
        seen = set()

        for alias in protein_aliases(protein_id):
            for psm in psms_by_protein.get(alias, []):
                psm_id = psm["psm_id"]

                if psm_id in seen:
                    continue

                seen.add(psm_id)
                linked_psms.append(psm)

        peptides = split_items(row.get("peptideIds", ""))

        rows.append({
            "protein_id": protein_id,
            "protein_group_id": row.get("ProteinGroupId", ""),
            "percolator_q": row.get("q-value", ""),
            "percolator_pep": row.get("posterior_error_prob", ""),
            "n_percolator_peptides": len(peptides),
            "n_linked_psms": len(linked_psms),
            "n_files": len(set(x["filename"] for x in linked_psms if x["filename"])),
            "files": join_unique([x["filename"] for x in linked_psms]),
            "scan_keys": join_unique([x["scan_key"] for x in linked_psms]),
            "psm_ids": join_unique([x["psm_id"] for x in linked_psms]),
            "peptides": ";".join(peptides),
        })

    rows.sort(key=lambda x: (f_float(str(x["percolator_q"])) or 1.0, x["protein_id"]))
    return rows


def build_file_peptides(
    filename: str,
    psms: list[dict[str, Any]],
    quant: list[dict[str, Any]],
    global_peptides: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    psms_by_peptide: dict[str, list[dict[str, Any]]] = defaultdict(list)
    quant_by_peptide: dict[str, list[dict[str, Any]]] = defaultdict(list)
    global_by_peptide = {row["peptide"]: row for row in global_peptides}

    for row in psms:
        if row["filename"] == filename:
            psms_by_peptide[row["peptide"]].append(row)

    for row in quant:
        if row["filename"] == filename:
            quant_by_peptide[row["peptide"]].append(row)

    peptides = sorted(set(psms_by_peptide) | set(quant_by_peptide))
    rows = []

    for peptide in peptides:
        linked_psms = psms_by_peptide.get(peptide, [])
        linked_quant = quant_by_peptide.get(peptide, [])
        global_row = global_by_peptide.get(peptide, {})

        rows.append({
            "peptide": peptide,
            "peptide_plain": strip_mods(peptide),
            "proteins": join_unique([x["proteins"] for x in linked_psms] + [x["proteins"] for x in linked_quant]),
            "n_psms": len(linked_psms),
            "scans": join_unique([x["scan"] for x in linked_psms]),
            "psm_ids": join_unique([x["psm_id"] for x in linked_psms]),
            "best_percolator_q": best_float(linked_psms, "percolator_q", global_row.get("percolator_q", "")),
            "best_percolator_score": best_score(linked_psms, "percolator_score", global_row.get("percolator_score", "")),
            "quantities": join_unique([x["quantity"] for x in linked_quant]),
            "charges": join_unique([x["charge"] for x in linked_quant]),
            "sage_lfq_q": best_float(linked_quant, "sage_lfq_q"),
            "sage_lfq_score": best_score(linked_quant, "sage_lfq_score"),
            "spectral_angle": best_score(linked_quant, "spectral_angle"),
        })

    rows.sort(key=lambda x: (x["peptide"], x["best_percolator_q"]))
    return rows


def build_file_proteins(
    filename: str,
    psms: list[dict[str, Any]],
    global_proteins: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    global_by_protein = {}

    for row in global_proteins:
        for alias in protein_aliases(row["protein_id"]):
            global_by_protein[alias] = row

    protein_psms: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for row in psms:
        if row["filename"] != filename:
            continue

        for protein in split_items(row["proteins"]):
            protein_psms[protein].append(row)

    rows = []

    for protein_id, linked_psms in sorted(protein_psms.items()):
        global_row = global_by_protein.get(protein_id, {})
        peptides = sorted(set(x["peptide"] for x in linked_psms if x["peptide"]))

        rows.append({
            "protein_id": protein_id,
            "protein_group_id": global_row.get("protein_group_id", ""),
            "protein_q": global_row.get("percolator_q", ""),
            "protein_pep": global_row.get("percolator_pep", ""),
            "n_peptides": len(peptides),
            "n_psms": len(linked_psms),
            "peptides": ";".join(peptides),
            "scans": join_unique([x["scan"] for x in linked_psms]),
            "psm_ids": join_unique([x["psm_id"] for x in linked_psms]),
        })

    rows.sort(key=lambda x: (f_float(str(x["protein_q"])) or 1.0, x["protein_id"]))
    return rows


def build_file_summary(
    file_row: dict[str, Any],
    file_psms: list[dict[str, Any]],
    file_peptides: list[dict[str, Any]],
    file_proteins: list[dict[str, Any]],
    file_quant: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "filename": file_row["filename"],
        "run_dir": file_row["run_dir"],
        "mzml_path": file_row["mzml_path"],
        "fraction": file_row["fraction"],
        "replicate": file_row["replicate"],
        "n_psms": len(file_psms),
        "n_peptides": len(file_peptides),
        "n_proteins": len(file_proteins),
        "n_quant_rows": len(file_quant),
        "n_scans": len(set(row["scan"] for row in file_psms if row["scan"])),
    }


def main() -> None:
    args = parse_args()
    paths = make_paths(args)

    require_files(paths)
    paths.out_dir.mkdir(parents=True, exist_ok=True)
    paths.by_file_dir.mkdir(parents=True, exist_ok=True)

    sage_psms = read_tsv(paths.sage_psms)
    sage_lfq = read_tsv(paths.sage_lfq)
    percolator_psms = read_tsv(paths.percolator_psms)
    percolator_peptides = read_tsv(paths.percolator_peptides)
    percolator_proteins = read_tsv(paths.percolator_proteins)

    filenames = sorted(
        set(row.get("filename", "") for row in sage_psms)
        | set(row.get("filename", "") for row in percolator_psms)
        | set(row.get("filename", "") for row in percolator_peptides)
    )
    filenames = [x for x in filenames if x]

    file_rows = build_file_rows(paths, filenames)
    psms = build_psms(paths, sage_psms, percolator_psms, args.q_max)
    quant = build_quant(paths, sage_lfq, args.q_max, args.include_zero_quant)
    global_peptides = build_global_peptides(percolator_peptides, psms, quant, args.q_max)
    global_proteins = build_global_proteins(percolator_proteins, psms, args.q_max)

    by_file_psm_fields = [
        "scan",
        "scan_native",
        "psm_id",
        "peptide",
        "peptide_plain",
        "peptide_flanked",
        "flank_left",
        "flank_right",
        "proteins",
        "percolator_score",
        "percolator_q",
        "percolator_pep",
        "sage_spectrum_q",
        "sage_peptide_q",
        "sage_protein_q",
        "sage_discriminant_score",
        "sage_posterior_error",
        "rank",
        "label",
        "charge",
        "exp_mass",
        "calc_mass",
        "precursor_ppm",
        "fragment_ppm",
        "rt",
        "aligned_rt",
        "predicted_rt",
        "delta_rt_model",
        "hyperscore",
        "delta_next",
        "delta_best",
        "matched_peaks",
        "matched_intensity_pct",
        "ms2_intensity",
    ]

    by_file_peptide_fields = [
        "peptide",
        "peptide_plain",
        "proteins",
        "n_psms",
        "scans",
        "psm_ids",
        "best_percolator_q",
        "best_percolator_score",
        "quantities",
        "charges",
        "sage_lfq_q",
        "sage_lfq_score",
        "spectral_angle",
    ]

    by_file_protein_fields = [
        "protein_id",
        "protein_group_id",
        "protein_q",
        "protein_pep",
        "n_peptides",
        "n_psms",
        "peptides",
        "scans",
        "psm_ids",
    ]

    by_file_quant_fields = [
        "peptide",
        "peptide_plain",
        "charge",
        "proteins",
        "quantity",
        "sage_lfq_q",
        "sage_lfq_score",
        "spectral_angle",
    ]

    file_summary_rows = []

    for file_row in file_rows:
        filename = file_row["filename"]
        run_dir = paths.by_file_dir / file_row["run_dir"]

        file_psms = [row for row in psms if row["filename"] == filename]
        file_quant = [row for row in quant if row["filename"] == filename]
        file_peptides = build_file_peptides(filename, psms, quant, global_peptides)
        file_proteins = build_file_proteins(filename, psms, global_proteins)
        file_summary = build_file_summary(file_row, file_psms, file_peptides, file_proteins, file_quant)

        reduced_psms = [
            {field: row.get(field, "") for field in by_file_psm_fields}
            for row in file_psms
        ]

        reduced_quant = [
            {field: row.get(field, "") for field in by_file_quant_fields}
            for row in file_quant
        ]

        write_tsv(run_dir / "scan_lookup.tsv", reduced_psms, by_file_psm_fields)
        write_tsv(run_dir / "psms.tsv", reduced_psms, by_file_psm_fields)
        write_tsv(run_dir / "peptides.tsv", file_peptides, by_file_peptide_fields)
        write_tsv(run_dir / "proteins.tsv", file_proteins, by_file_protein_fields)
        write_tsv(run_dir / "peptide_quant.tsv", reduced_quant, by_file_quant_fields)
        write_json(run_dir / "file.json", file_summary)

        file_summary_rows.append(file_summary)

    write_tsv(paths.out_dir / "files.tsv", file_summary_rows, [
        "filename",
        "run_dir",
        "mzml_path",
        "fraction",
        "replicate",
        "n_psms",
        "n_peptides",
        "n_proteins",
        "n_quant_rows",
        "n_scans",
    ])

    write_tsv(paths.out_dir / "peptides.tsv", global_peptides, [
        "peptide",
        "peptide_plain",
        "example_flanked_peptide",
        "proteins",
        "best_psm_id",
        "best_filename",
        "percolator_score",
        "percolator_q",
        "percolator_pep",
        "n_psms",
        "n_files",
        "files",
        "scan_keys",
        "psm_ids",
        "has_lfq",
    ])

    write_tsv(paths.out_dir / "proteins.tsv", global_proteins, [
        "protein_id",
        "protein_group_id",
        "percolator_q",
        "percolator_pep",
        "n_percolator_peptides",
        "n_linked_psms",
        "n_files",
        "files",
        "scan_keys",
        "psm_ids",
        "peptides",
    ])

    manifest = {
        "created": datetime.now().isoformat(timespec="seconds"),
        "mzml_dir": str(paths.mzml_dir),
        "sage_dir": str(paths.sage_dir),
        "percolator_dir": str(paths.percolator_dir),
        "out_dir": str(paths.out_dir),
        "by_file_dir": str(paths.by_file_dir),
        "q_max": args.q_max,
        "include_zero_quant": args.include_zero_quant,
        "inputs": {
            "sage_psms": str(paths.sage_psms),
            "sage_lfq": str(paths.sage_lfq),
            "percolator_psms": str(paths.percolator_psms),
            "percolator_peptides": str(paths.percolator_peptides),
            "percolator_proteins": str(paths.percolator_proteins),
        },
        "outputs": {
            "files": str(paths.out_dir / "files.tsv"),
            "peptides": str(paths.out_dir / "peptides.tsv"),
            "proteins": str(paths.out_dir / "proteins.tsv"),
            "by_file": str(paths.by_file_dir),
        },
        "rows": {
            "files": len(file_summary_rows),
            "global_peptides": len(global_peptides),
            "global_proteins": len(global_proteins),
            "psms": len(psms),
            "quant": len(quant),
        },
    }

    write_json(paths.out_dir / "manifest.json", manifest)


if __name__ == "__main__":
    main()
