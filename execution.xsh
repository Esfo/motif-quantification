#!/usr/bin/env xonsh

from pathlib import Path
from datetime import datetime
import json
import re


# ----------------------------
# edit these
# ----------------------------

project_dir = Path("/home/sfo/store/proteomics/PXD015057")
mzml_dir = project_dir / "mzML"
searches_dir = project_dir / "searches"

proteome_fasta = Path("/home/sfo/data/proteomics/fastas/human_uniprot_reviewed_plus_contaminants.fasta")

reorganize_py = Path("/home/sfo/motif-quantification/reorganize-results.py")
distributions_py = Path("/home/sfo/motif-quantification/index-distributions.py")

decoy_tag = "decoy_"
batch_size = 4

run_sage = True
patch_pins = True
run_percolator = True
run_reorganize = True
run_distributions = True


# ----------------------------
# output paths
# ----------------------------

sage_dir = searches_dir / "sage"
percolator_dir = searches_dir / "percolator"

decoy_fasta = searches_dir / ".decoy.fasta"
sage_json = searches_dir / "sage.json"
manifest_json = searches_dir / "manifest.json"

reorganized_dir = searches_dir / "reorganized"

searches_dir.mkdir(parents=True, exist_ok=True)
sage_dir.mkdir(parents=True, exist_ok=True)
percolator_dir.mkdir(parents=True, exist_ok=True)


# ----------------------------
# helpers
# ----------------------------

def read_fasta(path):
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


def write_fasta(path, entries):
    with path.open("w") as out:
        for header, seq in entries:
            out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")


def protein_keys(header):
    first = header.split()[0]
    keys = {header, first}

    if "|" in first:
        parts = first.split("|")
        keys.update(parts)

        if len(parts) >= 2:
            keys.add(parts[1])

    return keys


def clean_peptide(peptide):
    x = peptide

    if len(x) >= 5 and x[1] == "." and x[-2] == ".":
        x = x[2:-2]

    x = re.sub(r"\[[^\]]*\]", "", x)
    x = re.sub(r"\([^\)]*\)", "", x)
    x = re.sub(r"\{[^\}]*\}", "", x)
    x = re.sub(r"[^A-Z]", "", x)

    return x


def split_proteins(values):
    proteins = []

    for value in values:
        for item in re.split(r"[;,]", value):
            item = item.strip()

            if item:
                proteins.append(item)

    return proteins


def cleavage_sites(seq):
    sites = [0]

    for i, aa in enumerate(seq[:-1]):
        if aa in "KR" and seq[i + 1] != "P":
            sites.append(i + 1)

    sites.append(len(seq))
    return sites


def read_pin_peptides(pins):
    peptides = set()

    for pin in pins:
        with pin.open(errors="replace") as f:
            header = f.readline().rstrip("\n").split("\t")

            if "Peptide" in header:
                peptide_i = header.index("Peptide")
            elif "peptide" in header:
                peptide_i = header.index("peptide")
            else:
                raise SystemExit(f"no Peptide column found in {pin}")

            for line in f:
                parts = line.rstrip("\n").split("\t")

                if len(parts) <= peptide_i:
                    continue

                peptide = clean_peptide(parts[peptide_i])

                if peptide:
                    peptides.add(peptide)

    return peptides


def make_flank_index(entries, peptides, missed_cleavages=2, min_len=7, max_len=50):
    flanks = {}
    flanks_by_protein = {}

    for header, seq in entries:
        keys = protein_keys(header)
        sites = cleavage_sites(seq)

        for i in range(len(sites) - 1):
            last = min(i + missed_cleavages + 1, len(sites) - 1)

            for j in range(i, last):
                start = sites[i]
                end = sites[j + 1]
                peptide = seq[start:end]

                if len(peptide) < min_len or len(peptide) > max_len:
                    continue

                if peptide not in peptides:
                    continue

                left = seq[start - 1] if start > 0 else "-"
                right = seq[end] if end < len(seq) else "-"
                hit = (left, right, header, start)

                if peptide not in flanks:
                    flanks[peptide] = []

                flanks[peptide].append(hit)

                for key in keys:
                    item = (key, peptide)

                    if item not in flanks_by_protein:
                        flanks_by_protein[item] = []

                    flanks_by_protein[item].append(hit)

    for peptide in flanks:
        flanks[peptide] = sorted(flanks[peptide], key=lambda x: (x[2], x[3]))

    for item in flanks_by_protein:
        flanks_by_protein[item] = sorted(flanks_by_protein[item], key=lambda x: (x[2], x[3]))

    return flanks, flanks_by_protein


def choose_flank(peptide, proteins, flanks, flanks_by_protein):
    plain = clean_peptide(peptide)

    if not plain:
        return None, "empty_peptide"

    for protein in proteins:
        hits = flanks_by_protein.get((protein, plain), [])

        if hits:
            left, right, _, _ = hits[0]
            global_pairs = set((x[0], x[1]) for x in flanks.get(plain, []))
            status = "ok" if len(global_pairs) <= 1 else "chosen_from_listed_protein"
            return f"{left}.{peptide}.{right}", status

    hits = flanks.get(plain, [])

    if not hits:
        return None, "not_found"

    pairs = sorted(set((x[0], x[1]) for x in hits))

    if len(pairs) == 1:
        left, right = pairs[0]
        return f"{left}.{peptide}.{right}", "ok_fallback"

    left, right, _, _ = hits[0]
    return f"{left}.{peptide}.{right}", "chosen_from_fasta"


def patch_pin(pin, flanks, flanks_by_protein):
    tmp = pin.with_suffix(pin.suffix + ".tmp")
    report_path = pin.with_suffix(pin.suffix + ".flank_report.tsv")

    bad = 0

    with pin.open(errors="replace") as inp, tmp.open("w") as out, report_path.open("w") as report:
        header_line = inp.readline().rstrip("\n")
        header = header_line.split("\t")

        if "Peptide" in header:
            peptide_i = header.index("Peptide")
        elif "peptide" in header:
            peptide_i = header.index("peptide")
        else:
            raise SystemExit(f"no Peptide column found in {pin}")

        out.write(header_line + "\n")
        report.write("row\tstatus\tpeptide\tplain_peptide\tpatched_peptide\n")

        for row_n, line in enumerate(inp, start=2):
            line = line.rstrip("\n")
            parts = line.split("\t")

            if len(parts) <= peptide_i:
                bad += 1
                report.write(f"{row_n}\tbad_row\t\t\t\n")
                continue

            peptide = parts[peptide_i]

            if len(peptide) >= 5 and peptide[1] == "." and peptide[-2] == ".":
                out.write(line + "\n")
                report.write(f"{row_n}\talready_flanked\t{peptide}\t{clean_peptide(peptide)}\t{peptide}\n")
                continue

            proteins = split_proteins(parts[peptide_i + 1:])
            patched_peptide, status = choose_flank(peptide, proteins, flanks, flanks_by_protein)

            if patched_peptide is None:
                bad += 1
                report.write(f"{row_n}\t{status}\t{peptide}\t{clean_peptide(peptide)}\t\n")
                continue

            parts[peptide_i] = patched_peptide
            out.write("\t".join(parts) + "\n")
            report.write(f"{row_n}\t{status}\t{peptide}\t{clean_peptide(peptide)}\t{patched_peptide}\n")

    if bad:
        tmp.unlink(missing_ok=True)
        raise SystemExit(f"refusing to overwrite {pin}: {bad} rows had missing flanks; see {report_path}")

    tmp.replace(pin)


# ----------------------------
# inputs
# ----------------------------

mzmls = sorted(mzml_dir.glob("*.centroid.mzML")) + sorted(mzml_dir.glob("*.centroid.mzml"))

if not mzmls:
    raise SystemExit(f"no centroided mzML files found in {mzml_dir}")

for mzml in mzmls:
    if ".centroid." not in mzml.name:
        raise SystemExit(f"refusing non-centroid mzML: {mzml}")

if not proteome_fasta.exists():
    raise SystemExit(f"missing proteome_fasta: {proteome_fasta}")


# ----------------------------
# decoy fasta
# ----------------------------

entries = read_fasta(proteome_fasta)
entries = [(h, s) for h, s in entries if not h.startswith(decoy_tag)]

if not entries:
    raise SystemExit(f"no FASTA entries found in {proteome_fasta}")

fasta_entries = []

for header, seq in entries:
    fasta_entries.append((header, seq))

for header, seq in entries:
    fasta_entries.append((f"{decoy_tag}{header}", seq[::-1]))

write_fasta(decoy_fasta, fasta_entries)


# ----------------------------
# sage config
# ----------------------------

sage_config = {
    "database": {
        "fasta": str(decoy_fasta),
        "generate_decoys": False,
        "decoy_tag": decoy_tag,
        "enzyme": {
            "missed_cleavages": 2,
            "min_len": 7,
            "max_len": 50,
            "cleave_at": "KR",
            "restrict": "P",
            "c_terminal": True,
            "semi_enzymatic": False,
        },
        "peptide_min_mass": 500.0,
        "peptide_max_mass": 5000.0,
        "ion_kinds": ["b", "y"],
        "static_mods": {
            "C": 57.021464,
        },
        "variable_mods": {
            "M": [15.994915],
            "[": [42.010565],
            "^Q": [-17.026549],
            "^E": [-18.010565],
        },
        "max_variable_mods": 2,
    },
    "precursor_tol": {
        "ppm": [-10, 10],
    },
    "fragment_tol": {
        "ppm": [-20, 20],
    },
    "precursor_charge": [2, 5],
    "isotope_errors": [-1, 2],
    "deisotope": False,
    "chimera": False,
    "wide_window": False,
    "report_psms": 1,
    "min_peaks": 15,
    "max_peaks": 150,
    "min_matched_peaks": 6,
    "quant": {
        "lfq": True,
        "lfq_settings": {
            "peak_scoring": "Hybrid",
            "integration": "Sum",
            "spectral_angle": 0.7,
            "ppm_tolerance": 5.0,
        },
    },
    "output_directory": str(sage_dir),
}

sage_json.write_text(json.dumps(sage_config, indent=2) + "\n")


# ----------------------------
# manifest
# ----------------------------

mzml_paths = [str(x) for x in mzmls]

manifest = {
    "created": datetime.now().isoformat(timespec="seconds"),
    "project_dir": str(project_dir),
    "mzml_dir": str(mzml_dir),
    "searches_dir": str(searches_dir),
    "proteome_fasta": str(proteome_fasta),
    "decoy_fasta": str(decoy_fasta),
    "decoy_tag": decoy_tag,
    "fasta_entries": len(entries),
    "decoy_entries": len(entries),
    "sage_json": str(sage_json),
    "sage_dir": str(sage_dir),
    "percolator_dir": str(percolator_dir),
    "batch_size": batch_size,
    "mzmls": mzml_paths,
    "sage_command": [
        "sage",
        "--write-pin",
        "--batch-size",
        str(batch_size),
        str(sage_json),
        *mzml_paths,
    ],
}

manifest_json.write_text(json.dumps(manifest, indent=2) + "\n")


# ----------------------------
# sage
# ----------------------------

if run_sage:
    sage --write-pin --batch-size @(str(batch_size)) @(str(sage_json)) @(mzml_paths)


# ----------------------------
# patch pins
# ----------------------------

pins = sorted(sage_dir.glob("*.pin"))

if not pins:
    raise SystemExit(f"no PIN files found in {sage_dir}")

if patch_pins:
    peptides = read_pin_peptides(pins)
    flanks, flanks_by_protein = make_flank_index(read_fasta(decoy_fasta), peptides)

    for pin in pins:
        patch_pin(pin, flanks, flanks_by_protein)


# ----------------------------
# percolator
# ----------------------------

pin_paths = [str(x) for x in pins]

manifest["pins"] = pin_paths
manifest["pin_patch"] = {
    "enabled": patch_pins,
    "overwritten_in_place": True,
    "fasta": str(decoy_fasta),
    "reports": [str(x.with_suffix(x.suffix + ".flank_report.tsv")) for x in pins],
}

manifest["percolator_command"] = [
    "percolator",
    "-Y",
    "-P",
    decoy_tag,
    "--picked-protein",
    str(decoy_fasta),
    "--protein-enzyme",
    "trypsin",
    "--results-psms",
    str(percolator_dir / "psms.tsv"),
    "--decoy-results-psms",
    str(percolator_dir / "decoy_psms.tsv"),
    "--results-peptides",
    str(percolator_dir / "peptides.tsv"),
    "--decoy-results-peptides",
    str(percolator_dir / "decoy_peptides.tsv"),
    "--results-proteins",
    str(percolator_dir / "proteins.tsv"),
    "--decoy-results-proteins",
    str(percolator_dir / "decoy_proteins.tsv"),
    *pin_paths,
]

manifest_json.write_text(json.dumps(manifest, indent=2) + "\n")

if run_percolator:
    percolator \
        -Y \
        -P @(decoy_tag) \
        --picked-protein @(str(decoy_fasta)) \
        --protein-enzyme trypsin \
        --results-psms @(str(percolator_dir / "psms.tsv")) \
        --decoy-results-psms @(str(percolator_dir / "decoy_psms.tsv")) \
        --results-peptides @(str(percolator_dir / "peptides.tsv")) \
        --decoy-results-peptides @(str(percolator_dir / "decoy_peptides.tsv")) \
        --results-proteins @(str(percolator_dir / "proteins.tsv")) \
        --decoy-results-proteins @(str(percolator_dir / "decoy_proteins.tsv")) \
        @(pin_paths)

# ----------------------------
# reorganize outputs
# ----------------------------


if run_reorganize:
    python @(str(reorganize_py)) \
        --mzml-dir @(str(mzml_dir)) \
        --sage-dir @(str(sage_dir)) \
        --percolator-dir @(str(percolator_dir)) \
        --out-dir @(str(reorganized_dir)) \
        --q-max 0.01


# ----------------------------
# ms1 distribution determination -> project/distributions/<file>.distributions.sqlite
# ----------------------------

distributions_dir = project_dir / "distributions"

if run_distributions:
    # The rust detector is the default engine: it builds (incremental, ~instant
    # if unchanged) and runs the native binary, which also writes scan_points for
    # the fast GUI. Add --overwrite to regenerate existing sqlites.
    python @(str(distributions_py)) \
        --mzml-dir @(str(mzml_dir)) \
        --out-dir @(str(distributions_dir))
