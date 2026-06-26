"""Reader for the motif skeleton index (``index-motifs.py`` / Rust indexer).

Directory layout (e.g. ``.../motifs/human-proteome-skeletons``):
    build_info.tsv   key/value build parameters
    motifs.tsv       motif_id, motif_text, protein_count, posting_offset, posting_bytes
    postings.bin     per-motif protein-id lists, LEB128 varint, delta-encoded
    proteins.tsv     protein_id, accession, source, representative_protein_id,
                     representative_accession, header

A motif_text like ``A.....G`` is a fixed-length skeleton: letters are anchored
residues, dots are wildcards. Tab 4 represents groups of proteins by the motif
they share, so the core lookups here are motif -> protein ids and accession
-> motifs. Posting lists are decoded lazily from byte offsets so the whole
``postings.bin`` never has to load at once.
"""

import csv
from pathlib import Path


def _decode_varints(buf):
    """Yield unsigned LEB128 varints from a bytes buffer."""
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


def decode_posting_list(buf):
    """Decode one motif's posting block (count-prefixed, delta-encoded) to ids."""
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


class MotifIndex:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.motifs = []                 # list of motif rows (dicts)
        self.by_id = {}                  # motif_id -> row
        self.by_text = {}                # motif_text -> row
        self.proteins = []               # list of protein rows (dicts)
        self.protein_by_id = {}          # protein_id -> row
        self.protein_by_accession = {}   # accession -> protein_id
        self.build_info = {}
        self._postings_path = self.directory / "postings.bin"

    @classmethod
    def load(cls, directory):
        index = cls(directory)
        index._load_build_info()
        index._load_proteins()
        index._load_motifs()
        return index

    def _read_tsv(self, name):
        path = self.directory / name
        if not path.exists():
            return []
        with path.open(newline="", errors="replace") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def _load_build_info(self):
        for row in self._read_tsv("build_info.tsv"):
            self.build_info[row.get("key", "")] = row.get("value", "")

    def _load_motifs(self):
        for row in self._read_tsv("motifs.tsv"):
            row["motif_id"] = int(row.get("motif_id", 0))
            row["protein_count"] = int(row.get("protein_count", 0) or 0)
            row["posting_offset"] = int(row.get("posting_offset", 0) or 0)
            row["posting_bytes"] = int(row.get("posting_bytes", 0) or 0)
            self.motifs.append(row)
            self.by_id[row["motif_id"]] = row
            self.by_text[row.get("motif_text", "")] = row

    def _load_proteins(self):
        for row in self._read_tsv("proteins.tsv"):
            pid = int(row.get("protein_id", 0))
            row["protein_id"] = pid
            self.proteins.append(row)
            self.protein_by_id[pid] = row
            accession = row.get("accession", "")
            if accession:
                self.protein_by_accession[accession] = pid

    def protein_ids_for_motif(self, motif_id):
        row = self.by_id.get(motif_id)
        if row is None or row["posting_bytes"] == 0:
            return []
        with self._postings_path.open("rb") as f:
            f.seek(row["posting_offset"])
            buf = f.read(row["posting_bytes"])
        return decode_posting_list(buf)

    def accessions_for_motif(self, motif_id):
        out = []
        for pid in self.protein_ids_for_motif(motif_id):
            row = self.protein_by_id.get(pid)
            if row:
                out.append(row.get("accession", str(pid)))
        return out

    def motifs_for_accession(self, accession, limit=None):
        """Motifs whose posting list contains this protein. Linear scan -- the
        index is motif->proteins, so reverse lookups read each posting block."""
        pid = self.protein_by_accession.get(accession)
        if pid is None:
            return []
        hits = []
        for row in self.motifs:
            if pid in self.protein_ids_for_motif(row["motif_id"]):
                hits.append(row)
                if limit and len(hits) >= limit:
                    break
        return hits
