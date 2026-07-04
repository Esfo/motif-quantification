"""Elemental and amino-acid chemistry constants for the viewer.

This is the single source of truth for masses, isotope abundances, amino-acid
compositions, and the derived lookup tables the isotope-distribution math in
``isotopes.py`` needs. The numbers are NIST values and match
``distributions/elementalcomponents.py`` in the indexing pipeline; they live
here so the viewer never has to import across the pipeline package.

Derived tables (built once at import) mirror the structures the original
``descending_partial_products`` exploration relied on, so the ported isotope
walk in ``isotopes.py`` can stay faithful to the validated implementation.
"""

PROTON = 1.007276554940804

# element: {isotope_label: (mass, natural_abundance)}, isotopes in ascending mass
ELEMENT_INFO = {
    "H": {"H1": (1.00782503223, 0.999885), "H2": (2.01410177812, 0.000115)},
    "C": {"C12": (12.0, 0.9893), "C13": (13.00335483507, 0.0107)},
    "N": {"N14": (14.00307400443, 0.99636), "N15": (15.00010889888, 0.00364)},
    "O": {
        "O16": (15.99491461957, 0.99757),
        "O17": (16.99913175650, 0.00038),
        "O18": (17.99915961286, 0.00205),
    },
    "S": {
        "S32": (31.9720711744, 0.9499),
        "S33": (32.9714589098, 0.0075),
        "S34": (33.967867004, 0.0425),
        "S36": (35.96708071, 0.0001),
    },
    "P": {"P31": (30.97376199842, 1.0)},
    "Se": {
        "Se74": (73.922475934, 0.0089),
        "Se76": (75.919213704, 0.0937),
        "Se77": (76.919914154, 0.0763),
        "Se78": (77.91730928, 0.2377),
        "Se80": (79.9165218, 0.4961),
        "Se82": (81.9166995, 0.0873),
    },
}

# residue: {element: count} for the +H2O-free residue mass (water added per peptide)
AMINO_ACID_COMPOSITION = {
    "A": {"C": 3, "H": 5, "N": 1, "O": 1},
    "R": {"C": 6, "H": 12, "N": 4, "O": 1},
    "N": {"C": 4, "H": 6, "N": 2, "O": 2},
    "D": {"C": 4, "H": 5, "N": 1, "O": 3},
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},
    "Q": {"C": 5, "H": 8, "N": 2, "O": 2},
    "E": {"C": 5, "H": 7, "N": 1, "O": 3},
    "G": {"C": 2, "H": 3, "N": 1, "O": 1},
    "H": {"C": 6, "H": 7, "N": 3, "O": 1},
    "I": {"C": 6, "H": 11, "N": 1, "O": 1},
    "L": {"C": 6, "H": 11, "N": 1, "O": 1},
    "K": {"C": 6, "H": 12, "N": 2, "O": 1},
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},
    "F": {"C": 9, "H": 9, "N": 1, "O": 1},
    "P": {"C": 5, "H": 7, "N": 1, "O": 1},
    "S": {"C": 3, "H": 5, "N": 1, "O": 2},
    "T": {"C": 4, "H": 7, "N": 1, "O": 2},
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2},
    "V": {"C": 5, "H": 9, "N": 1, "O": 1},
    "U": {"C": 3, "H": 5, "N": 1, "O": 1, "Se": 1},  # selenocysteine
    "O": {"C": 12, "H": 19, "N": 3, "O": 2},  # pyrrolysine
    # J = Leu/Ile ambiguity; both have the SAME atomic composition, so this is
    # exact (some search engines emit J for the indistinguishable pair).
    "J": {"C": 6, "H": 11, "N": 1, "O": 1},
}

# n-/c-terminal fragment composition deltas relative to the cumulative residue sum
N_FRAGMENT_COMPOSITIONS = {
    "a": {"C": -1, "O": -1},
    "b": {},
    "c": {"N": 1, "H": 3},
}
C_FRAGMENT_COMPOSITIONS = {
    "x": {"C": 1, "O": 2},
    "y": {"H": 2, "O": 1},
    "z": {"N": -1, "H": -1, "O": 1},
}


# ---- derived tables (built once) -------------------------------------------

ELEMENTAL_MASSES = {}        # isotope label -> mass
ELEMENTAL_PROBABILITIES = {}  # isotope label -> abundance
ISOTOPES_BY_ELEMENT = {}     # element -> (isotope labels, ascending mass)
MONOISOTOPIC_KEYS = {}       # element -> lowest-mass isotope label
NONMONOISOTOPIC_GROUPS = {}  # element -> tuple of the remaining isotope labels
NONMONOISOTOPIC_ELEMENTS = set()  # all non-monoisotopic isotope labels
ELEMENT_VECTORS = {}         # element -> [0, 0, ...] one slot per isotope
VECTOR_POSITIONS = {}        # element -> {isotope label: index}
ELEMENT_POSITIONS = {}       # element -> {index: isotope label}

for _element, _isos in ELEMENT_INFO.items():
    _labels = list(_isos.keys())
    ISOTOPES_BY_ELEMENT[_element] = tuple(_labels)
    MONOISOTOPIC_KEYS[_element] = _labels[0]
    NONMONOISOTOPIC_GROUPS[_element] = tuple(_labels[1:])
    ELEMENT_VECTORS[_element] = [0 for _ in _labels]
    VECTOR_POSITIONS[_element] = {label: i for i, label in enumerate(_labels)}
    ELEMENT_POSITIONS[_element] = {i: label for i, label in enumerate(_labels)}

    for _i, (_label, (_mass, _abundance)) in enumerate(_isos.items()):
        ELEMENTAL_MASSES[_label] = _mass
        ELEMENTAL_PROBABILITIES[_label] = _abundance
        if _i > 0:
            NONMONOISOTOPIC_ELEMENTS.add(_label)


def peptide_atomic_composition(sequence):
    """Atomic composition of a peptide (residue sums + one water)."""
    from collections import Counter

    composition = Counter()
    for residue in sequence:
        try:
            composition += Counter(AMINO_ACID_COMPOSITION[residue])
        except KeyError:
            raise ValueError(
                f"unknown residue {residue!r} in peptide {sequence!r} "
                f"(no atomic composition)") from None

    composition["H"] += 2
    composition["O"] += 1
    return composition


def monoisotopic_mass(composition):
    """Neutral monoisotopic mass for an element->count composition."""
    return sum(ELEMENTAL_MASSES[MONOISOTOPIC_KEYS[element]] * count
               for element, count in composition.items())
