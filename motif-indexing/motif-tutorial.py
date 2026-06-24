#!/usr/bin/env python3

from math import ceil
from random import Random

specificity = 0.35
min_proteins = 2
seed = 7

# ---------------------------------------------------------------------
# 1. Make raw sequence windows.
#
# These mimic the windows Rust has already collected under one endpoint
# group. The parent skeleton is not hard-coded yet; it will be derived
# from the shared first and last amino acids.
#
# For tutorial clarity, this uses one window per protein.
# Rust can see many windows per protein, then counts unique proteins.
# ---------------------------------------------------------------------

rng = Random(seed)

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
windows = []

for i in range(1, 21):
    seq = [
        "A",
        rng.choice(amino_acids),
        rng.choice(["P", "P", "P", "P", "K", "K", "S", "T"]),
        rng.choice(amino_acids),
        rng.choice(["L", "L", "L", "V", "V", "I"]),
        rng.choice(amino_acids),
        "G",
    ]

    windows.append((f"protein_{i:02d}", "".join(seq)))


# ---------------------------------------------------------------------
# An optional fixed example.
# ---------------------------------------------------------------------

# windows = [
#     ("protein_01", "AAPLLQG"),
#     ("protein_02", "AAPLLQG"),
#     ("protein_03", "AAPVVQG"),
#     ("protein_04", "AAKVVQG"),
#     ("protein_05", "AAKVVQG"),
#     ("protein_06", "AASVVQG"),
#     ("protein_07", "AASIIQG"),
#     ("protein_08", "AASIIQG"),
#     ("protein_09", "AATIIQG"),
#     ("protein_10", "AATIIQG"),
# ]


# ---------------------------------------------------------------------
# 2. Derive the starting parent skeleton from the sequence endpoints.
# ---------------------------------------------------------------------

k = len(windows[0][1])
root_skeleton = windows[0][1][0] + "." * (k - 2) + windows[0][1][-1]

internal_positions = k - 2
internal_budget = round(internal_positions * specificity)

print("RAW WINDOWS")
for protein, seq in windows:
    print(f"  {protein}: {seq}")

print()
print("START")
print(f"  parent skeleton: {root_skeleton}")
print(f"  k: {k}")
print(f"  internal positions: {internal_positions}")
print(f"  specificity: {specificity}")
print(f"  internal budget: {internal_budget}")
print()


# ---------------------------------------------------------------------
# 3. Recursive tutorial function.
#
# This is the same logical shape as Rust:
#
#   current skeleton owns a set of windows
#   test every still-open internal position
#   see what child skeletons would be born from that position
#   pick the best position
#   recurse into the children from that winning position
# ---------------------------------------------------------------------

def walk(parent_skeleton, node_windows, fixed_slots, depth=0):
    indent = "  " * depth

    parent_proteins = {protein for protein, seq in node_windows}
    parent_support = len(parent_proteins)

    fixed_internal = len(fixed_slots)

    print(f"{indent}NODE {parent_skeleton}")
    print(f"{indent}  parent support: {parent_support}")
    print(f"{indent}  fixed internal positions: {fixed_internal}/{internal_budget}")

    if fixed_internal >= internal_budget:
        print(f"{indent}  STOP: specificity budget is used")
        print()
        return

    # This is the dynamic part.
    #
    # Early in the recursion, child_min is low.
    # Later in the recursion, child_min gets stricter.
    next_fixed_internal = fixed_internal + 1
    progress = next_fixed_internal / internal_budget
    required_fraction = specificity * progress
    child_min = max(min_proteins, ceil(parent_support * required_fraction))

    print(
        f"{indent}  child_min: max({min_proteins}, "
        f"ceil({parent_support} * {required_fraction:.3f})) = {child_min}"
    )
    print()

    slot_results = []

    # This is the lynch pin:
    # each open slot is tested as the next possible position to fill.
    for slot in range(1, k - 1):
        if slot in fixed_slots:
            continue

        print(f"{indent}  TEST SLOT {slot}")

        residue_to_windows = {}

        for protein, seq in node_windows:
            residue = seq[slot]
            residue_to_windows.setdefault(residue, []).append((protein, seq))

        surviving_children = {}
        new_coverage = set()
        same_parent_count = 0

        for residue, child_windows in sorted(residue_to_windows.items()):
            child_proteins = {protein for protein, seq in child_windows}
            child_support = len(child_proteins)

            child_skeleton = list(parent_skeleton)
            child_skeleton[slot] = residue
            child_skeleton = "".join(child_skeleton)

            if child_support >= child_min:
                surviving_children[child_skeleton] = child_windows

                if child_proteins == parent_proteins:
                    same_parent_count += 1
                    kind = "same-as-parent"
                else:
                    new_coverage.update(child_proteins)
                    kind = "new-child"

                print(
                    f"{indent}    KEEP   {child_skeleton} "
                    f"support={child_support}  {kind}"
                )
            else:
                print(
                    f"{indent}    reject {child_skeleton} "
                    f"support={child_support}"
                )

        new_child_count = len(surviving_children) - same_parent_count
        total_surviving = len(surviving_children)

        # This matches the current Rust ranking idea:
        #
        #   1. more new child skeletons
        #   2. more total protein coverage among those new children
        #   3. more same-as-parent children
        #   4. more total surviving children
        #
        # Python tuple comparison, like Rust tuple/struct ordering here,
        # compares left to right.
        score = (
            new_child_count,
            len(new_coverage),
            same_parent_count,
            total_surviving,
        )

        print(f"{indent}    slot score = {score}")
        print()

        if total_surviving > 0:
            slot_results.append((score, slot, surviving_children))

    if not slot_results:
        print(f"{indent}  STOP: no slot can birth a child")
        print()
        return

    winning_score, winning_slot, born_children = max(slot_results, key=lambda x: x[0])

    print(f"{indent}  WINNING SLOT: {winning_slot}")
    print(f"{indent}  WINNING SCORE: {winning_score}")
    print(f"{indent}  CHILDREN BORN:")
    for child in born_children:
        print(f"{indent}    {child}")
    print()

    # Rust recurses into every child born from the winning slot.
    for child_skeleton, child_windows in born_children.items():
        next_fixed_slots = set(fixed_slots)
        next_fixed_slots.add(winning_slot)

        walk(
            child_skeleton,
            child_windows,
            next_fixed_slots,
            depth + 1,
        )


walk(root_skeleton, windows, fixed_slots=set())
