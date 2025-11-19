#!/usr/bin/env python3
import argparse, sys, os
from typing import List, Tuple, Optional
# Convert RNA FASTA + dot-bracket secondary structure into:
#   (1) QRNAS-style BASEPAIR restraints
#   (2) Phenix nucleic_acid secondary-structure restraint blocks (.eff)
#
# This script is designed for RNA structural biology workflows, where base-pair
# restraints can be explicitly provided. It extracts base pairs from
# dot-bracket notation (including pseudoknots using (), [], {}, <>), aligns
# them to the FASTA sequence, and assigns Saenger classes assuming all base
# pairs are cis Watson–Crick / Watson–Crick:
#
#       G–C → class 19
#       A–U (A–T) → class 20
#       G–U (G–T) → class 28
#
#
# INPUT (positional):
#   1) FASTA file containing the RNA sequence
#   2) Dot-bracket file describing secondary structure
#
# KEY OPTIONS:
#   --chain <ID>       Set PDB chain ID (default: A)
#   --start <N>        First residue number (default: 1)
#   --format <mode>    qrnas | phenix | both (default: phenix)
#   --saenger <N>      Override all inferred Saenger classes
#
# OUTPUT:
#   - Printed QRNAS "BASEPAIR" lines (if requested)
#   - Printed Phenix .eff block defining base_pair {...} restraints
#
# LIMITATIONS:
#   - Only cis WC–WC classes (19/20/28) are inferred automatically.
#   - Does not annotate trans, Hoogsteen, or noncanonical interactions.
#   - FASTA and dot-bracket lengths must match; script warns otherwise.
#
# Example:
#   python dot_bracket_to_eff.py 6hbc.fasta struct.db \
#       --chain B --format phenix > SS.eff
#
# Author: Ewan K. S. McRae (2025)
# Center for RNA Therapeutics, Houston Methodist Research Institute
# -----------------------------------------------------------------------------
# ---------------------------
# Dot-bracket parsing
# ---------------------------

def parse_dotbracket(dotbracket: str, start_resnum: int = 1) -> List[Tuple[int, int]]:
    """
    Convert Vienna dot-bracket into list of (i, j) pairs.
    Supports (), [], {}, <> for pseudoknots.
    """
    pairs = []
    stacks = {'(': [], '[': [], '{': [], '<': []}
    match = {')': '(', ']': '[', '}': '{', '>': '<'}

    s = "".join(c for c in dotbracket.strip() if c in ".()[]{}<>")

    for idx, c in enumerate(s, start=start_resnum):
        if c in stacks:
            stacks[c].append(idx)
        elif c in match:
            left = match[c]
            if not stacks[left]:
                raise ValueError(f"Unmatched closing '{c}' at {idx}")
            j = stacks[left].pop()
            pairs.append((j, idx))
        elif c != '.':
            raise ValueError(f"Unexpected character '{c}' at {idx}")

    for k, st in stacks.items():
        if st:
            raise ValueError(f"Unmatched opening '{k}' at {st}")

    pairs.sort(key=lambda x:(x[0], x[1]))
    return pairs

# ---------------------------
# FASTA loading + Saenger inference
# ---------------------------

def load_fasta_file(path: str) -> str:
    """Load FASTA, returning uppercase A/C/G/U/T string."""
    with open(path) as fh:
        seq_chars = []
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            for c in line:
                cu = c.upper()
                if cu in ("A","C","G","U","T"):
                    seq_chars.append(cu)
        return "".join(seq_chars)

def count_dot_length(dotbracket: str) -> int:
    return len("".join(c for c in dotbracket.strip() if c in ".()[]{}<>"))

def infer_saenger_pair(b1: str, b2: str) -> Optional[int]:
    """Assign Saenger class for cis Watson–Crick base pairs."""
    b1 = b1.upper().replace("T","U")
    b2 = b2.upper().replace("T","U")
    pair = {b1, b2}

    if pair == {"G","C"}: return 19
    if pair == {"A","U"}: return 20
    if pair == {"G","U"}: return 28
    return None

def infer_saenger_for_pairs(pairs, seq, start_resnum):
    classes = []
    n = len(seq)
    for (i,j) in pairs:
        i0 = i - start_resnum
        j0 = j - start_resnum
        if not (0 <= i0 < n and 0 <= j0 < n):
            print(f"Warning: residues {i}, {j} out of range for sequence length {n}", file=sys.stderr)
            classes.append(None)
            continue

        sclass = infer_saenger_pair(seq[i0], seq[j0])
        if sclass is None:
            print(f"Warning: cannot infer Saenger class for {seq[i0]}-{seq[j0]}", file=sys.stderr)
        classes.append(sclass)
    return classes

# ---------------------------
# Output formatting
# ---------------------------

def format_qrnas(pairs, chain):
    return [f"BASEPAIR   {chain}/{i}   {chain}/{j}" for i,j in pairs]

def format_phenix(pairs, chain, saenger_list, override=None):
    out = []
    out.append("pdb_interpretation {")
    out.append("  secondary_structure {")
    out.append("    nucleic_acid {")

    for idx, (i,j) in enumerate(pairs):
        out.append("      base_pair {")
        out.append(f"        base1 = chain '{chain}' and resid {i}")
        out.append(f"        base2 = chain '{chain}' and resid {j}")

        sclass = saenger_list[idx]
        if override is not None:
            sclass = override
        if sclass is not None:
            out.append(f"        saenger_class = {sclass}")

        out.append("      }")

    out.append("    }")
    out.append("  }")
    out.append("}")
    return "\n".join(out)

# ---------------------------
# Main
# ---------------------------

if __name__ == "__main__":

    ap = argparse.ArgumentParser(
        description="Generate QRNAS and/or Phenix restraints from FASTA + dot-bracket"
    )

    ap.add_argument("fasta", help="FASTA file")
    ap.add_argument("secstruct", help="Dot-bracket secondary structure file")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--format", choices=["qrnas","phenix","both"], default="phenix")
    ap.add_argument("--saenger", help="Override Saenger class for all pairs (e.g. 19)")

    args = ap.parse_args()

    # Load FASTA
    if not os.path.exists(args.fasta):
        sys.exit(f"Error: FASTA file '{args.fasta}' not found.")
    seq = load_fasta_file(args.fasta)

    # Load dot-bracket
    if not os.path.exists(args.secstruct):
        sys.exit(f"Error: structure file '{args.secstruct}' not found.")
    dot = open(args.secstruct).read()

    # Parse
    pairs = parse_dotbracket(dot, start_resnum=args.start)

    # Infer Saenger classes
    saenger_list = infer_saenger_for_pairs(pairs, seq, args.start)

    # Output
    if args.format in ("qrnas","both"):
        for line in format_qrnas(pairs, args.chain):
            print(line)

    if args.format in ("phenix","both"):
        if args.format == "both":
            print()
        print(format_phenix(
            pairs,
            args.chain,
            saenger_list,
            override=args.saenger
        ))
