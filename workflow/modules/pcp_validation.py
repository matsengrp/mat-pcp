"""PCP validation functions for DASM-ready filtering."""

import pandas as pd

STOP_CODONS = ["TAA", "TAG", "TGA"]


def iter_codons(nt_seq):
    """Iterate over codons in a nucleotide sequence."""
    for i in range(0, (len(nt_seq) // 3) * 3, 3):
        yield nt_seq[i : i + 3]


def contains_stop_codon(sequence):
    """Check if sequence contains a stop codon."""
    sequence = sequence.upper()
    if len(sequence) % 3 != 0:
        raise ValueError("Sequence length must be divisible by 3")
    return any(codon in STOP_CODONS for codon in iter_codons(sequence))


def codon_is_valid(codon):
    """Return True if codon should be kept (no ambiguous bases)."""
    return "N" not in codon


def get_codon_mask(parent, child):
    """Return list of bools: True for codons to keep, False for masked codons."""
    return [
        codon_is_valid(p_codon) and codon_is_valid(c_codon)
        for p_codon, c_codon in zip(iter_codons(parent), iter_codons(child))
    ]


def apply_codon_mask(nt_seq, mask):
    """Return sequence with masked codons removed."""
    codons = list(iter_codons(nt_seq))
    return "".join(codon for codon, keep in zip(codons, mask) if keep)


def is_pcp_valid(parent, child):
    """
    Check if a parent-child pair is valid for DASM training.

    Requirements:
    1. Same length
    2. At least one unmasked codon (no 'N' in either sequence)
    3. Sequences differ after masking ambiguous codons
    """
    if len(parent) != len(child):
        return False

    mask = get_codon_mask(parent, child)

    # Must have at least one valid codon
    if not any(mask):
        return False

    # After masking, sequences must still differ
    masked_parent = apply_codon_mask(parent, mask)
    masked_child = apply_codon_mask(child, mask)

    return masked_parent != masked_child


def validate_pcp_df(pcp_df):
    """
    Validate a DataFrame of PCPs and return filtering results.

    Parameters:
        pcp_df: DataFrame with 'parent' and 'child' columns

    Returns:
        valid_mask: boolean Series (True = valid)
        invalid_reasons: dict mapping reason -> list of indices
    """
    invalid_reasons = {
        "identical": [],
        "parent_stop_codon": [],
        "child_stop_codon": [],
        "invalid_pcp": [],
    }

    valid_mask = pd.Series(True, index=pcp_df.index)

    for idx, row in pcp_df.iterrows():
        parent, child = row["parent"], row["child"]

        # Check identical
        if parent == child:
            invalid_reasons["identical"].append(idx)
            valid_mask[idx] = False
            continue

        # Check stop codons
        if contains_stop_codon(parent):
            invalid_reasons["parent_stop_codon"].append(idx)
            valid_mask[idx] = False
            continue

        if contains_stop_codon(child):
            invalid_reasons["child_stop_codon"].append(idx)
            valid_mask[idx] = False
            continue

        # Check PCP validity (masking, etc.)
        if not is_pcp_valid(parent, child):
            invalid_reasons["invalid_pcp"].append(idx)
            valid_mask[idx] = False

    return valid_mask, invalid_reasons
