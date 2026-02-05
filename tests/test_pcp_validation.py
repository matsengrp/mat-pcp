"""Unit tests for PCP validation module."""

import pytest
import pandas as pd
import sys
from pathlib import Path

# Add workflow modules to path
workflow_modules = Path(__file__).parent.parent / "workflow" / "modules"
sys.path.insert(0, str(workflow_modules))

from pcp_validation import validate_pcp_df


def test_identical_sequences_filtered():
    """PCPs with identical sequences should be filtered."""
    df = pd.DataFrame({
        "parent": ["ATGATG", "ATGATG"],
        "child": ["ATGATG", "ATGCTG"]  # First identical, second differs
    })
    valid_mask, reasons = validate_pcp_df(df)
    assert valid_mask[0] == False  # Identical filtered
    assert valid_mask[1] == True   # Different kept
    assert len(reasons["identical"]) == 1


def test_all_different_sequences_pass():
    """All different sequences should pass validation."""
    df = pd.DataFrame({
        "parent": ["ATGATG", "GCGCGC", "TTTAAA"],
        "child": ["ATGCTG", "GCGGGC", "TTTCAA"]
    })
    valid_mask, reasons = validate_pcp_df(df)
    assert all(valid_mask)
    assert len(reasons["identical"]) == 0


def test_stop_codons_allowed():
    """Stop codons should not cause filtering."""
    df = pd.DataFrame({
        "parent": ["ATGTAA", "ATGTAG"],  # Contains stop codons TAA, TAG
        "child": ["ATGCAA", "ATGCAG"]    # Different from parent
    })
    valid_mask, reasons = validate_pcp_df(df)
    assert valid_mask[0] == True  # Not filtered
    assert valid_mask[1] == True  # Not filtered
    assert len(reasons["identical"]) == 0


def test_ambiguous_bases_allowed():
    """Ambiguous bases should not cause filtering."""
    df = pd.DataFrame({
        "parent": ["ATGNAT", "NNNCAT"],  # Contains N
        "child": ["ATGCAT", "NNNAAT"]    # Different from parent
    })
    valid_mask, reasons = validate_pcp_df(df)
    assert valid_mask[0] == True  # Not filtered
    assert valid_mask[1] == True  # Not filtered
    assert len(reasons["identical"]) == 0


def test_empty_dataframe():
    """Empty dataframe should return empty results."""
    df = pd.DataFrame(columns=["parent", "child"])
    valid_mask, reasons = validate_pcp_df(df)
    assert len(valid_mask) == 0
    assert len(reasons["identical"]) == 0


def test_mixed_valid_invalid():
    """Mixed valid and invalid PCPs should be correctly classified."""
    df = pd.DataFrame({
        "parent": ["ATGATG", "GCGCGC", "TTTAAA", "CCCGGG"],
        "child": ["ATGATG", "GCGGGC", "TTTAAA", "CCCAAA"]
    })
    valid_mask, reasons = validate_pcp_df(df)
    assert valid_mask[0] == False  # Identical
    assert valid_mask[1] == True   # Different
    assert valid_mask[2] == False  # Identical
    assert valid_mask[3] == True   # Different
    assert len(reasons["identical"]) == 2


def test_sequences_already_sliced():
    """
    Test that validation works correctly when sequences are already sliced
    to gene region (no gene_region parameter needed).
    """
    # Simulate sequences already sliced to Spike region
    df = pd.DataFrame({
        "parent": ["ATGATGATG", "GCGCGCGCG"],
        "child": ["ATGATGATG", "GCGGGCGCG"]
    })
    valid_mask, reasons = validate_pcp_df(df)
    # First is identical in sliced region -> filtered
    # Second differs in sliced region -> kept
    assert valid_mask[0] == False
    assert valid_mask[1] == True
