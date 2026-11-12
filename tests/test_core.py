"""Basic tests for RNA-Insight core modules."""

import pandas as pd

from core.io import validate_rna_records
from core.rg4 import scan_rg4, composition_window_score
from core.apa import standardize_apa_sites, summarize_apa_sites


def test_rna_validation():
    records = pd.DataFrame({
        "sequence_id": ["test1", "test2"],
        "sequence": ["AUGCUU", "AUGCXU"],
    })

    result = validate_rna_records(records)

    assert result.loc[0, "valid_sequence"]
    assert not result.loc[1, "valid_sequence"]


def test_composition_score():
    assert composition_window_score("GGGG") == 1.0
    assert composition_window_score("CCCC") == -1.0


def test_rg4_canonical_candidate():
    records = pd.DataFrame({
        "sequence_id": ["test1"],
        "sequence": ["GGGAAAGGGAAAGGGAAAGGG"],
    })

    result = scan_rg4(records)

    assert not result.empty
    assert "canonical_motif" in set(result["candidate_type"])


def test_apa_standardization():
    frame = pd.DataFrame({
        "sequence_id": ["test1"],
        "chrom": ["chr1"],
        "start": [100],
        "end": [101],
        "strand": ["+"],
        "site_type": ["proximal"],
    })

    result = standardize_apa_sites(frame)
    summary = summarize_apa_sites(result)

    assert summary["total_sites"] == 1
    assert summary["proximal_sites"] == 1
