"""Integration of rG4 candidates with genomic APA annotations."""

from __future__ import annotations

import pandas as pd

from core.apa import standardize_apa_sites


MAPPING_COLUMNS = {
    "sequence_id", "chrom", "start", "end", "strand"
}


def integrate_rg4_apa(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    sequence_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """Map rG4 regions to genomic coordinates and test APA overlap."""
    mapping = sequence_mapping.copy()
    mapping.columns = [str(c).strip().lower() for c in mapping.columns]

    missing = MAPPING_COLUMNS - set(mapping.columns)
    if missing:
        raise ValueError(
            "Mapping table is missing: " + ", ".join(sorted(missing))
        )

    mapping = mapping.rename(columns={
        "start": "mapping_start",
        "end": "mapping_end",
        "strand": "mapping_strand",
        "chrom": "mapping_chrom",
    })

    mapping["sequence_id"] = mapping["sequence_id"].astype(str)
    mapping["mapping_start"] = pd.to_numeric(
        mapping["mapping_start"], errors="raise"
    )
    mapping["mapping_end"] = pd.to_numeric(
        mapping["mapping_end"], errors="raise"
    )

    if mapping["sequence_id"].duplicated().any():
        raise ValueError("Each sequence_id must have exactly one mapping.")

    if (mapping["mapping_start"] < 0).any():
        raise ValueError("Mapping start coordinates cannot be negative.")

    if (mapping["mapping_end"] <= mapping["mapping_start"]).any():
        raise ValueError("Mapping end must be greater than mapping start.")

    apa = standardize_apa_sites(apa_sites)

    columns = [
        "sequence_id", "candidate_start", "candidate_end",
        "chrom", "genomic_start", "genomic_end", "strand",
        "candidate_type", "score", "overlaps_apa",
        "overlapping_apa_sites",
    ]

    if rg4_candidates.empty:
        return pd.DataFrame(columns=columns)

    required = {
        "sequence_id", "start", "end", "candidate_type", "score"
    }
    if not required.issubset(rg4_candidates.columns):
        raise ValueError("rG4 candidate table is missing required columns.")

    merged = rg4_candidates.merge(
        mapping,
        on="sequence_id",
        how="left",
        validate="many_to_one",
    )

    results = []

    for _, candidate in merged.iterrows():
        if pd.isna(candidate["mapping_chrom"]):
            continue

        local_start = int(candidate["start"])
        local_end = int(candidate["end"])
        mapping_start = int(candidate["mapping_start"])
        mapping_end = int(candidate["mapping_end"])

        if local_start < 0 or local_end <= local_start:
            continue

        if local_end > mapping_end - mapping_start:
            continue

        chrom = str(candidate["mapping_chrom"])
        strand = str(candidate["mapping_strand"])
        genomic_start = mapping_start + local_start
        genomic_end = mapping_start + local_end

        overlapping = apa[
            (apa["chrom"] == chrom)
            & (apa["start"] < genomic_end)
            & (apa["end"] > genomic_start)
        ]

        overlapping = overlapping[
            (overlapping["strand"] == ".")
            | (overlapping["strand"] == strand)
        ]

        results.append({
            "sequence_id": candidate["sequence_id"],
            "candidate_start": local_start,
            "candidate_end": local_end,
            "chrom": chrom,
            "genomic_start": genomic_start,
            "genomic_end": genomic_end,
            "strand": strand,
            "candidate_type": candidate["candidate_type"],
            "score": candidate["score"],
            "overlaps_apa": not overlapping.empty,
            "overlapping_apa_sites": int(len(overlapping)),
        })

    return pd.DataFrame(results, columns=columns)
