
"""Genomic integration of RNA G4 candidates and APA annotations.

Coordinates:
- Sequence candidate coordinates: 0-based, half-open.
- Mapping start/end: 0-based, half-open genomic coordinates.
- APA start/end: expected to be 0-based, half-open.

All inputs must use the same genome assembly and compatible
transcript-to-genome mapping conventions.
"""

from __future__ import annotations

import pandas as pd


MAPPING_REQUIRED = {
    "sequence_id",
    "chrom",
    "start",
    "end",
    "strand",
}

RG4_REQUIRED = {
    "sequence_id",
    "start",
    "end",
    "candidate_type",
}

APA_REQUIRED = {
    "sequence_id",
    "chrom",
    "start",
    "end",
    "strand",
    "site_type",
}


def _validate_columns(
    frame: pd.DataFrame,
    required: set[str],
    label: str,
) -> None:
    missing = required - set(frame.columns)

    if missing:
        raise ValueError(
            f"{label} is missing required columns: "
            + ", ".join(sorted(missing))
        )


def _prepare_coordinates(
    frame: pd.DataFrame,
    label: str,
) -> pd.DataFrame:
    """Convert coordinate fields to numeric values and validate them."""
    result = frame.copy()

    for column in ("start", "end"):
        result[column] = pd.to_numeric(
            result[column],
            errors="coerce",
        )

    if result[["start", "end"]].isna().any().any():
        raise ValueError(
            f"{label} contains non-numeric or missing coordinates."
        )

    if (result["start"] < 0).any():
        raise ValueError(
            f"{label} contains negative start coordinates."
        )

    if (result["end"] <= result["start"]).any():
        raise ValueError(
            f"{label} contains invalid intervals. "
            "End must be greater than start."
        )

    result["start"] = result["start"].astype(int)
    result["end"] = result["end"].astype(int)

    return result


def integrate_rg4_apa(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    sequence_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """
    Map RNA candidate intervals to genomic coordinates and
    identify overlaps with APA annotations.

    The mapping table must contain one genomic interval per
    input sequence. Candidate local coordinates are added to
    the mapping start coordinate.

    Overlap is calculated using half-open intervals:
        candidate_start < apa_end
        candidate_end > apa_start

    Chromosome, strand and sequence identifier must match.
    """
    _validate_columns(
        rg4_candidates,
        RG4_REQUIRED,
        "rG4 candidate table",
    )

    _validate_columns(
        apa_sites,
        APA_REQUIRED,
        "APA annotation table",
    )

    _validate_columns(
        sequence_mapping,
        MAPPING_REQUIRED,
        "Sequence mapping table",
    )

    if rg4_candidates.empty:
        return pd.DataFrame(columns=[
            "sequence_id",
            "candidate_start",
            "candidate_end",
            "chrom",
            "genomic_start",
            "genomic_end",
            "strand",
            "candidate_type",
            "score",
            "overlaps_apa",
            "overlapping_apa_sites",
            "overlapping_apa_types",
            "overlapping_apa_coordinates",
        ])

    candidates = rg4_candidates.copy()
    mappings = sequence_mapping.copy()
    apa = apa_sites.copy()

    # Validate candidate coordinates.
    candidates["start"] = pd.to_numeric(
        candidates["start"],
        errors="coerce",
    )
    candidates["end"] = pd.to_numeric(
        candidates["end"],
        errors="coerce",
    )

    if candidates[["start", "end"]].isna().any().any():
        raise ValueError(
            "rG4 candidates contain invalid coordinates."
        )

    if (
        (candidates["start"] < 0)
        | (candidates["end"] <= candidates["start"])
    ).any():
        raise ValueError(
            "rG4 candidates must have valid 0-based half-open intervals."
        )

    candidates["start"] = candidates["start"].astype(int)
    candidates["end"] = candidates["end"].astype(int)

    mappings = _prepare_coordinates(
        mappings,
        "Sequence mapping table",
    )

    apa = _prepare_coordinates(
        apa,
        "APA annotation table",
    )

    # Normalize chromosome and strand labels.
    for frame in (mappings, apa):
        frame["chrom"] = frame["chrom"].astype(str).str.strip()
        frame["strand"] = frame["strand"].astype(str).str.strip()

    if not mappings["strand"].isin(["+", "-"]).all():
        raise ValueError(
            "Mapping strand must be '+' or '-'."
        )

    if not apa["strand"].isin(["+", "-"]).all():
        raise ValueError(
            "APA strand must be '+' or '-'."
        )

    if mappings["sequence_id"].duplicated().any():
        raise ValueError(
            "Mapping table must contain one row per sequence_id."
        )

    mapping_lookup = mappings.set_index("sequence_id").to_dict(
        orient="index"
    )

    results = []

    for _, candidate in candidates.iterrows():
        sequence_id = candidate["sequence_id"]

        if sequence_id not in mapping_lookup:
            raise ValueError(
                f"No genomic mapping found for sequence_id "
                f"'{sequence_id}'."
            )

        mapping = mapping_lookup[sequence_id]

        local_start = int(candidate["start"])
        local_end = int(candidate["end"])

        mapping_start = int(mapping["start"])
        mapping_end = int(mapping["end"])

        sequence_length = mapping_end - mapping_start

        if local_end > sequence_length:
            raise ValueError(
                f"Candidate {sequence_id}:{local_start}-{local_end} "
                "extends beyond its mapped sequence interval."
            )

        genomic_start = mapping_start + local_start
        genomic_end = mapping_start + local_end

        chromosome = str(mapping["chrom"])
        strand = str(mapping["strand"])

        matching_apa = apa[
            (apa["sequence_id"] == sequence_id)
            & (apa["chrom"] == chromosome)
            & (apa["strand"] == strand)
            & (apa["start"] < genomic_end)
            & (apa["end"] > genomic_start)
        ].copy()

        overlapping_types = sorted(
            matching_apa["site_type"]
            .astype(str)
            .unique()
            .tolist()
        )

        overlapping_coordinates = [
            f"{row.chrom}:{row.start}-{row.end} ({row.site_type})"
            for row in matching_apa.itertuples()
        ]

        score = candidate.get("score", None)

        results.append({
            "sequence_id": sequence_id,
            "candidate_start": local_start,
            "candidate_end": local_end,
            "chrom": chromosome,
            "genomic_start": genomic_start,
            "genomic_end": genomic_end,
            "strand": strand,
            "candidate_type": candidate["candidate_type"],
            "score": score,
            "overlaps_apa": not matching_apa.empty,
            "overlapping_apa_sites": len(matching_apa),
            "overlapping_apa_types": ", ".join(overlapping_types),
            "overlapping_apa_coordinates": "; ".join(
                overlapping_coordinates
            ),
        })

    return pd.DataFrame(results).sort_values(
        ["sequence_id", "genomic_start", "genomic_end"]
    ).reset_index(drop=True)
