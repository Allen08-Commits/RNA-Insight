
"""Exon-aware genomic integration of RNA G4 candidates and APA sites.

Coordinate conventions:
- rG4 candidate coordinates: 0-based, half-open transcript coordinates.
- Exon genomic coordinates: 0-based, half-open.
- Exon transcript coordinates: 0-based, half-open.
- APA genomic coordinates: 0-based, half-open.

The mapping table must contain one row per exon, not one row per
whole transcript.

Required mapping columns:
sequence_id, chrom, start, end, strand, transcript_start,
transcript_end

The mapping and APA annotation files must use the same genome
assembly and chromosome naming convention.
"""

from __future__ import annotations

import pandas as pd


MAPPING_REQUIRED = {
    "sequence_id",
    "chrom",
    "start",
    "end",
    "strand",
    "transcript_start",
    "transcript_end",
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


def _validate_columns(frame, required, label):
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"{label} is missing required columns: "
            + ", ".join(sorted(missing))
        )


def _numeric_coordinates(frame, columns, label):
    result = frame.copy()

    for column in columns:
        result[column] = pd.to_numeric(
            result[column], errors="coerce"
        )

    if result[columns].isna().any().any():
        raise ValueError(
            f"{label} contains missing or non-numeric coordinates."
        )

    for column in columns:
        if (result[column] < 0).any():
            raise ValueError(
                f"{label} contains negative coordinates."
            )

        if (result[column] % 1 != 0).any():
            raise ValueError(
                f"{label} coordinates must be integers."
            )

        result[column] = result[column].astype(int)

    return result


def _prepare_mapping(mapping):
    mapping = _numeric_coordinates(
        mapping,
        ["start", "end", "transcript_start", "transcript_end"],
        "Exon mapping",
    )

    if (mapping["end"] <= mapping["start"]).any():
        raise ValueError("Exon genomic end must exceed start.")

    if (
        mapping["transcript_end"] <= mapping["transcript_start"]
    ).any():
        raise ValueError("Exon transcript end must exceed start.")

    mapping["sequence_id"] = mapping["sequence_id"].astype(str).str.strip()
    mapping["chrom"] = mapping["chrom"].astype(str).str.strip()
    mapping["strand"] = mapping["strand"].astype(str).str.strip()

    if not mapping["strand"].isin(["+", "-"]).all():
        raise ValueError("Mapping strand must be '+' or '-'.")

    genomic_lengths = mapping["end"] - mapping["start"]
    transcript_lengths = (
        mapping["transcript_end"] - mapping["transcript_start"]
    )

    if not (genomic_lengths == transcript_lengths).all():
        raise ValueError(
            "Each exon must have equal genomic and transcript lengths."
        )

    mapping = mapping.sort_values(
        ["sequence_id", "transcript_start"]
    ).reset_index(drop=True)

    for sequence_id, group in mapping.groupby("sequence_id"):
        if group["chrom"].nunique() != 1:
            raise ValueError(
                f"Transcript {sequence_id} maps to multiple chromosomes."
            )

        if group["strand"].nunique() != 1:
            raise ValueError(
                f"Transcript {sequence_id} has inconsistent strands."
            )

        if group["transcript_start"].iloc[0] != 0:
            raise ValueError(
                f"Transcript {sequence_id} mapping does not start at 0."
            )

        if (
            group["transcript_end"].iloc[:-1].to_numpy()
            != group["transcript_start"].iloc[1:].to_numpy()
        ).any():
            raise ValueError(
                f"Transcript {sequence_id} has gaps or overlaps "
                "in transcript-coordinate exon mapping."
            )

        if group["strand"].iloc[0] == "+":
            genomic_order = group.sort_values("start")
        else:
            genomic_order = group.sort_values(
                "start", ascending=False
            )

        genomic_starts = genomic_order["start"].to_numpy()
        genomic_ends = genomic_order["end"].to_numpy()

        if len(genomic_starts) > 1:
            if (
                group["strand"].iloc[0] == "+"
                and (genomic_ends[:-1] > genomic_starts[1:]).any()
            ):
                raise ValueError(
                    f"Overlapping genomic exons found for {sequence_id}."
                )

            if (
                group["strand"].iloc[0] == "-"
                and (genomic_starts[:-1] < genomic_ends[1:]).any()
            ):
                raise ValueError(
                    f"Overlapping genomic exons found for {sequence_id}."
                )

    return mapping


def _map_candidate_to_exons(candidate, exon_rows):
    """Map a transcript candidate to one or more genomic exon segments."""
    local_start = int(candidate["start"])
    local_end = int(candidate["end"])

    transcript_length = int(
        exon_rows["transcript_end"].max()
    )

    if local_end > transcript_length:
        raise ValueError(
            f"Candidate {candidate['sequence_id']}:"
            f"{local_start}-{local_end} exceeds transcript length "
            f"{transcript_length}."
        )

    segments = []
    covered_bases = 0

    for exon in exon_rows.itertuples(index=False):
        overlap_start = max(
            local_start, int(exon.transcript_start)
        )
        overlap_end = min(
            local_end, int(exon.transcript_end)
        )

        if overlap_start >= overlap_end:
            continue

        covered_bases += overlap_end - overlap_start

        if exon.strand == "+":
            genomic_start = (
                int(exon.start)
                + overlap_start
                - int(exon.transcript_start)
            )
            genomic_end = (
                int(exon.start)
                + overlap_end
                - int(exon.transcript_start)
            )
        else:
            genomic_start = (
                int(exon.end)
                - (overlap_end - int(exon.transcript_start))
            )
            genomic_end = (
                int(exon.end)
                - (overlap_start - int(exon.transcript_start))
            )

        segments.append({
            "chrom": str(exon.chrom),
            "start": genomic_start,
            "end": genomic_end,
            "strand": str(exon.strand),
            "transcript_start": overlap_start,
            "transcript_end": overlap_end,
        })

    if covered_bases != local_end - local_start:
        raise ValueError(
            f"Candidate {candidate['sequence_id']}:"
            f"{local_start}-{local_end} is not completely covered "
            "by the supplied exon mapping."
        )

    return segments


def integrate_rg4_apa(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    sequence_mapping: pd.DataFrame,
) -> pd.DataFrame:
    """
    Map rG4 candidates through exon coordinates and detect APA overlaps.

    APA sequence_id must identify the same transcript as the
    rG4 candidate and exon mapping.

    Candidates spanning exon junctions are represented by multiple
    genomic segments. An overlap is reported if any segment overlaps
    an APA interval on the same chromosome and strand.
    """
    _validate_columns(
        rg4_candidates, RG4_REQUIRED, "rG4 candidate table"
    )
    _validate_columns(
        apa_sites, APA_REQUIRED, "APA annotation table"
    )
    _validate_columns(
        sequence_mapping, MAPPING_REQUIRED, "Exon mapping table"
    )

    output_columns = [
        "sequence_id",
        "candidate_start",
        "candidate_end",
        "chrom",
        "genomic_segments",
        "strand",
        "candidate_type",
        "score",
        "overlaps_apa",
        "overlapping_apa_sites",
        "overlapping_apa_types",
        "overlapping_apa_coordinates",
    ]

    if rg4_candidates.empty:
        return pd.DataFrame(columns=output_columns)

    candidates = _numeric_coordinates(
        rg4_candidates,
        ["start", "end"],
        "rG4 candidates",
    )

    if (candidates["end"] <= candidates["start"]).any():
        raise ValueError("rG4 candidate intervals are invalid.")

    candidates["sequence_id"] = (
        candidates["sequence_id"].astype(str).str.strip()
    )

    mapping = _prepare_mapping(sequence_mapping)

    apa = _numeric_coordinates(
        apa_sites,
        ["start", "end"],
        "APA annotation",
    )

    if (apa["end"] <= apa["start"]).any():
        raise ValueError("APA intervals are invalid.")

    apa["sequence_id"] = apa["sequence_id"].astype(str).str.strip()
    apa["chrom"] = apa["chrom"].astype(str).str.strip()
    apa["strand"] = apa["strand"].astype(str).str.strip()
    apa["site_type"] = apa["site_type"].fillna(
        "unspecified"
    ).astype(str)

    if not apa["strand"].isin(["+", "-"]).all():
        raise ValueError("APA strand must be '+' or '-'.")

    mapping_groups = {
        sequence_id: group.copy()
        for sequence_id, group in mapping.groupby("sequence_id")
    }

    apa_groups = {
        sequence_id: group.copy()
        for sequence_id, group in apa.groupby("sequence_id")
    }

    results = []

    for candidate in candidates.itertuples(index=False):
        sequence_id = str(candidate.sequence_id)

        if sequence_id not in mapping_groups:
            raise ValueError(
                f"No exon mapping found for transcript {sequence_id}."
            )

        exon_rows = mapping_groups[sequence_id]

        segments = _map_candidate_to_exons(
            {
                "sequence_id": sequence_id,
                "start": candidate.start,
                "end": candidate.end,
            },
            exon_rows,
        )

        transcript_chrom = str(exon_rows["chrom"].iloc[0])
        transcript_strand = str(exon_rows["strand"].iloc[0])

        transcript_apa = apa_groups.get(
            sequence_id,
            pd.DataFrame(columns=apa.columns),
        )

        matching_rows = []

        for segment in segments:
            matches = transcript_apa[
                (transcript_apa["chrom"] == segment["chrom"])
                & (transcript_apa["strand"] == segment["strand"])
                & (transcript_apa["start"] < segment["end"])
                & (transcript_apa["end"] > segment["start"])
            ]

            matching_rows.append(matches)

        if matching_rows:
            matched = pd.concat(
                matching_rows, ignore_index=True
            ).drop_duplicates()
        else:
            matched = pd.DataFrame(columns=apa.columns)

        segment_labels = [
            f"{segment['chrom']}:{segment['start']}-"
            f"{segment['end']}({segment['strand']})"
            for segment in segments
        ]

        apa_labels = [
            f"{row.chrom}:{row.start}-{row.end} ({row.site_type})"
            for row in matched.itertuples()
        ]

        apa_types = sorted(
            matched["site_type"].astype(str).unique().tolist()
        ) if not matched.empty else []

        results.append({
            "sequence_id": sequence_id,
            "candidate_start": int(candidate.start),
            "candidate_end": int(candidate.end),
            "chrom": transcript_chrom,
            "genomic_segments": "; ".join(segment_labels),
            "strand": transcript_strand,
            "candidate_type": candidate.candidate_type,
            "score": getattr(candidate, "score", None),
            "overlaps_apa": not matched.empty,
            "overlapping_apa_sites": int(len(matched)),
            "overlapping_apa_types": ", ".join(apa_types),
            "overlapping_apa_coordinates": "; ".join(apa_labels),
        })

    return pd.DataFrame(results, columns=output_columns).sort_values(
        ["sequence_id", "candidate_start", "candidate_end"]
    ).reset_index(drop=True)
