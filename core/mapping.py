
"""Exon-aware transcript-to-genome mapping using Ensembl REST API."""

from __future__ import annotations

import time
from typing import Any

import pandas as pd
import requests


ENSEMBL_BASE = "https://rest.ensembl.org"
HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
    "User-Agent": "RNA-Insight/1.0",
}

REQUIRED_TRANSCRIPT_COLUMNS = {
    "transcript_id",
    "sequence_length",
}


class MappingError(Exception):
    """Raised when Ensembl transcript mapping cannot be completed."""


def _get_transcript_annotation(
    transcript_id: str,
    timeout: int = 30,
    retries: int = 3,
) -> dict[str, Any]:
    """Retrieve an expanded Ensembl transcript annotation."""

    url = f"{ENSEMBL_BASE}/lookup/id/{transcript_id}"

    last_error = None

    for attempt in range(retries):

        try:
            response = requests.get(
                url,
                headers=HEADERS,
                params={"expand": 1},
                timeout=timeout,
            )

            if response.status_code == 429 or response.status_code >= 500:
                time.sleep(1.5 * (attempt + 1))
                continue

            response.raise_for_status()

            data = response.json()

            if not isinstance(data, dict):
                raise MappingError(
                    f"Unexpected Ensembl response for {transcript_id}"
                )

            return data

        except requests.RequestException as exc:
            last_error = exc
            time.sleep(1.5 * (attempt + 1))

    raise MappingError(
        f"Could not retrieve {transcript_id}: {last_error}"
    )


def _get_value(row: Any, column: str, default=None):
    """Safely obtain a value from a pandas row."""

    if column not in row.index:
        return default

    value = row[column]

    if pd.isna(value):
        return default

    return value


def build_transcript_exon_mapping(
    transcripts: pd.DataFrame,
    delay_seconds: float = 0.15,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Construct a 0-based, half-open exon mapping table.

    Input:
        transcript_id
        sequence_length

    Optional metadata:
        gene_symbol
        gene_id
        assembly_name

    Output mapping columns:
        sequence_id
        chrom
        start
        end
        strand
        transcript_start
        transcript_end
        exon_id
        assembly_name
        mapping_status

    Transcript coordinates are 0-based, half-open.
    Genomic coordinates are converted from Ensembl's
    1-based, closed coordinates to 0-based, half-open.

    Returns:
        (mapping_dataframe, status_dataframe)
    """

    missing = REQUIRED_TRANSCRIPT_COLUMNS - set(transcripts.columns)

    if missing:
        raise MappingError(
            "Transcript table is missing columns: "
            + ", ".join(sorted(missing))
        )

    mapping_rows = []
    status_rows = []

    selected = transcripts.copy()

    selected = selected.dropna(
        subset=["transcript_id"]
    )

    selected["transcript_id"] = (
        selected["transcript_id"].astype(str).str.strip()
    )

    selected = selected[
        selected["transcript_id"] != ""
    ]

    selected = selected.drop_duplicates(
        subset=["transcript_id"]
    )

    for index, transcript in selected.iterrows():

        transcript_id = str(transcript["transcript_id"])

        expected_length = _get_value(
            transcript,
            "sequence_length",
        )

        gene_symbol = _get_value(
            transcript,
            "gene_symbol",
            "",
        )

        gene_id = _get_value(
            transcript,
            "gene_id",
            "",
        )

        try:

            expected_length = int(expected_length)

            if expected_length <= 0:
                raise MappingError(
                    "Transcript sequence length must be positive."
                )

            annotation = _get_transcript_annotation(
                transcript_id
            )

            if annotation.get("object_type") != "Transcript":
                raise MappingError(
                    "Ensembl did not return a transcript object."
                )

            exons = annotation.get("Exon", [])

            if not exons:
                raise MappingError(
                    "No exon features returned by Ensembl."
                )

            transcript_strand = int(
                annotation.get("strand", 0)
            )

            if transcript_strand not in (1, -1):
                raise MappingError(
                    "Invalid transcript strand."
                )

            chromosome = str(
                annotation.get("seq_region_name", "")
            )

            assembly = str(
                annotation.get("assembly_name", "")
            )

            if not chromosome:
                raise MappingError(
                    "Chromosome information is missing."
                )

            if not assembly:
                raise MappingError(
                    "Genome assembly information is missing."
                )

            parsed_exons = []

            for exon in exons:

                if not isinstance(exon, dict):
                    continue

                exon_start = int(exon["start"])
                exon_end = int(exon["end"])

                exon_chrom = str(
                    exon.get(
                        "seq_region_name",
                        chromosome,
                    )
                )

                exon_strand = int(
                    exon.get(
                        "strand",
                        transcript_strand,
                    )
                )

                if exon_start < 1 or exon_end < exon_start:
                    raise MappingError(
                        "Invalid genomic exon coordinates."
                    )

                if exon_chrom != chromosome:
                    raise MappingError(
                        "Exons map to different chromosomes."
                    )

                if exon_strand != transcript_strand:
                    raise MappingError(
                        "Exon strand differs from transcript strand."
                    )

                if exon_chrom != chromosome:
                    raise MappingError(
                        "Exon chromosome mismatch."
                    )

                parsed_exons.append({
                    "exon_id": exon.get("id", ""),
                    "genomic_start_1based": exon_start,
                    "genomic_end_1based": exon_end,
                    "length": exon_end - exon_start + 1,
                })

            if not parsed_exons:
                raise MappingError(
                    "No valid exon coordinates were found."
                )

            # Transcript order follows transcription direction.
            if transcript_strand == 1:

                parsed_exons.sort(
                    key=lambda item: item["genomic_start_1based"]
                )

            else:

                parsed_exons.sort(
                    key=lambda item: item["genomic_start_1based"],
                    reverse=True,
                )

            total_exon_length = sum(
                exon["length"]
                for exon in parsed_exons
            )

            if total_exon_length != expected_length:

                raise MappingError(
                    "Exon length mismatch: "
                    f"Ensembl exon total={total_exon_length}, "
                    f"retrieved cDNA length={expected_length}. "
                    "Mapping was not accepted."
                )

            transcript_position = 0

            for exon in parsed_exons:

                exon_length = exon["length"]

                # Convert Ensembl 1-based closed genomic
                # coordinates to 0-based half-open coordinates.
                genomic_start = (
                    exon["genomic_start_1based"] - 1
                )

                genomic_end = (
                    exon["genomic_end_1based"]
                )

                transcript_start = transcript_position

                transcript_end = (
                    transcript_position + exon_length
                )

                mapping_rows.append({
                    "sequence_id": transcript_id,
                    "gene_symbol": gene_symbol,
                    "gene_id": gene_id,
                    "chrom": chromosome,
                    "start": genomic_start,
                    "end": genomic_end,
                    "strand": (
                        "+" if transcript_strand == 1 else "-"
                    ),
                    "transcript_start": transcript_start,
                    "transcript_end": transcript_end,
                    "exon_id": exon["exon_id"],
                    "assembly_name": assembly,
                    "mapping_status": "validated",
                })

                transcript_position = transcript_end

            status_rows.append({
                "transcript_id": transcript_id,
                "gene_symbol": gene_symbol,
                "assembly_name": assembly,
                "exon_count": len(parsed_exons),
                "cdna_length": expected_length,
                "mapped_exon_length": total_exon_length,
                "mapping_status": "validated",
                "message": "Exon lengths match cDNA length.",
            })

        except Exception as exc:

            status_rows.append({
                "transcript_id": transcript_id,
                "gene_symbol": gene_symbol,
                "assembly_name": "",
                "exon_count": 0,
                "cdna_length": expected_length,
                "mapped_exon_length": 0,
                "mapping_status": "failed",
                "message": str(exc),
            })

        if delay_seconds > 0:
            time.sleep(delay_seconds)

    mapping_columns = [
        "sequence_id",
        "gene_symbol",
        "gene_id",
        "chrom",
        "start",
        "end",
        "strand",
        "transcript_start",
        "transcript_end",
        "exon_id",
        "assembly_name",
        "mapping_status",
    ]

    status_columns = [
        "transcript_id",
        "gene_symbol",
        "assembly_name",
        "exon_count",
        "cdna_length",
        "mapped_exon_length",
        "mapping_status",
        "message",
    ]

    mapping_df = pd.DataFrame(
        mapping_rows,
        columns=mapping_columns,
    )

    status_df = pd.DataFrame(
        status_rows,
        columns=status_columns,
    )

    return mapping_df, status_df
