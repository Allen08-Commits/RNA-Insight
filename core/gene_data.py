
"""Retrieve real human gene and transcript data from Ensembl REST.

Data retrieved:
- Ensembl gene identifiers
- Gene names and biotypes
- Transcript identifiers
- Transcript sequences (cDNA)
- Chromosome, genomic coordinates, strand and assembly

Important:
Transcript cDNA coordinates are not interchangeable with genomic
coordinates because mature transcripts can contain spliced exons.
This module does not create a false one-interval transcript mapping.
"""

from __future__ import annotations

import re
import time
from typing import Any

import pandas as pd
import requests


ENSEMBL_BASE = "https://rest.ensembl.org"

HEADERS_JSON = {
    "Accept": "application/json",
    "Content-Type": "application/json",
    "User-Agent": "RNA-Insight/1.0 research application",
}

HEADERS_FASTA = {
    "Accept": "text/x-fasta",
    "User-Agent": "RNA-Insight/1.0 research application",
}

DEFAULT_GENES = [
    "HSPA1A",
    "HSP90AA1",
    "DDIT3",
    "ATF4",
]

GENE_COLUMNS = [
    "gene_symbol",
    "gene_id",
    "gene_name",
    "biotype",
    "chromosome",
    "genomic_start",
    "genomic_end",
    "strand",
    "assembly",
]

TRANSCRIPT_COLUMNS = [
    "gene_symbol",
    "gene_id",
    "transcript_id",
    "transcript_name",
    "biotype",
    "is_canonical",
    "chromosome",
    "genomic_start",
    "genomic_end",
    "strand",
    "assembly",
    "sequence_length",
    "sequence_status",
]


class EnsemblError(RuntimeError):
    """Raised when Ensembl data cannot be retrieved or validated."""


def _get_json(
    endpoint: str,
    params: dict[str, Any] | None = None,
) -> Any:
    """Request a JSON response from Ensembl."""
    url = f"{ENSEMBL_BASE}{endpoint}"

    try:
        response = requests.get(
            url,
            headers=HEADERS_JSON,
            params=params,
            timeout=30,
        )

        if response.status_code == 429:
            time.sleep(2)
            response = requests.get(
                url,
                headers=HEADERS_JSON,
                params=params,
                timeout=30,
            )

        response.raise_for_status()
        return response.json()

    except requests.RequestException as exc:
        raise EnsemblError(
            f"Ensembl request failed for {endpoint}: {exc}"
        ) from exc

    except ValueError as exc:
        raise EnsemblError(
            f"Ensembl returned invalid JSON for {endpoint}."
        ) from exc


def _get_sequence(transcript_id: str) -> str:
    """Retrieve a transcript cDNA sequence in FASTA format."""
    url = f"{ENSEMBL_BASE}/sequence/id/{transcript_id}"

    try:
        response = requests.get(
            url,
            headers=HEADERS_FASTA,
            params={"type": "cdna"},
            timeout=45,
        )
        response.raise_for_status()

    except requests.RequestException as exc:
        raise EnsemblError(
            f"Could not retrieve transcript {transcript_id}: {exc}"
        ) from exc

    lines = response.text.strip().splitlines()

    if not lines or not lines[0].startswith(">"):
        raise EnsemblError(
            f"Unexpected FASTA response for {transcript_id}."
        )

    sequence = "".join(
        line.strip()
        for line in lines[1:]
        if line.strip()
    ).upper()

    if not sequence:
        raise EnsemblError(
            f"Empty cDNA sequence returned for {transcript_id}."
        )

    allowed = set("ACGTUNRYKMSWBDHV")

    if set(sequence) - allowed:
        raise EnsemblError(
            f"Unexpected nucleotide characters in {transcript_id}."
        )

    return sequence.replace("T", "U")


def _find_gene(gene_symbol: str) -> dict[str, Any]:
    """Resolve a human gene symbol to an Ensembl gene record."""
    symbol = gene_symbol.strip().upper()

    if not re.fullmatch(r"[A-Z0-9_.-]+", symbol):
        raise ValueError(
            f"Invalid gene symbol: {gene_symbol}"
        )

    matches = _get_json(
        f"/xrefs/symbol/homo_sapiens/{symbol}"
    )

    gene_matches = [
        item
        for item in matches
        if item.get("type") == "gene"
        and item.get("id")
    ]

    if not gene_matches:
        raise EnsemblError(
            f"No Ensembl gene record found for {symbol}."
        )

    # Prefer an exact display-name match when multiple xrefs exist.
    exact_matches = [
        item
        for item in gene_matches
        if str(item.get("display_id", "")).upper() == symbol
    ]

    selected = exact_matches or gene_matches

    unique_ids = sorted({
        item["id"] for item in selected
    })

    if len(unique_ids) > 1:
        raise EnsemblError(
            f"Multiple Ensembl gene records found for {symbol}: "
            + ", ".join(unique_ids)
            + ". Resolve the gene identifier before proceeding."
        )

    gene_id = unique_ids[0]

    gene_data = _get_json(
        f"/lookup/id/{gene_id}",
        params={"expand": 1},
    )

    if gene_data.get("object_type") != "Gene":
        raise EnsemblError(
            f"Resolved identifier {gene_id} is not a gene."
        )

    return gene_data


def fetch_human_gene(
    gene_symbol: str,
    max_transcripts: int = 20,
) -> dict[str, Any]:
    """
    Retrieve one human gene and its transcript cDNA sequences.

    The canonical transcript is prioritized. Other transcripts
    are retrieved up to max_transcripts.

    Returns a dictionary containing:
        gene: gene annotation dictionary
        transcripts: transcript metadata DataFrame
        fasta: multi-record RNA FASTA string
    """
    if max_transcripts < 1:
        raise ValueError(
            "max_transcripts must be at least 1."
        )

    gene = _find_gene(gene_symbol)

    transcript_records = gene.get("Transcript", [])

    if not transcript_records:
        raise EnsemblError(
            f"No transcripts found for {gene_symbol}."
        )

    transcript_records = sorted(
        transcript_records,
        key=lambda item: (
            not bool(item.get("is_canonical", 0)),
            item.get("id", ""),
        ),
    )

    transcript_records = transcript_records[:max_transcripts]

    gene_id = gene["id"]
    resolved_symbol = (
        gene.get("display_name")
        or gene_symbol.upper()
    )

    transcript_rows = []
    fasta_records = []

    for index, transcript in enumerate(transcript_records):

        transcript_id = transcript.get("id")

        if not transcript_id:
            continue

        # Avoid rapid bursts of requests.
        if index:
            time.sleep(0.12)

        sequence_status = "retrieved"
        sequence = ""

        try:
            sequence = _get_sequence(transcript_id)

        except EnsemblError:
            sequence_status = "sequence_retrieval_failed"

        transcript_rows.append({
            "gene_symbol": resolved_symbol,
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "transcript_name": transcript.get(
                "display_name",
                transcript_id,
            ),
            "biotype": transcript.get("biotype", "unknown"),
            "is_canonical": bool(
                transcript.get("is_canonical", 0)
            ),
            "chromosome": str(
                transcript.get("seq_region_name", "")
            ),
            "genomic_start": transcript.get("start"),
            "genomic_end": transcript.get("end"),
            "strand": transcript.get("strand"),
            "assembly": transcript.get(
                "assembly_name",
                gene.get("assembly_name", "unknown"),
            ),
            "sequence_length": len(sequence),
            "sequence_status": sequence_status,
        })

        if sequence:
            fasta_records.append(
                f">{transcript_id}|{resolved_symbol}|"
                f"gene_id={gene_id}|"
                f"assembly={transcript.get('assembly_name', 'unknown')}\n"
                f"{sequence}"
            )

    transcript_df = pd.DataFrame(
        transcript_rows,
        columns=TRANSCRIPT_COLUMNS,
    )

    if transcript_df.empty:
        raise EnsemblError(
            f"No transcript records could be processed for "
            f"{gene_symbol}."
        )

    gene_record = {
        "gene_symbol": resolved_symbol,
        "gene_id": gene_id,
        "gene_name": gene.get("description", ""),
        "biotype": gene.get("biotype", "unknown"),
        "chromosome": str(gene.get("seq_region_name", "")),
        "genomic_start": gene.get("start"),
        "genomic_end": gene.get("end"),
        "strand": gene.get("strand"),
        "assembly": gene.get("assembly_name", "unknown"),
    }

    return {
        "gene": gene_record,
        "transcripts": transcript_df,
        "fasta": "\n".join(fasta_records),
    }


def fetch_gene_panel(
    gene_symbols: list[str] | None = None,
    max_transcripts_per_gene: int = 20,
) -> dict[str, Any]:
    """
    Retrieve a panel of human genes.

    A failure for one gene is recorded without stopping the
    retrieval of the remaining genes.
    """
    symbols = gene_symbols or DEFAULT_GENES

    gene_rows = []
    transcript_frames = []
    fasta_records = []
    errors = []

    for symbol in symbols:

        try:
            result = fetch_human_gene(
                symbol,
                max_transcripts=max_transcripts_per_gene,
            )

            gene_rows.append(result["gene"])
            transcript_frames.append(result["transcripts"])

            if result["fasta"]:
                fasta_records.append(result["fasta"])

        except (EnsemblError, ValueError) as exc:
            errors.append({
                "gene_symbol": symbol,
                "error": str(exc),
            })

    genes_df = pd.DataFrame(
        gene_rows,
        columns=GENE_COLUMNS,
    )

    if transcript_frames:
        transcripts_df = pd.concat(
            transcript_frames,
            ignore_index=True,
        )
    else:
        transcripts_df = pd.DataFrame(
            columns=TRANSCRIPT_COLUMNS
        )

    return {
        "genes": genes_df,
        "transcripts": transcripts_df,
        "fasta": "\n".join(fasta_records),
        "errors": pd.DataFrame(errors),
        "source": "Ensembl REST API",
        "species": "Homo sapiens",
    }
