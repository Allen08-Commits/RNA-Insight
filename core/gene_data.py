
"""Human gene and transcript retrieval using the Ensembl REST API.

Returns gene annotations, transcript metadata and transcript cDNA FASTA.
Transcript cDNA coordinates are not genomic coordinates.
"""

from __future__ import annotations

import json
import time
import urllib.parse
from typing import Any

import pandas as pd
import requests


GENE_COLUMNS = [
    "gene_symbol",
    "gene_id",
    "gene_name",
    "biotype",
    "chromosome",
    "genomic_start",
    "genomic_end",
    "strand",
    "assembly_name",
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
    "assembly_name",
    "sequence_status",
    "sequence_length",
    "sequence_error",
]

ERROR_COLUMNS = [
    "gene_symbol",
    "error",
]


class EnsemblError(Exception):
    """Raised when an Ensembl REST request cannot be completed."""


def _get_json(endpoint: str) -> Any:
    """Fetch JSON, failing over across Ensembl REST mirrors on HTTP errors."""
    last_error = None
    forbidden_hosts = []

    for base_url in ENSEMBL_BASES:
        url = base_url + endpoint
        for attempt in range(2):
            try:
                response = requests.get(url, headers=HEADERS, timeout=35)
                if response.status_code == 403:
                    forbidden_hosts.append(base_url)
                    last_error = f"HTTP 403 Forbidden from {base_url}"
                    break  # do not repeatedly hit a host that has refused access
                if response.status_code == 429:
                    wait = response.headers.get("Retry-After", "3")
                    try:
                        wait_seconds = min(max(float(wait), 1.0), 30.0)
                    except (TypeError, ValueError):
                        wait_seconds = 3.0
                    last_error = f"HTTP 429 rate limit from {base_url}"
                    if attempt == 0:
                        time.sleep(wait_seconds)
                        continue
                    break
                if response.status_code >= 500:
                    last_error = f"HTTP {response.status_code} from {base_url}"
                    if attempt == 0:
                        time.sleep(1.5)
                        continue
                    break
                response.raise_for_status()
                try:
                    return response.json()
                except ValueError as exc:
                    last_error = f"Invalid JSON from {base_url}: {exc}"
                    break
            except requests.RequestException as exc:
                last_error = f"{base_url}: {exc}"
                if attempt == 0:
                    time.sleep(1.0)
                    continue
                break

    if forbidden_hosts:
        raise EnsemblError(
            "Ensembl refused requests (HTTP 403) on: "
            + ", ".join(forbidden_hosts)
            + ". The hosting IP may be temporarily restricted by Ensembl. "
            + "Wait before retrying or use the manual Ensembl export workflow."
        )
    raise EnsemblError(f"All Ensembl REST mirrors failed. Last error: {last_error}")


def _lookup_gene(symbol: str) -> dict:
    """Resolve a gene symbol to one Ensembl gene record.

    First attempts Ensembl's direct symbol lookup. If that fails,
    searches symbol cross-references and checks candidate gene records.
    """

    encoded_symbol = urllib.parse.quote(symbol, safe="")

    direct_endpoint = (
        f"/lookup/symbol/{SPECIES}/{encoded_symbol}"
        "?expand=1"
    )

    try:
        result = _get_json(direct_endpoint)

        if (
            isinstance(result, dict)
            and result.get("id")
            and result.get("object_type") == "Gene"
        ):
            return result

    except EnsemblError as exc:
        if "HTTP 403" in str(exc) or "All Ensembl REST mirrors failed" in str(exc):
            raise
        # Continue to the cross-reference lookup for genuine lookup failures.
        pass

    # Fallback: search external cross-references.
    xref_endpoint = (
        f"/xrefs/symbol/{SPECIES}/{encoded_symbol}"
        "?object_type=gene"
    )

    try:
        xrefs = _get_json(xref_endpoint)
    except EnsemblError as exc:
        raise EnsemblError(
            f"Could not resolve {symbol}: {exc}"
        ) from exc

    if not isinstance(xrefs, list):
        xrefs = []

    candidate_ids = []

    for item in xrefs:
        if not isinstance(item, dict):
            continue

        identifier = item.get("id")

        if identifier and identifier not in candidate_ids:
            candidate_ids.append(identifier)

    exact_matches = []

    for gene_id in candidate_ids:
        try:
            gene = _get_json(
                f"/lookup/id/{urllib.parse.quote(gene_id)}"
                "?expand=1"
            )
        except EnsemblError:
            continue

        if not isinstance(gene, dict):
            continue

        if gene.get("object_type") != "Gene":
            continue

        display_name = str(
            gene.get("display_name", "")
        ).upper()

        external_name = str(
            gene.get("external_name", "")
        ).upper()

        if (
            display_name == symbol.upper()
            or external_name == symbol.upper()
        ):
            exact_matches.append(gene)

    if not exact_matches:
        raise EnsemblError(
            f"No exact gene-symbol match found for {symbol}. "
            f"Candidate Ensembl IDs: "
            f"{', '.join(candidate_ids[:10]) or 'none'}"
        )

    if len(exact_matches) > 1:
        ids = [
            item.get("id", "unknown")
            for item in exact_matches
        ]

        raise EnsemblError(
            f"Multiple exact Ensembl gene records found for "
            f"{symbol}: {', '.join(ids)}. "
            "Resolve the gene identifier before proceeding."
        )

    # Expand transcripts for the selected gene.
    selected_id = exact_matches[0]["id"]

    return _get_json(
        f"/lookup/id/{urllib.parse.quote(selected_id)}"
        "?expand=1"
    )


def _fetch_cdna(transcript_id: str) -> str:
    """Retrieve a transcript's spliced cDNA sequence."""

    encoded_id = urllib.parse.quote(transcript_id, safe="")

    result = _get_json(
        f"/sequence/id/{encoded_id}?type=cdna"
    )

    if not isinstance(result, dict):
        raise EnsemblError(
            "Unexpected sequence response from Ensembl."
        )

    sequence = str(result.get("seq", "")).upper()

    if not sequence:
        raise EnsemblError(
            "Ensembl returned an empty cDNA sequence."
        )

    # Store RNA alphabet consistently.
    return sequence.replace("T", "U")


def _safe_int(value: Any):
    """Convert coordinate values to integers when possible."""

    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def fetch_gene_panel(
    gene_symbols: list[str],
    max_transcripts_per_gene: int = 5,
) -> dict:
    """Retrieve genes, transcript annotations and cDNA FASTA.

    Parameters
    ----------
    gene_symbols:
        Human gene symbols such as HSPA1A or ATF4.
    max_transcripts_per_gene:
        Maximum number of transcripts to retrieve per gene.

    Returns
    -------
    dict containing:
        genes       : pandas DataFrame
        transcripts : pandas DataFrame
        errors      : pandas DataFrame
        fasta       : FASTA-formatted string
    """

    if not gene_symbols:
        raise ValueError(
            "At least one human gene symbol is required."
        )

    if not 1 <= int(max_transcripts_per_gene) <= 20:
        raise ValueError(
            "max_transcripts_per_gene must be between 1 and 20."
        )

    genes_data = []
    transcripts_data = []
    errors_data = []
    fasta_records = []

    for raw_symbol in gene_symbols:

        symbol = str(raw_symbol).strip().upper()

        if not symbol:
            continue

        try:
            gene = _lookup_gene(symbol)

            gene_id = gene.get("id", "")
            chromosome = gene.get("seq_region_name", "")
            assembly = gene.get("assembly_name", "")

            gene_row = {
                "gene_symbol": symbol,
                "gene_id": gene_id,
                "gene_name": (
                    gene.get("description")
                    or gene.get("display_name")
                    or symbol
                ),
                "biotype": gene.get("biotype", ""),
                "chromosome": chromosome,
                "genomic_start": _safe_int(
                    gene.get("start")
                ),
                "genomic_end": _safe_int(
                    gene.get("end")
                ),
                "strand": gene.get("strand"),
                "assembly_name": assembly,
            }

            genes_data.append(gene_row)

            transcript_list = gene.get("Transcript", [])

            if not isinstance(transcript_list, list):
                transcript_list = []

            # Prefer canonical transcripts, then other transcripts.
            transcript_list = sorted(
                transcript_list,
                key=lambda item: (
                    not bool(item.get("is_canonical", 0)),
                    str(item.get("id", "")),
                ),
            )

            selected_transcripts = transcript_list[
                :int(max_transcripts_per_gene)
            ]

            if not selected_transcripts:
                errors_data.append({
                    "gene_symbol": symbol,
                    "error": (
                        "Gene annotation was retrieved, but no "
                        "transcripts were returned by Ensembl."
                    ),
                })

            for transcript in selected_transcripts:

                transcript_id = transcript.get("id", "")

                if not transcript_id:
                    continue

                transcript_row = {
                    "gene_symbol": symbol,
                    "gene_id": gene_id,
                    "transcript_id": transcript_id,
                    "transcript_name": (
                        transcript.get("display_name")
                        or transcript_id
                    ),
                    "biotype": transcript.get("biotype", ""),
                    "is_canonical": bool(
                        transcript.get("is_canonical", 0)
                    ),
                    "chromosome": transcript.get(
                        "seq_region_name",
                        chromosome,
                    ),
                    "genomic_start": _safe_int(
                        transcript.get("start")
                    ),
                    "genomic_end": _safe_int(
                        transcript.get("end")
                    ),
                    "strand": transcript.get(
                        "strand",
                        gene.get("strand"),
                    ),
                    "assembly_name": assembly,
                    "sequence_status": "pending",
                    "sequence_length": 0,
                    "sequence_error": "",
                }

                try:
                    sequence = _fetch_cdna(transcript_id)

                    transcript_row["sequence_status"] = (
                        "retrieved"
                    )
                    transcript_row["sequence_length"] = len(
                        sequence
                    )

                    fasta_records.append(
                        f">{transcript_id}|gene={symbol}|"
                        f"gene_id={gene_id}\n{sequence}"
                    )

                except EnsemblError as exc:
                    transcript_row["sequence_status"] = "failed"
                    transcript_row["sequence_error"] = str(exc)

                    errors_data.append({
                        "gene_symbol": symbol,
                        "error": (
                            f"Transcript {transcript_id}: {exc}"
                        ),
                    })

                transcripts_data.append(transcript_row)

                # Be polite to the public REST service.
                time.sleep(0.1)

        except Exception as exc:
            errors_data.append({
                "gene_symbol": symbol,
                "error": str(exc),
            })

    genes_df = pd.DataFrame(
        genes_data,
        columns=GENE_COLUMNS,
    )

    transcripts_df = pd.DataFrame(
        transcripts_data,
        columns=TRANSCRIPT_COLUMNS,
    )

    errors_df = pd.DataFrame(
        errors_data,
        columns=ERROR_COLUMNS,
    )

    fasta_text = (
        "\n".join(fasta_records) + "\n"
        if fasta_records
        else ""
    )

    return {
        "genes": genes_df,
        "transcripts": transcripts_df,
        "errors": errors_df,
        "fasta": fasta_text,
    }
