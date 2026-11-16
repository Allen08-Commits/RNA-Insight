
"""Human gene and transcript retrieval using the Ensembl REST API."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


BASE_URL = "https://rest.ensembl.org"
SPECIES = "homo_sapiens"

DEFAULT_GENES = [
    "HSPA1A",
    "HSP90AA1",
    "DDIT3",
    "ATF4",
]

HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/json",
    "User-Agent": "RNA-Insight/1.0 research application",
}


class EnsemblError(Exception):
    """Raised when Ensembl retrieval or gene resolution fails."""


def _request_json(endpoint: str, retries: int = 3):
    """Request JSON from Ensembl with basic retry handling."""

    url = BASE_URL + endpoint

    for attempt in range(retries):
        request = urllib.request.Request(
            url,
            headers=HEADERS,
            method="GET",
        )

        try:
            with urllib.request.urlopen(
                request,
                timeout=35,
            ) as response:
                return json.loads(
                    response.read().decode("utf-8")
                )

        except urllib.error.HTTPError as exc:

            if exc.code in (429, 500, 502, 503, 504):
                if attempt < retries - 1:
                    time.sleep(2 ** attempt)
                    continue

            detail = exc.read().decode(
                "utf-8",
                errors="replace",
            )

            raise EnsemblError(
                f"Ensembl returned HTTP {exc.code}: {detail[:300]}"
            ) from exc

        except urllib.error.URLError as exc:

            if attempt < retries - 1:
                time.sleep(2 ** attempt)
                continue

            raise EnsemblError(
                f"Could not connect to Ensembl: {exc}"
            ) from exc

        except json.JSONDecodeError as exc:
            raise EnsemblError(
                "Ensembl returned an invalid JSON response."
            ) from exc

    raise EnsemblError("Ensembl request failed after retries.")


def _request_sequence(transcript_id: str) -> str:
    """Retrieve transcript cDNA and normalize it to RNA."""

    endpoint = (
        f"/sequence/id/{urllib.parse.quote(transcript_id)}"
        "?type=cdna"
    )

    url = BASE_URL + endpoint

    headers = {
        "Accept": "text/plain",
        "User-Agent": HEADERS["User-Agent"],
    }

    request = urllib.request.Request(
        url,
        headers=headers,
        method="GET",
    )

    for attempt in range(3):
        try:
            with urllib.request.urlopen(
                request,
                timeout=35,
            ) as response:
                sequence = response.read().decode("utf-8")

            sequence = "".join(
                sequence.split()
            ).upper()

            if sequence.startswith(">"):
                sequence = "".join(
                    sequence.splitlines()[1:]
                ).upper()

            sequence = sequence.replace("T", "U")

            if not sequence:
                raise EnsemblError(
                    f"Empty cDNA sequence for {transcript_id}."
                )

            return sequence

        except urllib.error.HTTPError as exc:

            if exc.code in (429, 500, 502, 503, 504):
                if attempt < 2:
                    time.sleep(2 ** attempt)
                    continue

            raise EnsemblError(
                f"Sequence retrieval failed for "
                f"{transcript_id}: HTTP {exc.code}"
            ) from exc

        except urllib.error.URLError as exc:

            if attempt < 2:
                time.sleep(2 ** attempt)
                continue

            raise EnsemblError(
                f"Could not retrieve sequence for "
                f"{transcript_id}: {exc}"
            ) from exc

    raise EnsemblError(
        f"Sequence retrieval failed for {transcript_id}."
    )


def _lookup_gene(gene_id: str) -> dict:
    """Retrieve expanded gene annotation."""

    endpoint = (
        f"/lookup/id/{urllib.parse.quote(gene_id)}"
        "?expand=1"
    )

    result = _request_json(endpoint)

    if not isinstance(result, dict):
        raise EnsemblError(
            f"Unexpected lookup response for {gene_id}."
        )

    if result.get("object_type") != "Gene":
        raise EnsemblError(
            f"{gene_id} is not an Ensembl gene record."
        )

    return result


def _find_gene(symbol: str) -> tuple[str, dict]:
    """
    Resolve a gene symbol to one exact Ensembl gene.

    Exact external_name matches are preferred.
    Protein-coding records are preferred over pseudogenes
    or other biotypes. If multiple equally suitable records
    remain, the function reports ambiguity instead of
    silently selecting an arbitrary gene.
    """

    encoded_symbol = urllib.parse.quote(symbol)

    endpoint = (
        f"/xrefs/symbol/{SPECIES}/{encoded_symbol}"
        "?object_type=gene"
    )

    xrefs = _request_json(endpoint)

    if not isinstance(xrefs, list):
        raise EnsemblError(
            f"Unexpected xref response for {symbol}."
        )

    gene_ids = sorted({
        item.get("id")
        for item in xrefs
        if item.get("type") == "gene"
        and item.get("id")
    })

    if not gene_ids:
        raise EnsemblError(
            f"No Ensembl gene records found for {symbol}."
        )

    candidates = []

    for gene_id in gene_ids:
        try:
            details = _lookup_gene(gene_id)

            external_name = str(
                details.get("external_name", "")
            ).upper()

            if external_name != symbol.upper():
                continue

            candidates.append(
                (gene_id, details)
            )

        except EnsemblError:
            continue

    if not candidates:
        raise EnsemblError(
            f"No exact gene-symbol match found for {symbol}. "
            f"Candidate Ensembl IDs: {', '.join(gene_ids)}"
        )

    # Prefer protein-coding genes.
    protein_coding = [
        item
        for item in candidates
        if item[1].get("biotype") == "protein_coding"
    ]

    if protein_coding:
        candidates = protein_coding

    # Prefer records with a canonical transcript annotation.
    canonical = [
        item
        for item in candidates
        if item[1].get("canonical_transcript")
    ]

    if len(canonical) == 1:
        candidates = canonical

    if len(candidates) > 1:
        candidate_summary = "; ".join(
            f"{gene_id} "
            f"(biotype={details.get('biotype')}, "
            f"assembly={details.get('assembly_name', 'unknown')})"
            for gene_id, details in candidates
        )

        raise EnsemblError(
            f"Multiple equally suitable Ensembl records remain "
            f"for {symbol}: {candidate_summary}. "
            "Please resolve the gene ID using an authoritative "
            "gene annotation source before analysis."
        )

    return candidates[0]


def fetch_human_gene(
    symbol: str,
    max_transcripts: int = 5,
) -> dict:
    """Retrieve one human gene, its transcripts and cDNA sequences."""

    symbol = symbol.strip().upper()

    if not symbol:
        raise EnsemblError("Gene symbol cannot be empty.")

    gene_id, gene = _find_gene(symbol)

    canonical_id = gene.get("canonical_transcript")

    transcript_records = gene.get("Transcript", [])

    if not isinstance(transcript_records, list):
        transcript_records = []

    # Rank canonical and protein-coding transcripts first.
    def transcript_rank(transcript: dict):
        transcript_id = transcript.get("id", "")

        is_canonical = (
            transcript_id == canonical_id
            or transcript_id.split(".")[0]
            == str(canonical_id).split(".")[0]
        )

        is_protein_coding = (
            transcript.get("biotype") == "protein_coding"
        )

        return (
            not is_canonical,
            not is_protein_coding,
            transcript_id,
        )

    transcript_records = sorted(
        transcript_records,
        key=transcript_rank,
    )

    transcript_records = transcript_records[
        :max_transcripts
    ]

    gene_record = {
        "gene_symbol": symbol,
        "gene_id": gene_id,
        "gene_name": gene.get("description", ""),
        "biotype": gene.get("biotype", ""),
        "chromosome": gene.get("seq_region_name", ""),
        "genomic_start": gene.get("start"),
        "genomic_end": gene.get("end"),
        "strand": gene.get("strand"),
        "assembly": gene.get("assembly_name", ""),
    }

    transcript_rows = []
    fasta_records = []

    for transcript in transcript_records:

        transcript_id = transcript.get("id")

        if not transcript_id:
            continue

        is_canonical = (
            transcript_id == canonical_id
            or transcript_id.split(".")[0]
            == str(canonical_id).split(".")[0]
        )

        row = {
            "gene_symbol": symbol,
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "transcript_name": transcript.get(
                "display_name", ""
            ),
            "biotype": transcript.get("biotype", ""),
            "is_canonical": is_canonical,
            "chromosome": transcript.get(
                "seq_region_name",
                gene.get("seq_region_name", ""),
            ),
            "genomic_start": transcript.get("start"),
            "genomic_end": transcript.get("end"),
            "strand": transcript.get(
                "strand",
                gene.get("strand"),
            ),
            "assembly": gene.get("assembly_name", ""),
            "sequence_status": "pending",
            "sequence_length": 0,
        }

        try:
            sequence = _request_sequence(transcript_id)

            row["sequence_status"] = "retrieved"
            row["sequence_length"] = len(sequence)

            fasta_records.append(
                f">{symbol}|{gene_id}|{transcript_id}\n"
                f"{sequence}\n"
            )

        except EnsemblError as exc:
            row["sequence_status"] = (
                f"failed: {str(exc)[:180]}"
            )

        transcript_rows.append(row)

    return {
        "gene": gene_record,
        "transcripts": transcript_rows,
        "fasta": "".join(fasta_records),
    }


def fetch_gene_panel(
    gene_symbols: list[str],
    max_transcripts_per_gene: int = 5,
) -> dict:
    """Retrieve a panel of human genes without stopping on one failure."""

    genes = []
    transcripts = []
    fasta_parts = []
    errors = []

    normalized_symbols = [
        str(symbol).strip().upper()
        for symbol in gene_symbols
        if str(symbol).strip()
    ]

    for symbol in normalized_symbols:

        try:
            result = fetch_human_gene(
                symbol=symbol,
                max_transcripts=max_transcripts_per_gene,
            )

            genes.append(result["gene"])
            transcripts.extend(result["transcripts"])

            if result["fasta"]:
                fasta_parts.append(result["fasta"])

        except Exception as exc:
            errors.append({
                "gene_symbol": symbol,
                "error": str(exc),
            })

    gene_columns = [
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

    transcript_columns = [
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
        "sequence_status",
        "sequence_length",
    ]

    error_columns = [
        "gene_symbol",
        "error",
    ]

    return {
        "genes": pd.DataFrame(
            genes,
            columns=gene_columns,
        ),
        "transcripts": pd.DataFrame(
            transcripts,
            columns=transcript_columns,
        ),
        "errors": pd.DataFrame(
            errors,
            columns=error_columns,
        ),
        "fasta": "".join(fasta_parts),
    }
