
"""Human gene and transcript retrieval using Ensembl REST API."""

import json
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd


BASE_URL = "https://rest.ensembl.org"
SPECIES = "homo_sapiens"

HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/json",
    "User-Agent": "RNA-Insight/1.0",
}


class EnsemblError(Exception):
    """Ensembl retrieval error."""


def _get(endpoint, accept="application/json"):
    url = BASE_URL + endpoint

    headers = {
        "Accept": accept,
        "User-Agent": "RNA-Insight/1.0",
    }

    for attempt in range(3):
        request = urllib.request.Request(
            url,
            headers=headers,
        )

        try:
            with urllib.request.urlopen(
                request, timeout=35
            ) as response:
                content = response.read().decode("utf-8")

            if accept == "application/json":
                return json.loads(content)

            return content

        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504):
                if attempt < 2:
                    time.sleep(2 ** attempt)
                    continue

            raise EnsemblError(
                f"Ensembl HTTP {exc.code}: {exc.reason}"
            ) from exc

        except urllib.error.URLError as exc:
            if attempt < 2:
                time.sleep(2 ** attempt)
                continue

            raise EnsemblError(
                f"Ensembl connection failed: {exc}"
            ) from exc

    raise EnsemblError("Ensembl request failed.")


def _lookup_gene(gene_id):
    endpoint = (
        f"/lookup/id/{urllib.parse.quote(gene_id)}?expand=1"
    )

    result = _get(endpoint)

    if not isinstance(result, dict):
        raise EnsemblError(
            f"Invalid gene annotation for {gene_id}"
        )

    return result


def _find_gene(symbol):
    """Resolve an exact human gene symbol without arbitrary selection."""

    symbol = symbol.strip().upper()

    endpoint = (
        f"/xrefs/symbol/{SPECIES}/"
        f"{urllib.parse.quote(symbol)}?object_type=gene"
    )

    xrefs = _get(endpoint)

    if not isinstance(xrefs, list):
        raise EnsemblError(
            f"Unexpected Ensembl response for {symbol}"
        )

    gene_ids = sorted({
        item.get("id")
        for item in xrefs
        if item.get("id")
        and item.get("type") == "gene"
    })

    if not gene_ids:
        raise EnsemblError(
            f"No Ensembl gene candidates found for {symbol}"
        )

    candidates = []

    for gene_id in gene_ids:
        try:
            details = _lookup_gene(gene_id)

            # Ensembl lookup records commonly expose the symbol
            # through display_name, not external_name.
            returned_symbol = (
                details.get("display_name")
                or details.get("external_name")
                or ""
            )

            if returned_symbol.strip().upper() == symbol:
                candidates.append((gene_id, details))

        except EnsemblError:
            continue

    if not candidates:
        raise EnsemblError(
            f"No exact display-name match found for {symbol}. "
            f"Candidate IDs: {', '.join(gene_ids)}"
        )

    # Prefer protein-coding gene records.
    coding = [
        item for item in candidates
        if item[1].get("biotype") == "protein_coding"
    ]

    if coding:
        candidates = coding

    # Prefer a unique record with a canonical transcript.
    canonical = [
        item for item in candidates
        if item[1].get("canonical_transcript")
    ]

    if len(canonical) == 1:
        candidates = canonical

    if len(candidates) != 1:
        descriptions = "; ".join(
            f"{gene_id} "
            f"(biotype={details.get('biotype')}, "
            f"assembly={details.get('assembly_name')})"
            for gene_id, details in candidates
        )

        raise EnsemblError(
            f"Ambiguous gene symbol {symbol}. "
            f"Matching records: {descriptions}. "
            "Please verify the Ensembl gene ID."
        )

    return candidates[0]


def _get_cdna(transcript_id):
    endpoint = (
        f"/sequence/id/{urllib.parse.quote(transcript_id)}"
        "?type=cdna"
    )

    sequence = _get(endpoint, accept="text/plain")

    sequence = "".join(sequence.split()).upper()
    sequence = sequence.replace("T", "U")

    if not sequence:
        raise EnsemblError(
            f"Empty cDNA sequence for {transcript_id}"
        )

    return sequence


def fetch_human_gene(symbol, max_transcripts=5):
    """Retrieve gene annotation, transcript metadata and cDNA."""

    symbol = symbol.strip().upper()

    if not symbol:
        raise EnsemblError("Gene symbol cannot be empty.")

    gene_id, gene = _find_gene(symbol)

    canonical_id = gene.get("canonical_transcript")

    transcript_records = gene.get("Transcript", [])

    if not isinstance(transcript_records, list):
        transcript_records = []

    def transcript_rank(transcript):
        transcript_id = transcript.get("id", "")

        is_canonical = (
            transcript_id == canonical_id
            or transcript_id.split(".")[0]
            == str(canonical_id).split(".")[0]
        )

        is_coding = (
            transcript.get("biotype") == "protein_coding"
        )

        return (
            not is_canonical,
            not is_coding,
            transcript_id,
        )

    transcript_records = sorted(
        transcript_records,
        key=transcript_rank,
    )[:max(1, int(max_transcripts))]

    gene_row = {
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
                "strand", gene.get("strand")
            ),
            "assembly": gene.get("assembly_name", ""),
            "sequence_status": "pending",
            "sequence_length": 0,
        }

        try:
            sequence = _get_cdna(transcript_id)

            row["sequence_status"] = "retrieved"
            row["sequence_length"] = len(sequence)

            fasta_records.append(
                f">{symbol}|{gene_id}|{transcript_id}\n"
                f"{sequence}\n"
            )

        except EnsemblError as exc:
            row["sequence_status"] = f"failed: {exc}"

        transcript_rows.append(row)

    return {
        "gene": gene_row,
        "transcripts": transcript_rows,
        "fasta": "".join(fasta_records),
    }


def fetch_gene_panel(
    gene_symbols,
    max_transcripts_per_gene=5,
):
    """Retrieve multiple human genes independently."""

    genes = []
    transcripts = []
    fasta_parts = []
    errors = []

    for symbol in gene_symbols:
        symbol = str(symbol).strip().upper()

        if not symbol:
            continue

        try:
            result = fetch_human_gene(
                symbol,
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
        "gene_symbol", "gene_id", "gene_name",
        "biotype", "chromosome", "genomic_start",
        "genomic_end", "strand", "assembly",
    ]

    transcript_columns = [
        "gene_symbol", "gene_id", "transcript_id",
        "transcript_name", "biotype", "is_canonical",
        "chromosome", "genomic_start", "genomic_end",
        "strand", "assembly", "sequence_status",
        "sequence_length",
    ]

    return {
        "genes": pd.DataFrame(genes, columns=gene_columns),
        "transcripts": pd.DataFrame(
            transcripts, columns=transcript_columns
        ),
        "errors": pd.DataFrame(
            errors, columns=["gene_symbol", "error"]
        ),
        "fasta": "".join(fasta_parts),
    }
