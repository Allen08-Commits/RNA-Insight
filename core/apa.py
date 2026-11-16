
"""APA annotation utilities with automatic PolyA_DB v4 retrieval.

PolyA_DB v4.1 human annotations are mapped to hg38.
The downloaded archive is parsed dynamically because its internal
file layout and column names may change between database releases.
"""

from __future__ import annotations

from io import BytesIO
from zipfile import ZipFile

import pandas as pd
import requests


POLYA_DB_URL = (
    "https://exon.apps.wistar.org/"
    "polya_db/v4/download/4.1/HumanPas.zip"
)

POLYA_DB_ASSEMBLY = "hg38"

REQUIRED_APA_COLUMNS = {
    "sequence_id",
    "chrom",
    "start",
    "end",
    "strand",
}

GENE_SYMBOLS = [
    "HSPA1A",
    "HSP90AA1",
    "DDIT3",
    "ATF4",
]


def _clean_column_name(value: str) -> str:
    """Normalize database column names for matching."""
    return (
        str(value)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
    )


def _read_database_archive(content: bytes) -> list[pd.DataFrame]:
    """Read supported tabular files from the downloaded ZIP archive."""
    frames = []

    with ZipFile(BytesIO(content)) as archive:
        files = [
            name for name in archive.namelist()
            if not name.endswith("/")
            and not name.startswith("__MACOSX/")
        ]

        for filename in files:
            lower_name = filename.lower()

            try:
                with archive.open(filename) as file:
                    if lower_name.endswith((".tsv", ".txt", ".bed")):
                        frame = pd.read_csv(
                            file,
                            sep="\t",
                            comment="#",
                            low_memory=False,
                        )

                    elif lower_name.endswith(".csv"):
                        frame = pd.read_csv(file, low_memory=False)

                    elif lower_name.endswith(".gz"):
                        frame = pd.read_csv(
                            file,
                            sep="\t",
                            compression="gzip",
                            comment="#",
                            low_memory=False,
                        )

                    else:
                        continue

                if not frame.empty:
                    frame.columns = [
                        _clean_column_name(c)
                        for c in frame.columns
                    ]
                    frames.append(frame)

            except Exception:
                continue

    return frames


def _find_column(frame: pd.DataFrame, aliases: list[str]):
    """Find a column using known naming alternatives."""
    for name in aliases:
        if name in frame.columns:
            return name
    return None


def _extract_gene_annotations(
    frames: list[pd.DataFrame],
    gene_symbols: list[str],
) -> pd.DataFrame:
    """Identify PAS records and standardize their genomic fields."""

    gene_aliases = [
        "gene_symbol",
        "genesymbol",
        "gene_name",
        "symbol",
        "gene",
        "refseq_gene_symbol",
    ]

    chromosome_aliases = [
        "chrom",
        "chr",
        "chromosome",
        "chromosome_name",
    ]

    start_aliases = [
        "start",
        "genomic_start",
        "position_start",
    ]

    end_aliases = [
        "end",
        "genomic_end",
        "position_end",
    ]

    position_aliases = [
        "position",
        "pas_position",
        "cleavage_site",
        "cs",
    ]

    strand_aliases = [
        "strand",
        "orientation",
    ]

    type_aliases = [
        "pas_type",
        "site_type",
        "type",
        "annotation",
    ]

    id_aliases = [
        "pas_id",
        "site_id",
        "id",
    ]

    requested = {
        str(symbol).upper()
        for symbol in gene_symbols
    }

    extracted = []

    for frame in frames:
        gene_col = _find_column(frame, gene_aliases)
        chrom_col = _find_column(frame, chromosome_aliases)
        strand_col = _find_column(frame, strand_aliases)

        start_col = _find_column(frame, start_aliases)
        end_col = _find_column(frame, end_aliases)
        position_col = _find_column(frame, position_aliases)

        if not all([gene_col, chrom_col, strand_col]):
            continue

        if not start_col and not position_col:
            continue

        working = frame.copy()

        working["_gene"] = (
            working[gene_col]
            .astype(str)
            .str.upper()
        )

        working = working[
            working["_gene"].isin(requested)
        ].copy()

        if working.empty:
            continue

        if start_col and end_col:
            working["_start"] = pd.to_numeric(
                working[start_col],
                errors="coerce",
            )
            working["_end"] = pd.to_numeric(
                working[end_col],
                errors="coerce",
            )

        elif position_col:
            working["_start"] = pd.to_numeric(
                working[position_col],
                errors="coerce",
            )
            working["_end"] = working["_start"] + 1

        else:
            continue

        type_col = _find_column(working, type_aliases)
        id_col = _find_column(working, id_aliases)

        result = pd.DataFrame({
            "sequence_id": working["_gene"],
            "gene_symbol": working["_gene"],
            "chrom": working[chrom_col].astype(str),
            "start": working["_start"],
            "end": working["_end"],
            "strand": working[strand_col].astype(str),
            "site_type": (
                working[type_col].astype(str)
                if type_col
                else "unspecified"
            ),
            "assembly": POLYA_DB_ASSEMBLY,
            "source": "PolyA_DB_v4.1",
        })

        if id_col:
            result["pas_id"] = working[id_col].astype(str)

        extracted.append(result)

    if not extracted:
        raise ValueError(
            "No compatible APA annotation table was found in "
            "the downloaded archive. The database file structure "
            "may have changed. Inspect the archive before proceeding."
        )

    return pd.concat(extracted, ignore_index=True)


def fetch_apa_sites(
    gene_symbols: list[str] | None = None,
    timeout: int = 90,
) -> pd.DataFrame:
    """Download PolyA_DB v4.1 human APA annotations and filter genes.

    Returns gene-level genomic APA annotations.
    sequence_id currently contains the gene symbol, NOT an Ensembl
    transcript ID. Transcript mapping must be performed separately.
    """

    if gene_symbols is None:
        gene_symbols = GENE_SYMBOLS

    gene_symbols = [
        str(symbol).strip().upper()
        for symbol in gene_symbols
        if str(symbol).strip()
    ]

    if not gene_symbols:
        raise ValueError("Provide at least one gene symbol.")

    try:
        response = requests.get(
            POLYA_DB_URL,
            timeout=timeout,
            headers={"User-Agent": "RNA-Insight/1.0"},
        )
        response.raise_for_status()

    except requests.RequestException as exc:
        raise RuntimeError(
            f"Could not download PolyA_DB annotations: {exc}"
        ) from exc

    if not response.content.startswith(b"PK"):
        raise RuntimeError(
            "The PolyA_DB download did not return a ZIP archive. "
            "Check whether the official download service is available."
        )

    frames = _read_database_archive(response.content)

    if not frames:
        raise RuntimeError(
            "The PolyA_DB archive contained no readable tabular files."
        )

    result = _extract_gene_annotations(
        frames,
        gene_symbols,
    )

    result = standardize_apa_sites(result)

    result = result.drop_duplicates(
        subset=["gene_symbol", "chrom", "start", "end", "strand"]
    ).reset_index(drop=True)

    result.attrs["source"] = "PolyA_DB v4.1"
    result.attrs["assembly"] = POLYA_DB_ASSEMBLY
    result.attrs["genes_requested"] = gene_symbols

    return result


def standardize_apa_sites(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and standardize an APA annotation table."""

    result = frame.copy()

    result.columns = [
        _clean_column_name(column)
        for column in result.columns
    ]

    missing = REQUIRED_APA_COLUMNS - set(result.columns)

    if missing:
        raise ValueError(
            "APA table is missing required columns: "
            + ", ".join(sorted(missing))
        )

    result["sequence_id"] = result["sequence_id"].astype(str)
    result["chrom"] = result["chrom"].astype(str)
    result["start"] = pd.to_numeric(
        result["start"], errors="raise"
    )
    result["end"] = pd.to_numeric(
        result["end"], errors="raise"
    )
    result["strand"] = result["strand"].astype(str)

    if (result["start"] < 0).any():
        raise ValueError(
            "APA start coordinates cannot be negative."
        )

    if (result["end"] <= result["start"]).any():
        raise ValueError(
            "APA end must be greater than start."
        )

    if not result["strand"].isin(["+", "-", "."]).all():
        raise ValueError(
            "Strand values must be +, - or ."
        )

    if "site_type" not in result.columns:
        result["site_type"] = "unspecified"

    result["site_type"] = (
        result["site_type"]
        .fillna("unspecified")
        .astype(str)
        .str.lower()
    )

    return result


def summarize_apa_sites(frame: pd.DataFrame) -> dict:
    """Summarize the number and distribution of annotated APA sites."""

    if frame.empty:
        return {
            "total_sites": 0,
            "sequence_count": 0,
            "proximal_sites": 0,
            "distal_sites": 0,
            "other_sites": 0,
        }

    site_types = frame["site_type"].astype(str).str.lower()

    proximal = site_types.str.contains(
        "proximal", na=False
    ).sum()

    distal = site_types.str.contains(
        "distal", na=False
    ).sum()

    return {
        "total_sites": int(len(frame)),
        "sequence_count": int(
            frame["sequence_id"].nunique()
        ),
        "proximal_sites": int(proximal),
        "distal_sites": int(distal),
        "other_sites": int(
            len(frame) - proximal - distal
        ),
    }
