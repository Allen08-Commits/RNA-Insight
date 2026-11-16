
"""RNA-Insight APA utilities for PolyA_DB v4.1 human hg38 data."""

from __future__ import annotations

from io import BytesIO
from zipfile import ZipFile
import re

import pandas as pd
import requests


POLYA_DB_URL = (
    "https://exon.apps.wistar.org/"
    "polya_db/v4/download/4.1/HumanPas.zip"
)

POLYA_DB_ASSEMBLY = "hg38"

GENE_SYMBOLS = [
    "HSPA1A",
    "HSP90AA1",
    "DDIT3",
    "ATF4",
]

REQUIRED_APA_COLUMNS = {
    "sequence_id",
    "chrom",
    "start",
    "end",
    "strand",
}


def _clean_column_name(value):
    """Normalize database column names."""
    return (
        str(value)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
        .replace(".", "_")
    )


def _read_database_archive(content):
    """Read the Main human PAS table from the ZIP archive."""

    frames = []
    diagnostics = []

    with ZipFile(BytesIO(content)) as archive:

        files = [
            name
            for name in archive.namelist()
            if not name.endswith("/")
            and not name.startswith("__MACOSX/")
        ]

        diagnostics.extend(files)

        preferred = [
            name
            for name in files
            if name.lower().endswith(
                "hg38.pas.main.tsv"
            )
        ]

        if not preferred:
            preferred = [
                name
                for name in files
                if name.lower().endswith(".tsv")
            ]

        for filename in preferred:

            with archive.open(filename) as file:

                frame = pd.read_csv(
                    file,
                    sep="\t",
                    low_memory=False,
                    dtype=str,
                )

            frame.columns = [
                _clean_column_name(column)
                for column in frame.columns
            ]

            frames.append(frame)

    return frames, diagnostics


def _parse_pas_id(value):
    """
    Decode PolyA_DB PAS identifiers.

    Actual PolyA_DB format:
        chr12:-:57517847
        chr22:+:39521193

    Returns:
        chromosome, 1-based position, strand
    """

    value = str(value).strip()

    pattern = (
        r"^(chr(?:[0-9]+|X|Y|M|MT)):"
        r"([+-]):"
        r"(\d+)$"
    )

    match = re.match(
        pattern,
        value,
        flags=re.IGNORECASE,
    )

    if not match:
        return None

    chrom = match.group(1)
    strand = match.group(2)
    position = int(match.group(3))

    if position <= 0:
        return None

    return chrom, position, strand


def _extract_gene_annotations(frames, gene_symbols):
    """Extract APA annotations for the selected genes."""

    gene_aliases = [
        "genesymbol",
        "gene_symbol",
        "gene_name",
        "symbol",
    ]

    pas_id_aliases = [
        "pas_id",
        "site_id",
    ]

    type_aliases = [
        "pas_type_refseq_label",
        "pas_lrtype_label",
        "pas_type",
        "site_type",
    ]

    requested = {
        str(gene).strip().upper()
        for gene in gene_symbols
    }

    extracted = []
    pas_id_examples = []

    for frame in frames:

        gene_col = next(
            (
                column
                for column in gene_aliases
                if column in frame.columns
            ),
            None,
        )

        pas_id_col = next(
            (
                column
                for column in pas_id_aliases
                if column in frame.columns
            ),
            None,
        )

        type_col = next(
            (
                column
                for column in type_aliases
                if column in frame.columns
            ),
            None,
        )

        if not gene_col or not pas_id_col:
            continue

        working = frame.copy()

        working["_gene"] = (
            working[gene_col]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.upper()
        )

        working = working[
            working["_gene"].isin(requested)
        ].copy()

        if working.empty:
            continue

        pas_id_examples.extend(
            working[pas_id_col]
            .dropna()
            .astype(str)
            .head(10)
            .tolist()
        )

        working["_parsed"] = (
            working[pas_id_col].apply(_parse_pas_id)
        )

        working = working[
            working["_parsed"].notna()
        ].copy()

        if working.empty:
            continue

        working["chrom"] = (
            working["_parsed"].apply(
                lambda value: value[0]
            )
        )

        working["_position"] = (
            working["_parsed"].apply(
                lambda value: value[1]
            )
        )

        working["strand"] = (
            working["_parsed"].apply(
                lambda value: value[2]
            )
        )

        # Convert 1-based cleavage position
        # to 0-based, half-open BED coordinates.
        working["start"] = (
            working["_position"] - 1
        )

        working["end"] = (
            working["_position"]
        )

        if type_col:

            working["site_type"] = (
                working[type_col]
                .fillna("unspecified")
                .astype(str)
                .str.lower()
            )

        else:
            working["site_type"] = "unspecified"

        result = pd.DataFrame({

            "sequence_id": working["_gene"],

            "gene_symbol": working["_gene"],

            "pas_id": working[pas_id_col].astype(str),

            "chrom": working["chrom"],

            "start": working["start"],

            "end": working["end"],

            "strand": working["strand"],

            "site_type": working["site_type"],

            "assembly": POLYA_DB_ASSEMBLY,

            "source": "PolyA_DB_v4.1",
        })

        additional_columns = [
            "pas_type_refseq_label",
            "pas_lrtype_label",
            "pas2gene",
            "pas2gene_refseq",
            "pse_3reads",
            "avgrpm_3reads",
            "polya_strength",
            "polya_strength_percentile",
            "conservation",
        ]

        for column in additional_columns:

            if column in working.columns:

                result[column] = (
                    working[column].values
                )

        extracted.append(result)

    if not extracted:

        examples = (
            " | ".join(pas_id_examples[:15])
            if pas_id_examples
            else "No PAS identifiers found."
        )

        raise ValueError(
            "Could not extract APA records.\n"
            "Sample PAS IDs: "
            + examples
        )

    return pd.concat(
        extracted,
        ignore_index=True,
    )


def fetch_apa_sites(
    gene_symbols=None,
    timeout=90,
):
    """Download and retrieve human APA annotations."""

    if gene_symbols is None:
        gene_symbols = GENE_SYMBOLS

    gene_symbols = [
        str(gene).strip().upper()
        for gene in gene_symbols
        if str(gene).strip()
    ]

    if not gene_symbols:
        raise ValueError(
            "Provide at least one gene symbol."
        )

    try:

        response = requests.get(
            POLYA_DB_URL,
            timeout=timeout,
            headers={
                "User-Agent": "RNA-Insight/1.0"
            },
        )

        response.raise_for_status()

    except requests.RequestException as exc:

        raise RuntimeError(
            f"PolyA_DB download failed: {exc}"
        ) from exc

    if not response.content.startswith(b"PK"):

        raise RuntimeError(
            "The downloaded file is not a valid ZIP archive."
        )

    frames, archive_files = (
        _read_database_archive(response.content)
    )

    if not frames:

        raise RuntimeError(
            "No compatible TSV files were found. "
            f"Archive contents: {archive_files}"
        )

    result = _extract_gene_annotations(
        frames,
        gene_symbols,
    )

    result = standardize_apa_sites(result)

    result = result.drop_duplicates(
        subset=[
            "gene_symbol",
            "pas_id",
        ]
    ).reset_index(drop=True)

    result.attrs["source"] = "PolyA_DB v4.1"
    result.attrs["assembly"] = POLYA_DB_ASSEMBLY
    result.attrs["genes_requested"] = gene_symbols

    return result


def standardize_apa_sites(frame):
    """Validate and standardize APA annotations."""

    result = frame.copy()

    result.columns = [
        _clean_column_name(column)
        for column in result.columns
    ]

    missing = (
        REQUIRED_APA_COLUMNS - set(result.columns)
    )

    if missing:

        raise ValueError(
            "Missing required APA columns: "
            + ", ".join(sorted(missing))
        )

    result["sequence_id"] = (
        result["sequence_id"].astype(str)
    )

    result["chrom"] = (
        result["chrom"].astype(str)
    )

    result["start"] = pd.to_numeric(
        result["start"],
        errors="raise",
    )

    result["end"] = pd.to_numeric(
        result["end"],
        errors="raise",
    )

    result["strand"] = (
        result["strand"]
        .astype(str)
        .str.strip()
    )

    if (result["start"] < 0).any():

        raise ValueError(
            "APA start coordinates cannot be negative."
        )

    if (result["end"] <= result["start"]).any():

        raise ValueError(
            "APA end must be greater than start."
        )

    if not result["strand"].isin(
        ["+", "-", "."]
    ).all():

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


def summarize_apa_sites(frame):
    """Summarize annotated APA sites."""

    if frame.empty:

        return {
            "total_sites": 0,
            "sequence_count": 0,
            "proximal_sites": 0,
            "distal_sites": 0,
            "other_sites": 0,
        }

    site_types = (
        frame["site_type"]
        .astype(str)
        .str.lower()
    )

    proximal = site_types.str.contains(
        "proximal",
        na=False,
    ).sum()

    distal = site_types.str.contains(
        "distal",
        na=False,
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
