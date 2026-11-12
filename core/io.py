"""Input parsing and validation utilities for RNA-Insight."""

from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def _get_bytes(source: Any) -> bytes:
    """Read bytes from an uploaded file or filesystem path."""
    if isinstance(source, bytes):
        return source

    if isinstance(source, bytearray):
        return bytes(source)

    if hasattr(source, "getvalue"):
        return source.getvalue()

    if hasattr(source, "read"):
        position = None
        try:
            position = source.tell()
        except Exception:
            pass

        data = source.read()

        if position is not None:
            try:
                source.seek(position)
            except Exception:
                pass

        return data.encode("utf-8") if isinstance(data, str) else data

    path = Path(source)
    return path.read_bytes()


def _get_text(source: Any) -> str:
    """Decode uploaded or local text files, including gzip files."""
    data = _get_bytes(source)

    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)

    return data.decode("utf-8-sig", errors="replace")


def detect_input(source: Any, filename: str = "") -> str:
    """Detect common sequence and annotation file types."""
    text = _get_text(source).lstrip()
    name = filename.lower()

    if text.startswith(">"):
        return "fasta"

    suffixes = name.split(".")

    if "bed" in suffixes:
        return "bed"

    if "gtf" in suffixes:
        return "gtf"

    if "fa" in suffixes or "fasta" in suffixes or "fna" in suffixes:
        return "fasta"

    if "csv" in suffixes:
        return "csv"

    if "tsv" in suffixes or "txt" in suffixes:
        return "tsv"

    return "unknown"


def read_fasta(source: Any) -> pd.DataFrame:
    """Read FASTA records and normalize DNA T bases to RNA U bases."""
    text = _get_text(source)
    records = []

    for record in SeqIO.parse(io.StringIO(text), "fasta"):
        sequence = str(record.seq).upper().replace("T", "U")

        records.append({
            "sequence_id": record.id,
            "description": record.description,
            "sequence": sequence,
            "length": len(sequence),
        })

    if not records:
        raise ValueError("No FASTA records were found.")

    return pd.DataFrame(records)


def read_tabular(source: Any) -> pd.DataFrame:
    """Read CSV or tab-separated annotation tables."""
    text = _get_text(source)

    try:
        frame = pd.read_csv(io.StringIO(text), sep=None, engine="python")
    except Exception:
        try:
            frame = pd.read_csv(io.StringIO(text))
        except Exception as exc:
            raise ValueError(f"Could not read tabular file: {exc}") from exc

    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def read_bed(source: Any) -> pd.DataFrame:
    """Read BED files using 0-based, half-open coordinates."""
    text = _get_text(source)

    frame = pd.read_csv(
        io.StringIO(text),
        sep="\t",
        comment="#",
        header=None,
    )

    if frame.shape[1] < 3:
        raise ValueError("BED input must contain at least 3 columns.")

    names = [
        "chrom", "start", "end", "name", "score", "strand"
    ]

    frame.columns = names[:frame.shape[1]]

    if "name" not in frame.columns:
        frame["name"] = ""

    if "strand" not in frame.columns:
        frame["strand"] = "."

    frame["start"] = pd.to_numeric(frame["start"], errors="raise")
    frame["end"] = pd.to_numeric(frame["end"], errors="raise")

    return frame


def read_gtf(source: Any) -> pd.DataFrame:
    """Read GTF records and convert coordinates to 0-based half-open."""
    text = _get_text(source)
    rows = []

    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.split("\t")

        if len(fields) != 9:
            continue

        attributes = {}

        for item in fields[8].strip().strip(";").split(";"):
            item = item.strip()

            if not item:
                continue

            parts = item.split(" ", 1)

            if len(parts) == 2:
                attributes[parts[0]] = parts[1].strip().strip('"')

        rows.append({
            "chrom": fields[0],
            "source": fields[1],
            "feature": fields[2],
            "start": int(fields[3]) - 1,
            "end": int(fields[4]),
            "score": fields[5],
            "strand": fields[6],
            "frame": fields[7],
            "attributes": attributes,
        })

    if not rows:
        raise ValueError("No valid GTF records were found.")

    return pd.DataFrame(rows)


def validate_rna_records(frame: pd.DataFrame) -> pd.DataFrame:
    """Add sequence validation information to a FASTA DataFrame."""
    if "sequence" not in frame.columns:
        raise ValueError("Input must contain a sequence column.")

    allowed = set("ACGUNRYKMSWBDHV")

    result = frame.copy()
    valid_flags = []
    invalid_characters = []

    for sequence in result["sequence"].fillna("").astype(str):
        sequence = sequence.upper().replace("T", "U")
        invalid = sorted(set(sequence) - allowed)

        valid_flags.append(bool(sequence) and not invalid)
        invalid_characters.append(",".join(invalid))

    result["valid_sequence"] = valid_flags
    result["invalid_characters"] = invalid_characters

    return result
