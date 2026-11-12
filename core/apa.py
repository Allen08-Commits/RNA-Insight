"""Alternative polyadenylation annotation utilities."""

from __future__ import annotations

import pandas as pd


REQUIRED_APA_COLUMNS = {
    "sequence_id",
    "chrom",
    "start",
    "end",
    "strand",
}


def standardize_apa_sites(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate and standardize an APA annotation table."""
    result = frame.copy()

    result.columns = [
        str(column).strip().lower()
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
    result["start"] = pd.to_numeric(result["start"], errors="raise")
    result["end"] = pd.to_numeric(result["end"], errors="raise")
    result["strand"] = result["strand"].astype(str)

    if (result["start"] < 0).any():
        raise ValueError("APA start coordinates cannot be negative.")

    if (result["end"] <= result["start"]).any():
        raise ValueError("APA end must be greater than start.")

    if not result["strand"].isin(["+", "-", "."]).all():
        raise ValueError("Strand values must be +, - or .")

    if "site_type" not in result.columns:
        result["site_type"] = "unspecified"

    result["site_type"] = (
        result["site_type"].fillna("unspecified").astype(str).str.lower()
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

    proximal = site_types.str.contains("proximal", na=False).sum()
    distal = site_types.str.contains("distal", na=False).sum()

    return {
        "total_sites": int(len(frame)),
        "sequence_count": int(frame["sequence_id"].nunique()),
        "proximal_sites": int(proximal),
        "distal_sites": int(distal),
        "other_sites": int(len(frame) - proximal - distal),
    }
