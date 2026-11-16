
"""
RNA-Insight Research Report Generator

Generates a structured Markdown report from:
1. Predicted RNA G-quadruplex candidates
2. Polyadenylation site annotations
3. Integrated gene-level genomic comparisons

Scientific note:
Gene-level genomic overlap is not equivalent to
transcript-specific APA usage or causal regulation.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

import pandas as pd


# --------------------------------------------------
# SAFE DATAFRAME HELPERS
# --------------------------------------------------

def _safe_frame(data: Any) -> pd.DataFrame:
    """Return a safe copy of a DataFrame."""

    if isinstance(data, pd.DataFrame):
        return data.copy()

    if data is None:
        return pd.DataFrame()

    try:
        return pd.DataFrame(data)
    except Exception:
        return pd.DataFrame()


def _count_boolean(series: pd.Series) -> int:
    """Count true-like values safely."""

    if series is None:
        return 0

    if pd.api.types.is_bool_dtype(series):
        return int(series.fillna(False).sum())

    values = (
        series.astype(str)
        .str.strip()
        .str.lower()
    )

    return int(values.isin(["true", "1", "yes", "y"]).sum())


def _get_overlap_count(frame: pd.DataFrame) -> int:
    """
    Support both the older and current integration schemas.

    Current RNA-Insight output:
        overlaps_gene_level_apa

    Older output:
        overlaps_apa
    """

    if frame.empty:
        return 0

    possible_columns = [
        "overlaps_gene_level_apa",
        "overlaps_apa",
    ]

    for column in possible_columns:
        if column in frame.columns:
            return _count_boolean(frame[column])

    return 0


def _unique_count(
    frame: pd.DataFrame,
    column: str,
) -> int:
    """Count unique non-empty values."""

    if frame.empty or column not in frame.columns:
        return 0

    values = (
        frame[column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    values = values[
        (values != "")
        & (values.str.lower() != "nan")
    ]

    return int(values.nunique())


def _markdown_table(frame: pd.DataFrame, limit: int = 15) -> str:
    """Convert a DataFrame into a safe Markdown table."""

    if frame.empty:
        return "_No records available._"

    preview = frame.head(limit).copy()

    return preview.to_markdown(index=False)


# --------------------------------------------------
# REPORT GENERATION
# --------------------------------------------------

def build_report(
    rg4_candidates: pd.DataFrame | None = None,
    apa_sites: pd.DataFrame | None = None,
    integrated_results: pd.DataFrame | None = None,
) -> str:
    """
    Build the complete RNA-Insight research report.

    Parameters
    ----------
    rg4_candidates:
        Predicted RNA G4 candidate intervals.

    apa_sites:
        Annotated polyadenylation site records.

    integrated_results:
        Mapped rG4 candidates and gene-level APA comparisons.

    Returns
    -------
    str
        Markdown research report.
    """

    rg4 = _safe_frame(rg4_candidates)
    apa = _safe_frame(apa_sites)
    integrated = _safe_frame(integrated_results)

    # --------------------------------------------------
    # BASIC STATISTICS
    # --------------------------------------------------

    total_rg4 = len(rg4)
    total_apa = len(apa)
    total_integrated = len(integrated)

    overlap_count = _get_overlap_count(integrated)

    mapped_genes = _unique_count(
        integrated,
        "gene_symbol",
    )

    if mapped_genes == 0:
        mapped_genes = _unique_count(
            integrated,
            "sequence_id",
        )

    mapped_transcripts = _unique_count(
        integrated,
        "sequence_id",
    )

    # --------------------------------------------------
    # rG4 CANDIDATE CLASSIFICATION
    # --------------------------------------------------

    candidate_counts = {}

    if "candidate_type" in rg4.columns:
        candidate_counts = (
            rg4["candidate_type"]
            .fillna("Unspecified")
            .astype(str)
            .value_counts()
            .to_dict()
        )

    candidate_summary = "\n".join(
        f"- {name}: {count}"
        for name, count in candidate_counts.items()
    )

    if not candidate_summary:
        candidate_summary = "- Candidate classification unavailable."

    # --------------------------------------------------
    # rG4 SCORE SUMMARY
    # --------------------------------------------------

    score_summary = "_No numerical score information available._"

    if "score" in rg4.columns:
        scores = pd.to_numeric(
            rg4["score"],
            errors="coerce",
        ).dropna()

        if not scores.empty:
            score_summary = (
                f"- Mean candidate score: {scores.mean():.4f}\n"
                f"- Minimum score: {scores.min():.4f}\n"
                f"- Maximum score: {scores.max():.4f}"
            )

    # --------------------------------------------------
    # APA ANNOTATION SUMMARY
    # --------------------------------------------------

    apa_summary = "_No APA annotation records available._"

    if not apa.empty:
        if "site_type" in apa.columns:
            site_counts = (
                apa["site_type"]
                .fillna("Unspecified")
                .astype(str)
                .value_counts()
            )

            apa_summary = "\n".join(
                f"- {site_type}: {count}"
                for site_type, count in site_counts.items()
            )

        else:
            apa_summary = (
                f"- Total annotated PAS records: {total_apa}"
            )

    apa_genes = _unique_count(
        apa,
        "gene_symbol",
    )

    if apa_genes == 0:
        apa_genes = _unique_count(
            apa,
            "sequence_id",
        )

    # --------------------------------------------------
    # INTEGRATION SUMMARY
    # --------------------------------------------------

    if total_integrated > 0:

        integration_status = (
            "The integration table contains mapped candidate records."
        )

        if overlap_count > 0:
            overlap_interpretation = (
                f"{overlap_count} mapped candidate records were "
                "labelled as overlapping annotated gene-level PAS intervals."
            )
        else:
            overlap_interpretation = (
                "No direct gene-level APA interval overlaps were "
                "identified under the current comparison criteria."
            )

    else:

        integration_status = (
            "Integrated results were not available when this report "
            "was generated."
        )

        overlap_interpretation = (
            "No overlap conclusion can be drawn because integrated "
            "records are unavailable."
        )

    # --------------------------------------------------
    # INTEGRATED TABLE PREVIEW
    # --------------------------------------------------

    integration_columns = [
        column
        for column in [
            "sequence_id",
            "gene_symbol",
            "candidate_start",
            "candidate_end",
            "chrom",
            "strand",
            "candidate_type",
            "score",
            "genomic_segments",
            "overlaps_gene_level_apa",
            "overlaps_apa",
        ]
        if column in integrated.columns
    ]

    if integration_columns:
        integration_preview = _markdown_table(
            integrated[integration_columns],
            limit=15,
        )
    else:
        integration_preview = (
            "_No compatible integrated result columns available._"
        )

    # --------------------------------------------------
    # APA TABLE PREVIEW
    # --------------------------------------------------

    apa_columns = [
        column
        for column in [
            "gene_symbol",
            "sequence_id",
            "chrom",
            "start",
            "end",
            "strand",
            "site_type",
            "pas_id",
        ]
        if column in apa.columns
    ]

    if apa_columns:
        apa_preview = _markdown_table(
            apa[apa_columns],
            limit=15,
        )
    else:
        apa_preview = "_No APA annotation preview available._"

    # --------------------------------------------------
    # rG4 TABLE PREVIEW
    # --------------------------------------------------

    rg4_columns = [
        column
        for column in [
            "sequence_id",
            "start",
            "end",
            "candidate_type",
            "score",
            "score_type",
            "length",
        ]
        if column in rg4.columns
    ]

    if rg4_columns:
        rg4_preview = _markdown_table(
            rg4[rg4_columns],
            limit=15,
        )
    else:
        rg4_preview = "_No rG4 candidate preview available._"

    # --------------------------------------------------
    # SCIENTIFIC INTERPRETATION
    # --------------------------------------------------

    interpretation = f"""
The computational analysis identified **{total_rg4} predicted
rG4 candidate intervals** and **{total_apa} annotated
polyadenylation site records**.

A total of **{total_integrated} mapped rG4 candidate records**
were included in the integrated comparison.

{overlap_interpretation}

These findings describe computational sequence predictions
and genomic positional comparisons. They do not establish
that an RNA G-quadruplex forms experimentally or that it
causes alternative polyadenylation.
"""

    # --------------------------------------------------
    # COMPLETE REPORT
    # --------------------------------------------------

    generated_at = datetime.now(
        timezone.utc
    ).strftime("%d %B %Y, %H:%M UTC")

    report = f"""# RNA-Insight Research Report

## In-Silico Analysis of RNA G-Quadruplex-Mediated
## Alternative Polyadenylation in Human Stress-Response Genes

**Platform:** RNA-Insight  
**Report generated:** {generated_at}  
**Analysis type:** Computational bioinformatics research

---

# 1. Abstract

RNA-Insight is a computational research workspace designed
to investigate predicted RNA G-quadruplex (rG4) candidate
regions and their genomic relationship with annotated
polyadenylation sites in human stress-response genes.

The workflow combines human transcript retrieval,
RNA sequence screening, polyadenylation site annotation,
transcript-to-genome exon mapping and gene-level genomic
comparison.

The analysis identifies candidate sequence regions and
examines their positional relationship with annotated
polyadenylation sites. The results provide computational
evidence for further investigation rather than experimental
confirmation of an rG4-mediated regulatory mechanism.

---

# 2. Research Objectives

1. Retrieve human stress-response gene and transcript
   information from Ensembl.
2. Obtain transcript cDNA sequences for computational analysis.
3. Identify candidate RNA G-quadruplex-forming regions.
4. Retrieve annotated polyadenylation sites from PolyA_DB.
5. Map transcript coordinates to genomic coordinates using
   exon-aware Ensembl mapping.
6. Compare mapped rG4 candidates with gene-level
   polyadenylation annotations.
7. Generate interpretable results and downloadable reports.

---

# 3. Methodology

## 3.1 Human Gene and Transcript Retrieval

Human gene information and transcript sequences are
retrieved using the Ensembl REST API.

The selected research panel includes stress-response genes
such as HSPA1A, HSP90AA1, DDIT3 and ATF4.

Transcript cDNA sequences are used as the input for
sequence-level rG4 screening.

## 3.2 RNA G-Quadruplex Candidate Screening

Retrieved RNA sequences are screened for potential
G-quadruplex-forming regions.

The workflow identifies:

- Canonical G-rich sequence motif candidates.
- G4Hunter-labelled candidate regions.
- Candidate coordinates within transcript sequences.
- Available numerical candidate scores.

These outputs represent computational predictions.
They do not confirm the physical formation of folded
RNA G-quadruplex structures.

## 3.3 Alternative Polyadenylation Annotation

Polyadenylation site annotations are retrieved from
PolyA_DB v4.1 using the hg38 reference assembly.

The retrieved annotations provide genomic positions,
chromosome information, strand orientation and
site-type labels where available.

The current PolyA_DB records are gene-level annotations.
They must not automatically be interpreted as
transcript-specific APA assignments.

## 3.4 Transcript-to-Genome Mapping

Transcript cDNA coordinates are mapped to genomic
coordinates using Ensembl exon annotations.

The mapping procedure:

- Retrieves exon coordinates.
- Preserves chromosome and strand information.
- Converts genomic coordinates into the required
  coordinate convention.
- Accounts for candidates spanning exon junctions.
- Validates the summed exon length against cDNA length.

Mappings are accepted only when their exon-length
validation succeeds.

## 3.5 Integrated Genomic Comparison

Mapped rG4 candidate intervals are compared with
annot
