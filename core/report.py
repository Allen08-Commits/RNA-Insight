
"""Research report generation for RNA-Insight."""

from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd


def _safe_dataframe(data):
    """Convert input into a safe DataFrame."""
    if isinstance(data, pd.DataFrame):
        return data.copy()

    if data is None:
        return pd.DataFrame()

    return pd.DataFrame(data)


def _count_true(series):
    """Count true values safely."""
    if series.empty:
        return 0

    if pd.api.types.is_bool_dtype(series):
        return int(series.fillna(False).sum())

    values = series.astype(str).str.strip().str.lower()
    return int(values.isin(["true", "1", "yes", "y"]).sum())


def _get_overlap_count(frame):
    """Read overlap count from either supported column name."""

    if frame.empty:
        return 0

    if "overlaps_gene_level_apa" in frame.columns:
        return _count_true(frame["overlaps_gene_level_apa"])

    if "overlaps_apa" in frame.columns:
        return _count_true(frame["overlaps_apa"])

    return 0


def _table_preview(frame, limit=100):
    """Create a Markdown table preview."""

    if frame.empty:
        return "_No records available._"

    return frame.head(limit).to_markdown(index=False)


def build_report(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    integrated_results: pd.DataFrame | None = None,
) -> str:
    """Build a complete Markdown report from RNA-Insight results."""

    rg4 = _safe_dataframe(rg4_candidates)
    apa = _safe_dataframe(apa_sites)
    integrated = _safe_dataframe(integrated_results)

    generated = datetime.now(timezone.utc).strftime(
        "%Y-%m-%d %H:%M UTC"
    )

    lines = [
        "# RNA-Insight Research Report",
        "",
        "**Study:** In-Silico Analysis of RNA G-Quadruplex-Mediated "
        "Alternative Polyadenylation in Human Stress-Response Genes",
        "",
        f"**Generated:** {generated}",
        "",
        "---",
        "",
        "## 1. Study scope",
        "",
        "RNA-Insight is a computational research workspace for "
        "identifying predicted RNA G-quadruplex candidate regions "
        "and comparing their genomic locations with annotated "
        "polyadenylation sites.",
        "",
        "The workflow combines human transcript retrieval, "
        "RNA sequence screening, APA annotation retrieval, "
        "exon-aware transcript-to-genome mapping and genomic "
        "comparison.",
        "",
        "## 2. rG4 candidate screening",
        "",
        f"Total predicted candidate intervals: **{len(rg4)}**",
        "",
    ]

    if not rg4.empty:

        if "candidate_type" in rg4.columns:
            counts = rg4["candidate_type"].value_counts()

            lines.append("### Candidate classification")
            lines.append("")

            for candidate_type, count in counts.items():
                lines.append(f"- {candidate_type}: {count}")

            lines.append("")

        if "score" in rg4.columns:
            scores = pd.to_numeric(
                rg4["score"],
                errors="coerce",
            ).dropna()

            if not scores.empty:
                lines.extend([
                    "### Candidate score summary",
                    "",
                    f"- Mean score: {scores.mean():.4f}",
                    f"- Minimum score: {scores.min():.4f}",
                    f"- Maximum score: {scores.max():.4f}",
                    "",
                ])

        lines.extend([
            "### rG4 candidate records",
            "",
            _table_preview(rg4),
            "",
        ])

    else:
        lines.extend([
            "No rG4 candidate records were available.",
            "",
        ])

    lines.extend([
        "## 3. Alternative polyadenylation annotation",
        "",
        f"Total annotated PAS records: **{len(apa)}**",
        "",
    ])

    if not apa.empty:

        if "site_type" in apa.columns:
            counts = apa["site_type"].value_counts()

            lines.append("### APA annotation categories")
            lines.append("")

            for site_type, count in counts.items():
                lines.append(f"- {site_type}: {count}")

            lines.append("")

        lines.extend([
            "### APA annotation records",
            "",
            _table_preview(apa),
            "",
        ])

    else:
        lines.extend([
            "No APA annotation records were available.",
            "",
        ])

    # Integrated results
    lines.extend([
        "## 4. Integrated rG4–APA analysis",
        "",
    ])

    if integrated_results is None:

        lines.extend([
            "Integrated analysis was not available when this report "
            "was generated.",
            "",
        ])

    else:

        overlap_count = _get_overlap_count(integrated)

        lines.extend([
            f"Mapped rG4 candidate records: **{len(integrated)}**",
            "",
            f"Gene-level APA-overlapping records: **{overlap_count}**",
            "",
        ])

        if not integrated.empty:

            if "gene_symbol" in integrated.columns:
                genes = integrated["gene_symbol"].dropna().nunique()

                lines.extend([
                    f"Genes represented: **{genes}**",
                    "",
                ])

            lines.extend([
                "### Integrated result records",
                "",
                _table_preview(integrated),
                "",
            ])

        else:
            lines.extend([
                "No integrated candidate records were available.",
                "",
            ])

        if overlap_count == 0 and not integrated.empty:
            lines.extend([
                "No direct gene-level APA interval overlaps were "
                "identified under the current comparison criteria.",
                "",
            ])

    # Interpretation
    lines.extend([
        "## 5. Interpretation",
        "",
        f"The analysis identified **{len(rg4)} predicted rG4 "
        "candidate intervals** and **"
        f"{len(apa)} annotated polyadenylation site records**.",
        "",
    ])

    if integrated_results is not None:
        lines.extend([
            f"A total of **{len(integrated)} mapped candidate records** "
            "were included in the genomic comparison.",
            "",
            f"The comparison identified **{_get_overlap_count(integrated)} "
            "gene-level APA-overlapping records**.",
            "",
        ])

    lines.extend([
        "These results describe computational predictions and "
        "genomic positional comparisons. They do not independently "
        "demonstrate that an RNA G-quadruplex forms experimentally "
        "or regulates alternative polyadenylation.",
        "",
        "## 6. Limitations",
        "",
        "- Predicted rG4 candidates are not experimentally confirmed "
        "RNA G-quadruplex structures.",
        "",
        "- Gene-level PAS annotations should not be interpreted as "
        "transcript-specific APA assignments.",
        "",
        "- Genomic overlap or proximity does not establish a causal "
        "relationship between rG4 formation and APA.",
        "",
        "- APA annotation counts do not measure differential APA usage.",
        "",
        "- Differential APA analysis requires sample-level measurements, "
        "biological replicates and appropriate statistical testing.",
        "",
        "- Genomic comparisons depend on compatible genome assemblies "
        "and coordinate conventions.",
        "",
        "## 7. Reproducibility",
        "",
        "For reproducibility, retain the input sequences, annotation "
        "source and version, genome assembly, analysis parameters "
        "and software versions.",
        "",
        "## 8. Conclusion",
        "",
        "RNA-Insight provides an integrated computational workflow "
        "for screening RNA G-quadruplex candidates, retrieving "
        "polyadenylation annotations and comparing genomic positions "
        "using exon-aware transcript mapping.",
        "",
        "The findings provide a basis for further investigation "
        "of potential relationships between RNA G-quadruplexes "
        "and alternative polyadenylation in human stress-response genes. "
        "Experimental validation and transcript-level APA measurements "
        "would be required to establish a regulatory mechanism.",
        "",
        "---",
        "",
        "*End of RNA-Insight research report.*",
    ])

    return "\n".join(lines)
