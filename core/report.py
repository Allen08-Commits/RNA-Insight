"""Research report generation for RNA-Insight."""

from __future__ import annotations

from datetime import datetime, timezone

import pandas as pd


def build_report(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    integrated_results: pd.DataFrame | None = None,
) -> str:
    """Build a Markdown report from RNA-Insight analysis results."""
    lines = [
        "# RNA-Insight Research Report",
        "",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "## Study scope",
        "",
        "In-silico screening of RNA G-quadruplex candidate regions "
        "and their genomic overlap with annotated alternative "
        "polyadenylation sites.",
        "",
        "## rG4 candidate screening",
        "",
        f"Total candidate regions: **{len(rg4_candidates)}**",
        "",
    ]

    if not rg4_candidates.empty:
        counts = rg4_candidates["candidate_type"].value_counts()

        for candidate_type, count in counts.items():
            lines.append(f"- {candidate_type}: {count}")

        lines.extend([
            "",
            "### Candidate table",
            "",
            rg4_candidates.head(100).to_markdown(index=False),
            "",
        ])

    lines.extend([
        "## APA annotation summary",
        "",
        f"Total annotated APA sites: **{len(apa_sites)}**",
        "",
    ])

    if not apa_sites.empty:
        if "site_type" in apa_sites.columns:
            counts = apa_sites["site_type"].value_counts()

            for site_type, count in counts.items():
                lines.append(f"- {site_type}: {count}")

        lines.extend([
            "",
            "### APA site table",
            "",
            apa_sites.head(100).to_markdown(index=False),
            "",
        ])

    if integrated_results is not None:
        lines.extend([
            "## Integrated rG4–APA analysis",
            "",
            f"Mapped candidate regions: **{len(integrated_results)}**",
            "",
        ])

        if not integrated_results.empty:
            overlap_count = int(
                integrated_results["overlaps_apa"].sum()
            )

            lines.extend([
                f"Candidate regions overlapping an APA annotation: "
                f"**{overlap_count}**",
                "",
                "### Integrated results",
                "",
                integrated_results.head(100).to_markdown(index=False),
                "",
            ])

    lines.extend([
        "## Interpretation and limitations",
        "",
        "- Sequence motif and composition screening identifies "
        "candidate regions, not experimentally confirmed rG4 structures.",
        "- Genomic overlap does not establish a causal regulatory "
        "relationship between rG4 formation and APA.",
        "- APA annotation counts do not establish differential APA usage.",
        "- Differential APA analysis requires sample-level measurements, "
        "biological replicates and suitable statistical testing.",
        "- Sequence mapping and APA annotations must use compatible "
        "genome assemblies and coordinate conventions.",
        "- Demonstration data supplied with this application are synthetic.",
        "",
        "## Reproducibility",
        "",
        "Record the input file names, genome assembly, annotation source, "
        "analysis parameters and software versions when preparing a "
        "scientific submission.",
        "",
    ])

    return "\n".join(lines)
