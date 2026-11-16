
"""Research report generation and PDF export for RNA-Insight."""

from __future__ import annotations

from datetime import datetime, timezone
from io import BytesIO
from xml.sax.saxutils import escape

import pandas as pd

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    PageBreak,
)


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
            overlap_column = (
                "overlaps_gene_level_apa"
                if "overlaps_gene_level_apa" in integrated_results.columns
                else "overlaps_apa"
            )

            if overlap_column in integrated_results.columns:
                overlap_count = int(
                    integrated_results[overlap_column].fillna(False).sum()
                )
            else:
                overlap_count = 0

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


def build_pdf_report(report_markdown: str) -> bytes:
    """Convert the generated Markdown report into a downloadable PDF."""

    buffer = BytesIO()

    document = SimpleDocTemplate(
        buffer,
        pagesize=landscape(A4),
        rightMargin=15 * mm,
        leftMargin=15 * mm,
        topMargin=17 * mm,
        bottomMargin=17 * mm,
        title="RNA-Insight Research Report",
        author="RNA-Insight",
    )

    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        name="ReportTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=22,
        leading=27,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#087E9B"),
        spaceAfter=14,
    ))

    styles.add(ParagraphStyle(
        name="ReportHeading",
        parent=styles["Heading1"],
        fontSize=16,
        leading=20,
        textColor=colors.HexColor("#087E9B"),
        spaceBefore=13,
        spaceAfter=8,
    ))

    styles.add(ParagraphStyle(
        name="ReportSubheading",
        parent=styles["Heading2"],
        fontSize=12,
        leading=15,
        textColor=colors.HexColor("#563B91"),
        spaceBefore=9,
        spaceAfter=5,
    ))

    styles.add(ParagraphStyle(
        name="ReportBody",
        parent=styles["BodyText"],
        fontSize=9,
        leading=13,
        spaceAfter=6,
    ))

    styles.add(ParagraphStyle(
        name="ReportBullet",
        parent=styles["BodyText"],
        fontSize=9,
        leading=12,
        leftIndent=14,
        firstLineIndent=-8,
        spaceAfter=4,
    ))

    styles.add(ParagraphStyle(
        name="ReportTableLine",
        parent=styles["Code"],
        fontName="Courier",
        fontSize=5.5,
        leading=7,
        spaceAfter=2,
        wordWrap="CJK",
    ))

    story = []

    for raw_line in report_markdown.splitlines():
        line = raw_line.strip()

        if not line:
            story.append(Spacer(1, 3))
            continue

        if line == "---":
            story.append(Spacer(1, 5))
            continue

        if line.startswith("# "):
            story.append(Paragraph(
                escape(line[2:]),
                styles["ReportTitle"]
            ))

        elif line.startswith("## "):
            story.append(Paragraph(
                escape(line[3:]),
                styles["ReportHeading"]
            ))

        elif line.startswith("### "):
            story.append(Paragraph(
                escape(line[4:]),
                styles["ReportSubheading"]
            ))

        elif line.startswith("- "):
            content = escape(line[2:])
            story.append(Paragraph(
                "• " + content,
                styles["ReportBullet"]
            ))

        elif line.startswith("|"):
            # Keep wide scientific table rows readable and wrapped.
            table_text = escape(line.replace("|", " | "))
            story.append(Paragraph(
                table_text,
                styles["ReportTableLine"]
            ))

        else:
            # Basic Markdown bold formatting.
            safe_line = escape(line)
            safe_line = safe_line.replace("**", "")
            story.append(Paragraph(
                safe_line,
                styles["ReportBody"]
            ))

    def add_page_number(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.grey)
        canvas.drawRightString(
            landscape(A4)[0] - 15 * mm,
            9 * mm,
            f"RNA-Insight | Page {doc.page}",
        )
        canvas.restoreState()

    document.build(
        story,
        onFirstPage=add_page_number,
        onLaterPages=add_page_number,
    )

    return buffer.getvalue()
