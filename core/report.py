
"""Research report generation and PDF export for RNA-Insight."""

from __future__ import annotations

from datetime import datetime, timezone
from io import BytesIO
from xml.sax.saxutils import escape
import re

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
    Table,
    TableStyle,
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
                    integrated_results[overlap_column]
                    .fillna(False)
                    .astype(bool)
                    .sum()
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
    """Convert RNA-Insight Markdown report into a formatted PDF."""

    buffer = BytesIO()

    document = SimpleDocTemplate(
        buffer,
        pagesize=landscape(A4),
        rightMargin=12 * mm,
        leftMargin=12 * mm,
        topMargin=15 * mm,
        bottomMargin=15 * mm,
        title="RNA-Insight Research Report",
        author="RNA-Insight",
    )

    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        name="RNAReportTitle",
        parent=styles["Title"],
        fontSize=22,
        leading=27,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#087F9B"),
        spaceAfter=16,
    ))

    styles.add(ParagraphStyle(
        name="RNAReportHeading",
        parent=styles["Heading1"],
        fontSize=16,
        leading=20,
        textColor=colors.HexColor("#087F9B"),
        spaceBefore=12,
        spaceAfter=7,
    ))

    styles.add(ParagraphStyle(
        name="RNAReportSubheading",
        parent=styles["Heading2"],
        fontSize=12,
        leading=15,
        textColor=colors.HexColor("#6045A5"),
        spaceBefore=8,
        spaceAfter=5,
    ))

    styles.add(ParagraphStyle(
        name="RNAReportBody",
        parent=styles["BodyText"],
        fontSize=9,
        leading=13,
        spaceAfter=6,
    ))

    styles.add(ParagraphStyle(
        name="RNAReportSmall",
        parent=styles["BodyText"],
        fontSize=6,
        leading=8,
        wordWrap="CJK",
    ))

    story = []
    table_lines = []

    def clean_markdown(value):
        value = re.sub(r"\*\*(.*?)\*\*", r"\1", value)
        value = re.sub(r"\*(.*?)\*", r"\1", value)
        value = re.sub(r"`(.*?)`", r"\1", value)
        return value.strip()

    def flush_table():
        nonlocal table_lines

        if not table_lines:
            return

        rows = []

        for line in table_lines:
            cells = [
                cell.strip()
                for cell in line.strip().strip("|").split("|")
            ]

            # Skip Markdown separator rows
            if all(
                re.fullmatch(r":?-{3,}:?", cell.replace(" ", ""))
                for cell in cells
            ):
                continue

            rows.append(cells)

        if rows:
            max_columns = max(len(row) for row in rows)

            normalized_rows = [
                row + [""] * (max_columns - len(row))
                for row in rows
            ]

            pdf_rows = []

            for row_index, row in enumerate(normalized_rows):
                formatted_row = []

                for cell in row:
                    formatted_row.append(
                        Paragraph(
                            escape(clean_markdown(cell)),
                            styles["RNAReportSmall"],
                        )
                    )

                pdf_rows.append(formatted_row)

            available_width = landscape(A4)[0] - 24 * mm
            col_width = available_width / max_columns

            table = Table(
                pdf_rows,
                colWidths=[col_width] * max_columns,
                repeatRows=1,
                hAlign="LEFT",
            )

            table.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0),
                 colors.HexColor("#DCEEF5")),
                ("TEXTCOLOR", (0, 0), (-1, 0),
                 colors.HexColor("#123A55")),
                ("GRID", (0, 0), (-1, -1), 0.35,
                 colors.HexColor("#B8C9D6")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]))

            story.append(Spacer(1, 5))
            story.append(table)
            story.append(Spacer(1, 8))

        table_lines = []

    for raw_line in report_markdown.splitlines():
        line = raw_line.strip()

        if line.startswith("|"):
            table_lines.append(line)
            continue

        flush_table()

        if not line:
            story.append(Spacer(1, 4))
            continue

        if line == "---":
            story.append(Spacer(1, 8))
            continue

        if line.startswith("# "):
            story.append(Paragraph(
                escape(clean_markdown(line[2:])),
                styles["RNAReportTitle"],
            ))

        elif line.startswith("## "):
            story.append(Paragraph(
                escape(clean_markdown(line[3:])),
                styles["RNAReportHeading"],
            ))

        elif line.startswith("### "):
            story.append(Paragraph(
                escape(clean_markdown(line[4:])),
                styles["RNAReportSubheading"],
            ))

        elif line.startswith("- "):
            story.append(Paragraph(
                "&#8226; " + escape(clean_markdown(line[2:])),
                styles["RNAReportBody"],
            ))

        else:
            story.append(Paragraph(
                escape(clean_markdown(line)),
                styles["RNAReportBody"],
            ))

    flush_table()

    def add_page_number(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.grey)

        canvas.drawRightString(
            landscape(A4)[0] - 12 * mm,
            8 * mm,
            f"RNA-Insight | Page {doc.page}",
        )

        canvas.restoreState()

    document.build(
        story,
        onFirstPage=add_page_number,
        onLaterPages=add_page_number,
    )

    return buffer.getvalue()
