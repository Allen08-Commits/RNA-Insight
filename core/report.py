
"""Professional Markdown and PDF report generation for RNA-Insight."""

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
    Table,
    LongTable,
    TableStyle,
    KeepTogether,
)


STUDY_TITLE = (
    "In-Silico Analysis of RNA G-Quadruplex-Mediated "
    "Alternative Polyadenylation in Human Stress-Response Genes"
)


def _safe_int(value) -> int:
    """Convert a value to an integer safely."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _overlap_column(df: pd.DataFrame) -> str | None:
    """Find the overlap column used by different RNA-Insight versions."""
    for column in (
        "overlaps_gene_level_apa",
        "overlaps_apa",
        "overlaps_gene_level_APA",
    ):
        if column in df.columns:
            return column
    return None


def _overlap_count(df: pd.DataFrame) -> int:
    """Count positive overlap records without assuming one column name."""
    column = _overlap_column(df)

    if column is None or df.empty:
        return 0

    values = df[column]

    if values.dtype == object:
        values = (
            values.astype(str)
            .str.strip()
            .str.lower()
            .isin(["true", "1", "yes"])
        )

    return int(pd.Series(values).fillna(False).astype(bool).sum())


def _markdown_table(df: pd.DataFrame, limit: int = 50) -> str:
    """Create a safe Markdown table."""
    if df is None or df.empty:
        return "_No records available._"

    return df.head(limit).to_markdown(index=False)


def build_report(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    integrated_results: pd.DataFrame | None = None,
) -> str:
    """Build a Markdown report for the in-app report preview."""

    rg4_candidates = (
        rg4_candidates
        if rg4_candidates is not None
        else pd.DataFrame()
    )
    apa_sites = (
        apa_sites
        if apa_sites is not None
        else pd.DataFrame()
    )

    lines = [
        "# RNA-Insight Research Report",
        "",
        f"**Study:** {STUDY_TITLE}",
        "",
        f"**Generated:** "
        f"{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "## 1. Study scope",
        "",
        "RNA-Insight is a computational research workspace designed "
        "to identify predicted RNA G-quadruplex candidate regions "
        "and compare their genomic locations with annotated "
        "polyadenylation sites in human stress-response genes.",
        "",
        "## 2. rG4 candidate screening",
        "",
        f"Total candidate intervals: **{len(rg4_candidates)}**",
        "",
    ]

    if not rg4_candidates.empty and "candidate_type" in rg4_candidates:
        counts = rg4_candidates["candidate_type"].value_counts()

        for name, count in counts.items():
            lines.append(f"- {name}: {count}")

    lines.extend([
        "",
        "### Candidate records",
        "",
        _markdown_table(rg4_candidates),
        "",
        "## 3. APA annotation summary",
        "",
        f"Total annotated APA records: **{len(apa_sites)}**",
        "",
    ])

    if not apa_sites.empty and "site_type" in apa_sites.columns:
        counts = apa_sites["site_type"].value_counts()

        for name, count in counts.items():
            lines.append(f"- {name}: {count}")

    lines.extend([
        "",
        "### APA annotation records",
        "",
        _markdown_table(apa_sites),
        "",
    ])

    if integrated_results is not None:
        overlap_count = _overlap_count(integrated_results)

        lines.extend([
            "## 4. Integrated rG4–APA analysis",
            "",
            f"Mapped candidate records: **{len(integrated_results)}**",
            "",
            f"Records flagged as overlapping APA annotations: "
            f"**{overlap_count}**",
            "",
            "### Integrated records",
            "",
            _markdown_table(integrated_results),
            "",
        ])

    lines.extend([
        "## 5. Interpretation",
        "",
        "The sequence-screening workflow identifies regions with "
        "G-rich sequence patterns or canonical G-rich motif patterns "
        "that may be candidates for RNA G-quadruplex formation.",
        "",
        "The APA dataset provides annotated genomic polyadenylation "
        "sites. The site labels in this dataset may include CDS, "
        "5'UTR, 3'UTR, intron and intergenic annotations. These labels "
        "should not automatically be interpreted as proximal or "
        "distal APA usage.",
        "",
        "The integrated analysis compares mapped candidate regions "
        "with the available gene-level APA annotations. A detected "
        "overlap indicates a positional relationship only. It does "
        "not demonstrate that an rG4 structure causes alternative "
        "polyadenylation.",
        "",
        "A zero-overlap result means that no overlap was detected "
        "under the current coordinate mapping and interval comparison "
        "procedure. It does not prove that rG4s have no biological "
        "relationship with APA.",
        "",
        "## 6. Scientific limitations",
        "",
        "- Predicted sequence candidates are not experimentally "
        "confirmed folded RNA G-quadruplex structures.",
        "- G4Hunter-style scores and canonical motif screening are "
        "computational indicators, not experimental validation.",
        "- Gene-level APA annotations are not equivalent to "
        "transcript-specific APA assignments.",
        "- APA annotation counts do not measure differential APA "
        "usage between conditions.",
        "- Differential APA analysis requires sample-level "
        "measurements, biological replicates and suitable statistical "
        "testing.",
        "- Genomic overlap alone cannot establish a causal "
        "rG4-mediated regulatory mechanism.",
        "- Coordinate comparisons depend on compatible genome "
        "assemblies, chromosome naming and coordinate conventions.",
        "",
        "## 7. Recommended future validation",
        "",
        "Future work should include transcript-specific APA annotation, "
        "sample-level 3' end sequencing data, differential APA testing, "
        "and experimental validation of selected rG4 candidates using "
        "appropriate biophysical or molecular techniques.",
        "",
        "## 8. Reproducibility",
        "",
        "For a scientific submission, record the input sequence files, "
        "gene identifiers, Ensembl release, genome assembly, PolyA_DB "
        "version, screening parameters, software versions and date "
        "of retrieval.",
        "",
        "---",
        "",
        "RNA-Insight | Computational research workspace",
    ])

    return "\n".join(lines)


def _paragraph_text(value, max_chars: int = 90) -> str:
    """Escape values before inserting them into PDF Paragraphs."""
    if pd.isna(value):
        return ""

    text = str(value).replace("\n", " ").strip()

    if len(text) > max_chars:
        text = text[:max_chars] + "..."

    return escape(text)


def _make_pdf_table(
    df: pd.DataFrame,
    columns: list[str],
    available_width: float,
    font_size: float = 6.5,
    max_chars: int = 70,
):
    """Build a readable, repeating, multi-page PDF table."""

    if df is None or df.empty:
        return Paragraph(
            "No records available.",
            ParagraphStyle(
                "Empty",
                fontSize=9,
                textColor=colors.grey,
            ),
        )

    selected = [c for c in columns if c in df.columns]

    if not selected:
        selected = list(df.columns[:6])

    data = [[
        Paragraph(
            f"<b>{escape(str(column))}</b>",
            ParagraphStyle(
                "TableHeader",
                fontSize=font_size,
                leading=font_size + 1,
                textColor=colors.white,
            ),
        )
        for column in selected
    ]]

    for _, row in df.iterrows():
        data.append([
            Paragraph(
                _paragraph_text(row[column], max_chars),
                ParagraphStyle(
                    "TableCell",
                    fontSize=font_size,
                    leading=font_size + 1.5,
                    wordWrap="CJK",
                ),
            )
            for column in selected
        ])

    col_width = available_width / len(selected)

    table = LongTable(
        data,
        colWidths=[col_width] * len(selected),
        repeatRows=1,
        hAlign="LEFT",
    )

    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#123456")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#F4F8FC")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [
            colors.white,
            colors.HexColor("#EAF1F8"),
        ]),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#B7C8D8")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))

    return table


def build_pdf_report(
    rg4_candidates: pd.DataFrame,
    apa_sites: pd.DataFrame,
    integrated_results: pd.DataFrame | None = None,
) -> bytes:
    """Generate a professionally formatted PDF and return its bytes."""

    rg4_candidates = (
        rg4_candidates
        if rg4_candidates is not None
        else pd.DataFrame()
    )
    apa_sites = (
        apa_sites
        if apa_sites is not None
        else pd.DataFrame()
    )

    buffer = BytesIO()

    doc = SimpleDocTemplate(
        buffer,
        pagesize=landscape(A4),
        rightMargin=15 * mm,
        leftMargin=15 * mm,
        topMargin=18 * mm,
        bottomMargin=17 * mm,
        title="RNA-Insight Research Report",
        author="RNA-Insight",
        subject=STUDY_TITLE,
    )

    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        name="ReportTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=23,
        leading=29,
        textColor=colors.HexColor("#087F9B"),
        alignment=TA_CENTER,
        spaceAfter=12,
    ))

    styles.add(ParagraphStyle(
        name="StudySubtitle",
        parent=styles["Normal"],
        fontSize=11,
        leading=16,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#263B50"),
        spaceAfter=10,
    ))

    styles.add(ParagraphStyle(
        name="SectionHeading",
        parent=styles["Heading2"],
        fontSize=14,
        leading=18,
        textColor=colors.HexColor("#087F9B"),
        spaceBefore=12,
        spaceAfter=7,
        keepWithNext=True,
    ))

    styles.add(ParagraphStyle(
        name="ReportBody",
        parent=styles["BodyText"],
        fontSize=9,
        leading=13,
        spaceAfter=6,
    ))

    styles.add(ParagraphStyle(
        name="SmallNote",
        parent=styles["BodyText"],
        fontSize=8,
        leading=11,
        textColor=colors.HexColor("#526477"),
        spaceAfter=5,
    ))

    story = []
    page_width = landscape(A4)[0] - 30 * mm

    generated = datetime.now(timezone.utc).strftime(
        "%d %B %Y, %H:%M UTC"
    )

    story.append(Spacer(1, 22 * mm))
    story.append(Paragraph("RNA-Insight", styles["ReportTitle"]))
    story.append(Paragraph(
        escape(STUDY_TITLE),
        styles["StudySubtitle"],
    ))
    story.append(Spacer(1, 8 * mm))
    story.append(Paragraph(
        f"<b>Research report</b><br/>Generated: {generated}",
        ParagraphStyle(
            "CoverInfo",
            parent=styles["Normal"],
            alignment=TA_CENTER,
            fontSize=10,
            leading=16,
        ),
    ))
    story.append(Spacer(1, 15 * mm))

    cover_summary = [
        ["Analysis component", "Records"],
        ["rG4 candidate intervals", str(len(rg4_candidates))],
        ["APA annotation records", str(len(apa_sites))],
    ]

    if integrated_results is not None:
        cover_summary.append([
            "Integrated candidate records",
            str(len(integrated_results)),
        ])
        cover_summary.append([
            "Records flagged as APA overlaps",
            str(_overlap_count(integrated_results)),
        ])

    cover_table = Table(
        cover_summary,
        colWidths=[100 * mm, 45 * mm],
        hAlign="CENTER",
    )
    cover_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#123456")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#EAF1F8")),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#B7C8D8")),
        ("ALIGN", (1, 0), (1, -1), "CENTER"),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 9),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    story.append(cover_table)
    story.append(Spacer(1, 10 * mm))
    story.append(Paragraph(
        "Prepared by the RNA-Insight computational research workspace.",
        styles["SmallNote"],
    ))

    story.append(Paragraph("1. Study scope", styles["SectionHeading"]))
    story.append(Paragraph(
        "This study uses computational screening to identify predicted "
        "RNA G-quadruplex candidate regions and compares their mapped "
        "genomic locations with annotated polyadenylation sites in "
        "human stress-response genes.",
        styles["ReportBody"],
    ))

    story.append(Paragraph(
        "2. rG4 candidate screening",
        styles["SectionHeading"],
    ))
    story.append(Paragraph(
        f"Total candidate intervals: <b>{len(rg4_candidates)}</b>",
        styles["ReportBody"],
    ))

    if not rg4_candidates.empty and "candidate_type" in rg4_candidates:
        for name, count in rg4_candidates["candidate_type"].value_counts().items():
            story.append(Paragraph(
                f"• {escape(str(name))}: {int(count)}",
                styles["ReportBody"],
            ))

    story.append(Spacer(1, 3 * mm))
    story.append(_make_pdf_table(
        rg4_candidates,
        [
            "sequence_id",
            "start",
            "end",
            "candidate_type",
            "score",
            "score_type",
            "length",
            "candidate_sequence",
        ],
        page_width,
        font_size=6,
        max_chars=42,
    ))

    story.append(Paragraph(
        "3. APA annotation summary",
        styles["SectionHeading"],
    ))
    story.append(Paragraph(
        f"Total annotated APA records: <b>{len(apa_sites)}</b>",
        styles["ReportBody"],
    ))

    if not apa_sites.empty and "site_type" in apa_sites.columns:
        for name, count in apa_sites["site_type"].value_counts().items():
            story.append(Paragraph(
                f"• {escape(str(name))}: {int(count)}",
                styles["ReportBody"],
            ))

    story.append(Spacer(1, 3 * mm))
    story.append(_make_pdf_table(
        apa_sites,
        [
            "gene_symbol",
            "sequence_id",
            "chrom",
            "start",
            "end",
            "strand",
            "site_type",
            "assembly",
            "source",
        ],
        page_width,
        font_size=6,
        max_chars=38,
    ))

    if integrated_results is not None:
        story.append(Paragraph(
            "4. Integrated rG4–APA analysis",
            styles["SectionHeading"],
        ))

        story.append(Paragraph(
            f"Mapped candidate records: <b>{len(integrated_results)}</b><br/>"
            f"Records flagged as overlapping APA annotations: "
            f"<b>{_overlap_count(integrated_results)}</b>",
            styles["ReportBody"],
        ))

        story.append(_make_pdf_table(
            integrated_results,
            [
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
                "overlapping_apa_sites",
            ],
            page_width,
            font_size=5.8,
            max_chars=38,
        ))

    story.append(Paragraph(
        "5. Interpretation",
        styles["SectionHeading"],
    ))
    story.append(Paragraph(
        "The sequence-screening workflow identifies regions with "
        "G-rich sequence patterns or canonical G-rich motif patterns "
        "that may be candidates for RNA G-quadruplex formation. "
        "These predictions are not experimental confirmation of "
        "folded RNA G-quadruplex structures.",
        styles["ReportBody"],
    ))
    story.append(Paragraph(
        "The APA dataset contains annotated genomic polyadenylation "
        "sites. Labels such as CDS, 5'UTR, 3'UTR, intron and intergenic "
        "describe annotation categories and should not automatically "
        "be interpreted as proximal or distal APA usage.",
        styles["ReportBody"],
    ))

    if integrated_results is not None:
        overlap_count = _overlap_count(integrated_results)

        if overlap_count == 0:
            interpretation = (
                "No exact overlaps were detected by the current "
                "integration procedure. This result does not prove "
                "that rG4s have no biological relationship with APA."
            )
        else:
            interpretation = (
                f"The current analysis flagged {overlap_count} "
                "candidate records as overlapping APA annotations. "
                "These are positional associations and do not establish "
                "a causal regulatory relationship."
            )

        story.append(Paragraph(
            interpretation,
            styles["ReportBody"],
        ))

    story.append(Paragraph(
        "6. Scientific limitations",
        styles["SectionHeading"],
    ))

    limitations = [
        "Predicted sequence candidates are not experimentally confirmed "
        "RNA G-quadruplex structures.",
        "Gene-level APA annotations are not equivalent to "
        "transcript-specific APA assignments.",
        "APA annotation counts do not measure differential APA usage.",
        "Differential APA analysis requires sample-level measurements, "
        "biological replicates and suitable statistical testing.",
        "Genomic overlap alone cannot establish a causal "
        "rG4-mediated regulatory mechanism.",
        "Coordinate comparisons depend on compatible genome assemblies "
        "and coordinate conventions.",
    ]

    for item in limitations:
        story.append(Paragraph(
            f"• {escape(item)}",
            styles["ReportBody"],
        ))

    story.append(Paragraph(
        "7. Recommended future validation",
        styles["SectionHeading"],
    ))
    story.append(Paragraph(
        "Future work should include transcript-specific APA annotation, "
        "sample-level 3' end sequencing data, differential APA testing, "
        "and experimental validation of selected rG4 candidates using "
        "appropriate biophysical or molecular techniques.",
        styles["ReportBody"],
    ))

    story.append(Paragraph(
        "8. Reproducibility",
        styles["SectionHeading"],
    ))
    story.append(Paragraph(
        "For scientific submission, record input sequence files, gene "
        "identifiers, Ensembl release, genome assembly, PolyA_DB version, "
        "screening parameters, software versions and retrieval date.",
        styles["ReportBody"],
    ))

    def add_page_number(canvas, document):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#526477"))
        canvas.drawString(
            15 * mm,
            9 * mm,
            "RNA-Insight | Computational research report",
        )
        canvas.drawRightString(
            landscape(A4)[0] - 15 * mm,
            9 * mm,
            f"Page {document.page}",
        )
        canvas.restoreState()

    doc.build(
        story,
        onFirstPage=add_page_number,
        onLaterPages=add_page_number,
    )

    pdf_bytes = buffer.getvalue()
    buffer.close()

    return pdf_bytes
