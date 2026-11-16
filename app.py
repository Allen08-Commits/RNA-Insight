
"""RNA-Insight: Human stress-response gene and rG4-APA research."""

from __future__ import annotations

import re

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from core.io import (
    read_fasta,
    read_tabular,
    read_bed,
    read_gtf,
    detect_input,
    validate_rna_records,
)

from core.rg4 import scan_rg4
from core.apa import standardize_apa_sites, summarize_apa_sites
from core.integration import integrate_rg4_apa
from core.report import build_report
from core.gene_data import DEFAULT_GENES, fetch_gene_panel


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="RNA-Insight",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ============================================================
# CUSTOM INTERFACE
# ============================================================

st.markdown("""
<style>
.stApp {
    background: linear-gradient(
        135deg, #07111f 0%, #0b1730 55%, #11152e 100%
    );
    color: #e8f2ff;
}

section[data-testid="stSidebar"] {
    background: #091426;
    border-right: 1px solid #1e4161;
}

h1, h2, h3 {
    color: #70e5f5 !important;
}

div[data-testid="stMetric"] {
    background: rgba(21, 45, 75, 0.72);
    border: 1px solid #245078;
    padding: 16px;
    border-radius: 12px;
}

div.stButton > button {
    background: linear-gradient(90deg, #087f9b, #6948ad);
    color: white;
    border: 0;
    border-radius: 9px;
    font-weight: 600;
}

div.stButton > button:hover {
    border: 1px solid #7cecff;
    color: white;
}

div[data-testid="stDownloadButton"] button {
    border-radius: 9px;
}
</style>
""", unsafe_allow_html=True)


# ============================================================
# GENOMIC VISUALIZATION
# ============================================================

def create_genomic_figure(
    integrated: pd.DataFrame,
    apa_sites: pd.DataFrame,
    chromosome: str,
) -> go.Figure:
    """Visualize exon-mapped rG4 segments and APA annotations."""

    fig = go.Figure()

    candidates = integrated[
        integrated["chrom"].astype(str) == str(chromosome)
    ].copy()

    apa = apa_sites[
        apa_sites["chrom"].astype(str) == str(chromosome)
    ].copy()

    for row in candidates.itertuples():

        segments_text = getattr(row, "genomic_segments", "")

        segments = re.findall(
            r"([^:;]+):(\d+)-(\d+)\(([+-])\)",
            str(segments_text),
        )

        for chrom, start, end, strand in segments:

            fig.add_trace(
                go.Scatter(
                    x=[int(start), int(end)],
                    y=["rG4 candidates", "rG4 candidates"],
                    mode="lines",
                    line=dict(width=10),
                    name=str(row.candidate_type),
                    text=[
                        f"{row.sequence_id}<br>"
                        f"Candidate: {row.candidate_type}<br>"
                        f"Coordinates: {chrom}:{start}-{end}<br>"
                        f"Strand: {strand}"
                        for _ in range(2)
                    ],
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False,
                )
            )

    if not apa.empty:

        apa["midpoint"] = (
            apa["start"] + apa["end"]
        ) / 2

        for site_type in sorted(
            apa["site_type"].astype(str).unique()
        ):

            subset = apa[
                apa["site_type"].astype(str) == site_type
            ]

            fig.add_trace(
                go.Scatter(
                    x=subset["midpoint"],
                    y=["APA annotations"] * len(subset),
                    mode="markers",
                    name=f"APA: {site_type}",
                    marker=dict(size=13, symbol="diamond"),
                    text=[
                        f"{row.sequence_id}<br>"
                        f"Type: {row.site_type}<br>"
                        f"Coordinates: {row.chrom}:"
                        f"{row.start}-{row.end}<br>"
                        f"Strand: {row.strand}"
                        for row in subset.itertuples()
                    ],
                    hovertemplate="%{text}<extra></extra>",
                )
            )

    fig.update_layout(
        title=f"Exon-mapped rG4 candidates and APA sites: {chromosome}",
        xaxis_title="Genomic coordinate (bp)",
        yaxis_title="Feature",
        template="plotly_dark",
        height=450,
        hovermode="closest",
    )

    return fig


# ============================================================
# SESSION STATE
# ============================================================

defaults = {
    "fasta_data": None,
    "rg4_results": None,
    "apa_results": None,
    "integrated_results": None,
    "gene_panel": None,
}

for key, value in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = value


# ============================================================
# HEADER
# ============================================================

st.markdown("# 🧬 RNA-Insight")

st.markdown(
    "**In-Silico Analysis of RNA G-Quadruplex-Mediated "
    "Alternative Polyadenylation in Human Stress-Response Genes**"
)

st.caption(
    "Computational research workspace for human gene retrieval, "
    "RNA G4 prediction and APA annotation analysis."
)

st.divider()


# ============================================================
# SIDEBAR
# ============================================================

with st.sidebar:

    st.title("Research workspace")

    page = st.radio(
        "Select analysis",
        [
            "Human gene data retrieval",
            "Data intake & quality control",
            "rG4 sequence analysis",
            "APA site analysis",
            "Integrated rG4–APA analysis",
            "Results & report",
        ],
    )

    st.divider()

    st.caption("RNA-Insight | Computational research tool")


# ============================================================
# 00. HUMAN GENE DATA RETRIEVAL
# ============================================================

if page == "Human gene data retrieval":

    st.header("00 | Human stress-response gene retrieval")

    st.write(
        "Retrieve human gene annotations and transcript cDNA "
        "sequences directly from Ensembl."
    )

    st.info(
        "This page uses live Ensembl REST requests. "
        "Internet access is required."
    )

    st.subheader("Research gene panel")

    gene_input = st.text_area(
        "Enter human gene symbols separated by commas",
        value=", ".join(DEFAULT_GENES),
        help="Example: HSPA1A, HSP90AA1, DDIT3, ATF4",
    )

    max_transcripts = st.slider(
        "Maximum transcripts retrieved per gene",
        min_value=1,
        max_value=20,
        value=5,
    )

    if st.button("Retrieve human gene data"):

        symbols = [
            symbol.strip().upper()
            for symbol in gene_input.split(",")
            if symbol.strip()
        ]

        if not symbols:
            st.error("Enter at least one gene symbol.")

        elif len(set(symbols)) != len(symbols):
            st.error("Remove duplicate gene symbols.")

        else:

            with st.spinner(
                "Connecting to Ensembl and retrieving transcripts..."
            ):

                try:

                    result = fetch_gene_panel(
                        gene_symbols=symbols,
                        max_transcripts_per_gene=max_transcripts,
                    )

                    st.session_state.gene_panel = result

                except Exception as exc:

                    st.error(
                        f"Gene retrieval failed: {exc}"
                    )

    result = st.session_state.gene_panel

    if result is not None:

        genes = result["genes"]
        transcripts = result["transcripts"]
        errors = result["errors"]

        st.divider()
        st.subheader("Retrieval summary")

        c1, c2, c3 = st.columns(3)

        c1.metric("Genes retrieved", len(genes))
        c2.metric("Transcript records", len(transcripts))

        successful_sequences = (
            int(
                transcripts["sequence_status"].eq("retrieved").sum()
            )
            if not transcripts.empty
            else 0
        )

        c3.metric(
            "Sequences retrieved",
            successful_sequences,
        )

        st.subheader("Human gene annotations")

        if not genes.empty:

            st.dataframe(
                genes,
                use_container_width=True,
                hide_index=True,
            )

            st.download_button(
                "Download gene annotations CSV",
                genes.to_csv(index=False).encode("utf-8"),
                file_name="RNA_Insight_human_gene_annotations.csv",
                mime="text/csv",
            )

        else:

            st.warning("No gene annotations were retrieved.")

        st.subheader("Transcript annotations")

        if not transcripts.empty:

            st.dataframe(
                transcripts,
                use_container_width=True,
                hide_index=True,
            )

            st.download_button(
                "Download transcript metadata CSV",
                transcripts.to_csv(index=False).encode("utf-8"),
                file_name="RNA_Insight_transcript_metadata.csv",
                mime="text/csv",
            )

        st.subheader("Retrieved RNA sequences")

        if result.get("fasta"):

            st.download_button(
                "Download transcript FASTA",
                result["fasta"].encode("utf-8"),
                file_name="RNA_Insight_human_transcripts.fasta",
                mime="text/plain",
            )

            st.success(
                "Transcript FASTA is ready for downstream rG4 screening."
            )

        if not errors.empty:

            st.warning(
                "Some gene or transcript retrieval operations failed."
            )

            st.dataframe(
                errors,
                use_container_width=True,
                hide_index=True,
            )

        st.caption(
            "Source: Ensembl REST API | Species: Homo sapiens."
        )

        st.warning(
            "Transcript sequences are mature cDNA sequences. "
            "Exon-aware genomic mapping is required before genomic "
            "APA integration."
        )


# ============================================================
# 01. DATA INTAKE AND QUALITY CONTROL
# ============================================================

elif page == "Data intake & quality control":

    st.header("01 | Data intake & quality control")

    st.write(
        "Upload RNA or DNA sequences in FASTA format. "
        "DNA sequences are normalized to RNA."
    )

    fasta_file = st.file_uploader(
        "Upload RNA or DNA FASTA",
        type=["fasta", "fa", "fna", "txt"],
        key="fasta_upload",
    )

    if fasta_file is not None:

        try:

            raw_records = read_fasta(fasta_file)

            checked_records = validate_rna_records(raw_records)

            st.session_state.fasta_data = checked_records
            st.session_state.rg4_results = None
            st.session_state.integrated_results = None

            st.success(
                f"Loaded {len(checked_records)} FASTA records."
            )

            st.dataframe(
                checked_records,
                use_container_width=True,
                hide_index=True,
            )

            valid_count = int(
                checked_records["valid_sequence"].sum()
            )

            c1, c2, c3 = st.columns(3)

            c1.metric(
                "Total sequences",
                len(checked_records),
            )

            c2.metric(
                "Valid sequences",
                valid_count,
            )

            c3.metric(
                "Invalid sequences",
                len(checked_records) - valid_count,
            )

        except Exception as exc:

            st.error(
                f"Could not process FASTA: {exc}"
            )

    st.divider()

    st.subheader("Example dataset")

    st.write(
        "Synthetic demonstration FASTA: "
        "`sample_data/example_stress_response_rna.fasta`"
    )


# ============================================================
# 02. rG4 SEQUENCE ANALYSIS
# ============================================================

elif page == "rG4 sequence analysis":

    st.header("02 | RNA G-quadruplex sequence analysis")

    records = st.session_state.fasta_data

    if records is None:

        st.info(
            "Upload a FASTA file in Data intake & quality control first."
        )

    else:

        c1, c2 = st.columns(2)

        with c1:

            window_size = st.slider(
                "G4Hunter sliding window size (nt)",
                10,
                100,
                25,
            )

        with c2:

            threshold = st.slider(
                "Minimum G4Hunter score",
                0.5,
                4.0,
                1.2,
                0.1,
            )

        st.caption(
            "G4Hunter-style sliding-window screening and "
            "canonical motif screening are performed separately."
        )

        if st.button("Run rG4 screening"):

            try:

                valid_records = records[
                    records["valid_sequence"]
                ].copy()

                results = scan_rg4(
                    valid_records,
                    window_size=window_size,
                    score_threshold=threshold,
                )

                st.session_state.rg4_results = results
                st.session_state.integrated_results = None

                st.success(
                    f"Screening completed. {len(results)} "
                    "candidate intervals identified."
                )

            except Exception as exc:

                st.error(
                    f"rG4 analysis failed: {exc}"
                )

        results = st.session_state.rg4_results

        if results is not None:

            if results.empty:

                st.warning(
                    "No candidate regions detected."
                )

            else:

                g4hunter = results[
                    results["candidate_type"] == "g4hunter"
                ]

                canonical = results[
                    results["candidate_type"] == "canonical_motif"
                ]

                c1, c2, c3 = st.columns(3)

                c1.metric(
                    "Total candidate intervals",
                    len(results),
                )

                c2.metric(
                    "G4Hunter-labelled regions",
                    len(g4hunter),
                )

                c3.metric(
                    "Canonical motifs",
                    len(canonical),
                )

                st.dataframe(
                    results,
                    use_container_width=True,
                    hide_index=True,
                )

                st.subheader("G4Hunter score distribution")

                if not g4hunter.empty:

                    fig = px.histogram(
                        g4hunter,
                        x="score",
                        nbins=15,
                        labels={
                            "score": "Candidate score",
                            "count": "Candidate count",
                        },
                    )

                    st.plotly_chart(
                        fig,
                        use_container_width=True,
                    )

                else:

                    st.info(
                        "No G4Hunter sliding-window candidates detected."
                    )

                st.subheader("Canonical motif results")

                if not canonical.empty:

                    display_columns = [
                        column
                        for column in [
                            "sequence_id",
                            "start",
                            "end",
                            "candidate_sequence",
                            "score",
                            "score_type",
                            "length",
                        ]
                        if column in canonical.columns
                    ]

                    st.dataframe(
                        canonical[display_columns],
                        use_container_width=True,
                        hide_index=True,
                    )

                else:

                    st.info(
                        "No canonical motifs detected."
                    )

                st.download_button(
                    "Download rG4 results CSV",
                    results.to_csv(index=False).encode("utf-8"),
                    file_name="RNA_Insight_rg4_candidates.csv",
                    mime="text/csv",
                )

        st.info(
            "Predicted candidates are not experimentally confirmed "
            "RNA G-quadruplex structures."
        )


# ============================================================
# 03. APA SITE ANALYSIS
# ============================================================

elif page == "APA site analysis":

    st.header("03 | Alternative polyadenylation site analysis")

    st.write(
        "Upload an APA annotation table with genomic coordinates."
    )

    apa_file = st.file_uploader(
        "Upload APA annotation table",
        type=["csv", "tsv", "txt", "bed", "gtf"],
        key="apa_upload",
    )

    if apa_file is not None:

        try:

            file_type = detect_input(
                apa_file,
                apa_file.name,
            )

            if file_type == "bed":

                raw_apa = read_bed(apa_file)

                if "name" in raw_apa.columns:

                    raw_apa = raw_apa.rename(
                        columns={"name": "sequence_id"}
                    )

                if "score" in raw_apa.columns:

                    raw_apa = raw_apa.drop(
                        columns=["score"]
                    )

            elif file_type == "gtf":

                raw_apa = read_gtf(apa_file)

                raw_apa["sequence_id"] = (
                    raw_apa["attributes"].apply(
                        lambda value: value.get(
                            "transcript_id",
                            value.get("gene_id", "unknown"),
                        )
                    )
                )

            else:

                raw_apa = read_tabular(apa_file)

            apa = standardize_apa_sites(raw_apa)

            st.session_state.apa_results = apa
            st.session_state.integrated_results = None

            summary = summarize_apa_sites(apa)

            c1, c2, c3 = st.columns(3)

            c1.metric(
                "Total APA sites",
                summary["total_sites"],
            )

            c2.metric(
                "Proximal sites",
                summary["proximal_sites"],
            )

            c3.metric(
                "Distal sites",
                summary["distal_sites"],
            )

            st.subheader("Standardized APA annotations")

            st.dataframe(
                apa,
                use_container_width=True,
                hide_index=True,
            )

            st.download_button(
                "Download standardized APA table",
                apa.to_csv(index=False).encode("utf-8"),
                file_name="RNA_Insight_apa_sites.csv",
                mime="text/csv",
            )

        except Exception as exc:

            st.error(
                f"APA processing failed: {exc}"
            )

    else:

        st.info(
            "Required columns: sequence_id, chrom, start, end, "
            "strand. site_type is optional."
        )

    st.caption(
        "APA annotation counts do not measure differential APA usage."
    )


# ============================================================
# 04. EXON-AWARE INTEGRATED ANALYSIS
# ============================================================

elif page == "Integrated rG4–APA analysis":

    st.header("04 | Integrated rG4–APA analysis")

    st.write(
        "Map predicted RNA G4 candidate intervals through transcript "
        "exons and compare their genomic segments with APA annotations."
    )

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results

    if rg4 is None:

        st.warning(
            "Run rG4 screening before starting integration."
        )

    elif apa is None:

        st.warning(
            "Upload and standardize an APA annotation table first."
        )

    else:

        st.subheader("Exon-level transcript mapping")

        st.info(
            "Upload an exon-level mapping CSV. Each transcript "
            "must have one row per exon. Do not upload a single "
            "whole-transcript genomic interval."
        )

        st.markdown("""
**Required mapping columns**

- `sequence_id`
- `chrom`
- `start`
- `end`
- `strand`
- `transcript_start`
- `transcript_end`

Coordinates must be 0-based and half-open.
""")

        mapping_file = st.file_uploader(
            "Upload exon-level transcript mapping CSV",
            type=["csv", "tsv", "txt"],
            key="exon_mapping_upload",
        )

        if mapping_file is not None:

            try:

                mapping = read_tabular(mapping_file)

                mapping.columns = [
                    str(column).strip().lower()
                    for column in mapping.columns
                ]

                required_mapping_columns = {
                    "sequence_id",
                    "chrom",
                    "start",
                    "end",
                    "strand",
                    "transcript_start",
                    "transcript_end",
                }

                missing_columns = (
                    required_mapping_columns - set(mapping.columns)
                )

                if missing_columns:

                    st.error(
                        "Mapping file is missing required columns: "
                        + ", ".join(sorted(missing_columns))
                    )

                else:

                    st.success(
                        f"Loaded {len(mapping)} exon mapping rows."
                    )

                    st.dataframe(
                        mapping.head(20),
                        use_container_width=True,
                        hide_index=True,
                    )

                    if st.button(
                        "Run exon-aware rG4–APA integration"
                    ):

                        try:

                            with st.spinner(
                                "Mapping candidates through exons "
                                "and checking APA overlap..."
                            ):

                                integrated = integrate_rg4_apa(
                                    rg4,
                                    apa,
                                    mapping,
                                )

                            st.session_state.integrated_results = (
                                integrated
                            )

                            st.success(
                                "Exon-aware integration completed."
                            )

                        except Exception as exc:

                            st.session_state.integrated_results = None

                            st.error(
                                f"Integration failed: {exc}"
                            )

            except Exception as exc:

                st.error(
                    f"Could not read mapping file: {exc}"
                )

        integrated = st.session_state.integrated_results

        if integrated is not None:

            st.divider()

            st.subheader("Integrated analysis results")

            if integrated.empty:

                st.warning(
                    "No mapped candidate intervals were returned."
                )

            else:

                overlap_count = int(
                    integrated["overlaps_apa"].sum()
                )

                c1, c2, c3 = st.columns(3)

                c1.metric(
                    "Mapped candidate intervals",
                    len(integrated),
                )

                c2.metric(
                    "Candidates overlapping APA sites",
                    overlap_count,
                )

                c3.metric(
                    "Candidates without APA overlap",
                    len(integrated) - overlap_count,
                )

                st.dataframe(
                    integrated,
                    use_container_width=True,
                    hide_index=True,
                )

                st.download_button(
                    "Download integrated results CSV",
                    integrated.to_csv(index=False).encode("utf-8"),
                    file_name=(
                        "RNA_Insight_exon_aware_integrated_results.csv"
                    ),
                    mime="text/csv",
                )

                st.divider()

                st.subheader("Genomic visualization")

                chromosomes = sorted(
                    integrated["chrom"].astype(str).unique()
                )

                chromosome = st.selectbox(
                    "Select chromosome",
                    chromosomes,
                )

                fig = create_genomic_figure(
                    integrated,
                    apa,
                    chromosome,
                )

                st.plotly_chart(
                    fig,
                    use_container_width=True,
                )

                st.caption(
                    "Genomic overlap indicates positional association "
                    "between predicted rG4 candidate segments and "
                    "annotated APA sites. It does not establish "
                    "causal regulation or differential APA usage."
                )

                st.warning(
                    "Verify that exon mapping and APA annotations use "
                    "the same genome assembly, chromosome naming "
                    "convention and coordinate system."
                )


# ============================================================
# 05. RESULTS AND REPORT
# ============================================================

elif page == "Results & report":

    st.header("05 | Results & research report")

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results
    integrated = st.session_state.integrated_results

    if rg4 is None:
        rg4 = pd.DataFrame()

    if apa is None:
        apa = pd.DataFrame()

    c1, c2, c3 = st.columns(3)

    c1.metric(
        "rG4 candidate intervals",
        len(rg4),
    )

    c2.metric(
        "APA annotation sites",
        len(apa),
    )

    c3.metric(
        "Integrated candidate intervals",
        len(integrated) if integrated is not None else 0,
    )

    report = build_report(
        rg4_candidates=rg4,
        apa_sites=apa,
        integrated_results=integrated,
    )

    st.download_button(
        "Download research report",
        report.encode("utf-8"),
        file_name="RNA_Insight_Research_Report.md",
        mime="text/markdown",
    )

    st.subheader("Report preview")

    st.markdown(report)
