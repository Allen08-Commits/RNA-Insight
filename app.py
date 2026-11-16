
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
from core.mapping import build_transcript_exon_mapping


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
# CUSTOM STYLING
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
                        (
                            f"Transcript: {row.sequence_id}<br>"
                            f"Candidate: {row.candidate_type}<br>"
                            f"Genomic region: {chrom}:{start}-{end}<br>"
                            f"Strand: {strand}"
                        )
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
                        (
                            f"Sequence: {row.sequence_id}<br>"
                            f"APA type: {row.site_type}<br>"
                            f"Region: {row.chrom}:"
                            f"{row.start}-{row.end}<br>"
                            f"Strand: {row.strand}"
                        )
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
    "exon_mapping": None,
    "mapping_status": None,
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
    "RNA G4 prediction, APA annotation and exon-aware genomic integration."
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
        "sequences from Ensembl."
    )

    st.info(
        "This module uses the Ensembl REST API. "
        "An internet connection is required."
    )

    st.subheader("Research gene panel")

    gene_input = st.text_area(
        "Enter human gene symbols separated by commas",
        value=", ".join(DEFAULT_GENES),
        help="Example: HSPA1A, HSP90AA1, DDIT3, ATF4",
    )

    max_transcripts = st.slider(
        "Maximum transcripts per gene",
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

                    # Invalidate downstream results when the
                    # selected gene panel changes.
                    st.session_state.exon_mapping = None
                    st.session_state.mapping_status = None
                    st.session_state.integrated_results = None

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

        successful_sequences = (
            int(
                transcripts["sequence_status"]
                .eq("retrieved").sum()
            )
            if not transcripts.empty
            else 0
        )

        c1, c2, c3 = st.columns(3)

        c1.metric("Genes retrieved", len(genes))
        c2.metric("Transcript records", len(transcripts))
        c3.metric("Sequences retrieved", successful_sequences)

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
                "Transcript FASTA is ready for rG4 screening."
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
            "Transcript cDNA coordinates are not continuous genomic "
            "coordinates. Exon-aware mapping is required for genomic "
            "integration."
        )


# ============================================================
# 01. DATA INTAKE AND QUALITY CONTROL
# ============================================================

elif page == "Data intake & quality control":

    st.header("01 | Data intake & quality control")

    st.write(
        "Upload RNA or DNA sequences in FASTA format."
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
            "Sliding-window screening and canonical motif screening "
            "are reported separately."
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
                    "Canonical motif candidates",
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
                        "No G4Hunter-labelled candidates detected."
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
                        "No canonical motif candidates detected."
                    )

                st.download_button(
                    "Download rG4 results CSV",
                    results.to_csv(index=False).encode("utf-8"),
                    file_name="RNA_Insight_rg4_candidates.csv",
                    mime="text/csv",
                )

        st.info(
            "Predicted candidates are computational predictions, "
            "not experimentally confirmed RNA G-quadruplex structures."
        )


# ============================================================
# 03. APA SITE ANALYSIS
# ============================================================

elif page == "APA site analysis":

    st.header("03 | Alternative polyadenylation site analysis")

    st.write(
        "Upload an APA annotation table containing genomic coordinates."
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
            "Required columns: sequence_id, chrom, start, end, strand. "
            "site_type is optional."
        )

    st.caption(
        "APA annotation counts do not measure differential APA usage."
    )


# ============================================================
# 04. AUTOMATED EXON-AWARE INTEGRATED ANALYSIS
# ============================================================

elif page == "Integrated rG4–APA analysis":

    st.header("04 | Integrated rG4–APA analysis")

    st.write(
        "Retrieve exon coordinates automatically, validate transcript "
        "mappings and integrate predicted rG4 candidates with APA sites."
    )

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results
    gene_panel = st.session_state.gene_panel

    if rg4 is None:

        st.warning(
            "Run rG4 screening before starting integration."
        )

    elif apa is None:

        st.warning(
            "Upload and standardize APA annotations first."
        )

    elif gene_panel is None:

        st.warning(
            "Retrieve human gene and transcript data first."
        )

    else:

        transcripts = gene_panel["transcripts"].copy()

        st.subheader("Step 1 | Transcript selection")

        if transcripts.empty:

            st.error(
                "No transcript records are available."
            )

        else:

            if "sequence_status" in transcripts.columns:

                successful = transcripts[
                    transcripts["sequence_status"].astype(str)
                    == "retrieved"
                ].copy()

            else:

                successful = transcripts.copy()

            st.metric(
                "Transcripts with retrieved cDNA",
                len(successful),
            )

            if successful.empty:

                st.warning(
                    "No successfully retrieved transcript sequences "
                    "are available for exon mapping."
                )

            elif "sequence_length" not in successful.columns:

                st.error(
                    "Transcript metadata does not contain sequence_length. "
                    "Check the gene_data.py output."
                )

            else:

                st.dataframe(
                    successful[
                        [
                            column
                            for column in [
                                "gene_symbol",
                                "transcript_id",
                                "sequence_length",
                                "assembly_name",
                            ]
                            if column in successful.columns
                        ]
                    ],
                    use_container_width=True,
                    hide_index=True,
                )

                st.divider()

                st.subheader(
                    "Step 2 | Automatic Ensembl exon retrieval"
                )

                st.write(
                    "Retrieve exon coordinates for the successfully "
                    "downloaded transcript sequences."
                )

                if st.button(
                    "Generate exon mapping automatically"
                ):

                    try:

                        with st.spinner(
                            "Retrieving exon annotations from Ensembl "
                            "and validating transcript coordinates..."
                        ):

                            mapping, status = (
                                build_transcript_exon_mapping(
                                    successful
                                )
                            )

                        st.session_state.exon_mapping = mapping
                        st.session_state.mapping_status = status
                        st.session_state.integrated_results = None

                        st.success(
                            "Exon mapping retrieval completed."
                        )

                    except Exception as exc:

                        st.error(
                            f"Exon mapping failed: {exc}"
                        )

                mapping = st.session_state.exon_mapping
                status = st.session_state.mapping_status

                if status is not None:

                    st.divider()

                    st.subheader(
                        "Step 3 | Transcript mapping validation"
                    )

                    validated = status[
                        status["mapping_status"] == "validated"
                    ]

                    failed = status[
                        status["mapping_status"] == "failed"
                    ]

                    c1, c2, c3 = st.columns(3)

                    c1.metric(
                        "Transcripts processed",
                        len(status),
                    )

                    c2.metric(
                        "Successfully validated",
                        len(validated),
                    )

                    c3.metric(
                        "Failed mappings",
                        len(failed),
                    )

                    st.dataframe(
                        status,
                        use_container_width=True,
                        hide_index=True,
                    )

                    if not failed.empty:

                        with st.expander(
                            "View failed transcript mappings"
                        ):

                            st.dataframe(
                                failed,
                                use_container_width=True,
                                hide_index=True,
                            )

                    if mapping is not None and not mapping.empty:

                        st.subheader(
                            "Validated exon mapping table"
                        )

                        st.dataframe(
                            mapping,
                            use_container_width=True,
                            hide_index=True,
                        )

                        st.download_button(
                            "Download validated exon mapping CSV",
                            mapping.to_csv(
                                index=False
                            ).encode("utf-8"),
                            file_name=(
                                "RNA_Insight_validated_exon_mapping.csv"
                            ),
                            mime="text/csv",
                        )

                        st.success(
                            f"{len(mapping)} validated exon records "
                            "are available."
                        )

                        st.divider()

                        st.subheader(
                            "Step 4 | rG4–APA integration"
                        )

                        if st.button(
                            "Run integrated rG4–APA analysis"
                        ):

                            try:

                                with st.spinner(
                                    "Mapping candidate intervals through "
                                    "transcript exons and comparing "
                                    "their genomic segments with APA sites..."
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
                                    "Integrated analysis completed."
                                )

                            except Exception as exc:

                                st.session_state.integrated_results = None

                                st.error(
                                    f"Integration failed: {exc}"
                                )

                    else:

                        st.warning(
                            "No validated exon mappings are available."
                        )

                integrated = st.session_state.integrated_results

                if integrated is not None:

                    st.divider()

                    st.subheader(
                        "Step 5 | Integrated research results"
                    )

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
                            "Mapped candidates",
                            len(integrated),
                        )

                        c2.metric(
                            "Candidates overlapping APA",
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
                            integrated.to_csv(
                                index=False
                            ).encode("utf-8"),
                            file_name=(
                                "RNA_Insight_integrated_results.csv"
                            ),
                            mime="text/csv",
                        )

                        if "chrom" in integrated.columns:

                            chromosomes = sorted(
                                integrated["chrom"]
                                .astype(str)
                                .unique()
                            )

                            if chromosomes:

                                st.subheader(
                                    "Genomic coordinate visualization"
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

                        st.warning(
                            "Genomic overlap indicates positional "
                            "association only. It does not establish "
                            "causal regulation or differential APA usage."
                        )

                        st.caption(
                            "Verify that Ensembl exon mappings and APA "
                            "annotations use the same genome assembly, "
                            "chromosome naming convention and coordinate "
                            "system."
                        )


# ============================================================
# 05. RESULTS AND RESEARCH REPORT
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

    if integrated is None:
        integrated_for_report = None
    else:
        integrated_for_report = integrated

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
        "Integrated candidates",
        (
            len(integrated)
            if integrated is not None
            else 0
        ),
    )

    report = build_report(
        rg4_candidates=rg4,
        apa_sites=apa,
        integrated_results=integrated_for_report,
    )

    st.download_button(
        "Download research report",
        report.encode("utf-8"),
        file_name="RNA_Insight_Research_Report.md",
        mime="text/markdown",
    )

    st.divider()

    st.subheader("Report preview")

    st.markdown(report)

    st.caption(
        "Interpret computational predictions as candidate associations. "
        "Experimental validation is required to establish biological "
        "structure and regulatory function."
    )
