
"""RNA-Insight: rG4 and alternative polyadenylation research workspace."""

from __future__ import annotations

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


st.set_page_config(
    page_title="RNA-Insight",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# --------------------------------------------------
# VISUAL STYLE
# --------------------------------------------------

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
</style>
""", unsafe_allow_html=True)


# --------------------------------------------------
# HELPER: GENOMIC VISUALIZATION
# --------------------------------------------------

def create_genomic_figure(
    integrated: pd.DataFrame,
    apa_sites: pd.DataFrame,
    chromosome: str,
) -> go.Figure:
    """Plot mapped candidate intervals and APA annotations."""

    fig = go.Figure()

    candidates = integrated[
        integrated["chrom"].astype(str) == str(chromosome)
    ].copy()

    apa = apa_sites[
        apa_sites["chrom"].astype(str) == str(chromosome)
    ].copy()

    candidate_types = sorted(
        candidates["candidate_type"].astype(str).unique()
    )

    # Plot candidate intervals as horizontal genomic segments.
    for candidate_type in candidate_types:

        subset = candidates[
            candidates["candidate_type"].astype(str)
            == candidate_type
        ]

        x_values = []
        y_values = []
        hover_values = []

        for row in subset.itertuples():

            label = (
                f"{row.sequence_id}<br>"
                f"{row.candidate_type}<br>"
                f"Coordinates: {row.genomic_start}-"
                f"{row.genomic_end}<br>"
                f"Score: {row.score}<br>"
                f"APA overlap: {row.overlaps_apa}"
            )

            x_values.extend([
                row.genomic_start,
                row.genomic_end,
                None,
            ])

            y_values.extend([
                "rG4 candidates",
                "rG4 candidates",
                None,
            ])

            hover_values.extend([
                label,
                label,
                None,
            ])

        fig.add_trace(
            go.Scatter(
                x=x_values,
                y=y_values,
                mode="lines",
                name=f"rG4: {candidate_type}",
                line=dict(width=9),
                text=hover_values,
                hovertemplate="%{text}<extra></extra>",
                connectgaps=False,
            )
        )

    # Plot APA annotation intervals as points at interval midpoints.
    if not apa.empty:

        apa["midpoint"] = (
            apa["start"] + apa["end"]
        ) / 2

        apa["site_label"] = (
            apa["site_type"].astype(str)
            + " APA site"
        )

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
                    marker=dict(
                        size=14,
                        symbol="diamond",
                        line=dict(width=1),
                    ),
                    text=[
                        (
                            f"{row.sequence_id}<br>"
                            f"{row.site_type} APA site<br>"
                            f"{row.chrom}:{row.start}-{row.end}<br>"
                            f"Strand: {row.strand}"
                        )
                        for row in subset.itertuples()
                    ],
                    hovertemplate="%{text}<extra></extra>",
                )
            )

    fig.update_layout(
        title=f"rG4 candidates and APA annotations on {chromosome}",
        xaxis_title="Genomic coordinate (bp)",
        yaxis_title="Feature",
        template="plotly_dark",
        height=430,
        hovermode="closest",
        legend_title="Feature type",
        margin=dict(l=20, r=20, t=60, b=30),
    )

    fig.update_yaxes(
        categoryorder="array",
        categoryarray=[
            "APA annotations",
            "rG4 candidates",
        ],
    )

    return fig


# --------------------------------------------------
# HEADER
# --------------------------------------------------

st.markdown("# 🧬 RNA-Insight")

st.markdown(
    "**In-Silico Analysis of RNA G-Quadruplex-Mediated "
    "Alternative Polyadenylation in Human Stress-Response Genes**"
)

st.caption(
    "A computational research workspace for RNA G4 prediction "
    "and alternative polyadenylation annotation analysis."
)

st.divider()


# --------------------------------------------------
# SESSION STATE
# --------------------------------------------------

defaults = {
    "fasta_data": None,
    "rg4_results": None,
    "apa_results": None,
    "integrated_results": None,
}

for key, value in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = value


# --------------------------------------------------
# NAVIGATION
# --------------------------------------------------

with st.sidebar:

    st.title("Research workspace")

    page = st.radio(
        "Select analysis",
        [
            "Data intake & quality control",
            "rG4 sequence analysis",
            "APA site analysis",
            "Integrated rG4–APA analysis",
            "Results & report",
        ],
    )

    st.divider()
    st.caption("RNA-Insight | Computational research tool")


# --------------------------------------------------
# 01 DATA INTAKE
# --------------------------------------------------

if page == "Data intake & quality control":

    st.header("01 | Data intake & quality control")

    st.write(
        "Upload RNA sequences in FASTA format. "
        "DNA FASTA sequences are automatically converted from T to U."
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

            # Invalidate analyses when the sequence input changes.
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

            invalid_count = len(checked_records) - valid_count

            c1, c2, c3 = st.columns(3)

            c1.metric("Total sequences", len(checked_records))
            c2.metric("Valid sequences", valid_count)
            c3.metric("Invalid sequences", invalid_count)

            if invalid_count:
                st.warning(
                    "Some records contain invalid or empty "
                    "sequence characters."
                )

        except Exception as exc:
            st.error(f"Could not process FASTA: {exc}")

    st.divider()

    st.subheader("Example dataset")

    st.write(
        "Synthetic demonstration FASTA: "
        "`sample_data/example_stress_response_rna.fasta`"
    )


# --------------------------------------------------
# 02 rG4 ANALYSIS
# --------------------------------------------------

elif page == "rG4 sequence analysis":

    st.header("02 | RNA G-quadruplex sequence analysis")

    records = st.session_state.fasta_data

    if records is None:

        st.info(
            "Upload a FASTA file in Data intake & quality control first."
        )

    else:

        st.subheader("Prediction parameters")

        c1, c2 = st.columns(2)

        with c1:
            window_size = st.slider(
                "G4Hunter sliding window size (nt)",
                min_value=10,
                max_value=100,
                value=25,
                step=1,
            )

        with c2:
            threshold = st.slider(
                "Minimum G4Hunter score",
                min_value=0.5,
                max_value=4.0,
                value=1.2,
                step=0.1,
            )

        st.caption(
            "G4Hunter-style screening evaluates G-rich sliding "
            "windows. Canonical four-tract motif screening is "
            "performed independently."
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
                    "candidate regions identified."
                )

            except Exception as exc:
                st.error(f"rG4 analysis failed: {exc}")

        results = st.session_state.rg4_results

        if results is not None:

            st.subheader("Candidate classification")

            total = len(results)

            if results.empty:

                st.warning(
                    "No candidates were detected with the "
                    "current parameters."
                )

            else:

                g4hunter_results = results[
                    results["candidate_type"] == "g4hunter"
                ]

                canonical_results = results[
                    results["candidate_type"] == "canonical_motif"
                ]

                c1, c2, c3 = st.columns(3)

                c1.metric("Total candidate intervals", total)
                c2.metric(
                    "G4Hunter regions",
                    len(g4hunter_results),
                )
                c3.metric(
                    "Canonical motifs",
                    len(canonical_results),
                )

                st.dataframe(
                    results,
                    use_container_width=True,
                    hide_index=True,
                )

                # G4Hunter score distribution only.
                st.subheader("G4Hunter score distribution")

                if not g4hunter_results.empty:

                    fig = px.histogram(
                        g4hunter_results,
                        x="score",
                        nbins=15,
                        title="G4Hunter sliding-window scores",
                        labels={
                            "score": "Mean G4Hunter window score",
                            "count": "Number of candidate regions",
                        },
                    )

                    st.plotly_chart(
                        fig,
                        use_container_width=True,
                    )

                else:

                    st.info(
                        "No G4Hunter sliding-window candidates "
                        "were detected. Canonical motif results "
                        "are shown separately below."
                    )

                # Canonical motif score display only.
                st.subheader("Canonical motif results")

                if not canonical_results.empty:

                    st.dataframe(
                        canonical_results[
                            [
                                "sequence_id",
                                "start",
                                "end",
                                "candidate_sequence",
                                "score",
                                "score_type",
                                "length",
                            ]
                        ],
                        use_container_width=True,
                        hide_index=True,
                    )

                    st.caption(
                        "Canonical motif scores are descriptive "
                        "motif composition scores, not G4Hunter "
                        "sliding-window scores."
                    )

                else:

                    st.info(
                        "No canonical four-tract motifs were detected."
                    )

                st.download_button(
                    "Download rG4 results CSV",
                    results.to_csv(index=False).encode("utf-8"),
                    file_name="rna_insight_rg4_candidates.csv",
                    mime="text/csv",
                )

        st.info(
            "Computational predictions identify candidate regions, "
            "not experimentally confirmed RNA G-quadruplex structures. "
            "Overlapping candidates are not necessarily independent "
            "folded structures."
        )


# --------------------------------------------------
# 03 APA ANALYSIS
# --------------------------------------------------

elif page == "APA site analysis":

    st.header("03 | Alternative polyadenylation site analysis")

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

            a, b, c = st.columns(3)

            a.metric("Total APA sites", summary["total_sites"])
            b.metric("Proximal sites", summary["proximal_sites"])
            c.metric("Distal sites", summary["distal_sites"])

            st.subheader("Standardized APA annotations")

            st.dataframe(
                apa,
                use_container_width=True,
                hide_index=True,
            )

            st.download_button(
                "Download standardized APA table",
                apa.to_csv(index=False).encode("utf-8"),
                file_name="rna_insight_apa_sites.csv",
                mime="text/csv",
            )

        except Exception as exc:
            st.error(
                f"APA annotation processing failed: {exc}"
            )

    else:

        st.info(
            "Upload an APA table with sequence_id, chrom, "
            "start, end and strand columns."
        )

    st.caption(
        "APA annotation counts do not measure differential APA usage."
    )


# --------------------------------------------------
# 04 INTEGRATED ANALYSIS
# --------------------------------------------------

elif page == "Integrated rG4–APA analysis":

    st.header("04 | Integrated rG4–APA analysis")

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results

    if rg4 is None or apa is None:

        st.info(
            "Run rG4 screening and upload an APA annotation table first."
        )

    else:

        mapping_file = st.file_uploader(
            "Upload sequence-to-genome mapping CSV",
            type=["csv", "tsv"],
            key="mapping_upload",
        )

        if mapping_file is not None:

            try:
                mapping = read_tabular(mapping_file)

                st.caption(
                    "Mapping coordinates must be 0-based, "
                    "half-open genomic intervals, using the same "
                    "genome assembly as the APA annotations."
                )

                if st.button("Run integrated analysis"):

                    integrated = integrate_rg4_apa(
                        rg4,
                        apa,
                        mapping,
                    )

                    st.session_state.integrated_results = integrated

                    st.success(
                        f"Integration completed: {len(integrated)} "
                        "mapped candidate regions."
                    )

            except Exception as exc:
                st.error(
                    f"Integration failed: {exc}"
                )

        integrated = st.session_state.integrated_results

        if integrated is not None:

            st.subheader("Integrated results")

            if integrated.empty:

                st.warning(
                    "No candidate regions could be mapped."
                )

            else:

                overlap_count = int(
                    integrated["overlaps_apa"].sum()
                )

                c1, c2 = st.columns(2)

                c1.metric(
                    "Mapped candidate regions",
                    len(integrated),
                )

                c2.metric(
                    "Candidate intervals overlapping APA",
                    overlap_count,
                )

                st.dataframe(
                    integrated,
                    use_container_width=True,
                    hide_index=True,
                )

                st.download_button(
                    "Download integrated results",
                    integrated.to_csv(index=False).encode("utf-8"),
                    file_name="rna_insight_integrated_results.csv",
                    mime="text/csv",
                )

                st.divider()
                st.subheader("Genomic coordinate visualization")

                chromosomes = sorted(
                    integrated["chrom"].astype(str).unique()
                )

                selected_chromosome = st.selectbox(
                    "Select chromosome",
                    chromosomes,
                    key="genomic_chromosome",
                )

                figure = create_genomic_figure(
                    integrated,
                    apa,
                    selected_chromosome,
                )

                st.plotly_chart(
                    figure,
                    use_container_width=True,
                )

                st.caption(
                    "Horizontal segments represent mapped candidate "
                    "intervals. Diamond markers represent APA "
                    "annotation positions. Hover over features for "
                    "their coordinates and annotation details."
                )

        st.info(
            "Genomic overlap is a positional relationship. It does "
            "not establish that an rG4 structure causes APA regulation."
        )


# --------------------------------------------------
# 05 RESULTS AND REPORT
# --------------------------------------------------

elif page == "Results & report":

    st.header("05 | Results & research report")

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results
    integrated = st.session_state.integrated_results

    if rg4 is None:
        rg4 = pd.DataFrame()

    if apa is None:
        apa = pd.DataFrame()

    if integrated is not None:

        st.metric(
            "Integrated candidate regions",
            len(integrated),
        )

    c1, c2 = st.columns(2)

    c1.metric(
        "rG4 candidate intervals",
        len(rg4),
    )

    c2.metric(
        "APA annotation sites",
        len(apa),
    )

    report = build_report(
        rg4_candidates=rg4,
        apa_sites=apa,
        integrated_results=integrated,
    )

    st.download_button(
        "Download research report",
        data=report.encode("utf-8"),
        file_name="RNA-Insight_Research_Report.md",
        mime="text/markdown",
    )

    st.markdown("### Report preview")
    st.markdown(report)
