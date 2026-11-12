"""RNA-Insight: rG4 and alternative polyadenylation research workspace."""

from __future__ import annotations

import io

import pandas as pd
import plotly.express as px
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


st.markdown("""
<style>
.stApp {
    background: linear-gradient(135deg, #07111f 0%, #0b1730 55%, #11152e 100%);
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
div[data-testid="stFileUploader"] {
    background: rgba(14, 35, 60, 0.6);
    border-radius: 10px;
    padding: 10px;
}
</style>
""", unsafe_allow_html=True)


st.markdown("# 🧬 RNA-Insight")
st.markdown(
    "**In-Silico Analysis of RNA G-Quadruplex-Mediated "
    "Alternative Polyadenylation in Human Stress-Response Genes**"
)
st.caption(
    "A sequence-screening and annotation-integration research workspace"
)

st.divider()


if "fasta_data" not in st.session_state:
    st.session_state.fasta_data = None

if "rg4_results" not in st.session_state:
    st.session_state.rg4_results = None

if "apa_results" not in st.session_state:
    st.session_state.apa_results = None

if "integrated_results" not in st.session_state:
    st.session_state.integrated_results = None


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


if page == "Data intake & quality control":
    st.header("01 | Data intake & quality control")

    st.write(
        "Upload RNA sequences in FASTA format. DNA FASTA sequences "
        "are automatically converted from T to U."
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
            st.success(f"Loaded {len(checked_records)} FASTA records.")

            st.dataframe(
                checked_records,
                use_container_width=True,
                hide_index=True,
            )

            valid_count = int(checked_records["valid_sequence"].sum())
            invalid_count = len(checked_records) - valid_count

            c1, c2, c3 = st.columns(3)
            c1.metric("Total sequences", len(checked_records))
            c2.metric("Valid sequences", valid_count)
            c3.metric("Invalid sequences", invalid_count)

            if invalid_count:
                st.warning(
                    "Some records contain invalid or empty sequence "
                    "characters. Review them before analysis."
                )

        except Exception as exc:
            st.error(f"Could not process FASTA: {exc}")

    st.divider()
    st.subheader("Example dataset")

    st.write(
        "You can also use the synthetic example FASTA file in "
        "`sample_data/example_stress_response_rna.fasta`."
    )


elif page == "rG4 sequence analysis":
    st.header("02 | RNA G-quadruplex sequence analysis")

    records = st.session_state.fasta_data

    if records is None:
        st.info("Upload a FASTA file in Data intake & quality control first.")
    else:
        window_size = st.slider(
            "Composition window size",
            min_value=10,
            max_value=100,
            value=25,
        )

        threshold = st.slider(
            "Minimum G/C composition-window score",
            min_value=0.0,
            max_value=1.0,
            value=0.4,
            step=0.05,
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
                    f"Screening completed. {len(results)} candidate "
                    "regions identified."
                )

            except Exception as exc:
                st.error(f"rG4 analysis failed: {exc}")

        results = st.session_state.rg4_results

        if results is not None:
            st.subheader("Candidate results")
            st.metric("Total candidate regions", len(results))

            if results.empty:
                st.warning(
                    "No candidates were detected with the current settings."
                )
            else:
                st.dataframe(
                    results,
                    use_container_width=True,
                    hide_index=True,
                )

                fig = px.histogram(
                    results,
                    x="score",
                    color="candidate_type",
                    title="Distribution of candidate scores",
                    labels={"score": "Composition-window score"},
                )
                st.plotly_chart(fig, use_container_width=True)

                st.download_button(
                    "Download rG4 results CSV",
                    results.to_csv(index=False).encode("utf-8"),
                    file_name="rna_insight_rg4_candidates.csv",
                    mime="text/csv",
                )

        st.caption(
            "The composition-window score is a simple heuristic. "
            "It is not a validated replacement for G4Hunter or "
            "experimental structural evidence."
        )


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
                    raw_apa = raw_apa.drop(columns=["score"])

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
            st.error(f"APA annotation processing failed: {exc}")

    else:
        st.info(
            "Upload an APA table with sequence_id, chrom, start, end "
            "and strand columns."
        )

    st.caption(
        "APA site counts are annotation summaries. They do not measure "
        "differential APA usage."
    )


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
                st.error(f"Integration failed: {exc}")

        integrated = st.session_state.integrated_results

        if integrated is not None:
            st.dataframe(
                integrated,
                use_container_width=True,
                hide_index=True,
            )

            if not integrated.empty:
                overlap_count = int(
                    integrated["overlaps_apa"].sum()
                )
                st.metric(
                    "rG4 candidate regions overlapping APA annotations",
                    overlap_count,
                )

            st.download_button(
                "Download integrated results",
                integrated.to_csv(index=False).encode("utf-8"),
                file_name="rna_insight_integrated_results.csv",
                mime="text/csv",
            )

    st.caption(
        "Overlap is a genomic proximity/overlap observation, not evidence "
        "of causal regulation. Confirm genome assembly and coordinate "
        "compatibility before interpreting results."
    )


elif page == "Results & report":
    st.header("05 | Results & research report")

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results
    integrated = st.session_state.integrated_results

    if rg4 is None:
        rg4 = pd.DataFrame()

    if apa is None:
        apa = pd.DataFrame()

    st.metric("rG4 candidate regions", len(rg4))
    st.metric("APA annotation sites", len(apa))

    if integrated is not None:
        st.metric("Integrated candidate regions", len(integrated))

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
