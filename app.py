
"""
RNA-Insight
In-Silico Analysis of RNA G-Quadruplex-Mediated
Alternative Polyadenylation in Human Stress-Response Genes
"""

from __future__ import annotations

import re
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from core.io import read_fasta, validate_rna_records
from core.rg4 import scan_rg4
from core.apa import fetch_apa_sites, summarize_apa_sites
from core.gene_data import DEFAULT_GENES, fetch_gene_panel
from core.mapping import build_transcript_exon_mapping
from core.report import build_report, build_pdf_report


# ==================================================
# PAGE CONFIGURATION
# ==================================================

st.set_page_config(
    page_title="RNA-Insight",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ==================================================
# CUSTOM UI
# ==================================================

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


# ==================================================
# HELPER FUNCTIONS
# ==================================================

def normalize_transcript_id(value):
    """Extract the Ensembl transcript ID from a FASTA header."""

    value = str(value).strip()

    if not value:
        return ""

    return value.split("|")[0].split()[0]


def normalize_chromosome(value):
    """Normalize chromosome labels for comparison."""

    value = str(value).strip().lower()

    if value.startswith("chr"):
        value = value[3:]

    return value


def check_hg38_assembly(mapping):
    """Check that Ensembl mapping uses the expected hg38 assembly."""

    if "assembly_name" not in mapping.columns:
        raise ValueError(
            "The Ensembl mapping does not contain assembly information."
        )

    assemblies = (
        mapping["assembly_name"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    if not assemblies:
        raise ValueError("Genome assembly information is unavailable.")

    incompatible = [
        assembly for assembly in assemblies
        if "GRCH38" not in assembly.upper()
    ]

    if incompatible:
        raise ValueError(
            "The mapping assembly may not be compatible with "
            f"PolyA_DB hg38. Retrieved assembly: {assemblies}"
        )


def integrate_gene_level_apa(rg4, apa, mapping):
    """
    Map transcript-level rG4 candidates to genomic exon segments
    and compare them with gene-level PolyA_DB annotations.

    This is a gene-level genomic overlap analysis, not a
    transcript-specific APA assignment or causal analysis.
    """

    required_rg4 = {
        "sequence_id", "start", "end", "candidate_type"
    }

    required_apa = {
        "sequence_id", "chrom", "start", "end",
        "strand", "site_type"
    }

    required_mapping = {
        "sequence_id", "gene_symbol", "chrom", "start", "end",
        "strand", "transcript_start", "transcript_end",
        "assembly_name"
    }

    for frame, required, label in [
        (rg4, required_rg4, "rG4"),
        (apa, required_apa, "APA"),
        (mapping, required_mapping, "Mapping")
    ]:
        missing = required - set(frame.columns)

        if missing:
            raise ValueError(
                f"{label} data is missing columns: "
                + ", ".join(sorted(missing))
            )

    candidates = rg4.copy()
    sites = apa.copy()
    exons = mapping.copy()

    candidates["sequence_id"] = (
        candidates["sequence_id"].apply(normalize_transcript_id)
    )

    exons["sequence_id"] = (
        exons["sequence_id"].astype(str).str.strip()
    )

    exons["gene_symbol"] = (
        exons["gene_symbol"].astype(str).str.upper().str.strip()
    )

    sites["sequence_id"] = (
        sites["sequence_id"].astype(str).str.upper().str.strip()
    )

    for frame, columns, label in [
        (candidates, ["start", "end"], "rG4"),
        (sites, ["start", "end"], "APA"),
        (exons, ["start", "end", "transcript_start", "transcript_end"],
         "Mapping")
    ]:
        for column in columns:
            frame[column] = pd.to_numeric(
                frame[column], errors="coerce"
            )

        if frame[columns].isna().any().any():
            raise ValueError(
                f"{label} data contains invalid coordinates."
            )

    if (candidates["end"] <= candidates["start"]).any():
        raise ValueError("Invalid rG4 candidate intervals.")

    if (sites["end"] <= sites["start"]).any():
        raise ValueError("Invalid APA intervals.")

    if (exons["end"] <= exons["start"]).any():
        raise ValueError("Invalid genomic exon intervals.")

    if (
        exons["transcript_end"] <= exons["transcript_start"]
    ).any():
        raise ValueError("Invalid transcript exon intervals.")

    results = []

    mapping_groups = {
        key: group.sort_values("transcript_start")
        for key, group in exons.groupby("sequence_id")
    }

    for candidate in candidates.itertuples(index=False):

        transcript_id = normalize_transcript_id(
            candidate.sequence_id
        )

        if transcript_id not in mapping_groups:
            continue

        exon_rows = mapping_groups[transcript_id]

        gene_symbol = str(
            exon_rows["gene_symbol"].iloc[0]
        ).upper()

        strand = str(exon_rows["strand"].iloc[0])
        chromosome = str(exon_rows["chrom"].iloc[0])

        local_start = int(candidate.start)
        local_end = int(candidate.end)

        transcript_length = int(
            exon_rows["transcript_end"].max()
        )

        if local_start < 0 or local_end > transcript_length:
            continue

        segments = []
        covered_bases = 0

        # Map transcript coordinates through individual exons.
        for exon in exon_rows.itertuples(index=False):

            overlap_start = max(
                local_start, int(exon.transcript_start)
            )

            overlap_end = min(
                local_end, int(exon.transcript_end)
            )

            if overlap_start >= overlap_end:
                continue

            covered_bases += overlap_end - overlap_start

            if strand == "+":

                genomic_start = (
                    int(exon.start)
                    + overlap_start
                    - int(exon.transcript_start)
                )

                genomic_end = (
                    int(exon.start)
                    + overlap_end
                    - int(exon.transcript_start)
                )

            elif strand == "-":

                genomic_start = (
                    int(exon.end)
                    - (
                        overlap_end
                        - int(exon.transcript_start)
                    )
                )

                genomic_end = (
                    int(exon.end)
                    - (
                        overlap_start
                        - int(exon.transcript_start)
                    )
                )

            else:
                raise ValueError(
                    f"Invalid strand for {transcript_id}."
                )

            segments.append({
                "chrom": str(exon.chrom),
                "start": genomic_start,
                "end": genomic_end,
                "strand": strand
            })

        # Reject incomplete mappings.
        if covered_bases != local_end - local_start:
            continue

        # Select gene-level PAS annotations matching gene,
        # chromosome and strand.
        gene_sites = sites[
            sites["sequence_id"].eq(gene_symbol)
            & sites["strand"].astype(str).eq(strand)
            & sites["chrom"].apply(normalize_chromosome).eq(
                normalize_chromosome(chromosome)
            )
        ].copy()

        matched_records = []

        for segment in segments:

            overlapping = gene_sites[
                (gene_sites["start"] < segment["end"])
                & (gene_sites["end"] > segment["start"])
            ]

            if not overlapping.empty:
                matched_records.extend(
                    overlapping.to_dict("records")
                )

        if matched_records:
            matched = (
                pd.DataFrame(matched_records)
                .drop_duplicates()
                .reset_index(drop=True)
            )
        else:
            matched = pd.DataFrame(columns=sites.columns)

        segment_labels = [
            (
                f"{segment['chrom']}:"
                f"{segment['start']}-{segment['end']}"
                f"({segment['strand']})"
            )
            for segment in segments
        ]

        if not matched.empty:

            apa_labels = [
                (
                    f"{row.chrom}:{row.start}-{row.end} "
                    f"({row.site_type})"
                )
                for row in matched.itertuples()
            ]

            apa_types = sorted(
                matched["site_type"]
                .fillna("unspecified")
                .astype(str)
                .unique()
                .tolist()
            )

        else:
            apa_labels = []
            apa_types = []

        results.append({
            "sequence_id": transcript_id,
            "gene_symbol": gene_symbol,
            "candidate_start": local_start,
            "candidate_end": local_end,
            "chrom": chromosome,
            "strand": strand,
            "candidate_type": candidate.candidate_type,
            "score": getattr(candidate, "score", None),
            "genomic_segments": "; ".join(segment_labels),
            "overlaps_gene_level_apa": not matched.empty,
            "overlapping_apa_sites": len(matched),
            "overlapping_apa_types": ", ".join(apa_types),
            "overlapping_apa_coordinates": "; ".join(apa_labels),
            "interpretation": "Gene-level genomic overlap"
        })

    output_columns = [
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
        "overlapping_apa_types",
        "overlapping_apa_coordinates",
        "interpretation"
    ]

    return pd.DataFrame(results, columns=output_columns)


def create_genomic_figure(integrated, apa_sites, chromosome):
    """Visualize mapped rG4 genomic segments and APA annotations."""

    fig = go.Figure()

    selected = integrated[
        integrated["chrom"].astype(str) == str(chromosome)
    ].copy()

    y_positions = {
        "rG4 candidates": 1,
        "APA annotations": 0
    }

    for row in selected.itertuples():

        segments = str(row.genomic_segments).split(";")

        for segment in segments:

            match = re.search(
                r":(\d+)-(\d+)\(",
                segment.strip()
            )

            if not match:
                continue

            start = int(match.group(1))
            end = int(match.group(2))

            fig.add_trace(
                go.Scatter(
                    x=[start, end],
                    y=[y_positions["rG4 candidates"]] * 2,
                    mode="lines",
                    line=dict(width=10),
                    name=f"rG4: {row.candidate_type}",
                    text=[
                        f"{row.sequence_id}<br>"
                        f"{row.gene_symbol}<br>"
                        f"{row.candidate_type}<br>"
                        f"{chromosome}:{start}-{end}",
                        f"{row.sequence_id}<br>"
                        f"{row.gene_symbol}<br>"
                        f"{row.candidate_type}<br>"
                        f"{chromosome}:{start}-{end}"
                    ],
                    hovertemplate="%{text}<extra></extra>",
                    showlegend=False
                )
            )

    apa = apa_sites[
        apa_sites["chrom"].apply(normalize_chromosome).eq(
            normalize_chromosome(chromosome)
        )
    ].copy()

    if not apa.empty:

        apa["midpoint"] = (
            apa["start"] + apa["end"]
        ) / 2

        fig.add_trace(
            go.Scatter(
                x=apa["midpoint"],
                y=[y_positions["APA annotations"]] * len(apa),
                mode="markers",
                marker=dict(size=10, symbol="diamond"),
                name="APA annotations",
                text=[
                    (
                        f"{row.sequence_id}<br>"
                        f"{row.site_type}<br>"
                        f"{row.chrom}:{row.start}-{row.end}"
                    )
                    for row in apa.itertuples()
                ],
                hovertemplate="%{text}<extra></extra>"
            )
        )

    fig.update_layout(
        title=f"rG4 candidates and APA annotations: {chromosome}",
        xaxis_title="Genomic coordinate (bp)",
        yaxis=dict(
            tickmode="array",
            tickvals=[0, 1],
            ticktext=["APA annotations", "rG4 candidates"]
        ),
        template="plotly_dark",
        height=430,
        hovermode="closest"
    )

    return fig


# ==================================================
# SESSION STATE
# ==================================================

defaults = {
    "fasta_data": None,
    "rg4_results": None,
    "apa_results": None,
    "integrated_results": None,
    "gene_panel": None,
    "exon_mapping": None,
    "mapping_status": None
}

for key, value in defaults.items():

    if key not in st.session_state:
        st.session_state[key] = value


# ==================================================
# APPLICATION HEADER
# ==================================================

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


# ==================================================
# SIDEBAR NAVIGATION
# ==================================================

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
            "Results & report"
        ]
    )

    st.divider()

    st.caption("RNA-Insight | Computational research tool")


# ==================================================
# 00 | HUMAN GENE DATA RETRIEVAL
# ==================================================

if page == "Human gene data retrieval":

    st.header("00 | Human stress-response gene retrieval")

    st.write(
        "Retrieve human gene annotations and transcript cDNA "
        "sequences from Ensembl."
    )

    st.info(
        "This page uses live Ensembl REST requests. "
        "Internet access is required."
    )

    st.subheader("Research gene panel")

    gene_input = st.text_area(
        "Enter human gene symbols separated by commas",
        value=", ".join(DEFAULT_GENES),
        help="Example: HSPA1A, HSP90AA1, DDIT3, ATF4"
    )

    max_transcripts = st.slider(
        "Maximum transcripts retrieved per gene",
        min_value=1,
        max_value=20,
        value=5
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
                        max_transcripts_per_gene=max_transcripts
                    )

                    st.session_state.gene_panel = result

                    # New gene retrieval requires fresh downstream analysis.
                    st.session_state.rg4_results = None
                    st.session_state.apa_results = None
                    st.session_state.integrated_results = None
                    st.session_state.exon_mapping = None
                    st.session_state.mapping_status = None

                except Exception as exc:
                    st.error(f"Gene retrieval failed: {exc}")

    result = st.session_state.gene_panel

    if result is not None:

        genes = result["genes"]
        transcripts = result["transcripts"]
        errors = result["errors"]

        st.divider()
        st.subheader("Retrieval summary")

        successful_sequences = (
            int(transcripts["sequence_status"].eq("retrieved").sum())
            if not transcripts.empty else 0
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
                hide_index=True
            )

            st.download_button(
                "Download gene annotations CSV",
                genes.to_csv(index=False).encode("utf-8"),
                file_name="RNA_Insight_human_gene_annotations.csv",
                mime="text/csv"
            )

        else:
            st.warning("No gene annotations were retrieved.")

        st.subheader("Transcript annotations")

        if not transcripts.empty:

            st.dataframe(
                transcripts,
                use_container_width=True,
                hide_index=True
            )

            st.download_button(
                "Download transcript metadata CSV",
                transcripts.to_csv(index=False).encode("utf-8"),
                file_name="RNA_Insight_transcript_metadata.csv",
                mime="text/csv"
            )

        st.subheader("Retrieved RNA sequences")

        if result["fasta"]:

            st.download_button(
                "Download transcript FASTA",
                result["fasta"].encode("utf-8"),
                file_name="RNA_Insight_human_transcripts.fasta",
                mime="text/plain"
            )

            st.success(
                "Transcript FASTA is ready for downstream rG4 screening."
            )

        if not errors.empty:

            st.warning(
                "Some gene or transcript lookups could not be completed."
            )

            st.dataframe(
                errors,
                use_container_width=True,
                hide_index=True
            )

        st.caption(
            "Source: Ensembl REST API | Species: Homo sapiens. "
            "Record the retrieval date, Ensembl identifiers and "
            "genome assembly when reporting these data."
        )

        st.warning(
            "Transcript sequences are mature cDNA sequences. "
            "Their coordinates must be mapped through individual "
            "exons before genomic analysis."
        )


# ==================================================
# 01 | DATA INTAKE AND QUALITY CONTROL
# ==================================================

elif page == "Data intake & quality control":

    st.header("01 | Data intake & quality control")

    st.write(
        "Upload RNA sequences in FASTA format. "
        "DNA FASTA sequences are automatically converted from T to U."
    )

    fasta_file = st.file_uploader(
        "Upload RNA or DNA FASTA",
        type=["fasta", "fa", "fna", "txt"],
        key="fasta_upload"
    )

    if fasta_file is not None:

        try:

            raw_records = read_fasta(fasta_file)
            checked_records = validate_rna_records(raw_records)

            st.session_state.fasta_data = checked_records
            st.session_state.rg4_results = None
            st.session_state.integrated_results = None
            st.session_state.exon_mapping = None
            st.session_state.mapping_status = None

            st.success(
                f"Loaded {len(checked_records)} FASTA records."
            )

            st.dataframe(
                checked_records,
                use_container_width=True,
                hide_index=True
            )

            valid_count = int(
                checked_records["valid_sequence"].sum()
            )

            c1, c2, c3 = st.columns(3)

            c1.metric("Total sequences", len(checked_records))
            c2.metric("Valid sequences", valid_count)
            c3.metric(
                "Invalid sequences",
                len(checked_records) - valid_count
            )

        except Exception as exc:
            st.error(f"Could not process FASTA: {exc}")

    st.divider()

    st.subheader("Example dataset")

    st.write(
        "Synthetic demonstration FASTA: "
        "`sample_data/example_stress_response_rna.fasta`"
    )

    if st.session_state.fasta_data is not None:

        st.info(
            "For automatic genomic integration, use the transcript "
            "FASTA downloaded from Human gene data retrieval. "
            "This ensures that candidate IDs can be matched to Ensembl."
        )


# ==================================================
# 02 | rG4 SEQUENCE ANALYSIS
# ==================================================

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
                10, 100, 25
            )

        with c2:

            threshold = st.slider(
                "Minimum G4Hunter score",
                0.5, 4.0, 1.2, 0.1
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
                    score_threshold=threshold
                )

                st.session_state.rg4_results = results
                st.session_state.integrated_results = None
                st.session_state.exon_mapping = None
                st.session_state.mapping_status = None

                st.success(
                    f"Screening completed. {len(results)} "
                    "candidate intervals identified."
                )

            except Exception as exc:
                st.error(f"rG4 analysis failed: {exc}")

        results = st.session_state.rg4_results

        if results is not None:

            if results.empty:

                st.warning("No candidate regions detected.")

            else:

                g4hunter = results[
                    results["candidate_type"] == "g4hunter"
                ]

                canonical = results[
                    results["candidate_type"] == "canonical_motif"
                ]

                c1, c2, c3 = st.columns(3)

                c1.metric("Total candidate intervals", len(results))
                c2.metric("G4Hunter regions", len(g4hunter))
                c3.metric("Canonical motifs", len(canonical))

                st.dataframe(
                    results,
                    use_container_width=True,
                    hide_index=True
                )

                st.subheader("G4Hunter score distribution")

                if not g4hunter.empty:

                    fig = px.histogram(
                        g4hunter,
                        x="score",
                        nbins=15,
                        labels={
                            "score": "Mean G4Hunter window score",
                            "count": "Candidate count"
                        }
                    )

                    st.plotly_chart(
                        fig,
                        use_container_width=True
                    )

                else:

                    st.info(
                        "No G4Hunter sliding-window candidates were detected."
                    )

                st.subheader("Canonical motif results")

                if not canonical.empty:

                    display_columns = [
                        column for column in [
                            "sequence_id",
                            "start",
                            "end",
                            "candidate_sequence",
                            "score",
                            "score_type",
                            "length"
                        ]
                        if column in canonical.columns
                    ]

                    st.dataframe(
                        canonical[display_columns],
                        use_container_width=True,
                        hide_index=True
                    )

                else:

                    st.info("No canonical motifs detected.")

                st.download_button(
                    "Download rG4 results CSV",
                    results.to_csv(index=False).encode("utf-8"),
                    file_name="RNA_Insight_rg4_candidates.csv",
                    mime="text/csv"
                )

        st.info(
            "Predicted candidates are not experimentally confirmed "
            "RNA G-quadruplex structures."
        )


# ==================================================
# 03 | AUTOMATIC APA ANALYSIS
# ==================================================

elif page == "APA site analysis":

    st.header("03 | Alternative polyadenylation site analysis")

    st.write(
        "Retrieve human polyadenylation site annotations "
        "automatically from PolyA_DB v4.1."
    )

    st.info(
        "Reference genome: hg38 | Species: Homo sapiens | "
        "Source: PolyA_DB v4.1"
    )

    gene_panel = st.session_state.gene_panel

    if gene_panel is not None and not gene_panel["genes"].empty:

        available_genes = (
            gene_panel["genes"]["gene_symbol"]
            .dropna()
            .astype(str)
            .str.upper()
            .unique()
            .tolist()
        )

    else:

        available_genes = DEFAULT_GENES

    selected_genes = st.multiselect(
        "Select genes for APA annotation retrieval",
        options=available_genes,
        default=available_genes
    )

    if st.button("Fetch APA annotations automatically"):

        if not selected_genes:

            st.warning("Select at least one gene.")

        else:

            with st.spinner(
                "Downloading and processing PolyA_DB annotations..."
            ):

                try:

                    apa = fetch_apa_sites(
                        gene_symbols=selected_genes
                    )

                    st.session_state.apa_results = apa
                    st.session_state.integrated_results = None

                    st.success(
                        f"Retrieved {len(apa)} APA site records."
                    )

                except Exception as exc:

                    st.error(
                        f"Automatic APA retrieval failed: {exc}"
                    )

    apa = st.session_state.apa_results

    if apa is not None:

        st.divider()
        st.subheader("APA annotation results")

        summary = summarize_apa_sites(apa)

        c1, c2, c3 = st.columns(3)

        c1.metric(
            "Total APA sites",
            summary["total_sites"]
        )

        c2.metric(
            "Genes represented",
            summary["sequence_count"]
        )

        c3.metric(
            "Other site annotations",
            summary["other_sites"]
        )

        st.dataframe(
            apa,
            use_container_width=True,
            hide_index=True
        )

        st.download_button(
            "Download APA annotations CSV",
            apa.to_csv(index=False).encode("utf-8"),
            file_name="RNA_Insight_apa_annotations.csv",
            mime="text/csv"
        )

        st.caption(
            "The current PolyA_DB retrieval returns gene-level "
            "genomic annotations, not transcript-specific assignments."
        )

    st.caption(
        "APA annotation counts do not measure differential APA usage."
    )


# ==================================================
# 04 | AUTOMATIC INTEGRATED ANALYSIS
# ==================================================

elif page == "Integrated rG4–APA analysis":

    st.header("04 | Integrated rG4–APA analysis")

    rg4 = st.session_state.rg4_results
    apa = st.session_state.apa_results
    gene_panel = st.session_state.gene_panel

    if rg4 is None or rg4.empty:

        st.info("Run rG4 screening first.")

    elif apa is None or apa.empty:

        st.info("Retrieve APA annotations first.")

    elif gene_panel is None:

        st.info(
            "Retrieve human gene and transcript data first. "
            "Ensembl transcript metadata is required for mapping."
        )

    else:

        st.markdown("### Automatic transcript-to-genome mapping")

        transcripts = gene_panel["transcripts"].copy()

        candidate_ids = set(
            rg4["sequence_id"]
            .apply(normalize_transcript_id)
            .tolist()
        )

        transcripts["transcript_id"] = (
            transcripts["transcript_id"].astype(str).str.strip()
        )

        matched_transcripts = transcripts[
            transcripts["transcript_id"].isin(candidate_ids)
            & transcripts["sequence_status"].eq("retrieved")
        ].copy()

        c1, c2 = st.columns(2)

        c1.metric(
            "rG4 transcript IDs",
            len(candidate_ids)
        )

        c2.metric(
            "Matching Ensembl transcripts",
            len(matched_transcripts)
        )

        unmatched = candidate_ids - set(
            matched_transcripts["transcript_id"]
        )

        if unmatched:

            st.warning(
                f"{len(unmatched)} candidate transcript IDs "
                "were not found in the retrieved Ensembl transcript panel."
            )

            with st.expander("View unmatched transcript IDs"):

                st.write(sorted(unmatched))

        if matched_transcripts.empty:

            st.error(
                "No matching transcript IDs were found. "
                "Use the FASTA downloaded directly from Human gene "
                "data retrieval and run rG4 screening on that FASTA."
            )

        else:

            st.caption(
                "Ensembl exon coordinates will be retrieved automatically. "
                "A mapping is accepted only when the summed exon length "
                "matches the retrieved cDNA length."
            )

            if st.button("Build exon mapping and integrate APA"):

                try:

                    with st.spinner(
                        "Retrieving and validating Ensembl exon coordinates..."
                    ):

                        mapping, mapping_status = (
                            build_transcript_exon_mapping(
                                matched_transcripts,
                                delay_seconds=0.15
                            )
                        )

                    st.session_state.exon_mapping = mapping
                    st.session_state.mapping_status = mapping_status

                    if mapping.empty:
                        raise ValueError(
                            "Ensembl did not return any validated exon mappings."
                        )

                    check_hg38_assembly(mapping)

                    validated_ids = set(
                        mapping_status.loc[
                            mapping_status["mapping_status"].eq("validated"),
                            "transcript_id"
                        ].astype(str)
                    )

                    validated_mapping = mapping[
                        mapping["sequence_id"].astype(str).isin(
                            validated_ids
                        )
                        & mapping["mapping_status"].eq("validated")
                    ].copy()

                    validated_rg4 = rg4[
                        rg4["sequence_id"]
                        .apply(normalize_transcript_id)
                        .isin(validated_ids)
                    ].copy()

                    if validated_rg4.empty:

                        raise ValueError(
                            "No rG4 candidates belong to transcripts "
                            "with validated exon mappings."
                        )

                    with st.spinner(
                        "Mapping rG4 candidates and comparing APA annotations..."
                    ):

                        integrated = integrate_gene_level_apa(
                            validated_rg4,
                            apa,
                            validated_mapping
                        )

                    st.session_state.integrated_results = integrated

                    st.success(
                        "Automatic mapping and gene-level APA comparison completed."
                    )

                except Exception as exc:

                    st.session_state.integrated_results = None

                    st.error(
                        f"Integrated analysis failed: {exc}"
                    )

            mapping_status = st.session_state.mapping_status

            if mapping_status is not None:

                st.divider()
                st.subheader("Ensembl exon mapping quality control")

                validated_count = int(
                    mapping_status["mapping_status"]
                    .eq("validated").sum()
                )

                failed_count = int(
                    mapping_status["mapping_status"]
                    .eq("failed").sum()
                )

                c1, c2 = st.columns(2)

                c1.metric(
                    "Validated transcripts",
                    validated_count
                )

                c2.metric(
                    "Failed mappings",
                    failed_count
                )

                st.dataframe(
                    mapping_status,
                    use_container_width=True,
                    hide_index=True
                )

                st.download_button(
                    "Download exon mapping status",
                    mapping_status.to_csv(index=False).encode("utf-8"),
                    file_name="RNA_Insight_exon_mapping_status.csv",
                    mime="text/csv"
                )

            integrated = st.session_state.integrated_results

            if integrated is not None:

                st.divider()
                st.subheader("Integrated results")

                if integrated.empty:

                    st.warning(
                        "No candidate regions could be integrated "
                        "using validated transcript mappings."
                    )

                else:

                    overlap_count = int(
                        integrated["overlaps_gene_level_apa"].sum()
                    )

                    c1, c2, c3 = st.columns(3)

                    c1.metric(
                        "Mapped rG4 candidates",
                        len(integrated)
                    )

                    c2.metric(
                        "Gene-level APA overlaps",
                        overlap_count
                    )

                    c3.metric(
                        "Genes represented",
                        integrated["gene_symbol"].nunique()
                    )

                    st.dataframe(
                        integrated,
                        use_container_width=True,
                        hide_index=True
                    )

                    st.download_button(
                        "Download integrated results CSV",
                        integrated.to_csv(index=False).encode("utf-8"),
                        file_name="RNA_Insight_integrated_results.csv",
                        mime="text/csv"
                    )

                    st.divider()
                    st.subheader("Genomic coordinate visualization")

                    chromosomes = sorted(
                        integrated["chrom"].astype(str).unique()
                    )

                    chromosome = st.selectbox(
                        "Select chromosome",
                        chromosomes
                    )

                    fig = create_genomic_figure(
                        integrated,
                        apa,
                        chromosome
                    )

                    st.plotly_chart(
                        fig,
                        use_container_width=True
                    )

                    st.warning(
                        "Interpretation: these results show gene-level "
                        "genomic overlap between predicted rG4 regions "
                        "and annotated PAS intervals. They do not establish "
                        "transcript-specific APA usage, differential APA "
                        "activity or a causal rG4-mediated mechanism."
                    )



# --------------------------------------------------
# 05 REPORT
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

    if integrated is None:
        integrated = pd.DataFrame()

    # Report summary metrics
    c1, c2, c3 = st.columns(3)

    c1.metric(
        "rG4 candidate intervals",
        len(rg4)
    )

    c2.metric(
        "APA annotation sites",
        len(apa)
    )

    c3.metric(
        "Integrated candidates",
        len(integrated)
    )

    # Generate report content
    report = build_report(
        rg4_candidates=rg4,
        apa_sites=apa,
        integrated_results=integrated,
    )

    # Generate PDF
    try:
        pdf_data = build_pdf_report(report)

        st.download_button(
            label="📄 Download research report (PDF)",
            data=pdf_data,
            file_name="RNA_Insight_Research_Report.pdf",
            mime="application/pdf",
            use_container_width=True,
        )

        st.success(
            "Your PDF research report is ready to download."
        )

    except Exception as exc:
        st.error(
            f"PDF generation failed: {exc}"
        )

    # Report preview
    st.divider()

    st.subheader("Report preview")

    st.markdown(report)
