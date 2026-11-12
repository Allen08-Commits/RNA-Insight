
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import re
import gzip
import io
import json
from datetime import datetime
from pathlib import Path

# ============================================================
# RNA-INSIGHT | UNIVERSAL RNA BIOINFORMATICS PLATFORM
# ============================================================

st.set_page_config(
    page_title="RNA-Insight | Molecular Research Lab",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================
# 1. FUTURISTIC WEBSITE DESIGN
# ============================================================

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@400;500&family=Manrope:wght@400;500;600;700;800&display=swap');

:root {
    --bg:#090d1b;
    --cyan:#35e7f2;
    --muted:#9aabc8;
}

html, body, [class*="css"] {
    font-family:'Manrope',sans-serif;
}

.stApp {
    background:
    radial-gradient(ellipse at 12% 4%, rgba(0,188,212,.16), transparent 32%),
    radial-gradient(ellipse at 88% 14%, rgba(126,70,210,.18), transparent 34%),
    linear-gradient(135deg,#080c19 0%,#10152b 52%,#17122c 100%);
    color:#eef5ff;
}

[data-testid="stHeader"] {
    background:rgba(7,11,24,.45);
}

[data-testid="stSidebar"] {
    background:linear-gradient(180deg,rgba(17,27,51,.98),rgba(24,20,48,.98));
    border-right:1px solid rgba(91,218,235,.18);
}

[data-testid="stSidebar"] * {
    color:#e8f4ff;
}

.block-container {
    padding-top:1.7rem;
    padding-bottom:3rem;
    max-width:1500px;
}

h1,h2,h3 {
    letter-spacing:-.035em;
}

h1 {
    font-weight:800 !important;
}

p,label,.stMarkdown {
    color:#dce8fa;
}

div[data-testid="stMetric"] {
    background:linear-gradient(145deg,rgba(31,48,78,.78),rgba(23,31,57,.72));
    border:1px solid rgba(87,218,235,.18);
    border-radius:16px;
    padding:15px 17px;
}

div[data-testid="stMetricLabel"] {
    color:#a9bddb;
}

div[data-testid="stMetricValue"] {
    color:#f4fbff;
}

div[data-testid="stVerticalBlockBorderWrapper"] {
    border:1px solid rgba(90,210,230,.18);
    border-radius:18px;
    background:rgba(17,26,49,.58);
}

.stButton>button,
.stDownloadButton>button {
    border-radius:12px;
    border:1px solid rgba(53,231,242,.48);
    background:linear-gradient(110deg,#0e8296,#2765a7);
    color:white;
    font-weight:700;
    padding:.55rem 1rem;
}

.stButton>button:hover,
.stDownloadButton>button:hover {
    border-color:#7af7ff;
    color:white;
    box-shadow:0 0 22px rgba(53,231,242,.2);
}

.stTextInput input,
.stTextArea textarea,
.stSelectbox div[data-baseweb="select"],
.stNumberInput input {
    background:#111a30 !important;
    color:#eef5ff !important;
}

[data-testid="stFileUploader"] {
    border:1px dashed rgba(53,231,242,.45);
    border-radius:16px;
    background:rgba(20,37,63,.42);
    padding:12px;
}

hr {
    border-color:rgba(113,174,210,.18);
}

.tag {
    display:inline-block;
    padding:5px 10px;
    border-radius:99px;
    background:rgba(44,219,231,.12);
    color:#65f3f6;
    border:1px solid rgba(44,219,231,.28);
    font-size:11px;
    letter-spacing:.12em;
    font-weight:800;
}

.hero {
    padding:28px 30px;
    border-radius:22px;
    border:1px solid rgba(74,222,239,.24);
    background:linear-gradient(110deg,rgba(25,58,87,.76),rgba(36,29,67,.72)),
    radial-gradient(circle at 75% 30%,rgba(34,213,234,.18),transparent 38%);
    box-shadow:0 18px 55px rgba(0,0,0,.22);
    position:relative;
    overflow:hidden;
}

.hero:after {
    content:"";
    position:absolute;
    width:310px;
    height:310px;
    right:-70px;
    top:-160px;
    border:1px solid rgba(67,225,242,.22);
    border-radius:50%;
    box-shadow:0 0 0 25px rgba(67,225,242,.025),
    0 0 0 55px rgba(67,225,242,.025),
    0 0 70px rgba(67,225,242,.08);
}

.kicker {
    color:#55eaf3;
    font-family:'DM Mono',monospace;
    letter-spacing:.18em;
    font-size:11px;
}

.hero h1 {
    font-size:clamp(32px,4vw,56px);
    line-height:1.02;
    margin:.45rem 0;
    color:#f5fbff;
}

.hero p {
    color:#b5c9e3;
    max-width:700px;
}

code {
    color:#63e9f2 !important;
}
</style>
""", unsafe_allow_html=True)


# ============================================================
# 2. UNIVERSAL DATA RECOGNITION
# ============================================================

ALIASES = {
    "gene_id": [
        "gene_id", "geneid", "gene", "gene_name",
        "gene_symbol", "symbol", "genes"
    ],

    "transcript_id": [
        "transcript_id", "transcriptid",
        "transcript", "tx_id", "mrna_id"
    ],

    "sequence": [
        "sequence", "seq", "rna_sequence",
        "dna_sequence", "nucleotide_sequence",
        "transcript_sequence", "utr_sequence"
    ],

    "chromosome": [
        "chromosome", "chrom", "chr", "seqname"
    ],

    "start": [
        "start", "start_position", "begin"
    ],

    "end": [
        "end", "end_position", "stop"
    ],

    "strand": [
        "strand", "orientation"
    ],

    "expression": [
        "expression", "tpm", "fpkm",
        "counts", "read_count", "abundance"
    ],

    "condition": [
        "condition", "treatment", "group",
        "phenotype", "experimental_condition"
    ],

    "sample_id": [
        "sample", "sample_id", "sampleid"
    ],

    "polyA_site": [
        "polya_site", "pas", "pas_position",
        "polyadenylation_site", "cleavage_site",
        "cleavage_position"
    ]
}


def normalize_column(value):
    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(value).strip().lower()
    )


def map_columns(columns):

    normalized = {
        normalize_column(c): c
        for c in columns
    }

    mapping = {}

    for canonical, aliases in ALIASES.items():

        for alias in aliases:

            if normalize_column(alias) in normalized:

                mapping[canonical] = normalized[
                    normalize_column(alias)
                ]

                break

    return mapping


# ============================================================
# 3. FILE FORMAT DETECTION
# ============================================================

def decompress_if_needed(filename, raw):

    if str(filename).lower().endswith(".gz"):
        return gzip.decompress(raw)

    return raw


def detect_format(filename, raw):

    filename = str(filename).lower()

    inner_name = (
        filename[:-3]
        if filename.endswith(".gz")
        else filename
    )

    extension = Path(inner_name).suffix.lower()

    formats = {
        ".csv": "CSV",
        ".tsv": "TSV",
        ".txt": "TEXT",
        ".xlsx": "XLSX",
        ".fa": "FASTA",
        ".fasta": "FASTA",
        ".fna": "FASTA",
        ".fas": "FASTA",
        ".faa": "FASTA",
        ".bed": "BED",
        ".gff": "GFF",
        ".gff3": "GFF3",
        ".gtf": "GTF"
    }

    if extension in formats:
        return formats[extension]

    sample = raw[:3000].decode(
        "utf-8",
        errors="ignore"
    ).lstrip("\ufeff \r\n")

    if sample.startswith(">"):
        return "FASTA"

    first_line = (
        sample.splitlines()[0]
        if sample.splitlines()
        else ""
    )

    if first_line.count("\t") >= 8:
        return "GFF3"

    if "\t" in first_line:
        return "TSV"

    if "," in first_line:
        return "CSV"

    return "UNKNOWN"


# ============================================================
# 4. FASTA PARSER
# ============================================================

def parse_fasta(raw):

    text = raw.decode(
        "utf-8",
        errors="replace"
    )

    records = []
    header = None
    sequence_lines = []

    for line in text.splitlines():

        line = line.strip()

        if not line:
            continue

        if line.startswith((";", "#", "!")) and header is None:
            continue

        if line.startswith(">"):

            if header is not None:

                sequence = "".join(
                    sequence_lines
                ).replace(" ", "").upper()

                records.append((
                    header.split()[0],
                    header,
                    sequence
                ))

            header = line[1:].strip()
            sequence_lines = []

        else:

            if header is None:

                raise ValueError(
                    "Sequence text found before the first FASTA header beginning with '>'."
                )

            if line.startswith((";", "#", "!")):
                continue

            sequence_lines.append(line)

    if header is not None:

        sequence = "".join(
            sequence_lines
        ).replace(" ", "").upper()

        records.append((
            header.split()[0],
            header,
            sequence
        ))

    if not records:
        raise ValueError("No FASTA records were found.")

    return pd.DataFrame(
        records,
        columns=[
            "record_id",
            "description",
            "sequence"
        ]
    ).assign(
        sequence_length=lambda d: d["sequence"].str.len()
    )


# ============================================================
# 5. GFF / GTF PARSER
# ============================================================

def parse_gff(raw, file_format):

    df = pd.read_csv(
        io.BytesIO(raw),
        sep="\t",
        comment="#",
        header=None,
        dtype=str
    )

    if df.shape[1] < 9:

        raise ValueError(
            f"{file_format} should contain 9 tab-separated columns. "
            f"Found {df.shape[1]}."
        )

    df = df.iloc[:, :9]

    df.columns = [
        "chromosome",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    def parse_attributes(value):

        output = {}

        for item in str(value).strip().strip(";").split(";"):

            if not item.strip():
                continue

            if "=" in item:

                key, val = item.split("=", 1)

                output[key.strip()] = val.strip().strip('"')

            elif " " in item:

                key, val = item.split(" ", 1)

                output[key.strip()] = val.strip().strip('"')

        return output

    attributes = df["attributes"].apply(parse_attributes)

    df["gene_id"] = attributes.apply(
        lambda x: x.get("gene_id")
        or x.get("gene")
        or x.get("Parent", "")
    )

    df["transcript_id"] = attributes.apply(
        lambda x: x.get("transcript_id")
        or x.get("transcript")
        or x.get("ID", "")
    )

    df["start"] = pd.to_numeric(
        df["start"],
        errors="coerce"
    )

    df["end"] = pd.to_numeric(
        df["end"],
        errors="coerce"
    )

    return df


# ============================================================
# 6. BED PARSER
# ============================================================

def parse_bed(raw):

    df = pd.read_csv(
        io.BytesIO(raw),
        sep="\t",
        comment="#",
        header=None,
        dtype=str
    )

    if df.shape[1] < 3:

        raise ValueError(
            "BED requires at least 3 columns: chromosome, start and end."
        )

    columns = [
        "chromosome",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thick_start",
        "thick_end",
        "item_rgb",
        "block_count",
        "block_sizes",
        "block_starts"
    ]

    df.columns = columns[:df.shape[1]]

    df["start"] = pd.to_numeric(
        df["start"],
        errors="coerce"
    )

    df["end"] = pd.to_numeric(
        df["end"],
        errors="coerce"
    )

    return df


# ============================================================
# 7. UNIVERSAL DATASET LOADER
# ============================================================

def load_dataset(upload):

    raw = decompress_if_needed(
        upload.name,
        upload.getvalue()
    )

    file_format = detect_format(
        upload.name,
        raw
    )

    if file_format == "FASTA":

        df = parse_fasta(raw)

    elif file_format in ("CSV", "TSV", "TEXT"):

        separator = (
            "," if file_format == "CSV"
            else "\t" if file_format == "TSV"
            else None
        )

        df = pd.read_csv(
            io.BytesIO(raw),
            sep=separator,
            engine="python"
        )

    elif file_format == "XLSX":

        df = pd.read_excel(
            io.BytesIO(raw)
        )

    elif file_format == "BED":

        df = parse_bed(raw)

    elif file_format in ("GFF", "GFF3", "GTF"):

        df = parse_gff(
            raw,
            file_format
        )

    else:

        raise ValueError(
            "File format could not be recognized. "
            "Try CSV, TSV, XLSX, FASTA, BED, GFF3 or GTF."
        )

    return df, file_format


# ============================================================
# 8. SEQUENCE CLASSIFICATION AND VALIDATION
# ============================================================

def classify_sequence(sequence):

    sequence = str(sequence).upper()

    if not sequence:
        return "Empty"

    valid = set("ACGTUNRYKMSWBDHV")

    if not set(sequence) <= valid:
        return "Protein / non-nucleotide"

    if "T" in sequence and "U" in sequence:
        return "Mixed DNA/RNA (check)"

    if "U" in sequence:
        return "RNA"

    if "T" in sequence:
        return "DNA"

    return "Ambiguous nucleotide"


def sequence_stats(sequence):

    sequence = re.sub(
        r"\s+",
        "",
        str(sequence).upper()
    )

    length = len(sequence)

    counts = {
        base: sequence.count(base)
        for base in "ACGTU"
    }

    gc = (
        (counts["G"] + counts["C"]) / length * 100
        if length
        else 0
    )

    return {
        "length": length,
        "A": counts["A"],
        "C": counts["C"],
        "G": counts["G"],
        "T": counts["T"],
        "U": counts["U"],
        "GC_percent": gc
    }


# ============================================================
# 9. DATA QUALITY REPORT
# ============================================================

def qc_report(df, mapping):

    report = {
        "Rows": len(df),
        "Columns": len(df.columns),
        "Duplicate rows": int(df.duplicated().sum()),
        "Missing cells": int(df.isna().sum().sum()),
        "Recognized fields": len(mapping)
    }

    if "sequence" in mapping:

        column = mapping["sequence"]

        report["Empty sequences"] = int(
            df[column]
            .fillna("")
            .astype(str)
            .str.strip()
            .eq("")
            .sum()
        )

    if "start" in mapping and "end" in mapping:

        start = pd.to_numeric(
            df[mapping["start"]],
            errors="coerce"
        )

        end = pd.to_numeric(
            df[mapping["end"]],
            errors="coerce"
        )

        report["Invalid coordinates"] = int(
            (
                start.isna()
                | end.isna()
                | (start < 0)
                | (end < start)
            ).sum()
        )

    return report


# ============================================================
# 10. RNA G-QUADRUPLEX MOTIF SCREENING
# ============================================================

def find_rg4(
    sequence,
    min_loop=1,
    max_loop=7,
    min_g=3
):

    sequence = str(sequence).upper().replace("T", "U")

    pattern = re.compile(
        rf"G{{{min_g},}}[ACGUN]{{{min_loop},{max_loop}}}"
        rf"G{{{min_g},}}[ACGUN]{{{min_loop},{max_loop}}}"
        rf"G{{{min_g},}}[ACGUN]{{{min_loop},{max_loop}}}"
        rf"G{{{min_g},}}"
    )

    results = []

    for match in pattern.finditer(sequence):

        motif = match.group()

        results.append({
            "start": match.start() + 1,
            "end": match.end(),
            "motif": motif,
            "length": len(motif),
            "score": round(
                motif.count("G") / len(motif) * 100,
                1
            )
        })

    return results


def make_csv(df):

    return df.to_csv(
        index=False
    ).encode("utf-8")


# ============================================================
# 11. SESSION STATE
# ============================================================

if "dataset" not in st.session_state:
    st.session_state.dataset = None

if "filename" not in st.session_state:
    st.session_state.filename = None

if "file_format" not in st.session_state:
    st.session_state.file_format = None


# ============================================================
# 12. SIDEBAR NAVIGATION
# ============================================================

with st.sidebar:

    st.markdown("""
    <div style='padding:8px 4px 18px'>
        <div class='kicker'>MOLECULAR RESEARCH PLATFORM</div>
        <h2 style='margin:5px 0;color:#f5fbff'>
            🧬 RNA-INSIGHT
        </h2>
        <div style='color:#91a9ca;font-size:12px'>
            Universal RNA Bioinformatics
        </div>
    </div>
    """, unsafe_allow_html=True)

    page = st.radio(
        "NAVIGATION",
        [
            "Overview",
            "Dataset Workspace",
            "Data Intelligence",
            "Sequence Analysis",
            "RNA G-Quadruplex",
            "APA Analysis",
            "Integrated Research",
            "Results Centre",
            "About & Validation"
        ],
        label_visibility="collapsed"
    )

    st.markdown("---")

    st.markdown(
        "<span class='tag'>SYSTEM READY</span>",
        unsafe_allow_html=True
    )

    if st.session_state.dataset is not None:

        st.caption(
            f"Loaded: {st.session_state.filename}"
        )

        st.caption(
            f"{len(st.session_state.dataset):,} records · "
            f"{st.session_state.file_format}"
        )

    else:

        st.caption("No dataset loaded")


# ============================================================
# 13. COMMON PAGE FUNCTIONS
# ============================================================

def page_header(kicker, title, subtitle):

    st.markdown(
        f"""
        <div class='kicker'>{kicker}</div>
        <h1 style='margin:.2rem 0 .3rem;color:#f4f8ff'>
            {title}
        </h1>
        <p style='color:#a9bdd8;margin-top:0'>
            {subtitle}
        </p>
        """,
        unsafe_allow_html=True
    )


def require_data():

    if st.session_state.dataset is None:

        st.warning(
            "Upload a dataset in Dataset Workspace to unlock this module."
        )

        st.stop()

    return st.session_state.dataset


# ============================================================
# 14. OVERVIEW DASHBOARD
# ============================================================

if page == "Overview":

    st.markdown("""
    <div class="hero">

        <div class="kicker">
            UNIVERSAL RNA BIOINFORMATICS ANALYSIS PLATFORM
        </div>

        <h1>
            Decode the language<br>of RNA.
        </h1>

        <p>
            Explore nucleotide sequences, screen for RNA
            G-quadruplex candidates and investigate alternative
            polyadenylation through one integrated research workspace.
        </p>

        <span class="tag">
            SEQUENCE • rG4 • APA • INTEGRATION
        </span>

    </div>
    """, unsafe_allow_html=True)

    st.write("")

    a, b, c, d = st.columns(4)

    a.metric(
        "Dataset status",
        "Loaded" if st.session_state.dataset is not None
        else "Awaiting upload"
    )

    b.metric(
        "Records",
        f"{len(st.session_state.dataset):,}"
        if st.session_state.dataset is not None
        else "—"
    )

    c.metric("Input formats", "8+")

    d.metric("Analysis modules", "7")

    st.write("")

    left, right = st.columns([1.25, 1])

    with left:

        with st.container(border=True):

            st.markdown("### Research workspace")

            st.write(
                "Start by uploading a compatible biological dataset. "
                "RNA-Insight will inspect its format, infer the data "
                "structure and identify analyses supported by the available fields."
            )

            st.markdown("**Supported formats**")

            st.markdown(
                "`CSV` `TSV` `XLSX` `FASTA` `BED` "
                "`GFF` `GFF3` `GTF` and supported `.gz` files."
            )

    with right:

        with st.container(border=True):

            st.markdown("### Analysis pipeline")

            pipeline = [
                (
                    "01 · Ingest",
                    "Recognize and parse biological files"
                ),
                (
                    "02 · Validate",
                    "Inspect columns, sequences and coordinates"
                ),
                (
                    "03 · Analyze",
                    "Run compatible sequence, rG4 and APA workflows"
                ),
                (
                    "04 · Interpret",
                    "Visualize and export reproducible results"
                )
            ]

            for i, (name, description) in enumerate(pipeline):

                st.markdown(f"**{name}**")

                st.caption(description)

                if i < 3:

                    st.markdown(
                        "<div style='color:#35e7f2;padding-left:8px'>↓</div>",
                        unsafe_allow_html=True
                    )

    st.caption(
        "Research note: computational motif screening and coordinate-based "
        "APA comparisons are exploratory analyses, not experimental validation."
    )


# ============================================================
# 15. DATASET WORKSPACE
# ============================================================

elif page == "Dataset Workspace":

    page_header(
        "01 / DATA INGESTION",
        "Dataset Workspace",
        "Upload a research file and let the platform inspect its format and biological structure."
    )

    uploaded = st.file_uploader(
        "DROP YOUR BIOLOGICAL DATASET HERE",
        type=None,
        help="Supported: CSV, TSV, XLSX, FASTA, BED, GFF, GFF3, GTF and supported gzip-compressed variants."
    )

    if uploaded:

        try:

            df, file_format = load_dataset(uploaded)

            st.session_state.dataset = df

            st.session_state.filename = uploaded.name

            st.session_state.file_format = file_format

            st.success(
                f"Successfully parsed {uploaded.name}"
            )

            mapping = map_columns(df.columns)

            x, y, z = st.columns(3)

            x.metric("Detected format", file_format)

            y.metric("Records", f"{len(df):,}")

            z.metric("Columns", len(df.columns))

            st.markdown(
                "### Automatic biological field recognition"
            )

            if mapping:

                mapping_df = pd.DataFrame([
                    {
                        "Standard field": key,
                        "Detected column": value
                    }
                    for key, value in mapping.items()
                ])

                st.dataframe(
                    mapping_df,
                    use_container_width=True,
                    hide_index=True
                )

            else:

                st.info(
                    "No standard biological fields were confidently recognized. "
                    "You can still inspect the table."
                )

            st.markdown("### Dataset preview")

            st.dataframe(
                df.head(20),
                use_container_width=True
            )

            st.download_button(
                "Download parsed dataset as CSV",
                make_csv(df),
                file_name="rna_insight_parsed_dataset.csv",
                mime="text/csv"
            )

        except Exception as error:

            st.error(
                f"Dataset processing failed: {error}"
            )

    elif st.session_state.dataset is not None:

        st.info(
            f"Current dataset: {st.session_state.filename}"
        )

        st.dataframe(
            st.session_state.dataset.head(20),
            use_container_width=True
        )

        if st.button("Clear current dataset"):

            st.session_state.dataset = None
            st.session_state.filename = None
            st.session_state.file_format = None

            st.rerun()

    else:

        st.markdown("""
        <div class="hero" style="margin-top:12px">

            <div class="kicker">INPUT GATEWAY</div>

            <h2 style="color:white">
                Bring your data into the lab.
            </h2>

            <p>
                Upload sequence files, annotation tables or genomic intervals.
                File extension and content are both inspected where possible.
            </p>

        </div>
        """, unsafe_allow_html=True)


# ============================================================
# 16. DATA INTELLIGENCE AND QUALITY CONTROL
# ============================================================

elif page == "Data Intelligence":

    page_header(
        "02 / QUALITY CONTROL",
        "Data Intelligence",
        "Inspect data structure, recognized fields, missing values and potential formatting problems."
    )

    df = require_data()

    mapping = map_columns(df.columns)

    quality = qc_report(df, mapping)

    columns = st.columns(4)

    for (label, value), column in zip(
        quality.items(),
        columns * 2
    ):

        column.metric(
            label,
            f"{value:,}"
            if isinstance(value, (int, float))
            else value
        )

    st.markdown("### Standard field mapping")

    if mapping:

        mapping_df = pd.DataFrame([
            {
                "Canonical field": key,
                "Input column": value
            }
            for key, value in mapping.items()
        ])

        st.dataframe(
            mapping_df,
            use_container_width=True,
            hide_index=True
        )

    else:

        st.warning(
            "No standard fields were detected automatically."
        )

    st.markdown("### Missing values by column")

    missing = (
        df.isna()
        .sum()
        .rename("Missing values")
        .reset_index()
        .rename(columns={"index": "Column"})
    )

    st.dataframe(
        missing,
        use_container_width=True,
        hide_index=True
    )

    if "sequence" in mapping:

        sequences = (
            df[mapping["sequence"]]
            .fillna("")
            .astype(str)
        )

        types = (
            sequences.apply(classify_sequence)
            .value_counts()
            .rename_axis("Sequence type")
            .reset_index(name="Count")
        )

        st.markdown("### Sequence type inspection")

        st.dataframe(
            types,
            use_container_width=True,
            hide_index=True
        )

    st.markdown("### Cleaning options")

    remove_duplicates = st.checkbox(
        "Remove exact duplicate rows"
    )

    remove_empty = st.checkbox(
        "Remove rows that are completely empty"
    )

    if st.button("Apply selected cleaning"):

        cleaned = df.copy()

        if remove_duplicates:
            cleaned = cleaned.drop_duplicates()

        if remove_empty:
            cleaned = cleaned.dropna(how="all")

        st.session_state.dataset = cleaned.reset_index(drop=True)

        st.success(
            f"Working dataset updated: {len(cleaned):,} rows remain."
        )


# ============================================================
# 17. SEQUENCE ANALYSIS
# ============================================================

elif page == "Sequence Analysis":

    page_header(
        "03 / SEQUENCE INTELLIGENCE",
        "Sequence Analysis",
        "Calculate nucleotide composition, sequence length and basic sequence quality metrics."
    )

    df = require_data()

    mapping = map_columns(df.columns)

    sequence_column = mapping.get("sequence")

    if not sequence_column:

        st.warning(
            "A sequence column was not recognized. "
            "Upload a FASTA file or a table containing a sequence column."
        )

        st.stop()

    work = df.copy()

    work["_seq"] = (
        work[sequence_column]
        .fillna("")
        .astype(str)
    )

    work["_type"] = work["_seq"].apply(
        classify_sequence
    )

    work["_length"] = work["_seq"].str.len()

    work["_gc"] = work["_seq"].apply(
        lambda sequence: sequence_stats(sequence)["GC_percent"]
    )

    m1, m2, m3 = st.columns(3)

    m1.metric(
        "Sequences",
        len(work)
    )

    m2.metric(
        "Mean length",
        f"{work['_length'].mean():,.1f} nt"
    )

    m3.metric(
        "Mean GC content",
        f"{work['_gc'].mean():.2f}%"
    )

    c1, c2 = st.columns(2)

    with c1:

        fig = px.histogram(
            work,
            x="_length",
            nbins=30,
            title="Sequence length distribution",
            template="plotly_dark"
        )

        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    with c2:

        fig = px.histogram(
            work,
            x="_gc",
            nbins=30,
            title="GC content distribution",
            template="plotly_dark"
        )

        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    st.markdown("### Per-sequence metrics")

    identifier_column = (
        mapping.get("transcript_id")
        or mapping.get("gene_id")
    )

    if identifier_column:

        identifiers = work[identifier_column].astype(str)

    else:

        identifiers = work.index.astype(str)

    result = pd.DataFrame({
        "Record": identifiers,
        "Sequence type": work["_type"],
        "Length": work["_length"],
        "GC %": work["_gc"]
    })

    st.dataframe(
        result,
        use_container_width=True,
        hide_index=True
    )

    st.download_button(
        "Export sequence metrics",
        make_csv(result),
        file_name="sequence_metrics.csv",
        mime="text/csv"
    )


# ============================================================
# 18. RNA G-QUADRUPLEX ANALYSIS
# ============================================================

elif page == "RNA G-Quadruplex":

    page_header(
        "04 / STRUCTURAL MOTIF SCREENING",
        "RNA G-Quadruplex Analysis",
        "Screen RNA sequences for canonical G-rich sequence patterns that may form G-quadruplex structures."
    )

    df = require_data()

    mapping = map_columns(df.columns)

    sequence_column = mapping.get("sequence")

    if not sequence_column:

        st.warning(
            "A sequence column is required for rG4 motif screening."
        )

        st.stop()

    col1, col2, col3 = st.columns(3)

    min_g = col1.slider(
        "Minimum G-tract length",
        2, 5, 3
    )

    min_loop = col2.slider(
        "Minimum loop length",
        1, 3, 1
    )

    max_loop = col3.slider(
        "Maximum loop length",
        2, 12, 7
    )

    if min_loop > max_loop:

        st.error(
            "Minimum loop length cannot exceed maximum."
        )

        st.stop()

    rows = []

    identifier_column = (
        mapping.get("transcript_id")
        or mapping.get("gene_id")
    )

    for index, row in df.iterrows():

        record_id = (
            str(row[identifier_column])
            if identifier_column
            else str(index + 1)
        )

        sequence = (
            str(row[sequence_column])
            if pd.notna(row[sequence_column])
            else ""
        )

        if classify_sequence(sequence) not in (
            "RNA",
            "DNA",
            "Ambiguous nucleotide"
        ):

            continue

        hits = find_rg4(
            sequence,
            min_loop,
            max_loop,
            min_g
        )

        for hit in hits:

            rows.append({
                "Record": record_id,
                "Start (1-based)": hit["start"],
                "End": hit["end"],
                "Motif": hit["motif"],
                "Motif length": hit["length"],
                "G content score (%)": hit["score"]
            })

    hits = pd.DataFrame(rows)

    a, b, c = st.columns(3)

    a.metric(
        "Sequences screened",
        len(df)
    )

    b.metric(
        "Candidate motifs",
        len(hits)
    )

    c.metric(
        "Sequences with candidates",
        hits["Record"].nunique()
        if len(hits)
        else 0
    )

    st.info(
        "These are canonical sequence-pattern candidates only. "
        "This is not a thermodynamic, machine-learning or experimental "
        "confirmation of an RNA G-quadruplex."
    )

    if len(hits):

        st.dataframe(
            hits,
            use_container_width=True,
            hide_index=True
        )

        st.download_button(
            "Download rG4 candidate table",
            make_csv(hits),
            file_name="rg4_candidates.csv",
            mime="text/csv"
        )

        fig = px.histogram(
            hits,
            x="Motif length",
            title="Candidate motif length distribution",
            template="plotly_dark"
        )

        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    else:

        st.warning(
            "No canonical motifs were detected with the current settings."
        )


# ============================================================
# 19. ALTERNATIVE POLYADENYLATION ANALYSIS
# ============================================================

elif page == "APA Analysis":

    page_header(
        "05 / TRANSCRIPT END ANALYSIS",
        "Alternative Polyadenylation",
        "Compare transcript end coordinates within genes when compatible genomic annotation is available."
    )

    df = require_data()

    mapping = map_columns(df.columns)

    required = [
        "gene_id",
        "transcript_id",
        "chromosome",
        "start",
        "end",
        "strand"
    ]

    missing = [
        field
        for field in required
        if field not in mapping
    ]

    if missing:

        st.warning(
            "This dataset does not contain all required fields "
            "for coordinate-based transcript-end comparison: "
            + ", ".join(missing)
        )

        st.info(
            "Upload a compatible transcript annotation table "
            "with gene, transcript, chromosome, start, end and strand fields. "
            "A poly(A)-site annotation can provide stronger evidence."
        )

        st.stop()

    work = pd.DataFrame({

        "gene_id": df[mapping["gene_id"]].astype(str),

        "transcript_id": df[mapping["transcript_id"]].astype(str),

        "chromosome": df[mapping["chromosome"]].astype(str),

        "start": pd.to_numeric(
            df[mapping["start"]],
            errors="coerce"
        ),

        "end": pd.to_numeric(
            df[mapping["end"]],
            errors="coerce"
        ),

        "strand": df[mapping["strand"]].astype(str)

    }).dropna(
        subset=["start", "end"]
    )

    work = work[
        (work["start"] >= 0)
        & (work["end"] >= work["start"])
        & work["strand"].isin(["+", "-"])
    ]

    work["transcript_end"] = np.where(
        work["strand"] == "+",
        work["end"],
        work["start"]
    )

    work["transcript_length"] = (
        work["end"] - work["start"] + 1
    )

    st.metric(
        "Valid transcript records",
        len(work)
    )

    st.caption(
        "Coordinate comparison is meaningful only when records share "
        "a compatible genome assembly and annotation convention."
    )

    apa_results = []

    for gene, group in work.groupby(
        ["gene_id", "chromosome", "strand"]
    ):

        unique = group.drop_duplicates(
            "transcript_id"
        ).copy()

        if len(unique) < 2:
            continue

        # Transcript orientation:
        # Plus strand: smaller genomic end is proximal.
        # Minus strand: larger genomic start is proximal.

        strand = unique["strand"].iloc[0]

        ordered = unique.sort_values(
            "transcript_end",
            ascending=(strand == "-"),
            kind="stable"
        )

        for rank, (_, row) in enumerate(
            ordered.iterrows()
        ):

            if rank == 0:
                label = "Proximal"

            elif rank == len(ordered) - 1:
                label = "Distal"

            else:
                label = "Intermediate"

            apa_results.append({

                "gene_id": gene[0],

                "chromosome": gene[1],

                "strand": gene[2],

                "transcript_id": row["transcript_id"],

                "transcript_end": row["transcript_end"],

                "relative_end_group": label

            })

    apa_df = pd.DataFrame(apa_results)

    if len(apa_df):

        st.success(
            f"Compared transcript ends for "
            f"{apa_df['gene_id'].nunique()} gene groups."
        )

        st.dataframe(
            apa_df,
            use_container_width=True,
            hide_index=True
        )

        st.download_button(
            "Download transcript-end comparison",
            make_csv(apa_df),
            file_name="apa_transcript_end_comparison.csv",
            mime="text/csv"
        )

        fig = px.histogram(
            apa_df,
            x="relative_end_group",
            color="strand",
            title="Relative transcript-end groups",
            template="plotly_dark"
        )

        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)"
        )

        st.plotly_chart(
            fig,
            use_container_width=True
        )

    else:

        st.warning(
            "No gene groups with two or more distinct transcripts "
            "were available for comparison."
        )

    st.warning(
        "This is an exploratory transcript-end comparison, "
        "not definitive APA-site calling. It does not independently "
        "establish cleavage/polyadenylation sites. Use curated PAS/APA "
        "annotations and a consistent genome build for stronger APA inference."
    )


# ============================================================
# 20. INTEGRATED rG4–APA RESEARCH
# ============================================================

elif page == "Integrated Research":

    page_header(
        "06 / CROSS-MODULE RESEARCH",
        "Integrated rG4–APA Research",
        "Connect sequence-level rG4 candidates with transcript-end and APA-related evidence."
    )

    df = require_data()

    mapping = map_columns(df.columns)

    if "sequence" not in mapping:

        st.warning(
            "A sequence column is required for integrated rG4 screening."
        )

        st.stop()

    st.markdown("### Sequence candidate screening")

    records = []

    identifier_column = (
        mapping.get("transcript_id")
        or mapping.get("gene_id")
    )

    for index, row in df.iterrows():

        record_id = (
            str(row[identifier_column])
            if identifier_column
            else str(index + 1)
        )

        sequence = (
            str(row[mapping["sequence"]])
            if pd.notna(row[mapping["sequence"]])
            else ""
        )

        for hit in find_rg4(sequence):

            records.append({

                "Record": record_id,

                "rG4_start": hit["start"],

                "rG4_end": hit["end"],

                "motif": hit["motif"]

            })

    rg4_df = pd.DataFrame(records)

    st.metric(
        "rG4 sequence candidates",
        len(rg4_df)
    )

    if len(rg4_df):

        st.dataframe(
            rg4_df,
            use_container_width=True,
            hide_index=True
        )

    else:

        st.info(
            "No canonical rG4 sequence candidates detected using default settings."
        )

    st.markdown("### Integration status")

    has_apa = all(
        field in mapping
        for field in [
            "gene_id",
            "transcript_id",
            "start",
            "end",
            "strand",
            "chromosome"
        ]
    )

    if has_apa:

        st.success(
            "The dataset contains fields that may support transcript-end comparison. "
            "Open APA Analysis for coordinate-based comparison."
        )

    else:

        st.warning(
            "The uploaded dataset does not contain the full coordinate "
            "annotation needed for transcript-end comparison."
        )

    st.info(
        "A true rG4–APA overlap analysis requires transcript sequences "
        "or genomic coordinates mapped to the same transcript/genome reference, "
        "plus supported APA/PAS annotations. This version does not infer "
        "that overlap from unrelated columns."
    )


# ============================================================
# 21. RESULTS CENTRE
# ============================================================

elif page == "Results Centre":

    page_header(
        "07 / EXPORT & REPORTING",
        "Results Centre",
        "Export the current dataset and analysis-ready tables."
    )

    df = require_data()

    st.markdown("### Current dataset")

    st.metric(
        "Records available",
        len(df)
    )

    st.download_button(
        "Download current dataset (CSV)",
        make_csv(df),
        file_name="rna_insight_dataset.csv",
        mime="text/csv"
    )

    if st.session_state.filename:

        summary = {

            "platform": "RNA-Insight",

            "source_file": st.session_state.filename,

            "detected_format": st.session_state.file_format,

            "records": len(df),

            "columns": [
                str(column)
                for column in df.columns
            ],

            "generated_at": datetime.now().isoformat(
                timespec="seconds"
            ),

            "recognized_fields": map_columns(df.columns)

        }

        st.download_button(
            "Download dataset QC metadata (JSON)",
            json.dumps(
                summary,
                indent=2
            ).encode(),

            file_name="rna_insight_qc_metadata.json",

            mime="application/json"
        )

    st.caption(
        "PDF report generation and persistent project history "
        "are planned extensions. This version exports validated tables and metadata."
    )


# ============================================================
# 22. ABOUT AND SCIENTIFIC VALIDATION
# ============================================================

elif page == "About & Validation":

    page_header(
        "PLATFORM INFORMATION",
        "About RNA-Insight",
        "A modular research interface for RNA sequence and transcript annotation analysis."
    )

    with st.container(border=True):

        st.markdown("### Current analytical capabilities")

        st.markdown("""
        - Multi-format file ingestion and basic content recognition.
        - Automatic mapping of common biological column names.
        - Missing-value, duplicate and coordinate quality checks.
        - DNA/RNA sequence statistics and motif screening.
        - Canonical G-rich rG4 sequence-pattern screening.
        - Exploratory strand-aware transcript-end comparisons.
        - CSV and JSON exports.
        """)

    with st.container(border=True):

        st.markdown("### Scientific limitations")

        st.markdown("""
        1. A sequence motif is a candidate pattern, not proof that an rG4 structure forms.

        2. Transcript-end differences alone do not prove alternative polyadenylation.

        3. APA comparisons require compatible genome assembly,
           strand-aware coordinates and preferably curated PAS annotations.

        4. The platform does not currently run wet-lab validation,
           thermodynamic folding, or external specialist predictors.

        5. Expression values alone cannot establish rG4 or APA events.
        """)

    with st.container(border=True):

        st.markdown("### Software information")

        st.write(
            "Built with Python, Streamlit, Pandas, NumPy and Plotly."
        )

        st.write(
            "Version: 1.0 research prototype"
        )

        st.write(
            "Generated:",
            datetime.now().strftime("%d %b %Y, %H:%M")
        )
