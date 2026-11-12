
import streamlit as st
import pandas as pd
import numpy as np
import csv
import gzip
import re
import io
from pathlib import Path
from Bio import SeqIO

# ==========================================
# 1. APPLICATION CONFIGURATION
# ==========================================

st.set_page_config(
    page_title="RNA-Insight",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 RNA-Insight")
st.caption("Universal RNA Bioinformatics Analysis Platform")

st.markdown("""
Upload biological datasets in different formats.
RNA-Insight automatically detects their structure,
identifies biological information and validates
the uploaded records.
""")

# ==========================================
# 2. COLUMN RECOGNITION DICTIONARY
# ==========================================

COLUMN_SYNONYMS = {
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
        "polyA_site", "pas", "pas_position",
        "polyadenylation_site", "cleavage_site",
        "cleavage_position"
    ]
}


def normalize_column(name):
    name = str(name).strip().lower()
    name = re.sub(r"[^a-z0-9]+", "_", name)
    return name.strip("_")


def recognize_columns(columns):
    mapping = {}

    for column in columns:
        normalized = normalize_column(column)

        for standard_name, synonyms in COLUMN_SYNONYMS.items():
            if normalized in synonyms:
                if standard_name not in mapping:
                    mapping[standard_name] = column
                break

    return mapping


# ==========================================
# 3. FILE FORMAT DETECTION
# ==========================================

def detect_format(filename, content):

    name = filename.lower()

    if name.endswith(".gz"):
        name = name[:-3]

    extension = Path(name).suffix.lower()

    formats = {
        ".csv": "CSV",
        ".tsv": "TSV",
        ".xlsx": "XLSX",
        ".fasta": "FASTA",
        ".fa": "FASTA",
        ".fna": "FASTA",
        ".fas": "FASTA",
        ".faa": "FASTA",
        ".bed": "BED",
        ".gff": "GFF3/GTF",
        ".gff3": "GFF3/GTF",
        ".gtf": "GFF3/GTF"
    }

    if extension in formats:
        return formats[extension]

    # Inspect the actual file content when extension
    # is missing or not recognized.

    try:
        sample = content[:10000].decode(
            "utf-8-sig", errors="replace"
        )

        lines = [
            line.strip() for line in sample.splitlines()
            if line.strip()
        ]

        if not lines:
            return "UNKNOWN"

        if lines[0].startswith(">"):
            return "FASTA"

        if "##gff-version" in sample[:1000]:
            return "GFF3/GTF"

        first = lines[0]

        if len(first.split("\t")) >= 9:
            return "GFF3/GTF"

        if "\t" in first:
            return "TSV"

        if "," in first:
            return "CSV"

    except Exception:
        pass

    return "UNKNOWN"


# ==========================================
# 4. FILE DECOMPRESSION
# ==========================================

def prepare_content(filename, content):

    if filename.lower().endswith(".gz"):
        content = gzip.decompress(content)

    return content


# ==========================================
# 5. DATASET PARSERS
# ==========================================

def parse_fasta(content):

    text = content.decode("utf-8-sig")

    records = list(
        SeqIO.parse(io.StringIO(text), "fasta")
    )

    if not records:
        raise ValueError("No valid FASTA records found.")

    data = []

    for record in records:
        sequence = str(record.seq).upper()

        data.append({
            "record_id": record.id,
            "description": record.description,
            "sequence": sequence,
            "sequence_length": len(sequence)
        })

    return pd.DataFrame(data)


def parse_tabular(content, file_format):

    text = content.decode("utf-8-sig")

    if file_format == "XLSX":
        return pd.read_excel(io.BytesIO(content))

    if file_format == "CSV":
        return pd.read_csv(io.StringIO(text))

    if file_format == "TSV":
        return pd.read_csv(
            io.StringIO(text), sep="\t"
        )

    # Automatic delimiter detection
    try:
        delimiter = csv.Sniffer().sniff(
            text[:10000],
            delimiters=",\t;"
        ).delimiter
    except csv.Error:
        delimiter = ","

    return pd.read_csv(
        io.StringIO(text),
        sep=delimiter,
        engine="python"
    )


def parse_bed(content):

    names = [
        "chromosome", "start", "end", "name",
        "score", "strand", "thick_start",
        "thick_end", "item_rgb", "block_count",
        "block_sizes", "block_starts"
    ]

    df = pd.read_csv(
        io.StringIO(content.decode("utf-8-sig")),
        sep="\t",
        comment="#",
        header=None
    )

    if df.empty or df.shape[1] < 3:
        raise ValueError("Invalid BED file.")

    if df.shape[1] > len(names):
        raise ValueError("Unsupported BED column count.")

    df.columns = names[:df.shape[1]]

    return df


def parse_annotation(content):

    names = [
        "chromosome", "source", "feature",
        "start", "end", "score", "strand",
        "frame", "attributes"
    ]

    df = pd.read_csv(
        io.StringIO(content.decode("utf-8-sig")),
        sep="\t",
        comment="#",
        header=None
    )

    if df.empty or df.shape[1] != 9:
        raise ValueError(
            "GFF3/GTF requires nine tab-separated fields."
        )

    df.columns = names

    return df


# ==========================================
# 6. BIOLOGICAL DATA CLASSIFICATION
# ==========================================

def classify_sequence(sequence):

    sequence = str(sequence).upper()

    dna = set("ACGTN")
    rna = set("ACGUN")
    protein = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

    if not sequence:
        return "Empty"

    if not set(sequence).issubset(dna | rna):
        if set(sequence).issubset(protein):
            return "Protein"
        return "Invalid/Unknown"

    if "U" in sequence and "T" in sequence:
        return "Mixed DNA/RNA"

    if "U" in sequence:
        return "RNA"

    if "T" in sequence:
        return "DNA"

    return "Nucleotide (DNA/RNA ambiguous)"


def detect_biological_type(df, file_format):

    if file_format == "FASTA":
        if "sequence" not in df.columns:
            return "Sequence Dataset"

        types = df["sequence"].apply(classify_sequence)
        return "Sequence Dataset: " + ", ".join(
            sorted(types.unique())
        )

    if file_format == "BED":
        return "Genomic Interval Dataset"

    if file_format == "GFF3/GTF":
        return "Gene / Transcript Annotation Dataset"

    mapping = recognize_columns(df.columns)

    if "sequence" in mapping:
        return "Sequence Annotation Dataset"

    if (
        "polyA_site" in mapping
        or (
            "start" in mapping
            and "end" in mapping
            and "chromosome" in mapping
        )
    ):
        return "Genomic / Potential APA Annotation Dataset"

    if "expression" in mapping:
        return "Gene Expression Dataset"

    if "gene_id" in mapping:
        return "Gene Identifier Dataset"

    return "General Biological Tabular Dataset"


# ==========================================
# 7. QUALITY CONTROL
# ==========================================

def validate_dataset(df, file_format):

    report = {}

    report["Total Records"] = len(df)
    report["Total Columns"] = len(df.columns)
    report["Duplicate Rows"] = int(
        df.duplicated().sum()
    )

    report["Missing Values"] = int(
        df.isna().sum().sum()
    )

    if "sequence" in df.columns:

        sequences = df["sequence"].fillna("").astype(str)

        report["Empty Sequences"] = int(
            (sequences.str.len() == 0).sum()
        )

        report["Invalid Sequences"] = int(
            sequences.apply(
                lambda x: classify_sequence(x)
                in ["Invalid/Unknown", "Mixed DNA/RNA"]
            ).sum()
        )

    if file_format in ["BED", "GFF3/GTF"]:

        if "start" in df.columns and "end" in df.columns:

            starts = pd.to_numeric(
                df["start"], errors="coerce"
            )

            ends = pd.to_numeric(
                df["end"], errors="coerce"
            )

            report["Invalid Coordinates"] = int(
                (
                    starts.isna()
                    | ends.isna()
                    | (starts < 0)
                    | (ends < starts)
                ).sum()
            )

    return report


# ==========================================
# 8. APPLICATION NAVIGATION
# ==========================================

st.sidebar.title("RNA-Insight")

page = st.sidebar.radio(
    "Navigation",
    [
        "Home",
        "Dataset Upload",
        "Data Validation",
        "Sequence Analysis",
        "RNA G4 Analysis",
        "APA Analysis",
        "Integrated Analysis",
        "Results & Reports"
    ]
)

# Store results across navigation
if "dataset" not in st.session_state:
    st.session_state.dataset = None

if "dataset_info" not in st.session_state:
    st.session_state.dataset_info = None

# ==========================================
# 9. HOME
# ==========================================

if page == "Home":

    st.header("Welcome to RNA-Insight")

    st.write(
        "A universal and modular RNA bioinformatics "
        "analysis platform."
    )

    col1, col2, col3 = st.columns(3)

    col1.metric("Analysis Modules", "5+")
    col2.metric("Supported Formats", "7+")
    col3.metric("Platform", "Universal")

    st.info(
        "Begin by uploading your biological dataset."
    )

# ==========================================
# 10. UNIVERSAL DATA UPLOAD
# ==========================================

elif page == "Dataset Upload":

    st.header("Universal Dataset Upload")

    st.write(
        "Upload a biological dataset. "
        "The platform will inspect its format and content."
    )

    uploaded_file = st.file_uploader(
        "Choose your dataset",
        type=None,
        max_upload_size=100,
        help=(
            "Supported: CSV, TSV, XLSX, FASTA, "
            "BED, GFF3, GTF and gzip-compressed variants."
        )
    )

    if uploaded_file:

        try:

            raw_content = uploaded_file.getvalue()

            if not raw_content:
                raise ValueError("The uploaded file is empty.")

            content = prepare_content(
                uploaded_file.name,
                raw_content
            )

            file_format = detect_format(
                uploaded_file.name,
                content
            )

            st.success(
                f"Detected format: {file_format}"
            )

            # Select appropriate parser
            if file_format == "FASTA":
                df = parse_fasta(content)

            elif file_format == "BED":
                df = parse_bed(content)

            elif file_format == "GFF3/GTF":
                df = parse_annotation(content)

            elif file_format in [
                "CSV", "TSV", "XLSX"
            ]:
                df = parse_tabular(
                    content, file_format
                )

            else:
                raise ValueError(
                    "Unsupported or unrecognized file format."
                )

            if df.empty:
                raise ValueError(
                    "The uploaded file contains no data records."
                )

            # Automatic recognition
            mapping = recognize_columns(df.columns)

            biological_type = detect_biological_type(
                df, file_format
            )

            quality = validate_dataset(
                df, file_format
            )

            # Store data
            st.session_state.dataset = df

            st.session_state.dataset_info = {
                "filename": uploaded_file.name,
                "format": file_format,
                "biological_type": biological_type,
                "mapping": mapping,
                "quality": quality
            }

            st.divider()

            st.subheader("Dataset Recognition")

            col1, col2 = st.columns(2)

            col1.metric(
                "Records",
                f"{len(df):,}"
            )

            col2.metric(
                "Columns",
                f"{len(df.columns):,}"
            )

            st.write(
                "**Biological classification:**",
                biological_type
            )

            st.write(
                "**Detected column mapping:**"
            )

            if mapping:
                st.dataframe(
                    pd.DataFrame(
                        list(mapping.items()),
                        columns=[
                            "Recognized Field",
                            "Original Column"
                        ]
                    ),
                    use_container_width=True
                )
            else:
                st.warning(
                    "No recognized biological columns. "
                    "You can still inspect the dataset."
                )

            st.divider()

            st.subheader("Dataset Preview")

            st.dataframe(
                df.head(20),
                use_container_width=True
            )

            st.divider()

            st.subheader("Data Quality Report")

            for key, value in quality.items():

                st.write(f"**{key}:** {value}")

            st.success(
                "Dataset successfully uploaded and processed."
            )

        except Exception as error:

            st.error(
                f"Dataset processing failed: {error}"
            )

# ==========================================
# 11. DATA VALIDATION
# ==========================================

elif page == "Data Validation":

    st.header("Data Validation")

    if st.session_state.dataset is None:

        st.info("Upload a dataset first.")

    else:

        info = st.session_state.dataset_info

        st.write("File:", info["filename"])
        st.write("Format:", info["format"])
        st.write("Data type:", info["biological_type"])

        st.subheader("Quality Metrics")

        for key, value in info["quality"].items():
            st.write(f"**{key}:** {value}")

        st.subheader("Column Mapping")

        st.json(info["mapping"])

        st.subheader("Missing Values Per Column")

        missing = (
            st.session_state.dataset
            .isna()
            .sum()
            .reset_index()
        )

        missing.columns = ["Column", "Missing Values"]

        st.dataframe(
            missing,
            use_container_width=True
        )

# ==========================================
# 12. FUTURE ANALYSIS MODULES
# ==========================================

else:

    st.header(page)

    if st.session_state.dataset is None:

        st.warning(
            "Please upload a dataset before starting analysis."
        )

    else:

        st.success(
            "Dataset available for analysis."
        )

        st.info(
            "This module will be implemented in the "
            "next development stages."
        )

        st.write(
            "Recognized data type:",
            st.session_state.dataset_info["biological_type"]
        )

        st.dataframe(
            st.session_state.dataset.head(10),
            use_container_width=True
        )
