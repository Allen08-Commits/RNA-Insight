# RNA-Insight

### In-Silico Analysis of RNA G-Quadruplex-Mediated
### Alternative Polyadenylation in Human Stress-Response Genes

RNA-Insight is a Streamlit research application for screening RNA
sequences for rG4 candidate regions and integrating those candidates
with genomic alternative polyadenylation annotations.

## Features

- FASTA sequence upload and validation
- DNA-to-RNA sequence normalization
- Canonical rG4 motif screening
- G/C composition-window heuristic
- APA annotation standardization
- Genomic coordinate mapping
- rG4–APA overlap analysis
- Downloadable CSV results
- Downloadable Markdown research report

## Installation

Install the dependencies:

    pip install -r requirements.txt

Run the application:

    streamlit run app.py

## Input requirements

### FASTA

Upload RNA or DNA FASTA sequences.

### APA annotation table

Required columns:

- sequence_id
- chrom
- start
- end
- strand

Optional column:

- site_type

### Sequence mapping table

Required columns:

- sequence_id
- chrom
- start
- end
- strand

Coordinates must be 0-based and half-open.

Use a consistent genome assembly and coordinate system across
the mapping and APA annotation files.

## Scientific limitations

The application identifies computational candidate regions.
It does not experimentally confirm rG4 folding.

Genomic overlap does not demonstrate causal regulation.

APA annotation counts do not establish differential APA usage.
Differential usage requires sample-level measurements, biological
replicates and appropriate statistical analysis.

The included example data are synthetic demonstration records.

## Project structure

    RNA-Insight/
    ├── app.py
    ├── requirements.txt
    ├── README.md
    ├── .gitignore
    ├── core/
    ├── tests/
    └── sample_data/
