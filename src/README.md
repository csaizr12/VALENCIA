## Table of contents
- [VALENCIA - Validating Annotation Levels Employing Nucleotide Comparison Identity Analysis](#valencia---validating-annotation-levels-employing-nucleotide-comparison-identity-analysis)
- [Requirements](#requirements)
    - [Python dependencies](#python-dependencies)
    - [Software dependencies](#software-dependencies)
- [Installation](#installation)
- [Usage](#usage)
  - [Generate quality panel](#generate-quality-panel)

# VALENCIA - Validating Annotation Levels Employing Nucleotide Comparison Identity Analysis
The main objective of this work is the design, development, and implementation of VALENCIA, a tool aimed at evaluating the quality of genome annotations by using experimental data such as transcripts and protein sequences.

## Requirements
### Python dependencies
- Python == 3.13.12
- pandas
- matplotlib
- seaborn
- numpy

### Software dependencies

- conda == 26.6.1
- Bioconda: agat == 1.6.1
- GFFread == 0.12.7 (https://github.com/gpertea/gffread)
- GFFcompare == 0.12.10 (https://github.com/gpertea/gffcompare)

[Back to Table of contents](#table-of-contents)

## Installation
Follow the steps below to set up the environment required to run VALENCIA.

Create and activate the conda environment:

```bash
conda create -n valencia_env python=3.13.12
conda activate valencia_env
```

Install Python dependencies:
```bash
pip install -r requirements.txt
```
Install external tools:
```bash
conda install -c bioconda agat gffread gffcompare
```
Verify installation
```bash
gffread --version
gffcompare --version
agat --help
```
[Back to Table of contents](#table-of-contents)

## Usage

VALENCIA evaluates genome annotations by integrating structural comparison and sequence similarity analysis.

To run the pipeline, four input files are required:

- Target annotation (-t)

- Evidence annotation (-p)

- Genome assembly (-g)

- Reference exon annotation (-x)

- Output directory (-o)

VALENCIA automatically extracts transcripts and proteins from both annotations and performs structural and sequence-level comparisons.

Example command

The following example illustrates a typical execution of VALENCIA using Arabidopsis thaliana datasets:

```bash
python3 VALENCIA/VALENCIA.py \
    -t INPUTS/Arabidopsis_thaliana/Artha_AllRNASeq.STAR.TAIR10.EVT_STv1.gff3 \
    -p INPUTS/Arabidopsis_thaliana/Arabidopsis_thaliana_GeMoMa_with_Oryza_sativa.gff\
    -g INPUTS/Arabidopsis_thaliana/Athaliana_447_TAIR10.fa \
    -x INPUTS/Arabidopsis_thaliana/Athaliana_447_Araport11.gene_exons.gff3 \
    -o test_valencia_araport1
```
### Parameter description

- [-t] Target annotation to be evaluated (GFF/GTF)

- [-p] Evidence annotation used for comparison (GFF/GTF)

- [-g] Genome assembly in FASTA format

- [-x] Reference exon annotation used for structural validation

- [-o] Output directory where all results will be stored

### What VALENCIA does
Running VALENCIA performs the following steps:

1. Sequence extraction
    
    Using GFFread, VALENCIA extracts:

    - Transcript  sequences

    - CDS sequences

    - Protein sequences

    These are stored in:
    ```bash
    target_sequences/
    evidence_sequences/
    ```
2. Structural comparison (GFFcompare)


