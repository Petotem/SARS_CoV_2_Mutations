# **COVID-19 Primer Analysis Pipeline**

This repository contains a Python-based workflow to clean COVID-19 genome FASTA files, merge metadata, and analyze primer binding performance across multiple genomes.

It is divided into three main modules, each handling a key part of the analysis pipeline.

---
## **Project Structure**
covid19-primer-pipeline/
├── covid_fasta_clean.py # Step 1 - Clean FASTA and merge metadata
├── primer_finder.py # Step 2 - Analyze primers against cleaned genomes
├── primer_loader.py # Helper - Load and organize primer information
├── data/
│ ├── raw/
│ │ ├── covid_raw_sequence.fasta
│ │ ├── covid_meta_data.tsv
│ │ ├── covid_variant.xlsx
│ │ ├── primers_info.xlsx
│ └── clean/ # Generated cleaned sequences
└── results/
├── result.xlsx # Primer matching statistics
└── report_reason.xlsx # Detailed mismatch reasons

---

## **Installation**
```bash
pip install pandas numpy biopython xlsxwriter openpyxl
```
---

# **How to Run the Pipeline**

Step 1: Clean and Prepare the FASTA Data
Run the following command to clean your raw COVID-19 FASTA and metadata:
python covid_fasta_clean.py
Reads the raw .fasta file (covid_raw_sequence.fasta)

Merges metadata from covid_meta_data.tsv and lineage info from covid_variant.xlsx

Validates sequences by length and content

Writes clean, structured CSV files in the folder:
👉 data/clean/split_*.csv

Also saves summary statistics in:
👉 data/clean/statistic.csv

Step 2: Analyze Primer-Genome Matching
After cleaning the sequences, run:
python primer_finder.py

Loads primer data from data/raw/primers_info.xlsx

Aligns each primer against all cleaned genome sequences

Computes matching probabilities and temperature differences
Writes detailed results to:
results/result.xlsx → statistical summaries by country, lineage, date, etc.
results/report_reason.xlsx → mismatch reasons and failed matches
