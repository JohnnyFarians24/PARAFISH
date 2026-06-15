# PARAFISH

**PARAFISH** is a scientific web platform for the analysis of fish parasite communities using amplicon sequencing (metabarcoding) data. It provides a complete, reproducible workflow running inside Docker — from raw FASTQ quality control through OTU/ASV clustering, taxonomic classification, and phylogenetic inference — with an interactive Plotly Dash interface accessible from any browser.

---

## Overview

The platform integrates the following tools and databases into a single unified interface:

| Stage | Tools |
|---|---|
| Quality Control | FastQC |
| Preprocessing / Merging | Trimmomatic, PEAR |
| Clustering / Denoising | VSEARCH (SINTAX), DADA2 |
| Taxonomic Classification | VSEARCH SINTAX, BLAST, BOLDigger3 |
| Reference Databases | SILVA 138.2, PR2 v5.1.1, NCBI/BOLD, custom |
| Phylogenetic Inference | MAFFT, FastTree, RAxML-NG |

All bioinformatics tools are bundled inside the Docker image via Conda. Users do not need to install any tool manually beyond Docker itself.

---

## Requirements

- [Docker](https://docs.docker.com/get-docker/) (version 20.10 or higher)
- [Docker Compose](https://docs.docker.com/compose/install/) (version 2.0 or higher)
- ~5 GB free disk space for the Docker image
- Additional disk space for reference databases (see below)

No Python environment, Conda installation, or bioinformatics tooling is required on the host machine.

---

## Quick Start

```bash
git clone https://github.com/JohnnyFarians24/PARAFISH.git
cd PARAFISH
cp .env.example .env
# Edit .env to set your database paths (see Configuration below)
docker compose up --build
```

Once the container is running, open your browser at:

```
http://localhost:8050
```

---

## Reference Databases

The reference databases are **not included in this repository** due to file size. You must download them separately and configure their paths in `.env`.

### SILVA 138.2 SSU NR99

```bash
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
```

### PR2 v5.1.1

```bash
wget https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.fasta.gz
wget https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.tax.gz
gunzip pr2_version_5.1.1_SSU_mothur.fasta.gz
gunzip pr2_version_5.1.1_SSU_mothur.tax.gz
```

Place the downloaded files in a directory of your choice (e.g. `~/parafish_dbs/silva/` and `~/parafish_dbs/pr2/`) and set the paths accordingly in your `.env` file.

---

## Configuration

Copy `.env.example` to `.env` and edit the paths to match your local database locations:

```bash
cp .env.example .env
```

The key variables to configure:

| Variable | Description |
|---|---|
| `DB_SILVA_PATH` | Full path to the SILVA FASTA file |
| `DB_PR2_PATH` | Full path to the PR2 mothur FASTA file |
| `DB_PR2_DADA2_PATH` | Full path to the PR2 DADA2 FASTA file (optional) |
| `RUNS_DIR` | Directory where analysis runs will be stored |
| `GLOBAL_DB_ROOT` | Root directory for all databases |
| `ENABLE_SILVA` | Enable SILVA classification (`True`/`False`) |
| `ENABLE_PR2` | Enable PR2 classification (`True`/`False`) |
| `ENABLE_BOLD` | Enable BOLD classification (`True`/`False`) |
| `ENABLE_NCBI` | Enable NCBI BLAST classification (`True`/`False`) |

BOLD and NCBI are disabled by default. BOLD requires credentials (`BOLD_USERNAME`, `BOLD_PASSWORD`) from a registered [BOLD Systems](https://boldsystems.org/) account.

---

## Directory Structure

```
PARAFISH/
├── app.py                  # Application entry point
├── core/                   # Config, models, storage, validation, utilities
├── ui/                     # Dash layout, pages, components, callbacks
├── databases/              # Database interface modules (SILVA, PR2, NCBI, BOLD, custom)
├── analysis/               # VSEARCH, DADA2, and analysis pipeline runners
├── preprocess/             # Trimmomatic and PEAR runners
├── qc/                     # FastQC runner and QC plot generation
├── phylo/                  # MAFFT, FastTree, RAxML-NG runners and phylogenetic plots
├── assets/                 # CSS and static assets for the Dash interface
├── scripts/                # Utility shell scripts
├── tests/                  # Test suite
├── Dockerfile
├── docker-compose.yml
├── environment.yml         # Conda environment specification (used inside Docker)
├── requirements.txt        # Python package dependencies
├── .env.example            # Environment variable template
└── docs/
    └── setup_docker.md     # Extended Docker setup notes
```

---

## Input Data

PARAFISH expects **paired-end Illumina FASTQ files** as primary input. Files can be uploaded directly through the web interface or placed in the `uploads/` directory before starting the container.

Supported formats: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`

A small test dataset (`parafish_test_fastq.zip`) is available in the repository for initial validation.

---

## Workflow

The platform guides the user through a modular pipeline:

1. **Project Setup** — name the project and select input samples
2. **Quality Control** — run FastQC and inspect per-sample QC metrics
3. **Preprocessing** — adapter trimming (Trimmomatic) and paired-end merging (PEAR)
4. **OTU/ASV Clustering** — VSEARCH SINTAX or DADA2
5. **Taxonomic Classification** — against SILVA, PR2, NCBI, BOLD, or a custom database
6. **Phylogenetic Analysis** — multiple sequence alignment (MAFFT) and tree inference (FastTree or RAxML-NG)
7. **Results and Visualisation** — interactive plots and downloadable outputs

Each step is independent. Users can enter the workflow at any stage if upstream results are already available.

---

## Stopping and Restarting

```bash
# Stop the container
docker compose down

# Restart without rebuilding
docker compose up

# Rebuild after code changes
docker compose up --build
```

Run data and project state persist in the `RUNS_DIR` directory defined in `.env`, so restarting the container does not lose analysis progress.

---

## Troubleshooting

**Container fails to start**
Check that Docker is running and that the paths in `.env` exist and are accessible.

**Database not found errors**
Verify that `DB_SILVA_PATH` and `DB_PR2_PATH` in `.env` point to the actual files (not directories) and that the files are uncompressed.

**Port 8050 already in use**
Change the port mapping in `docker-compose.yml` from `"8050:8050"` to e.g. `"8051:8050"` and access the app at `http://localhost:8051`.

**Pipeline does not advance after sample selection**
This is a known issue under active development. Check the browser console for errors and report via GitHub Issues.

---

## Development Status

PARAFISH is under active development. The current release supports the full pipeline described above, with the following known limitations:

- BOLD parameter visibility in the UI has a CSS rendering issue (fix in progress)
- Output visualisation within the app is not yet fully functional; outputs remain accessible as files in `RUNS_DIR`
- Custom database support is not yet available in this release
- Not all pipeline outputs are currently visible within the app interface; full in-app result rendering is under development







