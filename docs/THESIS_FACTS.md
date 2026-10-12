# PARAFISH — Thesis fact sheet

Reference document for fact-checking the MSc dissertation. Compiled on **2026-10-12** from the
working tree at commit `25281bf` (branch `main`, in sync with `origin/main`; no code change since
`90dc89e`, 2026-08-03).

## How to read this document

| Tag | Meaning |
|---|---|
| **[CODE]** | Verified by reading a file in the working tree on 2026-10-12: source code, a config file, or an output file of a stored run. The source in parentheses says which. |
| **[DOC]** | Stated only in `README.md`, `DEVLOG.md`, a code comment or a docstring. Not independently verified. |
| **[INFER]** | My deduction from [CODE]/[DOC] facts. Not stated anywhere. |
| **[UNKNOWN]** | Could not be determined from the repository. |

Rules followed: a [DOC] or [INFER] claim is never presented as [CODE]. Where code and documentation
disagree, both are recorded and the row is marked **CONTRADICTION**; all of them are collected in
[Section 12](#12-contradictions-and-unknowns-index).

**Three things to know before citing anything:**

1. **The stored runs were not produced in Docker.** They were produced on the host (WSL2, 3.0 GB RAM,
   12 cores) with tool versions that differ from the ones the Dockerfile pins (see §1.3 and §2.0).
   Cite the host versions for results, the Docker versions for the distributed platform.
2. **The test suite could not be run today.** The project virtualenv points at an interpreter that no
   longer exists (`.venv/bin/python -> /home/joaofaria/miniconda/bin/python`, missing). The last
   recorded result is 581 passed / 2 skipped on 2026-08-03. See §5.3.
3. **Line numbers are those of commit `25281bf`.** They will drift if the code changes.

---

## 1. Architecture

### 1.1 Stack

| Component | Declared requirement | Version installed in the dev virtualenv | Source |
|---|---|---|---|
| Python | 3.11 (Docker) | 3.13.9 (dev venv) | [CODE] `environment.yml`; `.venv/pyvenv.cfg` |
| Dash | `>=2.17.0` | 4.1.0 | [CODE] `requirements.txt`; `.venv/lib/python3.13/site-packages/dash-4.1.0.dist-info` |
| dash-bootstrap-components | `>=1.6.0` | 2.0.4 | [CODE] same |
| Flask | not declared (Dash dependency) | 3.1.3 | [CODE] `.venv` dist-info; used directly in `app.py:36-37` |
| pandas | `>=2.2.0` | 3.0.3 | [CODE] |
| NumPy | not declared (transitive) | 2.4.6 | [CODE] |
| SciPy | not declared (transitive) | 1.18.0 | [CODE] used in `analysis/ordination.py`, `databases/database_metrics.py:123` |
| scikit-bio | `>=0.7.0` | 0.7.3 | [CODE] |
| Plotly | `>=5.20.0` | 6.7.0 | [CODE] |
| Biopython | **not declared** | 1.87 | [CODE] imported in `ui/callbacks_progress.py:12` and `databases/ncbi_lineage.py:123`; absent from `requirements.txt`. [INFER] it arrives transitively (boldigger3). |
| boldigger3 | `>=3.0.0,<4.0.0` | 3.0.0 | [CODE] |
| cutadapt | `>=4.0` | 5.2 | [CODE] `requirements.txt`; version from `runs/_real/current/preprocess/aptl64b_02/cutadapt.log` |
| pyarrow / openpyxl | `>=16.0.0` / `>=3.1.0` | 24.0.0 / 3.1.5 | [CODE] used to read BOLDigger outputs (`databases/taxonomy_runner.py:99-126`) |
| pydantic, python-dotenv | `>=2.7.0`, `>=1.0.0` | 2.13.4, 1.2.2 | [CODE] declared; no import of either found in `core/`, `ui/`, `databases/`, `analysis/` (`core/config.py:94-95` has its own `.env` loader). [INFER] unused. |
| duckdb | not declared | 1.5.3 | [CODE] present in venv; used by BOLDigger3's local store |

There is **no `pyproject.toml`**; dependencies are declared only as lower bounds in `requirements.txt`
and as exact pins for the bioinformatics tools in `environment.yml`. [CODE]

The web application is a single-process Dash app on a Flask server, started by `python app.py` with
`debug=True, threaded=True, use_reloader=False` ([CODE] `app.py:171-176`). No WSGI server
(gunicorn/uwsgi) is configured. [CODE] `Dockerfile` `CMD ["python", "app.py"]`.

### 1.2 Repository structure

| Path | Content | Tracked in git? |
|---|---|---|
| `app.py` | Entry point: creates the Dash app, Flask download routes, startup reconciliation of orphan runs (183 lines) | yes |
| `core/` | `config.py` (env vars, tool paths), `models.py` (dataclasses, defaults), `preflight.py` (memory checks), `run_registry.py` (run state from disk), `export.py` (ZIP + PROVENANCE), `process.py` (`run_tool`, failure signatures), `progress_writer.py`, `storage.py`, `validation.py`, `utils.py` (14 files, 2 330 lines) | yes |
| `qc/` | FastQC runner, parser, metrics (4 files) | yes |
| `preprocess/` | `cutadapt_runner.py`, `trimmomatic_runner.py`, `pear_runner.py`, `preprocess_metrics.py` | yes |
| `analysis/` | `vsearch_runner.py`, `dada2_runner.py` (embeds the R script), `analysis_metrics.py`, `ordination.py`, `community_analysis.py`, `feature_set_comparison.py`, `custom_db_compare.py`, `sample_parsing.py` | yes |
| `databases/` | `taxonomy_runner.py` (SINTAX, consensus, reconciliation; 1 754 lines), `silva_db.py`, `pr2_db.py`, `ncbi_db.py`, `ncbi_lineage.py`, `bold_db.py`, `custom_db.py`, `reference_filter.py`, `database_metrics.py` | yes |
| `phylo/` | `mafft_runner.py`, `fasttree_runner.py`, `raxmlng_runner.py`, `phylo_runner.py`, `phylo_metrics.py`, `phylo_diversity.py` | yes |
| `ui/` | `layout.py`, `pages/` (10 page layouts), `components/`, 12 `callbacks_*.py`, `plots.py`, `runtime_state.py`, `feature_sets.py` (44 files, ~12 200 lines) | yes |
| `scripts/` | `run_pipeline.py` (headless), `golden.py`, `verify.py`, `make_fixture.py`, `run_coi_bold_demo.py` | yes |
| `tests/` | 53 test modules | yes |
| `golden/fixture/`, `golden/fixture_asv/` | Frozen stage contracts for the two fixture datasets (64 tracked files) | yes |
| `golden/real/` | Frozen contracts for the full-depth 18S run | **no** — 0 tracked files |
| `golden_baseline/` | Snapshot of run `20260716T080902Z` taken before the fix plan | [UNKNOWN] not checked |
| `docs/` | `setup_docker.md`, this file | yes |
| `assets/` | `parafish_theme.css` and static assets | yes |
| `runs/`, `projects/`, `fixtures/`, `uploads/` | Run outputs, project store, fixture FASTQ, upload dir | **no** (`.gitignore:12`, `:51`) |
| `dados stora teresa/` | Novogene metabarcoding data from the supervisor (untracked, 356 MB) | no |
| repo root `*.fasta` | SILVA and PR2 references and their derived caches (~7 GB) | no |

[CODE] `git ls-files`, `.gitignore`, `wc -l`.

### 1.3 Docker

[CODE] `docker-compose.yml`, `Dockerfile`, `environment.yml`.

| Item | Value |
|---|---|
| Compose services | **one**: `parafish` (image `parafish:latest`, container `parafish_app`, `restart: unless-stopped`) |
| Compose file version | `"3.9"` |
| Port | `8050:8050` |
| `env_file` | `.env` (required) |
| Inline environment | `PYTHONDONTWRITEBYTECODE=1`, `PYTHONUNBUFFERED=1` |
| Base image | `continuumio/miniconda3:24.1.2-0` |
| Volumes | `${DB_ROOT:-./data/dbs}/silva → /data/silva:ro`, `…/pr2 → /data/pr2:ro`, `…/bold_db → /data/bold_db:ro`, `…/blast_db → /data/blast_db:ro`, `${RUNS_DIR:-./data/runs} → /data/runs:rw`, `${CONFIG_DIR:-./data/config} → /data/config:ro` |
| apt packages | `wget curl git unzip build-essential default-jre-headless fasttree libarchive-tools` |
| conda (`environment.yml`) | `python=3.11`, `fastqc=0.12.1`, `trimmomatic=0.39`, `mafft=7.526`, `raxml-ng=1.2.2`, `blast=2.15.0`, `fasttree=2.1.11`, `pear=0.9.11`, `r-base=4.3.3`, `bioconductor-dada2=1.30.0` |
| VSEARCH | static binary **2.28.1** downloaded from GitHub (`ARG VSEARCH_VERSION=2.28.1`) |
| pip | everything in `requirements.txt` (includes cutadapt and boldigger3) |
| Entry point | `CMD ["python", "app.py"]`, `EXPOSE 8050` |

**The Docker build has never been verified end to end.** [CODE] Docker is not installed in the WSL
distribution used for development (`docker: command not found`); [DOC] `README.md` "Installation"
states the same. [UNKNOWN] whether the image builds and whether the pinned versions resolve.

### 1.4 Execution model

| Fact | Source |
|---|---|
| "Run Pipeline" on the Setup page only **configures** the run (builds a `RunConfig`, stores it in `store-run-config`) and navigates to Progress. It does not start anything. | [CODE] `ui/callbacks_setup.py:602-611` |
| The run starts only on "Start Pipeline" on the Progress page, which calls `_launch_pipeline_run`. | [CODE] `ui/callbacks_progress.py:1032-1040` |
| The pipeline runs in **one background daemon thread** of the web process (`threading.Thread(target=_run_pipeline_thread, daemon=True)`). No task queue, no separate worker process. | [CODE] `ui/callbacks_progress.py:1011-1017` |
| Only one run at a time: the launcher refuses if any pipeline thread is alive, checked under a lock. | [CODE] `ui/callbacks_progress.py:994-1018` |
| Samples are processed **sequentially** in a `for` loop. The Setup "Workers" field is stored (`vsearch_params["workers"]`) and never read. | [CODE] `ui/callbacks_progress.py:190,252`; `ui/callbacks_setup.py:582` (only occurrence) |
| Stage order: QC → Preprocess → Analysis → Databases → Phylo. | [CODE] `ui/callbacks_progress.py:181-935`; `scripts/run_pipeline.py:39` |
| State is written to `<run_dir>/progress.json` by `core/progress_writer.py` (thread lock, read-merge-write of JSON). | [CODE] `core/progress_writer.py:12-41` |
| The Progress page polls every **2 000 ms** (`dcc.Interval`). | [CODE] `ui/pages/page_progress.py:13` |
| Every page rebuilds its view from `progress.json` on disk, not from the browser session. | [CODE] `core/run_registry.py:179-286`; [DOC] `DEVLOG.md` "UI Eixo B" entries (2026-07-17) |
| At startup, runs whose `progress.json` still says `running` with no live thread are marked `interrupted` (not `done`). | [CODE] `app.py:160-162`; `core/run_registry.py:136-171` |
| The same code path is driven headless by `scripts/run_pipeline.py`, which calls `_run_pipeline_thread` synchronously. | [CODE] `scripts/run_pipeline.py:1-22` (docstring) and body |
| Each stage page also has its own "Run <stage>" button with a separate callback implementation. | [CODE] `btn-run-qc`, `btn-run-preprocess`, `btn-run-analysis`, `btn-run-databases`, `btn-run-phylo` in `ui/pages/*.py`; `ui/callbacks_databases.py:1260-1370` |

**Run directory layout** ([CODE] observed in `runs/_real/current`, 127 files):

```
<RUNS_DIR>/<project>/<run_id>/            run_id = UTC timestamp slug, e.g. 20260723T072541Z
  progress.json                           state, per-stage records, log lines, run_config block
  input/                                  uploaded FASTQ (written by the Samples page)
  qc/<sample>/*_fastqc.{html,zip}         qc/metrics/qc_{overview,per_base_quality,length_distribution,module_status}.csv
  preprocess/<sample>/                    *_R{1,2}_primer.fastq.gz, cutadapt.{json,log},
                                          *_R{1,2}_{trimmed,unpaired}.fastq.gz, trimlog.txt, pear/*.assembled.fastq
  preprocess/metrics/                     read_tracking.csv, assembled_length_distribution.csv
  analysis/<sample>/filtered.fasta        per-sample quality-filtered reads
  analysis/pooled/                        pooled_filtered.fasta, derep.fasta, otus_raw.fasta, otus.fasta,
                                          chimeras.fasta, chimera.uchime, clusters.uc, map.uc,
                                          otu_table.tsv, otu_id_map.csv, step{2..5}_*.log
  analysis/otu_matrix.csv                 CSV copy of otu_table.tsv
  analysis/dada2/                         asv_sequences.fasta, asv_abundance.csv, asv_table.{tsv,rds},
                                          asv_id_map.csv, read_tracking.csv, taxonomy.csv, dada2_script.R
  analysis/metrics/<otu|asv>/             alpha_diversity, rarefaction, rarefaction_depth, chimera_stats,
                                          rank_abundance, feature_occupancy, feature_prevalence, ...
  analysis/comparison/                    OTU-vs-ASV tables (only when both sets exist)
  databases/<otu|asv>/                    taxonomy_consensus.csv, taxonomy_summary.csv, taxonomy_hits.csv,
                                          source_assignments.csv, source_status.csv, conflict_report.csv,
                                          annotated.fasta, table.biom, taxa_tables/, phyloseq/,
                                          silva/, pr2/, ncbi/, bold/, custom_blast/, metrics/
  phylo/<otu|asv>/mafft/                  sanitized.fasta, aligned.fasta, reversed_asvs.csv, mafft.log
  phylo/<otu|asv>/fasttree/tree.nwk
  phylo/<otu|asv>/raxml/                  raxml.raxml.{bestTree,support,bootstraps,bestModel,log,...}
```

### 1.5 Components that access the network

| Component | Destination | Purpose | Default | Source |
|---|---|---|---|---|
| `blastn -db nt -remote` | NCBI BLAST servers | Alignment of all representatives against `nt`, one job for the whole multi-FASTA | on when "NCBI" is ticked in Setup (unticked by default) | [CODE] `databases/ncbi_db.py:164-171` |
| `Bio.Entrez.efetch(db="taxonomy")` | NCBI E-utilities | Resolve BLAST `staxids` to ranked lineages | only if `ENTREZ_EMAIL` is set | [CODE] `databases/ncbi_lineage.py:113-146` |
| `boldigger3 identify` | BOLD v5 identification engine | COI identification | on when "BOLD" is ticked (ticked by default) | [CODE] `databases/bold_db.py:304-333,532-600` |
| `boldigger3 download_db` | BOLD | One-off download of the local metadata DuckDB; needs credentials | manual setup step; no caller in `ui/` | [CODE] `databases/bold_db.py:625-653` |
| `update_blastdb.pl --decompress nt` | NCBI FTP | Download a local `nt` | manual setup step; no caller in `ui/` | [CODE] `databases/ncbi_db.py:270-279` |
| `urllib.request.urlretrieve` | arbitrary URL | Generic download helper `_download_file` | [CODE] defined in `core/utils.py:207-210`; no caller found in `ui/`, `databases/`, `analysis/`, `scripts/`. `ALLOW_DB_AUTO_DOWNLOAD` is read in `core/config.py:244` and referenced nowhere else. [INFER] automatic database download is **not implemented**. |
| Docker build | conda channels, GitHub (VSEARCH), PyPI | image construction | — | [CODE] `Dockerfile` |

Everything else (FastQC, cutadapt, Trimmomatic, PEAR, VSEARCH incl. SINTAX, DADA2, MAFFT, FastTree,
RAxML-NG, custom BLAST) runs locally. [CODE]

---

## 2. Pipeline stage by stage

### 2.0 Tool versions: Docker pins versus the versions that produced the stored runs

| Tool | Docker pin ([CODE] `environment.yml` / `Dockerfile`) | Version that produced the stored runs | Source of the right-hand column |
|---|---|---|---|
| FastQC | 0.12.1 | **0.11.5** | [CODE] `##FastQC 0.11.5` in `runs/_real/current/qc/aptl64b_02/APTL64B_R1_fastqc.zip` |
| cutadapt | `>=4.0` (pip) | **5.2** | [CODE] `runs/_real/current/preprocess/aptl64b_02/cutadapt.log` |
| Trimmomatic | 0.39 | 0.39 | [CODE] `.env` `TRIMMOMATIC_JAR_PATH=…/Trimmomatic-0.39/trimmomatic-0.39.jar`; the log does not print a version |
| PEAR | 0.9.11 | 0.9.11 | [CODE] `PEAR v0.9.11 [Nov 5, 2017]` in `runs/_real/current/preprocess/aptl64b_02/pear/pear.log` |
| VSEARCH | 2.28.1 | **2.30.0** | [CODE] `vsearch v2.30.0_linux_x86_64` in `runs/_real/current/analysis/pooled/step3_cluster.log` |
| R / DADA2 | 4.3.3 / 1.30.0 | **R 4.5.0 / dada2 1.36.0** on the host today | [CODE] `Rscript -e packageVersion("dada2")` today; [INFER] the same versions produced `runs/_fixture_asv` (not logged in the run) |
| BLAST+ | 2.15.0 | **2.17.0+** | [CODE] `.env` `BLASTN_PATH=/home/joaofaria/tools/ncbi-blast-2.17.0+/bin/blastn`. The system `/usr/bin/blastn` is 2.9.0 and cannot reach NCBI ([CODE] `databases/ncbi_db.py:232-239`). |
| MAFFT | 7.526 | **7.453** on the host today | [CODE] `mafft --version` today; [INFER] same for stored runs |
| FastTree | 2.1.11 | 2.1.11 (double precision, no SSE3) | [CODE] `runs/_real/current/phylo/otu/fasttree/fasttree.log` |
| RAxML-NG | 1.2.2 | 1.2.2-master (30.04.2024) | [CODE] `runs/_real/current/phylo/otu/raxml/raxml.raxml.log` |
| BOLDigger3 | `>=3.0.0,<4.0.0` | 3.0.0 | [CODE] `.venv` dist-info |

`java` is not on `PATH` on the host today ([CODE] `java: command not found`), so Trimmomatic could
not be run in the current shell. [INFER] it was available when the stored runs were made.

### 2.1 Implementation status of each named tool

| Tool | Status | Evidence |
|---|---|---|
| FastQC | **Implemented** | [CODE] `qc/fastqc_runner.py:57-125` |
| FQC dashboard | **Absent** | [CODE] no occurrence of `fqc` or `multiqc` in any `.py`, `Dockerfile`, `environment.yml`, `requirements.txt`. QC is shown in PARAFISH's own Dash page from parsed FastQC zips. |
| Trimmomatic | **Implemented** | [CODE] `preprocess/trimmomatic_runner.py:53-158` |
| PEAR | **Implemented** (VSEARCH/OTU branch, paired-end only) | [CODE] `preprocess/pear_runner.py:23-144` |
| cutadapt (primers) | **Implemented** | [CODE] `preprocess/cutadapt_runner.py:113-174` |
| VSEARCH | **Implemented** (filter, dereplication, clustering, chimera removal, OTU table, SINTAX) | [CODE] `analysis/vsearch_runner.py`, `databases/taxonomy_runner.py:945-993` |
| DADA2 | **Implemented** (PE and SE; R script embedded) | [CODE] `analysis/dada2_runner.py:116-397` |
| UNOISE | **Absent** | [CODE] no `unoise`/`cluster_unoise` anywhere |
| Kraken2 | **Absent** | [CODE] no `kraken` anywhere in code, Dockerfile or environment |
| SINTAX | **Implemented** (via VSEARCH `--sintax`) | [CODE] `databases/taxonomy_runner.py:954-976` |
| BLAST / megablast | **Implemented** — remote `nt` (megablast), optional local db (blastn, word_size 11), and custom-database BLAST | [CODE] `databases/ncbi_db.py:130-267`; `databases/custom_db.py:330-430` |
| BOLD / BOLDigger | **Implemented** (BOLDigger3 3.x CLI) | [CODE] `databases/bold_db.py`. Positive path proven on a subsampled COI run: `runs/_coi/current` — BOLD `ok`, 20/55 features. |
| MAFFT | **Implemented** | [CODE] `phylo/mafft_runner.py:16-102` |
| FastTree 2 | **Implemented** | [CODE] `phylo/fasttree_runner.py:15-80` |
| RAxML-NG | **Implemented** | [CODE] `phylo/raxmlng_runner.py:14-102` |
| Demultiplexing | **Absent** | [CODE] no `demultiplex`/`demux` anywhere; input is per-sample FASTQ |
| Phylogenetic placement (EPA-ng, pplacer) | **Absent** | [CODE] no occurrence; trees are built de novo |

### 2.2 Order of operations inside Preprocess and Analysis

[CODE] `ui/callbacks_progress.py:252-328, 399-512`.

1. cutadapt on the **raw** reads (primer removal).
2. Trimmomatic on cutadapt's output (or the raw reads if cutadapt passed through).
3. PEAR on Trimmomatic's paired output — **only** if the method includes VSEARCH, the sample is PE,
   and the sample id/label does not contain `extendedfrag`.
4. **OTU branch:** VSEARCH `--fastq_filter` per sample (on PEAR's `assembled.fastq` for PE, on the
   trimmed R1 for SE) → pool → dereplicate → cluster → de novo chimera removal **on the cluster
   centroids** → map all pooled reads back to the OTUs.
5. **ASV branch:** DADA2 receives Trimmomatic's `*_R1_trimmed` / `*_R2_trimmed` files (i.e. reads
   already primer-trimmed and quality-trimmed) and does its own filtering, error learning,
   denoising, merging and chimera removal. PEAR output is not used.

[INFER] A jury may ask about (5): DADA2 is fed reads that were already sliding-window trimmed by
Trimmomatic, which the DADA2 authors generally advise against because it alters the quality
profile the error model learns from. The code does this; no justification is recorded.

### 2.3 QC — FastQC

| Item | Value | Source |
|---|---|---|
| Command | `fastqc <R1> [<R2>] --outdir <run>/qc/<sample> --threads <n> --quiet` | [CODE] `qc/fastqc_runner.py:66-68` |
| Threads | 1 for a PE sample; `min(threads, 4)` for SE | [CODE] `qc/fastqc_runner.py:64` |
| Inputs | raw FASTQ | [CODE] |
| Outputs | `<stem>_fastqc.html`, `<stem>_fastqc.zip`, `<sample>_fastqc.log`; parsed into `qc/metrics/*.csv` | [CODE] `qc/fastqc_runner.py:21-28`; `qc/qc_metrics.py` |
| UI parameters | none ("FastQC has no manual parameters.") | [CODE] `ui/pages/page_setup.py:80` |
| PE vs SE | R2 added to the same invocation | [CODE] `qc/fastqc_runner.py:61-63` |

### 2.4 Primer trimming — cutadapt

**Presets** ([CODE] `analysis/vsearch_runner.py:60-76`):

| Preset key | Marker | Forward (5'→3') | Reverse (5'→3') | UI label |
|---|---|---|---|---|
| `18S_V4_TAReuk` | 18S | `CCAGCASCYGCGGTAATTCC` | `ACTTTCGTTCTTGATYRA` | "18S V4 — TAReuk454FWD1 / TAReukREV3" |
| `COI_Folmer` | COI | `GGTCAACAAATCATAAAGATATTGG` | `TAAACTTCAGGGTGACCAAAAAATCA` | "COI — Folmer (LCO1490 / HCO2198)" |
| `COI_Leray` | COI | `GGWACWGGWTGAACWGTWTAYCCYCC` | `TANACYTCNGGRTGNCCRAARAAYCA` | "COI — Leray (mlCOIintF / jgHCO2198)" |

UI labels: [CODE] `ui/pages/page_setup.py:207-213`. Two further options: `auto` (default) and `""`
("None — do not strip primers").

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `2a88619` (F2) **adds a fourth preset**; the three
above are unchanged, and `auto` still has no COI default.

| Preset key | Marker | Forward (5'→3') | Reverse (5'→3') | UI label |
|---|---|---|---|---|
| `COI_Leray_LoboR1` | COI | `GGWACWGGWTGAACWGTWTAYCCYCC` (mlCOIintF) | `TAAACYTCWGGRTGWCCRAARAAYCA` (LoboR1) | "COI — mlCOIintF / LoboR1 (lab)" |

[CODE] `analysis/vsearch_runner.py` `PRIMER_PRESETS`, `PRESET_MARKER`; `ui/pages/page_setup.py`
dropdown `setup-primer-preset`; `tests/test_cutadapt_runner.py::LabCoiPresetTests` (6 tests). The
sequences are those of `dados stora teresa/Primers utilizados metabarcoding.txt`, line "COI". No
stored run has used this preset yet.

| Item | Value | Source |
|---|---|---|
| Resolution order | explicit `primer_fwd`/`primer_rev` → preset → if preset is `auto`, the marker's standard pair | [CODE] `preprocess/cutadapt_runner.py:57-80` |
| `auto` mapping | only `{"18S": "18S_V4_TAReuk"}`. COI has **no** auto default (Folmer and Leray amplify different regions). | [CODE] `preprocess/cutadapt_runner.py:41` |
| PE command | `cutadapt -j <threads> --json cutadapt.json [--discard-untrimmed] -g <fwd> -G <rev> [--pair-filter=any] -o <id>_R1_primer.fastq.gz -p <id>_R2_primer.fastq.gz <R1> <R2>` | [CODE] `preprocess/cutadapt_runner.py:136-155` |
| SE command | `cutadapt -j <threads> --json cutadapt.json [--discard-untrimmed] -g <fwd> -o <id>_primer.fastq.gz <R1>` (reverse primer not used) | [CODE] `preprocess/cutadapt_runner.py:157-162` |
| `--discard-untrimmed` | on by default (`PARAFISH_DISCARD_UNTRIMMED`, default `True`); `--pair-filter=any` added in PE only | [CODE] `core/config.py` `DISCARD_UNTRIMMED`; `cutadapt_runner.py:137-154` |
| Error rate, overlap | cutadapt defaults (no `-e`, no `-O` passed) | [CODE] `cutadapt_runner.py:136-161` |
| Anchoring | not anchored (`-g SEQ`, not `-g ^SEQ`) | [CODE] same |
| Threads | `vsearch_params["threads"]` or 4 | [CODE] `cutadapt_runner.py:134` |
| Failure behaviour | never fails the stage: no primers, cutadapt missing, or a cutadapt error → reads pass through untouched with `primers_removed = 0.0` | [CODE] `cutadapt_runner.py:102-128,167-170` |
| Fallback | if cutadapt did not run, VSEARCH strips a **fixed number of bases** equal to the primer lengths (`--fastq_stripleft/--fastq_stripright`); if cutadapt ran, the flag `primers_removed_upstream` sets both to 0 | [CODE] `analysis/vsearch_runner.py:79-104,140-144`; `ui/callbacks_progress.py:274-277` |
| UI parameters | preset dropdown; free-text forward and reverse primers (override the preset) | [CODE] `ui/pages/page_setup.py:205-228` |
| Fixed | error rate, overlap, `--discard-untrimmed` (env var only) | [CODE] |

**CONTRADICTION (preset name vs laboratory primers).** The file
`dados stora teresa/Primers utilizados metabarcoding.txt` lists the COI pair actually used in the
laboratory as **mlCOIintF `GGWACWGGWTGAACWGTWTAYCCYCC` / LoboR1 `TAAACYTCWGGRTGWCCRAARAAYCA`**
[CODE]. The platform's `COI_Leray` preset has the same forward primer but reverse **jgHCO2198
`TANACYTCNGGRTGNCCRAARAAYCA`** [CODE]. The two reverse primers differ at three positions (3, 9 and 15 from the 5' end: A, W, W in LoboR1 versus
N, N, N in jgHCO2198). [INFER] jgHCO2198 is the more degenerate of
the two and its IUPAC codes include LoboR1's bases, so cutadapt would still match LoboR1-primed
reads with the preset; but the thesis must not say the data were amplified with jgHCO2198 unless
the laboratory confirms it. The same file lists a **16S Metazoa** pair
(`AGTTACYYTAGGGATAACAGCG` / `CCGGTCTGAACTCAGATCAYGT`) for which PARAFISH has no preset.

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** The laboratory's pair now has its own preset (`COI_Leray_LoboR1`,
above), so a new run can name the primers actually ordered. The `COI_Leray` preset and its label
still say jgHCO2198, by design (existing presets were not altered). **The stored COI run
`20260724T082931Z` was trimmed with the `COI_Leray` sequences** (jgHCO2198 as reverse) — that
remains what must be reported for that run.

### 2.5 Quality trimming — Trimmomatic

| Item | Value | Source |
|---|---|---|
| PE command | `java -jar <jar> PE -phred33 <R1> <R2> <R1_trimmed> <R1_unpaired> <R2_trimmed> <R2_unpaired> ILLUMINACLIP:<adapters>:2:30:10 LEADING:<l> TRAILING:<t> SLIDINGWINDOW:<w>:<q> MINLEN:<m>` | [CODE] `preprocess/trimmomatic_runner.py:83-106` |
| SE command | `java -jar <jar> SE -phred33 <R1> <trimmed> ILLUMINACLIP:… LEADING:… TRAILING:… SLIDINGWINDOW:… MINLEN:…` | [CODE] `trimmomatic_runner.py:110-123` |
| Adapters | `TruSeq3-PE.fa` (PE) / `TruSeq3-SE.fa` (SE); `ILLUMINACLIP` settings `2:30:10` are hard-coded | [CODE] `core/config.py:325-336`; `trimmomatic_runner.py:101` |
| Outputs | `<id>_R1_trimmed.fastq.gz`, `<id>_R2_trimmed.fastq.gz`, `<id>_R{1,2}_unpaired.fastq.gz`, `trimlog.txt` | [CODE] |
| UI parameters | LEADING, TRAILING, SLIDINGWINDOW (text `w:q`), MINLEN | [CODE] `ui/pages/page_setup.py:107-114` |

**Defaults — three sources that do not all agree:**

| Parameter | Setup page (what a UI run uses) | `RunConfig` default (what a headless run uses) | Runner's own fallback |
|---|---|---|---|
| LEADING | 3 | 3 | 3 |
| TRAILING | 3 | 3 | 3 |
| SLIDINGWINDOW | **4:15** | **4:20** | 4:20 |
| MINLEN | 50 | 50 | **200** |

[CODE] `ui/pages/page_setup.py:108-114`; `core/models.py:116-134`; `preprocess/trimmomatic_runner.py:102-105`.
**CONTRADICTION:** the sliding-window quality threshold is 15 from the UI and 20 from
`scripts/run_pipeline.py`. `tests/test_setup_defaults_single_source.py` does not cover it, because
it only inspects numeric `dcc.Input` fields and `setup-trim-slidingwindow` is a text field
([CODE] `tests/test_setup_defaults_single_source.py:27-39,47-53`). The golden datasets were produced
headless, i.e. with 4:20. [INFER]

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `8be350f` (F3): the Setup field now **reads the `RunConfig`
default**, so the page shows and uses **4:20**; `RunConfig` itself was not changed. The tooltip, which
said "Default: 4:15", now states the same value. [CODE] `ui/pages/page_setup.py` (`_RUN_DEFAULTS`,
`_DEFAULT_SLIDINGWINDOW`); `tests/test_setup_defaults_single_source.py::NonNumericSetupDefaultsTests`,
which reads the rendered layout and therefore covers the text field. **UI runs made before this
commit used 4:15**; the runner's own fallback of `MINLEN:200` (last column above) is unchanged and
still unreachable in practice, because `RunConfig` always supplies 50.

### 2.6 Merging — PEAR

| Item | Value | Source |
|---|---|---|
| Command | `pear -f <R1_trimmed> -r <R2_trimmed> -o <prefix> -v <min_overlap> -n <min_length> -m <max_length> -p <p_value> -u <max_uncalled_base> -j <threads> -g <test>` | [CODE] `preprocess/pear_runner.py:49-71` |
| Defaults | `-v 10`, `-n 200`, `-m 500`, `-p 0.01`, `-u 0.0`, `-j 4`, `-g 1` | [CODE] `core/models.py:107-115` |
| Outputs | `<id>.assembled.fastq`, `<id>.unassembled.{forward,reverse}.fastq`, `<id>.discarded.fastq`, `pear.log` | [CODE] `pear_runner.py:41-46` |
| UI parameters | Min overlap, Min length, P-value | [CODE] `ui/pages/page_setup.py:117-122` |
| Fixed | `max_length` 500, `max_uncalled_base` 0.0, `test` 1, threads 4 (the Setup "Threads" value is not passed to PEAR) | [CODE] `ui/callbacks_setup.py:512-516` updates only three keys |
| When it runs | VSEARCH branch, PE samples only; raises if `branch != VSEARCH` | [CODE] `pear_runner.py:24-25`; `ui/callbacks_progress.py:311` |
| Warning | merge rate < 50 % produces a warning string, not a failure | [CODE] `pear_runner.py:122-124` |
| `min_overlap` 10 pinned | [CODE] `tests/test_scientific_parameter_pins.py:91` |

[DOC] `DEVLOG.md:525-532`: raising `min_overlap` from 10 to ~50 is an open, undecided item.

### 2.7 Quality filtering — VSEARCH (OTU branch)

| Item | Value | Source |
|---|---|---|
| Command | `vsearch --fastq_filter <assembled.fastq> --sample <id> --fastq_maxee <ee> --fastq_minlen <len> --fastaout filtered.fasta --log step1_filter.log [--fastq_stripleft N] [--fastq_stripright N]` | [CODE] `analysis/vsearch_runner.py:118-144` |
| Defaults | `fastq_maxee` 1.0, `fastq_minlen` 200 | [CODE] `core/models.py:136-137`; UI the same (`page_setup.py:179-181`) |
| Post-processing | FASTA headers rewritten to `>{read_id};sample={sample_id}` so the sample tag survives VSEARCH's header truncation | [CODE] `vsearch_runner.py:34-53,158` |
| UI parameters | Max expected errors, Min length | [CODE] `page_setup.py:178-181` |
| No maximum length filter, no `--fastq_maxns` | [CODE] `vsearch_runner.py:118-132` |

### 2.8 Dereplication, OTU clustering, chimera removal, OTU table — VSEARCH

All samples are pooled and clustered **once** ([CODE] `analysis/vsearch_runner.py:286-437`).

| Step | Command | Source |
|---|---|---|
| Pool | concatenation of every sample's `filtered.fasta` into `pooled_filtered.fasta` | [CODE] `:310-314` |
| Dereplicate | `vsearch --derep_fulllength pooled_filtered.fasta --sizeout --output derep.fasta --log step2_pool_derep.log` plus `--minuniquesize N` **only if N > 1** | [CODE] `:326-328` |
| Cluster | `vsearch --cluster_size derep.fasta --id 0.97 --sizein --sizeout --centroids otus_raw.fasta --uc clusters.uc --log step3_cluster.log` | [CODE] `:337` |
| Chimeras | `vsearch --uchime_denovo otus_raw.fasta --sizein --sizeout --relabel OTU --relabel_keep --nonchimeras otus.fasta --chimeras chimeras.fasta --uchimeout chimera.uchime --log step4_chimera.log` | [CODE] `:344-349` |
| OTU table | `vsearch --usearch_global pooled_filtered.fasta --db otus.fasta --id 0.97 --uc map.uc --otutabout otu_table.tsv --threads <n> --log step5_otutable.log` | [CODE] `:372-388` |

| Fact | Source |
|---|---|
| Identity default **0.97**, allowed range [0.8, 1.0] (UI slider 0.80–1.00, step 0.01) | [CODE] `core/models.py:138`; `vsearch_runner.py:294-296`; `page_setup.py:167-175`; pinned in `tests/test_scientific_parameter_pins.py:71` |
| `min_unique_size` default **1** (`PARAFISH_MIN_UNIQUE_SIZE`); UI field 1–10 | [CODE] `core/config.py` `MIN_UNIQUE_SIZE`; `page_setup.py:196`; pinned `test_scientific_parameter_pins.py:108` |
| Chimera removal is **de novo, on OTU centroids, after clustering** — not on dereplicated reads, and with no reference-based step | [CODE] `vsearch_runner.py:337-349` |
| Clustering and `uchime_denovo` are not given `--threads` (only the mapping step is) | [CODE] `vsearch_runner.py:337-349,384-385` |
| Abundance is taken from the **mapping** (`--otutabout`), which maps *all* pooled reads, including those excluded from seeding; this is the run's canonical read count | [CODE] `vsearch_runner.py:212-235` docstring and code |
| `otu_id_map.csv` carries `size` (mapping) and `cluster_size` (clustering) separately | [CODE] `vsearch_runner.py:238-283` |
| OTU ids are `OTU1…OTUn`, assigned by `--relabel` after chimera removal | [CODE] `vsearch_runner.py:345-346` |
| Strand: no `--strand both` on clustering or mapping (VSEARCH default: plus only) | [CODE] `vsearch_runner.py:337,372-388` |

### 2.9 ASV denoising — DADA2

The R script is embedded in `analysis/dada2_runner.py:116-397` and run as
`Rscript --vanilla dada2_script.R <input_dir> <output_dir> <ref_fasta> <PE|SE> <truncLen_F> <truncLen_R> <maxEE_F> <maxEE_R> <truncQ> false <err_nbases>` ([CODE] `:483-498`).

| Step | Call | Source |
|---|---|---|
| Filter | `filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(F,R), maxEE=c(F,R), truncQ=truncQ, compress=TRUE)` | [CODE] `:269-278` |
| Error model | `learnErrors(filt, nbases=err_nbases, multithread=FALSE, randomize=TRUE)`, separately for F and R | [CODE] `:293-294` |
| Denoise | per sample: `derepFastq` → `dada(derep, err=…, multithread=FALSE)` | [CODE] `:307-323` |
| Merge | `mergePairs(dada_f, derep_f, dada_r, derep_r)` with DADA2 defaults (minOverlap 12, maxMismatch 0 — not overridden) | [CODE] `:314` |
| Chimeras | `removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE)` | [CODE] `:328` |
| Taxonomy | `assignTaxonomy()` only if a `ref_fasta` is configured **and** fits in RAM (6× file size rule); otherwise skipped | [CODE] `:37,78-100,235-243,429-437` |

| Fact | Source |
|---|---|
| Defaults: `truncLen_F` 0, `truncLen_R` 0, `maxEE_F` 2, `maxEE_R` 2, `truncQ` 2 | [CODE] `core/models.py:147-155` |
| UI parameters: TruncLen F/R, MaxEE F/R. `truncQ` is fixed at 2. | [CODE] `ui/pages/page_setup.py:238-245` |
| Pooling: samples are denoised independently (no `pool=`) | [CODE] `:307-323` |
| `err_nbases` adapts to free RAM: between 1e7 and 1e8 bases, at most 25 % of `MemAvailable` at 12 bytes/base | [CODE] `:41-63` |
| `truncLen` forced to 0,0 when there is no real pair or any sample is ITS | [CODE] `:449-460` |
| Multithreading is passed as `"false"` regardless of configuration | [CODE] `:496` |
| DADA2 refuses to start if less than ~850 MB is free | [CODE] `core/preflight.py:35,102-110`; `dada2_runner.py:409-413` |
| SE mode: single `learnErrors`/`dada`, no merging | [CODE] `:337-396` |
| Outputs: `asv_sequences.fasta` (headers `ASVn;size=N`), `asv_abundance.csv`, `asv_table.tsv`, `asv_table.rds`, `asv_id_map.csv`, `read_tracking.csv`, `taxonomy.csv` | [CODE] `:165-228` |
| Input files are selected by the regex `_R1_trimmed…fastq(.gz)` / `_R2_trimmed…` | [CODE] `:246-247` |
| **No `set.seed()`** in the R script, while `learnErrors(randomize=TRUE)` is used | [CODE] the string `set.seed` does not occur in `analysis/`. |

**CONTRADICTION (determinism of DADA2).** [DOC] `DEVLOG.md:1378` ("DADA2 is deterministic, so it was
safe to freeze") versus [CODE] no seed is set and `randomize=TRUE` is passed. [INFER] at fixture
scale all reads fit under `nbases`, so the random sample order may not change the result; this has
not been shown for larger inputs.

### 2.10 OTU and ASV in parallel

The analysis method is `vsearch`, `dada2` or `both`; the Setup page defaults to **`both`**
([CODE] `ui/pages/page_setup.py:158`). With `both`, the two branches run one after the other in the
same thread (not concurrently), and each later stage runs once per feature set into its own
directory (`databases/otu`, `databases/asv`, `phylo/otu`, `phylo/asv`)
([CODE] `ui/callbacks_progress.py:399-512,731-755,851-873`). If one branch fails the other is kept
([CODE] `:505-512,546-566`). A comparison joined on taxon, never on feature id, is written to
`analysis/comparison/` ([CODE] `:789-801`; `analysis/feature_set_comparison.py`).

**No stored run contains both feature sets.** [CODE] run inventory, §7.

### 2.11 Taxonomic assignment commands

**SINTAX** ([CODE] `databases/taxonomy_runner.py:954-976`):

```
vsearch --sintax <features.fasta> --db <reference.sintax.fasta> --sintax_cutoff <c>
        --tabbedout <prefix>_result.csv --threads <n> --randseed 42 --log <prefix>.log
```

Default cutoff 0.8 ([CODE] `core/models.py:139`; pinned `tests/test_scientific_parameter_pins.py:83`).
No `--strand both` ([CODE]).

**Remote NCBI BLAST** ([CODE] `databases/ncbi_db.py:142-173`):

```
blastn -query <features.fasta> -out blast_result.tsv
       -outfmt "6 qseqid sseqid pident qcovhsp length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames"
       -task megablast -evalue 1e-10 -max_target_seqs <n> -perc_identity 90.0 -qcov_hsp_perc 80.0
       -db nt -remote
```

| Fact | Source |
|---|---|
| Remote: `task=megablast`, no `-word_size`. Local db (`PARAFISH_NCBI_LOCAL_DB=1`): `task=blastn`, `word_size=11`, `-num_threads 8`. | [CODE] `ui/callbacks_progress.py:696-705` |
| `perc_identity` 90.0 and `qcov_hsp_perc` 80.0 are fixed (no UI control) | [CODE] `databases/ncbi_db.py:63-64` |
| `max_target_seqs`: UI default 5; the full-pipeline fallback is **20**, the standalone Databases callback fallback is 5 | [CODE] `page_setup.py:342`; `callbacks_progress.py:704`; `callbacks_databases.py:1321` |
| Refusals are retried with waits of 30 s and 120 s, then reported as a failure | [CODE] `ncbi_db.py:44,189-229` |
| blastn exits 0 when NCBI refuses; refusal is detected from stderr signatures (`cpu usage limit`, `sigxcpu`, `search failed`, …) | [CODE] `ncbi_db.py:29-37,215-223` |
| Timeout: `NCBI_TIMEOUT_SEC` (code default 1800; `.env` on the dev machine 3600) | [CODE] `core/config.py`; `.env:29` |
| Only the **best hit by bitscore** per query enters the consensus; all hits go to `taxonomy_hits.csv` | [CODE] `taxonomy_runner.py:1155-1167` |

**BOLDigger3** ([CODE] `databases/bold_db.py:304-333`):

```
boldigger3 identify <sanitized.fasta> [<local.ddb>] --db <1-8> --mode <1-3> --thresholds 97 95 90 85
```

`--db` UI options 1, 2, 5, 6 (default 1, "Animal Library (Public)"); `--mode` UI default **1**
("Rapid Species Search") while `RunConfig`/`BOLDConfig` default **2** ([CODE] `page_setup.py:344-365`;
`core/models.py:170-174`; `bold_db.py:186-188`). Up to 3 attempts with a 30 s wait on HTTP 503
([CODE] `bold_db.py:376-377,563-590`). Thresholds 97/95/90/85 map to species/genus/family/order
([CODE] `bold_db.py:433,454-493`).

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `8be350f` (F3): the Setup **BOLD mode dropdown now reads
`RunConfig`**, so its default is **2** ("Genus and Species Search"), the same as a headless run. UI runs
made before this commit used mode 1. **NCBI max hits was deliberately not changed:** `RunConfig`
defines no default for it, so there is no single source to align to. The three numbers remain —
Setup field 5, full-pipeline fallback 20, `NCBIConfig.max_target_seqs` 5 — and a test pins that
`RunConfig` has no `ncbi_maxhits` key ([CODE]
`tests/test_setup_defaults_single_source.py::test_ncbi_max_hits_has_no_run_config_default_to_align_to`).
In a UI run the field's value (5) is what reaches BLAST; the fallback of 20 applies only to a
headless run that sets nothing.

**Custom database** ([CODE] `databases/custom_db.py:303-409`): `makeblastdb -in <ref> -dbtype nucl -out custom_ref_db`, then
`blastn -query … -db … -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp" -evalue 1e-5 -perc_identity 90.0 -max_target_seqs 5 -num_threads <n>`.
No `-task` is passed ([INFER] blastn's default task, megablast, applies). Upload validation rejects
empty bodies, non-nucleotide characters, all-gap records, < 50 % A/C/G/T, and sequences shorter
than 50 bp ([CODE] `custom_db.py:106-150`).

### 2.12 Alignment and phylogeny commands

**MAFFT** ([CODE] `phylo/mafft_runner.py:39-54`):

- primary: `mafft --auto [--adjustdirection] --thread <n> sanitized.fasta` (or `--maxiterate 1000` instead of `--auto` when `mode == "accurate"`)
- fallback if the primary exits non-zero: `mafft --retree 1 --nofft [--adjustdirection] --thread <n> sanitized.fasta`

Input is the taxonomy-annotated FASTA when Databases ran, else the bare representatives
([CODE] `ui/callbacks_progress.py:855-856`). Tips are renamed `OTU…`/`ASV…`; the mapping is written
to `reversed_asvs.csv` ([CODE] `mafft_runner.py:19-21,86-91`). On the stored real run `--auto`
selected L-INS-i ([CODE] `runs/_real/current/phylo/otu/mafft/mafft.log`).

**FastTree** ([CODE] `phylo/fasttree_runner.py:41-50`): `FastTree -nt -gtr -gamma <aligned.fasta>`.
Refuses fewer than 4 sequences ([CODE] `phylo/__init__.py:14`).

**RAxML-NG** ([CODE] `phylo/raxmlng_runner.py:35-54`):

```
raxml-ng --parse --msa aligned.fasta --model GTR+G --prefix raxml_PARSE
raxml-ng --all   --msa aligned.fasta --model GTR+G --prefix raxml --threads auto{8}
         --seed 42 --bs-metric fbp --bs-trees <N> --redo [--tree <fasttree.nwk>]
```

**Modes** ([CODE] `core/models.py:27-32`; `phylo/phylo_runner.py:79-106`): `fast` = MAFFT + FastTree;
`robust` = MAFFT + RAxML-NG; `hybrid` (default) = MAFFT + FastTree, then RAxML-NG started from the
FastTree tree (`--tree`).

**Dead or inconsistent phylo controls:**

| Finding | Source |
|---|---|
| Setup "MAFFT: Auto / Accurate" (`setup-mafft-mode`) and "Adjust direction" (`setup-adjust-direction`) are **not read by any callback**. MAFFT always runs with `--auto --adjustdirection`. | [CODE] the two ids occur only in `ui/pages/page_setup.py:425,437` |
| FastTree reads `params.get("gtr")` / `params.get("gamma")`, but `RunConfig.phylo_params` defines `fasttree_gtr` / `fasttree_gamma`. The keys never match, so both flags are always on. | [CODE] `phylo/fasttree_runner.py:42-45`; `core/models.py:167-168` |
| RAxML bootstrap replicates: UI default **100**, `RunConfig` default **10**. The golden real run was built with 10. | [CODE] `page_setup.py:435`; `core/models.py:164`; `runs/_real/current/phylo/otu/raxml/raxml.raxml.log` ("bootstrap replicates: parsimony (10)") |
| The runner returns `support_tree = raxml.raxml.supportFBP`, but RAxML-NG 1.2.2 writes `raxml.raxml.support`. The returned path does not exist; readers find the file through a separate on-disk lookup. | [CODE] `phylo/raxmlng_runner.py:76-78`; `runs/_real/current/phylo/otu/raxml/`; `ui/feature_sets.py:88-92`; `phylo/phylo_metrics.py:274-281` |
| Model (`GTR+G`), bootstrap metric (`fbp`), seed (42), outgroup: no UI control | [CODE] `page_setup.py:394-449` |

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `8be350f` (F3): the Setup "RAxML Bootstraps" field now **reads
`RunConfig`**, so its default is **10**, the same as a headless run (the field still accepts 10–1000).
UI runs made before this commit defaulted to 100. The other findings in the table above (dead MAFFT
controls, FastTree key mismatch, support-tree file name) are **unchanged**.

### 2.13 Paired-end versus single-end — summary

| Stage | Paired-end | Single-end | Source |
|---|---|---|---|
| Detection | file stem ending in `_R1`/`_1`/`_forward`/`_fwd` (mate 1) or `_R2`/`_2`/`_reverse`/`_rev` (mate 2), separators `_ . -` | anything else; also any name containing `extendedfrag`, `assembled` or `merged` | [CODE] `analysis/sample_parsing.py:60-71` (`_detect_mode_and_mate`) |
| FastQC | R1 and R2 | R1 | [CODE] |
| cutadapt | `-g fwd -G rev`, `--pair-filter=any` | `-g fwd` only | [CODE] |
| Trimmomatic | `PE`, `TruSeq3-PE.fa` | `SE`, `TruSeq3-SE.fa` | [CODE] |
| PEAR | runs | skipped | [CODE] |
| VSEARCH filter input | PEAR `assembled.fastq` | trimmed R1 | [CODE] `callbacks_progress.py:404-411` |
| DADA2 | F/R error models, `mergePairs` | single error model, no merge | [CODE] |

[UNKNOWN] No stored run with `progress.json` evidence of a completed single-end analysis was
identified; SE is covered by unit tests only ([CODE] e.g. `tests/test_cutadapt_runner.py`).

---

## 3. Taxonomic assignment

### 3.1 Reference sources

| Source | Method | Version in the repository | How obtained | Where it lives | Source |
|---|---|---|---|---|---|
| SILVA | VSEARCH SINTAX | `SILVA_138.2_SSURef_NR99_tax_silva.fasta` (824 MB) | manual download (README gives the `wget`); not in git | `DB_SILVA_PATH`; converted copy `*.sintax.v4.fasta` cached beside it | [CODE] `.env:1`; `databases/silva_db.py:60-67` |
| PR2 | VSEARCH SINTAX | `pr2_version_5.1.1_SSU_mothur.fasta` + `.tax` (334 MB) | manual download | `DB_PR2_PATH`, `DB_PR2_TAX_PATH`; converted copy `*.sintax.v4.fasta` | [CODE] `.env:2`; `databases/pr2_db.py:77-93` |
| NCBI | blastn | remote `nt`, **version/date not recorded** | queried live | — | [CODE] `ncbi_db.py:171` |
| BOLD | BOLDigger3 | remote BOLD v5, **version/date not recorded** | queried live | optional local metadata `.ddb` under `DB_BOLD_PATH` | [CODE] `bold_db.py` |
| Custom | blastn or SINTAX (auto-detected from `;tax=` headers) | user upload, up to 500 MB | Setup or Databases page upload | run directory | [CODE] `page_setup.py:310-321`; `taxonomy_runner.py:1559-1617` |

The dataclass defaults `SILVAConfig.version = "138.2"` and `PR2Config.version = "5.0"` are labels
only; nothing reads them ([CODE] `silva_db.py:38`; `pr2_db.py:61`). **CONTRADICTION (minor):** the
PR2 label says 5.0, the file in use is 5.1.1.

**No checksum of any reference is computed or recorded.** [CODE] no hashing of database files in
`databases/` or `core/`; cache validity is decided by mtime and a format-version string
(`silva_db.py:64`, `reference_filter.py:314-323`).

**Marker-filtered references.** Before SINTAX, SILVA and PR2 are reduced to the run's amplicon
([CODE] `databases/reference_filter.py:336-447`): (1) keep headers containing `d:Eukaryota` for
markers 18S/28S/ITS, converting U→T; (2) `cutadapt -g "<fwd>;required...<revcomp(rev)>;required" -e 0.15 --match-read-wildcards --discard-untrimmed -m 150 -M 800`;
(3) if the primer step keeps < 20 % of records, fall back to lineage-only. The result is cached as
`<stem>.<marker>.<FWD8>_<REV8>.filtered.v1.fasta` with a JSON sidecar.

Measured sizes recorded in a stored run ([CODE] `golden/fixture/taxonomy/databases/otu/source_status.csv`):
SILVA 788 MB (510 495 sequences) → 58 940 Eukaryota → **25 MB, 54 207 sequences**; PR2 343 MB
(240 201) → 232 042 → **76 MB, 170 678 sequences**, with the TAReuk primers.

### 3.2 Marker × source matrix

[CODE] `databases/taxonomy_runner.py:741-779`.

| Source | 18S | 16S | COI | ITS / 28S | unknown / `OTHER` |
|---|---|---|---|---|---|
| SILVA | applies | applies | **not_applicable** | not_applicable | applies (never marked incompatible) |
| PR2 | applies | not_applicable | **not_applicable** | not_applicable | applies |
| BOLD | not_applicable | not_applicable | applies | not_applicable | applies |
| NCBI | applies | applies | applies | applies | applies (no restriction) |
| Custom | applies | applies | applies | applies | applies |

`SOURCE_SUPPORTED_MARKERS = {"BOLD": {"COI"}, "SILVA": {"16S","18S"}, "PR2": {"18S"}}`.
`_marker_incompatible` returns False for an empty or `OTHER` marker. `not_applicable` is decided
**before** the source runs (no filtering, no loading) and again in `_reconcile_source_status` for
sources that ran and assigned nothing ([CODE] `:1491-1497,1526-1530,1656-1663,865-870`).

`Marker` enum values are `COI`, `ITS`, `18S`, `28S`, `OTHER` — **there is no `16S` member**
([CODE] `core/models.py:35-42`), although `_infer_marker_from_name` can return `"16S"`
([CODE] `analysis/sample_parsing.py:28-29`) and `SOURCE_SUPPORTED_MARKERS` mentions it. [INFER] a
16S sample falls back to `OTHER`.

**CONTRADICTION (where the run's marker comes from).**

- [DOC] `DEVLOG.md:665-738` (Task 1, 2026-07-24): the marker is derived from the primer preset via
  `resolve_run_marker()`, so a COI/Leray run is recognised as COI and SILVA/PR2 are marked
  `not_applicable`.
- [CODE] That is true in the Setup page, in the standalone Databases callback
  (`ui/callbacks_databases.py:1328`) and for `RunConfig.marker` (`ui/runtime_state.py:278-282`).
- [CODE] But the **full-pipeline thread** — the path used by "Start Pipeline" and by
  `scripts/run_pipeline.py` — computes `run_marker = resolve_marker(sample.marker for sample in run_config.samples)`
  (`ui/callbacks_progress.py:669`), which looks only at the marker inferred from the **file names**
  and never reads `run_config.marker`.
- [CODE] Verified by execution on 2026-10-12: for two samples with marker `OTHER` and preset
  `COI_Leray`, `resolve_marker` returns `'OTHER'` and `resolve_run_marker` returns `'COI'`; with
  `'OTHER'`, `_marker_incompatible("SILVA", …)` is `False`.
- Consequence [INFER]: in a full-pipeline COI run whose file names do not contain `coi`, SILVA and
  PR2 are **not** marked `not_applicable`; they are attempted whole and either skipped by the memory
  gate or run and reported `no_hits`. The stored COI runs show marker `OTHER` in `progress.json`
  ([CODE] §7.2), consistent with this.

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `64c8446` (F1) **resolves this on the branch**. The full-pipeline
thread now calls `_pipeline_run_marker(run_config)`, which is `resolve_run_marker(sample markers,
vsearch_params)` — preset first, file names as fallback — the same resolution as Setup and the
standalone Databases callback. [CODE] `ui/callbacks_progress.py` (`_pipeline_run_marker` and its call
in the Databases stage); `tests/test_pipeline_run_marker.py` (4 tests), including one that drives the
real `_run_pipeline_thread` and captures the marker handed to the taxonomy stage (it failed with
`'OTHER' != 'COI'` before the change). Verified cases: files `OTHER` + preset `COI_Leray` → `COI`, and
SILVA/PR2 incompatible; 18S files + `auto` → `18S`; no preset and no marker → unchanged.
Two things are **not** changed: the `markers` key of the `run_config` block still lists the
file-name markers (so it can still read `OTHER` for such a run — the resolved marker is now recorded
separately as `parameters.databases.marker`, §5.1), and the stored COI runs were made before the fix.

### 3.3 Consensus rule

[CODE] `databases/taxonomy_runner.py:265-477, 1474-1483, 1717`.

**Input records.** One record per source per feature: SINTAX sources give one call each; NCBI gives
its best hit by bitscore; BOLD gives its top hit. Records are appended in the fixed order
SILVA → PR2 → custom SINTAX → NCBI → BOLD.

**Default method `lca`** (`_lca_consensus`) — the exact algorithm:

1. Walk the ranks from deepest to shallowest: species, genus, family, order, class, phylum, kingdom.
2. At each rank collect the values from the records that **have a value at that rank**; records
   silent at that rank are ignored.
3. The first rank at which those values are all identical (case-insensitive) is the consensus rank.
4. Every rank **above** it is filled with the **first non-empty value in record order** — not with a
   consensus of the ancestors.
5. If no rank qualifies, kingdom is set to the shared kingdom or to `unclassified`.

Two consequences that must be stated accurately in the thesis:

- **It is not a strict lowest common ancestor.** A rank asserted by a single source, with the other
  sources silent at that rank, is accepted. Since SILVA contributes only domain + species and PR2
  only domain + genus + species (§3.5), the middle ranks are in practice whatever NCBI or BOLD says.
- **The upper ranks can be internally inconsistent with the consensus rank.** Example in a stored
  run [CODE] `runs/_real/current/databases/otu/taxonomy_consensus.csv`, OTU9: `Consensus_Genus =
  Mycosphaerella`, `Consensus_Species = Nothopassalora personata`. [INFER] the genus is PR2's (first
  in record order), the species is the agreed/only species value.

**No thresholds and no tie-breaks are applied inside the LCA.** Identity, coverage and bootstrap do
not gate participation; a 90 % BLAST hit weighs the same as a 100 % one. [CODE] `:265-299`.

**Alternative methods** selectable in Setup ("Consensus method"): `strict` (walk from kingdom down,
keep ranks while all non-empty values agree, stop at first disagreement; falls back to LCA if
nothing agrees) and `majority` (deepest rank where one value is held by **more than half of all
records**; ties broken by `max` on `(count, name)`, i.e. alphabetically last; upper ranks filled as
in LCA). [CODE] `:302-348`; `ui/pages/page_setup.py:366-376`.

**Evidence columns** on each consensus row ([CODE] `:407-477`):

| Column | Rule |
|---|---|
| `Best_Source`, `Best_Identity`, `Best_Coverage` | the alignment record with the highest (identity, coverage); records with identity ≤ 0 or no taxon are excluded; SINTAX records never qualify |
| `Consensus_Sources` | every source that named a taxon at any rank for this feature, joined by `\|` |
| `Reliability` | graded on the best alignment record against the consensus rank; if there is no alignment, on the strongest SINTAX bootstrap |
| `Source_of_consensus` | the method name (`lca`, `strict`, `majority`) |

**Reliability grading** ([CODE] `:556-614`):

| Evidence | high | moderate | low |
|---|---|---|---|
| Alignment (NCBI, BOLD, custom BLAST) | identity ≥ rank threshold and coverage ≥ 80 % | identity ≥ threshold, coverage < 80 % | identity < threshold |
| SINTAX bootstrap | ≥ 0.95 | ≥ 0.85 | < 0.85 |

Rank identity thresholds: species 97, genus 95, family 90, order 85, class 80, phylum 75, kingdom 0.

**Dead control.** The Setup field "Species threshold (%)" (`setup-species-threshold`, default 97) is
stored as `vsearch_params["species_threshold"]` and **never read** ([CODE] the key is written at
`ui/callbacks_setup.py:580` and occurs nowhere else in the code). The 97 % used for grading is the hard-coded constant above.

**Dead control.** The Setup field "SINTAX cutoff" inside the VSEARCH advanced panel
(`setup-vsearch-sintax`) is passed into the callback as `v_sintax` and never used; the value that
takes effect is the slider `setup-sintax-cutoff` in the Databases panel
([CODE] `ui/callbacks_setup.py:306,350,561`). `tests/test_setup_defaults_single_source.py:36` maps
the dead control, not the live one.

### 3.4 NCBI taxonomy resolution

| Fact | Source |
|---|---|
| Method: `Bio.Entrez.efetch(db="taxonomy", id=<≤200 taxids>, retmode="xml")`, parsing `LineageEx` | [CODE] `databases/ncbi_lineage.py:93-163` |
| **No local taxdump is used.** `taxdump`/`taxonkit` occur nowhere in the code. | [CODE] grep |
| Requires `ENTREZ_EMAIL`; without it the resolver returns the cache plus a warning and ranks above genus stay blank | [CODE] `:113-120` |
| Optional `NCBI_API_KEY` / `ENTREZ_API_KEY`; `Entrez.tool = "PARAFISH"` | [CODE] `:128-132`; `core/config.py:210-211` |
| Cache: one JSON keyed by taxid, at `<DB_ROOT>/ncbi/lineage_cache.json`, never invalidated | [CODE] `:57-74`; path [DOC] `DEVLOG.md` 2026-07-15 entry |
| Rank mapping: `superkingdom` and `kingdom` both map to the `kingdom` column, `kingdom` winning when both exist | [CODE] `:23-36` |
| Merged taxids: `AkaTaxIds` registered so the queried id resolves | [CODE] `:151-159` |
| Only **blank** fields are filled; a name BLAST already reported is not overwritten | [CODE] `:166-178` |
| Genus is derived from the BLAST binomial (first word), with placeholders such as "uncultured", "environmental sample" rejected | [CODE] `taxonomy_runner.py:193-220` |
| blastn needs NCBI's `taxdb` (`$BLASTDB`) to print `sscinames` | [CODE] `ncbi_db.py:175-180`; `.env:26` |

[INFER] Because lineages are fetched live and cached without a date, the lineage attached to a
taxid reflects NCBI Taxonomy on the day it was first fetched; this is not recorded per run.

### 3.5 Linnaean backbone

Seven ranks: kingdom, phylum, class, order, family, genus, species
([CODE] `databases/taxonomy_runner.py:64`; internal key for class is `class_`, displayed as `class`).

**SILVA and PR2 are not excluded from the consensus.** They are excluded from supplying the
**middle ranks**:

| Source | Ranks it contributes | Implementation | Source |
|---|---|---|---|
| PR2 (9 levels: Domain; Supergroup; Division; Subdivision; Class; Order; Family; Genus; Species) | `d` (→ kingdom column), `g`, `s` | `PR2_RANK_LABELS = ("d", None, None, None, None, None, None, "g", "s")`; placeholders `*_X`, `*_XX`, `*_sp.` are blanked, not dropped | [CODE] `core/utils.py:255-289`; `databases/pr2_db.py:14-55` |
| SILVA (`tax_silva` lineage of variable depth, no declared ranks) | first field as `d`, last field (organism name) as `s` | `silva_rank_labels(depth)` | [CODE] `core/utils.py:300-330`; `databases/silva_db.py:13-33` |
| NCBI | all seven, via Entrez | §3.4 | [CODE] |
| BOLD | phylum … species (BOLD's own columns), limited by the identity thresholds | [CODE] `bold_db.py:433-493` |

Why ([DOC] `core/utils.py:255-277,300-321` comments; `DEVLOG.md` 2026-07-16 "Harmonised taxonomic
backbone"): mapping PR2's Division/Class and SILVA's positional levels onto Linnaean ranks produced
`Phylum = Opisthokonta`, `Class = Craniata`, `Phylum = Amorphea`, `Class = Obazoa`,
`Genus = Teleostei` for a fish.

Stated cost ([DOC] same comments): PR2's genuine protist classes (e.g. Dinophyceae) are lost; SILVA's
ranks could be recovered from its `tax_slv_*.txt` rank file, which is not used.

Note that the SINTAX prefix `d:` (domain) is stored in the **kingdom** column
([CODE] `taxonomy_runner.py:1014-1023`), so `Consensus_Kingdom` typically reads `Eukaryota`.

`scripts/verify.py` check (f) uses a separate `LINNEAN_SOURCES` set that excludes SILVA and PR2 and
fails a run in which no Linnaean source fed any consensus ([CODE] `scripts/verify.py`, check (f);
[DOC] `DEVLOG.md:2335-2346`).

### 3.6 Alignment identity versus SINTAX confidence

They are kept in **separate columns** in every deliverable.

| File | Alignment evidence | SINTAX evidence | Source |
|---|---|---|---|
| `taxonomy_summary.csv` | `align_source`, `align_taxon`, `align_identity` | `sintax_source`, `sintax_taxon`, `sintax_confidence`, `sintax_status` | [CODE] `taxonomy_runner.py:1295-1314` |
| `source_assignments.csv` | `identity`, `coverage` | `confidence` | [CODE] `:579` |
| `taxonomy_hits.csv` | `identity`, `coverage` (2 decimals) | `confidence` | [CODE] `HITS_COLUMNS` `:1211` |
| `taxonomy_consensus.csv` | `Best_Identity`, `Best_Coverage`, `Best_Source`, `BOLD_similarity` | `SINTAX_species`, `SINTAX_bootstrap` | [CODE] `:494-515` |

`confidence` no longer falls back to the BLAST bitscore ([CODE] `:238-245`).
**CONTRADICTION (stale docstring):** the docstring of `best_sintax_record` (`:386-389`) still says
the confidence field "falls back to an aligner's BLAST bitscore"; the code at `:245` does not.

`sintax_status` values: `ok`, `no SINTAX assignment`, `SINTAX not run — see source_status.csv`
([CODE] `:1368-1375`).

`source_status.csv` has columns `source, status, detail`; status ∈ `ok`, `no_hits`, `failed`,
`skipped`, `not_applicable` ([CODE] `:806-882`).

### 3.7 Seeds and determinism

| Step | Seed / determinism | Source |
|---|---|---|
| VSEARCH SINTAX | `--randseed 42` (`SINTAX_RANDSEED`) | [CODE] `taxonomy_runner.py:803,972-973` |
| RAxML-NG | `--seed 42` | [CODE] `core/models.py:166`; `raxmlng_runner.py:30,48` |
| Rarefaction, alpha diversity | `np.random.default_rng(seed=0)` (default argument of `build_analysis_metrics`) | [CODE] `analysis/analysis_metrics.py:443,474` |
| Fixture subsampling | fixed seed (default 1) | [CODE] `scripts/make_fixture.py` docstring |
| `conflict_report.csv` | rows sorted by feature id | [CODE] `taxonomy_runner.py:649` |
| Consensus rows | `sorted(records_by_asv.items())` | [CODE] `:1717` |
| DADA2 | **no seed** (`randomize=TRUE`) | [CODE] §2.9 |
| PERMANOVA | **no seed** passed to `skbio.stats.distance.permanova` in either implementation | [CODE] `analysis/ordination.py:268`; `databases/database_metrics.py:263` |
| VSEARCH clustering | no randomisation; `--cluster_size` orders by abundance, so ties between equal-abundance clusters can permute OTU ids | [DOC] `DEVLOG.md:2483-2489` |
| NCBI BLAST, BOLD, Entrez | remote and time-dependent; not reproducible by construction | [INFER] |
| BIOM | contains a wall-clock `date`; excluded from golden comparison | [CODE] `scripts/golden.py:172` |

[DOC] `DEVLOG.md:2249-2271`: before the seed was added, two runs of identical fixture data gave
SILVA bootstrap 0.99 and 0.98 for the same OTU.

### 3.8 Status of the multi-source SILVA + PR2 + NCBI run

**CONTRADICTION between documentation and stored data.**

- [DOC] `DEVLOG.md:374-388` (2026-07-24): "the multi-source consensus SILVA+PR2+NCBI has NEVER run
  with enough RAM … has never executed on real data at full depth".
- [CODE] `runs/_real/current/databases/otu/source_status.csv` (identical to
  `golden/real/taxonomy/…`): `SILVA ok 9/9`, `PR2 ok 9/9`, `NCBI ok 9/9, 1046 hits`; and
  `taxonomy_consensus.csv` has `Consensus_Sources = NCBI|PR2|SILVA` on all nine OTUs.
  `scripts/verify.py runs/_real/current` today: 7 passed, 0 failed, check (e) "NCBI=9 · PR2=9 ·
  SILVA=9". [DOC] `DEVLOG.md:2348-2378` records this run as "Phase 3 CLOSED and frozen" on
  2026-07-16.

What is therefore true, stated precisely:

| Configuration | Ever completed? | Evidence |
|---|---|---|
| SILVA + PR2 (marker-filtered) + NCBI remote megablast, 18S V4, **2 samples, 9 OTUs**, full read depth | **Yes** | [CODE] `runs/_real/current` |
| Same, on the 5 000-read fixture (4 OTUs) | Yes | [CODE] `runs/_fixture/current` |
| SILVA + PR2 **unfiltered** (full references) | **No** | [CODE] every stored run with the full reference shows `skipped` (e.g. `runs/projeto_1/20260723T072541Z/…/source_status.csv`: needs ~2.4 GB, 0.9 GB available) |
| All four sources on the COI sample | SILVA/PR2 skipped for RAM; NCBI + BOLD completed (one run) | [CODE] same file |
| Any multi-source run with more than 2 samples | **No** | [CODE] run inventory §7 |
| Multi-source run on the ASV feature set | Only SILVA + PR2 (no NCBI) on the fixture | [CODE] `runs/_fixture_asv/current` |

Where it fails: VSEARCH SINTAX needs about 3× the reference FASTA in RAM
([CODE] `core/preflight.py:24,73`), i.e. ~2.5 GB for SILVA and ~1.0 GB for PR2, on a machine with
3.0 GB total. The memory gate skips the source before it starts
([CODE] `taxonomy_runner.py:893-912,1512-1516`). Two later COI runs stopped inside the Databases
stage without finishing ([CODE] §7.2); [DOC] `DEVLOG.md:621-622` attributes one to the OOM killer.

---

## 4. Downstream analysis

### 4.1 Phylogeny

| Item | Value | Source |
|---|---|---|
| Alignment | MAFFT `--auto --adjustdirection` | §2.12 |
| Tree, fast | FastTree `-nt -gtr -gamma`; SH-like local supports (FastTree default, 1 000 resamples) | [CODE] `fasttree_runner.py`; `runs/_real/current/phylo/otu/fasttree/fasttree.log` |
| Tree, robust/hybrid | RAxML-NG `--all` (ML search + bootstrap), model `GTR+G`, Felsenstein bootstrap (`fbp`), seed 42 | [CODE] `raxmlng_runner.py` |
| Bootstrap replicates | UI default 100 (min 10, max 1000, step 10); headless default 10 | [CODE] `page_setup.py:435`; `core/models.py:164` |
| Minimum tips | 4 (`MIN_TREE_SEQUENCES`), enforced before the stage and inside FastTree | [CODE] `phylo/__init__.py:14`; `callbacks_progress.py:862` |
| Model selection | none (no ModelTest-NG); model is fixed | [CODE] |
| Rooting / outgroup | none from the UI; `outgroup` parameter exists in the runner but is never set | [CODE] `raxmlng_runner.py:32,53` |
| Support threshold in the viewer | 70 | [CODE] `phylo/phylo_metrics.py:25` |
| Outputs | `aligned.fasta`, `tree.nwk`, `raxml.raxml.{bestTree,support,bootstraps,bestModel}` | [CODE] |
| Derived metrics | alignment summary, per-column composition, branch support table, weighted and unweighted UniFrac matrices, Faith's PD | [CODE] `phylo/phylo_metrics.py:226-327`; `phylo/phylo_diversity.py:67` |
| Visualisation | tree rendered in the Phylo page for up to 250 tips | [CODE] `phylo/phylo_diversity.py:125` |
| Trees are byte-reproducible between runs | [DOC] `DEVLOG.md:1076` (measured on the OTU sets); enforced for the fixture by `golden/fixture/phylo` |

### 4.2 Statistics and visualisation

**Alpha diversity** ([CODE] `analysis/analysis_metrics.py:256-283,474-493`): computed per sample on a
matrix **rarefied to the depth of the shallowest non-empty sample** (hypergeometric subsampling
without replacement, seed 0). Columns: `reads, observed, chao1, shannon, simpson, pielou, hill_q0,
hill_q1, hill_q2` (scikit-bio `chao1`, `shannon`, `simpson`, `pielou_e`). The depth used is written
to `rarefaction_depth.csv`. Rarefaction curves: 20 depths × 10 iterations ([CODE] `:31-32,304-316`).
Faith's PD from the tree ([CODE] `phylo/phylo_diversity.py:67`).

**Beta diversity / ordination** ([CODE] `analysis/ordination.py`):

| Item | Value |
|---|---|
| Default | **Aitchison PCA**: prevalence filter → pseudocount → CLR → PCA by SVD, centred, **not** scaled |
| Prevalence filter | feature must be present in ≥ 2 samples (`MIN_PREVALENCE = 2`; UI adjustable) |
| Zero replacement | fixed pseudocount **0.5** added to counts (`ZERO_REPLACEMENT`) |
| CLR | `log(x + 0.5) − mean(log(x + 0.5))` per sample |
| Alternatives | PCoA (scikit-bio) on Bray-Curtis, Jaccard, weighted UniFrac, unweighted UniFrac (the last two need a tree) |
| Guard, n < 3 | `MIN_SAMPLES = 3`: raises `OrdinationError`; the UI shows an "Ordination unavailable" state, **no plot** |
| Guard, 3 ≤ n < 5 | plot drawn with a low-confidence warning (`LOW_CONFIDENCE_SAMPLES = 5`) |
| Guard, features | fewer than 3 features after filtering → error (`MIN_FEATURES = 3`) |
| Degenerate distances | all-zero distance matrix or undefined axes → error |
| Outputs | coordinates TSV, distance matrix TSV, top-10 loadings (biplot, Aitchison only), scree plot |
| Host exclusion | optional checkbox "Exclude host and unclassified" (default on) |

Sources: [CODE] `analysis/ordination.py:48-58,76-106,160-215`; `ui/components/ordination_panel.py:92,123,138-139`.

**PERMANOVA** — two separate implementations:

| Where | Distance | Requirement | Reports | Source |
|---|---|---|---|---|
| Ordination panel | the ordination's own distance (Aitchison = Euclidean on CLR coordinates) | ≥ 2 groups, each with ≥ 3 samples | pseudo-F, p, R² (derived), p-value floor for 2 groups; 999 permutations | [CODE] `analysis/ordination.py:238-289` |
| Databases metrics | Bray-Curtis on raw counts | ≥ 3 labelled samples, ≥ 2 groups, each with ≥ 2 | pseudo-F, p; 999 permutations | [CODE] `databases/database_metrics.py:239-271` |

[CODE] For the Aitchison case the distance matrix handed to PERMANOVA is `pdist(result["coordinates"])`
(`analysis/ordination.py:265`), and `coordinates` has been truncated to its first two columns
(`:189,205`). [INFER] the test is therefore run on the 2-D projection, not on the full Aitchison
distance; with exactly 3 samples the two coincide, with more they do not.

**Other statistics:** contaminant flagging by prevalence with a one-sided Fisher exact test against
samples marked as negative controls (`is_contaminant` when p < 0.1) ([CODE]
`databases/database_metrics.py:44,130-172`); taxon prevalence with Wilson 95 % intervals
([CODE] `:177-236`); taxonomic resolution per rank ([CODE] `:84-117`).

**No stored run has ≥ 3 samples with completed taxonomy**, so ordination, PERMANOVA and the
contaminant test have been exercised only by unit tests. [CODE] run inventory §7.

---

## 5. Provenance and reproducibility

### 5.1 What is recorded per run

| Recorded | Where | Source |
|---|---|---|
| `run_id` (UTC timestamp) | `progress.json` | [CODE] `ui/callbacks_progress.py:113-122` |
| `run_config` block: `db_sources` (as requested), `min_unique_size`, `analysis_method`, `branch`, `markers`, `primers {fwd, rev}`, `ncbi_mode` | `progress.json`, written at run start | [CODE] `:74-105` |
| Per-stage records: status, `updated_at`, output file paths, metrics file paths | `progress.json` | [CODE] |
| Time-stamped log lines (`ts`, `level`, `msg`) | `progress.json` `log_lines` | [CODE] |
| Each tool's own log, including its command line and version where the tool prints them | run directory (`step*_*.log`, `cutadapt.log`, `pear.log`, `raxml.raxml.log`, …) | [CODE] §2.0 |
| The marker-filter description (sizes, record counts, primers) | `source_status.csv` `detail` column | [CODE] `taxonomy_runner.py:1522-1523` |
| The exact DADA2 script | `analysis/dada2/dada2_script.R` | [CODE] `dada2_runner.py:470` |

**Not recorded per run** [CODE]: clustering identity, SINTAX cutoff, Trimmomatic/PEAR/DADA2
parameters, BLAST e-value, bootstrap count, phylo mode (they exist only in the browser session store
and, indirectly, in tool logs); tool versions; reference database versions, dates or checksums; the
date of the remote NCBI/BOLD query other than log timestamps; the PARAFISH git commit.

**Real example** — `runs/projeto_1/20260723T072541Z/progress.json`, `run_config` block [CODE]:

```json
{"db_sources": ["silva", "pr2", "bold", "ncbi"], "min_unique_size": 2,
 "analysis_method": "vsearch", "branch": "vsearch", "markers": ["OTHER"],
 "primers": {"fwd": "", "rev": ""}, "ncbi_mode": "remote"}
```

Runs made before 2026-07-17 have no `run_config` block ([CODE] inventory §7; [DOC] `DEVLOG.md:1665`).
`runs/_real/current` and `runs/_fixture/current` — the golden sources — are among them.

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `1cf052f` (F4): a **new** run now records four additional keys
in the `run_config` block, written when the run starts. The seven existing keys are unchanged and
still come first. No stored run has these keys; they appear only in runs made with this branch.

| New key | Content | Failure value |
|---|---|---|
| `parameters` | `enabled_stages`, `trimmomatic`, `pear`, `vsearch`, `dada2`, `databases` (`db_sources`, resolved `marker`, `sintax_cutoff`, `sintax_randseed`, `consensus_method`, and `ncbi` / `bold` sub-blocks when requested), `phylo` (`mode`, `mafft`, `tree`) — the dictionaries as the run will use them | `"unavailable"` |
| `tool_versions` | tool → version, resolved at run start from the **configured** binaries: vsearch, blastn, fastqc, mafft, FastTree, RAxML-NG, Rscript, PEAR, cutadapt, Trimmomatic, boldigger3, and DADA2 when the ASV branch runs | per tool `"unresolved"` |
| `references` | for each **local** reference the run will use (SILVA, PR2 FASTA, PR2 taxonomy, custom FASTA, local BLAST prefix): `source`, `file`, `size_bytes`, `mtime` (UTC), `sha256`. Recorded only when the Databases stage is enabled and the source requested | `sha256: "unavailable"` |
| `parafish_commit` | `git rev-parse HEAD` | `"unknown"` |

[CODE] `ui/callbacks_progress.py` (`_run_provenance`, `_recorded_parameters`, `_recorded_tool_versions`,
`_recorded_references`); `core/export.py` (`reference_fingerprint`, `git_commit`);
`tests/test_run_provenance_recorded.py` (23 tests).

- The sha256 is cached in `<reference>.sha256.json` beside the file and reused while size and mtime
  are unchanged; the cache is git-ignored. In a read-only directory it is not written and the file
  is hashed again next run. [CODE] `core/export.py`.
- Each new key is built in its own `try/except`; none can fail a run. [CODE] `_run_provenance`.
- Measured cost at run start with no reference to hash: 0.56 s. [CODE] measured 2026-10-12. Hashing
  the 824 MB SILVA the first time was not timed.
- Example of `tool_versions` recorded on the development host today [CODE]: vsearch 2.30.0, blastn
  2.17.0+, mafft 7.453, RAxML-NG 1.2.2, Rscript 4.5.0, PEAR 0.9.11, cutadapt 5.2, Trimmomatic 0.39,
  boldigger3 3.0.3, fastqc `unresolved`, FastTree `unresolved`.
- **Still not recorded:** the version or date of the remote NCBI `nt`, of BOLD and of NCBI Taxonomy;
  whether the working tree had uncommitted changes (only the commit is recorded); a checksum of a
  local BLAST database (it is a set of files, so its entry has `sha256: "unavailable"`); the derived
  marker-filtered reference actually loaded by SINTAX (the **source** file is what is fingerprinted).
- The `databases` values are computed by the same expressions as the Databases stage, written in two
  places rather than shared, by decision (no restructuring of the working stage). [CODE] comment in
  `_recorded_parameters`.
- The "Not recorded per run" paragraph above describes runs made **before** this commit.

### 5.2 Exportable run report

| Item | Fact | Source |
|---|---|---|
| Format | a ZIP built in memory, served at `/download/run_zip` and by "Download all (zip)" on Results | [CODE] `core/export.py:156-215`; `app.py:59` |
| Content | `PROVENANCE.txt` + every stage's output files + the persisted metric tables, under `<run_id>/<stage>/…` | [CODE] `core/export.py:163-212` |
| `PROVENANCE.txt` | plain text rendered from the `run_config` block: Reference sources, Marker(s), Analysis method, Branch, `min_unique_size`, primers, NCBI mode, then tool versions | [CODE] `core/export.py:110-153` |
| Tool versions | resolved **at export time** from the configured binaries, not recorded at run time; when none resolve the file says "known gap" | [CODE] `core/export.py:144-152`; `ui/callbacks_results_summary.py:440-460` |
| On-screen reproducibility block | Results page: run id, directory, method, databases, tools + versions, parameters by section, ordination settings, rarefaction depth | [CODE] `ui/callbacks_results_summary.py:462-507` |
| HTML or PDF report | **none** | [CODE] no report generator in the repository |

Example rendered today from the stored COI run's `run_config` (tool versions unresolved in this
shell) [CODE] `core.export.render_provenance`:

```
PARAFISH — run provenance
==============================================
Recipe
----------------------------------------------
  Reference sources : SILVA, PR2, BOLD, NCBI
  Marker(s)         : OTHER
  Analysis method   : vsearch
  Branch            : vsearch
  min_unique_size   : 2
  Primers  fwd      : —
  Primers  rev      : —
  NCBI              : remote
Tool versions
----------------------------------------------
  (not resolved in this environment — known gap)
```

[INFER] Because versions are read at export time, a ZIP exported after a tool upgrade reports the
new version, not the one that produced the data.

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Commit `1cf052f` (F4): for a run that recorded them (§5.1),
`PROVENANCE.txt` now uses the **recorded** tool versions, marked "(recorded when the run started)",
and ignores the export-time ones; it also gains three sections — "PARAFISH commit", "References"
(file, size, mtime, sha256) and "Parameters" (one `key = value` line per setting, nested keys
dotted). For a run **without** the new keys the output is unchanged, and a test pins that against
the previous rendering byte for byte. [CODE] `core/export.py` (`render_provenance`,
`_recorded_provenance_lines`); `tests/test_run_provenance_recorded.py::RenderProvenancePrefersWhatTheRunRecordedTests`.
The on-screen reproducibility block of the Results page is **unchanged**: it still reads parameters
from the browser session and resolves versions when rendered. The `/download/run_zip` Flask route
calls `build_run_zip` without a provenance text ([CODE] `app.py:86`), so a ZIP fetched from that URL
directly, rather than with the "Download all (zip)" button, contains no `PROVENANCE.txt` — as before.

### 5.3 Golden layers

**What they are** ([CODE] `scripts/golden.py:1-21,36-157`): for each closed stage, the contract files
are copied to `golden/<dataset>/<stage>/` and bulk outputs are recorded as MD5 checksums in
`CHECKSUMS.txt`. `python scripts/golden.py check <stage|all> --from <run_dir> --dataset <name>`
diffs a **stored run directory** against the frozen copy.

| Stage | Files frozen |
|---|---|
| qc | `qc/metrics/qc_{length_distribution,overview,per_base_quality,module_status}.csv` |
| preprocess | `read_tracking.csv`, `assembled_length_distribution.csv`; checksum of `*.assembled.fastq` |
| analysis | `otu_table.tsv`, `otu_id_map.csv`, `otu_matrix.csv`, metrics tables (and the `asv` equivalents) |
| taxonomy | the taxonomy deliverables incl. `source_status.csv` |
| export | `phyloseq/{tax_table,otu_table,sample_data}.tsv`, `taxa_tables/taxa_table_{class,genus,species}.csv`; checksum of `table.biom` (ignoring `date`) |
| phylo | `mafft/aligned.fasta`, `mafft/reversed_asvs.csv`, `fasttree/tree.nwk`, `raxml/raxml.raxml.{support,bestTree,bestModel}` |

| Dataset | Source run | Layers frozen | In git? | Last frozen |
|---|---|---|---|---|
| `fixture` | `runs/_fixture/current` (2 × 5 000 read pairs, 18S, OTU) | qc, preprocess, analysis, taxonomy, export, phylo | yes | last commits touching each: preprocess/export 2026-07-17; qc/analysis/taxonomy 2026-07-21; phylo 2026-07-22 [CODE] `git log` |
| `fixture_asv` | `runs/_fixture_asv/current` (same reads, DADA2) | analysis, taxonomy, export | yes | analysis/taxonomy 2026-07-21; export 2026-07-22 [CODE] `git log` |
| `real` | `runs/_real/current` (full depth, 18S, OTU) | qc, preprocess, analysis, taxonomy, export, phylo | **no (untracked)** | directory mtimes: preprocess 2026-07-16; export 07-17; qc, analysis, taxonomy 07-21; phylo 07-22 [CODE] `ls -la golden/real` |

**Status today** [CODE] run on 2026-10-12: `golden.py check all` → "no regressions" for `fixture`
(6 PASS), `real` (6 PASS) and `fixture_asv` (3 PASS, 3 SKIP "not frozen yet").

**Limits that must be stated:**

- `golden.py check` **never re-runs the pipeline**; it compares a stored run to frozen files. "No
  regressions" means those files are unchanged, not that the current code reproduces them.
  [CODE] `scripts/golden.py:57-95`; [DOC] `DEVLOG.md:460-493`.
- The taxonomy golden is **known-stale** with respect to the current code: the frozen files still
  have the id column `ASV_ID`, unrounded identities, no abundance columns, no `sintax_status`
  column, and Portuguese text in `source_status.csv` `detail`. [CODE] the header of
  `runs/_real/current/databases/otu/taxonomy_consensus.csv` begins `ASV_ID,…` while the code writes
  `feature_id` (`taxonomy_runner.py:495`). [DOC] `DEVLOG.md:170-299` and the 2026-08-03 entry.
- Two tests are skipped pending that re-freeze ([CODE] `tests/test_taxonomy_summary.py:145,152`).
- **CONTRADICTION:** [DOC] `README.md` ("Reproducibility") says contracts are frozen "across three
  datasets (`fixture`, `fixture_asv`, `real`)"; [CODE] `golden/real` is not in the repository, and
  none of the three source runs (`runs/`) or the fixture FASTQ (`fixtures/`) is tracked, so a fresh
  clone cannot execute any golden check.

### 5.4 `scripts/verify.py`

Usage `python scripts/verify.py <run_dir> [--quiet]`; standard library only; exits non-zero only on
FAIL; reports SKIP when an input is missing ([CODE] `scripts/verify.py:1-33`). Checks, run per feature
set present:

| Check | What it asserts |
|---|---|
| (a) | read totals agree across `otu_table`, `otu_matrix`, `chimera_stats`, `taxonomic_resolution` and every `taxa_table_*` |
| (b) | `otu_id_map.size` equals the feature's row sum in the abundance table |
| (c) | `taxa_table_class` agrees with the class column of `tax_table` |
| (d) | no invented middle ranks (Opisthokonta, Amorphea, Obazoa, Craniata) |
| (e) | `source_status` agrees with the assignments each source made and with raw outputs on disk |
| (f) | a Linnaean source reached the consensus (SKIP if none was requested) |
| (g) | the Results page can read this run |

Results today [CODE] 2026-10-12:

| Run | Result | Notable lines |
|---|---|---|
| `runs/_real/current` | 7 pass, 0 fail, 0 skip | (a) all 11 sources = 154 223; (e) NCBI=9 · PR2=9 · SILVA=9; (f) NCBI 9/9 |
| `runs/_fixture/current` | 7 pass | (a) = 8 107; 4/4 features |
| `runs/_fixture_asv/current` | 5 pass, 2 skip | (a/asv) = 8 458; 3/3 features; (f) skipped: no Linnaean source requested |
| `runs/_coi/current` | 7 pass | (a) = 3 226; 55/55; (e) BOLD=20; (f) BOLD 20/55 |

### 5.5 Test suite

**Current state: could not be executed.** [CODE] `.venv/bin/python` is a symlink to
`/home/joaofaria/miniconda/bin/python`, which does not exist on 2026-10-12; no other environment on
the machine has the dependencies (checked `miniconda3` base and its four environments). Installing
them would have changed the environment, which the task forbade.

**Last recorded result:** 581 passed, 2 skipped, 19 subtests (2026-08-03, commit `90dc89e`)
[DOC] `DEVLOG.md` 2026-08-03 entry and the commit message of `90dc89e`. [INFER] still valid, since
`git diff 90dc89e HEAD` is empty.

53 modules, 583 test functions ([CODE] AST count of `def test_*`):

| Module | Tests | What it verifies (module docstring) |
|---|---|---|
| `test_analysis_asv_residual.py` | 2 | selector-free Analysis callbacks render the resolved feature set, not "otu" |
| `test_analysis_metrics.py` | 22 | chimeras, rarefaction, alpha diversity, UpSet data |
| `test_analysis_occupancy_asv.py` | 2 | occupancy callback restore |
| `test_analysis_restore.py` | 3 | Analysis page reconstructs from disk in a fresh session |
| `test_analysis_stats_restore.py` | 5 | Analysis stat cards survive a restored session |
| `test_bold_positive_parsing.py` | 2 | BOLD positive-result parsing against the real BOLDigger3 3.x schema |
| `test_callbacks_smoke.py` | 2 | every server callback runs with empty inputs; no output claimed twice |
| `test_community_phylo_diversity.py` | 9 | community ecology and phylogenetic diversity numerics |
| `test_custom_db_compare.py` | 34 | custom-database BLAST comparison and upload validation |
| `test_cutadapt_runner.py` | 19 | primer removal, presets, degenerate bases, PE/SE, passthrough |
| `test_database_metrics.py` | 40 | resolution, confidence, decontamination, prevalence, BIOM, phyloseq |
| `test_databases_hits_consensus_asv.py` | 3 | hits and consensus callbacks on the ASV golden |
| `test_databases_qc.py` | 39 | database and QC runner contracts, source status, marker applicability |
| `test_databases_restore.py` | 3 | Databases page restore from disk |
| `test_dropdown_theme_tokens.py` | 3 | dropdown theme matches the installed Dash bundle |
| `test_export_provenance.py` | 5 | the run ZIP carries the recipe (`PROVENANCE.txt`) |
| `test_export_zip.py` | 4 | ZIP export never follows a stray non-path string |
| `test_feature_sets.py` | 14 | OTU and ASV sets generated alike, reported apart |
| `test_ncbi_lineage.py` | 7 | taxid → lineage resolver, cache, merged taxids |
| `test_ordination.py` | 15 | compositional correctness; refusal to plot with n < 3 |
| `test_partial_analysis.py` | 6 | a half-successful analysis still delivers its half |
| `test_phase2_contracts.py` | 5 | core data model and validation helpers |
| `test_phase3_security_phylo.py` | 17 | download path traversal, BOLD credentials, phylo keys |
| `test_phylo_golden.py` | 5 | phylo golden contract, per-set tree resolution |
| `test_phylo_input_guard.py` | 4 | the phylo guard refuses at the same floor the aligner enforces |
| `test_phylo_metrics.py` | 27 | alignment QC, branch support scaling, UniFrac |
| `test_pipeline_id_fixes.py` | 16 | regression tests for id handling bugs |
| `test_pipeline_launch_guard.py` | 7 | a pipeline starts only on an explicit click, never twice |
| `test_preprocess_metrics.py` | 19 | read-tracking waterfall parsed from runner logs |
| `test_preprocess_restore.py` | 3 | Preprocess page restore |
| `test_progress_stage_ordering.py` | 5 | stage badges do not lag behind the pipeline |
| `test_project_flow.py` | 3 | Project page states what happened and what to do next |
| `test_qc_metrics.py` | 17 | FastQC module parsing, persisted metrics |
| `test_qc_restore.py` | 2 | QC page restore |
| `test_qc_stats_restore.py` | 3 | QC KPI cards restore |
| `test_reference_filter.py` | 16 | marker-filtered references |
| `test_results_summary.py` | 22 | cross-stage read fate, sample dashboard, target report |
| `test_run_dir_single_source.py` | 5 | one selector for "the newest run to show" |
| `test_run_identity.py` | 12 | a stage record belongs only to the run that produced it |
| `test_run_provenance.py` | 2 | a run's config is recorded in its own `progress.json` |
| `test_run_registry.py` | 8 | UI facts come from disk, not session state |
| `test_sample_comparison.py` | 7 | cross-sample comparison metrics |
| `test_scientific_parameter_pins.py` | 5 | literal pins: identity 0.97, SINTAX cutoff 0.8, PEAR `min_overlap` 10, Trimmomatic `minlen` 50, `min_unique_size` 1 |
| `test_server_unresponsive_diagnosis.py` | 5 | "server did not respond" diagnosis banner |
| `test_setup_defaults_single_source.py` | 5 | Setup numeric fields agree with `RunConfig` defaults (11 mapped fields only) |
| `test_stage_runners.py` | 12 | runner contracts for preprocess, analysis, taxonomy, phylo |
| `test_startup_reconciliation.py` | 12 | orphan runs do not masquerade as the active run |
| `test_taxonomy_phylo.py` | 4 | taxonomy and phylogenetic contracts |
| `test_taxonomy_rank_fixes.py` | 48 | regression tests for rank mapping (PR2, SILVA, class alias, SINTAX parsing) |
| `test_taxonomy_summary.py` | 27 | per-feature summary, both evidence families side by side (2 skipped) |
| `test_ui_language.py` | 2 | **AST test**: no Portuguese word in any string constant of `ui/ core/ analysis/ databases/ preprocess/ phylo/ qc/` (docstrings excluded); Downloads table column id, cell value and `filter_query` agree |
| `test_ui_presentation.py` | 12 | table contrast (WCAG), small-n plots, id hygiene |
| `test_verify_asv.py` | 7 | `verify.py` sees an ASV run |

**Limit of the AST language test.** It matches a fixed word list, so Portuguese strings built from
words outside the list still pass. Found by reading the code today [CODE]:
`"VSEARCH SINTAX foi terminado pelo sistema (sinal …)"` (`taxonomy_runner.py:984`),
`"erro a preparar a base: …"` (`:1504,1537`), `"Base customizada (BLAST)"` / `"(SINTAX)"` as the
source name written to `source_status.csv` (`:1597,1611,1617`), `"PCA de Aitchison (CLR) —
recomendado"` (`analysis/ordination.py:63`), `"primers removidos"` (`callbacks_progress.py:284`),
`"(filtro)"` (`:442`), `"recortada aos primers …"` (`reference_filter.py:414`),
`"PEAR requer r1_trimmed e r2_trimmed…"` (`pear_runner.py:34`), `"Modo PE requer R2"`
(`core/models.py:56`), BOLD validation messages (`bold_db.py:209-213,570`),
`"(offline ou erro de rede?)"` (`ncbi_lineage.py:144`), `"economica"` (`mafft_runner.py:68`).

---

## 6. User interface

Ten pages, all rendered at once and shown/hidden by route ([CODE] `ui/layout.py:55-64`). Sidebar order
([CODE] `ui/components/sidebar.py:18-29`): Project, Samples, Setup, QC, Preprocess, Analysis,
Databases, Phylo, Progress, Results. Pages other than Project, Samples and Setup are disabled in the
sidebar while their step is `pending`.

Status was assessed **statically**: an element is counted as wired if its id appears as an `Output`,
`Input` or `State` of a registered callback ([CODE] script over `ui/pages/*.py` and
`ui/callbacks_*.py`), plus the smoke test's coverage. The application was **not opened in a browser**
for this document. [UNKNOWN] visual and interactive correctness today.

| # | Page (route) | Shows | Controls | Status |
|---|---|---|---|---|
| 1 | Project (`/project`) | counts of projects/samples/modes, recent projects, active-project state, "next step" link | project name, notes, "Create project", "Load example" | **FUNCTIONAL** |
| 2 | Samples (`/samples`) | active-project strip, upload status, sample table (marker, mode, selected, control flag, group) | upload (`.fastq, .fastq.gz, .fq, .fq.gz, .zip, .tar.gz, .rar`, max 2 GB), "Clear samples", per-sample checkboxes, control and group fields | **FUNCTIONAL**. There is **no marker selector**: the marker shown is inferred from file name or primer preset ([DOC] `DEVLOG.md:731-735`). |
| 3 | Setup (`/setup`) | workflow mode, five stage cards with on/off switches and parameters, memory preflight, sample table, run-target note | see §2 and §3; "Threads", "Workers", "Run Pipeline" | **PARTIAL** — five controls have no effect: `setup-mafft-mode`, `setup-adjust-direction` (no callback reads them); `setup-vsearch-sintax`, `setup-species-threshold`, `setup-workers` (read, then unused). `dada2-its-warning` is never shown by a callback. Workflow modes: "Complete" forces all five stages on and locks the switches; "One" and "Special" constrain how many may be on ([CODE] `ui/callbacks_setup.py:43-60,495`). [CODE] §2.12, §3.3, §1.4 |
| 4 | QC (`/qc`) | run banner, four KPI cards (samples, pass rate, reads, GC), run overview table, read-depth plot, per-base quality/length/module plots | "Run QC" | **FUNCTIONAL** |
| 5 | Preprocess (`/preprocess`) | run banner, three KPI cards, read-tracking table, read-fate plot, parameters used, pipeline log | "Run Preprocess" | **PARTIAL** — `preprocess-pipeline-bar` is an empty container (no callback) |
| 6 | Analysis (`/analysis`) | KPI cards (OTUs, ASVs, consensus), per-sample yield, alpha diversity with rarefaction note, ordination panel, chimera stats, DADA2 read tracking, abundance table, diversity/composition plots, shared/exclusive features (UpSet), feature occupancy, downloads | "Run Analysis", abundance mode, occupancy set and filter, ordination controls | **PARTIAL** — `analysis-pipeline-bar` is an empty container |
| 7 | Databases (`/databases`) | feature-set switch, active databases, custom-database upload, consensus table, custom comparison (table, heatmap, donut), ranked hits table, per-feature summary, OTU-vs-ASV comparison, taxonomic resolution, confidence, conflicts, contaminants, prevalence, composition | "Run Databases", feature-set radio, min identity / min coverage / top-N filters, species filter, two CSV downloads, upload (max 500 MB) | **FUNCTIONAL** (all 27 ids wired) |
| 8 | Phylo (`/phylo`) | KPI cards (tips, alignment, support), annotated tree viewer, alignment QC, branch support, phylogenetic alpha diversity | "Run Phylo" | **FUNCTIONAL** (all 12 ids wired; `phylo-unifrac` is a hidden container) |
| 9 | Progress (`/progress`) | overall status, per-stage rows, live log | "Start Pipeline" (the only launcher); 2 s polling | **FUNCTIONAL** |
| 10 | Results (`/results`) | run header with aggregate badge, read funnel, per-stage KPI cards, community composition (rank selector, host exclusion), compact ordination, between-sample comparison, target-taxon report, downloads table + "Download all (zip)", reproducibility block, run log, per-stage/per-sample viewer | stage and sample selectors, rank dropdown, exclusion checklist, target input, "New Run" | **PARTIAL** — three collapsible sections are **empty containers with no callback**: "Alpha diversity per sample" (`results-alpha-table`), "Taxonomy table (per-feature consensus)" (`results-taxonomy-table`), "Detailed read tracking" (`results-tracking-table`). `results-sankey`, `results-dashboard`, `results-stats` are hidden (`display: none`) though they have callbacks. The Downloads table hard-codes the `otu` paths, so a DADA2-only run shows those rows as missing ([CODE] `ui/callbacks_results_summary.py:300-321`). |

### 6.1 Known UI bugs and their state

| Bug | State | Source |
|---|---|---|
| `results-stage-selector` / `results-sample-selector` were empty after a restart, a new tab, or when a stage had failed | **Fixed**: options now come from `progress.json` via `core/run_registry.py`; failed/skipped stages are listed with their status; stages without per-sample output say so | [DOC] `DEVLOG.md:3058-3080` (2026-07-14); [CODE] `ui/callbacks_results.py:288-357` |
| `store-samples`: samples attributed to a different, still-active project | **Not a code bug** per the audit; mitigated by an active-project strip on Samples and Setup | [DOC] `DEVLOG.md:3003-3018` (2026-07-15) |
| Results and other pages could pick different runs | **Fixed** (one primitive, `newest_complete_run`) | [DOC] `DEVLOG.md:785`; [CODE] `core/run_registry.py:72-121` |
| An OOM-killed run blocked or blanked later sessions | **Fixed** by startup reconciliation | [CODE] `core/run_registry.py:136-171` |
| Dropdowns illegible on the dark theme | **Fixed** (theme tokens) | [DOC] `DEVLOG.md:740` |
| Analysis "Consensus" card could never show "ready" | **Fixed** | [DOC] `DEVLOG.md:107-132` |
| Empty Results sections, dead Setup controls, empty pipeline bars | **Open** | [CODE] table above |
| Session stores (`store-*-outputs`, `store-run-config`) are `storage_type="session"`; `store-project`, `store-samples` are `local` | by design | [CODE] `ui/layout.py:38-46` |

---

## 7. Data and runs available

### 7.1 Inventory of `runs/`

[CODE] generated from every `progress.json` on 2026-10-12. 107 run directories in 14 project folders
(plus an `exports/` folder). The per-feature-set layout (`databases/<set>/taxonomy_consensus.csv`)
exists in 8 runs; 5 older runs in `pr/` have the flat pre-feature-set layout
(`databases/taxonomy_consensus.csv`).

| Folder | Runs | Period | What it holds | Useful for the thesis? |
|---|---|---|---|---|
| `_real/current` | 1 | 2026-07-17 | **Full-depth 18S V4**, 2 samples (APTL64B, APTL64F), OTU, primers removed, SILVA+PR2+NCBI, phylo. 9 OTUs. Source of `golden/real`. | **Yes — the reference 18S result** |
| `_fixture/current` | 1 | 2026-07 | same samples subsampled to 5 000 pairs each; 4 OTUs; SILVA+PR2+NCBI; phylo | Yes — reproducibility demonstration |
| `_fixture_asv/current` | 1 | 2026-07 | fixture through DADA2; 3 ASVs; SILVA+PR2 | Yes — only completed ASV run |
| `_coi/current` | 1 | 2026-07-22 | 2 × 2 000 pairs, marker forced to COI, **BOLD only**, no primer removal; 55 OTUs, BOLD assigned 20 | Yes — only proof of BOLD's positive path. Caveat §7.3 |
| `projeto_1/` | 15 | 2026-07-15 → 08-03 | the 18S development runs and the **three COI runs** (§7.2) | Yes — COI case |
| `pr/` | 17 | 2026-07-08 → 07-14 | 18S development runs, 2 samples, 7–12 OTUs; several with phylo; taxonomy in the flat layout in 5 runs and in the per-set layout in `20260714T152744Z` | Historical only — they predate the taxonomy, primer and abundance fixes of 2026-07-15/16 ([DOC] `DEVLOG.md`) |
| `new/` | 4 | 2026-07-15 → 07-17 | one complete 2-sample run (7 OTUs); one interrupted **3-sample** run with 138 OTUs (`20260717T144357Z`, stopped in Analysis) | No |
| `novo/` | 24 | 2026-06-01 → 06-09 | early runs, Databases skipped; up to 6 samples | No |
| `project/` | 41 | 2026-05-23 → 07-08 | earliest runs; one 10-sample run (QC→Analysis only, no feature files kept); `20260708T091748Z`: 3 samples, 1 107 OTUs, Databases failed | No |
| `manual_run_validation`, `qc_validation_run`, `visual_validation_run`, `test_project`, `upload_test` | 5 | 2026-05 | validation scaffolding, no `progress.json` | No |

No stored run has a complete taxonomy for **more than two samples**. No stored run used both OTU
and ASV. [CODE]

### 7.2 The real COI runs

All three use the **same single sample**, `aptl64b_01` (input files `APTL64B_R1.fastq.gz` /
`APTL64B_R2.fastq.gz`, 207 528 read pairs), in `runs/projeto_1/`. [CODE] each run's `progress.json`
and `preprocess/metrics/read_tracking.csv`.

| Run id | Date | Primers removed? | `min_unique_size` | OTUs | Completed? | Taxonomy? |
|---|---|---|---|---|---|---|
| `20260723T072541Z` | 2026-07-23 | **No** | 2 | **90** | Yes (all five stages) | Yes: NCBI 83/90, BOLD 17/90; SILVA and PR2 skipped (RAM); custom BLAST failed |
| `20260724T082931Z` | 2026-07-24 | **Yes** (cutadapt, `COI_Leray` sequences; 92.32 % of reads carried the forward primer) | 2 | **112** | **No** — interrupted in Databases | **No** |
| `20260803T111144Z` | 2026-08-03 | **No** | 1 | **1 464** | **No** — Databases never finished (still marked `running`) | **No** |

- Marker recorded in all three: `OTHER` ([CODE] `run_config.markers`) — see the contradiction in §3.2.
- **Before/after primer trimming outputs both exist** for the analysis stage: `analysis/pooled/*`
  (pooled reads, dereplicated reads, raw centroids, chimeras, OTUs, tables) is present in both
  `20260723T072541Z` (before) and `20260724T082931Z` (after), and the after-run keeps
  `preprocess/aptl64b_01/*_primer.fastq.gz` and `cutadapt.json`. **Taxonomy exists only for the
  before-run.** [CODE] directory listings.
- The after-run predates `--discard-untrimmed` (its `primer_trimmed_reads` equals the input,
  207 528). [CODE] `read_tracking.csv`; [DOC] `DEVLOG.md:648-655`.
- **A COI run with primers removed and `min_unique_size = 1` — the platform's current defaults — has
  never been made.** [CODE] inventory.
- Origin of the sample: [UNKNOWN]. The repository does not document what APTL64B/APTL64F are (host
  species, tissue, locality, project, sequencing provider). [DOC] `scripts/run_coi_bold_demo.py`
  docstring refers to a "Micromesistius host" known from the 18S runs.

### 7.3 Test datasets and fixtures

| Dataset | Location | Content | Purpose | Tracked? |
|---|---|---|---|---|
| 18S full depth | repo root `APTL64B_R{1,2}.fastq.gz`, `APTL64F_R{1,2}.fastq.gz`; copies in `runs/projeto_1/20260715T151700Z/input/` | 101 190 and 90 220 read pairs | source of `_real` and of the fixture | no |
| Fixture | `fixtures/APTL64{B,F}_R{1,2}.fastq.gz` | 5 000 pairs each, drawn uniformly with a fixed seed, mates kept in lockstep | golden `fixture` and `fixture_asv`; `scripts/run_pipeline.py` default input | no (`.gitignore:51`) |
| COI demo input | `runs/_coi/current/input/` | 2 000 pairs each, named APTL64B/APTL64F | BOLD demonstration | no |
| Custom DB demos | `data/demo_custom_db_blast.fasta`, `data/demo_custom_db_sintax.fasta` | ~3 kB each | custom-database feature | [UNKNOWN] not checked |
| `parafish_test_fastq.zip` | repo root | two 20-byte files | upload test stub, **not** a dataset | no (`.gitignore:43`) |

**CONTRADICTION / caveat on `_coi`.** The COI demo's input files carry the same names as the 18S
samples (APTL64B, APTL64F) and its sample ids are `aptl64b_02` / `aptl64f_01`, the 18S ids
[CODE] `scripts/run_coi_bold_demo.py:26-30`. [UNKNOWN] whether those files are COI amplicons from the
same specimens or something else; the script's docstring says they are COI reads subsampled
beforehand, and BOLD did return identifications for 20 of 55 OTUs.

**Supervisor's data** — `dados stora teresa/SP_Novogene dados/SP_Novogene Dados metabarcoding Workshop/`
[CODE]: three marker folders (`18S`, `COI`, `16S Euk`), each with `01.RawData/` holding **11 samples**
(e.g. `OF.216F.i`, `OF.216G.i`, `OF.216M.i`, `OF.224F.i`, …) × 5 files (`.raw_1`, `.raw_2`, `_1`,
`_2`, `.extendedFrags`); 55 files per marker; 356 MB total; untracked. The accompanying primer list
is titled "Primers utilizados para metabarcoding dos tecidos de peixe SP". The directory was created
after the last code commit. [UNKNOWN] whether any of it has been run through PARAFISH — no run
directory references these sample names.

### 7.4 Guadiana dataset (PRJNA1334324)

**Absent.** [CODE] the strings `PRJNA`, `Guadiana` and `guadiana` occur nowhere in the working tree
outside `.venv/` and `.git/` (code, docs, DEVLOG, scripts, configuration). There is no download
script, no accession list, no sample sheet and no run. [UNKNOWN] whether the data exist elsewhere.

---

## 8. Decisions log

Dates are those of the DEVLOG entry. "Origin" gives the DEVLOG heading line in the current file.

| Decision | Justification as recorded | Date | Origin |
|---|---|---|---|
| OTU clustering identity **0.97** | No measured justification is recorded; it is the conventional species-level threshold and is described in the UI as "0.97 = species-level OTU". Pinned by test on 2026-07-29 because a mutation sweep showed the value could be changed to 0.90 with the suite still green. | 2026-07-24 / 07-29 | [DOC] `DEVLOG.md:299-344`; [CODE] `page_setup.py:166`; commit `a12c476` |
| SINTAX cutoff **0.8** | No measured justification recorded; pinned for the same reason. Reliability bands (0.95 / 0.85) were set above it so the grade is not true by construction. | 07-24 / 07-29 | [DOC] `DEVLOG.md:299-344`; [CODE] `taxonomy_runner.py:566-571` |
| `min_unique_size` **1** (singletons may seed OTUs) | A floor of 2 does not remove sequencing errors (those are within 3 % of a parent and are absorbed by clustering); it removes reads that are globally unique and > 3 % divergent — the profile of a rare taxon. On the 18S fixture a floor of 2 lost one read that at floor 1 founds an OTU. Cost on the 18S real run: clustering 57 s → 247 s (4 232 uniques instead of 1 475). | 2026-07-16 | [DOC] `DEVLOG.md:2439-2498` |
| Abundance comes from the mapping table, not from cluster sizes | Two accountings (147 010 vs 157 527 reads) existed because dereplication dropped singletons while mapping counted them. | 2026-07-16 | [DOC] `DEVLOG.md:2795-2852` |
| SINTAX `--randseed 42` | Unseeded SINTAX gave different bootstraps on identical input (0.99 vs 0.98), which can move a call across the cutoff. Any fixed non-zero value would do. | 2026-07-16 | [DOC] `DEVLOG.md:2249-2271` |
| RAxML-NG seed 42 | Recorded as what makes the trees byte-reproducible. | 2026-07-22 | [DOC] `DEVLOG.md:1076-1132` |
| Remote NCBI uses **megablast**, not blastn or dc-megablast | Measured on 9 OTUs against remote `nt`, same day: `blastn -word_size 11` and `dc-megablast` were refused (SIGXCPU) on all 9; megablast returned 123 hits for 9/9 in ~3 min. Accepted cost: megablast's 28 bp exact seed can miss a hit at ~90 % identity. | 2026-07-16 | [DOC] `DEVLOG.md:2230-2247,2292-2333` |
| Sensitive search only with a local database (`blastn`, `word_size 11`), opt-in | The public quota refuses it; a local db is the user's own CPU. Local BLAST is off by default so the source is one a reviewer can reproduce. | 2026-07-16 | [DOC] `DEVLOG.md:2244-2247,2330-2333` |
| NCBI lineage via **Entrez `efetch`**, cached per taxid | `blastn` reports only a binomial and a taxid, leaving kingdom…family blank. **No comparison with a local taxdump is recorded** — the alternative is not discussed in the DEVLOG. | 2026-07-15 | [DOC] `DEVLOG.md:3020-3054`; [UNKNOWN] why not taxdump |
| SILVA and PR2 contribute only domain and genus/species, not the middle ranks | Their backbones are not Linnaean; mapping them positionally invented ranks (`Phylum = Opisthokonta`, `Class = Craniata`, `Genus = Teleostei`) and created a false conflict on every OTU. Cost: PR2's real protist classes are lost. | 2026-07-16 | [DOC] `DEVLOG.md:2736-2793` |
| The Linnaean backbone depends on NCBI (or BOLD) | Consequence of the row above, acknowledged as an open weakness: without NCBI reachable there is no phylum…family. `verify.py` check (f) was added to fail such runs. | 2026-07-16 | [DOC] `DEVLOG.md:2283-2290,2335-2346` |
| Source status is decided by actual assignments, with `not_applicable` distinct from `no_hits` | BOLD reported `ok` with zero identifications; a refused BLAST reported `no_hits`. A source that cannot carry the marker must not read as "found nothing". | 2026-07-16 / 07-24 | [DOC] `DEVLOG.md:2690-2734,665-738` |
| Marker-filtered references for SINTAX | Full SILVA needs ~2.5 GB and was always skipped on a 2.8 GB machine. Cutting to Eukaryota + the amplicon gave SILVA 824 → 27 MB and PR2 334 → 85 MB; a real run classified 12/12 OTUs with both in 7.8 s and 258 MB peak RSS. | 2026-07-16 | [DOC] `DEVLOG.md:2898-2943` |
| Primer removal with cutadapt on raw reads, before Trimmomatic | The earlier fixed-length strip was blind and a no-op unless fields were filled. On the 18S run 100 % of reads still carried the primer and one species was split across OTU1 and OTU2 by the 5' spacer. | 2026-07-15 / 07-16 | [DOC] `DEVLOG.md:2972-3001,2645-2688` |
| `--discard-untrimmed` on by default | Provenance: a read with no primer cannot be shown to be an amplicon of the target region. Cost on the COI sample ~7.7 % of input reads, 0.19 % of reads reaching clustering. | 2026-07-24 | [DOC] `DEVLOG.md:648-655` |
| Trimmomatic `MINLEN` **50**, not 200 | On the COI sample 200 would discard 14 536 merges (8.00 %), 92.82 % of them valid 300–330 bp amplicons, to remove 2 306 frail merges (1.27 %). | 2026-07-24 | [DOC] `DEVLOG.md:495-536` |
| PEAR `min_overlap` left at 10 | Raising it to ~50 is identified as the correct lever but would re-baseline the golden; left open. | 2026-07-24 | [DOC] `DEVLOG.md:525-532` |
| Aitchison PCA (CLR) as the default ordination; pseudocount 0.5 | Counts are compositional; PCA on raw counts makes PC1 track sequencing depth. Fixed pseudocount chosen over Bayesian multiplicative replacement (zCompositions) to avoid an R dependency; measured distortion < 1 % of coordinates for a 50× depth difference. | 2026-07-14 | [DOC] `DEVLOG.md:3082-3104`; `analysis/ordination.py:1-48` |
| No ordination below 3 samples | With two samples there is one axis and the separation is an artefact. | 2026-07-14 | [DOC] `DEVLOG.md:3094-3104` |
| PCoA on Bray-Curtis/Jaccard/UniFrac offered as alternatives | "What the literature mostly reports", and UniFrac is the only option using the phylogeny. | — | [DOC] `analysis/ordination.py:23-25` |
| DADA2 does not assign taxonomy by default | Parity with OTUs (same evidence columns) and memory (`assignTaxonomy` loads the whole reference). | — | [DOC] `analysis/dada2_runner.py:1-20` |
| OTU and ASV never share a table | An OTU and an ASV are different units; comparison is by taxon only. | — | [DOC] `ui/callbacks_progress.py:526-529,786-788` |
| Pipeline starts only on an explicit "Start Pipeline" click; one run at a time | Remote queries cost hours, and two concurrent runs on a 3 GB machine are a guaranteed OOM. | 2026-07-23 | [DOC] `ui/callbacks_setup.py:605-611`; `callbacks_progress.py:994-1006` |
| UI state is read from disk | Session stores were empty after a restart, new tab or finished run. | 2026-07-14 / 07-17 | [DOC] `DEVLOG.md:3058-3080` |
| Golden contracts per stage; a stage is frozen only when all its deliverables are closed | An early freeze was misread as a regression. | 2026-07-16 | [DOC] `DEVLOG.md:2491-2496` |
| Rule: judge an external tool by its output, not its exit code; and an empty output is not automatically a failure | blastn exits 0 when NCBI refuses the search. | 2026-07-16 | [DOC] `DEVLOG.md:2193-2219`; [CODE] `core/process.py` |
| Rule: a scientific parameter must be pinned by a test that fails when the value changes | Mutation sweep: 3 of 5 parameters could be changed with the suite green. | 2026-07-24 | [DOC] `DEVLOG.md:299-344,390-441` |
| User-facing text in English | Consistency with an English README and thesis. | 2026-08-03 | [DOC] `DEVLOG.md:3-69` |

---

## 9. Gap analysis against the pre-thesis

The pre-thesis text itself is not in the repository; the items below are the ones listed in the
request. [UNKNOWN] the exact wording that was promised.

| Promised | Status | Evidence |
|---|---|---|
| End-to-end cross-platform GUI | **PARTIAL** | End-to-end GUI: implemented (10 pages, §6). Cross-platform: delivered only through Docker, and the Docker image has never been built or run ([CODE] §1.3). Development and all stored runs were on Linux/WSL2. No native Windows or macOS support in the code (`java`, POSIX paths, `/proc/meminfo` in `dada2_runner.py:66-75`). |
| Demultiplexing | **NOT IMPLEMENTED** | [CODE] no occurrence of demultiplexing anywhere; input must be per-sample FASTQ. |
| Quality filtering | **IMPLEMENTED** | Trimmomatic, VSEARCH `--fastq_filter`, DADA2 `filterAndTrim` (§2.5, §2.7, §2.9). |
| Denoising | **IMPLEMENTED** (DADA2); UNOISE absent | §2.9. Exercised only on the 5 000-read fixture (`runs/_fixture_asv`, 3 ASVs). |
| Chimera removal | **IMPLEMENTED** | `uchime_denovo` on OTU centroids; `removeBimeraDenovo` for ASVs (§2.8, §2.9). No reference-based chimera check. |
| OTU and ASV in parallel | **PARTIAL** | Both branches and the `both` mode are implemented and unit-tested (`tests/test_feature_sets.py`, `test_partial_analysis.py`); they run sequentially, not concurrently; **no stored run contains both sets** ([CODE] §7.1). |
| Curated, controlled and versioned reference libraries for fish and fish parasites | **NOT IMPLEMENTED** | The platform uses general public references (SILVA 138.2, PR2 5.1.1, NCBI `nt`, BOLD) plus a user-uploaded custom FASTA. No fish- or parasite-specific library ships with it; no reference is checksummed or version-recorded per run ([CODE] §3.1, §5.1). The marker filter narrows SILVA/PR2 to Eukaryota + amplicon, which is a size reduction, not curation. |
| Automatic taxonomic assignment | **IMPLEMENTED** | SINTAX, BLAST, BOLDigger3, consensus (§3). Multi-source consensus completed on 2 samples (§3.8). |
| Phylogenetic placement | **NOT IMPLEMENTED** | [CODE] no placement tool (EPA-ng, pplacer, SEPP); no reference tree. |
| Tree construction and visualisation | **IMPLEMENTED** | MAFFT, FastTree, RAxML-NG; tree viewer on the Phylo page (§2.12, §4.1). OTU trees only in stored runs. |
| Benchmark with real data (sensitivity, false positives, negative controls) | **NOT IMPLEMENTED** | No mock community, no spike-in, no truth set, no sensitivity or false-positive figure anywhere ([CODE] no such script, test or result). A prevalence-based contaminant test exists in code (§4.2) but **no stored run has a sample flagged as a negative control**. The only real-data comparison is the before/after primer-trimming count on one COI sample. |
| User documentation | **PARTIAL** | `README.md` (installation, configuration, outputs, limitations) and `docs/setup_docker.md`. No user manual, no screenshots, no tutorial dataset in the repository. Several README statements are contradicted by the code (§12). |
| Recommended default parameters | **PARTIAL** | Defaults exist and five are pinned by tests; but the UI and headless defaults disagree for sliding-window quality, bootstrap replicates and BOLD mode (§2.5, §2.12, §2.11), and no document justifies 0.97 or 0.8 from data (§8). |
| Run reports with versions and full provenance | **PARTIAL** | `PROVENANCE.txt` in the ZIP and an on-screen block exist; they omit most parameters, read tool versions at export time rather than run time, and record no database version or checksum (§5.1, §5.2). No HTML/PDF report. |
| Kraken2 | **NOT IMPLEMENTED** | [CODE] absent. |
| FastTree 2 | **IMPLEMENTED** | §2.12. |
| BOLDigger | **IMPLEMENTED** (BOLDigger3) | §2.11; one small COI run with positive identifications. |
| FQC | **NOT IMPLEMENTED** | [CODE] absent; FastQC results are shown in PARAFISH's own QC page instead. |

---

## 10. Numbers that can be cited

Each number is followed by its exact source and the conditions. **Conditions common to all runs
below:** host machine WSL2, 3.0 GB RAM, 12 cores ([CODE] `vsearch … 2.8GB RAM, 12 cores` in
`step3_cluster.log`; preflight log lines); host tool versions of §2.0, not Docker.

### 10.1 18S V4 reference run — `runs/_real/current`

Conditions: 2 samples (APTL64B, APTL64F), paired-end, TAReuk primers removed with cutadapt,
Trimmomatic → PEAR → VSEARCH 2.30.0 at 0.97, `min_unique_size` 1, headless defaults
(SLIDINGWINDOW 4:20), SILVA + PR2 (marker-filtered) + NCBI remote megablast, hybrid phylo with 10
bootstraps. No `run_config` block (predates it). [CODE] `preprocess/metrics/read_tracking.csv`,
`analysis/metrics/otu/*.csv`, `databases/otu/*.csv`.

| Quantity | aptl64b_02 | aptl64f_01 | Total |
|---|---|---|---|
| Input read pairs | 101 190 | 90 220 | 191 410 |
| Reads with primer | 101 178 (99.99 %) | 90 218 (100.0 %) | |
| After Trimmomatic (both mates kept) | 84 416 (83.42 %) | 78 040 (86.50 %) | 162 456 |
| Merged by PEAR | 84 411 (99.99 % of trimmed) | 78 035 (99.99 %) | 162 446 |
| After VSEARCH filter | 79 830 (94.57 % of merged) | 74 393 (95.33 %) | 154 223 |
| Mapped to OTUs | 79 830 | 74 393 | **154 223** |
| Unmapped | 0 | 0 | 0 |
| Mapped as % of input | 78.89 % | 82.46 % | 80.57 % (computed) |
| Observed OTUs (rarefied to 74 393) | 9 | 2 | 9 OTUs in the run |
| Shannon | 0.0083 | 0.0003 | |

Other figures: 9 raw centroids, **0 chimeras**; abundance multiset
`[2, 2, 2, 2, 4, 4, 12, 43, 154 152]` ([DOC] `DEVLOG.md:2487`); OTU1 = 99.95 % of reads
([DOC] `DEVLOG.md:230`); NCBI returned 1 046 hits; all 9 OTUs have `Consensus_Sources =
NCBI|PR2|SILVA` and `Best_Source = NCBI`; best identities 95.92–100 %; five classes resolved
(Actinopteri, Eurotiomycetes, Dothideomycetes, Insecta, Mammalia). Consensus: OTU1, OTU2 *Gadus*
(genus only), OTU6 *Melanogrammus*, OTU3 *Aspergillus*, OTU4 *Penicillium coffeae*, OTU5
*Pseudosplanchnonema phorcioides*, OTU7 *Saturnia pavonia*, OTU8 Gliridae, OTU9 *Nothopassalora
personata* [CODE] `taxonomy_consensus.csv`.

### 10.2 Fixture runs

| Quantity | `_fixture` (OTU) | `_fixture_asv` (DADA2) | Source |
|---|---|---|---|
| Input | 2 × 5 000 pairs | 2 × 5 000 pairs | [CODE] `read_tracking.csv` |
| After Trimmomatic | 4 182 / 4 349 | 4 182 / 4 349 | [CODE] |
| Final reads | 3 969 / 4 138 mapped (total 8 107) | 4 145 / 4 313 non-chimeric (total 8 458) | [CODE]; `verify.py` |
| Features | 4 OTUs | 3 ASVs | [CODE] |
| DADA2 tracking (b_02) | — | input 4 182 → filtered 4 148 → denoised F 4 145 / R 4 148 → merged 4 145 → non-chimeric 4 145 | [CODE] `analysis/dada2/read_tracking.csv` |
| Taxonomy | SILVA 4/4, PR2 4/4, NCBI 4/4 (649 hits) | SILVA 3/3, PR2 3/3 | [CODE] `golden/*/taxonomy/…/source_status.csv` |

### 10.3 COI case — one sample, three runs

Sample `aptl64b_01`, 207 528 input read pairs. [CODE] each run's `read_tracking.csv`,
`analysis/pooled/*.fasta`, `analysis/metrics/otu/*.csv`.

| Quantity | `…0723T072541Z` (no primer removal, floor 2) | `…0724T082931Z` (primers removed, floor 2) | `…0803T111144Z` (no primer removal, floor 1) |
|---|---|---|---|
| Reads with forward primer | not measured | 191 581 (92.32 %) | not measured |
| After Trimmomatic | 195 772 (94.34 %) | 192 981 (92.99 %) | 195 772 |
| Merged (PEAR) | 182 551 (93.25 %) | 181 715 (94.16 %) | 182 551 |
| After VSEARCH filter | 179 061 (98.09 %) | 179 054 (98.54 %) | 179 061 |
| Unique sequences after dereplication | 7 120 | 1 653 | 52 485 |
| Raw centroids | 99 | 125 | 1 659 |
| Chimeric centroids (reads) | 8 (10 reads, 0.006 %) | 13 (18 reads, 0.010 %) | 180 (114 reads, 0.064 %) |
| **OTUs** | **90** | **112** | **1 464** |
| Mapped reads | 174 900 | 176 492 | 178 578 |
| Unmapped reads | 4 161 | 2 562 | 483 |
| Mapped as % of input | 84.28 % | 85.04 % | 86.05 % |
| Shannon / Chao1 | 0.8954 / 90.0 | 0.3929 / 112.0 | 1.9143 / 3 402.0 |
| Taxonomy | yes | no | no |

**Taxonomy of the 90-OTU run** [CODE] `runs/projeto_1/20260723T072541Z/databases/otu/taxonomy_consensus.csv`,
`source_status.csv`, `analysis/otu_matrix.csv` (computed today):

| Quantity | Value |
|---|---|
| OTUs assigned to *Trisopterus luscus* | **61 of 90** |
| Reads in those 61 OTUs | 168 169 of 174 900 mapped (**96.15 %**) |
| Largest two OTUs | OTU1 128 921 reads, OTU2 36 253 reads |
| OTUs with no species | 11 |
| Other species (OTU count) | *Ramularia vizellae* 3, *Goneplax rhomboides* 2, *Micromesistius poutassou* 2, and one each of *Metarhizium* sp., *Penicillium cvjetkovicii*, *Labioleanira yhleni*, *Jassa herdmani*, *Homo sapiens*, Arthropoda sp., *Geotrichum candidum*, … |
| Classes (OTU count) | Actinopteri 63, Malacostraca 5, Polychaeta 4, Dothideomycetes 3, Eurotiomycetes 2, none 7 |
| NCBI | `ok`, 83/90 features, 1 339 BLAST hits |
| BOLD | `ok`, 17/90 features |
| SILVA / PR2 | `skipped` (needed ~2.4 GB / ~1.0 GB; 0.9 GB available) |
| Custom BLAST | `failed` |
| `Best_Source` | NCBI 73, BOLD 11, none 6 |
| `Consensus_Sources` | NCBI only 67, BOLD\|NCBI 16, BOLD only 1, none 6 |
| Reliability | high 70, low 13, moderate 1, unknown 6 |

**Interpretation recorded for 90 → 112** [DOC] `DEVLOG.md:593-622`: OTU1 and OTU2 of the 90-OTU run
differ at 11 positions, all inside the primers (global identity 96.99 %, bodies 100 % identical);
they merged after primer removal (−1). The net +22 is attributed to dereplication: with primer
noise removed, unique sequences fell from 7 120 to 1 653, the largest unique rose from 284 to
139 938 reads, and singletons discarded by the floor of 2 fell from 45 365 (25.3 % of filtered
reads) to 8 087 (4.5 %), so low-abundance variants reached the threshold to seed OTUs. Modal
amplicon length band 90.3 % → 96.40 %; 365 − 313 = 52 bp = the two primers. An estimate of host
features in the 112-OTU run was **withdrawn** in the same entry.

**Caveats for citing the COI case:**

- The 61/90 figure comes from the run **without** primer removal and with a singleton floor of 2 —
  neither is the platform's current default.
- The 1 464-OTU run shows that with singletons kept and primers **not** removed, the OTU count
  inflates sixteen-fold on this sample. [INFER] it is the primer-region noise (52 485 uniques) that
  drives this; it contradicts any general statement that floors 1 and 2 give the same output.
  **CONTRADICTION:** [DOC] `core/config.py` comment on `MIN_UNIQUE_SIZE` ("1 and 2 produce IDENTICAL
  output — same 9 OTUs…") is true only of the 18S reference run with primers removed, and
  [DOC] `DEVLOG.md:2483-2489` already corrects "identical" to "same sequences and abundances, ids
  permuted".
- [UNKNOWN] the OTU count and taxonomy for this sample under the current defaults (primers removed,
  floor 1).
- One sample only; no replicate, no negative control.

### 10.4 COI + BOLD demonstration — `runs/_coi/current`

2 × 2 000 read pairs, no primer removal, floor 1, BOLD only [CODE]: trimmed 1 578 / 1 672; merged
1 568 / 1 671; filtered 1 566 / 1 669; mapped 1 559 / 1 667 (total 3 226); **55 OTUs**; BOLD assigned
**20 of 55**; `verify.py` (f): "BOLD · 20/55 features with Phylum..Family".

### 10.5 Run times

| Measurement | Value | Conditions | Source |
|---|---|---|---|
| Whole pipeline, COI 90-OTU run | **18 min 14 s** | 1 sample, 207 528 pairs, 4 sources requested, hybrid phylo | [CODE] `progress.json` log timestamps 08:25:41 → 08:43:55 |
| — QC | 14 s | | [CODE] same |
| — Preprocess (Trimmomatic + PEAR) | 51 s | | [CODE] |
| — Analysis | ~23 s | | [CODE] |
| — Databases (NCBI remote megablast + BOLD, 90 OTUs) | **15 min 01 s** | | [CODE] |
| — Phylo (90 tips, MAFFT + FastTree + RAxML-NG) | 1 min 45 s | bootstrap count not recorded in the run | [CODE] |
| Databases, BOLD only, 55 OTUs | 2 min 15 s | `_coi` | [CODE] `progress.json` 16:23:50 → 16:26:05 |
| Analysis, 1 464 OTUs (floor 1, no primer removal) | 41 s | | [CODE] `…0803T111144Z` 12:12:47 → 12:13:28 |
| Remote megablast, 9 OTUs | ~3 min, 123 hits | 18S | [DOC] `DEVLOG.md:2242` |
| VSEARCH clustering, floor 2 vs floor 1 | 57 s vs 247 s | 18S real run, 1 475 vs 4 232 uniques | [DOC] `DEVLOG.md:2473-2474` |
| SILVA + PR2 SINTAX, filtered references, 12 OTUs | 7.8 s end to end | run `pr/20260709T083137Z`, 300 MB free | [DOC] `DEVLOG.md:2921-2923` |
| Building a filtered reference | ~4 s ([DOC] `DEVLOG.md:2924`) / "~2 minute build" ([DOC] `reference_filter.py:34`) | **the two figures disagree** | both [DOC] |
| Test suite | 64–116 s | 581 tests, 2026-08-03 | [DOC] session record / DEVLOG |

### 10.6 Memory

| Measurement | Value | Source |
|---|---|---|
| Machine | 3.0 GB total (2.8 GB reported by VSEARCH), swap 1 GB | [CODE] log lines in `progress.json`; `step3_cluster.log` |
| SINTAX peak RSS, filtered SILVA + PR2 | 258 MB | [DOC] `DEVLOG.md:2923` |
| Reference-filter build | 29 MB RSS | [DOC] `DEVLOG.md:2924` |
| SINTAX memory rule | 3.0 × reference FASTA size | [CODE] `core/preflight.py:24` (a rule, not a measurement) |
| Base pipeline / per sample / DADA2 baseline | 600 MB / 120 MB / 850 MB | [CODE] `core/preflight.py:28,39,35` (rules) |
| `assignTaxonomy` rule | 6 × reference size | [CODE] `dada2_runner.py:37` |
| Sample cap | 96 (`ABSOLUTE_MAX_SAMPLES`) | [CODE] `core/preflight.py:44` |

[UNKNOWN] measured peak RAM per stage for any complete run; no run records it.

### 10.7 Code and tests

~26 000 lines of Python in 102 non-test modules ([CODE] `wc -l` = 25 990 over `qc/ preprocess/
analysis/ databases/ phylo/ core/ ui/ scripts/ app.py`); 53 test modules, 583 test functions
([CODE] AST count); 107 run directories.

---

## 11. Open issues

There are **no `TODO`, `FIXME`, `XXX` or `HACK` markers** in the Python code ([CODE] grep). Open
items are recorded in the DEVLOG and found by reading the code. Ordered by how likely a jury is to
raise them.

### 11.1 Scientific

1. **No validation against ground truth.** No mock community, no sensitivity or false-positive
   estimate, no negative control in any run. [CODE] §9.
2. **All results rest on two specimens.** Every completed taxonomy is on APTL64B and/or APTL64F
   (one or two samples). No biological replicate; ordination, PERMANOVA and contaminant detection
   have never run on real data. [CODE] §7.
3. **Sample provenance is undocumented** in the repository. [UNKNOWN] §7.2.
4. **The COI headline number (61/90) predates primer removal and uses a non-default singleton
   floor**; the run under current defaults does not exist, and the run with primers removed has no
   taxonomy. [CODE] §10.3.
5. **Singleton policy on COI.** Floor 1 without primer removal gives 1 464 OTUs on one sample.
   [CODE] §10.3.
6. **The consensus is not a strict LCA**, applies no identity threshold, and can output a lineage
   whose genus does not match its species (OTU9). [CODE] §3.3.
7. **The Linnaean backbone depends on a live remote service** (NCBI, or BOLD for COI). SILVA and PR2
   supply only domain and genus/species. [CODE] §3.5; [DOC] `DEVLOG.md:2283-2290`.
8. **Remote megablast can miss hits near 90 % identity**; the sensitive search needs a local `nt`
   (~250 GB per [DOC] `.env:6`). [DOC] §8.
9. **Remote sources are unversioned**: NCBI `nt`, BOLD and NCBI Taxonomy state are not recorded per
   run. [CODE] §5.1.
10. **Marker resolution inconsistency** in the full-pipeline path (file-name marker, not preset).
    [CODE] §3.2.
11. **COI preset reverse primer (jgHCO2198) is not the laboratory's (LoboR1).** [CODE] §2.4.
12. **No 16S preset, and no `16S` marker value**, although 16S data exist in the supervisor's folder.
    [CODE] §3.2, §7.3.
13. **Fixed MINLEN and `min_overlap` 10** — a short mate can merge on ~11 bp; raising `min_overlap`
    is open. [DOC] `DEVLOG.md:525-536`.
14. **Chimera removal after clustering, de novo only.** [CODE] §2.8.
15. **DADA2 is fed Trimmomatic-trimmed reads; no seed.** [CODE] §2.2, §2.9.
16. **PERMANOVA unseeded; Aitchison PERMANOVA computed on the 2-D projection.** [CODE] §3.7, §4.2.
17. **Fixed substitution model (GTR+G), no model selection, no outgroup; 10 or 100 bootstraps
    depending on entry point.** [CODE] §4.1.
18. **No strand handling in clustering, mapping or SINTAX** (`--strand both` absent); MAFFT
    `--adjustdirection` is the only orientation step. [CODE] §2.8, §2.11.
19. **Class-level information from PR2 is discarded** (a loss for protist parasites). [DOC] §3.5.

### 11.2 Engineering and reproducibility

20. **Docker image never built or run.** [CODE] §1.3.
21. **Versions that produced the results differ from the pinned ones** (VSEARCH 2.30.0 vs 2.28.1,
    FastQC 0.11.5 vs 0.12.1, BLAST 2.17.0 vs 2.15.0, R/DADA2 4.5.0/1.36.0 vs 4.3.3/1.30.0, MAFFT
    7.453 vs 7.526). [CODE] §2.0.
22. **The development environment is currently broken** (virtualenv interpreter missing; `java` not
    on PATH), so neither the tests nor a pipeline run can be executed today. [CODE] §5.5.
23. **Per-run provenance is incomplete**: most parameters, tool versions, database versions and the
    code commit are not written to the run directory. [CODE] §5.1.
24. **UI and headless defaults disagree** (sliding window 15 vs 20; bootstraps 100 vs 10; BOLD mode
    1 vs 2; NCBI max hits 5 vs fallback 20). [CODE] §2.5, §2.11, §2.12.
25. **Dead controls** in Setup (MAFFT mode, adjust direction, VSEARCH-panel SINTAX cutoff, species
    threshold, workers) and FastTree flags that cannot be turned off. [CODE] §2.12, §3.3, §6.
26. **Empty UI sections** on Results; empty pipeline bars on Preprocess and Analysis. [CODE] §6.
27. **Golden taxonomy layer is stale** and `golden/real` is not in the repository; `golden.py check`
    does not re-run the pipeline. [CODE] §5.3.
28. **A fresh clone has no data**: no fixture, no runs, no references. [CODE] §1.2.
29. **Single-process, single-thread execution in the web server** with `debug=True`; one run at a
    time; a server restart kills the run. [CODE] §1.4.
30. **Sequential sample processing**; the "Workers" setting does nothing. [CODE] §1.4.
31. **Memory**: full SILVA/PR2 do not fit on the development machine; a run stopped in Databases
    twice. [CODE] §3.8, §7.2.
32. **No demultiplexing, no Kraken2, no UNOISE, no FQC, no phylogenetic placement, no curated
    fish/parasite library.** [CODE] §9.
33. **Residual Portuguese strings** in user-visible logs and in the `source` column of
    `source_status.csv`. [CODE] §5.5.
34. **RAxML support-tree key points at a non-existent file name.** [CODE] §2.12.
35. **ASV branch**: phylo and ordination never exercised (the ASV fixture has 3 features; the tree
    floor is 4). [DOC] `DEVLOG.md:963-1025`; [CODE] `runs/_fixture_asv`.
36. **Unlicensed**: no `LICENSE` file. [CODE]

---

## 12. Contradictions and unknowns index

### 12.1 Contradictions (code or data versus documentation)

| # | Topic | Documentation says | Code / data shows | Section |
|---|---|---|---|---|
| C1 | Run marker | [DOC] DEVLOG Task 1: derived from the primer preset | [CODE] full-pipeline thread uses file-name marker only (`callbacks_progress.py:669`) | §3.2 |
| C2 | Multi-source consensus | [DOC] DEVLOG 2026-07-24: "has NEVER run" | [CODE] `runs/_real/current`: SILVA, PR2 and NCBI all `ok`, 9/9, with filtered references | §3.8 |
| C3 | `ENABLE_SILVA/PR2/BOLD/NCBI` | [DOC] README "Which sources run" | [CODE] the flags are read only when choosing a DADA2 reference (`ui/runtime_state.py:304-308`); `enable_bold` and `enable_ncbi` are used nowhere. Sources are chosen by the Setup checklist. | §3.1 |
| C4 | BOLD credentials | [DOC] README: BOLD "requires a registered account … without credentials it is skipped" | [CODE] `identify` does not log in; credentials are optional (`bold_db.py:556-560`); only `download_db` needs them | §2.11 |
| C5 | Golden datasets | [DOC] README: three datasets frozen | [CODE] `golden/real` untracked; source runs and fixtures untracked | §5.3 |
| C6 | Tool versions | [DOC] README lists the Docker pins | [CODE] stored runs used different versions | §2.0 |
| C7 | Trimmomatic sliding window | [DOC] UI tooltip "Default: 4:15" | [CODE] `RunConfig` default 4:20; headless runs use 20 | §2.5 |
| C8 | Trimmomatic MINLEN | [DOC] comment "50, matching the Setup field" | [CODE] the runner's own fallback is still 200 (`trimmomatic_runner.py:105,122`) | §2.5 |
| C9 | RAxML bootstraps | [DOC] UI default 100 | [CODE] `RunConfig` 10; golden run used 10 | §2.12 |
| C10 | BOLD mode | UI default 1 | [CODE] `RunConfig`/`BOLDConfig` default 2 | §2.11 |
| C11 | COI primers | [CODE] preset label "mlCOIintF / jgHCO2198" | [CODE] laboratory list: mlCOIintF / LoboR1 | §2.4 |
| C12 | Singleton floor | [DOC] `core/config.py` comment: floors 1 and 2 give identical output | [CODE] true only for the 18S run; COI without primer removal: 90 vs 1 464 OTUs (different floors) | §10.3 |
| C13 | DADA2 determinism | [DOC] DEVLOG: deterministic | [CODE] no seed, `randomize=TRUE` | §2.9 |
| C14 | RAxML support file | [CODE] runner returns `raxml.raxml.supportFBP` | [CODE] file on disk is `raxml.raxml.support` | §2.12 |
| C15 | FastTree options | [CODE] `RunConfig` has `fasttree_gtr`/`fasttree_gamma` | [CODE] runner reads `gtr`/`gamma` | §2.12 |
| C16 | SINTAX cutoff control | [CODE] `tests/test_setup_defaults_single_source.py` maps `setup-vsearch-sintax` | [CODE] the effective control is `setup-sintax-cutoff` | §3.3 |
| C17 | `best_sintax_record` docstring | [DOC] says confidence falls back to bitscore | [CODE] it does not | §3.6 |
| C18 | PR2 version label | [CODE] `PR2Config.version = "5.0"` | [CODE] file used is 5.1.1 | §3.1 |
| C19 | Reference-filter build time | [DOC] "~4 s" (DEVLOG) | [DOC] "~2 minute build" (module docstring) | §10.5 |
| C20 | Golden taxonomy files | [CODE] code writes `feature_id`, rounded identities, `sintax_status` | [CODE] frozen files have `ASV_ID`, raw decimals, no such column (acknowledged as stale) | §5.3 |
| C21 | NCBI max hits | UI default 5 | [CODE] full-pipeline fallback 20 | §2.11 |
| C22 | "Setup fields agree with defaults" | [DOC] DEVLOG: `KNOWN_DIVERGENCES` empty, class closed | [CODE] the test covers 11 numeric fields; C7, C9, C10 are outside it | §2.5 |

**Update 2026-10-12 — branch `thesis-fixes`, not merged into `main`.** Status of the contradictions above after the four fixes (commits
`64c8446`, `2a88619`, `8be350f`, `1cf052f`). Everything not listed is unchanged.

| # | Status on `thesis-fixes` | Commit |
|---|---|---|
| C1 | **Resolved**: the full-pipeline thread uses the preset-first resolution (§3.2) | `64c8446` |
| C7 | **Resolved**: Setup reads `RunConfig` → 4:20; tooltip corrected (§2.5) | `8be350f` |
| C9 | **Resolved**: Setup reads `RunConfig` → 10 bootstraps (§2.12) | `8be350f` |
| C10 | **Resolved**: Setup reads `RunConfig` → BOLD mode 2 (§2.11) | `8be350f` |
| C11 | **Addressed, not removed**: a new preset `COI_Leray_LoboR1` carries the laboratory's pair; `COI_Leray` keeps its jgHCO2198 label and sequences (§2.4) | `2a88619` |
| C21 | **Open, by decision**: `RunConfig` has no NCBI max-hits default to align to (§2.11) | — |
| C22 | **Narrowed**: the single-source test now also covers the sliding-window text field, the bootstrap field and the BOLD dropdown; NCBI max hits stays outside it | `8be350f` |
| C8, C14, C15, C16 | Unchanged (out of scope of these fixes) | — |

Related open items of §11 that these fixes change: item 10 (marker resolution) resolved; item 11
(COI preset) addressed; item 23 (per-run provenance) largely addressed for new runs, with the
exceptions listed in §5.1; item 24 (UI vs headless defaults) resolved except NCBI max hits.

Verification on the branch [CODE] 2026-10-12, with the project virtualenv active: `pytest` 621 passed,
2 skipped (baseline 581 / 2; 40 new tests); `verify.py` and `golden.py check all` identical to the
baseline for all four runs and all three datasets. Note that U1 is thereby answered — the suite
runs again; the virtualenv was repaired on 2026-10-12 and now carries Dash 4.4.1, pandas 3.0.6 and
scikit-bio 0.7.4, newer than the versions listed in §1.1.

### 12.2 Unknowns

| # | Unknown |
|---|---|
| U1 | Whether the test suite passes today (environment broken; last recorded 581 pass / 2 skip on 2026-08-03) |
| U2 | Whether the Docker image builds and runs |
| U3 | Origin and biology of samples APTL64B / APTL64F (host, tissue, locality, project, provider) |
| U4 | Whether the `_coi` demo inputs are true COI amplicons of the same specimens |
| U5 | Whether the supervisor's Novogene datasets (11 samples × 18S, COI, 16S) have been analysed at all |
| U6 | Any trace of the Guadiana dataset (PRJNA1334324) — none in the repository |
| U7 | COI feature count and taxonomy under current defaults (primers removed, floor 1) |
| U8 | Taxonomy of the 112-OTU COI run |
| U9 | Measured peak RAM per stage |
| U10 | Version/date of NCBI `nt`, BOLD and NCBI Taxonomy at the time of each run |
| U11 | Checksums of the SILVA and PR2 files used |
| U12 | MAFFT, Trimmomatic and Java versions at the time of the stored runs (not logged) |
| U13 | Why Entrez was chosen over a local taxdump (no comparison recorded) |
| U14 | Data-based justification for 0.97 and 0.8 (none recorded) |
| U15 | Exact wording of the pre-thesis commitments |
| U16 | Whether DADA2 output is reproducible beyond fixture scale |
| U17 | Whether a single-end run has ever completed end to end |
| U18 | Current visual/interactive state of the UI (not opened in a browser for this document) |
| U19 | Content and tracking status of `golden_baseline/` and `data/demo_custom_db_*.fasta` |
| U20 | Bootstrap count used in the COI run's tree |

---

## Summary for the thesis

PARAFISH is a browser-based platform (Python, Plotly Dash) that takes per-sample Illumina amplicon
FASTQ files from raw reads to taxonomy and phylogeny in five automated stages: quality control
(FastQC); primer removal, quality trimming and merging (cutadapt, Trimmomatic, PEAR); feature
inference as 97 % OTUs (VSEARCH, de novo chimera removal) and/or ASVs (DADA2); taxonomic assignment;
and tree building (MAFFT, FastTree, RAxML-NG). Taxonomy combines SILVA and PR2 (SINTAX), NCBI `nt`
(remote megablast, lineages via Entrez), BOLD (BOLDigger3) and an optional user reference into one
call per feature, keeping alignment identity and k-mer confidence in separate columns and stating
per source whether it ran, found nothing, failed or did not apply to the marker. Large references
are cut to the amplified region before loading, which allows local classification on a 3 GB
machine. Downstream it computes rarefied alpha diversity and CLR/Aitchison or distance-based
ordinations with PERMANOVA, and exports abundance tables, BIOM and a phyloseq bundle. Seeded
classification and tree search, test-pinned parameters, frozen per-stage outputs, consistency checks
and a provenance file per run archive support reproducibility. It has so far been run end to end on
two specimens (18S: 9 OTUs; COI: 61 of 90 OTUs, 96 % of reads, assigned to the host *Trisopterus
luscus*); it has no demultiplexing, Kraken2, UNOISE, phylogenetic placement, curated fish-parasite
library or ground-truth benchmark, and its container image and multi-sample statistics are untested.
