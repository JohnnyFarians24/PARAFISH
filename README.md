# PARAFISH

**PARAFISH** is a containerised web platform for analysing fish-parasite communities from
amplicon sequencing (metabarcoding) data. It takes raw Illumina FASTQ files through quality
control, primer removal, merging, OTU/ASV inference, taxonomic classification and
phylogenetic inference, and presents the whole run in an interactive Plotly Dash interface
in the browser.

What distinguishes it from a script collection is that a run is a **contract**, not a
by-product: several reference databases are reconciled into one consensus per feature, OTUs
and ASVs can be produced and compared side by side in the same run, and every deliverable
travels with the recipe that made it. Scientific parameters are pinned by tests, the
randomised steps are seeded, and each stage's output is frozen against a golden baseline, so
a number that moves is a regression rather than a surprise.

---

## Pipeline overview

Ten pages in the UI, seven of them pipeline stages. Each stage is independent — you can
enter at any point where the upstream results already exist on disk.

| # | Stage | Tools | What it produces |
|---|---|---|---|
| 1 | **Project** | — | A named project holding samples and runs |
| 2 | **Samples** | — | Paired-end detection from `_R1`/`_R2`, sample groups, control flags |
| 3 | **Set up your run** | — | Marker, primers, method, database selection, and a memory preflight |
| 4 | **Quality Control** | FastQC | Per-sample QC reports plus `qc_overview`, per-base quality, length distribution, module status |
| 5 | **Preprocess** | cutadapt → Trimmomatic → PEAR | Primer-trimmed, quality-trimmed, merged reads and `read_tracking.csv` |
| 6 | **Analysis** | VSEARCH (OTU) and/or DADA2 (ASV) | `otu_matrix.csv`, representative FASTA, chimera stats, alpha diversity, rarefaction, occupancy |
| 7 | **Databases** | VSEARCH SINTAX, BLAST, BOLDigger3 | Per-source assignments reconciled into one consensus; hits, conflicts, per-rank taxa tables, BIOM, phyloseq |
| 8 | **Phylo** | MAFFT → FastTree and/or RAxML-NG | Alignment and tree (Newick, with branch support) |
| 9 | **Progress** | — | Live per-stage log while the run executes |
| 10 | **Results** | — | Run summary, read funnel, ordination, downloads, reproducibility block |

The execution order the pipeline actually runs is
`qc → preprocess → analysis → databases → phylo` (see `scripts/run_pipeline.py`); Project,
Samples and Set up are configuration, and Progress and Results are views.

**Ordination** (Aitchison PCA by default, plus PCoA on Bray-Curtis, Jaccard and weighted /
unweighted UniFrac, with PERMANOVA) is available on the Analysis and Results pages once a
run has at least 3 samples.

---

## Requirements

- **Docker** 24.0+ and **Docker Compose** 2.24+ (`docker compose`, not `docker-compose`) —
  check with `docker --version && docker compose version`
- **~5 GB** of disk for the image
- **Disk for reference databases** — they are *not* in this repository. Uncompressed:
  SILVA ~824 MB, PR2 ~334 MB, and a local NCBI `nt` (optional) is far larger.
- **RAM — read this section before blaming the app.**

### RAM, honestly

The memory cost is dominated by one step: VSEARCH SINTAX builds a k-mer index of the
reference and holds it in RAM. Measured multiplier in `core/preflight.py`:
**about 3× the reference FASTA size.**

| Configuration | Approximate peak RAM |
|---|---|
| Base pipeline (QC → Preprocess → Analysis), per run | ~600 MB + ~120 MB per sample |
| DADA2 (ASV branch) | +850 MB just to load the R libraries, regardless of sample size |
| SINTAX against **full** SILVA (~824 MB) | ~2.5 GB |
| SINTAX against **full** PR2 (~334 MB) | ~1.0 GB |
| SINTAX against a **marker-filtered** reference | ~75–230 MB |

**Minimum: 4 GB**, and only with marker filtering on (see below). **Recommended: 8–16 GB**,
which is what a multi-source run against full references or a large sample set needs.

PARAFISH mitigates this rather than just failing: when the run declares a marker and
primers, each reference is cut down to the amplified region before it is ever loaded.
Measured on a real 18S V4 run — SILVA **788 MB → 25 MB** (54 207 of 510 495 sequences) and
PR2 **343 MB → 76 MB** — which is what makes local classification possible on a small
machine, and also makes SINTAX *more* accurate, since a 380 bp query is then scored against
380 bp references rather than 1 800 bp whole-SSU sequences.

The Set up page runs a preflight against your actual free memory and says, before you start,
which sources will be filtered, which will be skipped, and how many samples are safe. A
source that will not fit is skipped with a stated reason rather than being killed mid-run by
the OS.

---

## Installation

```bash
git clone https://github.com/JohnnyFarians24/PARAFISH.git
cd PARAFISH
cp .env.example .env
# edit .env — at minimum DB_SILVA_PATH, DB_PR2_PATH and RUNS_DIR (see Configuration)
docker compose up --build
```

Then open <http://localhost:8050>.

The image is built from `continuumio/miniconda3` and installs every bioinformatics tool
itself (`environment.yml` for FastQC 0.12.1, Trimmomatic 0.39, MAFFT 7.526, RAxML-NG 1.2.2,
BLAST+ 2.15.0, FastTree 2.1.11, PEAR 0.9.11, R 4.3.3 + DADA2 1.30.0; VSEARCH 2.28.1 is
fetched as a static binary). Nothing beyond Docker is needed on the host.

> **Not yet verified on a clean machine.** The `git clone` → `docker compose up --build`
> path has *not* been executed end to end from a fresh checkout, because Docker is not
> available in the development environment this was last worked in. The commands above are
> transcribed from the real `docker-compose.yml` and `Dockerfile`, not from memory, but
> treat a first build as unproven and see [Known limitations](#known-limitations).

To stop, restart, or rebuild:

```bash
docker compose down          # stop
docker compose up            # restart without rebuilding
docker compose up --build    # rebuild after code changes
```

Run state lives in `RUNS_DIR` on the host, so restarting the container does not lose a run.

`docs/setup_docker.md` carries the longer setup notes (directory preparation, volume
layout, troubleshooting).

---

## Reference databases

The databases are **not distributed with this repository**. Download them, place them under
your `DB_ROOT`, and point `.env` at them. `docker-compose.yml` mounts
`${DB_ROOT:-./data/dbs}/silva`, `/pr2`, `/bold_db` and `/blast_db` read-only into `/data/…`
inside the container, which is what the default paths in `.env.example` assume.

### SILVA 138.2 SSU NR99 — local, fast, no network

```bash
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
# optional, recommended: SILVA's rank file, so SILVA also gives phylum..genus
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/taxonomy/tax_slv_ssu_138.2.txt.gz
# place both under $DB_ROOT/silva/
```

The SILVA FASTA lineages do not say which level is a phylum, a class or a genus. Without
the rank file PARAFISH therefore takes only the domain and the organism name from SILVA.
With `DB_SILVA_RANKS_PATH` pointing at `tax_slv_ssu_<release>.txt[.gz]`, each level gets
the rank SILVA itself declares for it (domain, phylum, class, order, family, genus; other
declared ranks such as `major_clade` are left out). The ranked reference is cached under
its own name (`*.sintax.v4r.<hash>.fasta`), so the unranked one is never overwritten. Note
that SILVA's ranks are its own and can differ from NCBI's (for fish, *Actinopterygii* in
SILVA against *Actinopteri* in NCBI), so a run that combines SILVA with NCBI or BOLD can
list rank conflicts that are naming, not biology.

### PR2 v5.1.1 — local, protist-focused

```bash
wget https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.fasta.gz
wget https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.tax.gz
# place under $DB_ROOT/pr2/
```

Both are filtered to your marker on first use and the filtered copy is cached next to the
source, keyed by the source's mtime and your primers — so the reduction cost (a few minutes)
is paid once.

### NCBI and BOLD — remote, and slow

These are **queried over the network, one feature at a time**, and that is by far the
slowest part of a run: **minutes to well over an hour**, depending on how many OTUs/ASVs you
have. The Databases stage can look frozen when it is only waiting on NCBI's queue. This is
not a bug, and the preflight warns about it before you start.

- **NCBI** is disabled by default (`ENABLE_NCBI=False`). It uses the public remote `nt`
  service. Set `ENTREZ_EMAIL` in `.env` so NCBI will resolve full Kingdom→Family lineages;
  without it, only genus and species are filled. `NCBI_API_KEY` raises the rate limit from
  3 to 10 requests/s. A local BLAST database is available as an emergency exit
  (`PARAFISH_NCBI_LOCAL_DB=1` + `DB_BLAST_PATH`).
- **BOLD** requires a registered [BOLD Systems](https://boldsystems.org/) account
  (`BOLD_USERNAME` / `BOLD_PASSWORD`). Without credentials it is skipped cleanly and the
  rest of the run continues. BOLD is COI-only, so on an 18S run it is correctly marked
  *not applicable* rather than reported as having found nothing.

If you choose `both` (OTU **and** ASV) as the analysis method, the remote sources run
**twice**, once per feature set. The preflight says so.

### Your own reference (custom database)

In the Databases stage you can upload your own FASTA and compare every OTU/ASV against it:

- a **plain reference FASTA** becomes a local BLAST database (`makeblastdb`) and the run's
  representative sequences are aligned to it with `blastn`, keeping the best hit per feature;
- a **SINTAX-formatted FASTA** (headers carrying `;tax=…`) is folded into the taxonomy
  consensus alongside SILVA and PR2.

Uploads are validated as DNA before use — an empty body, a protein FASTA, all-gap records,
or sequences below 50 bp are rejected with a specific reason rather than silently producing
zero hits. The comparison yields a per-feature table (unmatched features flagged `None` in
red), a taxon × sample detection heatmap with an explicit `No match` row, and a
matched/unmatched donut, with filters for minimum percent identity and species.

---

## Configuration (`.env`)

Copy `cp .env.example .env` and edit. Every variable below is read by `core/config.py`.

### Database paths

| Variable | Meaning | Example |
|---|---|---|
| `DB_SILVA_PATH` | SILVA FASTA (`.fasta` or `.fasta.gz`) | `/data/silva/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz` |
| `DB_SILVA_RANKS_PATH` | Optional. SILVA rank file; when set, SILVA contributes phylum..genus as SILVA declares them | `/data/silva/tax_slv_ssu_138.2.txt.gz` |
| `DB_PR2_PATH` | PR2 mothur FASTA | `/data/pr2/pr2_version_5.1.1_SSU_mothur.fasta.gz` |
| `DB_PR2_TAX_PATH` | PR2 mothur taxonomy file | `/data/pr2/pr2_version_5.1.1_SSU_mothur.tax` |
| `DB_PR2_DADA2_PATH` | PR2 in DADA2 format, for `assignTaxonomy()` | `/data/pr2/pr2_..._dada2.fasta` |
| `DB_BOLD_PATH` | BOLD local database directory | `/data/bold_db` |
| `DB_BLAST_PATH` | Local BLAST `nt` prefix | `/data/blast_db/nt` |
| `GLOBAL_DB_ROOT` / `DB_ROOT` | Root under which all of the above live | `/data/dbs` |
| `BLASTDB` | Directory holding NCBI's `taxdb.btd`/`.bti` — without it `blastn` cannot turn a taxid into a species name | `/data/blast_db` |

### Which sources run

| Variable | Default | Meaning |
|---|---|---|
| `ENABLE_SILVA` | `True` | Local SINTAX against SILVA |
| `ENABLE_PR2` | `True` | Local SINTAX against PR2 |
| `ENABLE_BOLD` | `True` | Remote BOLD (needs credentials, COI only) |
| `ENABLE_NCBI` | `False` | Remote NCBI BLAST (slow) |
| `BOLD_USERNAME`, `BOLD_PASSWORD` | empty | BOLD Systems account |
| `ENTREZ_EMAIL` | empty | Required by NCBI for automated access; also unlocks full lineage resolution |
| `NCBI_API_KEY` / `ENTREZ_API_KEY` | empty | Raises NCBI's rate limit 3 → 10 req/s |

### Scientific parameters

These are deliberate choices, each pinned by a test in
`tests/test_scientific_parameter_pins.py` — changing the value turns the suite red on
purpose.

| Variable | Default | Meaning |
|---|---|---|
| `PARAFISH_MIN_UNIQUE_SIZE` | `1` | `vsearch --minuniquesize`: the smallest cluster allowed to *seed* an OTU. `1` keeps singletons — on a parasite-detection platform a globally unique, >3 %-divergent read is a rare taxon, not an error. `2` buys clustering time on a large host-dominated study, and its cost shows up as `unmapped_reads` in `read_tracking.csv` rather than vanishing. |
| `PARAFISH_DISCARD_UNTRIMMED` | `1` | Drop reads in which cutadapt found no primer. On by default for provenance: a read with no primer match cannot be shown to be an amplicon of the targeted region. |
| `PARAFISH_NCBI_LOCAL_DB` | `0` | Use a local BLAST database instead of the public remote `nt`. Off by default — the public remote is what a reviewer can reproduce. |

Not exposed as environment variables, but pinned the same way: clustering identity `0.97`,
SINTAX cutoff `0.8`, PEAR `min_overlap` `10`, Trimmomatic `MINLEN` `50`.

### Runtime and housekeeping

| Variable | Default | Meaning |
|---|---|---|
| `RUNS_DIR` | `/data/runs` | Where runs are written; this is what persists across container restarts |
| `KEEP_INTERMEDIATES` | `0` | Keep intermediate FASTQ/FASTA after a stage completes |
| `DELETE_AFTER_STEPS` | empty | Comma-separated stages whose intermediates are deleted as soon as they finish |
| `ALLOW_DB_AUTO_DOWNLOAD` | `0` | Allow the app to fetch a missing reference itself |
| `NCBI_REMOTE_FALLBACK` | `1` | Fall back to remote NCBI when the local database is unusable |
| `NCBI_TIMEOUT_SEC` | `1800` | Remote BLAST timeout — raise it for large feature sets |
| `BLAST_TIMEOUT_SEC` | `900` | Local BLAST timeout |
| `MAKEBLASTDB_TIMEOUT_SEC` | `1800` | `makeblastdb` timeout when building a custom database |
| `BOLDIGGER_TIMEOUT_SEC` | `900` | BOLDigger3 timeout |
| `BOLDIGGER_PYTHON` | empty | Python interpreter for BOLDigger3, when it lives in its own virtualenv |
| `BOLD_CACHE_ROOT` | `$DB_ROOT/bold/cache` | BOLD download cache |

### Tool paths

Every tool path can be overridden — `FASTQC_PATH`, `VSEARCH_PATH`, `MAFFT_PATH`,
`FASTTREE_PATH`, `RAXML_PATH`, `PEAR_PATH`, `CUTADAPT_PATH`, `BLASTN_PATH`,
`MAKEBLASTDB_PATH`, `RSCRIPT_PATH`, `TRIMMOMATIC_JAR_PATH`, `TRIMMOMATIC_ADAPTERS`,
`TRIMMOMATIC_ADAPTERS_PE`. Inside the image these resolve automatically; the overrides exist
for running outside Docker. The top bar shows a `Tools n/m` badge naming anything missing.

---

## Input data

Single-end and paired-end Illumina FASTQ. Paired files are detected automatically from
`_R1`/`_R2` in the filename. Accepted: `.fastq`, `.fastq.gz`, `.fq`, `.fq.gz`.

Upload through the web interface, or drop files into `uploads/` before starting the
container.

---

## Running your first analysis

1. **Project** — create a project. Everything else attaches to the active project, which is
   named in the top bar and again on Samples and Set up, so you cannot silently upload into
   the wrong one.
2. **Samples** — upload or select FASTQ files. Pairs are detected automatically; assign
   metadata groups and mark negative controls here if you have them (groups are what
   PERMANOVA and the ordination colouring use later).
3. **Set up your run** — choose the marker and primer preset (or enter primers manually),
   the analysis method (`vsearch` for OTUs, `dada2` for ASVs, `both`), and which reference
   databases to use. Read the **Machine resources** panel: it states, before you commit,
   which sources fit in your free memory and roughly how long the remote ones will take.
4. **Start the run** — the pipeline runs *only* when you press the button. There is no
   auto-start.
5. **Progress** — watch the per-stage log live.
6. **Results** — the run summary: an aggregate status badge, an end-to-end read funnel
   (input → trimmed → merged → clustered → classified, with any single step losing more than
   20 % flagged red), per-stage KPI cards, community composition with a rank selector and a
   host-exclusion toggle, ordination, downloads, and a reproducibility block you can paste
   into a methods section.

**Closing the browser or reloading does not lose the run.** Every page reconstructs its
state from `progress.json` on disk rather than from the browser session, so F5, a new tab,
and a server restart all recover the run as it stands. If a run is interrupted (an
out-of-memory kill, for instance), the app reconciles it on the next startup instead of
presenting the orphan as the active run.

---

## Outputs

Everything lands under `RUNS_DIR/<project>/<run_id>/`, and the **Download all** button on
Results packages it as a ZIP. A real run directory holds ~127 files; these are the ones a
reader actually wants.

| Artefact | Path (feature set = `otu` or `asv`) |
|---|---|
| **`PROVENANCE.txt`** — the recipe: sources, marker, primers, method, branch, `min_unique_size`, NCBI mode, and resolved tool versions | ZIP root (rendered from `progress.json`) |
| Abundance matrix (feature × sample) | `analysis/otu_matrix.csv` |
| Representative sequences | `analysis/pooled/otus.fasta`, `analysis/dada2/asv_sequences.fasta` |
| Chimera statistics, cluster map | `analysis/pooled/`, `analysis/metrics/<set>/chimera_stats.csv` |
| Alpha diversity, rarefaction, occupancy, rank abundance | `analysis/metrics/<set>/` |
| **Taxonomy consensus** — one call per feature, reconciled across sources | `databases/<set>/taxonomy_consensus.csv` |
| **Per-feature summary** — consensus plus the best hit of each evidence family, with read counts, ordered by abundance | `databases/<set>/taxonomy_summary.csv` |
| **All hits**, ranked, with identity/coverage/confidence | `databases/<set>/taxonomy_hits.csv` |
| Per-source assignments, and where sources disagree | `databases/<set>/source_assignments.csv`, `conflict_report.csv` |
| **Per-source status** — which source ran, which found nothing, which crashed, and why | `databases/<set>/source_status.csv` |
| Per-rank taxa tables, absolute and relative (7 ranks) | `databases/<set>/taxa_tables/` |
| BIOM table | `databases/<set>/table.biom` |
| **phyloseq bundle** — `otu_table.tsv`, `tax_table.tsv`, `sample_data.tsv` + `load_phyloseq.R` | `databases/<set>/phyloseq/` |
| Taxonomy-annotated FASTA | `databases/<set>/annotated.fasta` |
| Alignment | `phylo/<set>/mafft/aligned.fasta` |
| Trees (Newick; RAxML's carries branch support) | `phylo/<set>/fasttree/tree.nwk`, `phylo/<set>/raxml/raxml.raxml.support` |
| Read tracking, per sample, per step | `preprocess/metrics/read_tracking.csv` |
| FastQC reports and parsed QC metrics | `qc/<sample>/`, `qc/metrics/` |

The phyloseq bundle is not redundant with BIOM: PARAFISH's BIOM export does not carry sample
groups or control flags, and those are exactly what a downstream comparison needs. Counts
are exported raw — rarefying is the analyst's call.

---

## Headless / scripted use

`scripts/run_pipeline.py` runs the pipeline without the UI, driving the *same* code the app
drives (`run_config_from_store` + `_run_pipeline_thread`), so a scripted run and a clicked
run cannot diverge. It is synchronous: the process exits when the pipeline is done.

```bash
# full pipeline on the bundled fixture inputs
python scripts/run_pipeline.py --out runs/_fixture/current

# full-depth inputs
python scripts/run_pipeline.py --real --out runs/_real/current

# a subset of stages, and a subset of sources
python scripts/run_pipeline.py --stages qc,preprocess,analysis --db-sources silva,pr2
```

| Flag | Default | Meaning |
|---|---|---|
| `--real` | off | Use the full-depth inputs instead of `fixtures/` |
| `--out` | `runs/_fixture/current` (or `runs/_real/current` with `--real`) | Output directory |
| `--stages` | `qc,preprocess,analysis,databases` | Comma-separated subset of `qc,preprocess,analysis,databases,phylo` |
| `--db-sources` | `silva,pr2,bold,ncbi,custom` | Comma-separated subset |
| `--keep` | off | Do not wipe the output directory first |
| `--min-unique-size` | from `PARAFISH_MIN_UNIQUE_SIZE` | Override the OTU seeding floor |

Note that the input FASTQ files are **not** in the repository (see Known limitations), so
this needs data you supply.

---

## Reproducibility

This is the part PARAFISH is built around, and each claim below is enforced by code rather
than by convention.

- **Pinned scientific parameters.** Clustering identity `0.97`, SINTAX cutoff `0.8`, PEAR
  `min_overlap` `10`, Trimmomatic `MINLEN` `50`, `min_unique_size` `1` — each has a test
  asserting its literal value, so changing a number cannot pass silently
  (`tests/test_scientific_parameter_pins.py`).
- **Seeded randomness.** VSEARCH SINTAX runs with `randseed=42` rather than VSEARCH's
  default of 0 ("seed from a random source"), and RAxML-NG with `seed=42`. A classification
  that changed between two runs of identical data would not be a result.
- **Golden contracts per layer.** Each closed stage's output is frozen under `golden/` and
  diffed by `scripts/golden.py` across three datasets (`fixture`, `fixture_asv`, `real`) and
  six layers (qc, preprocess, analysis, taxonomy, export, phylo). Artefacts that embed
  wall-clock metadata rather than content (FastQC's zip banner, BIOM's `date`) are
  deliberately excluded, so a real change cannot hide in that churn.
- **Cross-file consistency checks.** `scripts/verify.py <run_dir>` asserts that files
  written by *different* functions agree — read totals across every table, cluster sizes
  against the abundance matrix, per-rank tables against the taxonomy table, no invented
  middle ranks, and `source_status` against the assignments each source actually made. It
  reports SKIP rather than PASS when an input is missing, so a partial run cannot look like
  a clean bill of health.
- **The recipe travels with the results.** `PROVENANCE.txt` is written into every run ZIP,
  rendered from the run's own `progress.json`, and reports tool versions resolved from the
  *configured* binaries rather than from `$PATH`. When versions cannot be resolved it says
  so explicitly instead of looking complete.

One limit worth stating plainly: `golden.py check` compares a stored run against the frozen
files — it does **not** re-run the pipeline. "No regressions" therefore means the frozen
files are unchanged, not that today's code still produces them.

---

## Running the tests

```bash
python -m pytest tests/ -q
```

**583 tests** (581 passing, 2 skipped pending a golden re-baseline). The suite covers the
stage runners, the taxonomy consensus and its per-source reconciliation, restore-from-disk
on every page, the export and ZIP contract, ordination and phylogenetic diversity, the
scientific parameter pins, and a scan asserting that no Portuguese remains in any
user-facing string.

The golden diffs are separate, and need a stored run to compare against:

```bash
python scripts/golden.py check all --from runs/_fixture/current --dataset fixture
python scripts/verify.py runs/_fixture/current
```

---

## Known limitations

Stated plainly, because knowing where the edges are is more useful than a clean-looking list.

- **The multi-source consensus has not yet run at full depth.** SILVA + PR2 + NCBI live
  simultaneously, reconciled into one consensus per feature, is the platform's central
  feature — and it has so far only been exercised on fixtures or single-source
  configurations, because the development machine (2.8 GB) cannot hold the references. It is
  gated on hardware, not on code.
- **Clean-clone startup is unverified.** `git clone` → `docker compose up --build` has not
  been executed from a fresh checkout on a clean machine. The instructions are transcribed
  from the real compose and Docker files, but the path itself is unproven.
- **No input data ships with the repository.** `fixtures/`, `runs/` and the sample archive
  are excluded by `.gitignore`, so a fresh clone has code and golden baselines but no FASTQ
  files. `scripts/run_pipeline.py` and the golden diffs need data you supply.
- **The ASV branch is less exercised than the OTU branch.** DADA2 inference, taxonomy and
  export are contracted against a frozen ASV dataset, but the ASV *phylo* and *ordination*
  paths are not — the ASV fixture has 3 features and the aligner floor is 4 tips.
- **Remote sources are slow and depend on third parties.** NCBI and BOLD are queried one
  feature at a time; both are disabled or credential-gated by default for that reason.
- **Large references need real memory.** Without marker filtering, SINTAX against full SILVA
  wants ~2.5 GB and will be skipped (with a stated reason) on a machine that cannot supply
  it. See [RAM, honestly](#ram-honestly).
- **BOLD has never returned a real attribution in this project.** Its positive path is
  implemented and unit-tested against the real BOLDigger3 3.x store schema, but every run so
  far has been 18S against a COI-only database, so the path is proven by test rather than by
  a production run.

---

## License

This project is currently **unlicensed** — there is no `LICENSE` file in the repository, so
default copyright applies and no permissions are granted. Contact the author for usage
permissions.

## Citation

There is no publication to cite yet. Until there is, please cite this repository, and cite
the underlying tools according to their own publications — VSEARCH, DADA2, SILVA, PR2, BOLD,
FastQC, Trimmomatic, PEAR, cutadapt, MAFFT, FastTree, RAxML-NG, scikit-bio.

---

Bug reports and feature requests:
[GitHub Issues](https://github.com/JohnnyFarians24/PARAFISH/issues).
