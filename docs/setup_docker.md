# PARAFISH Docker Setup

## Prerequisites

- Docker 24.0 or newer
- Docker Compose 2.24 or newer

Check the installation:

```bash
docker --version && docker compose version
```

## Prepare the data directories

```bash
mkdir -p /data/silva
mkdir -p /data/pr2
mkdir -p /data/bold_db
mkdir -p /data/blast_db
mkdir -p /data/runs
mkdir -p /data/config
```

SILVA and PR2 must already exist on the host and are mounted read-only inside the container. BOLD and BLAST DBs are also mounted as volumes if available.

## Configure environment

```bash
cp .env.example .env
```

Edit `.env` so the container-visible paths match the mounted volumes.

## Build and run

```bash
docker compose build
docker compose up -d
docker compose logs -f parafish
```

The UI is exposed on port 8050.

## Verify tools inside the container

```bash
docker compose exec parafish python -c "import shutil; tools = ['fastqc', 'vsearch', 'pear', 'mafft', 'FastTree', 'raxml-ng', 'blastn']; print({t: bool(shutil.which(t)) for t in tools})"
```

## Updating databases

- Replace the SILVA or PR2 file on the host volume and restart the service.
- Re-download the BOLD DuckDB on the host volume if needed.
- Use the BLAST remote mode if no local NCBI DB is mounted.
