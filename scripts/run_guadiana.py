#!/usr/bin/env python
"""Run PARAFISH headless on the Guadiana COI reads (thesis results).

Driver only: it reuses scripts/run_pipeline.py (build_stores) and the same
run_config_from_store + _run_pipeline_thread the UI uses. No platform code is changed.
"""
import argparse, glob, json, os, sys, time
from pathlib import Path

# A BOLD-sized custom reference needs more than the platform's default timeouts
# (900 s blastn / 1800 s makeblastdb); environment variables are its own override.
os.environ.setdefault("BLAST_TIMEOUT_SEC", "172800")
os.environ.setdefault("MAKEBLASTDB_TIMEOUT_SEC", "172800")

ROOT = Path(__file__).resolve().parent.parent  # scripts/ -> project root
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
import run_pipeline as rp

FWD = "GCCCCCGAYATRGCNTTYCC"
REV = "TANACYTCNGGRTGNCCRAARAAYCA"

p = argparse.ArgumentParser()
p.add_argument("--fastq-dir", default=os.path.expanduser("~/guadiana/fastq"))
p.add_argument("--out", required=True)
p.add_argument("--stages", default="qc,preprocess,analysis,databases,phylo")
p.add_argument("--db-sources", default="bold")
p.add_argument("--method", default="vsearch", choices=["vsearch", "dada2", "both"])
p.add_argument("--threads", type=int, default=16)
p.add_argument("--limit", type=int, default=0, help="only the first N samples (quick test)")
p.add_argument("--fwd", default=FWD, help="forward primer (default: FishATLC2 assay)")
p.add_argument("--rev", default=REV, help="reverse primer (default: jgHCO2198)")
p.add_argument("--no-primers", action="store_true", help="RQ3: do not remove primers")
p.add_argument("--min-unique-size", type=int, default=None)
p.add_argument("--custom-db", default=None, help="reference FASTA for the custom (local BLAST) source")
p.add_argument("--dry-run", action="store_true")
a = p.parse_args()

r1s = sorted(glob.glob(os.path.join(a.fastq_dir, "*_COI_R1.fastq.gz")))
if a.limit:
    r1s = r1s[: a.limit]
inputs = {}
for r1 in r1s:
    run = os.path.basename(r1).replace("_COI_R1.fastq.gz", "")
    inputs[run.lower()] = (run + "_COI", Path(r1), Path(r1.replace("_COI_R1", "_COI_R2")))
if not inputs:
    sys.exit("FAIL no *_COI_R1.fastq.gz files in " + a.fastq_dir)

out_dir = Path(a.out).resolve()
stages = [s for s in a.stages.split(",") if s]
dbs = [s for s in a.db_sources.split(",") if s]
vp = {"threads": a.threads}
if a.no_primers:
    vp["primer_preset"] = ""          # explicit opt-out, obeyed by resolve_primers
else:
    vp["primer_fwd"], vp["primer_rev"] = a.fwd, a.rev
if a.min_unique_size is not None:
    vp["min_unique_size"] = a.min_unique_size
if a.custom_db:
    custom = Path(os.path.expanduser(a.custom_db)).resolve()
    if not custom.exists():
        sys.exit(f"FAIL custom database not found: {custom}")
    vp["custom_db_path"], vp["custom_db_kind"] = str(custom), "blast"

project, samples, run_cfg = rp.build_stores(inputs, out_dir, stages, dbs, out_dir.name, vp)
for s in samples["samples"]:
    s["marker"] = "COI"
project["project_name"] = "guadiana_COI"
project["analysis_method"] = a.method
run_cfg["analysis_method"] = a.method

from ui.runtime_state import run_config_from_store
from ui.callbacks_progress import _run_pipeline_thread, _pipeline_run_marker
from preprocess.cutadapt_runner import resolve_primers

cfg = run_config_from_store(project, samples, run_cfg)
cfg.run_id, cfg.output_dir = out_dir.name, out_dir
print(f"samples    : {len(inputs)}  (first: {r1s[0]})")
print(f"out_dir    : {out_dir}")
print(f"stages     : {stages}   db_sources: {dbs}   method: {cfg.analysis_method}")
print(f"marker     : {_pipeline_run_marker(cfg)}")
print(f"primers    : {resolve_primers(cfg)}")
print(f"custom db  : {vp.get('custom_db_path', '-')}")
print("-" * 60, flush=True)
if a.dry_run:
    sys.exit(0)

if out_dir.exists() and any(out_dir.iterdir()):
    sys.exit(f"FAIL {out_dir} already exists and is not empty; choose another --out")
out_dir.mkdir(parents=True, exist_ok=True)
t0 = time.time()
_run_pipeline_thread(cfg, run_cfg, out_dir / "progress.json")
prog = json.loads((out_dir / "progress.json").read_text())
print("-" * 60)
for st, status in (prog.get("stages") or {}).items():
    print(f"  {st:12} {status}")
print(f"overall: {prog.get('overall')}  ({time.time() - t0:.1f}s)")
for line in [l for l in (prog.get("log_lines") or []) if l.get("level") == "ERR"][:15]:
    print("  ERR", line.get("msg"))
