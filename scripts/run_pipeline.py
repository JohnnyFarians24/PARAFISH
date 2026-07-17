#!/usr/bin/env python
"""Run the PARAFISH pipeline end to end, headless, from one command.

Why this exists: the pipeline had no entry point outside the Dash UI, so the only
way to produce a run was to click through the app. That makes a run impossible to
reproduce exactly, and reproducibility is the precondition for the golden/ diffs
the fix plan is built on.

This deliberately drives the *same* code the UI drives — `run_config_from_store`
to build the config and `_run_pipeline_thread` to execute it — rather than
re-implementing the stage order. A runner that took its own path could pass while
the app stayed broken.

`_run_pipeline_thread` is called directly (not via `_launch_pipeline_run`) so the
run is synchronous: the process exits when the pipeline is done, which is what
makes it usable from a script.

Usage:
    python scripts/run_pipeline.py --out runs/_fixture/current            # fixture
    python scripts/run_pipeline.py --real --out runs/_real/current        # real data
    python scripts/run_pipeline.py --stages qc,preprocess,analysis        # subset
    python scripts/run_pipeline.py --db-sources silva,pr2                 # subset
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

ALL_STAGES = ["qc", "preprocess", "analysis", "databases", "phylo"]

# The fixture and the real inputs describe the SAME two samples; only the read
# depth differs. Sample ids match the reference run so golden/ diffs line up.
FIXTURE_INPUTS = {
    "aptl64b_02": ("APTL64B", ROOT / "fixtures/APTL64B_R1.fastq.gz", ROOT / "fixtures/APTL64B_R2.fastq.gz"),
    "aptl64f_01": ("APTL64F", ROOT / "fixtures/APTL64F_R1.fastq.gz", ROOT / "fixtures/APTL64F_R2.fastq.gz"),
}
REAL_INPUTS = {
    "aptl64b_02": ("APTL64B", ROOT / "runs/projeto_1/20260715T151700Z/input/APTL64B_R1.fastq.gz",
                              ROOT / "runs/projeto_1/20260715T151700Z/input/APTL64B_R2.fastq.gz"),
    "aptl64f_01": ("APTL64F", ROOT / "runs/projeto_1/20260715T151700Z/input/APTL64F_R1.fastq.gz",
                              ROOT / "runs/projeto_1/20260715T151700Z/input/APTL64F_R2.fastq.gz"),
}


def build_stores(inputs: dict, out_dir: Path, stages: list, db_sources: list, run_id: str, vsearch_params: dict | None = None) -> tuple:
    samples_data = {
        "samples": [
            {
                "internal_id": internal_id,
                "display_label": label,
                "r1_path": str(r1),
                "r2_path": str(r2),
                "mode": "PE",
                "marker": "18S",
                "selected": True,
                "is_control": False,
                "group": None,
            }
            for internal_id, (label, r1, r2) in inputs.items()
        ]
    }
    project_data = {"project_name": "projeto_1", "analysis_method": "vsearch"}
    run_config_data = {
        "run_id": run_id,
        "output_dir": str(out_dir),
        "enabled_stages": {stage: (stage in stages) for stage in ALL_STAGES},
        "analysis_method": "vsearch",
        "db_sources": db_sources,
        # primer_preset is deliberately left unset so the run exercises the real
        # default ("auto" -> resolved from the 18S marker). Pinning it here would
        # hide exactly the bug that made the reference run keep its primers.
        "vsearch_params": dict(vsearch_params or {"threads": 4}),
    }
    return project_data, samples_data, run_config_data


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--real", action="store_true", help="use the full-depth inputs instead of fixtures/")
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--stages", default="qc,preprocess,analysis,databases")
    parser.add_argument("--db-sources", default="silva,pr2,bold,ncbi,custom")
    parser.add_argument("--keep", action="store_true", help="do not wipe the output dir first")
    parser.add_argument(
        "--min-unique-size",
        type=int,
        default=None,
        help="vsearch --minuniquesize: 1 keeps singletons as OTU seeds, 2 (default) does not",
    )
    args = parser.parse_args()

    stages = [s.strip() for s in args.stages.split(",") if s.strip()]
    db_sources = [s.strip() for s in args.db_sources.split(",") if s.strip()]
    inputs = REAL_INPUTS if args.real else FIXTURE_INPUTS
    out_dir = args.out or (ROOT / ("runs/_real/current" if args.real else "runs/_fixture/current"))
    out_dir = Path(out_dir)

    for internal_id, (_, r1, r2) in inputs.items():
        for path in (r1, r2):
            if not path.exists():
                print(f"FAIL  missing input for {internal_id}: {path}")
                print("      (build the fixture first: python scripts/make_fixture.py)")
                return 1

    if out_dir.exists() and not args.keep:
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # A FIXED run_id, not a timestamp: the output path must be stable so that
    # diffing this run against golden/ is a plain directory diff.
    run_id = out_dir.name

    from ui.callbacks_progress import _run_pipeline_thread
    from ui.runtime_state import run_config_from_store

    vsearch_params = {"threads": 4}
    if args.min_unique_size is not None:
        vsearch_params["min_unique_size"] = args.min_unique_size

    project_data, samples_data, run_config_data = build_stores(inputs, out_dir, stages, db_sources, run_id, vsearch_params)
    run_config = run_config_from_store(project_data, samples_data, run_config_data)
    run_config.run_id = run_id
    run_config.output_dir = out_dir

    print(f"run_id     : {run_id}")
    print(f"inputs     : {'REAL' if args.real else 'FIXTURE'}")
    print(f"out_dir    : {out_dir}")
    print(f"stages     : {','.join(stages)}")
    print(f"db_sources : {','.join(db_sources)}")
    print("-" * 60, flush=True)

    started = time.time()
    _run_pipeline_thread(run_config, run_config_data, out_dir / "progress.json")
    elapsed = time.time() - started

    progress = json.loads((out_dir / "progress.json").read_text())
    print("-" * 60)
    for stage, status in (progress.get("stages") or {}).items():
        print(f"  {stage:12} {status}")
    overall = progress.get("overall", "?")
    print(f"\noverall: {overall}  ({elapsed:.1f}s)")

    # Surface errors here rather than making the caller open progress.json.
    errors = [line for line in (progress.get("log_lines") or []) if line.get("level") == "ERR"]
    if errors:
        print(f"\n{len(errors)} error line(s):")
        for line in errors[:15]:
            print(f"  ERR {line.get('msg')}")

    return 0 if overall not in ("failure", "failed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
