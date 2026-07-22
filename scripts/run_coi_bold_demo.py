"""Drive a subsampled COI + BOLD run through the SAME pipeline the UI drives.

Mirrors scripts/run_pipeline.py's build_stores, but sets marker=COI and db_sources=[bold]:
the demonstration is whether the (now case-fixed) BOLD parse identifies COI features on the
APTL64B/APTL64F samples whose biology we already know from the 18S runs (Micromesistius host).

Subsampled inputs (2k reads/sample) are produced by make_fixture.subsample_pair beforehand.
No primer preset is set — a first path-confirming pass leaves the COI reads untrimmed
(cutadapt is a safe no-op without primers; the analysis strip stays at 0,0), so this proves
the plumbing and that BOLD returns identifications. The full run should name the COI primer
pair (Folmer vs Leray) actually used.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT = ROOT / "runs/_coi/current"
STAGES = ["qc", "preprocess", "analysis", "databases"]
ALL_STAGES = ["qc", "preprocess", "analysis", "databases", "phylo"]
INPUTS = {
    "aptl64b_02": ("APTL64B", OUT / "input/APTL64B_R1.fastq.gz", OUT / "input/APTL64B_R2.fastq.gz"),
    "aptl64f_01": ("APTL64F", OUT / "input/APTL64F_R1.fastq.gz", OUT / "input/APTL64F_R2.fastq.gz"),
}


def main() -> int:
    for internal_id, (_, r1, r2) in INPUTS.items():
        for path in (r1, r2):
            if not path.exists():
                print(f"FAIL  missing subsampled input for {internal_id}: {path}")
                return 1

    run_id = OUT.name
    samples_data = {
        "samples": [
            {
                "internal_id": internal_id,
                "display_label": label,
                "r1_path": str(r1),
                "r2_path": str(r2),
                "mode": "PE",
                "marker": "COI",
                "selected": True,
                "is_control": False,
                "group": None,
            }
            for internal_id, (label, r1, r2) in INPUTS.items()
        ]
    }
    project_data = {"project_name": "coi_demo", "analysis_method": "vsearch"}
    run_config_data = {
        "run_id": run_id,
        "output_dir": str(OUT),
        "enabled_stages": {stage: (stage in STAGES) for stage in ALL_STAGES},
        "analysis_method": "vsearch",
        "db_sources": ["bold"],
        "vsearch_params": {"threads": 4},
    }

    from ui.callbacks_progress import _run_pipeline_thread
    from ui.runtime_state import run_config_from_store

    run_config = run_config_from_store(project_data, samples_data, run_config_data)
    run_config.run_id = run_id
    run_config.output_dir = OUT

    print(f"run_id     : {run_id}")
    print(f"marker     : COI   db_sources: bold   stages: {','.join(STAGES)}")
    print(f"out_dir    : {OUT}")
    print("-" * 60, flush=True)

    started = time.time()
    _run_pipeline_thread(run_config, run_config_data, OUT / "progress.json")
    elapsed = time.time() - started

    progress = json.loads((OUT / "progress.json").read_text())
    print("-" * 60)
    for stage, status in (progress.get("stages") or {}).items():
        print(f"  {stage:12} {status}")
    print(f"\noverall: {progress.get('overall', '?')}  ({elapsed:.1f}s)")
    errors = [line for line in (progress.get("log_lines") or []) if line.get("level") == "ERR"]
    if errors:
        print(f"\n{len(errors)} error line(s):")
        for line in errors[:20]:
            print(f"  ERR {line.get('msg')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
