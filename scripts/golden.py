#!/usr/bin/env python
"""Freeze a closed stage's output, and diff later runs against it.

This is the inter-stage contract from the fix plan. Once a stage is closed its output
is frozen here; every later change re-runs the pipeline and diffs against these files.
If a number from a CLOSED stage moves, that is a regression and the work stops until
it is explained.

What gets frozen is the stage's *contract* — the tables downstream stages and the
reader actually consume — plus checksums for the bulk outputs (FASTQ/FASTA), which are
too large to copy but whose identity still has to be pinned. A changed checksum is a
changed input to the next stage, and that is exactly what must not happen silently.

Golden sets are namespaced by dataset: the fixture and the real data legitimately
produce different numbers, and comparing one to the other would be meaningless.

Usage:
    python scripts/golden.py freeze preprocess --from runs/_real/current --dataset real
    python scripts/golden.py check  preprocess --from runs/_real/current --dataset real
    python scripts/golden.py check  all        --from runs/_real/current --dataset real
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GOLDEN = ROOT / "golden"

# Per stage: the exact files to copy (the contract) and the bulk files to checksum.
STAGES = {
    "preprocess": {
        "copy": ["preprocess/metrics/read_tracking.csv", "preprocess/metrics/assembled_length_distribution.csv"],
        "checksum": ["preprocess/*/pear/*.assembled.fastq"],
    },
    "analysis": {
        "copy": [
            "analysis/pooled/otu_table.tsv",
            "analysis/pooled/otu_id_map.csv",
            "analysis/otu_matrix.csv",
            "analysis/metrics/otu/chimera_stats.csv",
            "analysis/metrics/otu/alpha_diversity.csv",
            "analysis/metrics/otu/rank_abundance.csv",
        ],
        "checksum": ["analysis/pooled/otus.fasta"],
    },
    "taxonomy": {
        "copy": [
            "databases/otu/taxonomy_consensus.csv",
            "databases/otu/source_assignments.csv",
            "databases/otu/source_status.csv",
            "databases/otu/conflict_report.csv",
            "databases/otu/metrics/taxonomic_resolution.csv",
        ],
        "checksum": [],
    },
    "export": {
        "copy": [
            "databases/otu/phyloseq/tax_table.tsv",
            "databases/otu/phyloseq/otu_table.tsv",
            "databases/otu/phyloseq/sample_data.tsv",
            "databases/otu/taxa_tables/taxa_table_class.csv",
            "databases/otu/taxa_tables/taxa_table_genus.csv",
            "databases/otu/taxa_tables/taxa_table_species.csv",
        ],
        "checksum": ["databases/otu/table.biom"],
    },
}
ORDER = ["preprocess", "analysis", "taxonomy", "export"]


# Artifacts that embed generation metadata which is NOT content. BIOM carries a
# wall-clock ``date``, so its bytes change on every run while the table it holds does
# not — which made ``table.biom`` permanently undiffable and its checksum guaranteed to
# fail on any re-run. That is the exact defect conflict_report.csv had: a file that
# changes without its content changing cannot be diffed, so a real change would hide in
# the churn. It went unnoticed because Phase 4's two-run proof covered the nine COPIED
# contract files; table.biom is checksummed, not copied, and was never in that diff.
#
# Fixed here rather than in write_biom on purpose: the date is TRUE and downstream tools
# expect it. Making the output lie about when it was written, to keep a checker quiet,
# is the wrong trade for a platform whose whole argument is that its outputs are honest.
# So the field stays in the file; it just does not count as content.
_VOLATILE_JSON_FIELDS = {".biom": ("date",)}


def _md5(path: Path) -> str:
    volatile = _VOLATILE_JSON_FIELDS.get(path.suffix)
    if volatile:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass  # not the JSON we expected — hash the raw bytes and let it differ loudly
        else:
            if isinstance(payload, dict):
                for field in volatile:
                    payload.pop(field, None)
                canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
                return hashlib.md5(canonical.encode("utf-8")).hexdigest()
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _checksum_lines(run: Path, patterns: list[str]) -> list[str]:
    lines = []
    for pattern in patterns:
        for path in sorted(run.glob(pattern)):
            if path.is_file():
                lines.append(f"{_md5(path)}  {path.relative_to(run)}")
    return lines


def freeze(stage: str, run: Path, dataset: str) -> int:
    spec = STAGES[stage]
    dest = GOLDEN / dataset / stage
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)

    copied = 0
    for rel in spec["copy"]:
        src = run / rel
        if not src.exists():
            print(f"  --    {rel} (absent, not frozen)")
            continue
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, target)
        copied += 1
        print(f"  ok    {rel}")

    lines = _checksum_lines(run, spec["checksum"])
    if lines:
        (dest / "CHECKSUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"  ok    CHECKSUMS.txt ({len(lines)} file(s))")
    (dest / "PROVENANCE.txt").write_text(f"stage: {stage}\ndataset: {dataset}\nsource: {run}\n", encoding="utf-8")

    if not copied:
        print(f"FAIL  nothing frozen for {stage} — did the stage run?")
        return 1
    print(f"\nfrozen: golden/{dataset}/{stage}/  ({copied} file(s))")
    return 0


def check(stage: str, run: Path, dataset: str) -> int:
    dest = GOLDEN / dataset / stage
    if not dest.exists():
        print(f"  SKIP  {stage:11} not frozen yet")
        return 0

    problems = []
    for rel in STAGES[stage]["copy"]:
        golden_file, new_file = dest / rel, run / rel
        if not golden_file.exists():
            continue
        if not new_file.exists():
            problems.append(f"{rel}: MISSING from the new run")
        elif golden_file.read_bytes() != new_file.read_bytes():
            problems.append(f"{rel}: CONTENT CHANGED")

    checksums = dest / "CHECKSUMS.txt"
    if checksums.exists():
        expected = dict(
            (line.split("  ", 1)[1], line.split("  ", 1)[0])
            for line in checksums.read_text(encoding="utf-8").splitlines()
            if "  " in line
        )
        for rel, digest in expected.items():
            path = run / rel
            if not path.exists():
                problems.append(f"{rel}: MISSING from the new run")
            elif _md5(path) != digest:
                problems.append(f"{rel}: CHECKSUM CHANGED")

    if problems:
        print(f"  FAIL  {stage:11} {len(problems)} regression(s) against golden/{dataset}/{stage}")
        for problem in problems:
            print(f"          {problem}")
        return 1
    print(f"  PASS  {stage:11} identical to golden/{dataset}/{stage}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["freeze", "check"])
    parser.add_argument("stage", choices=[*ORDER, "all"])
    parser.add_argument("--from", dest="run", type=Path, required=True)
    parser.add_argument("--dataset", default="real")
    args = parser.parse_args()

    if not args.run.exists():
        print(f"FAIL  run dir does not exist: {args.run}")
        return 2

    if args.action == "freeze":
        if args.stage == "all":
            print("FAIL  freeze one stage at a time — that is the point of the plan")
            return 2
        return freeze(args.stage, args.run, args.dataset)

    stages = ORDER if args.stage == "all" else [args.stage]
    print(f"non-regression: {args.run} vs golden/{args.dataset}/")
    failed = sum(check(stage, args.run, args.dataset) for stage in stages)
    print(f"\n{'REGRESSION' if failed else 'no regressions'}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
