#!/usr/bin/env python
"""Build the small paired-end fixture the phased fix plan iterates on.

seqtk is not installed on this box, and installing it would add a dependency for
what is three lines of reservoir sampling. This does the same job with two
properties seqtk gives us that actually matter here:

  * R1 and R2 stay in lockstep — the same read indices are kept from both mates,
    so PEAR still sees real pairs. Sampling the files independently would produce
    a fixture that silently fails to assemble.
  * The sample is deterministic — a fixed seed, so the fixture is a stable input.
    A fixture that changed between runs would make golden/ diffs meaningless,
    which is the whole point of the exercise.

Reads are drawn uniformly across the file rather than taken from the head: the
first N reads of an Illumina lane are one corner of one tile and are not
representative of the run's error profile.

Usage:
    python scripts/make_fixture.py [--n 5000] [--seed 1] [--out fixtures/]
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path

# The reference run's raw inputs (runs/projeto_1/20260716T080902Z read these).
DEFAULT_SOURCES = [
    ("APTL64B", Path("runs/projeto_1/20260715T151700Z/input/APTL64B_R1.fastq.gz"),
                Path("runs/projeto_1/20260715T151700Z/input/APTL64B_R2.fastq.gz")),
    ("APTL64F", Path("runs/projeto_1/20260715T151700Z/input/APTL64F_R1.fastq.gz"),
                Path("runs/projeto_1/20260715T151700Z/input/APTL64F_R2.fastq.gz")),
]


def _records(path: Path):
    """Yield 4-line FASTQ records as single strings."""
    with gzip.open(path, "rt") as handle:
        while True:
            block = [handle.readline() for _ in range(4)]
            if not block[0]:
                return
            yield "".join(block)


def _count(path: Path) -> int:
    total = 0
    with gzip.open(path, "rt") as handle:
        for _ in handle:
            total += 1
    return total // 4


def subsample_pair(name: str, r1: Path, r2: Path, out_dir: Path, n: int, seed: int) -> dict:
    total = _count(r1)
    keep_n = min(n, total)
    # Choose the indices ONCE and apply them to both mates: this is what keeps
    # the pairing intact.
    rng = random.Random(seed)
    keep = set(rng.sample(range(total), keep_n))

    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    for mate, src in (("R1", r1), ("R2", r2)):
        dest = out_dir / f"{name}_{mate}.fastq.gz"
        count = 0
        # mtime=0 so the gzip container is byte-stable across rebuilds.
        with gzip.GzipFile(filename="", mode="wb", fileobj=dest.open("wb"), mtime=0) as gz:
            for idx, record in enumerate(_records(src)):
                if idx in keep:
                    gz.write(record.encode())
                    count += 1
        written[mate] = count
    return {"sample": name, "source_reads": total, "kept": written["R1"], "r1": written["R1"], "r2": written["R2"]}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--out", type=Path, default=Path("fixtures"))
    args = parser.parse_args()

    rows = []
    for name, r1, r2 in DEFAULT_SOURCES:
        if not r1.exists() or not r2.exists():
            print(f"FAIL  missing source for {name}: {r1} / {r2}")
            return 1
        row = subsample_pair(name, r1, r2, args.out, args.n, args.seed)
        rows.append(row)
        ok = row["r1"] == row["r2"] == min(args.n, row["source_reads"])
        print(f"{'OK  ' if ok else 'FAIL'}  {name}: {row['source_reads']} -> R1={row['r1']} R2={row['r2']}")
        if not ok:
            return 1

    print(f"\nFixture written to {args.out}/ (seed={args.seed}, n={args.n})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
