"""Persisted preprocess metrics: the read-tracking waterfall.

Trimmomatic and PEAR both print their read accounting to their console logs, but
nothing ever read it back, so the single most important preprocessing output —
how many reads survive each step — did not exist. This module parses those logs
once, in the runner, and writes the tables the UI reads.

Counts for paired-end samples are **read pairs**, matching how Trimmomatic and
PEAR report them; a "read" in a PE waterfall is one pair, not one FASTQ record.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import Counter
from pathlib import Path
from typing import Iterator, Optional

# A step that keeps less than this fraction of what it was given is a red flag:
# outside the initial quality filter, no step should discard most of the data.
MIN_STEP_RETENTION = 0.5

# The waterfall, in pipeline order. Every stage that can drop a read has a column, so
# a read that leaves the pipeline leaves a trace: input -> primer -> quality -> merge
# -> filter -> mapping.
#
# `filtered_*` and `mapped_*` are produced by the ANALYSIS stage, not preprocess, and
# are filled in only when build_preprocess_metrics() is given an analysis_dir. They live
# here rather than in a second table because splitting one waterfall across two files is
# what let 22k reads go missing: the preprocess table ended at "assembled" and nothing
# downstream ever reconciled against it.
#
# Dereplication and chimera removal are deliberately NOT columns. They are pooled
# operations, and --minuniquesize excludes reads from *clustering* while the mapping
# step counts them anyway — so a per-sample "dereplicated" column would show a loss
# that does not exist. That accounting is pooled and belongs in the vsearch pool log.
TRACKING_COLUMNS = [
    "sample",
    "mode",
    "input_reads",
    "reads_with_primer",
    "pct_primer",
    "primer_trimmed_reads",
    "trimmed_reads",
    "forward_only",
    "reverse_only",
    "dropped_reads",
    "assembled_reads",
    "discarded_reads",
    "not_assembled_reads",
    "filtered_reads",
    "filter_discarded_reads",
    "mapped_reads",
    "unmapped_reads",
    "surviving_reads",
    "pct_trimmed",
    "pct_assembled",
    "pct_filtered",
    "pct_mapped",
    "pct_surviving",
]

_TRIM_PE = re.compile(
    r"Input Read Pairs:\s*(\d+)\s+Both Surviving:\s*(\d+)\s*\([\d.]+%\)"
    r"\s+Forward Only Surviving:\s*(\d+)\s*\([\d.]+%\)"
    r"\s+Reverse Only Surviving:\s*(\d+)\s*\([\d.]+%\)"
    r"\s+Dropped:\s*(\d+)"
)
_TRIM_SE = re.compile(r"Input Reads:\s*(\d+)\s+Surviving:\s*(\d+)\s*\([\d.]+%\)\s+Dropped:\s*(\d+)")

# vsearch --fastq_filter: "82113 sequences kept (of which 0 truncated), 12544 sequences discarded."
_FILTER_KEPT = re.compile(r"(\d+)\s+sequences kept\s*\(of which \d+ truncated\),\s*(\d+)\s+sequences discarded")


def _pear_line(log_text: str, label: str) -> Optional[int]:
    """PEAR prints e.g. ``Assembled reads ...........: 94,657 / 99,618 (95.020%)``."""
    match = re.search(rf"{label}\s*\.*\s*:\s*([\d,]+)\s*/", log_text)
    return int(match.group(1).replace(",", "")) if match else None


def parse_trimmomatic_log(log_path) -> dict:
    """Return the read accounting Trimmomatic printed, or ``{}``."""
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {}
    pe = _TRIM_PE.search(text)
    if pe:
        total, both, forward, reverse, dropped = (int(value) for value in pe.groups())
        return {
            "mode": "PE",
            "input_reads": total,
            "trimmed_reads": both,
            "forward_only": forward,
            "reverse_only": reverse,
            "dropped_reads": dropped,
        }
    se = _TRIM_SE.search(text)
    if se:
        total, surviving, dropped = (int(value) for value in se.groups())
        return {
            "mode": "SE",
            "input_reads": total,
            "trimmed_reads": surviving,
            "forward_only": 0,
            "reverse_only": 0,
            "dropped_reads": dropped,
        }
    return {}


def parse_pear_log(log_path) -> dict:
    """Return PEAR's assembled/discarded/not-assembled counts, or ``{}``."""
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {}
    assembled = _pear_line(text, "Assembled reads")
    if assembled is None:
        return {}
    return {
        "assembled_reads": assembled,
        "discarded_reads": _pear_line(text, "Discarded reads") or 0,
        "not_assembled_reads": _pear_line(text, "Not assembled reads") or 0,
    }


def parse_cutadapt_report(json_path) -> dict:
    """Return ``{input, output_reads, reads_with_primer}`` from cutadapt's --json report.

    ``reads_with_primer`` is how many reads (pairs, for PE) actually carried the primer
    that was removed — the honest "primers found" number for the read-tracking table.
    ``output_reads`` is how many survived the step, which is a different question: a read
    can carry no primer and still pass through.
    """
    import json as _json

    try:
        data = _json.loads(Path(json_path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    read_counts = data.get("read_counts") or {}
    result: dict = {}
    for key, json_key in (("input", "input"), ("output_reads", "output"), ("reads_with_primer", "read1_with_adapter")):
        value = read_counts.get(json_key)
        if isinstance(value, (int, float)):
            result[key] = int(value)
    return result


def parse_vsearch_filter_log(log_path) -> dict:
    """Return ``{filtered_reads, filter_discarded_reads}`` from a vsearch --fastq_filter log.

    This is the step that actually swallowed the "~22k reads with no trace": the quality
    filter runs in the Analysis stage, so the preprocess waterfall never saw it and the
    reads simply stopped being mentioned between "assembled" and the OTU table.
    """
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {}
    match = _FILTER_KEPT.search(text)
    if not match:
        return {}
    return {"filtered_reads": int(match.group(1)), "filter_discarded_reads": int(match.group(2))}


def mapped_reads_per_sample(otu_table_path) -> dict[str, int]:
    """``{sample: reads assigned to some OTU}`` — each column's sum in the OTU table.

    The OTU table is the canonical abundance basis (see DEVLOG, "Single source of
    abundance truth"), so the waterfall's last step is read from it rather than
    recomputed, and cannot drift from it.
    """
    path = Path(otu_table_path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    rows = [line.split("\t") for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.strip()]
    if len(rows) < 1:
        return {}
    samples = [name.strip() for name in rows[0][1:]]
    totals = {sample: 0 for sample in samples}
    for row in rows[1:]:
        for index, value in enumerate(row[1:]):
            if index < len(samples):
                try:
                    totals[samples[index]] += int(float(value or 0))
                except ValueError:
                    continue
    return totals


def _open_fastq(path: Path):
    return gzip.open(path, "rt", errors="replace") if path.suffix == ".gz" else path.open("r", errors="replace")


def _sequence_lengths(path: Path) -> Iterator[int]:
    """Yield the length of every record in a FASTQ, streaming."""
    with _open_fastq(path) as handle:
        for index, line in enumerate(handle):
            if index % 4 == 1:
                yield len(line.strip())


def assembled_length_counts(path) -> Counter:
    """Length histogram of a merged FASTQ. Streams: these files run to gigabytes."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return Counter()
    return Counter(_sequence_lengths(path))


def _pct(numerator: Optional[int], denominator: Optional[int]) -> str:
    if not denominator or numerator is None:
        return ""
    return str(round(100.0 * numerator / denominator, 2))


def _find_assembled(pear_dir: Path) -> Optional[Path]:
    for candidate in sorted(pear_dir.glob("*.assembled.fastq*")):
        if candidate.stat().st_size > 0:
            return candidate
    return None


def _write_csv(path: Path, columns: list[str], rows: list[dict]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    return path


def build_preprocess_metrics(preprocess_dir, analysis_dir=None) -> dict[str, Path]:
    """Parse every sample's Trimmomatic/PEAR logs and persist the waterfall tables.

    ``analysis_dir`` is optional because this runs at the end of the Preprocess stage,
    when the Analysis outputs do not exist yet. The runner calls it a second time once
    Analysis has finished, passing the directory, and the filter/mapping columns are
    filled in then. Same writer both times — read_tracking.csv has exactly one author.
    """
    preprocess_dir = Path(preprocess_dir)
    if not preprocess_dir.exists():
        return {}

    analysis_dir = Path(analysis_dir) if analysis_dir else None
    mapped_by_sample = mapped_reads_per_sample(analysis_dir / "pooled" / "otu_table.tsv") if analysis_dir else {}

    tracking: list[dict] = []
    lengths: list[dict] = []

    for sample_dir in sorted(path for path in preprocess_dir.iterdir() if path.is_dir()):
        if sample_dir.name == "metrics":
            continue
        trim = parse_trimmomatic_log(sample_dir / "trimlog.txt")
        if not trim:
            continue

        pear_dir = sample_dir / "pear"
        pear = parse_pear_log(pear_dir / "pear.log") if pear_dir.is_dir() else {}

        # cutadapt ran before Trimmomatic; report how many reads carried the primer.
        cut = parse_cutadapt_report(sample_dir / "cutadapt.json")
        reads_with_primer = cut.get("reads_with_primer")

        # The Analysis stage's per-sample quality filter, and the OTU mapping that ends
        # the waterfall. Absent until Analysis has run.
        filt = parse_vsearch_filter_log(analysis_dir / sample_dir.name / "step1_filter.log") if analysis_dir else {}
        filtered_reads = filt.get("filtered_reads")
        mapped = mapped_by_sample.get(sample_dir.name)
        # Reads that passed the filter and still matched no OTU at the clustering
        # identity. Small, but it is the difference between a closed account and a
        # rounding excuse.
        unmapped = filtered_reads - mapped if (filtered_reads is not None and mapped is not None) else None

        # Without a merge step the trimmed pairs are what carries forward. This stays the
        # END OF PREPROCESS count; `mapped_reads` is the end of the pipeline. They differ
        # by the filter, and conflating them is what hid the loss.
        surviving = pear.get("assembled_reads", trim["trimmed_reads"])
        row = {
            "sample": sample_dir.name,
            **trim,
            "assembled_reads": pear.get("assembled_reads", ""),
            "discarded_reads": pear.get("discarded_reads", ""),
            "not_assembled_reads": pear.get("not_assembled_reads", ""),
            "surviving_reads": surviving,
            "reads_with_primer": reads_with_primer if reads_with_primer is not None else "",
            "primer_trimmed_reads": cut.get("output_reads", ""),
            "pct_primer": _pct(reads_with_primer, cut.get("input") or trim["input_reads"]) if reads_with_primer is not None else "",
            "filtered_reads": filtered_reads if filtered_reads is not None else "",
            "filter_discarded_reads": filt.get("filter_discarded_reads", ""),
            "mapped_reads": mapped if mapped is not None else "",
            "unmapped_reads": unmapped if unmapped is not None else "",
            "pct_trimmed": _pct(trim["trimmed_reads"], trim["input_reads"]),
            "pct_assembled": _pct(pear.get("assembled_reads"), trim["trimmed_reads"]),
            "pct_filtered": _pct(filtered_reads, pear.get("assembled_reads") or trim["trimmed_reads"]),
            "pct_mapped": _pct(mapped, trim["input_reads"]),
            "pct_surviving": _pct(surviving, trim["input_reads"]),
        }
        tracking.append(row)

        assembled = _find_assembled(pear_dir) if pear_dir.is_dir() else None
        if assembled is not None:
            for length, count in sorted(assembled_length_counts(assembled).items()):
                lengths.append({"sample": sample_dir.name, "length": length, "count": count})

    if not tracking:
        return {}

    metrics_dir = preprocess_dir / "metrics"
    written = {
        "read_tracking_csv": _write_csv(metrics_dir / "read_tracking.csv", TRACKING_COLUMNS, tracking),
    }
    if lengths:
        written["assembled_length_csv"] = _write_csv(
            metrics_dir / "assembled_length_distribution.csv", ["sample", "length", "count"], lengths
        )
    return written


def preprocess_metrics_paths(preprocess_dir) -> dict[str, Path]:
    """Return the persisted preprocess tables that exist under ``preprocess_dir``."""
    metrics_dir = Path(preprocess_dir) / "metrics"
    names = {
        "read_tracking_csv": "read_tracking.csv",
        "assembled_length_csv": "assembled_length_distribution.csv",
    }
    return {key: metrics_dir / name for key, name in names.items() if (metrics_dir / name).exists()}


def failing_steps(row: dict) -> list[str]:
    """Names of the steps in ``row`` that kept less than ``MIN_STEP_RETENTION``."""
    failing = []
    # "filtering" belongs here as much as the other two: on the reference run it is the
    # step that discarded 22k reads, and it was the only one nobody was watching.
    for label, pct_key in (("trimming", "pct_trimmed"), ("merging", "pct_assembled"), ("filtering", "pct_filtered")):
        value = row.get(pct_key)
        try:
            if value not in ("", None) and float(value) < MIN_STEP_RETENTION * 100:
                failing.append(label)
        except (TypeError, ValueError):
            continue
    return failing
