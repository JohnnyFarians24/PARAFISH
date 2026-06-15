"""Plot data preparation for QC visualizations."""

from __future__ import annotations

from qc.fastqc_parser import collect_fastqc_reports


def build_qc_plot_data(output_dir: str) -> dict:
    """
    Build a structured plotting payload from a QC output directory.

    Returns:
        {
            "per_sample": [
                {
                    "sample": str,
                    "total_sequences_r1": float | None,
                    "total_sequences_r2": float | None,
                    "pct_gc_r1": float | None,
                    "pct_gc_r2": float | None,
                    "failed_modules_r1": float,
                    "failed_modules_r2": float,
                }
            ],
            "totals": {
                "total_samples": int,
                "samples_with_failures": int,
            }
        }
    """
    reports = collect_fastqc_reports(output_dir)
    per_sample = []

    for sample_id, strands in sorted(reports.items()):
        r1 = strands.get("r1", {})
        r2 = strands.get("r2", {})

        def _float(data: dict, key: str):
            try:
                return float(data[key])
            except (KeyError, TypeError, ValueError):
                return None

        entry = {
            "sample": sample_id,
            "total_sequences_r1": _float(r1, "Total Sequences"),
            "total_sequences_r2": _float(r2, "Total Sequences"),
            "pct_gc_r1": _float(r1, "%GC"),
            "pct_gc_r2": _float(r2, "%GC"),
            "failed_modules_r1": r1.get("failed_modules_r1", 0.0),
            "failed_modules_r2": r2.get("failed_modules_r2", 0.0),
        }
        per_sample.append(entry)

    samples_with_failures = sum(
        1
        for sample in per_sample
        if (sample["failed_modules_r1"] or 0) > 0 or (sample["failed_modules_r2"] or 0) > 0
    )

    return {
        "per_sample": per_sample,
        "totals": {
            "total_samples": len(per_sample),
            "samples_with_failures": samples_with_failures,
        },
    }
