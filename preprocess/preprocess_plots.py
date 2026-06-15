"""Plot payload builders for preprocessing diagnostics."""

from __future__ import annotations

from preprocess.preprocess_parser import parse_preprocess_outputs


def build_preprocess_plot_data(output_dir: str) -> dict:
    """
    Build a plotting payload from preprocess summaries.

    Returns a dict with:
      - per_sample_counts: list of {sample, trimmed_r1, trimmed_r2, merged}
      - totals: {total_trimmed, total_merged}
    """
    summary = parse_preprocess_outputs(output_dir)
    per_sample = []

    for sample_id, files in summary.get("samples", {}).items():
        trimmed = files.get("trimmed", [])
        merged = files.get("merged", [])
        per_sample.append(
            {
                "sample": sample_id,
                "trimmed_r1": next((f for f in trimmed if "_r1_" in f.lower()), None),
                "trimmed_r2": next((f for f in trimmed if "_r2_" in f.lower()), None),
                "merged": merged[0] if merged else None,
                "n_trimmed": len(trimmed),
                "n_merged": len(merged),
            }
        )

    return {
        "per_sample": per_sample,
        "totals": {
            "total_samples": len(per_sample),
            "samples_with_merged": sum(1 for sample in per_sample if sample["n_merged"] > 0),
            "samples_with_trimmed": sum(1 for sample in per_sample if sample["n_trimmed"] > 0),
        },
    }
