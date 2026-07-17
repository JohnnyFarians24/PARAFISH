"""Cross-sample comparison against a user-uploaded (custom) reference database.

Given the per-OTU best-hit table produced by ``databases.custom_db.run_custom_blast``
and the OTU/ASV × sample abundance matrix, this builds the outputs the
supervisors asked for:

* a **per-OTU match table** (which OTU matched which species, at what identity,
  and how many reads it carries) where unmatched OTUs are flagged ``"None"``;
* a **taxon × sample detection matrix** collapsing OTUs by their matched species,
  with an explicit ``"No match"`` row for the reads that hit nothing — so a user
  comparing dozens to hundreds of samples immediately sees, per sample, which
  taxa were detected and how much of the community stayed unassigned.

Pure csv/dict (no pandas/plotly) so it stays unit-testable; the UI layer turns
the returned dict into DataTables and heatmaps.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

from analysis.sample_comparison import _load_matrix

NO_MATCH_LABEL = "No match"


def _read_besthit(besthit_csv: Path) -> list[dict]:
    if not besthit_csv.exists() or besthit_csv.stat().st_size == 0:
        return []
    with besthit_csv.open("r", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle))


def _to_float(value: str) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _load_consensus_species(consensus_csv: Optional[Path]) -> dict[str, str]:
    """Map OTU id → public-database consensus species (empty when absent)."""
    if not consensus_csv:
        return {}
    consensus_csv = Path(consensus_csv)
    if not consensus_csv.exists() or consensus_csv.stat().st_size == 0:
        return {}
    try:
        from analysis.community_analysis import _load_taxonomy

        return _load_taxonomy(consensus_csv, "species")
    except Exception:
        return {}


def build_custom_comparison(
    besthit_csv,
    abundance_matrix_path,
    *,
    min_identity: float = 0.0,
    species_filter: Optional[str] = None,
    consensus_csv=None,
) -> dict:
    """Combine the custom BLAST best-hit table with the abundance matrix.

    ``min_identity`` demotes any hit below the given percent identity to
    "unmatched" (so the user can tighten the comparison interactively).
    ``species_filter`` (case-insensitive substring) restricts the per-OTU table
    and detection matrix to matching taxa. ``consensus_csv`` (the public-database
    taxonomy consensus) adds a ``public_species`` column so each OTU shows the
    custom match and the public-database call side by side. Returns
    ``{"error": ...}`` on failure with empty collections so the UI degrades
    gracefully.
    """
    besthit_csv = Path(besthit_csv)
    rows = _read_besthit(besthit_csv)
    if not rows:
        return {"error": "Tabela de melhores hits do BLAST não encontrada ou vazia.", "per_otu": [], "detection": None, "samples": [], "summary": {}}

    feature_ids, samples, values = _load_matrix(Path(abundance_matrix_path))
    if not feature_ids or not samples:
        return {"error": "Matriz de abundância não encontrada ou sem amostras.", "per_otu": [], "detection": None, "samples": [], "summary": {}}

    public_species = _load_consensus_species(consensus_csv)

    index_of = {fid: i for i, fid in enumerate(feature_ids)}
    n_samples = len(samples)
    zero = [0.0] * n_samples

    def _row_for(otu_id: str) -> list[float]:
        idx = index_of.get(otu_id)
        return values[idx] if idx is not None else zero

    # Per-OTU records and taxon collapse.
    per_otu: list[dict] = []
    taxon_reads: dict[str, list[float]] = {}
    matched = 0
    filt = (species_filter or "").strip().lower()

    for row in rows:
        otu_id = (row.get("OTU_ID") or "").strip()
        if not otu_id:
            continue
        status = (row.get("status") or "").strip()
        pident = _to_float(row.get("pident"))
        species = (row.get("species") or "").strip()
        is_match = status == "Match" and species and species != "None"
        if is_match and min_identity > 0 and (pident is None or pident < min_identity):
            is_match = False

        label = species if is_match else NO_MATCH_LABEL
        counts = _row_for(otu_id)
        total_reads = sum(counts)
        n_detected = sum(1 for c in counts if c > 0)

        if is_match:
            matched += 1

        record = {
            "OTU_ID": otu_id,
            "species": label if is_match else "None",
            "public_species": public_species.get(otu_id) or "None",
            "subject_id": (row.get("subject_id") or "").strip(),
            "pident": round(pident, 2) if pident is not None else None,
            "status": "Match" if is_match else "None",
            "total_reads": int(round(total_reads)),
            "n_samples_detected": n_detected,
        }
        if not filt or (filt in label.lower()):
            per_otu.append(record)

        if not filt or (is_match and filt in label.lower()) or (not is_match and filt in NO_MATCH_LABEL.lower()):
            acc = taxon_reads.setdefault(label, [0.0] * n_samples)
            for j in range(n_samples):
                acc[j] += counts[j]

    total = sum(1 for r in rows if (r.get("OTU_ID") or "").strip())

    # Order taxa by total reads (descending), keep "No match" last for readability.
    ordered = sorted(
        (name for name in taxon_reads if name != NO_MATCH_LABEL),
        key=lambda n: sum(taxon_reads[n]),
        reverse=True,
    )
    taxa = ordered + ([NO_MATCH_LABEL] if NO_MATCH_LABEL in taxon_reads else [])
    detection = {
        "taxa": taxa,
        "samples": samples,
        "values": [[int(round(v)) for v in taxon_reads[name]] for name in taxa],
    }

    # Sort per-OTU: matches first (by reads desc), then None (by reads desc).
    per_otu.sort(key=lambda r: (r["status"] != "Match", -r["total_reads"]))

    summary = {
        "n_otus": total,
        "n_matched": matched,
        "n_unmatched": total - matched,
        "pct_matched": round(100.0 * matched / total, 1) if total else 0.0,
        "n_samples": n_samples,
        "n_taxa": len(ordered),
    }
    return {"per_otu": per_otu, "detection": detection, "samples": samples, "taxa": taxa, "summary": summary}


def write_detection_table(detection: dict, out_path) -> Optional[Path]:
    """Write the taxon × sample detection matrix as a downloadable CSV."""
    if not detection or not detection.get("taxa"):
        return None
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon"] + list(detection["samples"]))
        for name, row in zip(detection["taxa"], detection["values"]):
            writer.writerow([name] + list(row))
    return out_path
