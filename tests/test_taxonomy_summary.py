"""The per-feature taxonomy summary: one row per feature, both evidence families side by side.

The deliverable table. One row per feature carrying the consensus call (LCA) and the best hit
of EACH evidence family — alignment (identity) and SINTAX (bootstrap) — in SEPARATE columns,
with NO cross-metric winner. It is written by the pipeline (``taxonomy_summary.csv``), frozen
in the golden, and rendered as-is by the UI, so the screen and the download cannot diverge.

Three proofs:
 * the builder, including the guard that an aligner is never mistaken for a SINTAX hit — its
   ``confidence`` field falls back to the BLAST bitscore, so "has confidence" alone would pick
   NCBI as the best SINTAX call;
 * the strong golden proof: records collected from the raw SINTAX results the pipeline's own
   collectors read reproduce ``golden/fixture_asv``'s summary byte-for-byte;
 * the UI is single-source: render and download resolve the SAME file, and both callbacks can
   see ``store-project`` so a reconstructed session does not blank them.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# Source from the COMMITTED golden tree, not runs/ (which is gitignored): the fixture_asv
# contract travels with the repo, so these proofs run on a clean clone rather than skipping —
# a skip is not a pass. The raw SINTAX results live only under runs/, so the collector-path
# proof below guards on their presence.
GOLDEN_ASV = ROOT / "golden" / "fixture_asv" / "taxonomy" / "databases"
FROZEN_ASV_SUMMARY = GOLDEN_ASV / "asv" / "taxonomy_summary.csv"
ASV_RAW_RUN = ROOT / "runs" / "_fixture_asv" / "current" / "databases" / "asv"


def _reconstruct_records(assignments_csv: Path) -> dict:
    """Invert source_assignments.csv back to records_by_asv — its 1:1 serialization.

    Proven byte-identical to the pipeline's live collector path (see the collector-path test),
    so feeding these records to the real builder is the pipeline's own computation, replayed
    from the committed contract.
    """
    records: dict = defaultdict(list)
    with assignments_csv.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rank = (row.get("rank") or "").strip()
            record = {
                "source": row.get("source", ""),
                "similarity": row.get("identity", ""),
                "coverage": row.get("coverage", ""),
                "confidence": row.get("confidence", ""),
            }
            taxon = (row.get("taxon") or "").strip()
            if rank and rank != "unassigned" and taxon:
                record[rank] = taxon
            records[row["ASV_ID"]].append(record)
    return records


class SummaryBuilderTests(unittest.TestCase):
    def _summary(self, records):
        from databases.taxonomy_runner import _build_consensus_row, _build_summary_rows

        consensus = [_build_consensus_row(fid, recs, "lca", "OTU") for fid, recs in sorted(records.items())]
        return {row["ASV_ID"]: row for row in _build_summary_rows(records, consensus, "OTU")}

    def test_both_families_land_in_separate_columns(self):
        # NCBI is an aligner: it carries an identity AND a confidence that is really the BLAST
        # bitscore (717). SILVA/PR2 are SINTAX: a bootstrap confidence, no identity.
        records = {
            "F1": [
                {"source": "NCBI", "species": "Gadus morhua", "genus": "Gadus", "similarity": "99.0", "coverage": "100", "confidence": "717"},
                {"source": "SILVA", "species": "Gadus morhua", "genus": "Gadus", "confidence": "0.95"},
                {"source": "PR2", "species": "Gadus ogac", "genus": "Gadus", "confidence": "0.80"},
            ],
        }
        row = self._summary(records)["F1"]
        # Alignment family: NCBI, its identity, its taxon.
        self.assertEqual(row["align_source"], "NCBI")
        self.assertEqual(row["align_identity"], "99.0")
        self.assertEqual(row["align_taxon"], "Gadus morhua")
        # SINTAX family: the higher bootstrap (SILVA 0.95 > PR2 0.80). NOT NCBI, whose 717 is a
        # bitscore, not a bootstrap — the cross-metric guard.
        self.assertEqual(row["sintax_source"], "SILVA")
        self.assertEqual(row["sintax_confidence"], "0.95")
        self.assertNotEqual(row["sintax_source"], "NCBI")
        # Consensus: SILVA and PR2 disagree at species (morhua vs ogac), NCBI agrees morhua, so
        # the LCA lands at genus.
        self.assertEqual(row["consensus_taxon"], "Gadus")
        self.assertEqual(row["consensus_rank"], "genus")

    def test_disagreement_is_shown_not_resolved(self):
        # The two families name different taxa. The summary shows BOTH — no winner is picked.
        records = {
            "F1": [
                {"source": "NCBI", "genus": "Aspergillus", "species": "Aspergillus flavus", "similarity": "96.0", "coverage": "100", "confidence": "628"},
                {"source": "SILVA", "genus": "Gadus", "species": "Gadus morhua", "confidence": "1.00"},
            ],
        }
        row = self._summary(records)["F1"]
        self.assertEqual(row["align_taxon"], "Aspergillus flavus")
        self.assertEqual(row["sintax_taxon"], "Gadus morhua")
        self.assertNotEqual(row["align_taxon"], row["sintax_taxon"])

    def test_single_family_leaves_the_other_empty(self):
        only_sintax = {"F1": [{"source": "SILVA", "species": "Gadus morhua", "genus": "Gadus", "confidence": "1.00"}]}
        row = self._summary(only_sintax)["F1"]
        self.assertEqual(row["align_source"], "")
        self.assertEqual(row["align_taxon"], "")
        self.assertEqual(row["align_identity"], "")
        self.assertEqual(row["sintax_source"], "SILVA")

        only_align = {"F1": [{"source": "NCBI", "species": "Gadus morhua", "genus": "Gadus", "similarity": "99.0", "coverage": "100", "confidence": "700"}]}
        row = self._summary(only_align)["F1"]
        self.assertEqual(row["sintax_source"], "")
        self.assertEqual(row["sintax_taxon"], "")
        self.assertEqual(row["sintax_confidence"], "")
        self.assertEqual(row["align_source"], "NCBI")


class SummaryGoldenTests(unittest.TestCase):
    def _write_summary(self, records, feature_type):
        from databases.taxonomy_runner import SUMMARY_COLUMNS, _build_consensus_row, _build_summary_rows, _write_csv

        consensus = [_build_consensus_row(a, r, "lca", feature_type) for a, r in sorted(records.items())]
        rows = _build_summary_rows(records, consensus, feature_type)
        out = Path(tempfile.mkdtemp()) / "taxonomy_summary.csv"
        _write_csv(out, rows, SUMMARY_COLUMNS)
        return out

    def test_builder_reproduces_frozen_summary_from_committed_contract(self):
        """The builder, over records replayed from the committed source_assignments.csv,
        reproduces the frozen fixture_asv summary byte-for-byte — runs on a clean clone."""
        records = _reconstruct_records(GOLDEN_ASV / "asv" / "source_assignments.csv")
        out = self._write_summary(records, "ASV")
        self.assertEqual(out.read_bytes(), FROZEN_ASV_SUMMARY.read_bytes())

    def test_pipeline_collector_path_matches(self):
        """Stronger, dev-machine-only: records the pipeline's OWN collectors read from the raw
        SINTAX results reproduce the same frozen summary. Skips when runs/ is absent (the raw
        results are gitignored), so a clean clone does not fail on their absence."""
        from databases.taxonomy_runner import _collect_sintax_records

        raws = [ASV_RAW_RUN / "silva" / "silva_sintax_result.csv", ASV_RAW_RUN / "pr2" / "pr2_sintax_result.csv"]
        if not all(p.exists() for p in raws):
            self.skipTest("raw SINTAX results not present (runs/ is gitignored)")
        records = defaultdict(list)
        for source, path in (("SILVA", raws[0]), ("PR2", raws[1])):
            for asv_id, rows in _collect_sintax_records(source, path).items():
                records[asv_id].extend(rows)
        out = self._write_summary(records, "ASV")
        self.assertEqual(out.read_bytes(), FROZEN_ASV_SUMMARY.read_bytes())


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    # Reconstruct a run dir from the COMMITTED golden tree, so the resolver proof runs on a
    # clean clone. taxonomy_file resolves databases/<set>/taxonomy_summary.csv, which the
    # golden databases/asv/ subtree provides.
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(GOLDEN_ASV, run_dir / "databases", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class SummaryUiSingleSourceTests(unittest.TestCase):
    def setUp(self):
        self._prev = os.environ.get("RUNS_DIR")
        self.tmp = Path(tempfile.mkdtemp())
        os.environ["RUNS_DIR"] = str(self.tmp)
        _make_run(self.tmp, "restore_proj", "20260715T151700Z")
        self.project_data = {"project_name": "restore_proj"}

    def tearDown(self):
        if self._prev is None:
            os.environ.pop("RUNS_DIR", None)
        else:
            os.environ["RUNS_DIR"] = self._prev
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_render_and_download_resolve_the_same_frozen_file(self):
        from ui.callbacks_databases import _resolve_output_dir
        from ui.feature_sets import taxonomy_file

        run_output_dir = _resolve_output_dir(None, self.project_data)
        self.assertIsNotNone(run_output_dir)
        # The one file both consumers resolve, byte-identical to golden.
        path = taxonomy_file({}, run_output_dir, "asv", "taxonomy_summary_csv")
        self.assertIsNotNone(path)
        self.assertEqual(Path(path).read_bytes(), FROZEN_ASV_SUMMARY.read_bytes())
        # otu on an asv-only run yields nothing (the set is honored).
        self.assertIsNone(taxonomy_file({}, run_output_dir, "otu", "taxonomy_summary_csv"))

    def test_summary_callbacks_depend_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        targets = {"databases-feature-summary.children", "download-summary-csv.data"}
        checked = set()
        for cb in deps:
            out = cb["output"]
            outs = [o.strip(".") for o in out.strip(".").split("...")] if out.startswith("..") else [out]
            for target in targets:
                if target in outs:
                    ids = {d["id"] for d in cb["inputs"]} | {d["id"] for d in cb["state"]}
                    self.assertIn("store-project", ids, f"{target}: cannot see store-project")
                    checked.add(target)
        self.assertEqual(checked, targets, f"target callbacks not found: {targets - checked}")


if __name__ == "__main__":
    unittest.main()
