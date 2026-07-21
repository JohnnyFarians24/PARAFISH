"""A fresh session must reconstruct the Analysis page from disk, not blank it.

Eixo B, page 3/4. The selector-FREE Analysis render callbacks read ``output_dir`` off the
per-session ``store-run-config`` and blank on a fresh tab. The OTU golden pins the numbers,
so this is the STRONG form on the callbacks the golden covers.

Two resolvers feed those callbacks and both get the disk fallback:
 * ``_analysis_metric_paths`` (the otu-assuming metrics helper) — chimera, alpha, plots, upset;
 * ``_find_abundance_matrix`` (otu-first abundance table) — otu-table, sample-table.

Kept clear of Eixo A. The ONE feature-set selector (``analysis-occupancy-set`` /
``_render_occupancy``), the ASV-only ``_render_dada2_tracking``, and ``_update_analysis_stats``
(which restores store-recorded FASTA paths, a different shape) are LEFT for the ASV golden
and recorded in the living registry. The proof path here has no feature-set selector in it;
``_render_otu_table``'s ``abundance-mode`` is a presentation axis (relative/absolute), not a
feature-set axis, so it is on the provable side.

Seen to fail before the fix: the resolvers had no ``project_data`` / ``latest_run_dir``
fallback (→ ``{}`` / None), and the six render callbacks did not declare ``store-project``.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FROZEN_CHIMERA = ROOT / "golden" / "real" / "analysis" / "analysis" / "metrics" / "otu" / "chimera_stats.csv"
FROZEN_ALPHA = ROOT / "golden" / "real" / "analysis" / "analysis" / "metrics" / "otu" / "alpha_diversity.csv"
FROZEN_MATRIX = ROOT / "golden" / "real" / "analysis" / "analysis" / "otu_matrix.csv"
SRC_ANALYSIS = ROOT / "runs" / "_real" / "current" / "analysis"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_ANALYSIS, run_dir / "analysis", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class AnalysisRestoreLogicTests(unittest.TestCase):
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

    def test_metric_paths_reconstructs_frozen_chimera_and_alpha(self):
        from ui.callbacks_analysis import _analysis_metric_paths

        paths = _analysis_metric_paths(None, None, self.project_data)
        self.assertIn("chimera_stats_csv", paths)
        self.assertIn("alpha_diversity_csv", paths)
        self.assertEqual(Path(paths["chimera_stats_csv"]).read_bytes(), FROZEN_CHIMERA.read_bytes())
        self.assertEqual(Path(paths["alpha_diversity_csv"]).read_bytes(), FROZEN_ALPHA.read_bytes())

    def test_abundance_matrix_reconstructs_frozen_otu_matrix_from_disk(self):
        from ui.callbacks_analysis import _resolve_output_dir
        from ui.callbacks_results import _find_abundance_matrix

        run_output_dir = _resolve_output_dir(None, self.project_data)
        self.assertIsNotNone(run_output_dir)
        matrix = _find_abundance_matrix({}, run_output_dir)
        self.assertIsNotNone(matrix)
        self.assertEqual(Path(matrix).read_bytes(), FROZEN_MATRIX.read_bytes())


class AnalysisRestoreWiringTests(unittest.TestCase):
    def test_selector_free_render_callbacks_depend_on_store_project(self):
        from app import create_app

        app = create_app()
        client = app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()
        targets = {
            "analysis-chimera.children",
            "analysis-alpha-table.children",
            "analysis-otu-table.children",
            "analysis-sample-table.children",
            "analysis-plots.children",
            "analysis-upset.children",
        }
        checked = set()
        for cb in deps:
            out = cb["output"]
            outs = [o.strip(".") for o in out.strip(".").split("...")] if out.startswith("..") else [out]
            for target in targets:
                if target in outs:
                    ids = {d["id"] for d in cb["inputs"]} | {d["id"] for d in cb["state"]}
                    self.assertIn("store-project", ids, f"{target}: render callback cannot see store-project")
                    checked.add(target)
        self.assertEqual(checked, targets, f"target callbacks not found: {targets - checked}")


if __name__ == "__main__":
    unittest.main()
