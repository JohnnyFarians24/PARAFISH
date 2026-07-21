"""The Analysis residual: selector-free callbacks must render the RESOLVED set, not "otu".

Eixo A harvest, step 1, proven against the fixture_asv golden frozen this cycle. The six
selector-free Analysis callbacks had their Eixo B restore done and green, but two resolvers
sat behind them and only ONE carried the "otu" residual:

  * Resolver A — ``_analysis_metric_paths`` — hardcoded "otu": on an ASV run it returned
    ``{}`` (metrics/otu absent), blanking alpha / plots / upset. THIS is the residual fixed
    here: resolve the feature set from what the run produced, then read metrics/<set>/.
  * Resolver B — ``_find_abundance_matrix`` — already set-adaptive by its fallback order
    (otu_matrix → otu_table → asv_abundance), so otu_table / sample_table already rendered
    ASV. Locked below, not fixed — it was never otu-hardcoded.

``chimera_stats`` is OTU-only (DADA2 removes chimeras internally); on an ASV run resolver A
correctly yields no ``chimera_stats_csv``, so ``_render_chimera_stats`` shows empty — right,
not a regression.

Resolver A seen to fail before the fix: ``{}`` on the ASV run (the "otu" assumption).
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GA = ROOT / "golden" / "fixture_asv" / "analysis" / "analysis"
FROZEN_ASV_ALPHA = GA / "metrics" / "asv" / "alpha_diversity.csv"
FROZEN_ASV_RANK = GA / "metrics" / "asv" / "rank_abundance.csv"
FROZEN_ASV_ABUNDANCE = GA / "dada2" / "asv_abundance.csv"
SRC_ANALYSIS = ROOT / "runs" / "_fixture_asv" / "current" / "analysis"


def _make_run(runs_root: Path, project_slug: str, run_id: str) -> Path:
    run_dir = runs_root / project_slug / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SRC_ANALYSIS, run_dir / "analysis", dirs_exist_ok=True)
    (run_dir / "progress.json").write_text(json.dumps({"run_id": run_id, "done": True}), encoding="utf-8")
    return run_dir


class AnalysisAsvResidualTests(unittest.TestCase):
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

    def test_resolver_a_reconstructs_asv_metrics_not_blank(self):
        from ui.callbacks_analysis import _analysis_metric_paths

        paths = _analysis_metric_paths(None, None, self.project_data)
        self.assertIn("alpha_diversity_csv", paths)
        self.assertIn("rank_abundance_csv", paths)
        self.assertEqual(Path(paths["alpha_diversity_csv"]).read_bytes(), FROZEN_ASV_ALPHA.read_bytes())
        self.assertEqual(Path(paths["rank_abundance_csv"]).read_bytes(), FROZEN_ASV_RANK.read_bytes())
        # chimera_stats is OTU-only; an ASV run must not surface one.
        self.assertNotIn("chimera_stats_csv", paths)

    def test_resolver_b_abundance_matrix_already_resolves_asv(self):
        # Lock, not a fix: _find_abundance_matrix was set-adaptive before this change.
        from ui.callbacks_analysis import _resolve_output_dir
        from ui.callbacks_results import _find_abundance_matrix

        run_output_dir = _resolve_output_dir(None, self.project_data)
        matrix = _find_abundance_matrix({}, run_output_dir)
        self.assertIsNotNone(matrix)
        self.assertEqual(Path(matrix).read_bytes(), FROZEN_ASV_ABUNDANCE.read_bytes())


if __name__ == "__main__":
    unittest.main()
