"""Tests for cross-sample comparison metrics, ZIP inclusion, and figure builders."""

from __future__ import annotations

import math
import tempfile
import zipfile
from pathlib import Path

import unittest

from analysis.sample_comparison import build_comparison
from core.export import build_run_zip


class SampleComparisonMetricsTests(unittest.TestCase):
    def _write(self, temp_path: Path, text: str, name: str = "otu_matrix.csv") -> Path:
        p = temp_path / name
        p.write_text(text, encoding="utf-8")
        return p

    def test_richness_shannon_shared_exclusive_csv(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            matrix = self._write(
                temp_path,
                "#OTU ID,S1,S2\nOTU1,6,0\nOTU2,3,3\nOTU3,0,5\n",
            )
            data = build_comparison(matrix)
            self.assertEqual(data["n_samples"], 2)
            self.assertEqual(data["n_features"], 3)
            self.assertEqual(data["samples"], ["S1", "S2"])
            self.assertEqual(data["shared"], 1)  # OTU2 present in both
            self.assertEqual(data["exclusive"], {"S1": 1, "S2": 1})  # OTU1 / OTU3
            by_sample = {d["sample"]: d for d in data["diversity"]}
            self.assertEqual(by_sample["S1"]["richness"], 2)
            self.assertEqual(by_sample["S2"]["richness"], 2)
            self.assertEqual(by_sample["S1"]["reads"], 9)
            # Shannon for S1 counts [6,3]: -(2/3 ln 2/3 + 1/3 ln 1/3)
            expected = -((6 / 9) * math.log(6 / 9) + (3 / 9) * math.log(3 / 9))
            self.assertAlmostEqual(by_sample["S1"]["shannon"], round(expected, 4), places=3)

    def test_tab_separated_matrix(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            matrix = self._write(
                temp_path,
                "#OTU ID\tSA\tSB\tSC\nOTU1\t1\t1\t1\nOTU2\t4\t0\t0\n",
                name="otu_table.tsv",
            )
            data = build_comparison(matrix)
            self.assertEqual(data["samples"], ["SA", "SB", "SC"])
            self.assertEqual(data["shared"], 1)  # OTU1 in all three
            self.assertEqual(data["exclusive"]["SA"], 1)  # OTU2 only in SA

    def test_missing_matrix_returns_error(self):
        data = build_comparison(Path("/nonexistent/matrix.csv"))
        self.assertIn("error", data)
        self.assertEqual(data["samples"], [])

    def test_heatmap_limited_to_top_features(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            rows = "\n".join(f"OTU{i},{i},{i}" for i in range(1, 51))
            matrix = self._write(temp_path, "#OTU ID,S1,S2\n" + rows + "\n")
            data = build_comparison(matrix, top_n=10)
            self.assertEqual(len(data["heatmap"]["features"]), 10)


class ZipInclusionTests(unittest.TestCase):
    def test_otu_matrix_included_in_zip(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            otu_matrix = temp_path / "otu_matrix.csv"
            otu_matrix.write_text("#OTU ID,S1\nOTU1,5\n", encoding="utf-8")
            stage_outputs = {
                "analysis": {"outputs": {"otu_matrix": str(otu_matrix)}},
            }
            buf = build_run_zip("run_x", stage_outputs)
            with zipfile.ZipFile(buf) as zf:
                names = zf.namelist()
            self.assertTrue(any(name.endswith("otu_matrix.csv") for name in names), names)


class FigureBuilderTests(unittest.TestCase):
    def test_comparison_figures_build(self):
        import plotly.graph_objects as go
        from ui.plots import comparison_figures

        with tempfile.TemporaryDirectory() as temp_dir:
            matrix = Path(temp_dir) / "otu_matrix.csv"
            matrix.write_text("#OTU ID,S1,S2\nOTU1,6,0\nOTU2,3,3\nOTU3,0,5\n", encoding="utf-8")
            figures = comparison_figures(build_comparison(matrix))
            self.assertGreaterEqual(len(figures), 3)
            self.assertTrue(all(isinstance(f, go.Figure) for f in figures))

    def test_comparison_figures_empty_on_error(self):
        from ui.plots import comparison_figures

        self.assertEqual(comparison_figures({"error": "x", "samples": []}), [])


if __name__ == "__main__":
    unittest.main()
