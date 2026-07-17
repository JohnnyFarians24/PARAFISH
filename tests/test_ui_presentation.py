"""Presentation regressions: table contrast, small-n plots, id hygiene, layering."""

from __future__ import annotations

import unittest

import pandas as pd

from ui.components.table import STYLE_CELL, STYLE_HEADER, data_table
from ui.plots import SMALL_N, alpha_diversity_figure, qc_depth_box_figure


class TableContrastTests(unittest.TestCase):
    def test_every_themed_table_sets_an_explicit_text_colour(self):
        """QC and Preprocess passed only padding, so cells inherited Dash's near-black
        default on a dark background and the numbers were unreadable."""
        self.assertIn("color", STYLE_CELL)
        self.assertIn("backgroundColor", STYLE_CELL)
        self.assertIn("color", STYLE_HEADER)

    def test_data_table_keeps_caller_conditionals_on_top_of_striping(self):
        table = data_table(
            [{"a": 1}],
            [{"name": "A", "id": "a"}],
            conditional=[{"if": {"column_id": "a"}, "color": "red"}],
        )
        rules = table.style_data_conditional
        # Zebra striping first, caller's rules after, so the caller's colour wins.
        self.assertEqual(rules[0]["if"], {"row_index": "odd"})
        self.assertEqual(rules[-1]["color"], "red")


class SmallNPlotTests(unittest.TestCase):
    """A box plot through two points draws a quartile range that does not exist."""

    def _depth(self, n: int) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "read": ["R1"] * n,
                "sample": [f"S{i}" for i in range(n)],
                "total_reads": [10_000 + i for i in range(n)],
                "is_depth_outlier": [0] * n,
            }
        )

    def _alpha(self, n: int) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "sample": [f"S{i}" for i in range(n)],
                "observed": list(range(n)),
                "shannon": [1.0] * n,
                "pielou": [0.5] * n,
            }
        )

    def test_few_samples_render_as_labelled_bars(self):
        self.assertEqual(qc_depth_box_figure(self._depth(SMALL_N)).data[0].type, "bar")
        self.assertEqual(alpha_diversity_figure(self._alpha(SMALL_N)).data[0].type, "bar")

    def test_many_samples_still_render_as_a_distribution(self):
        self.assertEqual(qc_depth_box_figure(self._depth(SMALL_N + 4)).data[0].type, "box")
        self.assertEqual(alpha_diversity_figure(self._alpha(SMALL_N + 4)).data[0].type, "box")

    def test_a_single_sample_is_still_plotted(self):
        """The old box path bailed out below 2 samples and rendered nothing at all."""
        figure = qc_depth_box_figure(self._depth(1))
        self.assertIsNotNone(figure)
        self.assertEqual(figure.data[0].type, "bar")


class FeatureIdHygieneTests(unittest.TestCase):
    def test_raw_illumina_headers_never_reach_a_user_facing_table(self):
        from ui.callbacks_databases import _clean_feature_ids

        frame = pd.DataFrame(
            {
                "qseqid": ["OTU1;size=132828 A01426:1215:HW7YLDRX5:2:2101:6289:3458", "OTU2;size=99"],
                "pident": [99.1, 98.0],
            }
        )
        self.assertEqual(list(_clean_feature_ids(frame)["qseqid"]), ["OTU1", "OTU2"])

    def test_cleaning_does_not_mutate_the_caller_frame(self):
        from ui.callbacks_databases import _clean_feature_ids

        frame = pd.DataFrame({"ASV_ID": ["ASV1;size=5"]})
        _clean_feature_ids(frame)
        self.assertEqual(frame["ASV_ID"].iloc[0], "ASV1;size=5")


class LayeringTests(unittest.TestCase):
    def test_kpis_and_hero_plots_are_not_hidden_behind_a_disclosure(self):
        """A scientific page that opens empty cannot be judged at a glance.

        Layer-1 ids must be reachable without opening a <details>.
        """
        from dash import html

        from ui.pages.page_analysis import layout as analysis_layout
        from ui.pages.page_qc import layout as qc_layout

        def visible_ids(node, inside_details=False, found=None):
            found = set() if found is None else found
            is_details = isinstance(node, html.Details)
            if not inside_details and not is_details:
                node_id = getattr(node, "id", None)
                if isinstance(node_id, str):
                    found.add(node_id)
            children = getattr(node, "children", None)
            if children is not None:
                for child in children if isinstance(children, (list, tuple)) else [children]:
                    if hasattr(child, "children") or hasattr(child, "id"):
                        visible_ids(child, inside_details or is_details, found)
            return found

        qc = visible_ids(qc_layout())
        for required in ("qc-stat-total", "qc-stat-pass", "qc-sample-table", "qc-plot-hero"):
            self.assertIn(required, qc, f"{required} must stay on Layer 1")
        # The secondary curves are detail and belong behind the disclosure.
        self.assertNotIn("qc-plots", qc)

        analysis = visible_ids(analysis_layout())
        for required in ("analysis-sample-table", "analysis-alpha-table", "analysis-plots"):
            self.assertIn(required, analysis, f"{required} must stay on Layer 1")
        self.assertNotIn("analysis-chimera", analysis)


class SingleTableFactoryTests(unittest.TestCase):
    def test_no_page_instantiates_a_DataTable_directly(self):
        """T0's real acceptance criterion. Twenty tables were styled in four different
        places; two of those inherited Dash's default dark-grey-on-dark and were
        unreadable. One factory is the only way that cannot recur."""
        import subprocess

        out = subprocess.run(
            ["grep", "-rn", "dash_table.DataTable(", "--include=*.py", "ui/"],
            capture_output=True,
            text=True,
        )
        offenders = [
            line
            for line in out.stdout.splitlines()
            if line.strip() and not line.startswith("ui/components/table.py")
        ]
        self.assertEqual(offenders, [], f"DataTable instantiated outside the factory:\n" + "\n".join(offenders))

    def test_the_factory_meets_the_contrast_floor(self):
        """#e5e7eb on the card background is ~13:1. The inherited default was ~2:1."""
        from ui.components.table import STYLE_CELL

        def luminance(hex_colour):
            r, g, b = (int(hex_colour[i : i + 2], 16) / 255 for i in (1, 3, 5))
            channels = [(c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4) for c in (r, g, b)]
            return 0.2126 * channels[0] + 0.7152 * channels[1] + 0.0722 * channels[2]

        # Cells are transparent, so the effective background is the card's.
        text = luminance(STYLE_CELL["color"])
        card = luminance("#121923")
        ratio = (max(text, card) + 0.05) / (min(text, card) + 0.05)
        self.assertGreater(ratio, 4.5, f"contrast {ratio:.1f}:1 is below the WCAG AA floor")


if __name__ == "__main__":
    unittest.main()
