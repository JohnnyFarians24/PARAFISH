"""No Setup field may silently outrank the configured default it stands for.

Measured on a real COI run: every UI run clustered with ``min_unique_size = 2`` while
the project's documented default — decided deliberately, and implemented as a
configurable parameter so a rare parasite's singleton survives — was 1. Two sources
defined one value and the one nobody was reading won. The same shape as the phylo
guard and the two ``latest_run_dir``.

The cost was not theoretical. On the pre-cutadapt run a floor of 2 discarded 45 365 of
179 061 filtered reads (25.3%) as singletons — reads that were unique only because of
the degenerate positions in the primer, not because of biology.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

from core.config import MIN_UNIQUE_SIZE
from core.models import AnalysisBranch, RunConfig

PAGE = Path(__file__).resolve().parent.parent / "ui" / "pages" / "page_setup.py"

# Setup control id -> (RunConfig group, key) it sets.
FIELD_TO_PARAM = {
    "setup-trim-leading": ("trimmomatic_params", "leading"),
    "setup-trim-trailing": ("trimmomatic_params", "trailing"),
    "setup-trim-minlen": ("trimmomatic_params", "minlen"),
    "setup-pear-min-overlap": ("pear_params", "min_overlap"),
    "setup-pear-min-length": ("pear_params", "min_length"),
    "setup-pear-pvalue": ("pear_params", "p_value"),
    "setup-vsearch-maxee": ("vsearch_params", "fastq_maxee"),
    "setup-vsearch-minlen": ("vsearch_params", "fastq_minlen"),
    "setup-vsearch-sintax": ("vsearch_params", "sintax_cutoff"),
    "setup-vsearch-minsize": ("vsearch_params", "min_unique_size"),
    "setup-threads": ("vsearch_params", "threads"),
}

# Divergences that are a REPORTED, open scientific decision rather than a bug. Anything
# not listed here must agree. Keep this list empty if you can; each entry is a run whose
# parameters depend on whether it was started from the UI or from scripts/run_pipeline.py.
KNOWN_DIVERGENCES: set[str] = set()


def _page_defaults() -> dict[str, str]:
    """``value=`` of every numeric dcc.Input on the Setup page, as written."""
    source = PAGE.read_text(encoding="utf-8")
    found = {}
    for match in re.finditer(r'dcc\.Input\(\s*id="(setup-[a-z0-9-]+)"[^)]*?\bvalue=([^,)]+)', source):
        found[match.group(1)] = match.group(2).strip()
    return found


class SetupDefaultsSingleSourceTests(unittest.TestCase):
    def setUp(self):
        self.page = _page_defaults()
        self.config = RunConfig(run_id="t", samples=[], branch=AnalysisBranch.VSEARCH)

    def test_min_unique_size_is_not_a_literal_on_the_page(self):
        """The field must READ the constant. A number typed here is a second source,
        and it is the one that reaches the run."""
        self.assertEqual(self.page.get("setup-vsearch-minsize"), "MIN_UNIQUE_SIZE")

    def test_the_effective_default_keeps_singletons(self):
        self.assertEqual(MIN_UNIQUE_SIZE, 1)
        self.assertEqual(self.config.vsearch_params["min_unique_size"], 1)

    def test_no_setup_field_silently_disagrees_with_the_run_config(self):
        disagreements = []
        for field, (group, key) in FIELD_TO_PARAM.items():
            if field in KNOWN_DIVERGENCES:
                continue
            written = self.page.get(field)
            if written is None:
                continue
            expected = getattr(self.config, group)[key]
            try:
                same = float(written) == float(expected)
            except ValueError:
                # A named constant rather than a literal — that IS the single source.
                same = True
            if not same:
                disagreements.append(f"{field}={written} but {group}.{key}={expected}")
        self.assertEqual(disagreements, [], "; ".join(disagreements))

    def test_trim_minlen_is_aligned_at_50(self):
        """Resolved 2026-07-24, DOWN rather than up: 200 would have discarded 8.00% of
        merges on the COI run, 92.82% of them valid amplicons, to remove the 1.27% with
        frail overlap. The quantity at risk is PEAR's min_overlap, not read length."""
        self.assertEqual(self.page.get("setup-trim-minlen"), "50")
        self.assertEqual(self.config.trimmomatic_params["minlen"], 50)

    def test_no_divergence_is_currently_tolerated(self):
        """KNOWN_DIVERGENCES is empty and should stay that way. Each entry would be a
        run whose parameters depend on whether it started from the UI or the CLI."""
        self.assertEqual(KNOWN_DIVERGENCES, set())


def _find(node, wanted):
    """The component with id ``wanted`` in a Dash layout tree, or ``None``."""
    if getattr(node, "id", None) == wanted:
        return node
    children = getattr(node, "children", None)
    for child in children if isinstance(children, (list, tuple)) else [children]:
        if child is not None and not isinstance(child, str):
            found = _find(child, wanted)
            if found is not None:
                return found
    return None


class NonNumericSetupDefaultsTests(unittest.TestCase):
    """The controls the regex above cannot see: a text field, a dropdown, an expression.

    ``_page_defaults`` reads the ``value=`` of numeric ``dcc.Input`` fields as written, so
    three controls sat outside the net and each disagreed with ``RunConfig``: the sliding
    window is a TEXT field ("4:15" against 4:20), the BOLD mode is a dropdown (1 against
    2), and the bootstrap count was simply not mapped (100 against 10). A UI run and a
    ``scripts/run_pipeline.py`` run therefore used different parameters. These read the
    rendered layout instead, so the kind of control no longer matters.
    """

    def setUp(self):
        from ui.pages.page_setup import layout

        self.layout = layout()
        self.config = RunConfig(run_id="t", samples=[], branch=AnalysisBranch.VSEARCH)
        self.source = PAGE.read_text(encoding="utf-8")

    def _value_expression(self, component_id: str) -> str:
        """The page source from this component's id up to the next component's id."""
        start = self.source.index(f'id="{component_id}"')
        end = self.source.find('id="', start + 4)
        return self.source[start: end if end != -1 else len(self.source)]

    def test_sliding_window_text_field_matches_the_run_config(self):
        params = self.config.trimmomatic_params
        expected = f"{params['window_size']}:{params['window_quality']}"
        self.assertEqual(_find(self.layout, "setup-trim-slidingwindow").value, expected)

    def test_sliding_window_tooltip_does_not_advertise_another_default(self):
        params = self.config.trimmomatic_params
        expected = f"Default: {params['window_size']}:{params['window_quality']}."
        self.assertNotIn("Default: 4:15", self.source)

        titles = []

        def collect(node):
            title = getattr(node, "title", None)
            if isinstance(title, str):
                titles.append(title)
            children = getattr(node, "children", None)
            for child in children if isinstance(children, (list, tuple)) else [children]:
                if child is not None and not isinstance(child, str):
                    collect(child)

        collect(self.layout)
        window_tips = [title for title in titles if title.startswith("Window:quality")]
        self.assertEqual(len(window_tips), 1)
        self.assertTrue(window_tips[0].endswith(expected), window_tips[0])

    def test_raxml_bootstraps_match_the_run_config(self):
        self.assertEqual(_find(self.layout, "setup-raxml-bootstraps").value, self.config.phylo_params["bs_trees"])

    def test_bold_mode_matches_the_run_config(self):
        self.assertEqual(_find(self.layout, "setup-bold-mode").value, self.config.bold_params["mode"])

    def test_these_fields_read_the_run_config_rather_than_repeat_a_literal(self):
        """Equal today is not the same as single-sourced: a literal that happens to match
        drifts the next time RunConfig changes, which is how these three got here."""
        # The sliding window goes through one named constant, itself built from
        # _RUN_DEFAULTS, because the same string also feeds the field's tooltip.
        self.assertIn("_DEFAULT_SLIDINGWINDOW = f\"{_RUN_DEFAULTS.trimmomatic_params", self.source)
        reads = {
            "setup-trim-slidingwindow": "value=_DEFAULT_SLIDINGWINDOW",
            "setup-raxml-bootstraps": "value=_RUN_DEFAULTS.phylo_params",
            "setup-bold-mode": "value=_RUN_DEFAULTS.bold_params",
        }
        for component_id, expression in reads.items():
            with self.subTest(component_id=component_id):
                self.assertIn(expression, self._value_expression(component_id))

    def test_the_defaults_these_fields_stand_for_are_unchanged(self):
        """The alignment moved the PAGE to RunConfig, never the other way round."""
        self.assertEqual(self.config.trimmomatic_params["window_size"], 4)
        self.assertEqual(self.config.trimmomatic_params["window_quality"], 20)
        self.assertEqual(self.config.phylo_params["bs_trees"], 10)
        self.assertEqual(self.config.bold_params["mode"], 2)

    def test_ncbi_max_hits_has_no_run_config_default_to_align_to(self):
        """Reported, not fixed: RunConfig defines no NCBI max-hits default, so there is no
        single source to point the Setup field at. The field (5), the full-pipeline
        fallback (20) and NCBIConfig.max_target_seqs (5) remain three separate numbers.
        If a default is ever added to RunConfig this fails, and the field must then read it."""
        self.assertNotIn("ncbi_maxhits", self.config.vsearch_params)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
