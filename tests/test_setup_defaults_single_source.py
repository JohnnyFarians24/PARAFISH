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


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
