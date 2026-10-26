"""Literal pins for the scientific parameters the paper's Methods section cites.

The 2026-07-24 mutation sweep found these values were NOT verifiably fixed: changing the
clustering identity or the singleton floor left the whole suite green, so any edit — by
anyone — would pass unnoticed and the Methods section would silently stop describing what
the code does. That is the gravest form of the pattern this project keeps meeting
(coherence is not correctness), because the reproducibility the paper claims rests on
exactly these numbers.

Two rules, both deliberate and both learned the hard way:

1. **The asserted value is a LITERAL, never imported from the constant it protects.**
   `assert x == 0.97`, never `assert x == SOME_CONSTANT`. A test that reads the constant
   it means to pin can only ever prove the code agrees with itself — which is precisely
   what made the first MIN_REFERENCE_LENGTH boundary test useless (it derived its inputs
   from the constant, so moving it to 100 kept the test green). The literal is the point.

2. **Each pin states WHY the value is what it is.** A pin does not exist to prevent
   change; it exists to force change to be deliberate and visible in the diff. A mute pin
   guards against the accident; a pin with the reason makes whoever changes it confront
   the why, not just the number. The reasoning in full lives in the DEVLOG entries these
   docstrings point at.

Every source-of-truth is pinned, including the UI control — the clustering-identity
Slider escaped the single-source divergence guard (which scans only dcc.Input), so it is
read here by walking the resolved layout tree instead.
"""

from __future__ import annotations

import inspect
import unittest

from analysis.vsearch_runner import run_vsearch_otus
from core.config import MIN_UNIQUE_SIZE
from core.models import AnalysisBranch, RunConfig
from ui.pages.page_setup import layout


def _control_value(control_id: str):
    """The ``value`` of a Setup control, read from the resolved layout tree.

    Not a regex over the source: the identity threshold is a dcc.Slider, and the
    Input-only regex in test_setup_defaults_single_source is exactly how its divergence
    went unnoticed. Walking the component tree sees every control type.
    """
    def find(node):
        if getattr(node, "id", None) == control_id:
            return node
        children = getattr(node, "children", None)
        if children is None:
            return None
        if not isinstance(children, (list, tuple)):
            children = [children]
        for child in children:
            hit = find(child)
            if hit is not None:
                return hit
        return None

    node = find(layout())
    if node is None:
        raise AssertionError(f"control {control_id!r} not found in the Setup layout")
    return node.value


class ScientificParameterPins(unittest.TestCase):
    def setUp(self):
        self.config = RunConfig(run_id="pin", samples=[], branch=AnalysisBranch.VSEARCH)

    def test_clustering_identity_is_0_97(self):
        """0.97 = species-level OTU threshold: two sequences within 3% cluster into one
        OTU. This is the single value the whole feature count rests on — it decides
        whether two OTUs collapse (see Task 5, OTU1/OTU2). It had three literal sources
        and zero test; the sweep changed it to 0.90 and nothing went red. DEVLOG
        2026-07-24 'DIAGNOSTIC SESSION'."""
        self.assertEqual(self.config.vsearch_params["id"], 0.97)
        self.assertEqual(_control_value("setup-vsearch-id"), 0.97)
        self.assertEqual(
            inspect.signature(run_vsearch_otus).parameters["identity"].default, 0.97
        )

    def test_sintax_cutoff_is_0_8(self):
        """0.8 = the minimum SINTAX bootstrap confidence to accept a taxonomic
        assignment. Below it a taxon is reported one rank up rather than guessed. Only
        the divergence guard caught a change, and a coordinated edit passed it green.
        DEVLOG 2026-07-24."""
        self.assertEqual(self.config.vsearch_params["sintax_cutoff"], 0.8)
        # The control that takes effect: the Run Pipeline callback reads this slider.
        self.assertEqual(_control_value("setup-sintax-cutoff"), 0.8)
        self.assertEqual(_control_value("setup-vsearch-sintax"), 0.8)

    def test_pear_min_overlap_is_10(self):
        """10 bp is PEAR's minimum overlap to merge a read pair. Measured too low on the
        COI run — it admits 11 bp junctions — but raising it re-baselines the golden, so
        it is an OPEN decision (DEVLOG 2026-07-24 trim-minlen), pinned here at its
        current value so any change is deliberate and visible."""
        self.assertEqual(self.config.pear_params["min_overlap"], 10)
        self.assertEqual(_control_value("setup-pear-min-overlap"), 10)

    def test_trimmomatic_minlen_is_50(self):
        """50 bp, aligned DOWN from 200 on 2026-07-24: 200 would have discarded 8.00% of
        merges, 92.82% of them valid amplicons, to remove the 1.27% with frail overlap.
        The quantity at risk is min_overlap above, not read length. DEVLOG 2026-07-24
        'trim-minlen resolved DOWN at 50'. (Also pinned in
        test_setup_defaults_single_source; kept here so all five live in one place.)"""
        self.assertEqual(self.config.trimmomatic_params["minlen"], 50)
        self.assertEqual(_control_value("setup-trim-minlen"), 50)

    def test_min_unique_size_is_1(self):
        """1 keeps singletons: a rare parasite is present as few reads, and a floor of 2
        deletes exactly that profile before clustering sees it. Measured identical output
        to 2 on the reference run, so 1 costs nothing and 2 risks the platform's whole
        purpose. DEVLOG 2026-07-24 'min_unique_size'. (Also pinned as
        keeps_singletons; kept here for the single-place invariant.)"""
        self.assertEqual(MIN_UNIQUE_SIZE, 1)
        self.assertEqual(self.config.vsearch_params["min_unique_size"], 1)
        self.assertEqual(_control_value("setup-vsearch-minsize"), 1)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
