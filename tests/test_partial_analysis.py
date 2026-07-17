"""A half-successful analysis must still deliver its half.

The run that motivated these: VSEARCH clustered 7 OTUs perfectly, DADA2 was killed by
the OOM killer, and the pipeline then declared the whole Analysis stage a failure —
which made the Databases stage refuse to start ("análise anterior falhou"). A good OTU
set was thrown away because the *other* branch ran out of memory.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from analysis.analysis_metrics import feature_sets_from_outputs
from core.preflight import DADA2_BASELINE_MEMORY_BYTES, dada2_memory_shortfall, preflight_notes


class PartialAnalysisTests(unittest.TestCase):
    def test_a_dead_dada2_branch_still_leaves_the_otu_feature_set(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "otus.fasta").write_text(">OTU1;size=5\nACGT\n", encoding="utf-8")
            (root / "otu_matrix.csv").write_text("OTU_ID,S1,S2\nOTU1,3,2\n", encoding="utf-8")
            outputs = {
                "combined_fasta_vsearch": str(root / "otus.fasta"),
                "otu_matrix": str(root / "otu_matrix.csv"),
                # DADA2 was SIGKILLed: a result record with no files.
                "dada2": {"success": False, "files": {}},
            }
            sets = feature_sets_from_outputs(outputs)

        self.assertEqual([entry["key"] for entry in sets], ["otu"])
        # A non-empty feature-set list is what lets the Databases stage proceed.
        self.assertTrue(sets)

    def test_a_dead_vsearch_branch_still_leaves_the_asv_feature_set(self):
        """Symmetry: the same must hold the other way round."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "asv.fasta").write_text(">ASV1;size=4\nACGT\n", encoding="utf-8")
            (root / "asv_abundance.csv").write_text("ASV_ID,S1\nASV1,4\n", encoding="utf-8")
            outputs = {
                "combined_fasta_vsearch": "",
                "dada2": {
                    "success": True,
                    "files": {
                        "asv_sequences_fasta": str(root / "asv.fasta"),
                        "asv_abundance_csv": str(root / "asv_abundance.csv"),
                    },
                },
            }
            sets = feature_sets_from_outputs(outputs)

        self.assertEqual([entry["key"] for entry in sets], ["asv"])

    def test_both_branches_dead_leaves_nothing_to_classify(self):
        sets = feature_sets_from_outputs({"combined_fasta_vsearch": "", "dada2": {"success": False, "files": {}}})
        self.assertEqual(sets, [])


class Dada2MemoryGateTests(unittest.TestCase):
    """R holds ~747 MB resident just to load dada2 + Biostrings, before any data."""

    def test_refuses_when_the_libraries_alone_will_not_fit(self):
        with patch("core.preflight.available_memory_bytes", lambda: 200 * 1024 * 1024):
            message = dada2_memory_shortfall()
        self.assertIsNotNone(message)
        self.assertIn("OTUs", message)  # must say the OTU branch is unaffected

    def test_allows_when_there_is_headroom(self):
        with patch("core.preflight.available_memory_bytes", lambda: DADA2_BASELINE_MEMORY_BYTES + 1):
            self.assertIsNone(dada2_memory_shortfall())

    def test_the_gate_only_applies_to_methods_that_run_dada2(self):
        with patch("core.preflight.available_memory_bytes", lambda: 200 * 1024 * 1024), \
            patch("core.preflight.total_memory_bytes", lambda: 8 * 1024**3):
            vsearch_only = preflight_notes([], object(), 1, analysis_method="vsearch")
            with_dada2 = preflight_notes([], object(), 1, analysis_method="both")

        self.assertFalse(any("DADA2" in note for note in vsearch_only))
        self.assertTrue(any("DADA2" in note for note in with_dada2))


if __name__ == "__main__":
    unittest.main()
