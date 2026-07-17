"""Run-ZIP export: a stray non-path string must never reach the filesystem."""

from __future__ import annotations

import tempfile
import unittest
import zipfile
import warnings
from pathlib import Path

from core.export import build_run_zip


class ExportRobustnessTests(unittest.TestCase):
    def _run_dir(self, root: Path) -> dict:
        (root / "otus.fasta").write_text(">OTU1;size=5\nACGT\n", encoding="utf-8")
        (root / "alpha.csv").write_text("sample,observed\nS1,3\n", encoding="utf-8")
        (root / "consensus.csv").write_text("ASV_ID,Consensus_Genus\nOTU1,Diplostomum\n", encoding="utf-8")
        return {
            "analysis": {
                "outputs": {
                    "combined_fasta_vsearch": str(root / "otus.fasta"),
                    "combined_fasta": str(root / "otus.fasta"),  # legacy alias, same file
                    "vsearch_pool": {"files": {"otus_fasta": str(root / "otus.fasta")}},
                    # The DADA2 branch was OOM-killed: a StepOutput whose error_message
                    # is a long sentence, sitting right where a path is expected.
                    "dada2": {
                        "step_name": "dada2",
                        "success": False,
                        "files": {},
                        "error_message": "DADA2 (ASVs) precisa de cerca de 0.9 GB livres " * 8,
                    },
                    "branch_errors": {"asv": "x" * 400},
                    "dropped_samples": [],
                },
                "metrics_files": {
                    "otu.alpha_diversity_csv": str(root / "alpha.csv"),
                    "alpha_diversity_csv": str(root / "alpha.csv"),  # alias, same file
                },
            },
            "databases": {
                "outputs": {
                    "taxonomy_by_set": {
                        "otu": {"files": {"taxonomy_consensus_csv": str(root / "consensus.csv")}},
                    },
                    "taxonomy": {"files": {"taxonomy_consensus_csv": str(root / "consensus.csv")}},
                },
            },
        }

    def test_a_long_error_message_does_not_crash_the_export(self):
        """``Path.exists()`` raises ENAMETOOLONG rather than returning False, so an
        error_message treated as a filename took the whole download down with
        "OSError: [Errno 36] File name too long"."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            buf = build_run_zip("run1", self._run_dir(root))  # must not raise
        names = zipfile.ZipFile(buf).namelist()
        self.assertTrue(names)
        self.assertFalse(any("precisa de cerca" in name for name in names))

    def test_no_duplicate_entries(self):
        """The same artifact is reachable by several routes; it must be stored once."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            with warnings.catch_warnings():
                warnings.simplefilter("error")  # zipfile warns on duplicate names
                buf = build_run_zip("run1", self._run_dir(root))
        names = zipfile.ZipFile(buf).namelist()
        self.assertEqual(len(names), len(set(names)))

    def test_each_feature_set_keeps_its_own_folder(self):
        """OTU and ASV tables share basenames; flattened, one would overwrite the other."""
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            (root / "c.csv").write_text("ASV_ID\nX\n", encoding="utf-8")
            outputs = {
                "databases": {
                    "outputs": {
                        "taxonomy_by_set": {
                            "otu": {"files": {"taxonomy_consensus_csv": str(root / "c.csv")}},
                            "asv": {"files": {"taxonomy_consensus_csv": str(root / "c.csv")}},
                        }
                    }
                }
            }
            names = zipfile.ZipFile(build_run_zip("run1", outputs)).namelist()

        self.assertIn("run1/databases/otu/c.csv", names)
        self.assertIn("run1/databases/asv/c.csv", names)

    def test_metric_tables_are_included(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            names = zipfile.ZipFile(build_run_zip("run1", self._run_dir(root))).namelist()
        self.assertIn("run1/analysis/metrics/otu/alpha.csv", names)


if __name__ == "__main__":
    unittest.main()
