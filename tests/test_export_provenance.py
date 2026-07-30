"""The run ZIP the supervisor downloads must carry the RECIPE, not only the results.

Measured (2026-07-30): the old premise ("docs promise PROVENANCE.txt, no run emits it") was
stale — the README no longer promises one and the machine-readable recipe now lives in
progress.json (``_run_provenance``). But the real gap the measurement found is here:
``build_run_zip`` bundled outputs and metrics and NOT the recipe, so a supervisor opening the
ZIP had the tables without which sources / marker / parameters made them. Same principle as
sintax_status — the explanation travels with the deliverable.

``render_provenance`` renders the run_config progress.json already carries (a legible view of
the single source, never a second collection), plus tool versions. These tests pin that the
file enters the ZIP AND carries the recipe axes — not merely that it exists — and that a
missing recipe or missing versions is declared, never hidden behind a complete-looking file.
"""

from __future__ import annotations

import unittest
import zipfile

from core.export import build_run_zip, render_provenance

_CONFIG = {
    "db_sources": ["silva", "pr2", "ncbi"],
    "markers": ["COI"],
    "analysis_method": "vsearch",
    "branch": "otu",
    "min_unique_size": 1,
    "primers": {"fwd": "GGWACWGGWTGAACWGTWTAYCCYCC", "rev": "TAAACTTCAGGGTGACCAAAAAATCA"},
    "ncbi_mode": "local",
}


class RenderProvenanceTests(unittest.TestCase):
    def test_recipe_axes_are_human_readable_not_json(self):
        text = render_provenance(_CONFIG, [("vsearch", "2.30.0")])
        self.assertIn("SILVA, PR2, NCBI", text)  # sources, upper-cased for a human
        self.assertIn("COI", text)
        self.assertIn("min_unique_size", text)
        self.assertIn("2.30.0", text)
        self.assertNotIn("{", text)  # not the raw JSON dump — that adds nothing over progress.json

    def test_missing_versions_are_declared_as_a_known_gap(self):
        # Parameters without versions is half a recipe; the file must say so, not look complete.
        text = render_provenance(_CONFIG, None)
        self.assertIn("lacuna conhecida", text)

    def test_a_run_without_a_recorded_recipe_says_so(self):
        text = render_provenance({}, [("vsearch", "2.30.0")])
        self.assertIn("não foi registada", text)


class ZipProvenanceTests(unittest.TestCase):
    def test_provenance_enters_the_zip_carrying_the_recipe(self):
        text = render_provenance(_CONFIG, [("vsearch", "2.30.0"), ("blastn", "2.17.0")])
        buf = build_run_zip("run1", {}, provenance_text=text)
        with zipfile.ZipFile(buf) as zf:
            self.assertIn("run1/PROVENANCE.txt", zf.namelist())
            content = zf.read("run1/PROVENANCE.txt").decode("utf-8")
        # Not merely present — it carries the recipe axes a diff must be read under.
        self.assertIn("SILVA, PR2, NCBI", content)
        self.assertIn("COI", content)
        self.assertIn("blastn", content)
        self.assertIn("2.17.0", content)

    def test_no_provenance_text_leaves_the_zip_as_it_was(self):
        buf = build_run_zip("run1", {})
        with zipfile.ZipFile(buf) as zf:
            self.assertNotIn("run1/PROVENANCE.txt", zf.namelist())


if __name__ == "__main__":
    unittest.main()
