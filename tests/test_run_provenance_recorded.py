"""Provenance is recorded WHEN THE RUN STARTS, not reconstructed when someone exports it.

The run_config block of progress.json carried seven keys and none of: the parameters each
stage ran with, the tool versions, the reference files, or the code commit. PROVENANCE.txt
then resolved tool versions at EXPORT time — so a ZIP exported after an upgrade reported a
version that never touched the data.

The change is additive. Four new keys join the block (``parameters``, ``tool_versions``,
``references``, ``parafish_commit``); the seven existing keys keep their names and values.
PROVENANCE.txt prefers the recorded values and falls back to the old behaviour for a run
that predates them — pinned here against the old rendering, byte for byte.

Nothing in the recording may fail a run: an unresolvable tool is ``"unresolved"``, an
unhashable reference is ``"unavailable"``, a missing git is ``"unknown"``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent.parent

EXISTING_KEYS = ["db_sources", "min_unique_size", "analysis_method", "branch", "markers", "primers", "ncbi_mode"]
NEW_KEYS = ["parameters", "tool_versions", "references", "parafish_commit"]

# What a run recorded before this change — the real block of runs/projeto_1/20260723T072541Z.
OLD_CONFIG = {
    "db_sources": ["silva", "pr2", "bold", "ncbi"],
    "min_unique_size": 2,
    "analysis_method": "vsearch",
    "branch": "vsearch",
    "markers": ["OTHER"],
    "primers": {"fwd": "", "rev": ""},
    "ncbi_mode": "remote",
}

# render_provenance(OLD_CONFIG, [("vsearch", "2.30.0")]) as it was before this change.
OLD_RENDERING = (
    "PARAFISH — run provenance\n"
    "==============================================\n"
    "\n"
    "Recipe\n"
    "----------------------------------------------\n"
    "  Reference sources : SILVA, PR2, BOLD, NCBI\n"
    "  Marker(s)         : OTHER\n"
    "  Analysis method   : vsearch\n"
    "  Branch            : vsearch\n"
    "  min_unique_size   : 2\n"
    "  Primers  fwd      : —\n"
    "  Primers  rev      : —\n"
    "  NCBI              : remote\n"
    "\n"
    "Tool versions\n"
    "----------------------------------------------\n"
    "  vsearch      2.30.0\n"
)


def _load_run_pipeline():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("run_pipeline", ROOT / "scripts" / "run_pipeline.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fixture_config(tmp: Path, db_sources: list, stages: list):
    run_pipeline = _load_run_pipeline()
    from ui.runtime_state import run_config_from_store

    project_data, samples_data, run_config_data = run_pipeline.build_stores(
        run_pipeline.FIXTURE_INPUTS, tmp, stages, db_sources, "testrun", {"threads": 1}
    )
    run_config = run_config_from_store(project_data, samples_data, run_config_data)
    run_config.run_id = "testrun"
    run_config.output_dir = tmp
    return run_config, run_config_data


class NewRunRecordsItsProvenanceTests(unittest.TestCase):
    """Through the real pipeline entry point, every stage disabled so no tool is needed."""

    @classmethod
    def setUpClass(cls):
        from ui.callbacks_progress import _run_pipeline_thread

        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        run_config, run_config_data = _fixture_config(tmp, ["silva", "pr2"], [])
        _run_pipeline_thread(run_config, run_config_data, tmp / "progress.json")
        cls.block = json.loads((tmp / "progress.json").read_text())["run_config"]

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_the_new_keys_are_present(self):
        for key in NEW_KEYS:
            with self.subTest(key=key):
                self.assertIn(key, self.block)

    def test_the_existing_keys_are_still_there_first_and_unchanged(self):
        self.assertEqual(list(self.block)[: len(EXISTING_KEYS)], EXISTING_KEYS)
        self.assertEqual(self.block["db_sources"], ["silva", "pr2"])
        self.assertEqual(self.block["markers"], ["18S"])
        self.assertEqual(self.block["min_unique_size"], 1)
        self.assertEqual(self.block["ncbi_mode"] in ("local", "remote"), True)

    def test_parameters_carry_every_stage_as_it_will_run(self):
        parameters = self.block["parameters"]
        for stage in ("trimmomatic", "pear", "vsearch", "dada2", "databases", "phylo"):
            with self.subTest(stage=stage):
                self.assertIn(stage, parameters)
        # The pinned scientific parameters, read back from what the run recorded.
        self.assertEqual(parameters["vsearch"]["id"], 0.97)
        self.assertEqual(parameters["vsearch"]["min_unique_size"], 1)
        self.assertEqual(parameters["pear"]["min_overlap"], 10)
        self.assertEqual(parameters["trimmomatic"]["minlen"], 50)
        self.assertEqual(parameters["databases"]["sintax_cutoff"], 0.8)
        self.assertEqual(parameters["databases"]["sintax_randseed"], 42)
        self.assertEqual(parameters["phylo"]["tree"]["seed"], 42)
        self.assertEqual(parameters["databases"]["consensus_method"], "lca")

    def test_tool_versions_name_each_configured_tool(self):
        versions = self.block["tool_versions"]
        self.assertIsInstance(versions, dict)
        self.assertIn("vsearch", versions)
        for tool, version in versions.items():
            with self.subTest(tool=tool):
                self.assertIsInstance(version, str)
                self.assertTrue(version)

    def test_no_reference_is_fingerprinted_when_the_databases_stage_is_off(self):
        # Every stage is disabled here, so no reference is USED — and hashing 800 MB of
        # SILVA for a run that never opens it would be cost without meaning.
        self.assertEqual(self.block["references"], [])

    def test_the_commit_is_a_sha_or_says_unknown(self):
        commit = self.block["parafish_commit"]
        self.assertTrue(commit == "unknown" or (len(commit) == 40 and set(commit) <= set("0123456789abcdef")), commit)

    def test_the_block_survives_a_json_round_trip(self):
        self.assertEqual(json.loads(json.dumps(self.block)), self.block)


class ReferencesAreRecordedForSourcesThatAreUsedTests(unittest.TestCase):
    def test_a_requested_local_reference_is_fingerprinted(self):
        import ui.callbacks_progress as progress

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            silva = tmp / "silva.fasta"
            silva.write_text(">a\nACGT\n", encoding="utf-8")
            run_config, run_config_data = _fixture_config(tmp / "run", ["silva", "ncbi"], ["databases"])
            fake_paths = SimpleNamespace(silva=silva, pr2=None, pr2_tax=None, pr2_dada2=None, bold_dir=None, blast=None)
            with patch.object(progress, "get_db_paths", return_value=fake_paths):
                block = progress._run_provenance(run_config, run_config_data)

        self.assertEqual([entry["source"] for entry in block["references"]], ["silva"])
        entry = block["references"][0]
        self.assertEqual(entry["file"], "silva.fasta")
        self.assertEqual(entry["size_bytes"], 8)
        self.assertEqual(entry["sha256"], hashlib.sha256(b">a\nACGT\n").hexdigest())
        # Remote NCBI has no local file to fingerprint; its settings are in parameters.
        self.assertEqual(block["parameters"]["databases"]["ncbi"]["task"], "megablast")

    def test_a_failing_version_resolver_does_not_fail_the_run(self):
        import ui.callbacks_progress as progress

        with tempfile.TemporaryDirectory() as d:
            run_config, run_config_data = _fixture_config(Path(d), ["silva"], [])
            with patch("ui.callbacks_results_summary._resolve_tool_versions", side_effect=RuntimeError("boom")):
                block = progress._run_provenance(run_config, run_config_data)
        self.assertEqual(block["tool_versions"], "unresolved")
        self.assertEqual(block["db_sources"], ["silva"])


class ReferenceFingerprintTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.ref = self.tmp / "ref.fasta"
        self.ref.write_bytes(b">a\nACGTACGT\n")
        self.sidecar = self.tmp / "ref.fasta.sha256.json"

    def tearDown(self):
        self._tmp.cleanup()

    def test_name_size_mtime_and_sha256(self):
        from core.export import reference_fingerprint

        result = reference_fingerprint(self.ref)
        self.assertEqual(result["file"], "ref.fasta")
        self.assertEqual(result["size_bytes"], 12)
        self.assertEqual(result["sha256"], hashlib.sha256(b">a\nACGTACGT\n").hexdigest())
        self.assertTrue(result["mtime"].endswith("Z") or "+" in result["mtime"])

    def test_the_hash_is_cached_beside_the_file_and_reused(self):
        from core.export import reference_fingerprint

        reference_fingerprint(self.ref)
        self.assertTrue(self.sidecar.exists())
        cached = json.loads(self.sidecar.read_text())
        # Prove the second call READS the cache: plant a sentinel no hash could produce.
        cached["sha256"] = "sentinel"
        self.sidecar.write_text(json.dumps(cached))
        self.assertEqual(reference_fingerprint(self.ref)["sha256"], "sentinel")

    def test_a_changed_file_invalidates_the_cache(self):
        from core.export import reference_fingerprint

        reference_fingerprint(self.ref)
        cached = json.loads(self.sidecar.read_text())
        cached["sha256"] = "sentinel"
        self.sidecar.write_text(json.dumps(cached))
        self.ref.write_bytes(b">a\nACGTACGTTT\n")  # size changes
        self.assertEqual(reference_fingerprint(self.ref)["sha256"], hashlib.sha256(b">a\nACGTACGTTT\n").hexdigest())

    def test_a_touched_file_of_the_same_size_invalidates_the_cache(self):
        from core.export import reference_fingerprint

        reference_fingerprint(self.ref)
        cached = json.loads(self.sidecar.read_text())
        cached["sha256"] = "sentinel"
        self.sidecar.write_text(json.dumps(cached))
        stat = self.ref.stat()
        os.utime(self.ref, ns=(stat.st_atime_ns, stat.st_mtime_ns + 5_000_000_000))
        self.assertNotEqual(reference_fingerprint(self.ref)["sha256"], "sentinel")

    def test_a_missing_file_is_unavailable_not_an_exception(self):
        from core.export import reference_fingerprint

        result = reference_fingerprint(self.tmp / "absent.fasta")
        self.assertEqual(result["sha256"], "unavailable")
        self.assertEqual(result["file"], "absent.fasta")

    def test_an_unwritable_cache_still_returns_the_hash(self):
        # The Docker mounts are read-only: the sidecar cannot be written there, and that
        # must cost a re-hash next run, not the hash itself.
        import core.export as export

        with patch.object(export, "_write_fingerprint_cache", side_effect=OSError("read-only")):
            result = export.reference_fingerprint(self.ref)
        self.assertEqual(result["sha256"], hashlib.sha256(b">a\nACGTACGT\n").hexdigest())
        self.assertFalse(self.sidecar.exists())

    def test_a_corrupt_cache_is_ignored(self):
        from core.export import reference_fingerprint

        self.sidecar.write_text("{not json")
        self.assertEqual(reference_fingerprint(self.ref)["sha256"], hashlib.sha256(b">a\nACGTACGT\n").hexdigest())


class GitCommitTests(unittest.TestCase):
    def test_reads_the_repository_head(self):
        from core.export import git_commit

        commit = git_commit()
        self.assertEqual(len(commit), 40)

    def test_unknown_when_git_cannot_answer(self):
        import core.export as export

        with patch.object(export.subprocess, "run", side_effect=FileNotFoundError("git")):
            self.assertEqual(export.git_commit(), "unknown")


class RenderProvenancePrefersWhatTheRunRecordedTests(unittest.TestCase):
    NEW_CONFIG = {
        **OLD_CONFIG,
        "parameters": {"vsearch": {"id": 0.97, "threads": 8}, "phylo": {"mode": "hybrid", "tree": {"seed": 42}}},
        "tool_versions": {"vsearch": "2.30.0", "fastqc": "unresolved"},
        "references": [
            {"source": "silva", "file": "SILVA.fasta", "size_bytes": 824019487, "mtime": "2026-03-10T19:16:00Z", "sha256": "abc123"},
        ],
        "parafish_commit": "0123456789abcdef0123456789abcdef01234567",
    }

    def test_an_old_run_renders_exactly_as_before(self):
        from core.export import render_provenance

        self.assertEqual(render_provenance(OLD_CONFIG, [("vsearch", "2.30.0")]), OLD_RENDERING)

    def test_an_old_run_without_versions_still_declares_the_gap(self):
        from core.export import render_provenance

        text = render_provenance(OLD_CONFIG, None)
        self.assertIn("known gap", text)
        self.assertNotIn("Parameters", text)
        self.assertNotIn("References", text)

    def test_recorded_versions_win_over_export_time_versions(self):
        from core.export import render_provenance

        # Exported after an upgrade: the binary now says 9.9.9, the run was made on 2.30.0.
        text = render_provenance(self.NEW_CONFIG, [("vsearch", "9.9.9")])
        self.assertIn("2.30.0", text)
        self.assertNotIn("9.9.9", text)
        self.assertIn("recorded when the run started", text)
        self.assertIn("unresolved", text)

    def test_recorded_parameters_references_and_commit_are_rendered(self):
        from core.export import render_provenance

        text = render_provenance(self.NEW_CONFIG, None)
        self.assertIn("0123456789abcdef0123456789abcdef01234567", text)
        self.assertIn("SILVA.fasta", text)
        self.assertIn("abc123", text)
        self.assertIn("[vsearch]", text)
        self.assertIn("id = 0.97", text)
        self.assertIn("tree.seed = 42", text)
        self.assertNotIn("{", text)  # still a person-readable file, not a JSON dump

    def test_the_recipe_block_is_unchanged_for_a_new_run(self):
        from core.export import render_provenance

        recipe = OLD_RENDERING.split("Tool versions")[0]
        self.assertTrue(render_provenance(self.NEW_CONFIG, None).startswith(recipe))


class RecordedBlockReachesTheProvenanceFileTests(unittest.TestCase):
    """What the run WRITES must be what PROVENANCE.txt READS — key for key.

    The rendering tests above feed ``render_provenance`` a hand-written block, so they
    prove the reader against the keys the test author typed, not against the keys the
    recorder emits. A key spelled one way by ``_run_provenance`` / ``reference_fingerprint``
    and another way by the renderer (``sha25`` for ``sha256``) would pass all of them and
    print "—" in every real run. These build the block with the real recorder and render
    it with the real renderer; nothing in between is typed by hand.
    """

    CONTENT = b">a\nACGTACGTACGT\n"

    @classmethod
    def setUpClass(cls):
        import ui.callbacks_progress as progress
        from core.export import render_provenance

        cls._tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls._tmp.name)
        silva = tmp / "silva_reference.fasta"
        silva.write_bytes(cls.CONTENT)
        run_config, run_config_data = _fixture_config(tmp / "run", ["silva"], ["databases"])
        fake_paths = SimpleNamespace(silva=silva, pr2=None, pr2_tax=None, pr2_dada2=None, bold_dir=None, blast=None)
        with patch.object(progress, "get_db_paths", return_value=fake_paths):
            cls.block = progress._run_provenance(run_config, run_config_data)
        # Export-time versions that must LOSE to the recorded ones.
        cls.text = render_provenance(cls.block, [("vsearch", "0.0.0-export-time")])

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_the_recorded_sha256_appears_in_the_file(self):
        digest = hashlib.sha256(self.CONTENT).hexdigest()
        self.assertEqual(self.block["references"][0]["sha256"], digest)
        self.assertIn(f"sha256 : {digest}", self.text)

    def test_the_recorded_file_name_and_source_appear(self):
        self.assertIn("SILVA    silva_reference.fasta", self.text)

    def test_the_recorded_size_appears(self):
        self.assertEqual(self.block["references"][0]["size_bytes"], len(self.CONTENT))
        self.assertIn(f"size   : {len(self.CONTENT)} bytes", self.text)

    def test_the_recorded_mtime_appears(self):
        mtime = self.block["references"][0]["mtime"]
        self.assertTrue(mtime)
        self.assertIn(f"mtime  : {mtime}", self.text)

    def test_the_recorded_commit_appears(self):
        commit = self.block["parafish_commit"]
        self.assertTrue(commit)
        self.assertIn(f"PARAFISH commit\n{'-' * 46}\n  {commit}", self.text)

    def test_the_recorded_tool_versions_appear_and_export_time_ones_do_not(self):
        versions = self.block["tool_versions"]
        self.assertIsInstance(versions, dict)
        for tool, version in versions.items():
            with self.subTest(tool=tool):
                self.assertIn(f"  {tool:12s} {version}", self.text)
        self.assertNotIn("0.0.0-export-time", self.text)

    def test_every_recorded_parameter_section_appears(self):
        for section in self.block["parameters"]:
            with self.subTest(section=section):
                self.assertIn(f"  [{section}]", self.text)
        self.assertIn("    id = 0.97", self.text)
        self.assertIn("    sintax_randseed = 42", self.text)

    def test_no_reference_field_renders_as_missing(self):
        # The renderer prints an em dash for a key it cannot find. For a block the recorder
        # just produced, any em dash in the References section is a writer/reader mismatch.
        section = self.text.split("References\n", 1)[1].split("\nParameters", 1)[0]
        self.assertNotIn("—", section)


if __name__ == "__main__":
    unittest.main()
