"""Stage 6 (Results): cross-stage read fate, sample dashboard, target report."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _run_tree(root: Path) -> Path:
    """A minimal on-disk run, at THE LAYOUT THE PIPELINE ACTUALLY WRITES.

    The feature-set level (``databases/otu/``, ``analysis/metrics/otu/``) is not
    decoration: this fixture used to write the flat layout — the one the *reader*
    expected and the pipeline never produced — so the test and the code shared one false
    belief and agreed with each other perfectly. 401 tests passed while the Results page
    was blank on a run that had classified 9/9 features.

    ``GoldenLayoutTests`` below pins this layout against the real frozen run, so this
    fixture cannot drift back into agreeing only with its author. See DEVLOG 2026-07-17.
    """
    _write(
        root / "preprocess" / "metrics" / "read_tracking.csv",
        "sample,mode,input_reads,trimmed_reads,forward_only,reverse_only,dropped_reads,"
        "assembled_reads,discarded_reads,not_assembled_reads,surviving_reads,pct_trimmed,pct_assembled,pct_surviving\n"
        "s1,PE,1000,980,10,5,5,900,80,0,900,98.0,91.8,90.0\n"
        "blank,PE,200,180,5,5,10,50,130,0,50,90.0,27.8,25.0\n",
    )
    _write(
        root / "analysis" / "metrics" / "otu" / "alpha_diversity.csv",
        "sample,reads,observed,chao1,shannon,simpson,pielou,hill_q0,hill_q1,hill_q2\n"
        "s1,50,3,3.0,0.8,0.6,0.7,3.0,2.2,1.9\n"
        "blank,50,1,1.0,0.0,0.0,0.0,1.0,1.0,1.0\n",
    )
    _write(
        root / "analysis" / "otu_matrix.csv",
        "#OTU ID,s1,blank\nf_para,80,0\nf_fish,20,0\nf_reagent,0,10\n",
    )
    _write(
        root / "databases" / "otu" / "taxonomy_consensus.csv",
        "ASV_ID,Consensus_Kingdom,Consensus_Phylum,Consensus_Class,Consensus_Order,Consensus_Family,Consensus_Genus,Consensus_Species\n"
        "f_para,Eukaryota,Platyhelminthes,Trematoda,Diplostomida,Diplostomidae,Diplostomum,Diplostomum spathaceum\n"
        "f_fish,Eukaryota,Chordata,Actinopteri,,,,\n"
        "f_reagent,Bacteria,Pseudomonadota,,,,Ralstonia,\n",
    )
    _write(
        root / "databases" / "otu" / "metrics" / "taxon_prevalence.csv",
        "taxon,n_detected,n_samples,prevalence,ci_low,ci_high,total_reads\n"
        "Diplostomum,1,1,1.0,0.2065,1.0,80\n",
    )
    _write(
        root / "databases" / "otu" / "source_assignments.csv",
        "ASV_ID,source,rank,taxon,identity,coverage,confidence,reliability\n"
        "f_para,NCBI,species,Diplostomum spathaceum,99.5,98,850,high\n"
        "f_reagent,PR2,genus,Ralstonia,,,0.6,low\n",
    )
    return root


class ReadFateSankeyTests(unittest.TestCase):
    def test_stages_chain_from_raw_to_classified(self):
        from analysis.results_summary import read_fate_sankey

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            sankey = read_fate_sankey(root, feature_reads=800, classified_reads=700)

        self.assertEqual(sankey["stages"], ["Raw", "Trimmed", "Merged", "Clustered", "Classified"])
        # Raw = 1000 + 200; trimmed = 980 + 180; merged = 900 + 50.
        self.assertEqual(sankey["values"][:3], [1200, 1160, 950])

    def test_a_downstream_flow_never_exceeds_its_upstream(self):
        """A merged count above the trimmed count would mean broken accounting."""
        from analysis.results_summary import read_fate_sankey

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            # Claim more clustered reads than exist upstream.
            sankey = read_fate_sankey(root, feature_reads=99999)

        merged_to_clustered = next(flow for flow in sankey["flows"] if flow["target"] == "Clustered")
        self.assertLessEqual(merged_to_clustered["value"], 950)

    def test_without_preprocess_metrics_there_is_no_diagram(self):
        from analysis.results_summary import read_fate_sankey

        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(read_fate_sankey(Path(d)), {})


class SampleDashboardTests(unittest.TestCase):
    def test_percent_classified_is_a_share_of_the_sample_not_the_feature(self):
        """counts is samples x features; summing the wrong axis gave 650% before."""
        from analysis.results_summary import sample_dashboard

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            rows = {row["sample"]: row for row in sample_dashboard(root, control_samples=["blank"])}

        # s1: 80 of 100 reads (f_para) are species-level; f_fish stops at class.
        self.assertEqual(rows["s1"]["pct_species"], 80.0)
        self.assertGreaterEqual(rows["s1"]["pct_species"], 0.0)
        self.assertLessEqual(rows["s1"]["pct_species"], 100.0)

    def test_controls_are_flagged_not_dropped(self):
        from analysis.results_summary import sample_dashboard

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            rows = {row["sample"]: row for row in sample_dashboard(root, control_samples=["blank"])}

        self.assertEqual(rows["blank"]["is_control"], 1)
        self.assertEqual(rows["s1"]["is_control"], 0)

    def test_target_taxa_are_counted_per_sample(self):
        from analysis.results_summary import sample_dashboard

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            rows = {row["sample"]: row for row in sample_dashboard(root, target_taxa=["Diplostomum"])}

        self.assertEqual(rows["s1"]["n_targets"], 1)
        self.assertEqual(rows["blank"]["n_targets"], 0)


class TargetTaxonReportTests(unittest.TestCase):
    def test_prevalence_and_reliability_come_from_the_taxonomy_tables(self):
        from analysis.results_summary import target_taxon_report

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            rows = {row["taxon"]: row for row in target_taxon_report(root, ["Diplostomum"])}

        self.assertEqual(rows["Diplostomum"]["prevalence"], 1.0)
        self.assertEqual(rows["Diplostomum"]["ci_low"], 0.2065)
        self.assertEqual(rows["Diplostomum"]["best_reliability"], "high")
        self.assertIn("s1", rows["Diplostomum"]["detected_in"])

    def test_an_absent_taxon_reports_not_detected(self):
        from analysis.results_summary import target_taxon_report

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            row = target_taxon_report(root, ["Anisakis"])[0]

        self.assertIsNone(row["prevalence"])
        self.assertEqual(row["detected_in"], "—")

    def test_controls_are_excluded_from_where_a_taxon_was_detected(self):
        from analysis.results_summary import target_taxon_report

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            # Ralstonia only appears in the blank; excluding it leaves no real detection.
            row = target_taxon_report(root, ["Ralstonia"], control_samples=["blank"])[0]

        self.assertEqual(row["detected_in"], "—")

    def test_empty_target_list_returns_nothing(self):
        from analysis.results_summary import target_taxon_report

        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(target_taxon_report(Path(d), []), [])


class ResultsFigureTests(unittest.TestCase):
    def test_sankey_figure_is_built_from_the_flows(self):
        from analysis.results_summary import read_fate_sankey
        from ui.plots import read_fate_sankey_figure

        with tempfile.TemporaryDirectory() as d:
            root = _run_tree(Path(d))
            figure = read_fate_sankey_figure(read_fate_sankey(root, feature_reads=800, classified_reads=700))

        self.assertEqual(figure.data[0].type, "sankey")

    def test_empty_inputs_render_nothing(self):
        from ui.plots import dashboard_richness_figure, read_fate_sankey_figure

        self.assertIsNone(read_fate_sankey_figure({}))
        self.assertIsNone(dashboard_richness_figure([]))

    def test_dashboard_colours_controls_apart(self):
        from ui.plots import dashboard_richness_figure

        rows = [
            {"sample": "s1", "is_control": 0, "reads": 900, "richness": 3},
            {"sample": "blank", "is_control": 1, "reads": 50, "richness": 1},
        ]
        colours = list(dashboard_richness_figure(rows).data[0].marker.color)
        self.assertEqual(colours[0], "#67a6ff")
        self.assertEqual(colours[1], "#f85149")


if __name__ == "__main__":
    unittest.main()


class SessionRestoreTests(unittest.TestCase):
    """Stage stores are per-session: a restarted app lost every result on disk."""

    def _project(self, name="proj"):
        return {"active_id": "p1", "projects": {"p1": {"id": "p1", "project_name": name}}}

    def _run_on_disk(self, runs_root: Path, name="proj", run_id="20260101T000000Z", done=True):
        run_dir = runs_root / name / run_id
        _write(
            run_dir / "progress.json",
            '{"run_id": "%s", "done": %s, "qc_record": {"status": "success"}, "preprocess_record": {"status": "success"},'
            '"analysis_record": {"status": "success"}, "db_record": {"status": "failure"},'
            '"phylo_record": {"status": "skipped"}}' % (run_id, "true" if done else "false"),
        )
        return run_dir

    def test_the_newest_run_of_the_active_project_is_restored(self):
        from unittest.mock import patch

        from ui import runtime_state

        with tempfile.TemporaryDirectory() as d:
            runs_root = Path(d)
            self._run_on_disk(runs_root, run_id="20260101T000000Z")
            newest = self._run_on_disk(runs_root, run_id="20260202T000000Z")
            with patch.object(runtime_state, "runs_root", lambda: runs_root):
                restored = runtime_state.restore_run_state(self._project())

        self.assertEqual(restored["run_config"]["run_id"], "20260202T000000Z")
        self.assertEqual(restored["run_config"]["output_dir"], str(newest))
        self.assertEqual(restored["qc"]["status"], "success")
        self.assertEqual(restored["databases"]["status"], "failure")

    def test_a_run_still_in_progress_restores_config_but_no_stage_records(self):
        """A run without ``done`` is still going: restoring its half-written stage
        records would unlock QC/Preprocess/Analysis mid-run, one at a time, which is
        exactly what the progress poller withholds until the whole run finishes."""
        from unittest.mock import patch

        from ui import runtime_state

        with tempfile.TemporaryDirectory() as d:
            runs_root = Path(d)
            running = self._run_on_disk(runs_root, run_id="20260303T000000Z", done=False)
            with patch.object(runtime_state, "runs_root", lambda: runs_root):
                restored = runtime_state.restore_run_state(self._project())

        # The run-config still comes back so the poller can resume and publish
        # everything together once the run reaches ``done``.
        self.assertEqual(restored["run_config"]["run_id"], "20260303T000000Z")
        self.assertEqual(restored["run_config"]["output_dir"], str(running))
        # But no stage record leaks through while the run is live — the tabs stay locked.
        for stage in ("qc", "preprocess", "analysis", "databases", "phylo"):
            self.assertEqual(restored[stage], {}, f"{stage} should be empty mid-run")

    def test_a_project_with_no_runs_restores_nothing(self):
        from unittest.mock import patch

        from ui import runtime_state

        with tempfile.TemporaryDirectory() as d:
            with patch.object(runtime_state, "runs_root", lambda: Path(d)):
                self.assertIsNone(runtime_state.restore_run_state(self._project()))
                self.assertIsNone(runtime_state.restore_run_state({}))


class SankeyCostTests(unittest.TestCase):
    """The Sankey used to re-scan clusters.uc on every render, which stalled the server."""

    def test_clustered_total_comes_from_the_persisted_chimera_table(self):
        from analysis.results_summary import clustered_read_total

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _write(
                root / "analysis" / "metrics" / "otu" / "chimera_stats.csv",
                "chimeric_variants,total_variants,pct_chimeric_variants,chimeric_reads,total_reads,pct_chimeric_reads,read_basis\n"
                "2,12,16.6,500,10000,5.0,cluster\n",
            )
            # Chimeric reads do not survive into the clustered feature set.
            self.assertEqual(clustered_read_total(root), 9500)

    def test_species_total_comes_from_the_persisted_resolution_table(self):
        from analysis.results_summary import species_classified_total

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _write(
                root / "databases" / "otu" / "metrics" / "taxonomic_resolution.csv",
                "rank,n_features,pct_features,n_reads,pct_reads\n"
                "genus,3,50.0,1000,20.0\nspecies,3,50.0,4000,80.0\n",
            )
            self.assertEqual(species_classified_total(root), 4000)

    def test_missing_tables_report_nothing_rather_than_recomputing(self):
        from analysis.results_summary import clustered_read_total, species_classified_total

        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(clustered_read_total(Path(d)))
            self.assertIsNone(species_classified_total(Path(d)))


class GoldenLayoutTests(unittest.TestCase):
    """Anchor the reader on the artifacts the PIPELINE wrote, not on a fixture.

    Everything above this class builds its own directory tree. That is fast and good for
    edge cases, but it can only ever prove the reader agrees with whoever wrote the test
    — which is precisely how the Results page came to read blank on a run that had
    classified 9/9 features while 401 tests passed (DEVLOG 2026-07-17). These take their
    tree from ``golden/``, which no reader's author typed.

    ``golden/fixture/`` is **committed** so these run on any clone: a contract that lives
    on one machine anchors nothing, and a skip everywhere else is the contract going
    unchecked while the suite reports green. ``golden/real/`` is not committed (size), so
    the real-data test skips off this machine — it is the stronger anchor, not the load
    -bearing one.
    """

    def _golden_run(self, root: Path, dataset: str = "fixture"):
        import shutil

        golden = Path(__file__).resolve().parent.parent / "golden" / dataset
        if not golden.is_dir():
            self.skipTest(f"golden/{dataset} is absent — the contract cannot be checked here")
        for stage in ("preprocess", "analysis", "taxonomy", "export"):
            if (golden / stage).is_dir():
                shutil.copytree(golden / stage, root, dirs_exist_ok=True)
        return root

    def test_reader_sees_the_frozen_fixture_run(self):
        """The numbers on screen must be the numbers the pipeline froze. Runs anywhere.

        Fixture-only limitation, stated rather than papered over: this run's taxonomy
        reaches GENUS and no further (``taxonomic_resolution`` says species = 0 reads),
        so it cannot exercise the species path. ``clustered`` and ``richness`` are the
        dimensions it does anchor — and they are the two that read blank on the bug.
        """
        from analysis.results_summary import clustered_read_total, sample_dashboard

        with tempfile.TemporaryDirectory() as d:
            root = self._golden_run(Path(d), "fixture")
            clustered = clustered_read_total(root)
            rows = {row["sample"]: row for row in sample_dashboard(root)}

        # golden/fixture: chimera_stats total_reads 8107, chimeric 0.
        self.assertEqual(clustered, 8107)
        # golden/fixture: alpha_diversity observed = 4 and 2.
        self.assertEqual(rows["aptl64b_02"]["richness"], 4)
        self.assertEqual(rows["aptl64f_01"]["richness"], 2)

    def test_reader_sees_the_frozen_real_run(self):
        """The stronger anchor: a run that reached species. Skips without golden/real."""
        from analysis.results_summary import clustered_read_total, sample_dashboard, species_classified_total

        with tempfile.TemporaryDirectory() as d:
            root = self._golden_run(Path(d), "real")
            clustered = clustered_read_total(root)
            classified = species_classified_total(root)
            rows = {row["sample"]: row for row in sample_dashboard(root)}

        # Read from golden/real/analysis + /taxonomy, not asserted from memory.
        self.assertEqual(clustered, 154223)
        self.assertIsNotNone(classified)
        self.assertGreater(classified, 0)
        self.assertTrue(rows, "the dashboard is empty on a run that classified 9/9 features")
        for sample, row in rows.items():
            self.assertGreater(row["richness"], 0, f"{sample} reports richness 0 on the frozen run")

    def test_the_fixture_layout_matches_the_pipeline_layout(self):
        """The hand-written fixture must not drift back to a layout nothing writes.

        This is the guard on the actual defect: the fixture put its tables where the
        READER looked, so it agreed with the bug. It now has to agree with golden — and
        with the COMMITTED golden, so the guard exists on every clone rather than here.
        """
        with tempfile.TemporaryDirectory() as d:
            pipeline = self._golden_run(Path(d) / "pipeline", "fixture")
            hand_written = _run_tree(Path(d) / "hand_written")

            def tables(root: Path) -> set:
                return {
                    p.relative_to(root).as_posix()
                    for p in root.rglob("*.csv")
                    if p.relative_to(root).as_posix().split("/")[0] in {"analysis", "databases"}
                }

            pipeline_tables = tables(pipeline)
            for path in tables(hand_written):
                name = path.rsplit("/", 1)[-1]
                twin = [candidate for candidate in pipeline_tables if candidate.endswith("/" + name) or candidate == name]
                if not twin:
                    continue  # golden does not freeze this table; nothing to anchor to
                self.assertIn(
                    path,
                    twin,
                    f"fixture writes {name} at '{path}', the pipeline writes it at {twin} — "
                    "a fixture at a path nothing produces can only confirm the reader's own belief",
                )
