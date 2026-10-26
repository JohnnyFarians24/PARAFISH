"""Databases page: per-database counts, and downloads that honour ticked or filtered rows.

The supervisors asked to see, after a run, which databases assigned the features, ranked
by what was found most, and to download only the features they pick. The two download
buttons used to ignore the table (ticks and the table's own filter): they exported the
whole table. The default, with nothing ticked and nothing filtered, must stay exactly what
it was: the hits table as filtered by the sliders, and the summary file byte for byte.
"""

from __future__ import annotations

import base64
import csv
import io
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from ui.table_export import export_positions, rows_to_export, source_assignment_summary, table_view_state


HITS = [
	# feature, reads, hit_rank, source, identity, coverage, taxon, rank
	("OTU1", 900, 1, "NCBI", 100.0, 100.0, "Gadus morhua", "species"),
	("OTU1", 900, 2, "SILVA", "", "", "Gadus morhua", "species"),
	("OTU1", 900, 3, "PR2", "", "", "Gadus", "genus"),
	("OTU2", 80, 1, "SILVA", "", "", "Eukaryota", "kingdom"),
	("OTU2", 80, 2, "NCBI", 98.5, 100.0, "Anisakis simplex", "species"),
	("OTU3", 20, 1, "PR2", "", "", "Kudoa", "genus"),
]
HITS_COLUMNS = ["feature_id", "reads", "hit_rank", "source", "identity", "coverage", "taxon", "rank"]


class RowsToExportTests(unittest.TestCase):
	def test_nothing_ticked_nothing_filtered_keeps_the_old_behaviour(self):
		self.assertIsNone(rows_to_export(4))
		self.assertIsNone(rows_to_export(4, [], None))
		# Sorting alone is a permutation of every row: still the whole table, as before.
		self.assertIsNone(rows_to_export(4, [], [3, 2, 1, 0]))

	def test_filtered_table_exports_the_rows_shown_in_their_order(self):
		self.assertEqual(rows_to_export(5, [], [4, 1]), [4, 1])

	def test_a_filter_that_leaves_nothing_exports_nothing(self):
		self.assertEqual(rows_to_export(5, [], []), [])

	def test_ticked_rows_win_and_follow_the_view_order(self):
		self.assertEqual(rows_to_export(5, [0, 3], [3, 2, 1, 0, 4]), [3, 0])

	def test_a_ticked_row_hidden_by_the_filter_is_still_exported(self):
		self.assertEqual(rows_to_export(5, [0, 3], [3]), [3, 0])

	def test_out_of_range_or_bogus_positions_are_ignored(self):
		self.assertEqual(rows_to_export(3, [7, -1, True, "1", 2], None), [2])
		self.assertIsNone(rows_to_export("x", [0], None))

	def test_a_state_from_a_table_of_another_size_is_not_trusted(self):
		state = table_view_state([{}] * 4, [0], [0, 1])
		self.assertEqual(export_positions(state, 4), [0])
		self.assertIsNone(export_positions(state, 5))
		self.assertIsNone(export_positions(None, 4))


class SourceSummaryTests(unittest.TestCase):
	def _frame(self):
		return pd.DataFrame(HITS, columns=HITS_COLUMNS)

	def test_counts_hits_and_features_named_below_the_domain_per_database(self):
		summary = source_assignment_summary(self._frame())
		rows = {row["Database"]: row for row in summary.to_dict("records")}
		# SILVA hit OTU1 and OTU2 but named only OTU1: OTU2 stayed at the domain.
		self.assertEqual((rows["SILVA"]["Features with a hit"], rows["SILVA"]["Named below the domain"]), (2, 1))
		self.assertEqual(rows["SILVA"]["Reads in named features"], 900)
		self.assertEqual((rows["NCBI"]["Features with a hit"], rows["NCBI"]["Named below the domain"]), (2, 2))
		self.assertEqual(rows["NCBI"]["Reads in named features"], 980)
		self.assertEqual(rows["PR2"]["Reads in named features"], 920)
		self.assertEqual(rows["PR2"]["% of reads"], 92.0)
		# Sorted by named features, then their reads.
		self.assertEqual(list(summary["Database"]), ["NCBI", "PR2", "SILVA"])

	def test_a_feature_with_many_hits_from_one_source_counts_once(self):
		frame = pd.DataFrame(
			[("OTU1", 10, 1, "NCBI", "species"), ("OTU1", 10, 2, "NCBI", "genus"), ("OTU1", 10, 3, "NCBI", "kingdom")],
			columns=["feature_id", "reads", "hit_rank", "source", "rank"],
		)
		summary = source_assignment_summary(frame)
		self.assertEqual(summary.iloc[0]["Features with a hit"], 1)
		self.assertEqual(summary.iloc[0]["Named below the domain"], 1)
		self.assertEqual(summary.iloc[0]["Reads in named features"], 10)

	def test_domain_only_hits_are_not_named(self):
		frame = pd.DataFrame(
			[("OTU1", 10, 1, "SILVA", "kingdom"), ("OTU2", 5, 1, "SILVA", "unassigned")],
			columns=["feature_id", "reads", "hit_rank", "source", "rank"],
		)
		row = source_assignment_summary(frame).iloc[0]
		self.assertEqual((row["Features with a hit"], row["Named below the domain"]), (2, 0))
		self.assertEqual(row["Reads in named features"], "")

	def test_without_a_rank_column_every_hit_counts_as_named(self):
		frame = pd.DataFrame([("OTU1", 7, "NCBI")], columns=["feature_id", "reads", "source"])
		row = source_assignment_summary(frame).iloc[0]
		self.assertEqual((row["Features with a hit"], row["Named below the domain"], row["Reads in named features"]), (1, 1, 7))

	def test_without_reads_the_counts_still_build(self):
		frame = pd.DataFrame([("OTU1", "NCBI", "species")], columns=["feature_id", "source", "rank"])
		summary = source_assignment_summary(frame)
		self.assertEqual(summary.iloc[0]["Named below the domain"], 1)
		self.assertEqual(summary.iloc[0]["Reads in named features"], "")

	def test_empty_or_shapeless_input_gives_an_empty_table(self):
		self.assertTrue(source_assignment_summary(None).empty)
		self.assertTrue(source_assignment_summary(pd.DataFrame({"x": [1]})).empty)


class HitsOrderTests(unittest.TestCase):
	"""The hits table opens with the most-found features first, each feature's hits together."""

	def test_most_found_features_first_and_hits_kept_together(self):
		from ui.callbacks_databases import _load_hits_frame

		run_dir = Path(tempfile.mkdtemp())
		self.addCleanup(shutil.rmtree, run_dir, True)
		(run_dir / "databases" / "otu").mkdir(parents=True)
		rows = [
			("OTU2", 80, 2, "NCBI", 98.5, 100.0, "Anisakis simplex", "species"),
			("OTU10", 5000, 1, "SILVA", "", "", "Gadus", "genus"),
			("OTU2", 80, 1, "SILVA", "", "", "Eukaryota", "kingdom"),
			("OTU100", 3, 1, "PR2", "", "", "Kudoa", "genus"),
			("OTU10", 5000, 2, "PR2", "", "", "Gadus", "genus"),
		]
		pd.DataFrame(rows, columns=HITS_COLUMNS).to_csv(run_dir / "databases" / "otu" / "taxonomy_hits.csv", index=False)
		frame, error = _load_hits_frame({}, run_dir, "otu", 0, 0, 0)
		self.assertIsNone(error)
		self.assertEqual(
			list(zip(frame["feature_id"], frame["hit_rank"])),
			[("OTU10", 1), ("OTU10", 2), ("OTU2", 1), ("OTU2", 2), ("OTU100", 1)],
		)


class DownloadCallbackTests(unittest.TestCase):
	"""The real callbacks, through Dash's own HTTP protocol, on a run directory on disk."""

	@classmethod
	def setUpClass(cls):
		from app import create_app

		cls.app = create_app()
		cls.client = cls.app.server.test_client()

	def setUp(self):
		self.run_dir = Path(tempfile.mkdtemp())
		self.addCleanup(shutil.rmtree, self.run_dir, True)
		db_dir = self.run_dir / "databases" / "otu"
		db_dir.mkdir(parents=True)
		with (db_dir / "taxonomy_hits.csv").open("w", newline="", encoding="utf-8") as handle:
			writer = csv.writer(handle)
			writer.writerow(HITS_COLUMNS)
			writer.writerows(HITS)
		self.summary_path = db_dir / "taxonomy_summary.csv"
		self.summary_path.write_text(
			"feature_id,feature_type,reads,sample_pct,consensus_taxon\n"
			"OTU1,otu,900,90.00,Gadus morhua\n"
			"OTU2,otu,80,8.00,Anisakis simplex\n"
			"OTU3,otu,20,2.00,Kudoa\n",
			encoding="utf-8",
		)

	def _post(self, output: str, n_clicks_id: str, state_values: dict) -> dict:
		deps = self.client.get("/_dash-dependencies").get_json()
		dep = next(d for d in deps if d["output"] == output)
		state = [{"id": s["id"], "property": s["property"], "value": state_values.get(s["id"])} for s in dep["state"]]
		body = {
			"output": output,
			"outputs": {"id": output.split(".")[0], "property": output.split(".")[1]},
			"inputs": [{"id": n_clicks_id, "property": "n_clicks", "value": 1}],
			"state": state,
			"changedPropIds": [f"{n_clicks_id}.n_clicks"],
		}
		response = self.client.post("/_dash-update-component", json=body)
		if response.status_code == 204:
			return {}
		self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[-500:])
		return response.get_json()["response"][output.split(".")[0]]["data"]

	@staticmethod
	def _text(payload: dict) -> str:
		if payload.get("base64"):
			return base64.b64decode(payload["content"]).decode("utf-8")
		return payload["content"]

	def _hits_states(self, table_state):
		return {
			"databases-feature-set": "otu",
			"hits-min-identity": 0,
			"hits-min-coverage": 0,
			"hits-top-n": 0,
			"store-run-config": {"output_dir": str(self.run_dir)},
			"hits-table-state": table_state,
		}

	def test_hits_download_without_selection_is_the_whole_filtered_table(self):
		payload = self._post("download-hits-csv.data", "btn-download-hits", self._hits_states(None))
		self.assertEqual(payload["filename"], "otu_best_hits.csv")
		rows = list(csv.DictReader(io.StringIO(self._text(payload))))
		self.assertEqual(len(rows), len(HITS))

	def test_hits_download_with_ticked_rows_holds_only_those(self):
		# Rows 0 and 3 of the table: OTU1's NCBI hit and OTU2's SILVA hit.
		state = table_view_state([{}] * len(HITS), [0, 3], None)
		payload = self._post("download-hits-csv.data", "btn-download-hits", self._hits_states(state))
		self.assertEqual(payload["filename"], "otu_best_hits_selected.csv")
		rows = list(csv.DictReader(io.StringIO(self._text(payload))))
		self.assertEqual([(r["feature_id"], r["source"]) for r in rows], [("OTU1", "NCBI"), ("OTU2", "SILVA")])

	def test_summary_download_without_selection_is_the_pipeline_file_byte_for_byte(self):
		states = {"databases-feature-set": "otu", "store-run-config": {"output_dir": str(self.run_dir)}, "feature-summary-table-state": None}
		payload = self._post("download-summary-csv.data", "btn-download-summary", states)
		self.assertEqual(payload["filename"], "otu_taxonomy_summary.csv")
		self.assertEqual(self._text(payload), self.summary_path.read_text(encoding="utf-8"))
		# A sorted but unfiltered, unticked table is still the file itself.
		states["feature-summary-table-state"] = table_view_state([{}] * 3, [], [2, 1, 0])
		payload = self._post("download-summary-csv.data", "btn-download-summary", states)
		self.assertEqual(self._text(payload), self.summary_path.read_text(encoding="utf-8"))

	def test_summary_download_with_ticked_rows_holds_only_those_with_the_same_columns(self):
		states = {
			"databases-feature-set": "otu",
			"store-run-config": {"output_dir": str(self.run_dir)},
			"feature-summary-table-state": table_view_state([{}] * 3, [1], None),
		}
		payload = self._post("download-summary-csv.data", "btn-download-summary", states)
		self.assertEqual(payload["filename"], "otu_taxonomy_summary_selected.csv")
		text = self._text(payload)
		self.assertEqual(text.splitlines()[0], "feature_id,feature_type,reads,sample_pct,consensus_taxon")
		self.assertEqual(text.splitlines()[1:], ["OTU2,otu,80,8.00,Anisakis simplex"])


if __name__ == "__main__":
	unittest.main()
