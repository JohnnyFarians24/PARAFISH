"""Exercise the BOLD positive-result parsing that no real run ever reached — and that,
tested against the REAL output format, turns out to be broken.

BOLD has produced zero attributions in this project — every run hit COI-vs-18S marker
mismatch, so ``_read_id_engine_hits`` (reads BOLDigger3's raw hit store) and
``select_bold_top_hits`` (threshold -> rank collapse) never ran, in production OR in a
test: the existing BOLD tests feed already-parsed dicts and stub files, skipping both.

The synthetic DuckDB below uses the REAL BOLDigger3 3.x schema — verified against the
``id_engine_results`` table of the actual ``golden_baseline`` run: an ``id``, a
``pct_identity``, and the Linnaean ranks in **lowercase** (``phylum``..``species``), plus
``bin_uri`` / ``process_id``.

Two findings:
  * ``_read_id_engine_hits`` reads the store correctly — TESTED, passes.
  * ``select_bold_top_hits`` looks up ranks by ``_BOLD_RANK_ORDER`` = ``["Phylum".."Species"]``
    (CAPITALISED), which are absent from the lowercase DuckDB, so it blanks every rank even
    for a 99.5% hit. A real positive BOLD result would arrive with identity but no taxonomy.
    Marked ``expectedFailure`` — a known bug, NOT fixed here (reported for a supervisor
    decision on whether "BOLD functional" can be promised). When the case mismatch is fixed,
    this flips to an unexpected success and forces the marker off.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

try:
    import duckdb
    import pandas as pd

    _DEPS = True
except Exception:  # pragma: no cover
    _DEPS = False

# Real BOLDigger3 3.x id_engine_results columns are lowercase (see golden_baseline duckdb).
_REAL_HITS = [
    {"id": "ASV1", "pct_identity": 99.5, "phylum": "Platyhelminthes", "class": "Trematoda",
     "order": "Diplostomida", "family": "Diplostomidae", "genus": "Diplostomum",
     "species": "Diplostomum pseudospathaceum", "bin_uri": "BOLD:AAA1234", "process_id": "PROC1"},
    {"id": "ASV1", "pct_identity": 99.1, "phylum": "Platyhelminthes", "class": "Trematoda",
     "order": "Diplostomida", "family": "Diplostomidae", "genus": "Diplostomum",
     "species": "Diplostomum pseudospathaceum", "bin_uri": "BOLD:AAA1234", "process_id": "PROC2"},
]


@unittest.skipUnless(_DEPS, "duckdb/pandas required")
class BoldPositiveParsingTests(unittest.TestCase):
    def _write_duckdb(self, tmp: Path) -> None:
        data_dir = tmp / "boldigger3_data"
        data_dir.mkdir(parents=True, exist_ok=True)
        rows = pd.DataFrame(_REAL_HITS)
        with duckdb.connect(str(data_dir / "query.duckdb")) as con:
            con.register("rows", rows)
            con.execute("CREATE TABLE id_engine_results AS SELECT * FROM rows")

    def test_read_id_engine_hits_reads_the_store(self):
        # The raw-store reader works — this half of the positive path is sound.
        from databases.bold_db import _read_id_engine_hits

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            self._write_duckdb(tmp)
            hits = _read_id_engine_hits(tmp, "query")
            self.assertIsNotNone(hits, "the raw DuckDB hit store was not read")
            self.assertEqual(set(hits["id"]), {"ASV1"})
            self.assertEqual(float(hits["pct_identity"].max()), 99.5)

    @unittest.expectedFailure  # KNOWN BUG: rank lookup is capitalised, real columns are lowercase.
    def test_select_bold_top_hits_fills_species_for_a_99pct_hit(self):
        from databases.bold_db import _read_id_engine_hits, select_bold_top_hits

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            self._write_duckdb(tmp)
            hits = _read_id_engine_hits(tmp, "query")
            top = select_bold_top_hits(hits, [97, 95, 90, 85]).to_dict("records")[0]
            # A 99.5% hit must earn a species — today it comes back blank.
            self.assertEqual(top.get("Species") or top.get("species"), "Diplostomum pseudospathaceum")


if __name__ == "__main__":
    unittest.main()
