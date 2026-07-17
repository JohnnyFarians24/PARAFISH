"""Regression tests for the screen-recording bug fixes:

* sample-tag loss (VSEARCH header truncation) → correct per-sample columns;
* tree ↔ abundance-table ID mismatch (Faith's PD ASV mapping);
* honest per-source status classification (OOM/failed vs no-hits);
* detection heatmap degenerate (all "No match") case.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


class SampleHeaderNormalizationTests(unittest.TestCase):
    def test_sample_tag_survives_whitespace_in_header(self):
        from analysis.vsearch_runner import _normalize_sample_headers

        with tempfile.TemporaryDirectory() as d:
            fasta = Path(d) / "filtered.fasta"
            # Illumina header: description after a space, sample tag at the very end.
            fasta.write_text(
                ">A01426:1215:HW7YLDRX5:2:2101:16794:1423 1:N:0:ACGT+CTCT;sample=aptl64b_02\nACGT\n",
                encoding="utf-8",
            )
            _normalize_sample_headers(fasta, "aptl64b_02")
            header = next(l.strip() for l in fasta.open() if l.startswith(">"))
            # No whitespace before the sample tag → VSEARCH keeps it.
            self.assertEqual(header, ">A01426:1215:HW7YLDRX5:2:2101:16794:1423;sample=aptl64b_02")
            self.assertNotIn(" ", header)


class FaithPDMappingTests(unittest.TestCase):
    def _tree(self, d: Path) -> Path:
        # Two tips named ASV1/ASV2 (as the phylo stage renames them).
        tree = d / "tree.nwk"
        tree.write_text("(ASV1:0.3,ASV2:0.4);\n", encoding="utf-8")
        return tree

    def _matrix(self, d: Path) -> Path:
        # Abundance table keyed by the ORIGINAL headers (with ;size=), not ASV ids.
        matrix = d / "otu_matrix.csv"
        matrix.write_text(
            "#OTU ID,S1,S2\n"
            "A01426:1:1;size=10,5,0\n"
            "A01426:2:2;size=3,0,7\n",
            encoding="utf-8",
        )
        return matrix

    def _mapping(self, d: Path) -> Path:
        m = d / "reversed_asvs.csv"
        m.write_text(
            "ASV_ID,Original_Header\n"
            "ASV1,A01426:1:1;size=10\n"
            "ASV2,A01426:2:2;size=3\n",
            encoding="utf-8",
        )
        return m

    def test_without_mapping_zero_match(self):
        from phylo.phylo_diversity import faith_pd

        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            res = faith_pd(self._tree(d), self._matrix(d))
            self.assertEqual(res["matched_features"], 0)
            self.assertEqual(res["pd"], [0, 0])

    def test_with_mapping_matches_and_computes_pd(self):
        from phylo.phylo_diversity import faith_pd

        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            res = faith_pd(self._tree(d), self._matrix(d), id_mapping_csv=self._mapping(d))
            self.assertEqual(res["matched_features"], 1)  # each sample has one present tip
            self.assertGreater(res["pd"][0], 0)
            self.assertGreater(res["pd"][1], 0)


class SourceStatusClassificationTests(unittest.TestCase):
    def test_oom_kill_classified_as_failed(self):
        from databases.taxonomy_runner import _classify_sintax_source

        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            result = d / "silva_sintax_result.csv"
            result.write_text("", encoding="utf-8")  # empty (killed mid-run)
            log = d / "silva_sintax.log"
            log.write_text("VSEARCH SINTAX foi terminado pelo sistema (sinal 9), OOM.", encoding="utf-8")
            status, detail = _classify_sintax_source(result, log)
            self.assertEqual(status, "failed")
            self.assertIn("OOM", detail)

    def test_empty_without_oom_is_no_hits(self):
        from databases.taxonomy_runner import _classify_sintax_source

        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            result = d / "r.csv"
            result.write_text("", encoding="utf-8")
            log = d / "l.log"
            log.write_text("ok, ran fine, matched 0", encoding="utf-8")
            status, _ = _classify_sintax_source(result, log)
            self.assertEqual(status, "no_hits")

    def test_assignments_classified_ok(self):
        from databases.taxonomy_runner import _classify_sintax_source

        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            result = d / "r.csv"
            result.write_text("ASV1\td:Eukaryota\t+\td:Eukaryota\nASV2\td:Eukaryota\t+\td:Eukaryota\n", encoding="utf-8")
            status, detail = _classify_sintax_source(result, None)
            self.assertEqual(status, "ok")


class DetectionHeatmapDegenerateTests(unittest.TestCase):
    def test_all_no_match_returns_none(self):
        from ui.plots import detection_heatmap_figure

        detection = {"taxa": ["No match"], "samples": ["S1", "S2"], "values": [[10, 20]]}
        self.assertIsNone(detection_heatmap_figure(detection))

    def test_with_matched_taxon_builds_figure(self):
        import plotly.graph_objects as go

        from ui.plots import detection_heatmap_figure

        detection = {"taxa": ["Diplostomum spathaceum", "No match"], "samples": ["S1", "S2"], "values": [[8, 0], [2, 20]]}
        fig = detection_heatmap_figure(detection)
        self.assertIsInstance(fig, go.Figure)


if __name__ == "__main__":
    unittest.main()


class SintaxParsingTests(unittest.TestCase):
    """VSEARCH --tabbedout is a headerless TSV; csv.DictReader ate the first hit
    as a header and no SINTAX source ever reached the taxonomy consensus."""

    ROW = (
        "A01426:1:HW7:2:2101:6289:3458;sample=s1;size=12618\t"
        "d:Eukaryota(1.00),p:Obazoa(0.99),c:Opisthokonta(0.98),o:Metazoa(0.95),"
        "f:Craniata(0.91),g:Craniata_X(0.88),s:Teleostei(0.64)\t"
        "+\t"
        "d:Eukaryota,p:Obazoa,c:Opisthokonta,o:Metazoa,f:Craniata,g:Craniata_X\n"
    )

    def _write(self, directory: Path) -> Path:
        path = directory / "pr2_sintax_result.csv"
        path.write_text(self.ROW, encoding="utf-8")
        return path

    def test_headerless_sintax_rows_are_parsed(self):
        from databases.taxonomy_runner import _collect_sintax_records

        with tempfile.TemporaryDirectory() as d:
            rows = _collect_sintax_records("PR2", self._write(Path(d)))

        self.assertEqual(list(rows), ["A01426:1:HW7:2:2101:6289:3458"])
        record = rows["A01426:1:HW7:2:2101:6289:3458"][0]
        self.assertEqual(record["kingdom"], "Eukaryota")
        self.assertEqual(record["genus"], "Craniata_X")
        # "s:Teleostei" fell below the cutoff, so it is absent from column 4.
        self.assertEqual(record["species"], "")
        # The bootstrap for the deepest surviving rank comes from the raw column.
        self.assertEqual(record["confidence"], "0.88")

    def test_sintax_records_reach_the_consensus(self):
        from databases.taxonomy_runner import _build_consensus_row, _collect_sintax_records

        with tempfile.TemporaryDirectory() as d:
            rows = _collect_sintax_records("PR2", self._write(Path(d)))
        asv_id, records = next(iter(rows.items()))
        consensus = _build_consensus_row(asv_id, records, "lca")
        self.assertEqual(consensus["Consensus_Genus"], "Craniata_X")
        self.assertEqual(consensus["SINTAX_bootstrap"], "0.88")


class BlastRecordTests(unittest.TestCase):
    def test_headerless_blast_table_keeps_only_the_best_hit(self):
        from databases.taxonomy_runner import _collect_blast_records

        outfmt = "6 qseqid sseqid pident length evalue bitscore sscinames"
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "blast.tsv"
            path.write_text(
                "ASV1\tX1\t99.0\t300\t1e-50\t120.0\tGadus morhua\n"
                "ASV1\tX2\t95.0\t300\t1e-40\t300.5\tSalmo salar\n",
                encoding="utf-8",
            )
            rows = _collect_blast_records("NCBI", path, outfmt)

        self.assertEqual(len(rows["ASV1"]), 1)
        self.assertEqual(rows["ASV1"][0]["species"], "Salmo salar")
        self.assertEqual(rows["ASV1"][0]["similarity"], "95.0")


class BoldColumnCaseTests(unittest.TestCase):
    def test_capitalised_bold_columns_are_read(self):
        """BOLDigger3 emits Species/Genus/pct_identity; the old exact/lower/upper
        probe missed them and silently dropped every BOLD assignment."""
        from databases.taxonomy_runner import _normalize_source_row

        row = {"id": "ASV1", "Phylum": "Chordata", "Genus": "Salmo", "Species": "Salmo salar", "pct_identity": "98.7", "flags": ""}
        record = _normalize_source_row("BOLD", row)
        self.assertEqual(record["asv_id"], "ASV1")
        self.assertEqual(record["genus"], "Salmo")
        self.assertEqual(record["species"], "Salmo salar")
        self.assertEqual(record["similarity"], "98.7")

    def test_bold_no_match_placeholders_are_dropped(self):
        from databases.taxonomy_runner import _normalize_source_row

        record = _normalize_source_row("BOLD", {"id": "ASV1", "Species": "no-match", "Genus": "IncompleteTaxonomy"})
        self.assertEqual(record["species"], "")
        self.assertEqual(record["genus"], "")


class AnnotatedFastaTests(unittest.TestCase):
    def test_annotation_matches_ids_carrying_vsearch_suffixes(self):
        from databases.taxonomy_runner import _annotate_fasta

        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "otus.fasta"
            src.write_text(">ASV1;sample=s1;size=12\nACGT\n", encoding="utf-8")
            out = Path(d) / "annotated.fasta"
            _annotate_fasta(src, {"ASV1": {"Consensus_Genus": "Salmo"}}, out)
            header = out.read_text(encoding="utf-8").splitlines()[0]

        self.assertEqual(header, ">ASV1;size=12;tax=;;;;;Salmo;")


class ResultsSampleFilterTests(unittest.TestCase):
    def test_only_real_sample_keys_are_offered(self):
        """The analysis record mixes per-sample entries with run-wide artifacts."""
        from ui.callbacks_results import _sample_ids

        outputs = {
            "s1": {"vsearch_filter": {"success": True}},
            "s2": {"vsearch_filter": {"success": True}},
            "vsearch_pool": {"success": True},
            "otu_matrix": "/runs/x/otu_matrix.csv",
            "combined_fasta": "/runs/x/otus.fasta",
            "dada2": {"files": {}},
        }
        self.assertEqual(sorted(_sample_ids("analysis", outputs)), ["s1", "s2"])

    def test_stages_without_per_sample_outputs_offer_none(self):
        from ui.callbacks_results import _sample_ids

        self.assertEqual(_sample_ids("databases", {"taxonomy": {"files": {}}}), [])
