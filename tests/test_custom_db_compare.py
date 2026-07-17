"""Tests for the custom-database BLAST comparison feature."""

from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

from analysis.custom_db_compare import NO_MATCH_LABEL, build_custom_comparison, write_detection_table
from databases.custom_db import CustomDBConfig, blast_tools_available, detect_fasta_kind, reference_label, run_custom_blast


class DetectKindTests(unittest.TestCase):
    def _write(self, path: Path, text: str) -> Path:
        path.write_text(text, encoding="utf-8")
        return path

    def test_plain_fasta_is_blast_kind(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.fasta", ">seq1 Genus species\nACGT\n>seq2 Other taxon\nTTTT\n")
            kind, error = detect_fasta_kind(p)
            self.assertEqual(kind, "blast")
            self.assertEqual(error, "")

    def test_sintax_fasta_is_sintax_kind(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.fasta", ">seq1;tax=d:Eukaryota,g:Diplostomum,s:Diplostomum spathaceum;\nACGT\n")
            kind, error = detect_fasta_kind(p)
            self.assertEqual(kind, "sintax")

    def test_non_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.txt", "not a fasta\n")
            kind, error = detect_fasta_kind(p)
            self.assertEqual(kind, "")
            self.assertTrue(error)

    def test_config_validate_accepts_plain_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._write(Path(d) / "ref.fasta", ">seq1 Genus species\nACGT\n")
            cfg = CustomDBConfig(p)
            valid, error = cfg.validate()
            self.assertTrue(valid, error)
            self.assertEqual(cfg.kind, "blast")


class ReferenceLabelTests(unittest.TestCase):
    def test_sintax_species_field(self):
        self.assertEqual(reference_label("ID1;tax=d:Eukaryota,g:Diplostomum,s:Diplostomum_spathaceum;"), "Diplostomum spathaceum")

    def test_ncbi_style_genus_species(self):
        self.assertEqual(reference_label("ACC123 Gyrodactylus salaris voucher X"), "Gyrodactylus salaris")

    def test_bare_id(self):
        self.assertEqual(reference_label("ACC123"), "ACC123")


class BuildComparisonTests(unittest.TestCase):
    def _setup(self, d: Path):
        besthit = d / "besthit.csv"
        besthit.write_text(
            "OTU_ID,subject_id,species,pident,length,evalue,bitscore,status\n"
            "OTU1,REF1,Diplostomum spathaceum,99.5,300,0.0,550,Match\n"
            "OTU2,REF2,Gyrodactylus salaris,88.0,290,1e-40,300,Match\n"
            "OTU3,,None,,,,,None\n",
            encoding="utf-8",
        )
        matrix = d / "otu_matrix.csv"
        matrix.write_text(
            "#OTU ID,S1,S2\n"
            "OTU1,100,0\n"
            "OTU2,10,40\n"
            "OTU3,0,25\n",
            encoding="utf-8",
        )
        return besthit, matrix

    def test_detection_matrix_and_summary(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix)
            self.assertNotIn("error", comp)
            self.assertEqual(comp["summary"]["n_otus"], 3)
            self.assertEqual(comp["summary"]["n_matched"], 2)
            self.assertEqual(comp["summary"]["n_unmatched"], 1)
            self.assertEqual(comp["samples"], ["S1", "S2"])
            self.assertIn(NO_MATCH_LABEL, comp["taxa"])
            # No match row holds OTU3's reads (25 in S2).
            det = {name: row for name, row in zip(comp["detection"]["taxa"], comp["detection"]["values"])}
            self.assertEqual(det[NO_MATCH_LABEL], [0, 25])
            self.assertEqual(det["Diplostomum spathaceum"], [100, 0])

    def test_min_identity_demotes_weak_hits(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix, min_identity=95.0)
            # OTU2 at 88% now counts as unmatched.
            self.assertEqual(comp["summary"]["n_matched"], 1)
            otu2 = next(r for r in comp["per_otu"] if r["OTU_ID"] == "OTU2")
            self.assertEqual(otu2["status"], "None")

    def test_species_filter(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix, species_filter="diplostomum")
            self.assertEqual([r["OTU_ID"] for r in comp["per_otu"]], ["OTU1"])
            self.assertEqual(comp["taxa"], ["Diplostomum spathaceum"])

    def test_public_species_from_consensus(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            consensus = d / "consensus.csv"
            consensus.write_text(
                "ASV_ID,Consensus_Species\nOTU1,Diplostomum spathaceum\nOTU2,\n",
                encoding="utf-8",
            )
            comp = build_custom_comparison(besthit, matrix, consensus_csv=consensus)
            by_id = {r["OTU_ID"]: r for r in comp["per_otu"]}
            self.assertEqual(by_id["OTU1"]["public_species"], "Diplostomum spathaceum")
            self.assertEqual(by_id["OTU2"]["public_species"], "None")

    def test_write_detection_table(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            besthit, matrix = self._setup(d)
            comp = build_custom_comparison(besthit, matrix)
            out = write_detection_table(comp["detection"], d / "detection.csv")
            self.assertTrue(out.exists())
            text = out.read_text(encoding="utf-8")
            self.assertIn("taxon,S1,S2", text)
            self.assertIn("No match", text)

    def test_missing_besthit_returns_error(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            _, matrix = self._setup(d)
            comp = build_custom_comparison(d / "nope.csv", matrix)
            self.assertIn("error", comp)


@unittest.skipUnless(blast_tools_available(), "blastn/makeblastdb not available")
class RealBlastTests(unittest.TestCase):
    def test_blast_reports_none_for_unmatched_query(self):
        with tempfile.TemporaryDirectory() as dd:
            d = Path(dd)
            ref_seq = (
                "TCGACCTGGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTATAAGCAA"
                "TTATACCGCGAATGGCTCATTATAACAGTCGCACGTACCTTGACTATCTGGATAACCGTG"
                "GTAATTCTAGAGCTAATACATGCCGACGAGCGCTGACCCCCGGGGATGCGTGCATTTATC"
            )
            other_seq = (
                "GAACCTTACCAGGTCCAGACATAGTAAGGATTGACAGATTGAGAGCTCTTTCTTGATTTT"
                "GTGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTAATTCCGAT"
                "AACGAACGAGACTCTGGCATGCTAACTAGTTACGCGACCCCCATTCGCTAGTTGGCATCG"
            )
            reference = d / "ref.fasta"
            reference.write_text(f">REF1 Testus specius\n{ref_seq}\n", encoding="utf-8")
            query = d / "query.fasta"
            # OTU1 == reference (matches), OTU2 is unrelated (no hit).
            query.write_text(
                f">OTU1;size=5\n{ref_seq}\n>OTU2;size=3\n{other_seq}\n",
                encoding="utf-8",
            )
            result = run_custom_blast(query, reference, d / "out", perc_identity=90.0, threads=1)
            self.assertTrue(result.success, result.error_message)
            self.assertEqual(result.metrics["matched"], 1.0)
            self.assertEqual(result.metrics["total"], 2.0)
            besthit = Path(result.files["custom_blast_besthit_csv"])
            text = besthit.read_text(encoding="utf-8")
            self.assertIn("OTU1", text)
            self.assertIn("Testus specius", text)
            # The unmatched OTU2 row is present with a None status.
            import csv as _csv

            rows = {r["OTU_ID"]: r for r in _csv.DictReader(text.splitlines())}
            self.assertEqual(rows["OTU2"]["status"], "None")
            self.assertEqual(rows["OTU2"]["species"], "None")


if __name__ == "__main__":
    unittest.main()
