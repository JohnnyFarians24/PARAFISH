"""Stage 1 (QC): FastQC module parsing, persisted metrics and control metadata."""

from __future__ import annotations

import csv
import tempfile
import unittest
import zipfile
from pathlib import Path


FASTQC_DATA = """##FastQC\t0.12.1
>>Basic Statistics\tpass
#Measure\tValue
Filename\tS_R1.fastq.gz
Total Sequences\t1000
Sequence length\t35-256
%GC\t54
>>END_MODULE
>>Per base sequence quality\tpass
#Base\tMean\tMedian
1\t35.5\t36
2-3\t34.0\t34
>>END_MODULE
>>Per sequence quality scores\tpass
#Quality\tCount
20\t250
35\t750
>>END_MODULE
>>Sequence Length Distribution\twarn
#Length\tCount
150\t400
200-201\t600
>>END_MODULE
>>Sequence Duplication Levels\tfail
#Total Deduplicated Percentage\t40.0
#Duplication Level\tPercentage of deduplicated
1\t80.0
>>END_MODULE
>>Adapter Content\tpass
#Position\tIllumina Universal Adapter\tNextera
1\t0.1\t0.0
2\t7.5\t0.2
>>END_MODULE
"""


def _write_fastqc_zip(path: Path, data: str = FASTQC_DATA) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{path.stem}/fastqc_data.txt", data)
    return path


class FastqcDataParserTests(unittest.TestCase):
    def _parse(self, data: str = FASTQC_DATA) -> dict:
        from qc.fastqc_parser import parse_fastqc_data

        with tempfile.TemporaryDirectory() as d:
            return parse_fastqc_data(_write_fastqc_zip(Path(d) / "S_R1_fastqc.zip", data))

    def test_module_statuses_are_captured(self):
        status = self._parse()["module_status"]
        self.assertEqual(status["Basic Statistics"], "pass")
        self.assertEqual(status["Sequence Duplication Levels"], "fail")
        self.assertEqual(status["Sequence Length Distribution"], "warn")

    def test_basic_statistics(self):
        payload = self._parse()
        self.assertEqual(payload["total_sequences"], 1000.0)
        self.assertEqual(payload["pct_gc"], 54.0)

    def test_binned_positions_use_the_midpoint(self):
        # FastQC bins positions as "2-3"; the midpoint is 2.5.
        self.assertEqual(self._parse()["per_base_quality"], [(1.0, 35.5), (2.5, 34.0)])

    def test_pct_reads_q30_is_share_of_reads_not_bases(self):
        # 750 of 1000 reads have a mean quality >= 30.
        self.assertEqual(self._parse()["pct_reads_q30"], 75.0)

    def test_mean_length_is_abundance_weighted(self):
        # (150*400 + 200.5*600) / 1000
        self.assertEqual(self._parse()["mean_length"], 180.3)

    def test_duplication_is_the_complement_of_the_deduplicated_percentage(self):
        self.assertEqual(self._parse()["pct_duplication"], 60.0)

    def test_adapter_content_takes_the_worst_position(self):
        self.assertEqual(self._parse()["pct_adapter_max"], 7.5)

    def test_unreadable_zip_yields_empty_payload(self):
        from qc.fastqc_parser import parse_fastqc_data

        with tempfile.TemporaryDirectory() as d:
            junk = Path(d) / "not_a_zip.zip"
            junk.write_text("definitely not a zip", encoding="utf-8")
            self.assertEqual(parse_fastqc_data(junk), {})
            self.assertEqual(parse_fastqc_data(Path(d) / "missing.zip"), {})


class QcMetricsPersistenceTests(unittest.TestCase):
    def _build(self, depths: dict[str, int]):
        from qc.qc_metrics import build_qc_metrics

        tmp = tempfile.mkdtemp()
        qc_dir = Path(tmp)
        for sample, depth in depths.items():
            data = FASTQC_DATA.replace("Total Sequences\t1000", f"Total Sequences\t{depth}")
            _write_fastqc_zip(qc_dir / sample / f"{sample.upper()}_R1_fastqc.zip", data)
        return qc_dir, build_qc_metrics(qc_dir)

    def test_sample_id_comes_from_the_directory_not_the_filename(self):
        qc_dir, files = self._build({"aptl64f_01": 1000})
        rows = list(csv.DictReader(open(files["qc_overview_csv"])))
        self.assertEqual(rows[0]["sample"], "aptl64f_01")
        self.assertEqual(rows[0]["read"], "R1")
        del qc_dir

    def test_overview_carries_the_enriched_metrics(self):
        _, files = self._build({"s1": 1000})
        row = list(csv.DictReader(open(files["qc_overview_csv"])))[0]
        self.assertEqual(row["total_reads"], "1000")
        self.assertEqual(row["pct_reads_q30"], "75.0")
        self.assertEqual(row["pct_duplication"], "60.0")
        self.assertEqual(row["pct_adapter_max"], "7.5")
        self.assertEqual(row["n_fail"], "1")
        self.assertEqual(row["n_warn"], "1")

    def test_low_depth_samples_are_flagged_against_the_median(self):
        # median depth is 1000; 100 is below half of it, 900 is not.
        _, files = self._build({"deep_a": 1000, "deep_b": 1200, "shallow": 100})
        flags = {row["sample"]: row["is_depth_outlier"] for row in csv.DictReader(open(files["qc_overview_csv"]))}
        self.assertEqual(flags["shallow"], "1")
        self.assertEqual(flags["deep_a"], "0")
        self.assertEqual(flags["deep_b"], "0")

    def test_single_sample_run_flags_nothing(self):
        _, files = self._build({"only": 10})
        rows = list(csv.DictReader(open(files["qc_overview_csv"])))
        self.assertEqual(rows[0]["is_depth_outlier"], "0")

    def test_module_status_table_is_long_form(self):
        _, files = self._build({"s1": 1000})
        rows = list(csv.DictReader(open(files["qc_module_status_csv"])))
        statuses = {row["module"]: row["status"] for row in rows}
        self.assertEqual(statuses["Sequence Duplication Levels"], "fail")

    def test_missing_qc_dir_returns_nothing(self):
        from qc.qc_metrics import build_qc_metrics

        self.assertEqual(build_qc_metrics(Path("/does/not/exist")), {})


class SampleControlMetadataTests(unittest.TestCase):
    def test_control_and_group_round_trip_through_the_store(self):
        from ui.runtime_state import build_sample_store_entry, sample_record_to_sample

        with tempfile.TemporaryDirectory() as d:
            r1 = Path(d) / "r1.fastq"
            r1.write_text("", encoding="utf-8")
            sample = sample_record_to_sample(
                {"internal_id": "blank_01", "display_label": "Blank", "r1_path": str(r1), "mode": "SE", "is_control": True, "group": " site A "}
            )
            self.assertTrue(sample.is_control)
            self.assertEqual(sample.group, "site A")

            entry = build_sample_store_entry(sample, r1)
            self.assertTrue(entry["is_control"])
            self.assertEqual(entry["group"], "site A")

    def test_defaults_keep_older_sample_records_working(self):
        from ui.runtime_state import sample_record_to_sample

        with tempfile.TemporaryDirectory() as d:
            r1 = Path(d) / "r1.fastq"
            r1.write_text("", encoding="utf-8")
            sample = sample_record_to_sample({"internal_id": "s1", "r1_path": str(r1), "mode": "SE"})
            self.assertFalse(sample.is_control)
            self.assertIsNone(sample.group)

    def test_blank_group_string_is_normalised_to_none(self):
        from ui.runtime_state import sample_record_to_sample

        with tempfile.TemporaryDirectory() as d:
            r1 = Path(d) / "r1.fastq"
            r1.write_text("", encoding="utf-8")
            sample = sample_record_to_sample({"internal_id": "s1", "r1_path": str(r1), "mode": "SE", "group": "   "})
            self.assertIsNone(sample.group)


if __name__ == "__main__":
    unittest.main()
