"""Stage 2 (Preprocess): read-tracking waterfall parsed from the runner logs."""

from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path


TRIM_PE_LOG = """TrimmomaticPE: Started with arguments:
Quality encoding detected as phred33
Input Read Pairs: 101190 Both Surviving: 99618 (98.45%) Forward Only Surviving: 935 (0.92%) Reverse Only Surviving: 477 (0.47%) Dropped: 160 (0.16%)
TrimmomaticPE: Completed successfully
"""

TRIM_SE_LOG = """TrimmomaticSE: Started with arguments:
Input Reads: 1000 Surviving: 900 (90.00%) Dropped: 100 (10.00%)
TrimmomaticSE: Completed successfully
"""

PEAR_LOG = """Assemblying reads: 100%

Assembled reads ...................: 94,657 / 99,618 (95.020%)
Discarded reads ...................: 4,961 / 99,618 (4.980%)
Not assembled reads ...............: 0 / 99,618 (0.000%)
Assembled reads file...............: /tmp/x.assembled.fastq
"""


def _fastq(records: list[str]) -> str:
    return "".join(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(records))


def _make_sample(root: Path, name: str, trim_log: str, pear_log: str | None = None, assembled: list[str] | None = None):
    sample_dir = root / name
    sample_dir.mkdir(parents=True)
    (sample_dir / "trimlog.txt").write_text(trim_log, encoding="utf-8")
    if pear_log is not None:
        pear_dir = sample_dir / "pear"
        pear_dir.mkdir()
        (pear_dir / "pear.log").write_text(pear_log, encoding="utf-8")
        if assembled:
            (pear_dir / f"{name}.assembled.fastq").write_text(_fastq(assembled), encoding="utf-8")
    return sample_dir


class TrimmomaticLogTests(unittest.TestCase):
    def test_paired_end_accounting(self):
        from preprocess.preprocess_metrics import parse_trimmomatic_log

        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "trimlog.txt"
            log.write_text(TRIM_PE_LOG, encoding="utf-8")
            parsed = parse_trimmomatic_log(log)

        self.assertEqual(parsed["mode"], "PE")
        self.assertEqual(parsed["input_reads"], 101190)
        self.assertEqual(parsed["trimmed_reads"], 99618)
        self.assertEqual(parsed["forward_only"], 935)
        self.assertEqual(parsed["reverse_only"], 477)
        self.assertEqual(parsed["dropped_reads"], 160)

    def test_single_end_accounting(self):
        from preprocess.preprocess_metrics import parse_trimmomatic_log

        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "trimlog.txt"
            log.write_text(TRIM_SE_LOG, encoding="utf-8")
            parsed = parse_trimmomatic_log(log)

        self.assertEqual(parsed["mode"], "SE")
        self.assertEqual(parsed["input_reads"], 1000)
        self.assertEqual(parsed["trimmed_reads"], 900)
        self.assertEqual(parsed["dropped_reads"], 100)

    def test_missing_or_unparseable_log(self):
        from preprocess.preprocess_metrics import parse_trimmomatic_log

        with tempfile.TemporaryDirectory() as d:
            noise = Path(d) / "trimlog.txt"
            noise.write_text("Exception in thread main\n", encoding="utf-8")
            self.assertEqual(parse_trimmomatic_log(noise), {})
            self.assertEqual(parse_trimmomatic_log(Path(d) / "absent.txt"), {})


class PearLogTests(unittest.TestCase):
    def test_thousand_separators_are_stripped(self):
        from preprocess.preprocess_metrics import parse_pear_log

        with tempfile.TemporaryDirectory() as d:
            log = Path(d) / "pear.log"
            log.write_text(PEAR_LOG, encoding="utf-8")
            parsed = parse_pear_log(log)

        self.assertEqual(parsed["assembled_reads"], 94657)
        self.assertEqual(parsed["discarded_reads"], 4961)
        self.assertEqual(parsed["not_assembled_reads"], 0)

    def test_absent_log_yields_empty(self):
        from preprocess.preprocess_metrics import parse_pear_log

        self.assertEqual(parse_pear_log(Path("/does/not/exist.log")), {})


class ReadTrackingTests(unittest.TestCase):
    def _build(self, root: Path):
        from preprocess.preprocess_metrics import build_preprocess_metrics

        return build_preprocess_metrics(root)

    def test_fate_categories_partition_the_input_exactly(self):
        """The whole point of the waterfall: no read may vanish unaccounted for."""
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "s1", TRIM_PE_LOG, PEAR_LOG)
            files = self._build(root)
            row = list(csv.DictReader(open(files["read_tracking_csv"])))[0]

        accounted = (
            int(row["assembled_reads"])
            + int(row["discarded_reads"])
            + int(row["not_assembled_reads"])
            + int(row["forward_only"])
            + int(row["reverse_only"])
            + int(row["dropped_reads"])
        )
        self.assertEqual(accounted, int(row["input_reads"]))

    def test_percentages_are_relative_to_the_right_denominator(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "s1", TRIM_PE_LOG, PEAR_LOG)
            files = self._build(root)
            row = list(csv.DictReader(open(files["read_tracking_csv"])))[0]

        self.assertEqual(row["pct_trimmed"], "98.45")  # of input
        self.assertEqual(row["pct_assembled"], "95.02")  # of trimmed, not of input
        self.assertEqual(row["pct_surviving"], "93.54")  # of input

    def test_single_end_sample_survives_as_trimmed_reads(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "se1", TRIM_SE_LOG)
            files = self._build(root)
            row = list(csv.DictReader(open(files["read_tracking_csv"])))[0]

        self.assertEqual(row["mode"], "SE")
        self.assertEqual(row["surviving_reads"], "900")
        self.assertEqual(row["assembled_reads"], "")
        self.assertEqual(row["pct_assembled"], "")

    def test_failing_steps_flags_a_step_that_loses_most_reads(self):
        from preprocess.preprocess_metrics import failing_steps

        bad_pear = PEAR_LOG.replace("94,657 / 99,618 (95.020%)", "20,000 / 99,618 (20.076%)")
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "s1", TRIM_PE_LOG, bad_pear)
            files = self._build(root)
            row = list(csv.DictReader(open(files["read_tracking_csv"])))[0]

        self.assertEqual(failing_steps(row), ["merging"])

    def test_healthy_sample_flags_nothing(self):
        from preprocess.preprocess_metrics import failing_steps

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "s1", TRIM_PE_LOG, PEAR_LOG)
            files = self._build(root)
            row = list(csv.DictReader(open(files["read_tracking_csv"])))[0]

        self.assertEqual(failing_steps(row), [])

    def test_fragment_lengths_are_counted_from_the_merged_fastq(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "s1", TRIM_PE_LOG, PEAR_LOG, assembled=["A" * 300, "A" * 300, "A" * 250])
            files = self._build(root)
            rows = {int(row["length"]): int(row["count"]) for row in csv.DictReader(open(files["assembled_length_csv"]))}

        self.assertEqual(rows, {250: 1, 300: 2})

    def test_gzipped_merged_fastq_is_read(self):
        from preprocess.preprocess_metrics import assembled_length_counts

        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "x.assembled.fastq.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(_fastq(["ACGT", "ACGTACGT"]))
            self.assertEqual(dict(assembled_length_counts(path)), {4: 1, 8: 1})

    def test_metrics_dir_is_not_mistaken_for_a_sample(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            _make_sample(root, "s1", TRIM_PE_LOG, PEAR_LOG)
            self._build(root)
            # Second pass now sees the metrics/ directory it just created.
            files = self._build(root)
            samples = [row["sample"] for row in csv.DictReader(open(files["read_tracking_csv"]))]

        self.assertEqual(samples, ["s1"])

    def test_empty_dir_returns_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(self._build(Path(d)), {})
        self.assertEqual(self._build(Path("/does/not/exist")), {})


if __name__ == "__main__":
    unittest.main()


class PreprocessFigureTests(unittest.TestCase):
    """Single-end rows leave the PEAR columns empty; "" is neither NaN nor 0."""

    def _tracking(self, rows):
        import pandas as pd

        return pd.DataFrame(rows)

    def _se_row(self):
        return dict(
            sample="se1", mode="SE", input_reads=1000, trimmed_reads=900, forward_only=0, reverse_only=0,
            dropped_reads=100, assembled_reads="", discarded_reads="", not_assembled_reads="",
            surviving_reads=900, pct_trimmed=90.0, pct_assembled="", pct_surviving=90.0,
        )

    def _pe_row(self, sample="pe1", pct_assembled=95.0, assembled=930, discarded=50):
        return dict(
            sample=sample, mode="PE", input_reads=1000, trimmed_reads=980, forward_only=10, reverse_only=5,
            dropped_reads=5, assembled_reads=assembled, discarded_reads=discarded, not_assembled_reads=0,
            surviving_reads=assembled, pct_trimmed=98.0, pct_assembled=pct_assembled, pct_surviving=93.0,
        )

    def test_single_end_fate_omits_the_merge_categories(self):
        from ui.plots import read_fate_figure

        names = [trace.name for trace in read_fate_figure(self._tracking([self._se_row()])).data]
        self.assertEqual(names, ["Dropped by trimming"])

    def test_single_end_has_no_assembly_rate_chart(self):
        from ui.plots import assembly_rate_figure

        self.assertIsNone(assembly_rate_figure(self._tracking([self._se_row()])))

    def test_assembly_rate_colours_failing_samples_red(self):
        from ui.plots import assembly_rate_figure

        tracking = self._tracking([self._pe_row("ok", 95.0), self._pe_row("bad", 20.0, assembled=200, discarded=780)])
        colours = assembly_rate_figure(tracking).data[0].marker.color
        self.assertEqual(tuple(colours), ("#3fb950", "#f85149"))

    def test_mixed_single_and_paired_run_plots_only_the_paired_rates(self):
        from ui.plots import assembly_rate_figure

        tracking = self._tracking([self._se_row(), self._pe_row("pe1", 95.0)])
        figure = assembly_rate_figure(tracking)
        self.assertEqual(list(figure.data[0].x), ["pe1"])

    def test_zero_categories_are_dropped_from_the_stack(self):
        from ui.plots import read_fate_figure

        names = [trace.name for trace in read_fate_figure(self._tracking([self._pe_row(discarded=0)])).data]
        self.assertNotIn("Discarded by PEAR", names)
        self.assertNotIn("Not merged", names)
        self.assertIn("Merged (kept)", names)
