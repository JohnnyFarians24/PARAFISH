"""Regression tests for the run-review findings.

Each test pins down a bug that silently produced wrong science rather than an
error: ranks shifted onto the wrong level, BLAST killed by a zero timeout, a
reliability flag that was true by construction, and singletons clustered as if
they were real biology.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


# Real PR2 v5 lineages: nine levels, and *not* Linnean-aligned.
PR2_FISH = "Eukaryota;Obazoa;Opisthokonta;Metazoa;Craniata;Craniata_X;Teleostei;Auxis;Auxis_thazard".split(";")
PR2_YEAST = "Eukaryota;Obazoa;Opisthokonta;Fungi;Ascomycota;Saccharomycotina;Saccharomycetales;Cyberlindnera;Cyberlindnera_subsufficiens".split(";")


class Pr2RankMappingTests(unittest.TestCase):
    def _ranks(self, fields):
        from core.utils import PR2_RANK_LABELS, taxonomy_fields_to_sintax

        return dict(part.split(":", 1) for part in taxonomy_fields_to_sintax(fields, PR2_RANK_LABELS).split(","))

    def test_the_real_species_reaches_the_species_rank(self):
        """Positional mapping dropped Auxis_thazard entirely and called a class the species."""
        ranks = self._ranks(PR2_FISH)
        self.assertEqual(ranks["s"], "Auxis_thazard")
        self.assertEqual(ranks["g"], "Auxis")

    def test_a_class_is_no_longer_labelled_a_species(self):
        ranks = self._ranks(PR2_FISH)
        self.assertNotEqual(ranks["s"], "Teleostei")

    def test_an_order_is_no_longer_labelled_a_species(self):
        ranks = self._ranks(PR2_YEAST)
        self.assertEqual(ranks["s"], "Cyberlindnera_subsufficiens")
        self.assertNotEqual(ranks["s"], "Saccharomycetales")

    def test_supergroup_does_not_overwrite_the_domain(self):
        """SINTAX treats d: and k: as one rank, so folding Supergroup onto k lost the domain."""
        ranks = self._ranks(PR2_FISH)
        self.assertEqual(ranks["d"], "Eukaryota")
        self.assertNotIn("k", ranks)

    def test_the_old_positional_mapping_was_wrong(self):
        """Guards the bug itself: the default mapping still mislabels PR2."""
        from core.utils import taxonomy_fields_to_sintax

        positional = dict(part.split(":", 1) for part in taxonomy_fields_to_sintax(PR2_FISH).split(","))
        self.assertEqual(positional["s"], "Teleostei")  # a class, on the species rank


class SilvaRankMappingTests(unittest.TestCase):
    def _ranks(self, fields):
        from core.utils import silva_rank_labels, taxonomy_fields_to_sintax

        return dict(part.split(":", 1) for part in taxonomy_fields_to_sintax(fields, silva_rank_labels(len(fields))).split(","))

    def test_the_organism_name_lands_on_species_whatever_the_depth(self):
        deep = "Eukaryota;Archaeplastida;Chloroplastida;Chlorophyta;Chlorophyceae;Sphaeropleales;Monoraphidium;Monoraphidium sp.".split(";")
        shallow = "Bacteria;Proteobacteria;Vibrio;Vibrio cholerae".split(";")

        self.assertEqual(self._ranks(deep)["s"], "Monoraphidium_sp.")
        self.assertEqual(self._ranks(shallow)["s"], "Vibrio_cholerae")

    def test_the_last_but_one_level_is_no_longer_claimed_as_the_genus(self):
        """It is the genus often enough to look right, and not often enough to be right.

        For these two lineages the level before the organism genuinely IS the genus
        (Monoraphidium, Vibrio), which is why claiming ``g`` there survived review. On
        SILVA's coarse metazoan paths the same position holds an infraclass: the
        reference run's cod came out as Genus=Teleostei. SILVA declares no ranks, so
        the position cannot tell the two cases apart and the rank is left unset —
        PR2/NCBI carry the genus instead.
        """
        deep = "Eukaryota;Archaeplastida;Chloroplastida;Chlorophyta;Chlorophyceae;Sphaeropleales;Monoraphidium;Monoraphidium sp.".split(";")
        fish = "Eukaryota;Amorphea;Obazoa;Opisthokonta;Holozoa;Teleostei;Gadus morhua".split(";")

        self.assertNotIn("g", self._ranks(deep))
        self.assertNotIn("g", self._ranks(fish))
        self.assertEqual(self._ranks(fish)["s"], "Gadus_morhua")

    def test_the_domain_stays_at_the_top(self):
        self.assertEqual(self._ranks("Bacteria;Proteobacteria;Vibrio;Vibrio cholerae".split(";"))["d"], "Bacteria")


class SintaxCacheVersionTests(unittest.TestCase):
    def test_the_cache_name_carries_the_format_version(self):
        """mtime cannot tell that the rank *mapping* changed, so the name must."""
        from core.utils import SINTAX_FORMAT_VERSION
        from databases.silva_db import SILVAConfig

        with tempfile.TemporaryDirectory() as d:
            fasta = Path(d) / "silva.fasta"
            fasta.write_text(">A.1.2 Bacteria;Proteobacteria\nACGT\n", encoding="utf-8")
            path = SILVAConfig(db_path=fasta).sintax_ready_path()

        self.assertIn(SINTAX_FORMAT_VERSION, path.name)


class BlastTimeoutTests(unittest.TestCase):
    def test_a_local_blast_does_not_inherit_the_zero_remote_timeout(self):
        """NCBI_TIMEOUT_SEC defaults to 0 to keep remote opt-in; that 0 was killing local runs."""
        from core.config import BLAST_TIMEOUT_SEC
        from databases.ncbi_db import NCBIConfig

        config = NCBIConfig(local_db_path=Path("/db/nt"), use_remote=False, timeout=0)
        self.assertEqual(config.effective_timeout(), BLAST_TIMEOUT_SEC)
        self.assertGreater(config.effective_timeout(), 0)

    def test_remote_blast_is_not_gated_behind_a_timeout_threshold(self):
        """A zero/low timeout must not veto remote BLAST outright.

        The old gate refused to run whenever ``NCBI_TIMEOUT_SEC < 1800`` and told the
        user to set it — but the variable was absent from .env, so NCBI was dead on
        arrival on every run. A missing timeout means "use the default", not "refuse".
        """
        from databases.ncbi_db import NCBIConfig

        config = NCBIConfig(use_remote=True, timeout=0)
        self.assertEqual(config.validate(), [])
        self.assertGreater(config.effective_timeout(), 0)

    def test_an_explicit_timeout_is_honoured(self):
        from databases.ncbi_db import NCBIConfig

        self.assertEqual(NCBIConfig(local_db_path=Path("/db/nt"), timeout=120).effective_timeout(), 120)


class SintaxReliabilityTests(unittest.TestCase):
    """VSEARCH already dropped everything below the cutoff, so grading against the
    cutoff marked every surviving call "high" — a flag true by construction."""

    def test_the_grade_varies_with_the_bootstrap(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertEqual(assignment_reliability({"confidence": "0.99"}, "genus"), "high")
        self.assertEqual(assignment_reliability({"confidence": "0.90"}, "genus"), "moderate")
        self.assertEqual(assignment_reliability({"confidence": "0.82"}, "genus"), "low")

    def test_a_call_that_only_just_cleared_the_cutoff_is_not_high(self):
        from databases.taxonomy_runner import assignment_reliability

        self.assertNotEqual(assignment_reliability({"confidence": "0.81"}, "species"), "high")


class SilvaMemoryGuardTests(unittest.TestCase):
    def test_a_reference_bigger_than_ram_is_skipped_not_run(self):
        """Running it anyway means a SIGKILL minutes later and an empty result file."""
        from databases import taxonomy_runner

        with tempfile.TemporaryDirectory() as d:
            db = Path(d) / "silva.fasta"
            db.write_text("x" * 1000, encoding="utf-8")
            with patch.object(taxonomy_runner, "available_memory_bytes", lambda: 500):
                message = taxonomy_runner._sintax_memory_shortfall(db)

        self.assertIsNotNone(message)
        self.assertIn("RAM", message)

    def test_a_reference_that_fits_is_allowed_to_run(self):
        from databases import taxonomy_runner

        with tempfile.TemporaryDirectory() as d:
            db = Path(d) / "pr2.fasta"
            db.write_text("x" * 100, encoding="utf-8")
            with patch.object(taxonomy_runner, "available_memory_bytes", lambda: 10_000):
                self.assertIsNone(taxonomy_runner._sintax_memory_shortfall(db))

    def test_requirement_scales_past_the_raw_file_size(self):
        """SINTAX needs a MULTIPLE of the FASTA, not its size.

        Comparing against the raw file size let a 786 MB SILVA start on a box with
        1.0 GB free; VSEARCH then built its k-mer index and was OOM-killed partway
        through, which surfaced to the user as an unexplained failure.
        """
        from databases import taxonomy_runner

        with tempfile.TemporaryDirectory() as d:
            db = Path(d) / "silva.fasta"
            db.write_text("x" * 1000, encoding="utf-8")
            # Fits on disk-size alone, but not once the index is accounted for.
            with patch.object(taxonomy_runner, "available_memory_bytes", lambda: 1500):
                self.assertIsNotNone(taxonomy_runner._sintax_memory_shortfall(db))


class ClusteringTests(unittest.TestCase):
    def _commands(self, params=None):
        """Capture the VSEARCH command lines the pooling step would run."""
        from core.models import AnalysisBranch, PhyloMode, RunConfig
        from analysis import vsearch_runner

        captured = []

        def fake_run(command, **kwargs):
            captured.append(command)

            class Completed:
                stdout = ""
                stderr = ""

            return Completed()

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            filtered = root / "s1.fasta"
            filtered.write_text(">r1;sample=s1\nACGT\n", encoding="utf-8")
            config = RunConfig(
                run_id="r", samples=[], branch=AnalysisBranch.VSEARCH, phylo_mode=PhyloMode.FAST, output_dir=root
            )
            if params:
                config.vsearch_params.update(params)
            with patch.object(vsearch_runner, "_run", side_effect=fake_run), \
                patch.object(vsearch_runner, "_count_fasta_seqs", return_value=1):
                vsearch_runner.run_vsearch_pool_and_cluster({"s1": filtered}, config, root / "pooled")
        return captured

    def test_clustering_carries_abundances_in_and_out(self):
        """Without --sizein the chimera step judges abundance skew on counts it cannot see."""
        cluster = next(command for command in self._commands() if "--cluster_size" in command)
        self.assertIn("--sizein", cluster)
        self.assertIn("--sizeout", cluster)

    def test_sintax_is_seeded_so_two_runs_of_the_same_data_agree(self):
        """vsearch --randseed defaults to 0 = "seed from a random data source", so SINTAX's
        bootstrap subsampling was unseeded: identical input gave SILVA 0.99 on one run and
        0.98 on the next. Those numbers set `Reliability` and are compared against
        --sintax_cutoff, so a taxon could appear or vanish between two runs of the same
        data, and no reviewer could reproduce a published number."""
        import subprocess
        from pathlib import Path
        from unittest.mock import patch

        import tempfile as _tempfile

        captured = {}

        def fake_run(command, **kwargs):
            captured["command"] = command
            return subprocess.CompletedProcess(args=command, returncode=0, stdout="", stderr="")

        with _tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "q.fasta").write_text(">OTU1\nACGT\n", encoding="utf-8")
            (root / "db.fasta").write_text(">R1;tax=d:Eukaryota\nACGT\n", encoding="utf-8")
            with patch("databases.taxonomy_runner.subprocess.run", side_effect=fake_run):
                from databases.taxonomy_runner import run_sintax

                run_sintax(root / "q.fasta", root / "db.fasta", 0.8, 1, output_dir=root)

        command = captured["command"]
        self.assertIn("--randseed", command)
        seed = command[command.index("--randseed") + 1]
        self.assertNotEqual(seed, "0", "randseed 0 means 'use a random source' — the bug")

    def test_conflict_report_row_order_does_not_depend_on_source_order(self):
        """Two runs over identical data must emit identical bytes. This walked the dict
        in insertion order — i.e. the order the sources happened to answer in — so the
        same three rows came out shuffled between runs, and a file that changes without
        its content changing cannot be diffed against a reference."""
        from databases.taxonomy_runner import _build_conflict_report

        def records(genus_a, genus_b):
            return [{"source": "NCBI", "genus": genus_a}, {"source": "SILVA", "genus": genus_b}]

        forward = {"OTU1": records("Gadus", "Melanogrammus"), "OTU3": records("Gadus", "Pollachius")}
        reversed_insertion = {"OTU3": records("Gadus", "Pollachius"), "OTU1": records("Gadus", "Melanogrammus")}

        self.assertEqual(
            [row["feature_id"] for row in _build_conflict_report(forward)],
            [row["feature_id"] for row in _build_conflict_report(reversed_insertion)],
        )

    def test_singletons_are_retained_by_default(self):
        """The default is to let a singleton seed an OTU.

        This is a parasite-detection platform: the read a floor of 2 uniquely discards is
        the globally-unique, >3%-divergent one — a rare novel taxon, not an error (errors
        sit within 3% of their parent and are absorbed by clustering either way). Losing
        it is the most expensive outcome the pipeline has. Measured on the reference run,
        floors of 1 and 2 give identical output, so the floor buys nothing but time.
        """
        derep = next(command for command in self._commands() if "--derep_fulllength" in command)
        self.assertNotIn("--minuniquesize", derep)

    def test_singletons_can_be_dropped_explicitly(self):
        """A floor of 2 stays available — for a large host-dominated study buying back
        clustering time. Its cost is no longer silent: it shows up as `unmapped_reads`."""
        derep = next(command for command in self._commands({"min_unique_size": 2}) if "--derep_fulllength" in command)
        self.assertIn("--minuniquesize", derep)
        self.assertEqual(derep[derep.index("--minuniquesize") + 1], "2")


class FeatureNamingTests(unittest.TestCase):
    """VSEARCH clusters are OTUs; DADA2 variants are ASVs. They are different objects."""

    def _tips(self, prefix):
        from phylo.fasta_utils import sanitize_fasta_for_phylo

        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "in.fasta"
            src.write_text(">long_header_1\n" + "A" * 120 + "\n", encoding="utf-8")
            out = Path(d) / "out.fasta"
            mapping = sanitize_fasta_for_phylo(src, out, prefix=prefix)
        return list(mapping)

    def test_vsearch_features_are_named_otu(self):
        self.assertEqual(self._tips("OTU"), ["OTU1"])

    def test_dada2_features_are_named_asv(self):
        self.assertEqual(self._tips("ASV"), ["ASV1"])

    def test_the_phylo_runner_picks_the_prefix_from_the_analysis_method(self):
        from core.models import AnalysisBranch, PhyloMode, RunConfig
        from phylo import phylo_runner

        captured = {}

        def fake_mafft(fasta_path, output_dir, params):
            captured["prefix"] = params.get("feature_prefix")
            raise RuntimeError("stop here")

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / "in.fasta"
            fasta.write_text(">a\nACGT\n", encoding="utf-8")
            for method, expected in (("vsearch", "OTU"), ("dada2", "ASV")):
                config = RunConfig(
                    run_id="r", samples=[], branch=AnalysisBranch.VSEARCH, phylo_mode=PhyloMode.FAST, output_dir=root
                )
                config.analysis_method = method
                with patch.object(phylo_runner, "run_mafft", side_effect=fake_mafft):
                    with self.assertRaises(RuntimeError):
                        phylo_runner.run_phylo(fasta, config, root / "phylo")
                self.assertEqual(captured["prefix"], expected)


if __name__ == "__main__":
    unittest.main()


class PrimerRemovalTests(unittest.TestCase):
    """ILLUMINACLIP strips sequencing adapters, not PCR primers. Left in, the
    primer's degenerate positions reach the classifier and the tree as signal."""

    def test_a_preset_yields_the_primer_lengths(self):
        from analysis.vsearch_runner import PRIMER_PRESETS, _primer_strip_lengths

        forward, reverse = PRIMER_PRESETS["18S_V4_TAReuk"]
        self.assertEqual(_primer_strip_lengths({"primer_preset": "18S_V4_TAReuk"}), (len(forward), len(reverse)))

    def test_explicit_sequences_win_over_a_preset(self):
        from analysis.vsearch_runner import _primer_strip_lengths

        self.assertEqual(_primer_strip_lengths({"primer_fwd": "ACGTACGT", "primer_rev": "ACGT"}), (8, 4))

    def test_explicit_lengths_are_used_verbatim(self):
        from analysis.vsearch_runner import _primer_strip_lengths

        self.assertEqual(_primer_strip_lengths({"primer_strip_left": 21, "primer_strip_right": 17}), (21, 17))

    def test_nothing_is_stripped_when_no_primers_are_configured(self):
        """Guessing a length would silently truncate real sequence."""
        from analysis.vsearch_runner import _primer_strip_lengths

        self.assertEqual(_primer_strip_lengths({}), (0, 0))


class OtuIdTests(unittest.TestCase):
    def test_the_chimera_step_relabels_and_keeps_the_original_header(self):
        from analysis import vsearch_runner

        captured = []
        with patch.object(vsearch_runner, "_run", side_effect=lambda command, **kw: captured.append(command) or type("C", (), {"stdout": "", "stderr": ""})()), \
            patch.object(vsearch_runner, "_count_fasta_seqs", return_value=1), \
            patch.object(vsearch_runner, "_write_otu_id_map", lambda *a: Path("/x.csv")):
            from core.models import AnalysisBranch, PhyloMode, RunConfig

            with tempfile.TemporaryDirectory() as d:
                root = Path(d)
                filtered = root / "s1.fasta"
                filtered.write_text(">r1;sample=s1\nACGT\n", encoding="utf-8")
                config = RunConfig(run_id="r", samples=[], branch=AnalysisBranch.VSEARCH, phylo_mode=PhyloMode.FAST, output_dir=root)
                vsearch_runner.run_vsearch_pool_and_cluster({"s1": filtered}, config, root / "pooled")

        uchime = next(command for command in captured if "--uchime_denovo" in command)
        self.assertIn("--relabel", uchime)
        self.assertEqual(uchime[uchime.index("--relabel") + 1], "OTU")
        self.assertIn("--relabel_keep", uchime)

    def _id_map_fixture(self, root: Path) -> tuple[Path, Path]:
        otus = root / "otus.fasta"
        otus.write_text(
            ">OTU1;size=100 A01426:1:HW7:2:2101:6289:3458;sample=s1;size=100\nACGT\n"
            ">OTU2;size=14144 A01426:1:HW7:2:2101:23375:23531;sample=s1;size=14144\nTTTT\n",
            encoding="utf-8",
        )
        otu_table = root / "otu_table.tsv"
        # What the mapping actually found: OTU2 gathers far more reads than it seeded
        # with, which is the divergence this map used to hide.
        otu_table.write_text("#OTU ID\ts1\ts2\nOTU1\t60\t34\nOTU2\t1244\t62016\n", encoding="utf-8")
        return otus, otu_table

    def test_the_id_map_ties_each_otu_to_its_representative_read(self):
        import csv as _csv

        from analysis.vsearch_runner import _write_otu_id_map

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            otus, otu_table = self._id_map_fixture(root)
            rows = list(_csv.DictReader(open(_write_otu_id_map(otus, otu_table, root / "map.csv"))))

        self.assertEqual(rows[0]["otu_id"], "OTU1")
        self.assertIn("A01426", rows[0]["representative_read"])

    def test_id_map_size_is_the_mapping_total_not_the_cluster_size(self):
        """The bug: ``size`` came from the FASTA's ``;size=``, i.e. the clustering.

        On the reference run that put OTU2 at 14 144 while otu_table said 63 260, and
        made the run total 147 010 here and 157 527 everywhere else. ``size`` must be
        the mapping total; the clustering number keeps its own column.
        """
        import csv as _csv

        from analysis.vsearch_runner import _write_otu_id_map

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            otus, otu_table = self._id_map_fixture(root)
            rows = {r["otu_id"]: r for r in _csv.DictReader(open(_write_otu_id_map(otus, otu_table, root / "map.csv")))}

        self.assertEqual(rows["OTU2"]["size"], "63260")  # 1244 + 62016, from otu_table
        self.assertEqual(rows["OTU2"]["cluster_size"], "14144")  # preserved, not canonical
        self.assertEqual(rows["OTU1"]["size"], "94")
        # Rule (b): every size equals its otu_table row sum, so the file totals agree.
        self.assertEqual(sum(int(r["size"]) for r in rows.values()), 63354)

    def test_an_otu_absent_from_the_mapping_table_gets_zero_not_blank(self):
        import csv as _csv

        from analysis.vsearch_runner import _write_otu_id_map

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            otus, otu_table = self._id_map_fixture(root)
            otu_table.write_text("#OTU ID\ts1\nOTU1\t94\n", encoding="utf-8")
            rows = {r["otu_id"]: r for r in _csv.DictReader(open(_write_otu_id_map(otus, otu_table, root / "map.csv")))}

        self.assertEqual(rows["OTU2"]["size"], "0")

    def test_chimeric_reads_are_matched_on_the_original_label(self):
        """--relabel restarts numbering in the chimeras file, so its "OTU1" is a
        different sequence from the non-chimera "OTU1" — and neither is in the .uc."""
        from analysis.analysis_metrics import chimera_stats, original_label

        self.assertEqual(original_label("OTU1;size=3 C1;size=3"), "C1;size=3")

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "otus_raw.fasta").write_text(">P1;size=500\nAC\n>C1;size=3\nGT\n", encoding="utf-8")
            (root / "chimeras.fasta").write_text(">OTU1;size=3 C1;size=3\nGT\n", encoding="utf-8")
            (root / "clusters.uc").write_text(
                "S\t0\t4\t*\t*\t*\t*\t*\tP1;size=500\t*\n"
                "H\t0\t4\t99.0\t+\t0\t0\t4M\tm1;size=497\tP1;size=500\n"
                "S\t1\t4\t*\t*\t*\t*\t*\tC1;size=3\t*\n",
                encoding="utf-8",
            )
            stats = chimera_stats(root / "chimeras.fasta", root / "otus_raw.fasta", root / "clusters.uc")

        self.assertEqual(stats["chimeric_reads"], 3)
        self.assertEqual(stats["read_basis"], "cluster")


class ExportTests(unittest.TestCase):
    def _fixture(self, root: Path):
        matrix = root / "m.csv"
        matrix.write_text("#OTU ID,s1,s2\nOTU1,80,0\nOTU2,20,50\n", encoding="utf-8")
        consensus = root / "c.csv"
        consensus.write_text(
            "ASV_ID,Consensus_Kingdom,Consensus_Phylum,Consensus_Class,Consensus_Order,Consensus_Family,Consensus_Genus,Consensus_Species\n"
            "OTU1,Eukaryota,Platyhelminthes,Trematoda,D,Dp,Diplostomum,Diplostomum spathaceum\n"
            "OTU2,Eukaryota,Chordata,Actinopteri,,,,\n",
            encoding="utf-8",
        )
        return matrix, consensus

    def test_a_table_is_written_for_every_rank(self):
        from databases.database_metrics import taxa_tables_all_ranks

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            matrix, consensus = self._fixture(root)
            written = taxa_tables_all_ranks(matrix, consensus, root / "out")

        for rank in ("kingdom", "phylum", "class", "order", "family", "genus", "species"):
            self.assertIn(f"taxa_table_{rank}_csv", written)
            self.assertIn(f"taxa_table_{rank}_relative_csv", written)

    def test_the_relative_table_sums_to_one_hundred_per_sample(self):
        import csv as _csv

        from databases.database_metrics import taxa_tables_all_ranks

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            matrix, consensus = self._fixture(root)
            written = taxa_tables_all_ranks(matrix, consensus, root / "out")
            rows = list(_csv.reader(open(written["taxa_table_phylum_relative_csv"])))

        header, *body = rows
        for column in range(1, len(header)):
            self.assertAlmostEqual(sum(float(row[column]) for row in body), 100.0, places=3)

    def test_biom_carries_a_taxonomy_lineage_per_feature(self):
        import json

        from databases.database_metrics import write_biom

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            matrix, consensus = self._fixture(root)
            table = json.loads(Path(write_biom(matrix, consensus, root / "t.biom")).read_text())

        self.assertEqual(table["shape"], [2, 2])
        self.assertEqual(table["matrix_type"], "sparse")
        lineage = {row["id"]: row["metadata"]["taxonomy"] for row in table["rows"]}
        self.assertEqual(lineage["OTU1"][-1], "s__Diplostomum spathaceum")
        self.assertEqual(lineage["OTU2"][-1], "s__")  # unassigned, not fabricated
        # Sparse: only the non-zero cells.
        self.assertEqual(len(table["data"]), 3)


class SintaxParentheticalNameTests(unittest.TestCase):
    """SILVA puts common names in parentheses; SINTAX puts bootstraps there too.

    7.7% of SILVA's eukaryote species names carry a parenthetical
    ("Gadus_morhua_(Atlantic_cod)"). The parser stripped whatever was in the last
    pair of parentheses, so on the cutoff-filtered column — which has no bootstraps
    at all — it amputated the name and promoted the common name to a bootstrap. It
    only became reachable once marker-filtered references let SILVA run at all.
    """

    RAW = "d:Eukaryota(1.00),g:Teleostei(1.00),s:Gadus_morhua_(Atlantic_cod)(0.97)"
    FILTERED = "d:Eukaryota,g:Teleostei,s:Gadus_morhua_(Atlantic_cod)"

    def test_a_common_name_is_not_mistaken_for_a_bootstrap(self):
        from databases.taxonomy_runner import _parse_sintax_assignment

        ranks, confidence = _parse_sintax_assignment(self.FILTERED)
        self.assertEqual(ranks["species"], "Gadus_morhua_(Atlantic_cod)")
        self.assertEqual(confidence, "")

    def test_a_real_bootstrap_is_still_read_past_a_parenthetical_name(self):
        from databases.taxonomy_runner import _parse_sintax_assignment

        ranks, confidence = _parse_sintax_assignment(self.RAW)
        self.assertEqual(ranks["species"], "Gadus_morhua_(Atlantic_cod)")
        self.assertEqual(confidence, "0.97")

    def test_names_without_parentheses_are_untouched(self):
        from databases.taxonomy_runner import _parse_sintax_assignment

        ranks, confidence = _parse_sintax_assignment("d:Eukaryota(1.00),s:Diplostomum_spathaceum(0.64)")
        self.assertEqual(ranks["species"], "Diplostomum_spathaceum")
        self.assertEqual(confidence, "0.64")

    def test_rank_bootstrap_lookup_skips_the_common_name(self):
        from databases.taxonomy_runner import _sintax_confidence_for

        self.assertEqual(_sintax_confidence_for(self.RAW, "species"), "0.97")
        # No numeric bootstrap present: report none rather than the word.
        self.assertEqual(_sintax_confidence_for(self.FILTERED, "species"), "")


class RankVocabularyTests(unittest.TestCase):
    """One rank spelled two ways in two modules, and the lookup that silently missed."""

    def test_collapsing_by_class_reads_the_class_column(self):
        # database_metrics.TAXONOMY_RANKS says "class_"; community_analysis keyed on
        # "class". The get returned {} and every read fell into "Unclassified" — while
        # tax_table.tsv, which uses a different lookup, showed the real class.
        from analysis.community_analysis import collapsed_taxa_table

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            matrix = root / "m.csv"
            matrix.write_text("#OTU ID,s1\nOTU1,100\nOTU2,50\n", encoding="utf-8")
            consensus = root / "c.csv"
            consensus.write_text(
                "ASV_ID,Consensus_Kingdom,Consensus_Phylum,Consensus_Class,Consensus_Order,"
                "Consensus_Family,Consensus_Genus,Consensus_Species\n"
                "OTU1,Eukaryota,Chordata,Actinopteri,,,,\n"
                "OTU2,Eukaryota,Chordata,Actinopteri,,,,\n",
                encoding="utf-8",
            )
            out = collapsed_taxa_table(matrix, consensus, "class_", root / "t.csv")
            body = Path(out).read_text(encoding="utf-8").strip().splitlines()

        self.assertEqual(body[1], "Actinopteri,150")
        self.assertNotIn("Unclassified", "\n".join(body))

    def test_both_spellings_of_the_class_rank_agree(self):
        from analysis.community_analysis import canonical_rank

        self.assertEqual(canonical_rank("class_"), "class")
        self.assertEqual(canonical_rank("class"), "class")
        self.assertEqual(canonical_rank("genus"), "genus")


class BackboneHarmonisationTests(unittest.TestCase):
    """Neither PR2 nor SILVA declares Linnean middle ranks. Neither may invent them."""

    def test_pr2_does_not_map_division_onto_phylum_or_its_class_onto_class(self):
        from core.utils import PR2_RANK_LABELS, taxonomy_fields_to_sintax

        # Domain;Supergroup;Division;Subdivision;Class;Order;Family;Genus;Species
        fish = ["Eukaryota", "Obazoa", "Opisthokonta", "Metazoa", "Craniata", "Craniata_X", "Teleostei", "Auxis", "Auxis_thazard"]
        tax = taxonomy_fields_to_sintax(fish, PR2_RANK_LABELS)

        self.assertNotIn("p:Opisthokonta", tax)  # a division is not a phylum
        self.assertNotIn("c:Craniata", tax)  # a subphylum is not a class
        self.assertIn("d:Eukaryota", tax)
        self.assertIn("g:Auxis", tax)
        self.assertIn("s:Auxis_thazard", tax)

    def test_silva_claims_only_the_domain_and_the_organism(self):
        from core.utils import silva_rank_labels, taxonomy_fields_to_sintax

        # SILVA declares no ranks; the old code claimed p/c/o/f then g on the
        # last-but-one, which made an infraclass (Teleostei) a genus.
        silva = ["Eukaryota", "Amorphea", "Obazoa", "Opisthokonta", "Holozoa", "Teleostei", "Gadus_morhua"]
        tax = taxonomy_fields_to_sintax(silva, silva_rank_labels(len(silva)))

        self.assertEqual(tax, "d:Eukaryota,s:Gadus_morhua")
        for invented in ("p:Amorphea", "c:Obazoa", "o:Opisthokonta", "f:Holozoa", "g:Teleostei"):
            self.assertNotIn(invented, tax)

    def test_silva_rank_labels_never_truncates_a_lineage(self):
        from core.utils import silva_rank_labels

        for depth in range(0, 12):
            self.assertEqual(len(silva_rank_labels(depth)), depth)

    def test_pr2_placeholders_are_not_taxa(self):
        from core.utils import is_pr2_placeholder

        for placeholder in ("Craniata_X", "Craniata_XX", "Craniata_XXX", "Craniata_XXX_sp.", "Chlorophyta_sp"):
            self.assertTrue(is_pr2_placeholder(placeholder), placeholder)
        for real in ("Craniata", "Gadus_morhua", "Dinophyceae", "Unruhdinium_kevei"):
            self.assertFalse(is_pr2_placeholder(real), real)

    def test_blanking_a_placeholder_does_not_shift_the_ranks_below_it(self):
        # PR2's nine levels are positional: dropping a field instead of blanking it
        # would slide Genus onto Family and Species onto Genus.
        from core.utils import PR2_RANK_LABELS, is_pr2_placeholder, taxonomy_fields_to_sintax

        raw = ["Eukaryota", "TSAR", "Alveolata", "Dinoflagellata", "Dinophyceae", "Peridiniales", "Peridiniales_X", "Unruhdinium", "Unruhdinium_kevei"]
        fields = ["" if is_pr2_placeholder(f) else f for f in raw]
        tax = taxonomy_fields_to_sintax(fields, PR2_RANK_LABELS)

        self.assertIn("g:Unruhdinium", tax)
        self.assertIn("s:Unruhdinium_kevei", tax)
        self.assertNotIn("Peridiniales_X", tax)
