"""SILVA's declared ranks (tax_slv rank file) reach the SINTAX reference, opt-in only.

Without the rank file SILVA contributes domain + organism name (``silva_rank_labels``),
and that behaviour, cache name included, must not move: every stored 18S result was made
with it. With ``DB_SILVA_RANKS_PATH`` the intermediate levels take the rank SILVA itself
declares for their path, under a separate cache name keyed by the rank file's content.
"""

from __future__ import annotations

import csv
import gzip
import os
import random
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from core.config import get_db_paths
from core.utils import (
	load_silva_rank_map,
	silva_declared_labels,
	silva_rank_labels,
	taxonomy_fields_to_sintax,
)
from databases.silva_db import SILVAConfig, _convert_silva_to_sintax, silva_rank_path


# A miniature tax_slv file: the real format is path<TAB>taxid<TAB>rank<TAB>remark<TAB>release.
# The ranks are those the real file uses, including the non-Linnean ones to be ignored.
RANK_LINES = [
	("Eukaryota;", "2", "domain"),
	("Eukaryota;Amorphea;", "10", "major_clade"),
	("Eukaryota;Amorphea;Obazoa;", "11", "major_clade"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;", "12", "major_clade"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);", "13", "kingdom"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Vertebrata;", "14", "phylum"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Vertebrata;Actinopterygii;", "15", "class"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Vertebrata;Actinopterygii;Gadiformes;", "16", "order"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Vertebrata;Actinopterygii;Gadiformes;Gadidae;", "17", "family"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Vertebrata;Actinopterygii;Gadiformes;Gadidae;Gadus;", "18", "genus"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Nematozoa;", "20", "phylum"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Nematozoa;Chromadorea;", "21", "class"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Nematozoa;Chromadorea;Ascaridida;", "22", "order"),
	("Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Nematozoa;Chromadorea;Ascaridida;Anisakidae;", "23", "family"),
	("Eukaryota;Archaeplastida;", "30", "major_clade"),
	("Eukaryota;Archaeplastida;Chloroplastida;", "31", "phylum"),
	("Eukaryota;Archaeplastida;Chloroplastida;Charophyta;", "32", "class"),
]

FISH = "Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Vertebrata;Actinopterygii;Gadiformes;Gadidae;Gadus;Gadus morhua (Atlantic cod)"
WORM = "Eukaryota;Amorphea;Obazoa;Opisthokonta;Metazoa (Animalia);Nematozoa;Chromadorea;Ascaridida;Anisakidae;Anisakis;Anisakis simplex"
PLANT = "Eukaryota;Archaeplastida;Chloroplastida;Charophyta;Phragmoplastophyta;Pinus taeda (loblolly pine)"


def _write_rank_file(path: Path, lines=RANK_LINES, gz: bool = False) -> Path:
	text = "".join(f"{lineage}\t{taxid}\t{rank}\t\t138.2\n" for lineage, taxid, rank in lines)
	if gz:
		with gzip.open(path, "wt", encoding="utf-8") as handle:
			handle.write(text)
	else:
		path.write_text(text, encoding="utf-8")
	return path


def _random_sequence(rng: random.Random, length: int = 300) -> str:
	return "".join(rng.choice("ACGT") for _ in range(length))


def _mutate(rng: random.Random, sequence: str, n: int) -> str:
	seq = list(sequence)
	for position in rng.sample(range(len(seq)), n):
		seq[position] = rng.choice([b for b in "ACGT" if b != seq[position]])
	return "".join(seq)


def _write_silva_fasta(path: Path, rng: random.Random) -> dict:
	"""Three taxa, five close copies each, SILVA-style headers. Returns one query per taxon."""
	queries = {}
	with path.open("w", encoding="utf-8") as handle:
		for name, lineage in (("fish", FISH), ("worm", WORM), ("plant", PLANT)):
			base = _random_sequence(rng)
			queries[name] = base
			for copy in range(5):
				handle.write(f">{name.upper()}{copy}.1.300 {lineage}\n{_mutate(rng, base, 3)}\n")
	return queries


class RankFileParsingTests(unittest.TestCase):
	def setUp(self):
		self.tmp = Path(tempfile.mkdtemp())

	def tearDown(self):
		shutil.rmtree(self.tmp, ignore_errors=True)

	def test_reads_plain_and_gzipped_files_identically(self):
		plain = load_silva_rank_map(_write_rank_file(self.tmp / "tax.txt"))
		gzipped = load_silva_rank_map(_write_rank_file(self.tmp / "tax.txt.gz", gz=True))
		self.assertEqual(plain, gzipped)
		self.assertEqual(plain["Eukaryota;"], "domain")
		self.assertEqual(len(plain), len(RANK_LINES))

	def test_skips_malformed_lines_and_terminates_paths(self):
		path = self.tmp / "tax.txt"
		path.write_text("garbage line\n\nEukaryota\t2\tDomain\t\t138.2\n\t3\tphylum\n", encoding="utf-8")
		self.assertEqual(load_silva_rank_map(path), {"Eukaryota;": "domain"})


class DeclaredLabelsTests(unittest.TestCase):
	def setUp(self):
		tmp = Path(tempfile.mkdtemp())
		self.addCleanup(shutil.rmtree, tmp, True)
		self.rank_map = load_silva_rank_map(_write_rank_file(tmp / "tax.txt"))

	def _sintax(self, lineage: str) -> str:
		fields = lineage.split(";")
		return taxonomy_fields_to_sintax(fields, silva_declared_labels(fields, self.rank_map))

	def test_fish_gets_every_declared_rank_and_no_clade(self):
		self.assertEqual(
			self._sintax(FISH),
			"d:Eukaryota,p:Vertebrata,c:Actinopterygii,o:Gadiformes,f:Gadidae,g:Gadus,s:Gadus_morhua_(Atlantic_cod)",
		)

	def test_names_the_verifier_forbids_are_never_given_a_rank(self):
		# major_clade and kingdom have no Linnean column: Amorphea/Obazoa/Opisthokonta and
		# "Metazoa" must not land in p/c/o/f/g (scripts/verify.py check d).
		sintax = self._sintax(FISH)
		for name in ("Amorphea", "Obazoa", "Opisthokonta", "Metazoa"):
			self.assertNotIn(name, sintax)

	def test_undeclared_genus_stays_empty(self):
		# Anisakis has no line in the rank file: no genus is invented for it.
		self.assertEqual(
			self._sintax(WORM),
			"d:Eukaryota,p:Nematozoa,c:Chromadorea,o:Ascaridida,f:Anisakidae,s:Anisakis_simplex",
		)

	def test_rank_order_is_never_inverted(self):
		rank_map = {"A;": "domain", "A;B;": "family", "A;B;C;": "order", "A;B;C;D;": "genus"}
		labels = silva_declared_labels(["A", "B", "C", "D", "E sp."], rank_map)
		self.assertEqual(labels, ["d", "f", None, "g", "s"])

	def test_same_name_on_another_branch_is_looked_up_on_its_own_path(self):
		rank_map = {"A;": "domain", "A;X;": "phylum", "A;Y;": "major_clade", "A;Y;X;": "class"}
		self.assertEqual(silva_declared_labels(["A", "Y", "X", "sp"], rank_map), ["d", None, "c", "s"])

	def test_length_always_matches_depth(self):
		for depth in range(0, 12):
			fields = FISH.split(";")[:depth]
			self.assertEqual(len(silva_declared_labels(fields, self.rank_map)), depth)
		self.assertEqual(silva_declared_labels(["Eukaryota"], self.rank_map), ["d"])
		self.assertEqual(silva_declared_labels(["Eukaryota", "x"], self.rank_map), ["d", "s"])

	def test_empty_map_gives_exactly_the_unranked_labels(self):
		for lineage in (FISH, WORM, PLANT):
			fields = lineage.split(";")
			self.assertEqual(silva_declared_labels(fields, {}), silva_rank_labels(len(fields)))


class ReferenceCacheTests(unittest.TestCase):
	def setUp(self):
		self.tmp = Path(tempfile.mkdtemp())
		self.addCleanup(shutil.rmtree, self.tmp, True)
		self.silva = self.tmp / "SILVA_test_tax_silva.fasta"
		_write_silva_fasta(self.silva, random.Random(7))
		self.ranks = _write_rank_file(self.tmp / "tax_slv_ssu_test.txt.gz", gz=True)

	def test_without_rank_file_the_cache_name_and_content_are_unchanged(self):
		path = SILVAConfig(db_path=self.silva).sintax_ready_path()
		self.assertEqual(path.name, "SILVA_test_tax_silva.sintax.v4.fasta")
		headers = [line for line in path.read_text().splitlines() if line.startswith(">")]
		self.assertIn(">FISH0.1.300;tax=d:Eukaryota,s:Gadus_morhua_(Atlantic_cod);", headers)
		self.assertTrue(all(",p:" not in h and ",g:" not in h for h in headers))

	def test_ranked_cache_has_its_own_name_and_leaves_the_unranked_one_alone(self):
		unranked = SILVAConfig(db_path=self.silva).sintax_ready_path()
		before = (unranked.read_bytes(), unranked.stat().st_mtime_ns)
		ranked = SILVAConfig(db_path=self.silva, rank_path=self.ranks).sintax_ready_path()
		self.assertNotEqual(ranked, unranked)
		self.assertRegex(ranked.name, r"^SILVA_test_tax_silva\.sintax\.v4r\.[0-9a-f]{8}\.fasta$")
		self.assertEqual((unranked.read_bytes(), unranked.stat().st_mtime_ns), before)
		text = ranked.read_text()
		self.assertIn(">FISH0.1.300;tax=d:Eukaryota,p:Vertebrata,c:Actinopterygii,o:Gadiformes,f:Gadidae,g:Gadus,s:Gadus_morhua_(Atlantic_cod);", text)
		# Sequences are untouched: only headers change.
		self.assertEqual(
			[l for l in text.splitlines() if not l.startswith(">")],
			[l for l in self.silva.read_text().splitlines() if not l.startswith(">")],
		)

	def test_a_different_rank_file_never_reuses_another_files_cache(self):
		first = SILVAConfig(db_path=self.silva, rank_path=self.ranks).sintax_ready_path()
		other = _write_rank_file(self.tmp / "other.txt", lines=RANK_LINES[:6])
		second = SILVAConfig(db_path=self.silva, rank_path=other).sintax_ready_path()
		self.assertNotEqual(first.name, second.name)

	def test_a_rank_file_from_another_release_is_refused_and_leaves_no_cache(self):
		wrong = _write_rank_file(self.tmp / "wrong.txt", lines=[("Bacteria;", "3", "domain"), ("Bacteria;Firmicutes;", "4", "phylum")])
		with self.assertRaises(ValueError):
			SILVAConfig(db_path=self.silva, rank_path=wrong).sintax_ready_path()
		self.assertEqual(sorted(p.name for p in self.tmp.glob("*.v4r.*")), [])

	def test_a_file_that_is_not_a_rank_file_is_refused(self):
		junk = self.tmp / "junk.txt"
		junk.write_text("not a rank file\n", encoding="utf-8")
		with self.assertRaises(ValueError):
			SILVAConfig(db_path=self.silva, rank_path=junk).sintax_ready_path()

	def test_convert_without_map_matches_the_previous_converter_line_for_line(self):
		dest = self.tmp / "out.fasta"
		_convert_silva_to_sintax(self.silva, dest)
		expected = []
		for line in self.silva.read_text().splitlines():
			if line.startswith(">"):
				accession, _, taxonomy = line[1:].partition(" ")
				fields = [f for f in taxonomy.split(";") if f.strip()]
				expected.append(f">{accession};tax={taxonomy_fields_to_sintax(fields, silva_rank_labels(len(fields)))};")
			else:
				expected.append(line)
		self.assertEqual(dest.read_text().splitlines(), expected)


class ConfigurationTests(unittest.TestCase):
	def test_unset_variable_means_no_rank_file(self):
		with mock.patch.dict(os.environ, {"DB_SILVA_RANKS_PATH": ""}):
			paths = get_db_paths()
		self.assertIsNone(paths.silva_ranks)
		self.assertEqual(silva_rank_path(paths), (None, None))

	def test_missing_file_is_reported_not_silently_ignored(self):
		with mock.patch.dict(os.environ, {"DB_SILVA_RANKS_PATH": "/nonexistent/tax_slv.txt.gz"}):
			paths = get_db_paths()
		path, note = silva_rank_path(paths)
		self.assertIsNone(path)
		self.assertIn("DB_SILVA_RANKS_PATH", note)

	def test_existing_file_is_used(self):
		tmp = Path(tempfile.mkdtemp())
		self.addCleanup(shutil.rmtree, tmp, True)
		ranks = _write_rank_file(tmp / "tax.txt")
		with mock.patch.dict(os.environ, {"DB_SILVA_RANKS_PATH": str(ranks)}):
			self.assertEqual(silva_rank_path(get_db_paths()), (ranks, None))


@unittest.skipUnless(shutil.which("vsearch"), "vsearch not installed")
class RankedSilvaThroughTheConsensusTests(unittest.TestCase):
	"""The whole Databases path, SILVA only: SINTAX, parsing, LCA, deliverables."""

	def _run(self, rank_path):
		from databases.taxonomy_runner import run_taxonomy_pipeline

		tmp = Path(tempfile.mkdtemp())
		self.addCleanup(shutil.rmtree, tmp, True)
		silva = tmp / "SILVA_e2e_tax_silva.fasta"
		queries = _write_silva_fasta(silva, random.Random(11))
		features = tmp / "otus.fasta"
		features.write_text("".join(f">OTU{i};size=10\n{seq}\n" for i, seq in enumerate(queries.values(), start=1)))
		ranks = _write_rank_file(tmp / "tax.txt.gz", gz=True) if rank_path else None
		run_taxonomy_pipeline(
			features,
			tmp / "out",
			silva_config=SILVAConfig(db_path=silva, rank_path=ranks),
			threads=1,
		)
		with (tmp / "out" / "taxonomy_consensus.csv").open(encoding="utf-8") as handle:
			rows = {row["feature_id"]: row for row in csv.DictReader(handle)}
		with (tmp / "out" / "source_status.csv").open(encoding="utf-8") as handle:
			status = {row["source"]: row for row in csv.DictReader(handle)}
		return rows, status

	def test_ranked_reference_fills_the_middle_ranks(self):
		rows, status = self._run(rank_path=True)
		self.assertEqual(status["SILVA"]["status"], "ok")
		self.assertIn("declared ranks from tax.txt.gz", status["SILVA"]["detail"])
		fish = rows["OTU1"]
		self.assertEqual(
			(fish["Consensus_Phylum"], fish["Consensus_Class"], fish["Consensus_Order"], fish["Consensus_Family"], fish["Consensus_Genus"]),
			("Vertebrata", "Actinopterygii", "Gadiformes", "Gadidae", "Gadus"),
		)
		self.assertEqual(rows["OTU2"]["Consensus_Family"], "Anisakidae")
		self.assertEqual(rows["OTU2"]["Consensus_Genus"], "")

	def test_default_reference_still_gives_domain_and_species_only(self):
		rows, status = self._run(rank_path=False)
		self.assertNotIn("declared ranks", status["SILVA"]["detail"])
		for row in rows.values():
			for column in ("Consensus_Phylum", "Consensus_Class", "Consensus_Order", "Consensus_Family"):
				self.assertEqual(row[column], "", column)


if __name__ == "__main__":
	unittest.main()
