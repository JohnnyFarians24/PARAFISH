"""Remote NCBI search through the BLAST URL API, and BLAST database detection.

No test here touches the network: the URL API call, Entrez and the lineage lookup are
replaced, and the BLAST XML is a fixed sample in the format NCBI returns.
"""

import io
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

QUERY = "AGCTCCAATAGCGTATCTTAAAGTTGCTGCAGTTAAAAAGCTCGTAGTTGGATCTCGGGG"  # 60 bp
assert len(QUERY) == 60


def _hsp(score, ident, length, qfrom, qto, qseq, hseq, gaps=0):
	return f"""<Hsp>
<Hsp_num>1</Hsp_num><Hsp_bit-score>{score}</Hsp_bit-score><Hsp_score>{int(score)}</Hsp_score>
<Hsp_evalue>1e-20</Hsp_evalue><Hsp_query-from>{qfrom}</Hsp_query-from><Hsp_query-to>{qto}</Hsp_query-to>
<Hsp_hit-from>101</Hsp_hit-from><Hsp_hit-to>{100 + length}</Hsp_hit-to><Hsp_query-frame>1</Hsp_query-frame>
<Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>{ident}</Hsp_identity><Hsp_positive>{ident}</Hsp_positive>
<Hsp_gaps>{gaps}</Hsp_gaps><Hsp_align-len>{length}</Hsp_align-len>
<Hsp_qseq>{qseq}</Hsp_qseq><Hsp_hseq>{hseq}</Hsp_hseq><Hsp_midline>{"|" * len(hseq)}</Hsp_midline>
</Hsp>"""


def _hit(num, hit_id, acc, hsp):
	return f"""<Hit>
<Hit_num>{num}</Hit_num><Hit_id>{hit_id}</Hit_id><Hit_def>Some fish 18S ribosomal RNA</Hit_def>
<Hit_accession>{acc}</Hit_accession><Hit_len>1800</Hit_len><Hit_hsps>{hsp}</Hit_hsps>
</Hit>"""


def _iteration(num, query_def, qlen, hits):
	return f"""<Iteration>
<Iteration_iter-num>{num}</Iteration_iter-num><Iteration_query-ID>Query_{num}</Iteration_query-ID>
<Iteration_query-def>{query_def}</Iteration_query-def><Iteration_query-len>{qlen}</Iteration_query-len>
<Iteration_hits>{hits}</Iteration_hits>
<Iteration_stat><Statistics><Statistics_db-num>1</Statistics_db-num><Statistics_db-len>1</Statistics_db-len>
<Statistics_hsp-len>0</Statistics_hsp-len><Statistics_eff-space>0</Statistics_eff-space>
<Statistics_kappa>0.41</Statistics_kappa><Statistics_lambda>0.625</Statistics_lambda><Statistics_entropy>0.78</Statistics_entropy>
</Statistics></Iteration_stat>
</Iteration>"""


def sample_xml():
	good = _hsp(110.0, 60, 60, 1, 60, QUERY, QUERY)
	gapped = _hsp(100.0, 58, 61, 1, 60, QUERY[:30] + "-" + QUERY[30:], QUERY[:30] + "A" + QUERY[30:59] + "-", gaps=2)
	low_identity = _hsp(60.0, 50, 60, 1, 60, QUERY, QUERY)  # 83% < 90%
	short = _hsp(50.0, 30, 30, 1, 30, QUERY[:30], QUERY[:30])  # covers 50% of the query < 80%
	hits_q1 = (
		_hit(1, "gi|3053|ref|XR_014640960.1|", "XR_014640960", good)
		+ _hit(2, "gb|KF986698.1|", "KF986698", gapped)
		+ _hit(3, "gb|AAA111.1|", "AAA111", low_identity)
		+ _hit(4, "gb|BBB222.1|", "BBB222", short)
	)
	return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
<BlastOutput_program>blastn</BlastOutput_program><BlastOutput_version>BLASTN 2.17.0+</BlastOutput_version>
<BlastOutput_reference>ref</BlastOutput_reference><BlastOutput_db>core_nt</BlastOutput_db>
<BlastOutput_query-ID>Query_1</BlastOutput_query-ID><BlastOutput_query-def>OTU1;size=156789 extra</BlastOutput_query-def>
<BlastOutput_query-len>60</BlastOutput_query-len>
<BlastOutput_param><Parameters><Parameters_expect>1e-10</Parameters_expect><Parameters_sc-match>1</Parameters_sc-match>
<Parameters_sc-mismatch>-2</Parameters_sc-mismatch><Parameters_gap-open>0</Parameters_gap-open>
<Parameters_gap-extend>0</Parameters_gap-extend><Parameters_filter>L;m;</Parameters_filter></Parameters></BlastOutput_param>
<BlastOutput_iterations>
{_iteration(1, "OTU1;size=156789 extra", 60, hits_q1)}
{_iteration(2, "OTU2;size=44", 60, "")}
</BlastOutput_iterations>
</BlastOutput>"""


class URLAPIRemoteTests(unittest.TestCase):
	def _run(self, xml=None, error=None, delay=0.0, timeout=60, **cfg_kwargs):
		from databases.ncbi_db import NCBIConfig, run_blast

		def fake_qblast(sequences, config):
			self.sent.append(sequences)
			if delay:
				time.sleep(delay)
			if error:
				raise error
			return xml

		self.sent = []
		tmp = Path(tempfile.mkdtemp())
		query = tmp / "otus.fasta"
		query.write_text(f">OTU1;size=156789 extra\n{QUERY}\n>OTU2;size=44\n{QUERY}\n", encoding="utf-8")
		cfg = NCBIConfig(use_remote=True, remote_method="urlapi", timeout=timeout, **cfg_kwargs)
		taxids = {"XR_014640960.1": "8061", "XR_014640960": "8061", "KF986698.1": "8063", "KF986698": "8063"}
		with patch("databases.ncbi_db._qblast_xml", side_effect=fake_qblast), \
			patch("databases.ncbi_db._accession_taxids", return_value=taxids), \
			patch("databases.ncbi_db._taxid_names", return_value={"8061": "Micromesistius poutassou", "8063": "Merluccius merluccius"}), \
			patch("databases.ncbi_db.shutil.which", return_value=None):
			result = run_blast(query, cfg, tmp / "ncbi")
		return result, tmp / "ncbi"

	def test_the_url_api_route_does_not_need_blastn(self):
		"""shutil.which finds no blastn here, and the search still runs."""
		result, _ = self._run(xml=sample_xml())
		self.assertTrue(result.success, result.error_message)

	def test_all_otus_go_in_one_request(self):
		result, _ = self._run(xml=sample_xml())
		self.assertEqual(len(self.sent), 1)
		self.assertEqual(self.sent[0].count(">"), 2)

	def test_output_has_the_columns_of_outfmt_and_blastn_thresholds(self):
		from databases.ncbi_db import NCBIConfig

		result, out = self._run(xml=sample_xml())
		lines = (out / "blast_result.tsv").read_text(encoding="utf-8").splitlines()
		fields = NCBIConfig().outfmt.split()[1:]
		rows = [dict(zip(fields, line.split("\t"))) for line in lines]
		self.assertEqual({len(line.split("\t")) for line in lines}, {len(fields)})
		# the 83 % identity hit and the hit covering half the query are dropped, as
		# -perc_identity 90 and -qcov_hsp_perc 80 drop them in blastn
		self.assertEqual([r["sseqid"] for r in rows], ["ref|XR_014640960.1|", "gb|KF986698.1|"])
		self.assertEqual(result.metrics["hit_count"], 2.0)
		first, second = rows
		self.assertEqual(first["qseqid"], "OTU1;size=156789")
		self.assertEqual(first["pident"], "100.000")
		self.assertEqual(first["qcovhsp"], "100")
		self.assertEqual(first["staxids"], "8061")
		self.assertEqual(first["sscinames"], "Micromesistius poutassou")
		self.assertEqual(second["gapopen"], "2")
		self.assertEqual(second["mismatch"], "1")

	def test_a_query_without_hits_is_listed_as_no_hit(self):
		_, out = self._run(xml=sample_xml())
		self.assertEqual((out / "blast_no_hit.txt").read_text(encoding="utf-8").split(), ["OTU2;size=44"])

	def test_an_ncbi_error_is_a_failure_not_an_empty_result(self):
		result, _ = self._run(error=ValueError("Error message from NCBI: CPU usage limit was exceeded"))
		self.assertFalse(result.success)
		self.assertIn("refused", (result.error_message or "").lower())

	def test_no_answer_before_the_time_limit_is_a_failure(self):
		result, _ = self._run(xml=sample_xml(), delay=2.0, timeout=1)
		self.assertFalse(result.success)
		self.assertIn("did not answer", result.error_message or "")

	def test_remote_defaults(self):
		from databases.ncbi_db import NCBIConfig

		cfg = NCBIConfig(use_remote=True)
		self.assertEqual(cfg.remote_method, "urlapi")
		self.assertEqual(cfg.remote_db, "core_nt")


class AccessionTaxidTests(unittest.TestCase):
	def test_taxids_are_numbers_not_the_repr_of_entrez_objects(self):
		"""Entrez.read returns IntegerElement, whose str() is "IntegerElement(8063, ...)".
		That string reached blast_result.tsv as the taxid, the lineage lookup discarded it,
		and every remote hit was named "N/A"."""
		from databases import ncbi_db

		class IntegerElement(int):
			def __str__(self):
				return f"IntegerElement({int(self)}, attributes={{}})"

		records = [{"AccessionVersion": "KF986698.1", "TaxId": IntegerElement(8063)}]
		with patch("core.config.get_entrez_config", return_value={"email": "a@b.c", "api_key": ""}), \
			patch("Bio.Entrez.esummary", return_value=io.StringIO("")), \
			patch("Bio.Entrez.read", return_value=records):
			found = ncbi_db._accession_taxids(["KF986698"])
		self.assertEqual(found["KF986698"], "8063")
		self.assertEqual(found["KF986698.1"], "8063")

	def test_an_unnamed_hit_leaves_the_species_blank_for_the_lineage_to_fill(self):
		from databases.ncbi_db import NCBIConfig, _xml_to_rows

		rows = _xml_to_rows(sample_xml(), NCBIConfig())
		self.assertEqual({r["sscinames"] for r in rows}, {""})


class BlastDbExistsTests(unittest.TestCase):
	def test_a_database_prefix_counts_as_present(self):
		from databases.ncbi_db import blast_db_exists

		with tempfile.TemporaryDirectory() as d:
			prefix = Path(d) / "SSU_eukaryote_rRNA"
			self.assertFalse(blast_db_exists(prefix))
			Path(str(prefix) + ".nin").touch()
			self.assertTrue(blast_db_exists(prefix))
			self.assertFalse(blast_db_exists(None))


if __name__ == "__main__":
	unittest.main()
