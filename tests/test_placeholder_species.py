"""Environmental placeholders are not species.

NCBI Taxonomy files sequences such as "uncultured eukaryote" with rank "species", and
the remote GenBank search on the 18S data named one OTU that way; the consensus then
counted it as resolved to species level.
"""

import unittest


class PlaceholderNameTests(unittest.TestCase):
	def test_placeholders_are_recognised_in_ncbi_and_silva_spelling(self):
		from databases.ncbi_lineage import is_placeholder_name

		for name in ("uncultured eukaryote", "uncultured_marine_eukaryote", "Eukaryota environmental sample", "unidentified fungus"):
			self.assertTrue(is_placeholder_name(name), name)
		for name in ("Micromesistius poutassou", "Gadus_morhua_(Atlantic_cod)", "Cercospora sp.", "", None):
			self.assertFalse(is_placeholder_name(name), name)

	def test_a_placeholder_is_not_kept_as_the_species_of_a_hit(self):
		from databases.taxonomy_runner import _normalize_source_row

		row = _normalize_source_row("NCBI", {"qseqid": "OTU6", "sscinames": "uncultured eukaryote", "staxids": "100272", "pident": "99.47"})
		self.assertEqual(row["species"], "")
		self.assertEqual(row["genus"], "")

	def test_a_real_species_is_kept_and_gives_its_genus(self):
		from databases.taxonomy_runner import _normalize_source_row

		row = _normalize_source_row("NCBI", {"qseqid": "OTU1", "sscinames": "Micromesistius poutassou", "staxids": "81636", "pident": "100"})
		self.assertEqual(row["species"], "Micromesistius poutassou")
		self.assertEqual(row["genus"], "Micromesistius")

	def test_the_lineage_does_not_put_the_placeholder_back(self):
		from databases.ncbi_lineage import fill_records

		record = {"taxid": "100272", "species": ""}
		fill_records([record], {"100272": {"kingdom": "Eukaryota", "species": "uncultured eukaryote"}})
		self.assertEqual(record["kingdom"], "Eukaryota")
		self.assertEqual(record.get("species", ""), "")


if __name__ == "__main__":
	unittest.main()
