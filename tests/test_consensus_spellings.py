import unittest

from databases.taxonomy_runner import _lca_consensus, _majority_consensus, _strict_consensus


def _key(value):
    return value.split("_(")[0].replace(" ", "_").lower()


SILVA = {"kingdom": "Eukaryota", "species": "Gadus_morhua_(Atlantic_cod)"}
PR2 = {"kingdom": "Eukaryota", "genus": "Gadus", "species": "Gadus_morhua"}


class TaxonKeyTests(unittest.TestCase):
    def test_common_name_and_separator_do_not_change_the_key(self):
        from databases.taxonomy_runner import _taxon_key

        self.assertEqual(_taxon_key("Gadus_morhua_(Atlantic_cod)"), "gadus_morhua")
        self.assertEqual(_taxon_key("Gadus morhua"), "gadus_morhua")
        self.assertEqual(_taxon_key("Gadus_morhua"), "gadus_morhua")

    def test_placeholders_stay_empty(self):
        from databases.taxonomy_runner import _taxon_key

        self.assertEqual(_taxon_key("unclassified"), "")
        self.assertEqual(_taxon_key("   "), "")


class ConsensusAcrossSpellingsTests(unittest.TestCase):
    def test_lca_keeps_the_species_both_references_name(self):
        consensus = _lca_consensus([PR2, SILVA])
        self.assertEqual(_key(consensus.get("species", "")), "gadus_morhua")

    def test_lca_still_falls_back_when_species_really_differ(self):
        other = {"kingdom": "Eukaryota", "genus": "Gadus", "species": "Gadus_chalcogrammus"}
        consensus = _lca_consensus([PR2, other])
        self.assertNotIn("species", consensus)
        self.assertEqual(consensus.get("genus"), "Gadus")

    def test_strict_and_majority_keep_the_same_species(self):
        self.assertEqual(_key(_strict_consensus([PR2, SILVA]).get("species", "")), "gadus_morhua")
        self.assertEqual(_key(_majority_consensus([PR2, SILVA]).get("species", "")), "gadus_morhua")


if __name__ == "__main__":
    unittest.main()
