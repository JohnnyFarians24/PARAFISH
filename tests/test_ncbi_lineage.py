"""Tests for the NCBI taxid -> lineage resolver used by the Databases stage."""

from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from databases import ncbi_lineage


def _fish_record():
    """An Entrez taxonomy record shaped like Bio.Entrez.read() output."""
    return {
        "TaxId": "58100",
        "ScientificName": "Micromesistius poutassou",
        "Rank": "species",
        "AkaTaxIds": [],
        "LineageEx": [
            {"Rank": "superkingdom", "ScientificName": "Eukaryota"},
            {"Rank": "kingdom", "ScientificName": "Metazoa"},
            {"Rank": "phylum", "ScientificName": "Chordata"},
            {"Rank": "class", "ScientificName": "Actinopteri"},
            {"Rank": "order", "ScientificName": "Gadiformes"},
            {"Rank": "family", "ScientificName": "Gadidae"},
            {"Rank": "genus", "ScientificName": "Micromesistius"},
        ],
    }


class CleanTaxidsTests(unittest.TestCase):
    def test_splits_multi_and_dedups_preserving_order(self):
        self.assertEqual(
            ncbi_lineage._clean_taxids(["58100;9606", "9606", "", "abc", "  42 "]),
            ["58100", "9606", "42"],
        )


class LineageFromRecordTests(unittest.TestCase):
    def test_kingdom_beats_superkingdom_and_species_is_own_name(self):
        ranks = ncbi_lineage._lineage_from_record(_fish_record())
        self.assertEqual(ranks["kingdom"], "Metazoa")  # not Eukaryota
        self.assertEqual(ranks["phylum"], "Chordata")
        self.assertEqual(ranks["class_"], "Actinopteri")
        self.assertEqual(ranks["family"], "Gadidae")
        self.assertEqual(ranks["genus"], "Micromesistius")
        self.assertEqual(ranks["species"], "Micromesistius poutassou")


class FillRecordsTests(unittest.TestCase):
    def test_only_blank_fields_are_filled(self):
        lineages = {"58100": {"kingdom": "Metazoa", "family": "Gadidae", "species": "X"}}
        records = [
            {"taxid": "58100", "kingdom": "", "family": "", "species": "Micromesistius poutassou"},
            {"taxid": "999", "kingdom": ""},  # unknown taxid: untouched
        ]
        ncbi_lineage.fill_records(records, lineages)
        self.assertEqual(records[0]["kingdom"], "Metazoa")
        self.assertEqual(records[0]["family"], "Gadidae")
        # species was already asserted by BLAST -> never overwritten
        self.assertEqual(records[0]["species"], "Micromesistius poutassou")
        self.assertEqual(records[1]["kingdom"], "")


class ResolveLineagesTests(unittest.TestCase):
    def test_no_email_warns_and_resolves_nothing(self):
        with patch.object(ncbi_lineage, "get_entrez_config", lambda: {"email": "", "api_key": ""}):
            lineages, warnings = ncbi_lineage.resolve_lineages(["58100"])
        self.assertEqual(lineages, {})
        self.assertTrue(warnings and "ENTREZ_EMAIL" in warnings[0])

    def test_cache_short_circuits_without_any_network(self):
        with TemporaryDirectory() as d:
            cache = Path(d) / "lineage_cache.json"
            cache.write_text('{"58100": {"kingdom": "Metazoa"}}', encoding="utf-8")
            # No email set: if the cache is honoured, we still get the hit and no warning.
            with patch.object(ncbi_lineage, "get_entrez_config", lambda: {"email": "", "api_key": ""}):
                lineages, warnings = ncbi_lineage.resolve_lineages(["58100"], cache_path=cache)
        self.assertEqual(lineages["58100"]["kingdom"], "Metazoa")
        self.assertEqual(warnings, [])

    def test_fetch_populates_cache_and_maps_akataxids(self):
        record = _fish_record()
        record["AkaTaxIds"] = ["58101"]  # an old id merged into 58100
        fake_entrez = types.SimpleNamespace(
            email=None,
            api_key=None,
            efetch=lambda **kw: types.SimpleNamespace(close=lambda: None),
            read=lambda handle: [record],
        )
        fake_bio = types.ModuleType("Bio")
        fake_bio.Entrez = fake_entrez
        with TemporaryDirectory() as d:
            cache = Path(d) / "sub" / "lineage_cache.json"
            with patch.object(ncbi_lineage, "get_entrez_config", lambda: {"email": "me@x.org", "api_key": "k"}), \
                 patch.dict(sys.modules, {"Bio": fake_bio, "Bio.Entrez": fake_entrez}):
                lineages, warnings = ncbi_lineage.resolve_lineages(["58101"], cache_path=cache)
            self.assertEqual(warnings, [])
            self.assertEqual(lineages["58101"]["family"], "Gadidae")
            self.assertTrue(cache.exists())  # cache written for reuse
            self.assertIn("58100", ncbi_lineage._load_cache(cache))  # primary id cached too

    def test_network_error_is_swallowed_with_a_warning(self):
        def _boom(**kw):
            raise OSError("no route to host")

        fake_entrez = types.SimpleNamespace(email=None, api_key=None, efetch=_boom, read=lambda h: [])
        fake_bio = types.ModuleType("Bio")
        fake_bio.Entrez = fake_entrez
        with patch.object(ncbi_lineage, "get_entrez_config", lambda: {"email": "me@x.org", "api_key": ""}), \
             patch.dict(sys.modules, {"Bio": fake_bio, "Bio.Entrez": fake_entrez}):
            lineages, warnings = ncbi_lineage.resolve_lineages(["58100"])
        self.assertEqual(lineages, {})
        self.assertTrue(warnings and "offline" in warnings[0].lower())


if __name__ == "__main__":
    unittest.main()
