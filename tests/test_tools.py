"""seqdate verdicts, relative-annotation transfer and UniSave date handling (offline)."""

import datetime as dt
import unittest

from cutoffaudit.homolog import as_entry
from cutoffaudit.match import supports
from cutoffaudit.seqdate import verdict
from cutoffaudit.uniprot import release_date


class VerdictTest(unittest.TestCase):
    def test_labels(self):
        cut = "2025-01-01"
        self.assertEqual(verdict({"first_seen": "2002-10-01", "cluster_members_before_cutoff": 20}, cut), "old_sequence")
        self.assertEqual(verdict({"first_seen": "2025-02-05", "cluster_members_before_cutoff": 3}, cut), "new_sequence_old_family")
        self.assertEqual(verdict({"first_seen": "2025-02-05", "cluster_members_before_cutoff": 0}, cut), "new")
        self.assertEqual(verdict({"first_seen": "", "cluster_members_before_cutoff": 0}, cut), "unknown")
        self.assertEqual(verdict({"first_seen": "2025-02-05", "cluster_members_before_cutoff": -1}, cut), "unknown")


class TransferTest(unittest.TestCase):
    def test_relative_annotations_support_options(self):
        e = as_entry({"go": [["GO:0004867", "F", "serine-type endopeptidase inhibitor activity"]], "ec": ["3.4.21.4"],
                      "cc": {"COFACTOR": ["Name=Zn(2+); Xref=ChEBI:CHEBI:29105;"],
                             "PATHWAY": ["Amino-acid degradation; L-lysine degradation."]}, "kw": ["protease inhibitor"]})
        self.assertTrue(supports("GO_molecular_function", "serine-type endopeptidase inhibitor activity", e))
        self.assertTrue(supports("EC_number", "3.4.21.4", e))
        self.assertTrue(supports("cofactor", "Zn(2+)", e))
        self.assertTrue(supports("pathway", "Amino-acid degradation", e))
        self.assertFalse(supports("cofactor", "FAD", e))


class UniSaveDateTest(unittest.TestCase):
    def test_release_dates_compare_as_dates_not_strings(self):
        # string comparison would put "02-Sep-2026" before "2025-01-01"
        self.assertGreater(release_date({"firstReleaseDate": "02-Sep-2026"}), dt.date(2025, 1, 1))
        self.assertLess(release_date({"firstReleaseDate": "27-Nov-2024"}), dt.date(2025, 1, 1))


if __name__ == "__main__":
    unittest.main()
