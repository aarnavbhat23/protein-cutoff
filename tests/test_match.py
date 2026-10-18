"""Flat-file parsing and per-task matching on small synthetic entries."""

import unittest

from cutoffaudit.flatfile import parse
from cutoffaudit.match import classify, supports

NEW = """ID   TEST_HUMAN              Reviewed;         330 AA.
AC   Q8MVC4; Q00001;
DE   RecName: Full=Thing;
DE            EC=3.1.3.16 {ECO:0000269|PubMed:1};
CC   -!- FUNCTION: Does things.
CC   -!- CATALYTIC ACTIVITY:
CC       Reaction=O-phospho-L-seryl-[protein] + H2O = L-seryl-[protein] +
CC         phosphate; Xref=Rhea:RHEA:20629; EC=3.1.3.16;
CC         Evidence={ECO:0000269|PubMed:1};
CC   -!- COFACTOR:
CC       Name=Fe(2+); Xref=ChEBI:CHEBI:29033; Evidence={ECO:0000250};
CC   -!- PATHWAY: Capsule biogenesis; capsule polysaccharide biosynthesis.
CC   -!- BIOPHYSICOCHEMICAL PROPERTIES:
CC       pH dependence:
CC         Optimum pH is 7.5. {ECO:0000269|PubMed:1};
CC       Temperature dependence:
CC         Optimum temperature is 30 degrees Celsius. {ECO:0000269|PubMed:1};
CC   ---------------------------------------------------------------------------
CC   Copyrighted by the UniProt Consortium, see https://www.uniprot.org/terms
DR   GO; GO:0004867; F:serine-type endopeptidase inhibitor activity; IEA:InterPro.
DR   GO; GO:0005576; C:extracellular region; IEA:UniProtKB-SubCell.
KW   Disulfide bond; Protease inhibitor.
RX   PubMed=15116248; DOI=10.1/x;
FT   ACT_SITE        154
FT                   /note="Acyl-thioester intermediate"
FT                   /evidence="ECO:0000255|PROSITE-ProRule:PRU0001"
FT   MOTIF           66..75
FT                   /note="Nuclear localization signal"
FT   TRANSMEM        237..255
FT                   /note="Helical"
FT   TRANSMEM        289..311
SQ   SEQUENCE   12 AA;  1000 MW;  0 CRC64;
     MKVLAAGIVG LL
//
"""

OLD_FT = """ID   OLD_HUMAN               Unreviewed;        50 AA.
AC   P99999;
FT   ACT_SITE    154    154       Acyl-thioester intermediate.
FT   MOTIF        66     75       NLS.
SQ   SEQUENCE   3 AA;
     MKV
//
"""


class ParseTest(unittest.TestCase):
    def setUp(self):
        self.e = parse(NEW)

    def test_fields(self):
        e = self.e
        self.assertEqual(e.accessions, ["Q8MVC4", "Q00001"])
        self.assertTrue(e.reviewed)
        self.assertIn("3.1.3.16", e.ec)
        self.assertIn("fe(2+)", e.cofactors())
        self.assertIn("o-phospho-l-seryl-[protein] + h2o = l-seryl-[protein] + phosphate", e.reactions())
        self.assertIn("serine-type endopeptidase inhibitor activity", e.go_names("F"))
        self.assertEqual(e.feature_locations("ACT_SITE"), ["154"])
        self.assertEqual(e.feature_locations("TRANSMEM"), ["237..255", "289..311"])
        self.assertEqual(e.pubmed, ["15116248"])
        self.assertEqual(e.sequence, "MKVLAAGIVGLL")
        self.assertEqual(len(e.optimum("pH")), 1)

    def test_old_feature_layout(self):
        e = parse(OLD_FT)
        self.assertEqual(e.feature_locations("ACT_SITE"), ["154"])
        self.assertEqual(e.feature_locations("MOTIF"), ["66..75"])
        self.assertFalse(e.reviewed)


class SupportsTest(unittest.TestCase):
    def setUp(self):
        self.e = parse(NEW)

    def test_each_task_rule(self):
        e = self.e
        cases = [
            ("GO_molecular_function", "serine-type endopeptidase inhibitor activity", True),
            ("GO_molecular_function", "extracellular region", False),      # right term, wrong aspect
            ("GO_cellular_component", "extracellular region", True),
            ("EC_number", "3.1.3.16", True), ("EC_number", "3.1.3.1", False),
            ("cofactor", "Fe(2+)", True), ("cofactor", "FAD", False),
            ("catalytic_activity", "O-phospho-L-seryl-[protein] + H2O = L-seryl-[protein] + phosphate", True),
            ("pathway", "Capsule biogenesis", True), ("pathway", "Amino-acid degradation", False),
            ("active_site", "Position 154 (Acyl-thioester intermediate)", True),
            ("active_site", "Position 15 (Proton donor)", False),
            ("motif_position", "66..75", True), ("motif_position", "66..7", False),
            ("transmembrane", "237..255; 289..311", True), ("transmembrane", "237..255", False),
            ("ph", "7.5", True), ("ph", "7", False), ("ph", "6.0-7.9", False),
            ("temperature", "30", True), ("temperature", "3", False),
        ]
        for task, option, want in cases:
            with self.subTest(task=task, option=option):
                self.assertEqual(supports(task, option, e), want)

    def test_classify(self):
        e = self.e
        ch = {"A": "FAD", "B": "Fe(2+)", "C": "heme", "D": "Mg(2+)"}
        self.assertEqual(classify("cofactor", ch, "B", e)["status"], "leaked")
        self.assertEqual(classify("cofactor", ch, "A", e)["status"], "misleading")
        self.assertEqual(classify("cofactor", ch, "B", None)["status"], "no_entry")
        r = classify("cofactor", {"A": "FAD", "B": "heme"}, "A", e)
        self.assertEqual((r["status"], r["lookup_acc"]), ("absent", 0.5))


if __name__ == "__main__":
    unittest.main()
