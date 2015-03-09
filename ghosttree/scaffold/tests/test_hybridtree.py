import unittest
from StringIO import StringIO

from ghosttree.scaffold.hybridtree import scaffold_tips_into_backbone


class TestScaffoldTipsIntoBackbone(unittest.TestCase):
    def setUp(self):
        self.otu_table = StringIO(otu_table)
        self.tips_taxonomy = StringIO(tips_taxonomy)
        self.tips_sequences = StringIO(tips_sequences)
        self.backbone_alignment = StringIO(backbone_alignment)
        self.ghost_tree_fp = StringIO(ghost_tree_fp)

    def test_scaffold_tips_into_backbone(self):
        result = scaffold_tips_into_backbone(self.otu_table,
                                             self.tips_taxonomy,
                                             self.tips_sequences,
                                             self.backbone_alignment,
                                             self.ghost_tree_fp)
        self.assertEqual(result, '(((P1:0.00027,P2:0.00027))PBB1:0.2289,' +
                                 '((M2:0.26368,M3:0.26368))MBB1:0.57748,' +
                                 '(((C4:0.0,C5:0.0):0.30767,(C3:0.13757,' +
                                 'M4:0.00055)0.817:0.13744,(C1:0.00055,' +
                                 '(C2:0.0,M1:0.0):0.13757)0.000:0.00055))' +
                                 'CBB1:0.23576);')


otu_table = """C1\tC1\tC2\tC3\tM1
M2\tM2\tM3
P1\tP1\tP2
C4\tC4\tC5\tM4
"""

tips_taxonomy = """P1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Phoma;s__El
P2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__phoma;s__El
C1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C4\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C5\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
M1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Mucor;s__El
M2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Mucor;s__El
M3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Mucor;s__El
M4\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Mucor;s__El
"""

tips_sequences = """>C1 RR_location
AAAAAAAA
>C2 RR_location
AAGAAAAA
>C3 RR_location
TTAAAAAA
>C4 RR_location
AAAATTAA
>C5 RR_location
AAAATTAA
>M1 RR_location
AAGAAAAA
>M2 RR_location
AAGAAAAA
>M3 RR_location
ATATAAAA
>M4 RR_location
TAAAAAAA
>P1 RR_location
TTAAAAAA
>P2 RR_location
AATTAAAA
"""

backbone_alignment = """>PBB1 Fungi;Phoma;
AAA---
>PBB2 Fungi;Phoma;
AAA-AA
>PBB3 Fungi;Phoma;
AAA-AT
>CBB1 Fungi;Cladosporium;
AAAAAA
>CBB2 Fungi;Aspergillus;
A-A---
>MBB1 Fungi;Mucor;
A-TAAA
>CBB1 Fungi;Candida;
AAG---
"""

ghost_tree_fp = "ghost_tree_from_UNIT_TEST.nwk"

if __name__ == "__main__":
    unittest.main()
