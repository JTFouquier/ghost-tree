import unittest
from StringIO import StringIO

from skbio import BiologicalSequence
from skbio import SequenceCollection

from ghosttree.scaffold.hybridtree import _make_nr_backbone_alignment
from ghosttree.scaffold.hybridtree import _create_taxonomy_dic
from ghosttree.scaffold.hybridtree import _make_mini_otu_files
from ghosttree.scaffold.hybridtree import _tips_genus_accession_dic


class TestScaffoldTipsIntoBackbone(unittest.TestCase):
    def setUp(self):
        self.otu_clusters = StringIO(otu_clusters)
        self.tips_taxonomy = StringIO(tips_taxonomy)
        self.tips_taxonomy_unids = StringIO(tips_taxonomy_unids)
        self.tips_taxonomy_none = StringIO(tips_taxonomy_none)
        self.tips_sequences = StringIO(tips_sequences)
        self.backbone_alignment = StringIO(backbone_alignment)
        self.ghost_tree_fp = StringIO(ghost_tree_fp)
        self.tips_genus_dic_few = tips_genus_dic_few
        self.tips_genus_dic_none = tips_genus_dic_none
        self.key_node = key_node

    def test_make_nr_backbone_alignment_few(self):
        result = _make_nr_backbone_alignment(self.backbone_alignment,
                                             self.tips_genus_dic_few)
        self.assertEqual(list(result), [
            BiologicalSequence("AAA---", id="PBB1", description="Phoma"),
            BiologicalSequence("AAG---", id="CBB1", description="Candida"),
        ])

    def test_make_nr_backbone_alignment_none(self):
        result = _make_nr_backbone_alignment(self.backbone_alignment,
                                             self.tips_genus_dic_none)
        self.assertEqual(list(result), [])

    def test_create_taxonomy_dic_many(self):
        test = {'P1': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Phoma;s__El',
                'P2': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__phoma;s__El',
                'C1': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El',
                'C2': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El',
                'C3': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El',
                'C4': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El',
                'C5': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El',
                'M1': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El',
                'M2': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El',
                'M3': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El',
                'M4': 'k__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El'}
        result = _create_taxonomy_dic(self.tips_taxonomy)
        self.assertDictEqual(result, test)

    def test_create_taxonomy_dic_none(self):
        with self.assertRaises(ValueError):
            list(_create_taxonomy_dic(self.tips_taxonomy_none))

    def test_make_mini_otu_files(self):
        self.tips_sequences = SequenceCollection.read(self.tips_sequences)
        result = _make_mini_otu_files(self.key_node, self.tips_genus_dic_few,
                                      self.tips_sequences)
        self.assertEqual(result, """>P1\nTTAAAAAA\n""")
        # only one loop, therefor only testing one returned sequence from list
        # with two sequence accession numbers

    def test_tips_genus_accession_dic(self):
        test = {'Candida': ['C1', 'C2', 'C3', 'M1', 'C4', 'C5', 'M4'],
                'Mucor': ['M2', 'M3'],
                'Phoma': ['P1', 'P2']}
        result = _tips_genus_accession_dic(self.otu_clusters,
                                           self.tips_taxonomy)
        self.assertDictEqual(result, test)

    def test_tips_genus_accession_dic_unidentifieds(self):
        test = {'Candida': ['C1', 'C2', 'C3', 'M1'],
                'Mucor': ['M2', 'M3', 'C4', 'C5', 'M4'],
                'Phoma': ['P1', 'P2']}
        result = _tips_genus_accession_dic(self.otu_clusters,
                                           self.tips_taxonomy_unids)
        self.assertDictEqual(result, test)

key_node = 'Phoma'

tips_genus_dic_few = {'Candida': ['C1', 'C2', 'C3', 'P2'], 'Phoma': ['P1']}
tips_genus_dic_none = {'Asperggg': ['C2', 'C3', 'P2'], 'Nope': ['P1', 'M1']}

otu_clusters = """C1\tC1\tC2\tC3\tM1
M2\tM2\tM3
P1\tP1\tP2
C4\tC4\tC5\tM4
"""

tips_taxonomy = """P1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Phoma;s__El
P2\tk__Fungi;p__As;c__Do;o__My;f__Els;g__phoma;s__El
C1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
C2\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
C3\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
C4\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
C5\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
M1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
M2\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
M3\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
M4\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
"""

tips_taxonomy_unids = """P1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Phoma;s__El
P2\tk__Fungi;p__As;c__Do;o__My;f__Els;g__phoma;s__El
C1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
C2\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Candida;s__El
C3\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Unidentified;s__El
C4\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Unidentified;s__El
C5\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Unidentified;s__El
M1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
M2\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
M3\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
M4\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Mucor;s__El
"""

tips_taxonomy_none = """P1\tk__Fungi;p__As;c__Do;o__My;f__Els;Phoma;s__El
P2\tk__Fungi;p__As;c__Do;o__My;f__Els;_phoma;s__El
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

ghost_tree_fp = "ghost_tree_from_unit_test.nwk"

if __name__ == "__main__":
    unittest.main()
