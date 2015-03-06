import unittest
from StringIO import StringIO

from ghosttree.tips.classifytips import find_rep_otu_genus


class TestFindRepOtuGenus(unittest.TestCase):
    def setUp(self):
        self.otus = StringIO(otus)
        self.otus_no_duplicates = StringIO(otus_no_duplicates)
        self.otus_mixed_genera_clusters = StringIO(otus_mixed_genera_clusters)
        self.otus_mixed_genera_tie = StringIO(otus_mixed_genera_tie)
        self.otus_same_genus_names = StringIO(otus_same_genus_names)
        self.taxonomy = StringIO(taxonomy)
        self.taxonomy_without_genus = StringIO(taxonomy_without_genus)
        self.modfile1 = StringIO(modfile1)
        self.modfile2 = StringIO(modfile2)

    def test_tip_seqs_and_taxonomy_correct(self):
        result = find_rep_otu_genus(self.otus, self.taxonomy, self.modfile1,
                                    self.modfile2)
        self.assertDictEqual(result, {'Candida': ['C1', 'C2'],
                                      'Murcor': ['M1', 'M2'],
                                      'Phoma': ['P1', 'P2']})

    def test_otus_no_duplicate_accessions(self):
        result = find_rep_otu_genus(self.otus_no_duplicates, self.taxonomy,
                                    self.modfile1, self.modfile2)
        self.assertDictEqual(result, {'Candida': ['C1', 'C2'],
                                      'Murcor': ['M1', 'M2'],
                                      'Phoma': ['P1', 'P2']})

    def test_tip_seqs_with_mixed_genera(self):
        result = find_rep_otu_genus(self.otus_mixed_genera_clusters,
                                    self.taxonomy, self.modfile1,
                                    self.modfile2)
        self.assertDictEqual(result, {'Candida': ['C1', 'C2', 'P3'],
                                      'Murcor': ['M1', 'M2'],
                                      'Phoma': ['P1', 'P2']})

    def test_tip_seqs_with_mixed_genera_tie(self):
        result = find_rep_otu_genus(self.otus_mixed_genera_tie,
                                    self.taxonomy, self.modfile1,
                                    self.modfile2)
        # Needs to address ties that conflict with dictionary key.
        self.assertDictEqual(result, {'Candida': ['C1', 'C2', 'M3', 'M4'],
                                      'Murcor': ['M1', 'M2'],
                                      'Phoma': ['P1', 'P2']})

    def test_otus_same_genus_names(self):
        result = find_rep_otu_genus(self.otus_same_genus_names,
                                    self.taxonomy, self.modfile1,
                                    self.modfile2)
        # Needs to address ties that conflict with dictionary key.
        self.assertDictEqual(result, {'Candida': ['C1', 'C2', 'M1', 'C3',
                                                  'C4', 'C5', 'P1']})

    def test_broken_taxonomy_file(self):
        with self.assertRaises(ValueError):
            list(find_rep_otu_genus(self.otus, self.taxonomy_without_genus,
                                    self.modfile1, self.modfile2))

otus = """C1\tC1\tC2
M1\tM1\tM2
P1\tP1\tP2
"""

otus_no_duplicates = """C1\tC2
M1\tM2
P1\tP2
"""

otus_mixed_genera_clusters = """C1\tC1\tC2\tP3
M1\tM1\tM2
P1\tP1\tP2
"""

otus_mixed_genera_tie = """C1\tC1\tC2\tM3\tM4
M1\tM1\tM2
P1\tP1\tP2
"""

otus_same_genus_names = """C1\tC1\tC2\tM1
C3\tC3\tC4\tC5\tP1
"""


taxonomy = """P1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Phoma;s__El
P2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__phoma;s__El
P3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__phoma;s__El
P4\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__phoma;s__El
P5\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Phoma;s__El
C1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C4\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
C5\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
M1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Murcor;s__El
M2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Murcor;s__El
M3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Murcor;s__El
M4\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Murcor;s__El
"""

taxonomy_without_genus = """P1\tk__Fungi;p__As;c__Do;o__M;f__E;__Pho;s__El
P2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Phoma;s__El
P3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Phoma;s__El
C1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Candida;s__El
C2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Candida;s__El
C3\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Candida;s__El
M1\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Murcor;s__El
M2\tk__Fungi;p__Asco;c__Do;o__My;f__Els;__Murcor;s__El
"""

modfile1 = "repgenusOTUfile.txt"
modfile2 = "repgenusOTUfile_non-redundant.txt"


if __name__ == "__main__":
    unittest.main()
