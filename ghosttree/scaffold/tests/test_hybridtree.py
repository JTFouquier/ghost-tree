# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import os
from io import StringIO

import skbio
from skbio import Sequence

from ghosttree.scaffold.hybridtree import _make_nr_foundation_alignment
from ghosttree.scaffold.hybridtree import _create_taxonomy_dict
from ghosttree.scaffold.hybridtree import _make_mini_otu_files
from ghosttree.scaffold.hybridtree import _extension_genus_accession_dict


class TestScaffoldExtensionsIntofoundation(unittest.TestCase):
    def setUp(self):
        self.otu_clusters = StringIO(otu_clusters)
        self.extension_taxonomy = StringIO(extension_taxonomy)
        self.extension_taxonomy_unids = StringIO(extension_taxonomy_unids)
        self.extension_taxonomy_none = StringIO(extension_taxonomy_none)
        self.extension_seqs = StringIO(extension_seqs)
        self.foundation_alignment = StringIO(foundation_alignment)
        self.ghost_tree_fp = StringIO(ghost_tree_fp)
        self.extension_genus_dic_few = extension_genus_dic_few
        self.extension_genus_dic_none = extension_genus_dic_none
        self.key_node = key_node
        self.graft_letter_g = 'g'
        self.graft_level_6 = 6
        self.graft_level_5 = 5
        self.foundation_taxonomy = StringIO(foundation_taxonomy)
        self.foundation_newick = StringIO(foundation_newick)

    def test_make_nr_foundation_alignment_few(self):
        result = _make_nr_foundation_alignment(self.foundation_alignment,
                                               self.extension_genus_dic_few,
                                               self.graft_letter_g)

        self.assertEqual(list(result), [
            Sequence("AAA---", metadata={"id": "PBB1",
                                         "description": "Phoma"}),
            Sequence("AAG---", metadata={"id": "CBB3",
                                         "description": "Candida"}),
        ])

    def test_make_nr_foundation_alignment_none(self):
        result = _make_nr_foundation_alignment(self.foundation_alignment,
                                               self.extension_genus_dic_none,
                                               self.graft_letter_g)
        self.assertEqual(list(result), [])

    # def test_newick_file_few_extensions(self):
    #     result = _make_nr_foundation_newick(self.foundation_newick,
    #                                         self.extension_genus_dic_few)
    #
    #     print('result', result)

    def test_create_taxonomy_dic_many_genus(self):

        test = {'P1': 'Phoma', 'P2': 'Phoma', 'C1': 'Candida', 'C2': 'Candida',
                'C3': 'Candida', 'C4': 'Candida', 'C5': 'Candida',
                'M1': 'Mucor', 'M2': 'Mucor', 'M3': 'Mucor', 'M4': 'Mucor'}

        result = _create_taxonomy_dict(self.extension_taxonomy,
                                       self.graft_level_6)
        print('\n\n\nlalala', result)
        self.assertDictEqual(result, test)

    def test_create_taxonomy_dic_none(self):
        with self.assertRaises(ValueError):
            list(_create_taxonomy_dict(self.extension_taxonomy_none,
                                       9))

    def test_make_mini_otu_files(self):
        os.system("mkdir tmp")
        extension_seqs = list(skbio.io.read(self.extension_seqs,
                                            format='fasta'))

        _make_mini_otu_files(self.key_node,
                             self.extension_genus_dic_few,
                             extension_seqs)

        result = str(list(skbio.io.read('tmp/mini_seq_gt.fasta',
                                        format='fasta'))[0])

        os.system("rm -r tmp")
        self.assertEqual(result, 'TTAAAAAA')

    def test_extension_genus_accession_dic(self):
        test = {'Candida': ['C1', 'C2', 'C3', 'M1', 'C4', 'C5', 'M4'],
                'Mucor': ['M2', 'M3'],
                'Phoma': ['P1', 'P2']}
        result = _extension_genus_accession_dict(self.otu_clusters,
                                                 self.extension_taxonomy,
                                                 self.graft_level_6)
        self.assertDictEqual(result, test)

    def test_extension_genus_accession_dic_unidentifieds_level_6_genus(self):
        test = {'Candida': ['C1', 'C2', 'C3', 'M1'],
                'Mucor': ['M2', 'M3', 'C4', 'C5', 'M4'],
                'Phoma': ['P1', 'P2']}
        result = _extension_genus_accession_dict(self.otu_clusters,
                                                 self.extension_taxonomy_unids,
                                                 self.graft_level_6)
        self.assertDictEqual(result, test)

    def test_extension_genus_accession_dic_unidentifieds_level_5_family(self):
        test = {'Saccharomycetaceae': ['C1', 'C2', 'C3', 'M1', 'C4', 'C5',
                'M4'], 'Didymellaceae': ['P1', 'P2'], 'Mucoraceae':
                ['M2', 'M3']}
        result = _extension_genus_accession_dict(self.otu_clusters,
                                                 self.extension_taxonomy_unids,
                                                 self.graft_level_5)

        self.assertDictEqual(result, test)


key_node = 'Phoma'

extension_genus_dic_few = {'Candida': ['C1', 'C2', 'C3', 'P2'],
                           'Phoma': ['P1']}
extension_genus_dic_none = {'Asperggg': ['C2', 'C3', 'P2'],
                            'Nope': ['P1', 'M1']}

otu_clusters = """C1\tC1\tC2\tC3\tM1
M2\tM2\tM3
P1\tP1\tP2
C4\tC4\tC5\tM4
"""

extension_taxonomy = """P1\tk__Fungi;p__As;c__Do;o__My;f__Els;g__Phoma;s__El
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

extension_taxonomy_unids = """P1\tk__Fungi;p__As;c__Do;o__My;f__Didymellaceae;g__Phoma;s__El
P2\tk__Fungi;p__As;c__Do;o__My;f__Didymellaceae;g__Phoma;s__El
C1\tk__Fungi;p__As;c__Do;o__My;f__Saccharomycetaceae;g__Candida;s__El
C2\tk__Fungi;p__As;c__Do;o__My;f__Saccharomycetaceae;g__Candida;s__El
C3\tk__Fungi;p__As;c__Do;o__My;f__Saccharomycetaceae;g__Unidentified;s__El
C4\tk__Fungi;p__As;c__Do;o__My;f__Saccharomycetaceae;g__Unidentified;s__El
C5\tk__Fungi;p__As;c__Do;o__My;f__Unidentified;g__Unidentified;s__El
M1\tk__Fungi;p__As;c__Do;o__My;f__Mucoraceae;g__Mucor;s__El
M2\tk__Fungi;p__As;c__Do;o__My;f__Mucoraceae;g__Mucor;s__El
M3\tk__Fungi;p__As;c__Do;o__My;f__Mucoraceae;g__Mucor;s__El
M4\tk__Fungi;p__As;c__Do;o__My;f__Mucoraceae;g__Mucor;s__El
"""

extension_taxonomy_none = """P1\tk__Fungi;p__As;c__Do;o__My;f__Els;Phoma;s__El
P2\tk__Fungi;p__As;c__Do;o__My;f__Els;_phoma;s__El
"""

extension_seqs = """>C1 RR_location
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

foundation_alignment = """>PBB1 Fungi;Phoma;
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
>CBB3 Fungi;Candida;
AAG---
"""

foundation_newick = "(CBB2:0.00055,MBB1:0.23625,(CBB1:0.00055,(PBB2:0.00055," \
                    "(CBB3:0.44621,(PBB3:0.00055,PBB1:0.00055)0.736:0.08656)" \
                    "0.801:0.14667)0.437:0.00054)0.316:0.00055);"

foundation_taxonomy = """PBB1\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Phoma;s__El
PBB2\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Phoma;s__El
PBB3\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Phoma;s__El
CBB1\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Cladosporium;s__El
CBB2\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Aspergillus;s__El
MBB1\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Mucor;s__El
CBB3\tk__Fungi;p__Asc;c__Do;o__My;f__El;g__Candida;s__El"""

ghost_tree_fp = "ghost_tree_from_unit_test.nwk"

if __name__ == "__main__":
    unittest.main()
