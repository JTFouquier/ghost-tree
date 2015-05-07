# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from StringIO import StringIO

from skbio import BiologicalSequence

from ghosttree.silva.filter import fungi_from_fasta


class TestFungiFromFasta(unittest.TestCase):
    def setUp(self):
        self.fasta_no_fungi = StringIO(fasta_no_fungi)
        self.fasta_with_fungi = StringIO(fasta_with_fungi)
        self.fasta_many_fungi = StringIO(fasta_many_fungi)
        self.accession = StringIO(accession)
        self.accession_with_duplicates = StringIO(accession_with_duplicates)
        self.taxonomy_no_fungi = StringIO(taxonomy_no_fungi)
        self.taxonomy_with_fungi = StringIO(taxonomy_with_fungi)
        self.taxonomy_with_duplicates = StringIO(taxonomy_with_duplicates)

    def test_fasta_no_fungi(self):
        result = fungi_from_fasta(self.fasta_no_fungi,
                                  self.accession,
                                  self.taxonomy_no_fungi)
        self.assertEqual(list(result), [])

    def test_fasta_with_fungi(self):
        result = fungi_from_fasta(self.fasta_with_fungi,
                                  self.accession,
                                  self.taxonomy_with_fungi)
        self.assertEqual(list(result),
                         [BiologicalSequence("ATCG", id="AB21",
                                             description="Fungi")])

    def test_fasta_with_many_fungi(self):
        result = fungi_from_fasta(self.fasta_many_fungi,
                                  self.accession,
                                  self.taxonomy_with_fungi)
        self.assertEqual(list(result), [
            BiologicalSequence("GGGG", id="AB123", description="Fungi"),
            BiologicalSequence("CCCC", id="AB125", description="Fungi"),
            BiologicalSequence("AAAA", id="AB126", description="Fungi"),
        ])

    def test_duplicate_accession_numbers(self):
        with self.assertRaises(ValueError):
            list(fungi_from_fasta(self.fasta_with_fungi,
                                  self.accession_with_duplicates,
                                  self.taxonomy_with_fungi))

    def test_duplicate_map_numbers(self):
        with self.assertRaises(ValueError):
            list(fungi_from_fasta(self.fasta_with_fungi,
                                  self.accession,
                                  self.taxonomy_with_duplicates))


fasta_many_fungi = """>AB123 Fungi
GGGG
>AB124 Bacteria
TTTT
>AB125 Fungi
CCCC
>AB126 Fungi
AAAA
>AB100 Fungi
ATAC
"""

fasta_no_fungi = """>AB21 Bacteria
ATCG
"""
fasta_with_fungi = """>AB21 Fungi
ATCG
"""
taxonomy_with_fungi = """Fungi\t123\tgenus\t\t
Eukaryote;Fungi;Fungal Species\t456\tgenus\t\t
Fungi;Fungal species\t789\tgenus\t\t
Bacteria;bacterial species\t111\tgenus\t\t
Fungi;cladosporium\t100\tdomain\t\t
"""
taxonomy_no_fungi = "Bacteria\t123\tdomain\t\t"
taxonomy_with_duplicates = """Fungi\t123\tgenus\t\t
Fungi\t123\tgenus\t\t
Fungi\t122\tgenus\t\t
"""
accession = """AB21\t123
AB123\t456
AB124\t111
AB125\t789
AB126\t123
AB100\t100
"""
accession_with_duplicates = """AB21\t123
AC33\t456
AC33\t789
"""


if __name__ == "__main__":
    unittest.main()
