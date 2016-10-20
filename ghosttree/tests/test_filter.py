# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from io import StringIO

from skbio import Sequence

from ghosttree.filter import filter_positions


class TestFilterPositions(unittest.TestCase):
    def setUp(self):
        self.alignment_with_gaps = StringIO(alignment_with_gaps)
        self.maximum_position_entropy_10 = maximum_position_entropy_10
        self.maximum_position_entropy_50 = maximum_position_entropy_50
        self.maximum_position_entropy_100 = maximum_position_entropy_100
        self.maximum_gap_frequency_0 = maximum_gap_frequency_0
        self.maximum_gap_frequency_50 = maximum_gap_frequency_50
        self.maximum_gap_frequency_100 = maximum_gap_frequency_100

    def test_filter_gap_med_entropy_high(self):
        result = filter_positions(self.alignment_with_gaps,
                                  self.maximum_gap_frequency_50,
                                  self.maximum_position_entropy_100)

        aln = [Sequence("ACC--G-GGTA..", metadata={'id': "seq1"}),
               Sequence("TCC--G-GGCA..", metadata={'id': "seq2"})]
        msa = TabularMSA(aln, minter='id')
        print(msa)




        self.assertEqual(result, aln)

    def test_filter_gap_low_entropy_med(self):
        result = filter_positions(self.alignment_with_gaps,
                                  self.maximum_gap_frequency_0,
                                  self.maximum_position_entropy_50)
        aln = Alignment([Sequence('A', id="seq1"),
                         Sequence('A', id="seq2"),
                         Sequence('A', id="seq3"),
                         Sequence('A', id="seq4")])
        self.assertEqual(result, aln)

    def test_filter_gap_high_entropy_low(self):
        result = filter_positions(self.alignment_with_gaps,
                                  self.maximum_gap_frequency_100,
                                  self.maximum_position_entropy_10)
        aln = Alignment([Sequence('A-', id="seq1"),
                         Sequence('A-', id="seq2"),
                         Sequence('A-', id="seq3"),
                         Sequence('A-', id="seq4")])
        self.assertEqual(result, aln)


alignment_with_gaps = """>seq1
ACC-T
>seq2
ACT--
>seq3
ATG--
>seq4
ATA--
"""


maximum_gap_frequency_0 = 0.0
maximum_gap_frequency_50 = 0.5
maximum_gap_frequency_100 = 1.0

maximum_position_entropy_10 = 0.1
maximum_position_entropy_50 = 0.5
maximum_position_entropy_100 = 1.0


if __name__ == "__main__":
    unittest.main()
