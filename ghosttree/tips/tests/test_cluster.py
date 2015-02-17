import unittest
from StringIO import StringIO

from skbio import BiologicalSequence

from ghosttree.tips.cluster import preprocess_tip_sequences

# Fix input OTU level sequence files (unnecessary characters and other
# requirements set forth by SWARM software)


class TestClusterTipSequences(unittest.TestCase):
    def setUp(self):
        self.tips_with_normal_format = StringIO(tips_with_normal_format)
        self.tips_random_characters = StringIO(tips_random_characters)
        self.taxonomy_file_genera = StringIO(taxonomy_file_genera)
        self.taxonomy_file_no_genera = StringIO(taxonomy_file_no_genera)

    def test_tip_seqs_and_taxonomy_correct(self):
        result = preprocess_tip_sequences(self.tips_with_normal_format,
                                          self.taxonomy_file_genera)
        self.assertEqual(list(result), [BiologicalSequence("ATC",
                                                           id="SSS456_1")])

    def test_tip_seqs_random_characters(self):
        result = preprocess_tip_sequences(self.tips_random_characters,
                                          self.taxonomy_file_genera)
        self.assertEqual(list(result),
                         [BiologicalSequence("AGGAAAAA", id="JJJ123_1")])

    def text_taxonomy_file_no_genera(self):
        with self.assertRaises(ValueError):
            list(preprocess_tip_sequences(self.tips_with_normal_format,
                                          self.taxonomy_file_no_genera))

tips_with_normal_format = """>SSS456
ATC
"""
tips_random_characters = """>JJJ123_1
AgGS S@P
"""

taxonomy_file_genera = "JJJ999\tk__Fungi;p__;c__;o__;f__Tca;g__Th;s__species"

taxonomy_file_no_genera = """JJJ999\tk__Fungi;p__Bas;c__Ag;o__Th;f__Tcae;s__fungalspecies
TTT000\tk__Fungi;p__Bas;c__Ag;o__Th;f__Tcae;s__uncultured
"""


if __name__ == "__main__":
    unittest.main()
