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

    def test_tip_seqs_normal(self):
        result = preprocess_tip_sequences(self.tips_with_normal_format)
        self.assertEqual(list(result), [BiologicalSequence("ATC",
                                                           id="SSS456_1")])

    def test_tip_seqs_random_characters(self):
        result = preprocess_tip_sequences(self.tips_random_characters)
        self.assertEqual(list(result),
                         [BiologicalSequence("AGGAAAAA", id="JJJ123_1")])

tips_with_normal_format = """>SSS456
ATC
"""
tips_random_characters = """>JJJ123_1
AgGS S@P
"""

if __name__ == "__main__":
    unittest.main()
