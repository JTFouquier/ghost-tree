import unittest
from StringIO import StringIO

from skbio import BiologicalSequence

from ghosttree.tips.cluster import preprocess_tip_sequences


# Fix input OTU level sequence files (unnecessary characters and other
# requirements set forth by SWARM software)
class TestClusterTipSequences(unittest.TestCase):
    def setUp(self):
        self.tips_with_returns = StringIO(tips_with_returns)

    def test_tip_sequences_with_returns(self):
        result = preprocess_tip_sequences(self.tips_with_returns)
        self.assertEqual(list(result), [BiologicalSequence("ATC",
                                                           id="SSS456_1")])

tips_with_returns = """>SSS456
ATC
"""

if __name__ == "__main__":
    unittest.main()
