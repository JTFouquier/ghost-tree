# ----------------------------------------------------------------------------
# Copyright (c) 2015--, ghost-tree development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
from io import StringIO

import numpy as np

from ghosttree.util import compare_tip_to_tip_distances


# name of class corresponds to the function we are testing
# we test library code NOT command line so we use the library code function
# names
class TestCompareTipToTipDistances(unittest.TestCase):
    def setUp(self):
        self.tree1 = StringIO(tree1)
        self.tree1_copy = StringIO(tree1)
        self.tree2 = StringIO(tree2)

    def test_same_trees(self):
        np.random.seed(1234)  # remove randomness for sake of testing p-value
        result = compare_tip_to_tip_distances(self.tree1, self.tree1_copy)
        self.assertEqual(result, (1.0, 0.042, 4))

    def test_different_trees(self):
        np.random.seed(1234)  # remove randomness for sake of testing p-value
        coeff, p_value, n = compare_tip_to_tip_distances(self.tree1,
                                                         self.tree2)
        self.assertAlmostEqual(coeff, 0.59603956067926978)
        self.assertAlmostEqual(p_value, 0.69599999999999995)
        self.assertEqual(n, 3)

tree1 = "[example](a:0.1, 'b_b''':0.2, (c:0.3, d_d:0.4)e:0.5)f:0.0;"
tree2 = "[example](a:0.1, 'b_b''':0.2, (g:0.3, d_d:0.1)e:0.1)f:0.0;"

# lets you run tests in this file individually, if you want to test one
# piece at a time
if __name__ == "__main__":
    unittest.main()
