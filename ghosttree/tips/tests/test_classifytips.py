import unittest
from StringIO import StringIO

from ghosttree.tips.classifytips import find_rep_otu_genus


class TestFindRepOtuGenus(unittest.TestCase):
    def setUp(self):
        self.otus = StringIO(otus)
        self.taxonomy = StringIO(taxonomy)
        self.modfile1 = StringIO(modfile1)
        self.modfile2 = StringIO(modfile2)

    def test_tip_seqs_and_taxonomy_correct(self):
        testotu = open("testotu.txt", "w")
        testotu.write(otus)
        testotu.close()
        testotu = open("testotu.txt", "U")
        test = open("test.txt", "w")
        test.write(taxonomy)
        test.close()
        test = open("test.txt", "U")
        result = find_rep_otu_genus(testotu, test, modfile1, modfile2)
        testotu.close()
        test.close()
        modfile1.close()
        modfile2.close()
        self.assertDictEqual(result, {'Candida': ['A222', 'A223'],
                                      'Murcor': ['A333', 'A334'],
                                      'Phoma': ['A111', 'A112']})

otus = """A111\tA111\tA112
A222\tA222\tA223
A333\tA333\tA334
"""

taxonomy = """A111\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Phoma;s__El
A112\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Phoma;s__El
A222\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
A223\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Candida;s__El
A333\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Murcor;s__El
A334\tk__Fungi;p__Asco;c__Do;o__My;f__Els;g__Murcor;s__El
"""

modfile1 = open("out1.txt", "w")
modfile2 = open("out2.txt", "w")


if __name__ == "__main__":
    unittest.main()
