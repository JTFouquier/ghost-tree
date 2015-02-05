
import unittest
from StringIO import StringIO

from ghosttree.silva.filter import fungi_from_fasta


class TestFungiFromFasta(unittest.TestCase):
    def test_fasta_no_fungi(self):
        fasta_file = StringIO(fasta_no_fungi)
        accession_file = StringIO(accession)
        taxonomy_file = StringIO(taxonomy)
        result = fungi_from_fasta(fasta_file, accession_file, taxonomy_file)
        self.assertEqual(list(result), [])

fasta_no_fungi = """>AB21 Bacteria
ATCG
"""
accession = "AB21\t123"
taxonomy = "Bacteria\t123\tdomain"


if __name__ == "__main__":
    unittest.main()
