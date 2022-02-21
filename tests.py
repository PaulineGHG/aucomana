import unittest
from reactions_loss import Reactions

FILE_TEST = "data/run2_reactions.tsv"


class Test(unittest.TestCase):

    def test_init(self):
        R = Reactions(FILE_TEST)
        self.assertEqual(type(R), Reactions)


if __name__ == '__main__':
    unittest.main()


