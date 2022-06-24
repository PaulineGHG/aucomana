import unittest
import pandas as pd
import typing
from analysis_runs.analysis import Analysis
from analysis_runs.reactions import Reactions

STUDY_PATH = 'Study_folder'
RUNS_PATH = 'Runs_aucome'

RUN = 'bact7'
A = Analysis(RUNS_PATH, STUDY_PATH)
R = A.reactions(RUN)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(R, Reactions)

    def test_init_species_default(self):
        species_obj = ['HS', 'IAI1', 'CFT073', 'UTI89', 'sf301', 'ec042', 'LF82']
        self.assertEqual(species_obj, R.species_list)

    def test_species_precised(self):
        sp_l = ['HS', 'UTI89', 'ec042']
        rs = A.reactions(RUN, species_list=sp_l)
        self.assertEqual(sp_l, rs.species_list)

    def test_name(self):
        self.assertEqual(R.name, 'bact7')

    def test_reactions(self):
        self.assertIs(type(R.reactions_list), list)
        self.assertEqual(2402, len(R.reactions_list))
        self.assertEqual('3.2.2.21-RXN', R.reactions_list[0])
        self.assertEqual('2.7.7.46-RXN', R.reactions_list[-1])

    def test_reactions_filtered1(self):
        ro1 = A.reactions(RUN, out=1)
        self.assertEqual(2196, len(ro1.reactions_list))
        self.assertIn('RXN-17824', ro1.reactions_list)
        self.assertIn('BETAGALACTOSID-RXN', ro1.reactions_list)
        self.assertNotIn('RXN-10016', ro1.reactions_list)
        ro2 = A.reactions(RUN, out=2)
        self.assertIn('RXN-10016', ro2.reactions_list)

    def test_data_reactions(self):
        self.assertIsInstance(R.data_reactions, pd.DataFrame)
        self.assertListEqual(R.reactions_list, list(R.data_reactions.index))
        self.assertListEqual(R.species_list, list(R.data_reactions.columns))
        self.assertEqual(R.data_reactions.loc['RXN0-962', 'CFT073'], 0)

    def test_data_genes_assoc(self):
        self.assertIsInstance(R.data_genes_assoc, pd.DataFrame)
        self.assertListEqual(R.reactions_list, list(R.data_genes_assoc.index))
        index_list = [x + R.STR_GENE_ASSOC for x in R.species_list]
        self.assertListEqual(index_list, list(R.data_genes_assoc.columns))


if __name__ == '__main__':
    unittest.main()


