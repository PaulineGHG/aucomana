import os.path
import unittest
import pandas as pd
from analysis_runs.analysis import Analysis
from analysis_runs.pathways import Pathways

STUDY_PATH = 'Study_folder'
RUNS_PATH = 'Runs_aucome'

RUN = 'bact7'
A = Analysis(RUNS_PATH, STUDY_PATH)
PW = A.pathways(RUN)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(PW, Pathways)

    def test_init_species_default(self):
        species_obj = ['HS', 'IAI1', 'CFT073', 'UTI89', 'sf301', 'ec042', 'LF82']
        self.assertEqual(species_obj, PW.species_list)

    def test_species_precised(self):
        # SPECIES_LIST = ['HS', 'UTI89', 'ec042']
        sp_l = ['HS', 'UTI89', 'ec042']
        pws = A.pathways(RUN, species_list=sp_l)
        self.assertEqual(sp_l, pws.species_list)

    def test_name(self):
        self.assertEqual(PW.name, 'bact7')

    def test_pathways(self):
        self.assertIs(type(PW.pathways_list), list)
        self.assertEqual(1281, len(PW.pathways_list))
        self.assertEqual('PWY-8010', PW.pathways_list[0])
        self.assertEqual('PWY-5517', PW.pathways_list[-1])

    def test_pathways_filtered1(self):
        # OUT = 1
        pwo1 = A.pathways(RUN, out=1)
        self.assertEqual(1225, len(pwo1.pathways_list))
        self.assertIn('PWY0-1533', pwo1.pathways_list)
        self.assertIn('PWY-6416', pwo1.pathways_list)
        self.assertNotIn('PWY0-1569', pwo1.pathways_list)

        # OUT = 2
        pwo2 = A.pathways(RUN, out=2)
        self.assertIn('PWY-7798', pwo2.pathways_list)
        self.assertNotIn('PWY0-1569', pwo2.pathways_list)

    def test_data_pathways_str(self):
        self.assertIsInstance(R.data_reactions, pd.DataFrame)
        self.assertListEqual(R.reactions_list, list(R.data_reactions.index))
        self.assertListEqual(R.species_list, list(R.data_reactions.columns))
        self.assertEqual(R.data_reactions.loc['RXN0-962', 'CFT073'], 0)

    def test_data_pathways_float(self):
        self.assertIsInstance(R.data_reactions, pd.DataFrame)
        self.assertListEqual(R.reactions_list, list(R.data_reactions.index))
        self.assertListEqual(R.species_list, list(R.data_reactions.columns))
        self.assertEqual(R.data_reactions.loc['RXN0-962', 'CFT073'], 0)
    #
    # def test_data_genes_assoc(self):
    #     self.assertIsInstance(R.data_genes_assoc, pd.DataFrame)
    #     self.assertListEqual(R.reactions_list, list(R.data_genes_assoc.index))
    #     index_list = [x + R.STR_GENE_ASSOC for x in R.species_list]
    #     self.assertListEqual(index_list, list(R.data_genes_assoc.columns))