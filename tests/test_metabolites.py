import os.path
import unittest
import pandas as pd
from analysis_runs.analysis import Analysis
from analysis_runs.metabolites import Metabolites

STUDY_PATH = 'Study_folder'
RUNS_PATH = 'Runs_aucome'

RUN = 'bact7'
A = Analysis(RUNS_PATH, STUDY_PATH)
M = A.metabolites(RUN)


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
        self.assertIsInstance(PW.data_pathways_str, pd.DataFrame)
        columns_list = [x + PW.STR_COMP for x in PW.species_list]
        self.assertListEqual(PW.pathways_list, list(PW.data_pathways_str.index))
        self.assertListEqual(columns_list, list(PW.data_pathways_str.columns))
        self.assertEqual(PW.data_pathways_str.loc['PWY-5515', 'UTI89_completion_rate'], '1/3')
        self.assertEqual(PW.data_pathways_str.loc['PWY-6550', 'UTI89_completion_rate'], '0/5')

    def test_data_pathways_float(self):
        self.assertIsInstance(PW.data_pathways_float, pd.DataFrame)
        columns_list = [x + PW.STR_COMP for x in PW.species_list]
        self.assertListEqual(PW.pathways_list, list(PW.data_pathways_float.index))
        self.assertListEqual(columns_list, list(PW.data_pathways_float.columns))
        self.assertEqual(PW.data_pathways_float.loc['PWY-5515', 'UTI89_completion_rate'], float(1/3))
        self.assertEqual(PW.data_pathways_float.loc['PWY-6550', 'UTI89_completion_rate'], 0)