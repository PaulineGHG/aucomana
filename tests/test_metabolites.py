import os.path
import unittest
import pandas as pd
from analysis_runs.analysis import Analysis
from analysis_runs.metabolites import Metabolites

STUDY_PATH = 'Study_folder'
RUNS_PATH = 'Runs_aucome'

RUN = 'bact7'
A = Analysis(RUNS_PATH, STUDY_PATH)
MB = A.metabolites(RUN)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(MB, Metabolites)

    def test_init_species_default(self):
        species_obj = ['HS', 'IAI1', 'CFT073', 'UTI89', 'sf301', 'ec042', 'LF82']
        self.assertEqual(species_obj, MB.species_list)

    def test_species_precised(self):
        # SPECIES_LIST = ['HS', 'UTI89', 'ec042']
        sp_l = ['HS', 'UTI89', 'ec042']
        mbs = A.metabolites(RUN, species_list=sp_l)
        self.assertEqual(sp_l, mbs.species_list)

    def test_name(self):
        self.assertEqual(MB.name, 'bact7')

    def test_metabolites(self):
        self.assertIs(type(MB.metabolites_list), list)
        self.assertEqual(2581, len(MB.metabolites_list))
        self.assertEqual('DETHIOBIOTIN', MB.metabolites_list[0])
        self.assertEqual('MESACONATE', MB.metabolites_list[-1])

    def test_data_metabolites_consumed(self):
        self.assertIsInstance(PW.data_pathways_str, pd.DataFrame)
        columns_list = [x + PW.STR_COMP for x in PW.species_list]
        self.assertListEqual(PW.pathways_list, list(PW.data_pathways_str.index))
        self.assertListEqual(columns_list, list(PW.data_pathways_str.columns))
        self.assertEqual(PW.data_pathways_str.loc['PWY-5515', 'UTI89_completion_rate'], '1/3')
        self.assertEqual(PW.data_pathways_str.loc['PWY-6550', 'UTI89_completion_rate'], '0/5')

    def test_data_metabolites_produced(self):
        self.assertIsInstance(PW.data_pathways_float, pd.DataFrame)
        columns_list = [x + PW.STR_COMP for x in PW.species_list]
        self.assertListEqual(PW.pathways_list, list(PW.data_pathways_float.index))
        self.assertListEqual(columns_list, list(PW.data_pathways_float.columns))
        self.assertEqual(PW.data_pathways_float.loc['PWY-5515', 'UTI89_completion_rate'], float(1/3))
        self.assertEqual(PW.data_pathways_float.loc['PWY-6550', 'UTI89_completion_rate'], 0)

    def test_data_metabolites(self):
        pass
