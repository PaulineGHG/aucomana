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
        self.assertIsInstance(MB.data_metabolites_consumed, pd.DataFrame)
        columns_list = [x + MB.STR_CONSUME for x in MB.species_list]
        self.assertListEqual(MB.metabolites_list, list(MB.data_metabolites_consumed.index))
        self.assertListEqual(columns_list, list(MB.data_metabolites_consumed.columns))
        self.assertEqual(MB.data_metabolites_consumed.loc['GLYCOLLATE', 'HS_rxn_consume'],
                         'GLYCOLATEDEHYDRO-RXN;RXN-969;RXN0-7229')
        self.assertEqual(MB.data_metabolites_consumed.loc['DNA-Ligase-L-lysine-adenylate', 'LF82_rxn_consume'],
                         'RXN-17918')

    def test_data_metabolites_produced(self):
        self.assertIsInstance(MB.data_metabolites_produced, pd.DataFrame)
        columns_list = [x + MB.STR_PRODUCE for x in MB.species_list]
        self.assertListEqual(MB.metabolites_list, list(MB.data_metabolites_produced.index))
        self.assertListEqual(columns_list, list(MB.data_metabolites_produced.columns))
        self.assertEqual(MB.data_metabolites_produced.loc['Alkylated-Bases', 'HS_rxn_produce'], '3.2.2.21-RXN')
        self.assertEqual(MB.data_metabolites_produced.loc['RNA-3prime-Cytidine-3prime-P', 'CFT073_rxn_produce'],
                         'RXN-19935;RXN-19932')

    def test_data_metabolites(self):
        self.assertIsInstance(MB.data_metabolites, pd.DataFrame)
        self.assertListEqual(MB.metabolites_list, list(MB.data_metabolites.index))
        self.assertListEqual(MB.species_list, list(MB.data_metabolites.columns))
        self.assertEqual(MB.data_metabolites.loc['GDP-4-DEHYDRO-6-L-DEOXYGALACTOSE', 'HS'], 0)
        self.assertEqual(MB.data_metabolites.loc['GLN', 'LF82'], 1)
