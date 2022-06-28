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

    def test_data_rxn_assoc(self):
        self.assertIsInstance(PW.data_rxn_assoc, pd.DataFrame)
        self.assertListEqual(PW.pathways_list, list(PW.data_rxn_assoc.index))
        columns_list = [x + PW.STR_RXN_ASSOC for x in PW.species_list]
        self.assertListEqual(columns_list, list(PW.data_rxn_assoc.columns))

    def test_is_present(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_present('HS', 'PWY-8259'))
        self.assertFalse(PW.is_present('UTI89', 'PWY-7716'))

        # UNIQUE = TRUE
        self.assertTrue(PW.is_present('sf301', 'PWY-5517', unique=True))
        self.assertFalse(PW.is_present('sf301', 'AEROBACTINSYN-PWY', unique=True))

    # def test_is_absent(self):
    #     # NO PARAMETER (UNIQUE = FALSE)
    #     self.assertTrue(R.is_absent('CFT073', 'RXN0-962'))
    #     self.assertFalse(R.is_absent('CFT073', 'RXN-12615'))
    #
    #     # UNIQUE = TRUE
    #     self.assertTrue(R.is_absent('ec042', '5.1.99.4-RXN', unique=True))
    #     self.assertFalse(R.is_absent('ec042', 'RXN-12615', unique=True))
    #
    # def test_get_rxn_present(self):
    #     # NO PARAMETER
    #     dict_rxn_pres = R.get_rxn_present()
    #     self.assertIs(type(dict_rxn_pres), dict)
    #     self.assertEqual(dict_rxn_pres['HS'][0], 2304)
    #     self.assertIs(type(dict_rxn_pres['HS'][1]), set)
    #     self.assertIn('RXN-15149', dict_rxn_pres['HS'][1])
    #     self.assertEqual(set(dict_rxn_pres.keys()), set(R.species_list))
    #
    #     # UNIQUE = TRUE
    #     dict_rxn_pres_uni = R.get_rxn_present(unique=True)
    #     self.assertEqual(dict_rxn_pres_uni['HS'][0], 0)
    #     self.assertEqual(dict_rxn_pres_uni['IAI1'][0], 5)
    #     self.assertIn('RXN-20706', dict_rxn_pres_uni['IAI1'][1])
    #
    #     # SPECIES STR
    #     dict_rxn_pres = R.get_rxn_present(species='HS')
    #     self.assertEqual(dict_rxn_pres.keys(), {'HS'})
    #
    #     # SPECIES LIST
    #     dict_rxn_pres = R.get_rxn_present(species=['HS', 'IAI1', 'ec042'])
    #     self.assertEqual(dict_rxn_pres.keys(), {'HS', 'IAI1', 'ec042'})
    #
    # def test_get_rxn_absent(self):
    #     # NO PARAMETER
    #     dict_rxn_abs = R.get_rxn_absent()
    #     self.assertIs(type(dict_rxn_abs), dict)
    #     self.assertEqual(dict_rxn_abs['HS'][0], 98)
    #     self.assertIs(type(dict_rxn_abs['HS'][1]), set)
    #     self.assertIn('TRANS-RXN0-268', dict_rxn_abs['HS'][1])
    #     self.assertEqual(set(dict_rxn_abs.keys()), set(R.species_list))
    #
    #     # UNIQUE = TRUE
    #     dict_rxn_abs_uni = R.get_rxn_absent(unique=True)
    #     self.assertEqual(dict_rxn_abs_uni['HS'][0], 15)
    #     self.assertEqual(dict_rxn_abs_uni['IAI1'][0], 0)
    #     self.assertIn('RXN0-7341', dict_rxn_abs_uni['HS'][1])
    #
    #     # SPECIES STR
    #     dict_rxn_pres = R.get_rxn_present(species='HS')
    #     self.assertEqual(dict_rxn_pres.keys(), {'HS'})
    #
    #     # SPECIES LIST
    #     dict_rxn_pres = R.get_rxn_present(species=['HS', 'IAI1', 'ec042'])
    #     self.assertEqual(dict_rxn_pres.keys(), {'HS', 'IAI1', 'ec042'})