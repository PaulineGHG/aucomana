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

    def test_nb_rnx_pw(self):
        self.assertEqual(PW.nb_rnx_pw["NADPHOS-DEPHOS-PWY-1"], 2)
        self.assertEqual(PW.nb_rnx_pw["PWY-7716"], 4)
        self.assertEqual(PW.nb_rnx_pw["PWY-5517"], 6)

    def test_nb_rnx_pw_min(self):
        PW3 = A.pathways(RUN, nb_rxn_pw_min=3)
        self.assertNotIn("NADPHOS-DEPHOS-PWY-1", PW3.pathways_list)  # 2
        self.assertIn("PWY-2722", PW3.pathways_list)  # 3
        self.assertIn("PWY-7716", PW3.pathways_list)  # 4

    # Gain

    def test_is_present(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_present('HS', 'PWY-8259'))
        self.assertFalse(PW.is_present('UTI89', 'PWY-7716'))

        # UNIQUE = TRUE
        self.assertTrue(PW.is_present('sf301', 'PWY-5517', unique=True))
        self.assertFalse(PW.is_present('sf301', 'AEROBACTINSYN-PWY', unique=True))

    def test_is_max(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_max('HS', 'PWY-6538'))
        self.assertFalse(PW.is_max('ec042', 'PWY-6538'))

        # UNIQUE = TRUE
        self.assertTrue(PW.is_max('IAI1', 'PWY-6854', unique=True))
        self.assertFalse(PW.is_max('HS', 'PWY-6538', unique=True))

    def test_is_complete(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_complete('HS', 'GLYCOLYSIS'))
        self.assertFalse(PW.is_complete('ec042', 'SUCUTIL-PWY'))

        # UNIQUE = TRUE
        self.assertFalse(PW.is_complete('HS', 'PWY-6538', unique=True))

    def test_get_pw_present(self):
        # NO PARAMETER
        dict_pw_pres = PW.get_pw_present()
        self.assertIs(type(dict_pw_pres), dict)
        self.assertEqual(dict_pw_pres['HS'][0], 1252)
        self.assertIs(type(dict_pw_pres['HS'][1]), set)
        self.assertIn('PWY-5087', dict_pw_pres['HS'][1])
        self.assertEqual(set(dict_pw_pres.keys()), set(PW.species_list))

        # UNIQUE = TRUE
        dict_pw_pres_uni = PW.get_pw_present(unique=True)
        self.assertEqual(dict_pw_pres_uni['HS'][0], 0)
        self.assertEqual(dict_pw_pres_uni['IAI1'][0], 1)
        self.assertIn('PWY-5934', dict_pw_pres_uni['IAI1'][1])

        # SPECIES STR
        dict_pw_pres = PW.get_pw_present(species='HS')
        self.assertEqual(dict_pw_pres.keys(), {'HS'})

        # SPECIES LIST
        dict_pw_pres = PW.get_pw_present(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_pw_pres.keys(), {'HS', 'IAI1', 'ec042'})

    def test_get_pw_max(self):
        # NO PARAMETER
        dict_pw_max = PW.get_pw_max()
        self.assertIs(type(dict_pw_max), dict)
        self.assertEqual(dict_pw_max['HS'][0], 1217)
        self.assertIs(type(dict_pw_max['HS'][1]), set)
        self.assertIn('PWY0-901', dict_pw_max['HS'][1])
        self.assertEqual(set(dict_pw_max.keys()), set(PW.species_list))

        # UNIQUE = TRUE
        dict_pw_max_uni = PW.get_pw_max(unique=True)
        self.assertEqual(dict_pw_max_uni['HS'][0], 0)
        self.assertEqual(dict_pw_max_uni['IAI1'][0], 3)
        self.assertIn('PWY-6854', dict_pw_max_uni['IAI1'][1])

        # SPECIES STR
        dict_pw_max = PW.get_pw_max(species='HS')
        self.assertEqual(dict_pw_max.keys(), {'HS'})

        # SPECIES LIST
        dict_pw_max = PW.get_pw_max(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_pw_max.keys(), {'HS', 'IAI1', 'ec042'})

    def test_get_pw_complete(self):
        # NO PARAMETER
        dict_pw_complete = PW.get_pw_complete()
        self.assertIs(type(dict_pw_complete), dict)
        self.assertEqual(dict_pw_complete['HS'][0], 391)
        self.assertIs(type(dict_pw_complete['HS'][1]), set)
        self.assertIn('PWY-6587', dict_pw_complete['HS'][1])
        self.assertEqual(set(dict_pw_complete.keys()), set(PW.species_list))

        # UNIQUE = TRUE
        dict_pw_complete_uni = PW.get_pw_complete(unique=True)
        self.assertEqual(dict_pw_complete_uni['HS'][0], 0)
        self.assertEqual(dict_pw_complete_uni['IAI1'][0], 0)

        # SPECIES STR
        dict_pw_complete = PW.get_pw_complete(species='HS')
        self.assertEqual(dict_pw_complete.keys(), {'HS'})

        # SPECIES LIST
        dict_pw_complete = PW.get_pw_complete(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_pw_complete.keys(), {'HS', 'IAI1', 'ec042'})

    # Loss

    def test_is_absent(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_absent('HS', 'PWY-5338'))
        self.assertFalse(PW.is_absent('UTI89', 'PWY-5338'))

        # UNIQUE = TRUE
        self.assertTrue(PW.is_absent('ec042', 'PWY-6784', unique=True))
        self.assertFalse(PW.is_absent('UTI89', 'PWY-5647', unique=True))

    def test_is_min(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_min('IAI1', 'PWY-6370'))
        self.assertFalse(PW.is_min('HS', 'PWY-7007'))

        # UNIQUE = TRUE
        self.assertTrue(PW.is_min('sf301', 'PWY-7799', unique=True))
        self.assertFalse(PW.is_min('ec042', 'PWY-6538', unique=True))

    def test_is_incomplete(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(PW.is_incomplete('HS', 'PWY-8171'))
        self.assertFalse(PW.is_incomplete('HS', 'ORNDEG-PWY'))

        # UNIQUE = TRUE
        self.assertTrue(PW.is_incomplete('sf301', 'PWY-1801', unique=True))
        self.assertFalse(PW.is_incomplete('CFT073', 'CYANCAT-PWY', unique=True))

    def test_get_pw_absent(self):
        # NO PARAMETER
        dict_pw_abs = PW.get_pw_absent()
        self.assertIs(type(dict_pw_abs), dict)
        self.assertEqual(dict_pw_abs['HS'][0], 29)
        self.assertIs(type(dict_pw_abs['HS'][1]), set)
        self.assertIn('GDPRHAMSYN-PWY', dict_pw_abs['HS'][1])
        self.assertEqual(set(dict_pw_abs.keys()), set(PW.species_list))

        # UNIQUE = TRUE
        dict_pw_abs_uni = PW.get_pw_absent(unique=True)
        self.assertEqual(dict_pw_abs_uni['HS'][0], 10)
        self.assertEqual(dict_pw_abs_uni['IAI1'][0], 0)
        self.assertIn('PWY-8225', dict_pw_abs_uni['HS'][1])

        # SPECIES STR
        dict_pw_abs = PW.get_pw_present(species='HS')
        self.assertEqual(dict_pw_abs.keys(), {'HS'})

        # SPECIES LIST
        dict_pw_abs = PW.get_pw_present(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_pw_abs.keys(), {'HS', 'IAI1', 'ec042'})

    def test_get_pw_min(self):
        # NO PARAMETER
        dict_pw_min = PW.get_pw_min()
        self.assertIs(type(dict_pw_min), dict)
        self.assertEqual(dict_pw_min['HS'][0], 1058)
        self.assertIs(type(dict_pw_min['HS'][1]), set)
        self.assertIn('TCA-1', dict_pw_min['HS'][1])
        self.assertEqual(set(dict_pw_min.keys()), set(PW.species_list))

        # UNIQUE = TRUE
        dict_pw_min_uni = PW.get_pw_min(unique=True)
        self.assertEqual(dict_pw_min_uni['HS'][0], 5)
        self.assertEqual(dict_pw_min_uni['IAI1'][0], 0)
        self.assertIn('PWY-5177', dict_pw_min_uni['HS'][1])

        # SPECIES STR
        dict_pw_min = PW.get_pw_min(species='HS')
        self.assertEqual(dict_pw_min.keys(), {'HS'})

        # SPECIES LIST
        dict_pw_min = PW.get_pw_min(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_pw_min.keys(), {'HS', 'IAI1', 'ec042'})

    def test_get_pw_incomplete(self):
        # NO PARAMETER
        dict_pw_incomplete = PW.get_pw_incomplete()
        self.assertIs(type(dict_pw_incomplete), dict)
        self.assertEqual(dict_pw_incomplete['HS'][0], 861)
        self.assertIs(type(dict_pw_incomplete['HS'][1]), set)
        self.assertIn('PWY0-1591', dict_pw_incomplete['HS'][1])
        self.assertEqual(set(dict_pw_incomplete.keys()), set(PW.species_list))

        # UNIQUE = TRUE
        dict_pw_incomplete_uni = PW.get_pw_incomplete(unique=True)
        self.assertEqual(dict_pw_incomplete_uni['HS'][0], 22)
        self.assertEqual(dict_pw_incomplete_uni['IAI1'][0], 20)
        self.assertIn('PWY-5690', dict_pw_incomplete_uni['HS'][1])

        # SPECIES STR
        dict_pw_incomplete = PW.get_pw_incomplete(species='HS')
        self.assertEqual(dict_pw_incomplete.keys(), {'HS'})

        # SPECIES LIST
        dict_pw_incomplete = PW.get_pw_incomplete(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_pw_incomplete.keys(), {'HS', 'IAI1', 'ec042'})

    # Other

    def test_get_pw_over_treshold(self):
        hs_ov50 = PW.get_pw_over_treshold('HS', 0.5)['HS']
        self.assertEqual(hs_ov50[0], 781)
        self.assertIn('PWY-5278', hs_ov50[1])  # 1/2
        self.assertIn('PWY0-41', hs_ov50[1])  # 5/6
        self.assertNotIn('PWY-6523', hs_ov50[1])  # 1/3
        hs_ov50 = PW.get_pw_over_treshold('HS', 0.5, strict=True)['HS']
        self.assertNotIn('PWY-5278', hs_ov50[1])  # 1/2
        self.assertIn('PWY0-41', hs_ov50[1])  # 5/6
        self.assertNotIn('PWY-6523', hs_ov50[1])  # 1/3

    def test_get_pw_under_treshold(self):
        hs_ov50 = PW.get_pw_under_treshold('HS', 0.5)['HS']
        self.assertEqual(hs_ov50[0], 650)
        self.assertIn('PWY-5278', hs_ov50[1])  # 1/2
        self.assertNotIn('PWY0-41', hs_ov50[1])  # 5/6
        self.assertIn('PWY-6523', hs_ov50[1])  # 1/3
        hs_ov50 = PW.get_pw_under_treshold('HS', 0.5, strict=True)['HS']
        self.assertNotIn('PWY-5278', hs_ov50[1])  # 1/2
        self.assertNotIn('PWY0-41', hs_ov50[1])  # 5/6
        self.assertIn('PWY-6523', hs_ov50[1])  # 1/3

    def test_get_pw_names(self):
        pw_names = PW.get_pw_names()
        self.assertEqual(pw_names['1CMET2-PWY'], 'folate transformations III (E. coli)')
        self.assertEqual(pw_names['PWY-5137'], 'fatty acid &beta;-oxidation III (unsaturated, odd number)')
        self.assertEqual(pw_names.keys(), set(PW.pathways_list))

    def test_convert_df_to_binary(self):
        # ONLY THRESHOLD
        df_bin = PW.convert_df_to_binary(0.5)
        self.assertEqual(df_bin.loc['PWY66-385', 'HS'], 1)
        self.assertEqual(df_bin.loc['PWY-8219', 'HS'], 0)
        self.assertEqual(list(df_bin.columns), PW.species_list)
        self.assertEqual(list(df_bin.index), PW.pathways_list)
        # TRSHOLD STRICT
        df_bin = PW.convert_df_to_binary(0.5, strict=True)
        self.assertEqual(df_bin.loc['PWY-6578', 'HS'], 0)  # 1/2
        # WITH COMMON NAME
        df_bin = PW.convert_df_to_binary(0.5, common_name=True)
        self.assertEqual(df_bin.columns[0], 'common name')
        # WITH NB RNX IN PW
        df_bin = PW.convert_df_to_binary(0.5, nb_rxn=True)
        self.assertEqual(df_bin.columns[-1], 'nb rnx pw')
        # WITH COMMON NAME AND NB RNX IN PW
        df_bin = PW.convert_df_to_binary(0.5, nb_rxn=True, common_name=True)
        self.assertEqual(df_bin.columns[0], 'common name')
        self.assertEqual(df_bin.columns[-1], 'nb rnx pw')

    def test_generate_pw_dendrogram(self):
        pass
        # FILE Study_folder/output_data/dendro_tanglegrams/dendro_groups.tsv FILLED BEFORE
        # NO PHYLO REF FILE
        # PW.generate_pw_dendrogram(0.8, name="test_no_phylo", n_boot=10)
        # self.assertEqual(os.path.exists(os.path.join(PW.path_study, 'output_data', 'dendro_tanglegrams', 'bact7',
        #                                              'pw_test_no_phylo', 'pw_test_no_phylo_dendextend_dend.png')),
        #                  True)
        # self.assertEqual(os.path.exists(os.path.join(PW.path_study, 'output_data', 'dendro_tanglegrams', 'bact7',
        #                                              'pw_test_no_phylo', 'pw_test_no_phylo_pvclust_dend.png')),
        #                  True)

        # PHYLO REF FILE
        # phylo_file = os.path.join(A.path_study, 'Phylo.nexus')
        # R.generate_rnx_dendrogram(name="test_phylo", n_boot=10, phylo_file=phylo_file)



