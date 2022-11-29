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

    # def test_data_metabolites_consumed(self):
    #     self.assertIsInstance(MB.data_metabolites_consumed, pd.DataFrame)
    #     columns_list = [x + MB.STR_CONSUME for x in MB.species_list]
    #     self.assertListEqual(MB.metabolites_list, list(MB.data_metabolites_consumed.index))
    #     self.assertListEqual(columns_list, list(MB.data_metabolites_consumed.columns))
    #     self.assertEqual(MB.data_metabolites_consumed.loc['GLYCOLLATE', 'HS_rxn_consume'],
    #                      'GLYCOLATEDEHYDRO-RXN;RXN-969;RXN0-7229')
    #     self.assertEqual(MB.data_metabolites_consumed.loc['DNA-Ligase-L-lysine-adenylate', 'LF82_rxn_consume'],
    #                      'RXN-17918')
    #
    # def test_data_metabolites_produced(self):
    #     self.assertIsInstance(MB.data_metabolites_produced, pd.DataFrame)
    #     columns_list = [x + MB.STR_PRODUCE for x in MB.species_list]
    #     self.assertListEqual(MB.metabolites_list, list(MB.data_metabolites_produced.index))
    #     self.assertListEqual(columns_list, list(MB.data_metabolites_produced.columns))
    #     self.assertEqual(MB.data_metabolites_produced.loc['Alkylated-Bases', 'HS_rxn_produce'], '3.2.2.21-RXN')
    #     self.assertEqual(MB.data_metabolites_produced.loc['RNA-3prime-Cytidine-3prime-P', 'CFT073_rxn_produce'],
    #                      'RXN-19935;RXN-19932')
    #
    # def test_data_metabolites(self):
    #     self.assertIsInstance(MB.data_metabolites, pd.DataFrame)
    #     self.assertListEqual(MB.metabolites_list, list(MB.data_metabolites.index))
    #     self.assertListEqual(MB.species_list, list(MB.data_metabolites.columns))
    #     self.assertEqual(MB.data_metabolites.loc['GDP-4-DEHYDRO-6-L-DEOXYGALACTOSE', 'HS'], 0)
    #     self.assertEqual(MB.data_metabolites.loc['GLN', 'LF82'], 1)
    #
    # def test_get_rxn_consuming(self):
    #     # NO PARAMETER
    #     rxn_consuming = MB.get_rxn_consuming()
    #     self.assertEqual(rxn_consuming.keys(), set(MB.metabolites_list))
    #     self.assertEqual(rxn_consuming['CPD-15158'].keys(), set(MB.species_list))
    #     self.assertIn('RXN-969', rxn_consuming['GLYCOLLATE']['HS'])
    #
    #     # METABOLITES_LIST STR
    #     rxn_consuming = MB.get_rxn_consuming(metabolites_list='GLYCOLLATE')
    #     self.assertEqual(rxn_consuming.keys(), {'GLYCOLLATE'})
    #     self.assertEqual(rxn_consuming['GLYCOLLATE'].keys(), set(MB.species_list))
    #     self.assertIn('RXN-969', rxn_consuming['GLYCOLLATE']['HS'])
    #
    #     # METABOLITES_LIST LIST
    #     rxn_consuming = MB.get_rxn_consuming(metabolites_list=['GLYCOLLATE', 'CPD-15158'])
    #     self.assertEqual(rxn_consuming.keys(), {'GLYCOLLATE', 'CPD-15158'})
    #     self.assertEqual(rxn_consuming['GLYCOLLATE'].keys(), set(MB.species_list))
    #     self.assertIn('RXN-969', rxn_consuming['GLYCOLLATE']['HS'])
    #
    # def test_get_rxn_producing(self):
    #     # NO PARAMETER
    #     rxn_producing = MB.get_rxn_producing()
    #     self.assertEqual(rxn_producing.keys(), set(MB.metabolites_list))
    #     self.assertEqual(rxn_producing['CPD-15158'].keys(), set(MB.species_list))
    #     self.assertIn('RXN-20148', rxn_producing['OXYGEN-MOLECULE']['CFT073'])
    #
    #     # METABOLITES_LIST STR
    #     rxn_producing = MB.get_rxn_producing(metabolites_list='OXYGEN-MOLECULE')
    #     self.assertEqual(rxn_producing.keys(), {'OXYGEN-MOLECULE'})
    #     self.assertEqual(rxn_producing['OXYGEN-MOLECULE'].keys(), set(MB.species_list))
    #     self.assertIn('RXN-20148', rxn_producing['OXYGEN-MOLECULE']['CFT073'])
    #
    #     # METABOLITES_LIST LIST
    #     rxn_producing = MB.get_rxn_producing(metabolites_list=['GLYCOLLATE', 'OXYGEN-MOLECULE'])
    #     self.assertEqual(rxn_producing.keys(), {'GLYCOLLATE', 'OXYGEN-MOLECULE'})
    #     self.assertEqual(rxn_producing['GLYCOLLATE'].keys(), set(MB.species_list))
    #     self.assertIn('RXN-20148', rxn_producing['OXYGEN-MOLECULE']['CFT073'])

    def test_get_metabolites_names(self):
        mb_names = MB.get_metabolites_names()
        self.assertEqual(mb_names['CPD-12117'], 'demethylmenaquinol-7')
        self.assertEqual(mb_names['CPD-12817'], '1,2-dipalmitoyl-phosphatidylserine')

    def test_generate_met_dendrogram(self):
        # TODO : FINISH TEST
        pass
        # FILE Study_folder/output_data/dendro_tanglegrams/dendro_groups.tsv FILLED BEFORE
        # NO PHYLO REF FILE
        MB.generate_met_dendrogram(name="test_no_phylo", n_boot=10)
        self.assertEqual(os.path.exists(os.path.join(MB.path_study, 'output_data', 'dendro_tanglegrams', 'bact7',
                                                     'met_test_no_phylo', 'met_test_no_phylo_dendextend_dend.png')),
                         True)
        self.assertEqual(os.path.exists(os.path.join(MB.path_study, 'output_data', 'dendro_tanglegrams', 'bact7',
                                                     'met_test_no_phylo', 'met_test_no_phylo_pvclust_dend.png')),
                         True)

        # PHYLO REF FILE
        # phylo_file = os.path.join(A.path_study, 'Phylo.nexus')
        # MB.generate_met_dendrogram(name="test_phylo", n_boot=10, phylo_file=phylo_file)

