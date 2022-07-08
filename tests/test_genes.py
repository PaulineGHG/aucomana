import unittest
import pandas as pd
from analysis_runs.analysis import Analysis
from analysis_runs.genes import Genes

STUDY_PATH = 'Study_folder'
RUNS_PATH = 'Runs_aucome'

RUN = 'bact7'
A = Analysis(RUNS_PATH, STUDY_PATH)
G = A.genes(RUN)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(G, Genes)

    def test_init_species_default(self):
        species_obj = ['HS', 'IAI1', 'CFT073', 'UTI89', 'sf301', 'ec042', 'LF82']
        self.assertEqual(species_obj, G.species_list)

    def test_species_precised(self):
        # SPECIES_LIST = ['HS', 'UTI89', 'ec042']
        sp_l = ['HS', 'UTI89', 'ec042']
        gs = A.genes(RUN, species_list=sp_l)
        self.assertEqual(sp_l, gs.species_list)

    def test_name(self):
        self.assertEqual(G.name, 'bact7')

    def test_genes(self):
        self.assertIs(type(G.genes_list), list)
        self.assertEqual(12115, len(G.genes_list))
        self.assertEqual('EcHS_A3606', G.genes_list[0])
        self.assertEqual('LF82_1397', G.genes_list[-1])

    def test_data_genes(self):
        self.assertIsInstance(G.data_genes, pd.DataFrame)
        self.assertListEqual(G.genes_list, list(G.data_genes.index))
        self.assertListEqual(G.species_list, list(G.data_genes.columns))
        self.assertEqual(G.data_genes.loc['EcHS_A3031', 'HS'], 1)
        self.assertEqual(G.data_genes.loc['EcHS_A3031', 'LF82'], 0)
        self.assertEqual(G.data_genes.loc['LF82_1690', 'LF82'], 1)

    def test_data_rxn_assoc(self):
        self.assertIsInstance(G.data_rxn_assoc, pd.DataFrame)
        columns_list = [x + G.STR_RXN_ASSOC for x in G.species_list]
        self.assertListEqual(G.genes_list, list(G.data_rxn_assoc.index))
        self.assertListEqual(columns_list, list(G.data_rxn_assoc.columns))
        self.assertEqual(G.data_rxn_assoc.loc['EcHS_A3031', 'HS' + G.STR_RXN_ASSOC], '4.3.1.15-RXN')
        self.assertEqual(G.data_rxn_assoc.loc['LF82_1690', 'LF82' + G.STR_RXN_ASSOC], 'ABC-25-RXN')

    def test_nb_genes_species(self):
        self.assertIs(type(G.nb_genes_species), dict)
        self.assertEqual(G.nb_genes_species.keys(), set(G.species_list))
        self.assertEqual(G.nb_genes_species['HS'], 1691)

    def test_get_genes_species(self):
        # NO PARAMETER
        genes_species = G.get_genes_species()
        self.assertIs(type(genes_species), dict)
        self.assertIs(type(genes_species['HS']), set)
        self.assertEqual(genes_species.keys(), set(G.species_list))
        self.assertIn('EcHS_A3369', genes_species['HS'])
        self.assertEqual(len(genes_species['HS']), 1691)

        # SPECIES_LIST STR
        self.assertEqual(G.get_genes_species(species_list='HS').keys(), {'HS'})

        # SPECIES_LIST STR
        self.assertEqual(G.get_genes_species(species_list=['HS', 'UTI89']).keys(), {'HS', 'UTI89'})

