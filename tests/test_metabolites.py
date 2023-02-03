import os.path
import unittest
import pandas as pd
from aucomana.aucomana import Analysis
from aucomana.metabolites import Metabolites

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

    def test_reactions_consuming(self):
        # Test types
        self.assertIs(type(MB.reactions_consuming), dict)
        self.assertEqual(MB.reactions_consuming.keys(), set(MB.metabolites_list))
        self.assertIs(type(MB.reactions_consuming['DETHIOBIOTIN']), dict)
        self.assertIs(type(MB.reactions_consuming['DETHIOBIOTIN']['HS']), list)
        self.assertIs(type(MB.reactions_consuming['DETHIOBIOTIN']['HS'][0]), str)
        # Test values
        self.assertEqual(MB.reactions_consuming['GMP']['HS'], ['GUANYL-KIN-RXN', 'GMP-REDUCT-RXN', 'RXN-7609'])
        self.assertEqual(MB.reactions_consuming['GLYCOLLATE']['HS'], ['GLYCOLATEDEHYDRO-RXN', 'RXN-969', 'RXN0-7229'])
        self.assertEqual(MB.reactions_consuming['DNA-Ligase-L-lysine-adenylate']['LF82'], ['RXN-17918'])

    def test_reactions_producing(self):
        print(MB.reactions_producing)
        # Test types
        self.assertIs(type(MB.reactions_producing), dict)
        self.assertEqual(MB.reactions_producing.keys(), set(MB.metabolites_list))
        self.assertIs(type(MB.reactions_producing['DETHIOBIOTIN']), dict)
        self.assertIs(type(MB.reactions_producing['DETHIOBIOTIN']['HS']), list)
        self.assertIs(type(MB.reactions_producing['DETHIOBIOTIN']['HS'][0]), str)
        # Test values
        self.assertIn('RXN-19302', MB.reactions_producing['GMP']['HS'])
        self.assertEqual(MB.reactions_producing['Alkylated-Bases']['HS'], ['3.2.2.21-RXN'])
        self.assertEqual(MB.reactions_producing['RNA-3prime-Cytidine-3prime-P']['CFT073'],
                         ['RXN-19935', 'RXN-19932'])

    def test_metabolites_produced(self):
        self.assertIsInstance(MB.metabolites_produced, pd.DataFrame)
        self.assertListEqual(list(MB.metabolites_produced.index), MB.metabolites_list)
        self.assertListEqual(list(MB.metabolites_produced.columns), MB.species_list)
        self.assertEqual(MB.metabolites_produced.loc['5-HYDROXYU34-TRNA', 'HS'], 1)
        self.assertEqual(MB.metabolites_produced.loc['H2CO3', 'IAI1'], 0)

    def test_metabolites_consumed(self):
        self.assertIsInstance(MB.metabolites_consumed, pd.DataFrame)
        self.assertListEqual(list(MB.metabolites_consumed.index), MB.metabolites_list)
        self.assertListEqual(list(MB.metabolites_consumed.columns), MB.species_list)
        self.assertEqual(MB.metabolites_consumed.loc['5-HYDROXYU34-TRNA', 'HS'], 0)
        self.assertEqual(MB.metabolites_consumed.loc['H2CO3', 'IAI1'], 1)

    def test_is_consumed(self):
        self.assertTrue(MB.is_consumed('HS', 'OLEATE-CPD'))
        self.assertFalse(MB.is_consumed('HS', 'Methylketones'))
        self.assertTrue(MB.is_consumed('IAI1', 'CPD-479', unique=True))
        self.assertFalse(MB.is_consumed('HS', 'PHENYLACETALDEHYDE', unique=True))

    def test_is_produced(self):
        self.assertTrue(MB.is_produced('HS', 'D-ALA-D-ALA'))
        self.assertFalse(MB.is_produced('HS', 'Lipoyl-ACPs'))
        self.assertTrue(MB.is_produced('ec042', 'MESACONATE', unique=True))
        self.assertFalse(MB.is_produced('HS', 'GMP', unique=True))

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

