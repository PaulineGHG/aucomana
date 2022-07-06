import os.path
import unittest
import pandas as pd
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
        # SPECIES_LIST = ['HS', 'UTI89', 'ec042']
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
        # OUT = 1
        ro1 = A.reactions(RUN, out=1)
        self.assertEqual(2196, len(ro1.reactions_list))
        self.assertIn('RXN-17824', ro1.reactions_list)
        self.assertIn('BETAGALACTOSID-RXN', ro1.reactions_list)
        self.assertNotIn('RXN-10016', ro1.reactions_list)

        # OUT = 2
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
        columns_list = [x + R.STR_GENE_ASSOC for x in R.species_list]
        self.assertListEqual(columns_list, list(R.data_genes_assoc.columns))

    def test_is_present(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(R.is_present('UTI89', 'RXN-16632'))
        self.assertFalse(R.is_present('UTI89', 'RXN0-963'))

        # UNIQUE = TRUE
        self.assertTrue(R.is_present('sf301', 'RXN-14500', unique=True))
        self.assertFalse(R.is_present('IAI1', 'RXNMETA-12672', unique=True))

    def test_is_absent(self):
        # NO PARAMETER (UNIQUE = FALSE)
        self.assertTrue(R.is_absent('CFT073', 'RXN0-962'))
        self.assertFalse(R.is_absent('CFT073', 'RXN-12615'))

        # UNIQUE = TRUE
        self.assertTrue(R.is_absent('ec042', '5.1.99.4-RXN', unique=True))
        self.assertFalse(R.is_absent('ec042', 'RXN-12615', unique=True))

    def test_get_rxn_present(self):
        # NO PARAMETER
        dict_rxn_pres = R.get_rxn_present()
        self.assertIs(type(dict_rxn_pres), dict)
        self.assertEqual(dict_rxn_pres['HS'][0], 2304)
        self.assertIs(type(dict_rxn_pres['HS'][1]), set)
        self.assertIn('RXN-15149', dict_rxn_pres['HS'][1])
        self.assertEqual(set(dict_rxn_pres.keys()), set(R.species_list))

        # UNIQUE = TRUE
        dict_rxn_pres_uni = R.get_rxn_present(unique=True)
        self.assertEqual(dict_rxn_pres_uni['HS'][0], 0)
        self.assertEqual(dict_rxn_pres_uni['IAI1'][0], 5)
        self.assertIn('RXN-20706', dict_rxn_pres_uni['IAI1'][1])

        # SPECIES STR
        dict_rxn_pres = R.get_rxn_present(species='HS')
        self.assertEqual(dict_rxn_pres.keys(), {'HS'})

        # SPECIES LIST
        dict_rxn_pres = R.get_rxn_present(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_rxn_pres.keys(), {'HS', 'IAI1', 'ec042'})

    def test_get_rxn_absent(self):
        # NO PARAMETER
        dict_rxn_abs = R.get_rxn_absent()
        self.assertIs(type(dict_rxn_abs), dict)
        self.assertEqual(dict_rxn_abs['HS'][0], 98)
        self.assertIs(type(dict_rxn_abs['HS'][1]), set)
        self.assertIn('TRANS-RXN0-268', dict_rxn_abs['HS'][1])
        self.assertEqual(set(dict_rxn_abs.keys()), set(R.species_list))

        # UNIQUE = TRUE
        dict_rxn_abs_uni = R.get_rxn_absent(unique=True)
        self.assertEqual(dict_rxn_abs_uni['HS'][0], 15)
        self.assertEqual(dict_rxn_abs_uni['IAI1'][0], 0)
        self.assertIn('RXN0-7341', dict_rxn_abs_uni['HS'][1])

        # SPECIES STR
        dict_rxn_pres = R.get_rxn_present(species='HS')
        self.assertEqual(dict_rxn_pres.keys(), {'HS'})

        # SPECIES LIST
        dict_rxn_pres = R.get_rxn_present(species=['HS', 'IAI1', 'ec042'])
        self.assertEqual(dict_rxn_pres.keys(), {'HS', 'IAI1', 'ec042'})

    def test_get_genes_assoc(self):
        # NO PARAMETER
        dict_assoc = R.get_genes_assoc()
        self.assertEqual(dict_assoc.keys(), set(R.reactions_list))
        self.assertIs(type(dict_assoc), dict)
        self.assertIs(type(dict_assoc['3.2.2.21-RXN']), dict)
        self.assertEqual(dict_assoc['3.2.2.21-RXN'].keys(), set(R.species_list))
        self.assertIs(type(dict_assoc['3.2.2.21-RXN']['HS']), list)
        self.assertIn('EcHS_A2209', dict_assoc['3.2.2.21-RXN']['HS'])

        # REACTIONS_LIST = ['3.2.2.21-RXN', 'RXN-16632'] (LIST)
        dict_assoc = R.get_genes_assoc(reactions_list=['3.2.2.21-RXN', 'RXN-16632'])
        self.assertEqual(dict_assoc.keys(), {'3.2.2.21-RXN', 'RXN-16632'})

        # REACTIONS_LIST = '3.2.2.21-RXN' (STR)
        dict_assoc = R.get_genes_assoc(reactions_list='3.2.2.21-RXN')
        self.assertEqual(dict_assoc.keys(), {'3.2.2.21-RXN'})

        # OUTPUT_FILE = TRUE
        R.get_genes_assoc(reactions_list=['3.2.2.21-RXN', 'RXN-16632'], output_file=True)
        self.assertEqual(os.path.exists(os.path.join(R.path_study, "output_data", "reactions_data", "genes_assoc",
                                        str(R.nb_genes_assoc), '3.2.2.21-RXN.fa')), True)
        self.assertEqual(os.path.exists(os.path.join(R.path_study, "output_data", "reactions_data", "genes_assoc",
                                        str(R.nb_genes_assoc), 'RXN-16632.fa')), True)
        self.assertEqual(os.path.exists(os.path.join(R.path_study, "output_data", "reactions_data", "genes_assoc",
                                        str(R.nb_genes_assoc), 'bact7_genes_assoc.json')), True)

    def test_get_reactions_names(self):
        rxn_names = R.get_reactions_names()
        self.assertEqual(rxn_names['GDPMANDEHYDRA-RXN'], 'GDP-D-mannose dehydratase, NAD(P)-binding')
        self.assertEqual(rxn_names['RXN-20756'], 'tRNA 2-selenouridine synthase')

    def test_generate_rnx_dendrogram(self):
        # TODO : FINISH TEST + CREATE PHYLO FILE
        # FILE Study_folder/output_data/dendro_tanglegrams/dendro_groups.tsv FILLED BEFORE
        # NO PHYLO REF FILE
        R.generate_rnx_dendrogram(name="test_no_phylo", n_boot=10)
        self.assertEqual(os.path.exists(os.path.join(R.path_study, 'output_data', 'dendro_tanglegrams', 'bact7',
                                                     'rnx_test_no_phylo', 'rnx_test_no_phylo_dendextend_dend.png')),
                         True)
        self.assertEqual(os.path.exists(os.path.join(R.path_study, 'output_data', 'dendro_tanglegrams', 'bact7',
                                                     'rnx_test_no_phylo', 'rnx_test_no_phylo_pvclust_dend.png')),
                         True)

        # PHYLO REF FILE
        # phylo_file = os.path.join(A.path_study, 'Phylo.nexus')
        # R.generate_rnx_dendrogram(name="test_phylo", n_boot=10, phylo_file=phylo_file)


if __name__ == '__main__':
    unittest.main()





