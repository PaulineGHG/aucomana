import unittest
from analysis_runs.analysis import Analysis
from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.genes import Genes
from analysis_runs.metabolites import Metabolites

STUDY_PATH = 'Study_folder'
RUNS_PATH = 'Runs_aucome'
RUN = 'bact7'

A = Analysis(RUNS_PATH, STUDY_PATH)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(A, Analysis)

    def test_get_grp_l(self):
        self.assertEqual(A.get_grp_set(RUN, 'group1'), {'UTI89', 'IAI1', 'HS', 'CFT073'})
        self.assertEqual(A.get_grp_set(RUN, 'group1', ['UTI89', 'HS', 'LF82']), {'UTI89', 'HS'})
        self.assertEqual(A.get_grp_set(RUN, 'group2'), {'UTI89', 'LF82', 'ec042', 'sf301'})
        self.assertEqual(A.get_grp_set(RUN, 'groupA'), {'IAI1', 'HS'})
        self.assertEqual(A.get_grp_set(RUN, 'groupB'), {'sf301', 'UTI89', 'CFT073'})
        self.assertEqual(A.get_grp_set(RUN, 'groupC'), {'LF82', 'ec042', 'sf301'})
        with self.assertRaises(ValueError):
            A.get_grp_set(RUN, 'group3')

    def test_reactions(self):
        RXN = A.reactions(RUN)
        self.assertIsInstance(RXN, Reactions)
        RXN = A.reactions(RUN, species_list=['IAI1', 'HS', 'sf301'])
        self.assertEqual(RXN.species_list, ['IAI1', 'HS', 'sf301'])
        RXN = A.reactions(RUN, species_list=['IAI1', 'HS', 'sf301'], group='groupA')
        self.assertEqual(set(RXN.species_list), {'IAI1', 'HS'})
        RXN = A.reactions(RUN, out=1)
        self.assertEqual(len(RXN.reactions_list), 2196)
        with self.assertRaises(OSError):
            A.reactions('run_not_existing')

    def test_pathways(self):
        PW = A.pathways(RUN)
        self.assertIsInstance(PW, Pathways)
        PW = A.pathways(RUN, species_list=['IAI1', 'HS', 'sf301'])
        self.assertEqual(PW.species_list, ['IAI1', 'HS', 'sf301'])
        PW = A.pathways(RUN, species_list=['IAI1', 'HS', 'sf301'], group='groupA')
        self.assertEqual(set(PW.species_list), {'IAI1', 'HS'})
        PW = A.pathways(RUN, out=1)
        self.assertEqual(len(PW.pathways_list), 1225)
        PW = A.pathways(RUN, nb_rxn_pw_min=3)
        self.assertEqual(len(PW.pathways_list), 993)
        with self.assertRaises(OSError):
            A.pathways('run_not_existing')

    def test_genes(self):
        GN = A.genes(RUN)
        self.assertIsInstance(GN, Genes)
        GN = A.genes(RUN, species_list=['IAI1', 'HS', 'sf301'])
        self.assertEqual(GN.species_list, ['IAI1', 'HS', 'sf301'])
        GN = A.genes(RUN, species_list=['IAI1', 'HS', 'sf301'], group='groupA')
        self.assertEqual(set(GN.species_list), {'IAI1', 'HS'})
        with self.assertRaises(OSError):
            A.genes('run_not_existing')

    def test_metabolites(self):
        MB = A.metabolites(RUN)
        self.assertIsInstance(MB, Metabolites)
        MB = A.metabolites(RUN, species_list=['IAI1', 'HS', 'sf301'])
        self.assertEqual(MB.species_list, ['IAI1', 'HS', 'sf301'])
        MB = A.metabolites(RUN, species_list=['IAI1', 'HS', 'sf301'], group='groupA')
        self.assertEqual(set(MB.species_list), {'IAI1', 'HS'})
        with self.assertRaises(OSError):
            A.metabolites('run_not_existing')

    def test_compare_groups(self):
        pass




