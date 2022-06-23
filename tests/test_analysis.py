import unittest
from analysis_runs.analysis import Analysis
from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways

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

    def test_pathways(self):
        PW = A.pathways(RUN)
        self.assertIsInstance(PW, Pathways)
        PW = A.pathways(RUN, species_list=['IAI1', 'HS', 'sf301'])
        self.assertEqual(PW.species_list, ['IAI1', 'HS', 'sf301'])
        PW = A.pathways(RUN, species_list=['IAI1', 'HS', 'sf301'], group='groupA')
        self.assertEqual(set(PW.species_list), {'IAI1', 'HS'})
        PW = A.pathways(RUN, out=1)

