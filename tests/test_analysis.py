import unittest
from analysis_runs.analysis import Analysis

STUDY_PATH = "Study_folder"
RUNS_PATH = "Runs_aucome"
RUN = "bact7"

A = Analysis(RUNS_PATH, STUDY_PATH)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(A, Analysis)

    def test_get_grp_l(self):
        self.assertEqual(A.get_grp_set(RUN, "group1"), {'UTI89', 'IAI1', 'HS', 'CFT073'})
        self.assertEqual(A.get_grp_set(RUN, "group1", ['UTI89', 'HS', 'LF82']), {'UTI89', 'HS'})
        self.assertEqual(A.get_grp_set(RUN, "group2"), {'UTI89', 'LF82', 'ec042', 'sf301'})
        self.assertEqual(A.get_grp_set(RUN, "groupA"), {'IAI1', 'HS'})
        self.assertEqual(A.get_grp_set(RUN, "groupB"), {'sf301', 'UTI89', 'CFT073'})
        self.assertEqual(A.get_grp_set(RUN, "groupC"), {'LF82', 'ec042', 'sf301'})
        with self.assertRaises(ValueError):
            A.get_grp_set(RUN, "group3")


