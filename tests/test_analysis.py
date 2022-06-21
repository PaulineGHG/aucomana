import unittest
from analysis_runs.analysis import Analysis

STUDY_PATH = "Study_folder"
RUNS_PATH = "Runs_aucome"
RUN = "bact7"

A = Analysis(RUNS_PATH, STUDY_PATH)
print(A.get_grp_l(RUN, "groupA"))


# class Test(unittest.TestCase):
#
#     def test_init(self):
#         self.assertIsInstance(R, Reactions)