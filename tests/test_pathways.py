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