import unittest
import os
import analysis_runs.init_analysis as ia

STUDY_PATH = "Study_folder"
RUNS_PATH = "Runs_aucome"

ia.create_folders(STUDY_PATH)

arbo = [['input_data'], ['output_data'],
        ['output_data', 'dendro_tanglegrams'],
        ['output_data', 'pathways_data'],
        ['output_data', 'reactions_data'],
        ['output_data', 'reactions_data', 'common_reac'],
        ['output_data', 'reactions_data', 'common_reac', 'union'],
        ['output_data', 'reactions_data', 'common_reac', 'intersection'],
        ['output_data', 'reactions_data', 'genes_assoc'],
        ['output_data', 'genes_data'],
        ['output_data', 'metabolites_data'],
        ['output_data', 'pages']]


class Test(unittest.TestCase):

    def test_arbo(self):
        for folder in arbo:
            folder_path = os.path.join(STUDY_PATH, *folder)
            self.assertTrue(os.path.exists(folder_path))


if __name__ == '__main__':
    unittest.main()
