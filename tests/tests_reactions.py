import unittest
import pandas as pd
from analysis_runs.reactions import Reactions
from analysis_runs.utils import get_cat_l

STUDY_PATH = "Study_folder"
RUNS_PATH = "Runs_aucome"

FILE_TEST = f"{RUNS_PATH}/runA1/analysis/all/reactions.tsv"
ORG_TSV = f"{STUDY_PATH}/species_group.tsv"

brown_algae = get_cat_l(FILE_TEST, ORG_TSV, "brown")

R = Reactions(FILE_TEST)
RS = Reactions(FILE_TEST, brown_algae)
RSO1 = Reactions(FILE_TEST, brown_algae, out=1)
RSO2 = Reactions(FILE_TEST, brown_algae, out=2)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertIsInstance(R, Reactions)

    def test_init_species_default(self):
        species_obj = ['Ectocarpus_crouaniorum_m', 'Laminarionema_elsbetiae',
                       'Scytosiphon_promiscuus_MALE', 'Undaria_pinnatifida_Kr2015',
                       'Desmarestia_dudresnayi', 'Thalassiosira_pseudonana',
                       'Fragilariopsis_cylindrus', 'Phaeodactylum_tricornutum',
                       'Saccharina_japonica', 'Ectocarpus_subulatus',
                       'Porterinema_fluviatile', 'Ectocarpus_siliculosus',
                       'Cladosiphon_okamuranus', 'Ectocarpus_fasciculatus_m',
                       'Nemacystus_decipiens', 'Nannochloropsis_gaditana']
        self.assertEqual(species_obj, R.species_list)

    def test_species_precised(self):
        self.assertEqual(brown_algae, RS.species_list)

    def test_name(self):
        self.assertEqual(R.name, 'runA1')

    def test_reactions(self):
        self.assertEqual(3737, len(RS.reactions_list))
        self.assertEqual('RXN-14728', RS.reactions_list[0])
        self.assertEqual('RXN-18669', RS.reactions_list[-1])

    def test_reactions_filtered1(self):
        self.assertEqual(2859, len(RSO1.reactions_list))
        self.assertIn('RXN-14728', RSO1.reactions_list)
        self.assertIn('XYLULOKIN-RXN', RSO1.reactions_list)
        self.assertNotIn('1.1.1.285-RXN', RSO1.reactions_list)

    def test_reactions_filtered2(self):
        self.assertEqual(3053, len(RSO2.reactions_list))
        self.assertIn('RXN-14728', RSO2.reactions_list)
        self.assertIn('1.1.1.285-RXN', RSO2.reactions_list)
        self.assertNotIn('1.1.4.2-RXN', RSO1.reactions_list)

    def test_data_reactions(self):
        self.assertIsInstance(RSO1.data_reactions, pd.DataFrame)
        self.assertListEqual(RSO1.reactions_list, list(RSO1.data_reactions.index))
        self.assertListEqual(RSO1.species_list, list(RSO1.data_reactions.columns))

    def test_data_genes_assoc(self):
        self.assertIsInstance(RSO1.data_genes_assoc, pd.DataFrame)
        self.assertListEqual(RSO1.reactions_list, list(RSO1.data_genes_assoc.index))
        index_list = [x + RSO1.STR_GENE_ASSOC for x in RSO1.species_list]
        self.assertListEqual(index_list, list(RSO1.data_genes_assoc.columns))


if __name__ == '__main__':
    unittest.main()


