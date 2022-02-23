import unittest
from reactions_loss import Reactions

FILE_TEST = "data/runA1_reactions.tsv"
R = Reactions(FILE_TEST)


class Test(unittest.TestCase):

    def test_init(self):
        self.assertEqual(type(R), Reactions)

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
        species_list = ['Ectocarpus_crouaniorum_m', 'Laminarionema_elsbetiae',
                       'Scytosiphon_promiscuus_MALE', 'Undaria_pinnatifida_Kr2015',
                       'Desmarestia_dudresnayi', 'Thalassiosira_pseudonana']
        RS = Reactions(FILE_TEST, species_list)
        self.assertEqual(species_list, RS.species_list)

    def test_name(self):
        self.assertEqual(R.name, 'runA1_reactions')


if __name__ == '__main__':
    unittest.main()


