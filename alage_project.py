from reactions_loss import Reactions

DATA_FILE_1 = "run1_reactions.tsv"

DATA_FILE_40 = 'run40_reactions.tsv'
BROWN_ALGAE_40 = ['Thalassiosira_pseudonana',
                  'Fragilariopsis_cylindrus',
                  'Phaeodactylum_tricornutum',
                  'Nannochloropsis_gaditana',
                  'Ectocarpus_siliculosus',
                  'Ectocarpus_crouaniorum',
                  'Ectocarpus_subulatus',
                  'Ectocarpus_fasciculatus',
                  'Scytosiphon_lomentaria',
                  'Porterinema_fluviatile',
                  'Nemacystus_decipiens',
                  'Cladosiphon_okamuranus',
                  'Laminarionema_elsbetiae',
                  'Saccharina_japonica',
                  'Undaria_pinnatifida']

INTEREST_SPECIES = 'Laminarionema_elsbetiae'

# R1 = Reactions(DATA_FILE_1)
R40 = Reactions(DATA_FILE_40, BROWN_ALGAE_40)


R40_loss_lami_e = R40.reactions_loss[INTEREST_SPECIES]
# R1_loss_lami_e = R1.reactions_loss[INTEREST_SPECIES]
print(R40_loss_lami_e)
# print(R1_loss_lami_e)
# print(R1.get_common_reactions(R40, INTEREST_SPECIES))
# print(R40.reactions_loss)
# print(R40.data_genes_assoc)
# R40.print_genes_assoc(R40.get_genes_assoc(R40_loss_lami_e[1]))
