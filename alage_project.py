from reactions_loss import Reactions
from pandas_ods_reader import read_ods


def reactions_from_file(file):
    reactions = []
    d = read_ods(file)
    for i in range(d.shape[0]):
        r = d.loc[i, "Identifiant Metacyc"]
        reactions.append(r)
    return reactions


DATA_FILE_1 = "run1_reactions.tsv"

DATA_FILE_40 = 'run40_reactions.tsv'

DATA_FILE_A = 'reactions_phaeo2.tsv'

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

LAMINARIONEMA_E = 'Laminarionema_elsbetiae'

DATA_LELSB_LOSSES = "Lelsb_losses.ods"

R1 = Reactions(DATA_FILE_1)
R40 = Reactions(DATA_FILE_40, BROWN_ALGAE_40)
reac_lostA = reactions_from_file(DATA_LELSB_LOSSES)
RA = Reactions(DATA_FILE_A)

R40_loss_lami_e = R40.reactions_loss[LAMINARIONEMA_E]
R1_loss_lami_e = R1.reactions_loss[LAMINARIONEMA_E]
print(R40_loss_lami_e)
print(R1_loss_lami_e)
print(RA.reactions_loss[LAMINARIONEMA_E])


# reac_comR1_R40 = R1.get_common_reactions([R40], LAMINARIONEMA_E)
#
# reac_comR1_R40_RA = []
# for r in reac_comR1_R40[1]:
#     if r in reac_lostA:
#         reac_comR1_R40_RA.append(r)
# print(len(reac_comR1_R40_RA), reac_comR1_R40_RA)

print(R1.get_common_reactions([R40, RA], LAMINARIONEMA_E))

# R40.print_genes_assoc(R40.get_genes_assoc(r_interest40))
# R1.print_genes_assoc(R1.get_genes_assoc(r_interest1))

