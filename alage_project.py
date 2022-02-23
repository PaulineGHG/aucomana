from reactions_loss import Reactions
from pandas_ods_reader import read_ods


def reactions_from_file(file):
    reactions = []
    d = read_ods(file)
    for i in range(d.shape[0]):
        r = d.loc[i, "Identifiant Metacyc"]
        reactions.append(r)
    return reactions


# ###Description ###

# 40 : run 40+7
# 01 : run Pauline all species
# A0 : run Alexandre
# A1 : run Pauline, Alexandre like exactly
# A2 : run Pauline, Alexandre like improved

# ### FILES ####

DATA_FILE_01 = "data/run01_reactions.tsv"
DATA_FILE_40 = 'data/run40_reactions.tsv'
DATA_FILE_A0 = 'data/runA0_reactions.tsv'
DATA_FILE_A1 = 'data/runA1_reactions.tsv'
DATA_FILE_A2 = 'data/runA2_reactions.tsv'
DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"

# ### Select species ###

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

# ### Class instances ###

R01 = Reactions(DATA_FILE_01)
R40 = Reactions(DATA_FILE_40, BROWN_ALGAE_40)
RA0 = Reactions(DATA_FILE_A0)
RA1 = Reactions(DATA_FILE_A1)
RA2 = Reactions(DATA_FILE_A2)

reac_lostA = reactions_from_file(DATA_LELSB_LOSSES)


# ### Laminarionema loss ###

print(R40.reactions_loss[LAMINARIONEMA_E])
print(R01.reactions_loss[LAMINARIONEMA_E])
print(RA0.reactions_loss[LAMINARIONEMA_E])
print(RA1.reactions_loss[LAMINARIONEMA_E])
print(RA2.reactions_loss[LAMINARIONEMA_E])


# ### Common Reactions ###

print(Reactions.get_common_reactions([R01, R40, RA0], LAMINARIONEMA_E, output_file=True))
print(Reactions.get_common_reactions([R01, R40, RA1], LAMINARIONEMA_E, output_file=True))
print(Reactions.get_common_reactions([R01, R40, RA2], LAMINARIONEMA_E, output_file=True))


