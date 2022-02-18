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
# 0A : run Alexandre
# 02 : run Pauline, Alexandre like exactly
# 03 : run Pauline, Alexandre like improved

# ### FILES ####

DATA_FILE_01 = "data/run1_reactions.tsv"
DATA_FILE_40 = 'data/run40_reactions.tsv'
DATA_FILE_0A = 'data/reactions_phaeo2.tsv'
DATA_FILE_02 = 'data/run2_reactions.tsv'
DATA_FILE_03 = 'data/run3_reactions.tsv'
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
R0A = Reactions(DATA_FILE_0A)
R02 = Reactions(DATA_FILE_02)
R03 = Reactions(DATA_FILE_03)

reac_lostA = reactions_from_file(DATA_LELSB_LOSSES)


# ### Laminarionema loss ###

print(R40.reactions_loss[LAMINARIONEMA_E])
print(R01.reactions_loss[LAMINARIONEMA_E])
print(R0A.reactions_loss[LAMINARIONEMA_E])
print(R02.reactions_loss[LAMINARIONEMA_E])
print(R03.reactions_loss[LAMINARIONEMA_E])


# ### Common Reactions ###

print(Reactions.get_common_reactions([R01, R40, R0A], LAMINARIONEMA_E, output_file=True))
print(Reactions.get_common_reactions([R01, R40, R02], LAMINARIONEMA_E, output_file=True))
print(Reactions.get_common_reactions([R01, R40, R03], LAMINARIONEMA_E, output_file=True))


