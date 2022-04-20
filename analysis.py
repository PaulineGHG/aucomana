from analysis_runs.utils import *
import matplotlib.pyplot as plt
import matplotlib_venn as mven

# ### FILES #######################################################################################

# DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"
ORG_TSV = "data/species_group.tsv"

# ### Select species ##############################################################################

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

LELS = 'Laminarionema_elsbetiae'
SLAT = "Saccharina_latissima_FEMALE"
PLAC = "Pleurocladia_lacustris"

# ### Class instances #############################################################################

# REACTIONS = get_reactions_inst(PATH_RUNS, ORG_TSV, ("brown", 1), 1)

# PATHWAYS = get_pathways_inst(PATH_RUNS, ORG_TSV, ("brown", 1))

R01 = "run01"
R02 = "run02"
R03 = "run03"
R04 = "run04"
R40 = "run40"

# ### Common Reactions ############################################################################

# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True)
# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True, union=True)

# ### genes assoc #################################################################################

# REACTIONS[R01].get_genes_assoc(LELS, {"12-OXOPHYTODIENOATE-REDUCTASE-RXN"}, output_file=True)

# print(PATHWAYS[R01].get_pw_complete([SLAT, PLAC], unique=True))

# compare_groups(R04, ("SR", 2), ("LR", 2), ORG_TSV)

file = f"{PATH_RUNS}{R04}/analysis/all/reactions.tsv"
sr = get_cat_l(file, ORG_TSV, ("SR", 2))
lr = get_cat_l(file, ORG_TSV, ("LR", 2))

rsr = Reactions(file, sr)
rlr = Reactions(file, lr)

reactions_sr = set(rsr.reactions_list)
reactions_lr = set(rlr.reactions_list)

print(len(reactions_sr), len(reactions_lr))
# mven.venn2([reactions_sr, reactions_lr], ("Short Read", "Long Read"))
# plt.show()
