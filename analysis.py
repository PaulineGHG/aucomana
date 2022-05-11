from analysis_runs.utils import *
from analysis_runs.init_analysis import *

# ### FILES #######################################################################################

# DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"
# ORG_TSV = "data/species_group.tsv"

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

R01 = "run01"
R02 = "run02"
R03 = "run03"
R04 = "run04"
R40 = "run40"

REACTIONS = get_reactions_inst(runs=[R04])
PATHWAYS = get_pathways_inst(runs=[R04], group=("brown", 1), nb_rnx_px_min=3)

# ### Common Reactions ############################################################################

# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True)
# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True, union=True)

# ### genes assoc #################################################################################

# REACTIONS[R01].get_genes_assoc(LELS, {"12-OXOPHYTODIENOATE-REDUCTASE-RXN"}, output_file=True)
# print(REACTIONS[R04].get_genes_assoc({"12-OXOPHYTODIENOATE-REDUCTASE-RXN",
#                                       "LEUKOTRIENE-C4-SYNTHASE-RXN",
#                                       "PROSTAGLANDIN-E-SYNTHASE-RXN"}, output_file=True))


# ### pathways ##################################################################################

# print(PATHWAYS[R01].get_pw_complete([SLAT, PLAC], unique=True))
# PATHWAYS[R04].convert_df_to_binary(1, output_file=True)
# print(PATHWAYS[R04].get_pw_names())
print(REACTIONS)
print(PATHWAYS)

