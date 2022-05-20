from analysis_runs.analysis import *

# ### FILES #######################################################################################

ORG_FILE = 'data/species_group.tsv'
PATH_STUDY = os.getcwd()
PATH_RUNS = 'data/runs/'
# PATH_RUNS = '/home/phamongi/Documents/Runs/'
# DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"
# ORG_TSV = "data/species_group.tsv"
phylo_f1 = "data/Phaeoexplorer_MLtree_rooted.nex"
phylo_f2 = "data/SpeciesTree_rooted.nex"

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

A = Analysis(PATH_RUNS, PATH_STUDY, ORG_FILE)

REACTIONS = A.get_reactions_inst(runs=[R04], out=10)
# PATHWAYS = A.get_pathways_inst(runs=[R04])
# GENES = A.get_genes_inst(runs=[R04], group=("brown", 1))
# METABOLITES = A.get_metabolites_inst(runs=[R04], group=("brown", 1))

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
# print(REACTIONS)
# print(PATHWAYS)
# PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_allrxn", phylo_f2, 1000)
# print(GENES)
# print(METABOLITES)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out10", phylo_f2, n_boot=1000)
