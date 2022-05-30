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

# REACTIONS = A.get_reactions_inst(runs=[R04], out=34)
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
# REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out10", phylo_f2, n_boot=1000)

print("REACTIONS : 1")
REACTIONS = A.get_reactions_inst(runs=[R04], out=34)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out05perc", phylo_f2, n_boot=10000)

print("REACTIONS : 2")
REACTIONS = A.get_reactions_inst(runs=[R04], out=32)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out10perc", phylo_f2, n_boot=10000)

print("REACTIONS : 3")
REACTIONS = A.get_reactions_inst(runs=[R04], out=28)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out20perc", phylo_f2, n_boot=10000)

print("PATHWAYS ALL RXN : 1")
PATHWAYS = A.get_pathways_inst(runs=[R04])
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_allrxn", phylo_f2, 10000)

PATHWAYS = A.get_pathways_inst(runs=[R04], nb_rnx_px_min=2)
print("PATHWAYS 2 MIN : 1")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_2min", phylo_f2, 10000)

PATHWAYS = A.get_pathways_inst(runs=[R04], nb_rnx_px_min=3)
print("PATHWAYS 3 MIN : 1")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_3min", phylo_f2, 10000)

PATHWAYS = A.get_pathways_inst(runs=[R04], nb_rnx_px_min=4)
print("PATHWAYS 4 MIN : 1")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_4min", phylo_f2, 10000)
