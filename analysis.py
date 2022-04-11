from analysis_runs.utils import *
from analysis_runs.init_analysis import PATH_RUNS

# create_folders(PATH_STUDY)

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


# with open("pathways_slat_plac.txt", "w") as f:
#     f.write("1) Pathways absents chez Slat et présents chez les autres : ")
#     p_l = PATHWAYS[R01].get_pw_absent(SLAT)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#     f.write("\nPathways absents chez Plac et présents chez les autres : ")
#     p_l = PATHWAYS[R01].get_pw_absent(PLAC)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#
#     f.write("\n\n2) Pathways de completion minimale chez Slat : ")
#     p_l = PATHWAYS[R01].get_pw_min(SLAT)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#     f.write("\nPathways de completion minimale chez Plac : ")
#     p_l = PATHWAYS[R01].get_pw_min(PLAC)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#
#     f.write("\n\n3) Pathways incomplets chez Slat et complets chez autres espèces : ")
#     p_l = PATHWAYS[R01].get_pw_incomplete(SLAT)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#     f.write("\nPathways incomplets chez Plac et complets chez autres espèces : ")
#     p_l = PATHWAYS[R01].get_pw_incomplete(PLAC)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")

SHORT_READS = get_cat_l("data/runs/run04/analysis/all/reactions.tsv", ORG_TSV, ("SR", 2))
LONG_READS = get_cat_l("data/runs/run04/analysis/all/reactions.tsv", ORG_TSV, ("LR", 2))
print(SHORT_READS)
print(LONG_READS)

R04_SR = Reactions("data/runs/run04/analysis/all/reactions.tsv", species_list=SHORT_READS)
R04_LR = Reactions("data/runs/run04/analysis/all/reactions.tsv", species_list=LONG_READS)

# P04_SR = Pathways("data/runs/run04/analysis/all/pathways.tsv", species_list=SHORT_READS)
# P04_LR = Pathways("data/runs/run04/analysis/all/pathways.tsv", species_list=LONG_READS)

G04_SR = Genes("data/runs/run04/analysis/all/genes.tsv", species_list=SHORT_READS)
print(G04_SR.nb_genes_species)


