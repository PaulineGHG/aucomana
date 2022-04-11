from analysis_runs.utils import *
from analysis_runs.init_analysis import PATH_RUNS
import pandas as pd

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


SHORT_READS = get_cat_l("data/runs/run04/analysis/all/reactions.tsv", ORG_TSV, ("SR", 2))
LONG_READS = get_cat_l("data/runs/run04/analysis/all/reactions.tsv", ORG_TSV, ("LR", 2))


def compare_groups(run, group1, group2):
    to_calculate = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")
    stat = (" mean", " med", " sd", " min", " max")
    path = f"{PATH_RUNS}/{run}/analysis/all/"
    R_G1 = Reactions(f"{path}reactions.tsv", group1)
    R_G2 = Reactions(f"{path}reactions.tsv", group2)
    P_G1 = Pathways(f"{path}pathways.tsv", group1)
    P_G2 = Pathways(f"{path}pathways.tsv", group2)
    G_G1 = Genes(f"{path}genes.tsv", group1)
    G_G2 = Genes(f"{path}genes.tsv", group2)

    g1_res_file = "output_data/group1.tsv"
    g2_res_file = "output_data/group2.tsv"
    res_file = "output_data/compare_groups.tsv"

    df_g1 = pd.DataFrame(columns=to_calculate, index=group1)
    df_g2 = pd.DataFrame(columns=to_calculate, index=group2)
    col_stat = []
    for x in to_calculate:
        for y in stat:
            col_stat.append(x + y)
    df_comp = pd.DataFrame(columns=col_stat, index=["group 1", "group 2"])

    for sp in group1:
        df_g1.loc[sp, to_calculate[0]] = G_G1.nb_genes_species[sp]
        df_g1.loc[sp, to_calculate[1]] = R_G1.nb_reactions_sp[sp]
        df_g1.loc[sp, to_calculate[2]] = P_G1.get_pw_over_treshold(sp, 0.8)[sp][0]

    print(df_g1)



compare_groups(R04, SHORT_READS, LONG_READS)
