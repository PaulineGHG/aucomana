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


def compare_groups(run, group1, group2, org_file):
    to_calculate = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")
    stat = (" mean", " med", " sd", " min", " max")
    path = f"{PATH_RUNS}/{run}/analysis/all/"

    group1_list = get_cat_l(f"{path}reactions.tsv", org_file, group1)
    group2_list = get_cat_l(f"{path}reactions.tsv", org_file, group2)

    r_g1 = Reactions(f"{path}reactions.tsv", group1_list)
    r_g2 = Reactions(f"{path}reactions.tsv", group2_list)
    p_g1 = Pathways(f"{path}pathways.tsv", group1_list)
    p_g2 = Pathways(f"{path}pathways.tsv", group2_list)
    g_g1 = Genes(f"{path}genes.tsv", group1_list)
    g_g2 = Genes(f"{path}genes.tsv", group2_list)

    g1_res_file = "output_data/group1.tsv"
    g2_res_file = "output_data/group2.tsv"
    res_file = "output_data/compare_groups.tsv"

    df_g1 = pd.DataFrame(columns=to_calculate, index=group1_list)
    df_g2 = pd.DataFrame(columns=to_calculate, index=group2_list)
    col_stat = []
    for x in to_calculate:
        for y in stat:
            col_stat.append(x + y)
    df_comp = pd.DataFrame(columns=col_stat, index=[group1[0], group2[0]])

    for sp in group1_list:
        df_g1.loc[sp, to_calculate[0]] = g_g1.nb_genes_species[sp]
        df_g1.loc[sp, to_calculate[1]] = r_g1.nb_reactions_sp[sp]
        df_g1.loc[sp, to_calculate[2]] = p_g1.get_pw_over_treshold(sp, 0.8)[sp][0]
        df_g1.loc[sp, to_calculate[3]] = p_g1.get_pw_complete(sp, False)[sp][0]

    for sp in group2_list:
        df_g2.loc[sp, to_calculate[0]] = g_g2.nb_genes_species[sp]
        df_g2.loc[sp, to_calculate[1]] = r_g2.nb_reactions_sp[sp]
        df_g2.loc[sp, to_calculate[2]] = p_g2.get_pw_over_treshold(sp, 0.8)[sp][0]
        df_g2.loc[sp, to_calculate[3]] = p_g2.get_pw_complete(sp, False)[sp][0]

    for grp in [(group1[0], df_g1), (group2[0], df_g2)]:
        for calc in to_calculate:
            df_comp.loc[grp[0], calc + stat[0]] = None

    df_g1.to_csv(g1_res_file, sep="\t", index=False)
    df_g2.to_csv(g2_res_file, sep="\t", index=False)
    print(df_comp)


compare_groups(R04, ("SR", 2), ("LR", 2), ORG_TSV)
