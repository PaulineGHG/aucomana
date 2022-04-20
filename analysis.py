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


def intersect_groups(run, group1, group2, org_file, venn_plot=False):
    file = f"{PATH_RUNS}{run}/analysis/all/reactions.tsv"
    group1_list = get_cat_l(file, org_file, group1)
    group2_list = get_cat_l(file, org_file, group2)

    r1 = Reactions(file, group1_list)
    r2 = Reactions(file, group2_list)

    reactions_g1 = set(r1.reactions_list)
    reactions_g2 = set(r2.reactions_list)

    g1_nb = len(reactions_g1)
    g2_nb = len(reactions_g2)
    intersect_nb = len(reactions_g1.intersection(reactions_g2))
    union_nb = len(reactions_g1.union(reactions_g2))

    print(f"Over {union_nb} reactions, {intersect_nb} reactions in common.\n"
          f"{g1_nb - intersect_nb} only present among the {group1[0]} group = "
          f"{round(((g1_nb - intersect_nb)/union_nb) * 100, 2)} %.\n"
          f"{g2_nb - intersect_nb} only present among the {group2[0]} group = "
          f"{round(((g2_nb - intersect_nb)/union_nb) * 100, 2)} %.\n")

    if venn_plot:
        mven.venn2([reactions_g1, reactions_g2], (group1[0], group2[0]))
        plt.show()


intersect_groups(R04, ("LR", 2), ("SR", 2), ORG_TSV, True)
