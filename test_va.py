import scipy.stats
import analysis_runs
from analysis_runs.utils import *
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


file = "data/runs/run04/analysis/all/reactions.tsv"
ORG_FILE = "data/species_group.tsv"
comp_dir_SR_LR = "output_data/compare_SR_LR/"
lr = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("LR", 2))
sr = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("SR", 2))
brown = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("brown", 1))
lr = set(lr).intersection(set(brown))
sr = set(sr).intersection(set(brown))


RM = analysis_runs.reactions.Reactions(file, species_list=brown, out=5)


def barplot_freq(rnx_obj):
    freq = []
    sp_l = []
    color_l = []
    for sp, nb in rnx_obj.nb_reactions_sp.items():
        if sp in lr:
            freq.append(round(nb / rnx_obj.nb_reactions, 4))
            sp_l.append(analysis_runs.utils.get_abr_name(sp))
            color_l.append("#889EBE")
        elif sp in sr:
            freq.append(round(nb / rnx_obj.nb_reactions, 4))
            sp_l.append(analysis_runs.utils.get_abr_name(sp))
            color_l.append("#BE88A2")

    freq, sp_l, color_l = zip(*sorted(zip(freq, sp_l, color_l)))
    fig, ax = plt.subplots(figsize=(15, 10))
    bars = ax.barh(sp_l, freq, color=color_l)
    red_patch = mpatches.Patch(color='#BE88A2', label='Short Read')
    blue_patch = mpatches.Patch(color='#889EBE', label='Long Read')
    ax.legend(handles=[red_patch, blue_patch], bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower right", ncol=2)
    ax.bar_label(bars)
    plt.show()


def boxplot_grp_completion_rnx(run, org_file, group1, group2):
    reaction_file = os.path.join(PATH_RUNS, run, "analysis", "all", "reactions.tsv")
    g1_l = get_cat_l(reaction_file, org_file, group1)
    g2_l = get_cat_l(reaction_file, org_file, group2)
    all_sp = g1_l + g2_l
    out = int(len(all_sp)*0.2)
    print(out)
    freq_g1 = []
    freq_g2 = []
    r = Reactions(reaction_file, species_list=all_sp, out=out)
    for sp, nb in r.nb_reactions_sp.items():
        if sp in g1_l:
            freq_g1.append(nb)
        elif sp in g2_l:
            freq_g2.append(nb)
    print(freq_g1, freq_g2)
    # fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_axes([0, 0, 1, 1])
    # ax.boxplot(g1_l)
    # plt.show()




# barplot_freq(RM)
boxplot_grp_completion_rnx("run04", ORG_FILE, ("LR", 2), ("SR", 2))
# print(scipy.stats.chisquare(get_freq_list(RM)))
