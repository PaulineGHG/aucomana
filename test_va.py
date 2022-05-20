import scipy.stats
import analysis_runs
from analysis_runs.analysis import Analysis
from analysis_runs.reactions import Reactions
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


run = "run04"
ORG_FILE = 'data/species_group.tsv'
PATH_STUDY = os.getcwd()
PATH_RUNS = 'data/runs/'
comp_dir_SR_LR = "output_data/compare_SR_LR/"
A = Analysis(PATH_RUNS, PATH_STUDY, ORG_FILE)


RM = A.get_reactions_inst(runs=[run], group=("brown", 1))[run]


def barplot_freq(rnx_obj):
    freq = []
    sp_l = []
    color_l = []
    for sp, nb in rnx_obj.nb_reactions_sp.items():
        if sp in lr:
            freq.append(round(nb / rnx_obj.nb_reactions, 4))
            sp_l.append(analysis_runs.analysis.get_abbr_name(sp))
            color_l.append("#889EBE")
        elif sp in sr:
            freq.append(round(nb / rnx_obj.nb_reactions, 4))
            sp_l.append(analysis_runs.analysis.get_abbr_name(sp))
            color_l.append("#BE88A2")

    freq, sp_l, color_l = zip(*sorted(zip(freq, sp_l, color_l)))
    fig, ax = plt.subplots(figsize=(15, 10))
    bars = ax.barh(sp_l, freq, color=color_l)
    red_patch = mpatches.Patch(color='#BE88A2', label='Short Read')
    blue_patch = mpatches.Patch(color='#889EBE', label='Long Read')
    ax.legend(handles=[red_patch, blue_patch], bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower right", ncol=2)
    ax.bar_label(bars)
    plt.show()


def boxplot_grp_completion_rnx(run, group1, group2, a):
    reaction_file = os.path.join(PATH_RUNS, run, "analysis", "all", "reactions.tsv")
    g1_l = a.get_grp_l(run, group1)
    g2_l = a.get_grp_l(run, group2)
    all_sp = g1_l + g2_l
    out = int(len(all_sp)*0.2)
    freq_g1 = []
    freq_g2 = []
    r = a.get_reactions_inst([run], species_list=all_sp, out=out)[run]
    for sp, nb in r.nb_reactions_sp.items():
        if sp in g1_l:
            freq_g1.append(nb)
        elif sp in g2_l:
            freq_g2.append(nb)
    print(min(freq_g1))
    print(np.median(freq_g2))
    plt.boxplot([freq_g1, freq_g2], labels=[group1[0], group2[0]])
    plt.ylabel("Nb reactions")
    plt.title(f"Boxplot of number of reactions among a conserved subset\n"
              f"of reactions for {group1[0]} and {group2[0]} groups.")
    plt.show()


# barplot_freq(RM)
boxplot_grp_completion_rnx(run, ("LR", 2), ("SR", 2), A)
# print(scipy.stats.chisquare(get_freq_list(RM)))
