import scipy.stats
import analysis_runs
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


file = "data/runs/run04/analysis/all/reactions.tsv"
comp_dir_SR_LR = "output_data/compare_SR_LR/"
lr = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("LR", 2))
sr = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("SR", 2))
brown = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("brown", 1))
lr = set(lr).intersection(set(brown))
sr = set(sr).intersection(set(brown))

def get_freq_list(r: analysis_runs.reactions.Reactions):
    frequency = []
    for sp in r.species_list:
        f = sum(r.data_reactions[sp])
        if f < 5:
            frequency.append(f)
    return frequency


def get_under_mean(comp_dir):
    sp_under_mean = []
    for dir in os.listdir(comp_dir):
        decomp = dir.split(".")[0].split("_")
        if decomp[0] == "compare":
            g1 = decomp[1]
            g2 = decomp[2]

            dict_comp = pd.read_csv(f"{comp_dir}{dir}", sep="\t", index_col=0).to_dict()
            dict_g1 = pd.read_csv(f"{comp_dir}{g1}.tsv", sep="\t", index_col=0).to_dict()
            dict_g2 = pd.read_csv(f"{comp_dir}{g2}.tsv", sep="\t", index_col=0).to_dict()

            mean_g1 = dict_comp["nb_genes mean"][g1]
            mean_g2 = dict_comp["nb_genes mean"][g2]
            for sp in dict_g1["nb_genes"].keys():
                if dict_g1["nb_genes"][sp] < mean_g1:
                    sp_under_mean.append(sp)
            return sp_under_mean


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


barplot_freq(RM)
# print(scipy.stats.chisquare(get_freq_list(RM)))



