import scipy.stats
import analysis_runs
import pandas as pd
import os
import matplotlib.pyplot as plt


file = "data/runs/run04/analysis/all/reactions.tsv"
comp_dir_SR_LR = "output_data/compare_SR_LR/"
brown = analysis_runs.utils.get_cat_l(file, "data/species_group.tsv", ("brown", 1))


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


# RM = analysis_runs.reactions.Reactions(file, get_under_mean(comp_dir_SR_LR))
# print(scipy.stats.chisquare(get_freq_list(RM)))



