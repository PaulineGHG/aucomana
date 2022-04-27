import os

from analysis_runs.utils import *
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

R04 = "run04"
ORG_TSV = "data/species_group.tsv"


# ### compare SR & LR #############################################################################

# compare_groups(R04, ("SR", 2), ("LR", 2), ORG_TSV)

# intersect_groups(R04, ("LR", 2), ("SR", 2), ORG_TSV, True)

# Corelation nb rnx & fragmentation


def reg_lin(file_tsv):
    df = pd.read_csv(file_tsv, sep="\t", index_col="sp")
    X = list(df["#contigs"])
    Y = list(df["#rnx"])

    res = stats.linregress(x=X, y=Y)
    reg_line = [res.slope * x + res.intercept for x in X]

    axes = plt.axes()
    axes.grid()
    plt.scatter(X, Y)
    plt.plot(X, reg_line, c='r')
    plt.xlabel("#contigs")
    plt.ylabel("#rnx")
    plt.title(f"r = {res.rvalue}")
    plt.show()


# reg_lin("data/SR_frag_nbrnx.tsv")

# ### illustration ###########


def save_figur_comp(folder, df_comp, group):
    stat = (" mean", " med", " sd", " min", " max")
    to_calculate = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")
    group_file = None
    for file in os.listdir(folder):
        if file.split(".")[0] == group:
            group_file = os.path.join(folder, file)
    if group_file is not None:
        df_group = pd.read_csv(group_file, sep="\t", index_col=0)

        fig = plt.figure(figsize=[12.8, 9.6])
        i = 1

        for calc in to_calculate:
            axs = fig.add_subplot(2, 2, i)
            i += 1
            X = df_group[calc]
            axs.hist(X, len(X) - len(X)//2, density=True, facecolor='g')
            axs.set_xlabel(calc)
            axs.set_ylabel('Density')
            mu = df_comp.loc[group, calc + stat[0]]
            med = df_comp.loc[group, calc + stat[1]]
            sd = df_comp.loc[group, calc + stat[2]]
            mini = df_comp.loc[group, calc + stat[3]]
            maxi = df_comp.loc[group, calc + stat[4]]
            axs.set_title(f"Histogram of {calc} for the {group} group\n"
                                          f"$\mu={mu},\ \sigma={sd}$")
            axs.axvline(x=mu, color='b', label=f"$\mu={mu}$")
            axs.axvline(x=med, color='c', label=f"med={med}")
            axs.axvline(x=mini, color='m', label=f"min={mini}")
            axs.axvline(x=maxi, color='r', label=f"max={maxi}")
            axs.grid(True)
            axs.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        fig.tight_layout()
        path_fig = os.path.join(folder, f"{group}_hist.png")
        plt.savefig(path_fig)
        # plt.show()


def illustrate_comp(folder):
    comp_file = None
    for file in os.listdir(folder):
        if file.split("_")[0] == "compare":
            comp_file = os.path.join(folder, file)
    if comp_file is not None:
        df_comp = pd.read_csv(comp_file, sep="\t", index_col=0)
        for group in df_comp.index:
            save_figur_comp(folder, df_comp, group)


illustrate_comp("output_data/compare_groups/run04_compare_SR_LR")
