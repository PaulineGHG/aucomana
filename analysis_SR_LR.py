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


def save_figur_comp(folder, df_comp, group, calc):
    stat = (" mean", " med", " sd", " min", " max")
    group_file = None
    for file in os.listdir(folder):
        if file.split(".")[0] == group:
            group_file = os.path.join(folder, file)
    if group_file is not None:
        df_group = pd.read_csv(group_file, sep="\t", index_col=0)

        X = df_group[calc]

        plt.hist(X, len(X) - len(X)//2, density=True, facecolor='g')
        plt.xlabel(calc)
        plt.ylabel('Density')
        mu = df_comp.loc[group, calc + stat[0]]
        med = df_comp.loc[group, calc + stat[1]]
        sd = df_comp.loc[group, calc + stat[2]]
        mini = df_comp.loc[group, calc + stat[3]]
        maxi = df_comp.loc[group, calc + stat[4]]
        plt.title(f"Histogram of {calc} for the {group} group\n"
                  f"$\mu={mu},\ \sigma={sd}$")
        plt.axvline(x=mu, color='b', label=f"$\mu={mu}$")
        plt.axvline(x=med, color='c', label=f"med={med}")
        plt.axvline(x=mini, color='m', label=f"min={mini}")
        plt.axvline(x=maxi, color='r', label=f"max={maxi}")
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        # plt.savefig("output.png", bbox_inches="tight")
        plt.show()


def illustrate_comp(folder):
    to_calculate = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")
    comp_file = None
    for file in os.listdir(folder):
        if file.split("_")[0] == "compare":
            comp_file = os.path.join(folder, file)
    if comp_file is not None:
        df_comp = pd.read_csv(comp_file, sep="\t", index_col=0)
        for group in df_comp.index:
            for calc in to_calculate:
                save_figur_comp(folder, df_comp, group, calc)


illustrate_comp("output_data/compare_SR_LR/")
