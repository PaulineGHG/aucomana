import os
from analysis_runs.analysis import *
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

R04 = "run04"
# ORG_TSV = "data/species_group.tsv"


# ### compare SR & LR #############################################################################

# compare_groups("run04", [("SR", 2), ("LR", 2), ("PUB", 2)], boxplot=True)

intersect_rnx_groups(R04, [("LR", 2), ("SR", 2), ("PUB", 2)], True)

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

# illustrate_comp("output_data/compare_groups/run04_compare_SR_LR")
