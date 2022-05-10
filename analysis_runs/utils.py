import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn

from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.genes import Genes
from analysis_runs.init_analysis import PATH_RUNS
from typing import Tuple, List, Dict


def get_grp_l(run: str, organisms_file: str, group: Tuple[str, int]) \
        -> List[str]:
    """ Select species according to the group they belong to. The groups must be specified in a
    TSV file (organisms-file as input).

    Parameters
    ----------
    run: str
        ID of the run
    organisms_file: str
        File providing groups information
    group: Tuple[str, int]
        The group to consider : Tuple(name of the group, column of the group in the organisms_file)

    Returns
    -------
    group_l: List[str]
        List of species corresponding to the category chosen
    """
    grp_template_f = os.path.join(PATH_RUNS, run, "analysis",
                                  "group_template.tsv")
    with open(organisms_file, "r") as org_f, open(grp_template_f, "r") as grp_f:
        species_l = set()
        group_list = []
        for line in grp_f:
            line = line.split("\t")
            if line[0] == "all":
                for species in line[1:]:
                    species_l.add(species)
                break
        for line in org_f:
            line = line.split()
            if line[group[1]] == group[0] and line[0] in species_l:
                group_list.append(line[0])
        return group_list


# def write_cut_reactions_file(original_file, cut_nb, reac_list):
#     name = original_file.split("/")[-1]
#     with open(original_file, "r") as f, open(f"outputs/cut_reactions_data/"
#                                              f"cut{cut_nb}_{name}", "w") as o:
#         for line in f:
#             l = line.split("\t")
#             if l[0] in reac_list or l[0] == "reaction":
#                 o.write(line)


def get_reactions_inst(organisms_file: str = None, runs: List[str] = None,
                       group: Tuple[str, int] = None, out: int = None) -> Dict[str, 'Reactions']:
    """ Create Reactions instances in a dictionary.

    Parameters
    ----------
    organisms_file: str, optional (default=None)
        File providing groups information
    runs: List[str], optional (default=None)
        List of the runs ID to consider
        If None, will be the list of all runs in the Runs path
    group: Tuple[str, int], optional (default=None)
        The group to consider : Tuple(name of the group, column of the
        group in the organisms_file)
        If None will be all the species of each run
    out: int, optional (default=None)
        Number of species maximum not having the reaction for the reaction to be kept
        If None, will not filter the reactions

    Returns
    -------
    reactions_dic: Dict[str, 'Reactions']
        Dictionary of Reactions instances : Dict[ID of the run, Reactions instance of the run]
    """
    rnx_dic = {}
    if runs is None:
        runs = os.listdir(PATH_RUNS)
    for run in runs:
        r_path = os.path.join(PATH_RUNS, run, "analysis", "all",
                              "reactions.tsv")
        if os.path.exists(r_path):
            if group is not None:
                if organisms_file is None:
                    print("If group is specified, organisms_file must also be specified for the "
                          "filter to be applied. Here no organisms_file has been specified so all "
                          "species has been kept.")
                    group = None
                    rnx_dic[run] = Reactions(r_path, out)
                else:
                    species_l = get_grp_l(run, organisms_file, group)
                    rnx_dic[run] = Reactions(r_path, species_l, out)
            else:
                rnx_dic[run] = Reactions(r_path, group, out)
    print(f"Reactions instances has been created for runs : {runs} with group = {group} and out = "
          f"{out}")
    return rnx_dic


def get_pathways_inst(path_runs, org_tsv, runs=None, cat=None, out=None, nb_rnx_px_min=0):
    p_dic = {}
    if runs is None:
        runs = os.listdir(path_runs)
    for run in runs:
        r_path = os.path.join(path_runs, run, "analysis", "all", "reactions.tsv")
        p_path = os.path.join(path_runs, run, "analysis", "all", "pathways.tsv")
        if os.path.exists(r_path) and os.path.exists(p_path):
            if cat is not None:
                species_l = get_grp_l(r_path, org_tsv, cat)
                p_dic[run] = Pathways(p_path, species_l, out, nb_rnx_px_min)
            else:
                p_dic[run] = Pathways(p_path, cat, out, nb_rnx_px_min)
    return p_dic


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


def compare_groups(run, group1, group2, org_file, hist=True):
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

    out_dir = os.path.join("output_data", "compare_groups", f"{run}_compare_{group1[0]}_{group2[0]}")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    g1_res_file = os.path.join(out_dir, f"{group1[0]}.tsv")
    g2_res_file = os.path.join(out_dir, f"{group2[0]}.tsv")
    res_file = os.path.join(out_dir, f"compare_{group1[0]}_{group2[0]}.tsv")

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
            df_comp.loc[grp[0], calc + stat[0]] = round(float(np.mean(grp[1][calc])), 1)
            df_comp.loc[grp[0], calc + stat[1]] = np.median(grp[1][calc])
            df_comp.loc[grp[0], calc + stat[2]] = round(np.sqrt(np.var(grp[1][calc])), 1)
            df_comp.loc[grp[0], calc + stat[3]] = min(grp[1][calc])
            df_comp.loc[grp[0], calc + stat[4]] = max(grp[1][calc])

    df_g1.to_csv(g1_res_file, sep="\t")
    df_g2.to_csv(g2_res_file, sep="\t")
    df_comp.to_csv(res_file, sep="\t")

    if hist:
        illustrate_comp(out_dir)


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

    print(f"Over {union_nb} reactions, {intersect_nb} reactions in common = "
          f"{round((intersect_nb/union_nb) * 100, 2)} % of the union.\n"
          f"{g1_nb - intersect_nb} only present among the {group1[0]} group = "
          f"{round(((g1_nb - intersect_nb)/union_nb) * 100, 2)} % of the union.\n"
          f"{g2_nb - intersect_nb} only present among the {group2[0]} group = "
          f"{round(((g2_nb - intersect_nb)/union_nb) * 100, 2)} % of the union.\n")

    if venn_plot:
        matplotlib_venn.venn2([reactions_g1, reactions_g2], (group1[0], group2[0]))
        plt.show()


def get_abr_name(name):
    name = name.split("_")
    return f"{name[0][0]}.{name[1][:4]}"
