import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib_venn

from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.genes import Genes
from analysis_runs.init_analysis import PATH_RUNS
from typing import Tuple


def get_cat_l(reactions_file, organisms_file, cat: Tuple[str, int]):
    with open(organisms_file, "r") as org_f, open(reactions_file, "r") as rea_f:
        species_l = set()
        cat_l = []
        for l in rea_f:
            l = l.split("\t")
            for x in l[1:]:
                if x[-7:] == '(sep=;)':
                    break
                species_l.add(x)
            break
        for l in org_f:
            l = l.split()
            if l[cat[1]] == cat[0] and l[0] in species_l:
                cat_l.append(l[0])
        return cat_l


def write_cut_reactions_file(original_file, cut_nb, reac_list):
    name = original_file.split("/")[-1]
    with open(original_file, "r") as f, open(f"outputs/cut_reactions_data/cut{cut_nb}_{name}", "w") as o:
        for line in f:
            l = line.split("\t")
            if l[0] in reac_list or l[0] == "reaction":
                o.write(line)


def get_reactions_inst(path_runs, org_tsv, cat=None, out=None):
    r_dic = {}
    for run in os.listdir(path_runs):
        r_path = os.path.join(path_runs, run, "analysis", "all", "reactions.tsv")
        if os.path.exists(r_path):
            if cat is not None:
                species_l = get_cat_l(r_path, org_tsv, cat)
                r_dic[run] = Reactions(r_path, species_l, out)
            else:
                r_dic[run] = Reactions(r_path, cat, out)
    return r_dic


def get_pathways_inst(path_runs, org_tsv, cat=None, out=None):
    p_dic = {}
    for run in os.listdir(path_runs):
        r_path = os.path.join(path_runs, run, "analysis", "all", "reactions.tsv")
        p_path = os.path.join(path_runs, run, "analysis", "all", "pathways.tsv")
        if os.path.exists(r_path) and os.path.exists(p_path):
            if cat is not None:
                species_l = get_cat_l(r_path, org_tsv, cat)
                p_dic[run] = Pathways(p_path, species_l, out)
            else:
                p_dic[run] = Pathways(p_path, cat, out)
    return p_dic


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

    out_dir = f"output_data/compare_{group1[0]}_{group2[0]}/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    g1_res_file = f"{out_dir}{group1[0]}.tsv"
    g2_res_file = f"{out_dir}{group2[0]}.tsv"
    res_file = f"{out_dir}compare_{group1[0]}_{group2[0]}.tsv"

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
