import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import venn

from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.genes import Genes
from analysis_runs.init_analysis import PATH_RUNS, PATH_STUDY, ORG_FILE
from typing import Tuple, List, Dict


# ## Various utils fnc

def get_abbr_name(name: str) -> str:
    """Abbreviate and returns the name of the species in format : X.xxxx.
    Ex : Escherichia coli returns E.coli

    Parameters
    ----------
    name: str
        name of the species

    Returns
    -------
    abbr_name: str
        Abbreviated name of the species
    """
    name = name.split("_")
    return f"{name[0][0]}.{name[1][:4]}"


# ## Extract species according to groups from organisms file

def get_grp_l(run: str, group: Tuple[str, int],
              species_list: List[str] = None) -> List[str]:
    """ Select species according to the group they belong to. The groups must be specified in a
    TSV file (organisms-file as input).

    Parameters
    ----------
    run: str
        ID of the run
    group: Tuple[str, int]
        The group to consider : Tuple(name of the group, column of the group in the organisms_file)
    species_list: List[str], optional (default=None)
        List of species to consider to be filtered according to their group belonging
        If None will be all the species of the run

    Returns
    -------
    group_l: List[str]
        List of species corresponding to the category chosen
    """
    if species_list is not None:
        species_l = set(species_list)
    else:
        grp_template_f = os.path.join(PATH_RUNS, run, "analysis",
                                      "group_template.tsv")
        with open(grp_template_f, "r") as grp_f:
            species_l = set()
            for line in grp_f:
                line = line.split("\t")
                if line[0] == "all":
                    for species in line[1:]:
                        species_l.add(species.strip())
    group_list = []
    df = pd.read_csv(ORG_FILE, sep="\t", index_col=0, header=None)
    if group[1] not in list(df.columns):
        raise ValueError(f"No column {group[1]} in {ORG_FILE} file. (Number of columns = "
                         f"{df.shape[1]})")
    if group[0] not in list(df[group[1]]):
        warnings.warn(f"No species of group {group[0]} are in the column {group[1]} of the "
                      f"{ORG_FILE}. All species has been kept.")
        return list(species_l)
    else:
        for sp in species_l:
            if df.loc[sp, group[1]] == group[0]:
                group_list.append(sp)
    return group_list


# ## Create instances

def get_reactions_inst(runs: List[str] = None, species_list: List[str] = None,
                       group: Tuple[str, int] = None, out: int = None) -> Dict[str, 'Reactions']:
    """ Create Reactions instances in a dictionary.

    Parameters
    ----------
    runs: List[str], optional (default=None)
        List of the runs ID to consider
        If None, will be the list of all runs in the Runs path
    species_list: List[str], optional (default=None)
        List of species to consider for instances creation
    group: Tuple[str, int], optional (default=None)
        The group to consider for species filtering : Tuple(name of the group, column of the
        group in the organisms_file)
        If None will be all the species of each run
    out: int, optional (default=None)
        Number of species maximum not having the reaction for the reaction to be kept
        If None, will not filter the reactions

    Returns
    -------
    reactions_dict: Dict[str, 'Reactions']
        Dictionary of Reactions instances : Dict[ID of the run, Reactions instance of the run]
    """
    reactions_dict = {}
    if runs is None:
        runs = os.listdir(PATH_RUNS)
    for run in runs:
        r_path = os.path.join(PATH_RUNS, run, "analysis", "all",
                              "reactions.tsv")
        if os.path.exists(r_path):
            if group is not None:
                grp_species_l = get_grp_l(run, group, species_list)
                reactions_dict[run] = Reactions(r_path, grp_species_l, out)
            else:
                reactions_dict[run] = Reactions(r_path, species_list, out)
    print(f"Reactions instances has been created for runs : {runs} with group = {group} and out = "
          f"{out}")
    return reactions_dict


def get_pathways_inst(runs: List[str] = None, species_list: List[str] = None,
                      group: Tuple[str, int] = None, out: int = None,
                      nb_rnx_px_min: int = 0) -> Dict[str, 'Pathways']:
    """ Create Pathways instances in a dictionary.

    Parameters
    ----------
    runs: List[str], optional (default=None)
        List of the runs ID to consider
        If None, will be the list of all runs in the Runs path
    species_list: List[str], optional (default=None)
        List of species to consider for instances creation
    group: Tuple[str, int], optional (default=None)
        The group to consider for species filtering : Tuple(name of the group, column of the
        group in the organisms_file)
        If None will be all the species of each run
    out: int, optional (default=None)
        Number of species maximum not having the pathway for the pathway to be kept
        If None, will not filter the reactions
    nb_rnx_px_min: int, optional (default=0)
        Minimal number of reactions in a pathway for the pathway to be kept
        If not specified, will be 0 (=> no filter)

    Returns
    -------
    pathways_dict: Dict[str, 'Pathways']
        Dictionary of Pathways instances : Dict[ID of the run, Pathways instance of the run]
    """
    pathways_dict = {}
    if runs is None:
        runs = os.listdir(PATH_RUNS)
    for run in runs:
        p_path = os.path.join(PATH_RUNS, run, "analysis", "all", "pathways.tsv")
        if os.path.exists(p_path):
            if group is not None:
                grp_species_l = get_grp_l(run, group, species_list)
                pathways_dict[run] = Pathways(p_path, grp_species_l, out, nb_rnx_px_min)
            else:
                pathways_dict[run] = Pathways(p_path, species_list, out, nb_rnx_px_min)
        print(f"Pathways instances has been created for runs : {runs} with group = {group}, out = "
              f"{out} and minimal number of rnx in pw = {nb_rnx_px_min}")
    return pathways_dict


# ## Group comparisons and figures generation

STAT = (" mean", " med", " sd", " min", " max")
TO_CALCULATE = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")


def hist_comp(out_dir: str, df_comp: 'pd.DataFrame',
              df_grp_dict: Dict[Tuple[str, int], 'pd.DataFrame']):
    """Generate histogram figures for each group in df_grp_dict keys.

    Parameters
    ----------
    out_dir: str
        Path of the directory to save figures
    df_comp: 'pd.DataFrame'
        Data Frame of the comparison between all the groups
    df_grp_dict: Dict[Tuple[str, int], 'pd.DataFrame']
        Dictionary associating for each group, its data_frame of calculated stats
        Dict[Tuple[name of group, column], DataFrame of the group]

    Save a total of the number of groups of histogram figures in out_dir path.
    """
    for group, df_group in df_grp_dict.items():
        fig = plt.figure(figsize=[12.8, 9.6])
        i = 1
        for calc in TO_CALCULATE:
            axs = fig.add_subplot(2, 2, i)
            i += 1
            X = df_group[calc]
            axs.hist(X, len(X) - len(X) // 2, density=True, facecolor='g')
            axs.set_xlabel(calc)
            axs.set_ylabel('Density')
            mu = df_comp.loc[group[0], calc + STAT[0]]
            med = df_comp.loc[group[0], calc + STAT[1]]
            sd = df_comp.loc[group[0], calc + STAT[2]]
            mini = df_comp.loc[group[0], calc + STAT[3]]
            maxi = df_comp.loc[group[0], calc + STAT[4]]
            axs.set_title(f"Histogram of {calc} for the {group[0]} group\n"
                          f"$\mu={mu},\ \sigma={sd}$")
            axs.axvline(x=mu, color='b', label=f"$\mu={mu}$")
            axs.axvline(x=med, color='c', label=f"med={med}")
            axs.axvline(x=mini, color='m', label=f"min={mini}")
            axs.axvline(x=maxi, color='r', label=f"max={maxi}")
            axs.grid(True)
            axs.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        fig.tight_layout()
        path_fig = os.path.join(out_dir, f"{group[0]}_hist.png")
        plt.savefig(path_fig)


def boxplot_comp(out_dir: str, df_grp_dict: Dict[Tuple[str, int], 'pd.DataFrame']):
    """Generate boxplot figure comparing each group of df_grp_dict keys.

       Parameters
       ----------
       out_dir: str
           Path of the directory to save figures
       df_grp_dict: Dict[Tuple[str, int], 'pd.DataFrame']
           Dictionary associating for each group, its data_frame of calculated stats
           Dict[Tuple[name of group, column], DataFrame of the group]

       Save a boxplot figure in out_dir path.
       """
    nb_grp = len(df_grp_dict)
    i = 1
    str_grp_names = ""
    for group in df_grp_dict:
        if i == nb_grp - 1:
            str_grp_names += f"{group[0]} and "
            i += 1
        else:
            str_grp_names += f"{group[0]}, "
            i += 1
    fig = plt.figure(figsize=[12.8, 9.6])
    i = 1
    for calc in TO_CALCULATE:
        axs = fig.add_subplot(2, 2, i)
        i += 1
        x_lst = []
        labels_lst = []
        for group, df_group in df_grp_dict.items():
            x_lst.append(list(df_group[calc]))
            labels_lst.append(group[0])
        axs.boxplot(x_lst, labels=labels_lst)
        axs.set_ylabel(calc)
        axs.set_title(f"Boxplot of {calc} for the {str_grp_names[:-2]} groups")
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, f"boxplot.png"))


def compare_groups(run: str, groups_list: List[Tuple[str, int]], org_file: str, hist: bool = False,
                   boxplot: bool = False):
    """ Compare groups given according to : their number of genes, number of reactions, number of
    pathways with completion > 80% and number of pathways with completion = 100%. For each of this
    elements are calculated the mean, median, standard deviation, minimum and maximum between
    species of the group. Results are stored in tables and can be illustrated in histograms and/or
    boxplots figures.

    Parameters
    ----------
    run: str
        ID of the run to consider
    groups_list: List[Tuple[str, int]]
        List of groups to compare, groups are selected from groups specified in org_file
    org_file: str
        TSV file providing groups information
    hist: bool, optional (default=False)
        Indicate whether histograms output (PNG files) must be created
    boxplot: bool, optional (default=False)
        Indicate whether boxplot output (PNG file) must be created
    """
    path = os.path.join(PATH_RUNS, run, "analysis", "all")

    grp_str = ""
    for group in groups_list:
        grp_str += group[0] + "_"

    out_dir = os.path.join("output_data", "compare_groups", f"{run}_compare_{grp_str[:-1]}")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    res_file = os.path.join(out_dir, f"compare_{grp_str[:-1]}.tsv")

    col_stat = []
    for x in TO_CALCULATE:
        for y in STAT:
            col_stat.append(x + y)
    df_comp = pd.DataFrame(columns=col_stat, index=[group[0] for group in groups_list])

    df_grp_dict = {}
    for group in groups_list:
        grp_sp_list = get_grp_l(run, org_file, group)

        r_g = Reactions(os.path.join(path, "reactions.tsv"), grp_sp_list)
        p_g = Pathways(os.path.join(path, "pathways.tsv"), grp_sp_list)
        g_g = Genes(os.path.join(path, "genes.tsv"), grp_sp_list)

        g_res_file = os.path.join(out_dir, f"{group[0]}.tsv")

        df_g = pd.DataFrame(columns=TO_CALCULATE, index=grp_sp_list)

        for sp in grp_sp_list:
            df_g.loc[sp, TO_CALCULATE[0]] = g_g.nb_genes_species[sp]
            df_g.loc[sp, TO_CALCULATE[1]] = r_g.nb_reactions_sp[sp]
            df_g.loc[sp, TO_CALCULATE[2]] = p_g.get_pw_over_treshold(sp, 0.8)[sp][0]
            df_g.loc[sp, TO_CALCULATE[3]] = p_g.get_pw_complete(sp, False)[sp][0]

        df_grp_dict[group] = df_g
        df_g.to_csv(g_res_file, sep="\t")

        for calc in TO_CALCULATE:
            df_comp.loc[group[0], calc + STAT[0]] = round(float(np.mean(df_g[calc])), 1)
            df_comp.loc[group[0], calc + STAT[1]] = np.median(df_g[calc])
            df_comp.loc[group[0], calc + STAT[2]] = round(np.sqrt(np.var(df_g[calc])), 1)
            df_comp.loc[group[0], calc + STAT[3]] = min(df_g[calc])
            df_comp.loc[group[0], calc + STAT[4]] = max(df_g[calc])

        df_comp.to_csv(res_file, sep="\t")

    if hist:
        hist_comp(out_dir, df_comp, df_grp_dict)
    if boxplot:
        boxplot_comp(out_dir, df_grp_dict)


def intersect_rnx_groups(run: str, groups_list, org_file, percentage=True, venn_plot=False):
    file = os.path.join(PATH_RUNS, run, "analysis", "all", "reactions.tsv")
    group_rnx_dict = {}
    grp_str = ""
    for group in groups_list:
        grp_str += group[0] + "_"
    for group in groups_list:
        group_list = get_grp_l(run, org_file, group)
        r = Reactions(file, group_list)
        reactions_g = set(r.reactions_list)
        group_rnx_dict[group[0]] = reactions_g

    if percentage:
        union = set.union(*group_rnx_dict.values())
        intersect = set.intersection(*group_rnx_dict.values())
        union_nb = len(union)
        intersect_nb = len(intersect)
        txt_file = os.path.join(PATH_STUDY, "output_data", "compare_groups",
                                f"{run}_intersect_{grp_str}.txt")
        with open(txt_file, "w") as f:
            s = f"Over {union_nb} reactions, {intersect_nb} reactions in common = " \
                f"{round((intersect_nb/union_nb) * 100, 2)} % of the union."
            print(s)
            f.write(s)
            for group, reactions in group_rnx_dict.items():
                g_rnx_nb = len(reactions)
                s = f"{g_rnx_nb - intersect_nb} only present among the {group} group = " \
                    f"{round(((g_rnx_nb - intersect_nb)/union_nb) * 100, 2)} % of the union."
                print(s)
                f.write(s)

    if venn_plot:
        venn.venn(group_rnx_dict, cmap="plasma")
        plt.savefig(os.path.join(PATH_STUDY, "output_data", "compare_groups",
                                 f"{run}_intersect_{grp_str}.png"))
