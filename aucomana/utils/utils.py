import os
import pandas as pd
from typing import List, Set, Dict, Iterable


def create_dendrogram_groups_file(output):
    file = os.path.join(output, "dendro_tanglegrams", "dendrogram_groups.tsv")
    if not os.path.exists(file):
        with open(file, "w") as f:
            f.write("\t".join(["group name", "color", "element (B for branch, L for leave)"]))
        print(f"dendrogram_groups.tsv file created in path : {file}")


def get_abbr_name(name: str) -> str:
    """Abbreviate and returns the name of the species in format : X.xxxx
    Ex : Escherichia-coli returns E.coli

    Parameters
    ----------
    name: str
        name of the species

    Returns
    -------
    abbr_name: str
        Abbreviated name of the species
    """
    name = name.split("_")[0].split("-")
    if len(name) > 1:
        if len(name[1]) > 3:
            return f"{name[0][0].upper()}.{name[1][:4].lower()}"
        else:
            return f"{name[0][0].upper()}.{name[1].lower()}"
    else:
        return f"{name[0][0].upper()}{name[0][1:5].lower()}"


# ## Extract species according to groups from organisms file


def extract_groups(group_file: str) -> Dict[str, Set[str]]:
    dic_group = dict()
    with open(group_file, 'r') as f:
        for l in f:
            l = l.strip().split('\t')
            group = l[0]
            species = set(l[1:])
            if '' in species:
                species.remove('')
            dic_group[group] = species
    return dic_group


def get_grp_set(group_file: str, group: str or Iterable[str], species_list: Iterable[str] = None) \
        -> Set[str]:
    """ Select species according to the group they belong to. The groups must be specified in the
    "group_template.tsv" file.

    Parameters
    ----------
    group_file: str
        group_file path
    group: str or Iterable[str]
        The group or list of groups to consider
    species_list: Iterable[str], optional (default=None)
        List of species to consider to be filtered according to their group belonging
        If None, it will be all the species of the run


    Returns
    -------
    Set[str]
        Set of species corresponding to the intersection of the groups chosen
    """
    dic_groups = extract_groups(group_file)
    if group in dic_groups['all']:
        return {group}

    if species_list is not None:
        species_set = set(species_list)
    else:
        species_set = dic_groups['all']

    if type(group) == str:
        group = [group]

    for g in group:
        if g not in dic_groups.keys():
            raise ValueError(f"No group {g} in {group_file} file. Groups are : {list(dic_groups.keys())}.")
        species_set = species_set.intersection(dic_groups[g])
    return species_set
