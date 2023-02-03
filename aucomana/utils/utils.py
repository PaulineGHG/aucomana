import os
import pandas as pd
from typing import List, Set


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

def get_grp_set(group_file: str, group: str, species_list: List[str] = None) -> Set[str]:
    """ Select species according to the group they belong to. The groups must be specified in the
    "group_template.tsv" file in path : <path runs>/<runID>/analysis/

    Parameters
    ----------
    group_file: str
        group_file path
    group: str
        The group to consider
    species_list: List[str], optional (default=None)
        List of species to consider to be filtered according to their group belonging
        If None, it will be all the species of the run

    Returns
    -------
    Set[str]
        Set of species corresponding to the group chosen
    """
    df = pd.read_csv(group_file, sep="\t", index_col=0)
    if group not in list(df.index):
        raise ValueError(f"No group {group} in {group_file} file. Groups are : {list(df.index)}.")
    group_set = set(df.loc[group])
    if species_list is not None:
        species_set = set(species_list)
    else:
        species_set = set(df.columns)
    return group_set.intersection(species_set)
