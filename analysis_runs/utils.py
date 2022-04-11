import os
from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.genes import Genes
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
