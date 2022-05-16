import os
import pandas as pd

from analysis_runs.analysis import Analysis
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri


packnames = ('dendextend', 'pvclust', 'grDevices', 'ape')


def install_packages():
    utils = importr('utils')
    utils.chooseCRANmirror(ind=1)
    for package in packnames:
        if package not in rpy2.robjects.r['installed.packages']():
            utils.install_packages(package)


dendextend = importr(packnames[0])
pvclust = importr(packnames[1])
grdevices = importr(packnames[2])
ape = importr(packnames[3])


def create_dendro_groups_file(path_study):
    file = os.path.join(path_study, "../output_data", "dendro_tanglegrams", "dendro_groups.tsv")
    if not os.path.exists(file):
        with open(file, "w") as f:
            f.write("\t".join(["group name", "column", "color",
                               "element (B for branch, L for leave)"]))
        print(f"dendro_groups.tsv file created in path : {file}")
    else:
        print(f"dendro_groups.tsv file already exists in path : {file}")


def get_dendro_pvclust(path_study, df_binary, name, run, phylo_file=None, n_boot=100000):
    # Make pandas dataframe compatible with R dataframe.
    pandas2ri.activate()
    # Launch pvclust on the data silently and in parallel.
    print(f"Running pvclust with nboot = {n_boot}, this step may take a while.")
    result = pvclust.pvclust(df_binary, method_dist="binary", method_hclust="complete",
                             nboot=n_boot, quiet=True, parallel=True)
    print("Running pvclust : Done.")
    # Create the dendrogram picture.
    out_dir = os.path.join(path_study, "output_data", "dendro_tanglegrams", run, name)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        print(f"Directory {out_dir} already exists")
    out_file = os.path.join(out_dir, f"{name}_pvclust_dend.png")
    grdevices.png(file=out_file, width=2048, height=2048, pointsize=24)
    pvclust.plot_pvclust(result)
    grdevices.dev_off()
    print(f"pvclust dendrogram has been saved to : {out_file}")
    # Dendextend dendrogram
    dendro_groups_file = os.path.join(path_study, "output_data", "dendro_tanglegrams",
                                      "dendro_groups.tsv")
    create_dendextend(result, name, dendro_groups_file, run, out_dir, phylo_file)


def get_groups(d_grp, run):
    groups_branch = {}
    groups_leaves = {}
    for i in d_grp.index:
        group = (d_grp.loc[i, "group name"], d_grp.loc[i, "column"])
        color = d_grp.loc[i, "color"]
        element = d_grp.loc[i, "element (B for branch, L for leave)"]
        group_list = get_grp_l(run, group)
        if element == "B":
            groups_branch[group] = {"list": group_list, "color": color}
        elif element == "L":
            groups_leaves[group] = {"list": group_list, "color": color}
    return groups_branch, groups_leaves


def set_branch_color(dend, groups_branch):
    for gval in groups_branch.values():
        value = rpy2.robjects.StrVector(gval["list"])
        color = rpy2.robjects.StrVector([gval["color"]])
        dend = dendextend.set(dend, "by_labels_branches_col", value=value, TF_values=color,
                              quiet=True)
    return dend


def set_labels_color(dend, groups_branch):
    lab = rpy2.robjects.r['labels'](dend)
    lab_color = []
    for sp in lab:
        in_list = False
        for gval in groups_branch.values():
            if sp in gval["list"]:
                in_list = True
                lab_color.append(gval["color"])
        if not in_list:
            lab_color.append("black")
    lab_color = rpy2.robjects.StrVector(lab_color)
    dend = dendextend.set(dend, "labels_col", lab_color, quiet=True)
    return dend


def set_leaves_color(dend, groups_leaves):
    lab = rpy2.robjects.r['labels'](dend)
    leaves_color = []
    leaves_pch = []
    for sp in lab:
        in_list = False
        for gval in groups_leaves.values():
            if sp in gval["list"]:
                in_list = True
                leaves_color.append(gval["color"])
                leaves_pch.append(19)
        if not in_list:
            leaves_color.append("0")
            leaves_pch.append(0)
    leaves_color = rpy2.robjects.StrVector(leaves_color)
    leaves_pch = rpy2.robjects.IntVector(leaves_pch)
    dend = dendextend.set(dend, "leaves_col", leaves_color, quiet=True)
    dend = dendextend.set(dend, "leaves_pch", leaves_pch, quiet=True)
    return dend


def set_abbr_labels(dend):
    abbr_label = []
    for l in rpy2.robjects.r['labels'](dend):
        abbr_label.append(Analysis.get_abbr_name(l))
    abbr_label = rpy2.robjects.StrVector(abbr_label)
    dend = dendextend.set(dend, "labels", abbr_label, quiet=True)
    return dend


def create_dendextend(pvclust_res, name, dendro_groups_file, run, out_dir, phylo_file):
    print(f"Creating shaped dendextend dendrogram from {dendro_groups_file} parameters.")
    d_grp = pd.read_csv(dendro_groups_file, sep="\t")
    # Get groups
    groups_branch, groups_leaves = get_groups(d_grp, run)
    # Create dendrograms
    dend = rpy2.robjects.r['as.dendrogram'](pvclust_res)
    dend = dendextend.set(dend, "branches_lwd", 4)
    # Branch Color
    dend = set_branch_color(dend, groups_branch)
    # Labels Color
    dend = set_labels_color(dend, groups_branch)
    # Leaves Color
    dend = set_leaves_color(dend, groups_leaves)
    # Sort & Rename labels
    dend = rpy2.robjects.r['sort'](dend, type="nodes")
    dend = set_abbr_labels(dend)
    # Save figure
    title = f"{name} original phylogeny dendrogram"
    out_file = os.path.join(out_dir, f"{name}_dendextend_dend.png")
    grdevices.png(file=out_file, width=2048, height=2048, pointsize=24)
    rpy2.robjects.r.plot(dend, horiz=True, main=title)
    grdevices.dev_off()
    print(f"Shaped dendextend dendrogram has been saved to : {out_file}")
    # Compare with original phylogeny
    if phylo_file is not None:
        phylo = get_original_phylo_dend(phylo_file, name, out_dir, groups_branch, groups_leaves)
        get_tanglegram(phylo, dend, name, out_dir)


def get_original_phylo_dend(phylo_file, name, out_dir, groups_branch, groups_leaves):
    print(f"Creating phylogeny dendrogram from {phylo_file} information.")
    phylo = rpy2.robjects.r['read.nexus'](phylo_file)
    phylo = ape.chronos(phylo, quiet=True)
    phylo = rpy2.robjects.r['as.dendrogram'](phylo)
    phylo = dendextend.set(phylo, "branches_lwd", 4)
    # Branch Color
    phylo = set_branch_color(phylo, groups_branch)
    # Labels Color
    phylo = set_labels_color(phylo, groups_branch)
    # Leaves Color
    phylo = set_leaves_color(phylo, groups_leaves)
    # Sort & Rename
    phylo = rpy2.robjects.r['sort'](phylo, type="nodes")
    phylo = set_abbr_labels(phylo)
    # Save figure
    title = f"{name} metabolic dendrogram"
    out_file = os.path.join(out_dir, f"{name}_phylo_dend.png")
    grdevices.png(file=out_file, width=2048, height=2048, pointsize=24)
    rpy2.robjects.r.plot(phylo, horiz=True, main=title)
    grdevices.dev_off()
    print(f"Phylogeny dendrogram has been saved to : {out_file}")
    return phylo


def get_tanglegram(phylo, dend, name, out_dir):
    print("Creating tanglegram comparing metabolic dendrogram to phylogeny dendrogram.")
    d1 = dendextend.intersect_trees(phylo, dend, quiet=True)
    col_lines = dendextend.labels_col(d1[0])

    out_file = os.path.join(out_dir, f"{name}_tanglegram.png")
    grdevices.png(file=out_file, width=1500, height=1000, pointsize=24)
    title_right = f"{name} Metabolic Dendrogram"
    dendextend.tanglegram(d1, lwd=2, color_lines=col_lines, main_left="Original phylogeny",
                          main_right=title_right, quiet=True)
    grdevices.dev_off()
    print(f"Tanglegram has been saved to : {out_file}")
    get_similarity_indicators(d1, name, out_dir)


def get_similarity_indicators(d1, name, out_dir):
    cor_coph = dendextend.cor_cophenetic(d1, method_coef="pearson")
    cor_bakers_gamma = dendextend.cor_bakers_gamma(d1)
    out_file = os.path.join(out_dir, f"{name}_similarity_indicators.tsv")
    with open(out_file, 'w') as f:
        f.write("Correlation cophenetic\tCorrelation bakers gamma\n")
        f.write(f"{cor_coph[0]}\t{cor_bakers_gamma[0]}")
    print(f"Correlation indicators table has been saved to : {out_file}")
