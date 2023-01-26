"""
Dendrogram class
"""
import os
import pandas as pd
import rpy2

from aucomana.analysis import Analysis
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

packnames = ('dendextend', 'pvclust', 'grDevices', 'ape')


class Dendrogram:
    dendextend = importr(packnames[0])
    pvclust = importr(packnames[1])
    grdevices = importr(packnames[2])
    ape = importr(packnames[3])

    def __init__(self, path_runs, path_study, df_binary, run, name, phylo_file):
        self.path_runs = path_runs
        self.path_study = path_study
        self.df_binary = df_binary
        self.run = run
        self.name = name
        self.phylo_file = phylo_file
        self.out_dir = os.path.join(self.path_study, "output_data", "dendro_tanglegrams", self.run,
                                    self.name)
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.dendro_groups_file = os.path.join(self.path_study, "output_data", "dendro_tanglegrams",
                                               "dendro_groups.tsv")

    def get_dendro_pvclust(self, n_boot=10000):
        # Make pandas dataframe compatible with R dataframe.
        pandas2ri.activate()
        # Launch pvclust on the data silently and in parallel.
        print(f"Running pvclust with nboot = {n_boot}, this step may take a while.")
        result = self.pvclust.pvclust(self.df_binary, method_dist="binary", method_hclust="complete",
                                 nboot=n_boot, quiet=True, parallel=True)
        print("Running pvclust : Done.")
        # Create the dendrogram picture.

        out_file = os.path.join(self.out_dir, f"{self.name}_pvclust_dend.png")
        self.grdevices.png(file=out_file, width=1500, height=1000, pointsize=24)
        self.pvclust.plot_pvclust(result)
        self.grdevices.dev_off()
        print(f"pvclust dendrogram has been saved to : {out_file}")
        # Dendextend dendrogram
        self.__create_dendextend(result)

    def __get_groups(self, d_grp):
        groups_branch = {}
        groups_leaves = {}
        for i in d_grp.index:
            group = d_grp.loc[i, "group name"]
            color = d_grp.loc[i, "color"]
            element = d_grp.loc[i, "element (B for branch, L for leave)"]
            a = Analysis(self.path_runs, self.path_study)
            group_list = a.get_grp_set(run=self.run, group=group)
            if element == "B":
                groups_branch[group] = {"list": group_list, "color": color}
            elif element == "L":
                groups_leaves[group] = {"list": group_list, "color": color}
        return groups_branch, groups_leaves

    def __set_branch_color(self, dend, groups_branch):
        for gval in groups_branch.values():
            value = rpy2.robjects.StrVector(gval["list"])
            color = rpy2.robjects.StrVector([gval["color"]])
            dend = self.dendextend.set(dend, "by_labels_branches_col", value=value, TF_values=color,
                                       quiet=True)
        return dend

    def __set_labels_color(self, dend, groups_branch):
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
        dend = self.dendextend.set(dend, "labels_col", lab_color, quiet=True)
        return dend

    def __set_leaves_color(self, dend, groups_leaves):
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
        dend = self.dendextend.set(dend, "leaves_col", leaves_color, quiet=True)
        dend = self.dendextend.set(dend, "leaves_pch", leaves_pch, quiet=True)
        return dend

    # !!! If the abbreviatio gives duplicate names, untangle and tanglegram will give error.
    def __set_abbr_labels(self, dend):
        abbr_label = []
        for l in rpy2.robjects.r['labels'](dend):
            abbr_label.append(Analysis.get_abbr_name(l))
        abbr_label = rpy2.robjects.StrVector(abbr_label)
        dend = self.dendextend.set(dend, "labels", abbr_label, quiet=True)
        return dend

    def __create_dendextend(self, pvclust_res):
        print(f"Creating shaped dendextend dendrogram from {self.dendro_groups_file} parameters.")
        d_grp = pd.read_csv(self.dendro_groups_file, sep="\t")
        # Get groups
        groups_branch, groups_leaves = self.__get_groups(d_grp)
        # Create dendrograms
        dend = rpy2.robjects.r['as.dendrogram'](pvclust_res)
        dend = self.dendextend.set(dend, "branches_lwd", 4)
        # Branch Color
        dend = self.__set_branch_color(dend, groups_branch)
        # Labels Color
        dend = self.__set_labels_color(dend, groups_branch)
        # Leaves Color
        dend = self.__set_leaves_color(dend, groups_leaves)
        # Sort & Rename labels
        dend = rpy2.robjects.r['sort'](dend, type="nodes")
        dend = self.__set_abbr_labels(dend)
        # Save figure
        title = f"{self.name} original phylogeny dendrogram"
        out_file = os.path.join(self.out_dir, f"{self.name}_dendextend_dend.png")
        self.grdevices.png(file=out_file, width=3000, height=2000, pointsize=24)
        rpy2.robjects.r.plot(dend, horiz=True, main=title)
        self.grdevices.dev_off()
        print(f"Shaped dendextend dendrogram has been saved to : {out_file}")
        # Compare with original phylogeny
        if self.phylo_file is not None:
            phylo = self.__get_original_phylo_dend(groups_branch, groups_leaves)
            self.__get_tanglegram(phylo, dend)

    def __get_original_phylo_dend(self, groups_branch, groups_leaves):
        print(f"Creating phylogeny dendrogram from {self.phylo_file} information.")
        phylo = rpy2.robjects.r['read.tree'](self.phylo_file)

        # Make the tree ultrametric
        phylo = self.ape.chronos(phylo, quiet=True)

        # Make the tree dichotomous
        phylo = self.ape.multi2di(phylo)
        # phylo = self.ape.unroot(phylo)
        # phylo = self.ape.root(phylo, outgroup="Heterosigma-akashiwo", resolve_root=True)
        # print(self.ape.is_rooted(phylo))

        # Convert tree to dendrogram
        phylo = rpy2.robjects.r['as.dendrogram'](phylo)

        # Color branch / leaves / labels
        phylo = self.dendextend.set(phylo, "branches_lwd", 4)  # Thickness
        phylo = self.__set_branch_color(phylo, groups_branch)
        phylo = self.__set_labels_color(phylo, groups_branch)
        phylo = self.__set_leaves_color(phylo, groups_leaves)

        # Sort & Rename
        phylo = rpy2.robjects.r['sort'](phylo, type="nodes")
        phylo = self.__set_abbr_labels(phylo)

        # Save figure
        title = f"{self.name} metabolic dendrogram"
        out_file = os.path.join(self.out_dir, f"{self.name}_phylo_dend.png")
        self.grdevices.png(file=out_file, width=3000, height=2000, pointsize=24)
        rpy2.robjects.r.plot(phylo, horiz=True, main=title)
        self.grdevices.dev_off()
        print(f"Phylogeny dendrogram has been saved to : {out_file}")
        return phylo

    def __get_tanglegram(self, phylo, dend):
        print("Creating tanglegram comparing metabolic dendrogram to phylogeny dendrogram.")
        d1 = self.dendextend.intersect_trees(phylo, dend, quiet=True)
        col_lines = self.dendextend.labels_col(d1[0])

        out_file = os.path.join(self.out_dir, f"{self.name}_tanglegram.png")
        self.grdevices.png(file=out_file, width=3000, height=2000, pointsize=24)
        title_right = f"{self.name} Metabolic Dendrogram"
        d1 = self.dendextend.untangle(d1)
        self.dendextend.tanglegram(d1, lwd=4, color_lines=col_lines, main_left="Original phylogeny",
                                   main_right=title_right, quiet=True)
        self.grdevices.dev_off()
        print(f"Tanglegram has been saved to : {out_file}")
        self.__get_similarity_indicators(d1)

    def __get_similarity_indicators(self, d1):
        cor_coph = self.dendextend.cor_cophenetic(d1, method_coef="pearson")
        cor_bakers_gamma = self.dendextend.cor_bakers_gamma(d1)
        out_file = os.path.join(self.out_dir, f"{self.name}_similarity_indicators.tsv")
        with open(out_file, 'w') as f:
            f.write("Correlation cophenetic\tCorrelation bakers gamma\n")
            f.write(f"{cor_coph[0]}\t{cor_bakers_gamma[0]}")
        print(f"Correlation indicators table has been saved to : {out_file}")
