import os
import pandas as pd
from typing import List, Set, Iterable
from utils.utils import *
from utils.pathways import Pathways


class AuCoMAna:
    """
    Attributes
    ----------
    group_template: str
    reactions_tsv: str
    pathways_tsv: str
    genes_tsv: str
    metabolites_tsv: str
    """

    STAT = (" mean", " med", " sd", " min", " max")
    TO_CALCULATE = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")

    def __int__(self, compare_path: str, group_template: str):
        """ Init the AuCoMAna class.

        Parameters
        ----------
        compare_path: str
            Path of the directory with compare padmets results
        group_template: str
            Groups template tsv file
        """
        self.group_template = group_template
        self.reactions_tsv = os.path.join(compare_path, 'reactions.tsv')
        self.pathways_tsv = os.path.join(compare_path, 'pathways.tsv')
        self.genes_tsv = os.path.join(compare_path, 'genes.tsv')
        self.metabolites_tsv = os.path.join(compare_path, 'metabolites.tsv')

    def group_pathway_completion_comparison(self, groups_comp: Iterable[str], group_analysis: str = 'all'):
        species_analysis = get_grp_set(self.group_template, group_analysis)

        pathway = Pathways(file_pathways_tsv=self.pathways_tsv,
                           species_list=species_analysis)
        dic_groups = dict()
        for group in groups_comp:
            dic_groups[group] = get_grp_set(self.group_template, group, species_analysis)

        abs = pathway.get_pw_absent(species=species_analysis, unique=True)
        print(abs)
        # abs_sr = []
        # abs_lr = []
        # abs_all = []
        # abs_lam = []
        #
        # for sp, val in abs.items():
        #     print(sp, val[0])
        #     abs_all.append(val[0])
        #     if sp in sr:
        #         abs_sr.append(val[0])
        #     if sp in lr:
        #         abs_lr.append(val[0])
        #     if sp == LELS:
        #         abs_lam.append(val[0])
        #
        # mini = PATHWAYS[R04].get_pw_min(species=species, unique=True)
        # min_sr = []
        # min_lr = []
        # min_all = []
        # min_lam = []
        #
        # for sp, val in mini.items():
        #     print(sp, val[0])
        #     min_all.append(val[0])
        #     if sp in sr:
        #         min_sr.append(val[0])
        #     if sp in lr:
        #         min_lr.append(val[0])
        #     if sp == LELS:
        #         min_lam.append(val[0])
        #
        # inc = PATHWAYS[R04].get_pw_incomplete(species=species, unique=True)
        # inc_sr = []
        # inc_lr = []
        # inc_all = []
        # inc_lam = []
        #
        # for sp, val in inc.items():
        #     print(sp, val[0])
        #     inc_all.append(val[0])
        #     if sp in sr:
        #         inc_sr.append(val[0])
        #     if sp in lr:
        #         inc_lr.append(val[0])
        #     if sp == LELS:
        #         inc_lam.append(val[0])
        #
        # # print(np.mean(loss_l))
        # # print(np.sqrt(np.var(loss_l)))
        #
        # col = ["absent", "minimal", "incomplete"]
        # ind = ["all", "LR", "SR", "L.elsb"]
        # df = pd.DataFrame(columns=ind, index=col)
        # df.loc["absent", "all"] = round(float(np.mean(abs_all)), 2)
        # df.loc["absent", "LR"] = round(float(np.mean(abs_lr)), 2)
        # df.loc["absent", "SR"] = round(float(np.mean(abs_sr)), 2)
        # df.loc["absent", "L.elsb"] = round(float(np.mean(abs_lam)), 2)
        #
        # df.loc["minimal", "all"] = round(float(np.mean(min_all)), 2)
        # df.loc["minimal", "LR"] = round(float(np.mean(min_lr)), 2)
        # df.loc["minimal", "SR"] = round(float(np.mean(min_sr)), 2)
        # df.loc["minimal", "L.elsb"] = round(float(np.mean(min_lam)), 2)
        #
        # df.loc["incomplete", "all"] = round(float(np.mean(inc_all)), 2)
        # df.loc["incomplete", "LR"] = round(float(np.mean(inc_lr)), 2)
        # df.loc["incomplete", "SR"] = round(float(np.mean(inc_sr)), 2)
        # df.loc["incomplete", "L.elsb"] = round(float(np.mean(inc_lam)), 2)
        #
        # print(df)
        # ax = df.plot.bar()
        # for container in ax.containers:
        #     ax.bar_label(container)
        # ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        # ax.set_ylabel("Number of pathways")
        # ax.set_title("Unique absent/minimal/incomplete pathway for L.eslb compared to all/LR/SR groups mean")
        # plt.xticks(rotation=0)
        # plt.tight_layout()
        # plt.show()

group_file = '../Runs/run62/analysis/group_template.tsv'
compare_dir = '../Runs/run62/analysis/'

A = AuCoMAna(compare_dir, group_file)
print(A.group_template)



