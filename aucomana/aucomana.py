"""
AuCoMAna class
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Iterable
from aucomana.utils.utils import get_grp_set
from aucomana.utils.pathways import Pathways
from aucomana.utils.reactions import Reactions


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

    def __init__(self, compare_path: str, group_template: str):
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


    def group_reactions_comparison(self, groups_comp: Iterable[str], group_analysis: str = 'all'):
        species_analysis = get_grp_set(self.group_template, group_analysis)

        reaction = Reactions(file_reactions_tsv=self.reactions_tsv,
                             species_list=list(species_analysis))
        dic_groups = dict()
        for group in list(groups_comp) + [group_analysis]:
            dic_groups[group] = get_grp_set(self.group_template, (group, group_analysis))

        absent = reaction.get_rxn_absent(species=species_analysis, unique=True)
        abs_groups = dict()

        for group in list(groups_comp) + [group_analysis]:
            abs_groups[group] = list()

            for sp in dic_groups[group]:
                abs_groups[group].append(absent[sp][0])

        ind = ["absent"]
        col = list(abs_groups.keys())
        df = pd.DataFrame(columns=col, index=ind)
        for group, numbers_list in abs_groups.items():
            df.loc[ind[0], group] = round(float(np.mean(numbers_list)), 2)

        ax = df.plot.bar()
        for container in ax.containers:
            ax.bar_label(container)
        ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        ax.set_ylabel("Number of pathways")
        ax.set_title(f'Unique absent reactions for {"/".join(groups_comp)} groups')
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.show()


    def group_pathway_completion_comparison(self, groups_comp: Iterable[str], group_analysis: str = 'all'):
        species_analysis = get_grp_set(self.group_template, group_analysis)

        pathway = Pathways(file_pathways_tsv=self.pathways_tsv,
                           species_list=list(species_analysis))
        dic_groups = dict()
        for group in list(groups_comp) + [group_analysis]:
            dic_groups[group] = get_grp_set(self.group_template, (group, group_analysis))

        absent = pathway.get_pw_absent(species=species_analysis, unique=True)
        mini = pathway.get_pw_min(species=species_analysis, unique=True)
        incomplete = pathway.get_pw_incomplete(species=species_analysis, unique=True)

        abs_groups = dict()
        min_groups = dict()
        inc_groups = dict()

        for group in list(groups_comp) + [group_analysis]:
            abs_groups[group] = list()
            min_groups[group] = list()
            inc_groups[group] = list()

            for sp in dic_groups[group]:
                abs_groups[group].append(absent[sp][0])
                min_groups[group].append(mini[sp][0])
                inc_groups[group].append(incomplete[sp][0])

        ind = ["absent", "minimal", "incomplete"]
        col = list(abs_groups.keys())
        df = pd.DataFrame(columns=col, index=ind)
        for group, numbers_list in abs_groups.items():
            df.loc[ind[0], group] = round(float(np.mean(numbers_list)), 2)
        for group, numbers_list in min_groups.items():
            df.loc[ind[1], group] = round(float(np.mean(numbers_list)), 2)
        for group, numbers_list in inc_groups.items():
            df.loc[ind[2], group] = round(float(np.mean(numbers_list)), 2)

        ax = df.plot.bar()
        for container in ax.containers:
            ax.bar_label(container)
        ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
        ax.set_ylabel("Number of pathways")
        ax.set_title(f'Unique absent/minimal/incomplete pathway for {"/".join(groups_comp)} groups')
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.show()



