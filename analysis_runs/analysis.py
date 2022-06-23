import os
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import venn

from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.genes import Genes
from analysis_runs.metabolites import Metabolites
from analysis_runs.rename_padmets_id import make_automaton, apply_automaton, get_dict
from typing import Tuple, List, Dict, Set


class Analysis:
    STAT = (" mean", " med", " sd", " min", " max")
    TO_CALCULATE = ("nb_genes", "nb_rnx", "nb_pw > 80%", "nb_pw 100%")
# ## Various utils fnc

    def __init__(self, path_runs, path_study):
        self.path_runs = path_runs
        self.path_study = path_study
        self.__create_folders()

    # Create folders to save output files

    def __create_folders(self):
        """ Create the folder hierarchy where the outputs will be stored. Folders created in <path_study> directory
        """
        arbo = [['output_data'],
                ['output_data', 'compare_groups'],
                ['output_data', 'dendro_tanglegrams'],
                ['output_data', 'pathways_data'],
                ['output_data', 'pathways_data', 'binary_df'],
                ['output_data', 'reactions_data'],
                ['output_data', 'reactions_data', 'common_reac'],
                ['output_data', 'reactions_data', 'common_reac', 'union'],
                ['output_data', 'reactions_data', 'common_reac', 'intersection'],
                ['output_data', 'reactions_data', 'genes_assoc'],
                ['output_data', 'genes_data'],
                ['output_data', 'metabolites_data'],
                ['output_data', 'renamed_id_padmet']]
        for folder in arbo:
            folder_path = os.path.join(self.path_study, *folder)
            if not os.path.exists(folder_path):
                os.mkdir(folder_path)
                print(f'Folder {folder_path} created')

    @staticmethod
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
        if len(name) > 1:
            if len(name[1]) > 3:
                return f"{name[0][0].upper()}.{name[1][:4].lower()}"
            else:
                return f"{name[0][0].upper()}.{name[1].lower()}"
        else:
            return f"{name[0][0].upper()}{name[0][1:5].lower()}"

    # ## Extract species according to groups from organisms file

    def get_grp_set(self, run: str, group: str,
                    species_list: List[str] = None) -> Set[str]:
        """ Select species according to the group they belong to. The groups must be specified in the
        "group_template.tsv" file in path : <path runs>/<runID>/analysis/

        Parameters
        ----------
        run: str
            ID of the run
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
        grp_template_f = os.path.join(self.path_runs, run, "analysis",
                                      "group_template.tsv")
        df = pd.read_csv(grp_template_f, sep="\t", index_col=0)
        if group not in list(df.index):
            raise ValueError(f"No group {group} in {grp_template_f} file. Groups are : {list(df.index)}.")
        group_set = set(df.loc[group])
        if species_list is not None:
            species_set = set(species_list)
        else:
            species_set = set(df.columns)
        return group_set.intersection(species_set)

    # ## Create instances

    def reactions(self, run: str, species_list: List[str] = None, group: str = None, out: int = None) -> 'Reactions':
        """ Create Reactions instance.

        Parameters
        ----------
        run: str
            The run ID to consider
        species_list: List[str], optional (default=None)
            List of species to consider for instance creation
        group: str, optional (default=None)
            The group to consider for species filtering
            If None will be all the species in species_list
        out: int, optional (default=None)
            Number of species maximum not having the reaction for the reaction to be kept
            If None, will not filter the reactions

        Returns
        -------
        RXN : Reactions
            Reactions instance
        """
        r_path = os.path.join(self.path_runs, run, "analysis", "all", "reactions.tsv")
        if os.path.exists(r_path):
            if group is not None:
                grp_species_l = list(self.get_grp_set(run, group, species_list))
                RXN = Reactions(self.path_runs, self.path_study, r_path, grp_species_l, out)
                print(f"Reactions instance has been created for run : {run} with group = {group} and out = {out}")
                return RXN
            else:
                RXN = Reactions(self.path_runs, self.path_study, r_path, species_list, out)
                print(f"Reactions instances has been created for run : {run} and out = {out}")
                return RXN
        else:
            raise OSError(f"No file reactions.tsv in path {r_path}")

    def pathways(self, run: str, species_list: List[str] = None, group: str = None, out: int = None,
                 nb_rxn_pw_min: int = 0) -> 'Pathways':
        """ Create Pathways instance.

        Parameters
        ----------
        run: str
            The run ID to consider
        species_list: List[str], optional (default=None)
            List of species to consider for instance creation
        group: str, optional (default=None)
            The group to consider for species filtering
            If None will be all the species in species_list
        out: int, optional (default=None)
            Number of species maximum not having the pathway for the pathway to be kept
            If None, will not filter the reactions
        nb_rxn_pw_min: int, optional (default=0)
            Minimal number of reactions in a pathway for the pathway to be kept
            If not specified, will be 0 (<==> no filter)

        Returns
        -------
        PW: Pathways
            Pathways instance
        """
        p_path = os.path.join(self.path_runs, run, "analysis", "all", "pathways.tsv")
        if os.path.exists(p_path):
            if group is not None:
                grp_species_l = list(self.get_grp_set(run, group, species_list))
                PW = Pathways(self.path_runs, self.path_study, p_path, grp_species_l, out, nb_rxn_pw_min)
                print(f"Pathways instance has been created for run : {run} with group = {group}, out = "
                      f"{out} and minimal number of reactions in pathway = {nb_rxn_pw_min}")
                return PW
            else:
                PW = Pathways(self.path_runs, self.path_study, p_path, species_list, out, nb_rxn_pw_min)
                print(f"Pathways instances has been created for run : {run}, out = "
                      f"{out} and minimal number of reactions in pathway = {nb_rxn_pw_min}")
                return PW
        else:
            raise OSError(f"No file reactions.tsv in path {p_path}")

    def genes(self, run: str, species_list: List[str] = None, group: str = None) -> 'Genes':
        """ Create Genes instance.

        Parameters
        ----------
        run: str
            The run ID to consider
        species_list: List[str], optional (default=None)
            List of species to consider for instance creation
        group: str, optional (default=None)
            The group to consider for species filtering
            If None will be all the species in species_list

        Returns
        -------
        GN: Genes
            Genes instance
        """
        g_path = os.path.join(self.path_runs, run, "analysis", "all", "genes.tsv")
        if os.path.exists(g_path):
            if group is not None:
                grp_species_l = list(self.get_grp_set(run, group, species_list))
                GN = Genes(g_path, grp_species_l)
                print(f"Genes instance has been created for run : {run} with group = {group}")
                return GN
            else:
                GN = Genes(g_path, species_list)
                print(f"Genes instance has been created for run : {run}")
                return GN
        else:
            raise OSError(f"No file reactions.tsv in path {g_path}")

    def metabolites(self, run: str, species_list: List[str] = None, group: str = None) -> 'Metabolites':
        """ Create Metabolites instance.

        Parameters
        ----------
        run: str
            The runs ID to consider
        species_list: List[str], optional (default=None)
            List of species to consider for instance creation
        group: str, optional (default=None)
            The group to consider for species filtering
            If None will be all the species in species_list

        Returns
        -------
        MB: Metabolites
            Metabolites instance
        """
        m_path = os.path.join(self.path_runs, run, "analysis", "all", "metabolites.tsv")
        if os.path.exists(m_path):
            if group is not None:
                grp_species_l = list(self.get_grp_set(run, group, species_list))
                MB = Metabolites(m_path, grp_species_l)
                print(f"Metabolites instances has been created for run : {run} with group = {group}")
                return MB
            else:
                MB = Metabolites(m_path, species_list)
                print(f"Metabolites instances has been created for run : {run}")
                return MB
        else:
            raise OSError(f"No file reactions.tsv in path {m_path}")

    # ## Group comparisons and figures generation

    def __hist_comp(self, out_dir: str, df_comp: 'pd.DataFrame', df_grp_dict: Dict[str, 'pd.DataFrame']):
        """Generate histogram figures for each group in df_grp_dict keys.

        Parameters
        ----------
        out_dir: str
            Path of the directory to save figures
        df_comp: pd.DataFrame
            Data Frame of the comparison between all the groups
        df_grp_dict: Dict[str, pd.DataFrame]
            Dictionary associating for each group, its data_frame of calculated stats
            Dict[name of group, DataFrame of the group]

        Save a total of the number of groups of histogram figures in out_dir path.
        """
        for group, df_group in df_grp_dict.items():
            fig = plt.figure(figsize=[12.8, 9.6])
            i = 1
            for calc in self.TO_CALCULATE:
                axs = fig.add_subplot(2, 2, i)
                i += 1
                X = df_group[calc]
                axs.hist(X, len(X) - len(X) // 2, density=True, facecolor='g')
                axs.set_xlabel(calc)
                axs.set_ylabel('Density')
                mu = df_comp.loc[group, calc + self.STAT[0]]
                med = df_comp.loc[group, calc + self.STAT[1]]
                sd = df_comp.loc[group, calc + self.STAT[2]]
                mini = df_comp.loc[group, calc + self.STAT[3]]
                maxi = df_comp.loc[group, calc + self.STAT[4]]
                axs.set_title(f"Histogram of {calc} for the {group[0]} group\n"
                              f"$\mu={mu},\ \sigma={sd}$")
                axs.axvline(x=mu, color='b', label=f"$\mu={mu}$")
                axs.axvline(x=med, color='c', label=f"med={med}")
                axs.axvline(x=mini, color='m', label=f"min={mini}")
                axs.axvline(x=maxi, color='r', label=f"max={maxi}")
                axs.grid(True)
                axs.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
            fig.tight_layout()
            path_fig = os.path.join(out_dir, f"{group}_hist.png")
            plt.savefig(path_fig)

    def __boxplot_comp(self, out_dir: str, df_grp_dict: Dict[str, 'pd.DataFrame']):
        """Generate boxplot figure comparing each group of df_grp_dict keys.

           Parameters
           ----------
           out_dir: str
               Path of the directory to save figures
           df_grp_dict: Dict[str, pd.DataFrame]
               Dictionary associating for each group, its data_frame of calculated stats
               Dict[name of group, DataFrame of the group]

           Save a boxplot figure in out_dir path.
           """
        nb_grp = len(df_grp_dict)
        i = 1
        str_grp_names = ""
        for group in df_grp_dict:
            if i == nb_grp - 1:
                str_grp_names += f"{group} and "
                i += 1
            else:
                str_grp_names += f"{group}, "
                i += 1
        fig = plt.figure(figsize=[12.8, 9.6])
        i = 1
        for calc in self.TO_CALCULATE:
            axs = fig.add_subplot(2, 2, i)
            i += 1
            x_lst = []
            labels_lst = []
            for group, df_group in df_grp_dict.items():
                x_lst.append(list(df_group[calc]))
                labels_lst.append(group)
            axs.boxplot(x_lst, labels=labels_lst)
            axs.set_ylabel(calc)
            axs.set_title(f"Boxplot of {calc} for the {str_grp_names[:-2]} groups")
        fig.tight_layout()
        plt.savefig(os.path.join(out_dir, f"boxplot.png"))

    def compare_groups(self, run: str, groups_list: List[str], hist: bool = False,
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
        groups_list: List[str]
            List of groups to compare
        hist: bool, optional (default=False)
            Indicate whether histograms output (PNG files) must be created
        boxplot: bool, optional (default=False)
            Indicate whether boxplot output (PNG file) must be created
        """
        # Input path of TSV files
        in_path = os.path.join(self.path_runs, run, "analysis", "all")

        # String concatenating named of each groups
        grp_str = ""
        for group in groups_list:
            grp_str += group + "_"

        # Create output directory
        out_dir = os.path.join(self.path_study, "output_data", "compare_groups", f"{run}_compare_{grp_str[:-1]}")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        # Path of the compare_gX.tsv file
        res_file = os.path.join(out_dir, f"compare_{grp_str[:-1]}.tsv")

        # Create main DataFrame with groups as index and stats to calculate as columns
        col_stat = []
        for x in self.TO_CALCULATE:
            for y in self.STAT:
                col_stat.append(x + y)
        df_comp = pd.DataFrame(columns=col_stat, index=groups_list)

        # Create dictionary storing for each group his DataFrame
        df_grp_dict = {}
        for group in groups_list:
            grp_sp_list = list(self.get_grp_set(run, group))
            # Create Reactions, Pathways and Genes instances for the group
            r_g = Reactions(self.path_runs, self.path_study, os.path.join(in_path, "reactions.tsv"), grp_sp_list)
            p_g = Pathways(self.path_runs, self.path_study, os.path.join(in_path, "pathways.tsv"), grp_sp_list)
            g_g = Genes(os.path.join(in_path, "genes.tsv"), grp_sp_list)
            # Path of group.tsv file (group stats)
            g_res_file = os.path.join(out_dir, f"{group}.tsv")
            # Create DataFrame with species of the group as index and stats to calculate as columns
            df_g = pd.DataFrame(columns=self.TO_CALCULATE, index=grp_sp_list)
            # Calculate stats for each species and fill the DataFrame
            for sp in grp_sp_list:
                df_g.loc[sp, self.TO_CALCULATE[0]] = g_g.nb_genes_species[sp]
                df_g.loc[sp, self.TO_CALCULATE[1]] = r_g.nb_reactions_sp[sp]
                df_g.loc[sp, self.TO_CALCULATE[2]] = p_g.get_pw_over_treshold(sp, 0.8)[sp][0]
                df_g.loc[sp, self.TO_CALCULATE[3]] = p_g.get_pw_complete(sp, False)[sp][0]
            # Store DataFrame in dict and save it in TSV file
            df_grp_dict[group] = df_g
            df_g.to_csv(g_res_file, sep="\t")
            print(f"{group} stats saved to path : {g_res_file}")

            # Fill the main DataFrame with stats for the group
            for calc in self.TO_CALCULATE:
                df_comp.loc[group, calc + self.STAT[0]] = round(float(np.mean(df_g[calc])), 1)
                df_comp.loc[group, calc + self.STAT[1]] = np.median(df_g[calc])
                df_comp.loc[group, calc + self.STAT[2]] = round(np.sqrt(np.var(df_g[calc])), 1)
                df_comp.loc[group, calc + self.STAT[3]] = min(df_g[calc])
                df_comp.loc[group, calc + self.STAT[4]] = max(df_g[calc])

        # Save main DataFrame in TSV file
        df_comp.to_csv(res_file, sep="\t")
        print(f"groups stats comparisons saved to path : {res_file}")

        # Generate histograms for each groups
        if hist:
            self.__hist_comp(out_dir, df_comp, df_grp_dict)
        # Generate boxplot comparing all groups
        if boxplot:
            self.__boxplot_comp(out_dir, df_grp_dict)

    def intersect_rnx_groups(self, run: str, groups_list: List[str], percentage: bool = True, venn_plot: bool = False):
        """ Intersects the set of reactions of groups to compare them. Can show percentages comparisons and/or venn
        plot.

        Parameters
        ----------
        run : str
            The ID of the run to consider
        groups_list : List[str]
            List of the groups to intersect (max 6 if venn plot)
        percentage : bool, optional (default=True)
            True to save txt file of the percentages
        venn_plot : bool, optional (default=True)
            True to save venn plot of the intersection
        """
        # Path to the reactions.tsv file
        rxn_file = os.path.join(self.path_runs, run, "analysis", "all", "reactions.tsv")
        group_rnx_dict = {}
        grp_str = ""
        for group in groups_list:
            # Create string concatenating groups names
            grp_str += group + "_"
            group_list = list(self.get_grp_set(run, group))
            r = Reactions(self.path_runs, self.path_study, rxn_file, group_list)
            reactions_g = set(r.reactions_list)
            # Add reactions list of the groups in dict
            group_rnx_dict[group] = reactions_g

        if percentage:
            union = set.union(*group_rnx_dict.values())
            intersect = set.intersection(*group_rnx_dict.values())
            union_nb = len(union)
            intersect_nb = len(intersect)
            txt_file = os.path.join(self.path_study, "output_data", "compare_groups",
                                    f"{run}_intersect_{grp_str[:-1]}.txt")
            with open(txt_file, "w") as f:
                s = f"Over {union_nb} reactions, {intersect_nb} reactions in common = " \
                    f"{round((intersect_nb/union_nb) * 100, 2)} % of the union."
                print(s)
                f.write(s + "\n")
                for group, reactions in group_rnx_dict.items():
                    g_rnx_nb = len(reactions)
                    s = f"{g_rnx_nb} present among the {group} group = " \
                        f"{round((g_rnx_nb/union_nb) * 100, 2)} % of the union."
                    print(s)
                    f.write(s + "\n")
                    s = f"{g_rnx_nb - intersect_nb} only present among the {group} group = " \
                        f"{round(((g_rnx_nb - intersect_nb)/union_nb) * 100, 2)} % of the union."
                    print(s)
                    f.write(s + "\n")

        if venn_plot:
            venn.venn(group_rnx_dict, cmap="plasma")
            plt.savefig(os.path.join(self.path_study, "output_data", "compare_groups",
                                     f"{run}_intersect_{grp_str[:-1]}.png"))

    # Rename padmet id

    def rename_padmet_id(self, run):
        """ Rename the genes IDs to their original IDs in the final network PADMETs if they had been renamed by AuCoMe.
        Useful to create WIKIs.

        Parameters
        ----------
        run : str
            The run to consider to rename IDs
        """
        asso_dict = {}
        # Path to studied_organisms
        path_species = os.path.join(self.path_runs, run, "studied_organisms")
        # Path to final neworks PADMETs
        path_padmet = os.path.join(self.path_runs, run, "networks", "PADMETs")
        # Enrich asso_dict for every species having been renamed
        for species_folder in os.listdir(path_species):
            asso_dict = get_dict(run, species_folder, asso_dict, self.path_runs)
        # Make python ahocorasick automaton from asso_dict
        automaton = make_automaton(asso_dict)
        out_path = os.path.join(self.path_study, "output_data", "renamed_id_padmet")
        # Replace renamed IDs to original IDs in PADMETs files
        for species in os.listdir(path_species):
            print(f"ID renamed for {species}.padmet")
            apply_automaton(automaton, os.path.join(path_padmet, f"{species}.padmet"),
                            os.path.join(out_path, f"{species}.padmet"))
        print(f"New PADMETs saved in : {out_path}")

