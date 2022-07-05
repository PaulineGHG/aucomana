"""
Pathways class
"""
import copy
from typing import Dict, List, Tuple, Set
import os
import analysis_runs.dendrograms
import pandas as pd
import numpy as np


class Pathways:
    """
    Attributes
    ----------
    path_runs: str
        path of AuCoMe runs results
    path_study: str
        path of outputs_data of the study
    name : str
        name of the file
    species_list : List[str]
        List of species studied
    data_pathways_float :
        Dataframe indicating completion rate of each pathway for each species (prop float)
    data_pathways_str :
        Dataframe indicating completion rate of each pathway for each species (prop str)
    data_rxn_assoc :
        Dataframe indicating reactions associated with each pathway for each species
    pathways_list : List[str]
        List of all pathways
    nb_pathways : int
        number of pathways
    nb_species : int
        number of species studied
    """
    nb_dend = 0
    STR_COMP = "_completion_rate"
    STR_RXN_ASSOC = "_rxn_assoc (sep=;)"

    def __init__(self, path_runs: str, path_study: str, file_pathways_tsv: str, species_list: List[str] = None,
                 out: int = None, nb_rnx_pw_min: int = 1):
        """ Init the Reactions class

        Parameters
        ----------
        path_runs: str
            path of AuCoMe runs results
        path_study: str
            path of outputs_data of the study
        file_pathways_tsv : str
            file pathways.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study.
            If not specified, will contain all the species from the run.
        out: int, optional (default=None)
            Number of species having at least a reaction in the pathway for the pathway to be kept
        nb_rnx_pw_min: int, optional (default=1)
            Minimal total number of reactions in the pathway for the pathway to be kept
        """
        self.path_runs = path_runs
        self.path_study = path_study
        self.name = file_pathways_tsv.split("/")[-4]
        self.species_list = species_list
        self.data_pathways_str, self.data_rxn_assoc = self.__init_data(file_pathways_tsv)
        self.data_pathways_float = self.data_pathways_str.copy(deep=True)
        self.nb_rnx_pw = self.__convert_data_pathways()
        self.pathways_list = self.__get_filtered_pathways(out, nb_rnx_pw_min)
        self.nb_pathways, self.nb_species = self.data_pathways_float.shape

    def __init_data(self, file_pathways_tsv: str) \
            -> Tuple['pd.DataFrame', 'pd.DataFrame']:
        """ Generate the data_reactions, data_genes_assoc and reactions_list attributes

        Parameters
        ----------
        file_pathways_tsv : str
            file pathways.tsv

        Returns
        -------
        data_pathways_str :
            data_pathways attribute
        data_reacs_assoc :
            data_reacs_assoc attribute
        """
        data = pd.read_csv(file_pathways_tsv, sep="\t", header=0, index_col='pathway')
        if self.species_list is None:
            self.__generate_species_list(data)
        comp_list = [x + self.STR_COMP for x in self.species_list]
        data_species_all_pathways = data[comp_list]
        rxn_assoc_list = [x + self.STR_RXN_ASSOC for x in self.species_list]
        data_rxn_assoc = data[rxn_assoc_list]
        return data_species_all_pathways, data_rxn_assoc

    def __generate_species_list(self, data: 'pd.DataFrame'):
        """ Generate the species_list attribute if is None

        Parameters
        ----------
        data :
            The dataframe created from pathways.tsv file
        """
        self.species_list = []
        for x in data.columns:
            if x[-7:] == '(sep=;)':
                break
            self.species_list.append(x[:-16])

    def __convert_data_pathways(self):
        """ Converts data_pathways_floats and data_pathways_str values

        data_pathways_floats values transformed in float
        data_pathways_str NA transforms in str
        """
        nb_rnx_pw = {}
        for sp in self.data_pathways_float.columns:
            for pw in self.data_pathways_float.index:
                val_str = self.data_pathways_float.loc[pw, sp]
                if type(val_str) == str:
                    val_l = val_str.split("/")
                    self.data_pathways_float.loc[pw, sp] = int(val_l[0]) / int(val_l[1])
                    nb_rnx_pw[pw] = int(val_l[1])
                else:
                    self.data_pathways_float.loc[pw, sp] = float(0)
                    nb_r = "?"
                    for v in self.data_pathways_str.loc[pw]:
                        if type(v) == str:
                            nb_r = v.split("/")[1]
                            break
                    self.data_pathways_str.loc[pw, sp] = f"0/{nb_r}"
                    if nb_r != "?":
                        nb_rnx_pw[pw] = int(nb_r)
        return nb_rnx_pw

    def __get_filtered_pathways(self, out: int, nb_rnx_pw_min: int) -> List[str]:
        """ Filter the pathways according to the number of species not having the pathway

        Parameters
        ----------
        out: int
            number of species maximum not having the pathway for the pathway to be kept

        Returns
        -------
        filtered_pw: List[str]
            List of pathway filtered
        """
        if out is None:
            out = len(self.species_list) - 1
        filtered_pw = []
        nb_species = len(self.species_list)
        for pw in self.data_pathways_float.index:
            count_pw = 0
            for sp in self.data_pathways_float.columns:
                if self.data_pathways_float.loc[pw, sp] != 0:
                    count_pw += 1
            if count_pw > nb_species - (out + 1) and self.nb_rnx_pw[pw] >= nb_rnx_pw_min:
                filtered_pw.append(pw)
        self.data_pathways_float = self.data_pathways_float.loc[filtered_pw]
        self.data_pathways_str = self.data_pathways_str.loc[filtered_pw]
        return filtered_pw

    # Loss

    # ## Min

    def is_min(self, species: str, pathway: str, unique: bool = False) -> bool:
        """ Indicate if the completion of the pathway for the species is minimal (unique or not) between all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool, optional (default=False)
            True if the minimum is unique, False otherwise

        Returns
        -------
        bool
            True if the completion of the pathway for the species is minimal (unique or not) between all species,
            False otherwise
        """
        species += self.STR_COMP
        row = list(self.data_pathways_float.loc[pathway])
        min_completion = min(row)
        if unique:
            return self.data_pathways_float.loc[pathway, species] == min_completion and min_completion > 0 and \
                   row.count(min_completion) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == min_completion and min_completion > 0

    def get_pw_min(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of pathways having the minimal completion (unique or not),
        the pathways returned are not absent (at least 1 reaction in the pathway)

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the minimum is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            pathways and its set having the minimal (unique or not) completion value
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        min_pw_dict = {}
        for sp in species:
            min_pw = set()
            for pw in self.pathways_list:
                if self.is_min(sp, pw, unique):
                    min_pw.add(pw)
            min_pw_dict[sp] = (len(min_pw), min_pw)
        return min_pw_dict

    # ## Absent

    def is_absent(self, species: str, pathway: str, unique: bool = False) -> bool:
        """ Indicate if the pathway is absent for the species (unique or not) : considered unique if only this species
        is not having the pathway among all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool, optional (default=False)
            True if the absence is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is absent (unique or not),
            False otherwise
        """
        species += self.STR_COMP
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return self.data_pathways_float.loc[pathway, species] == 0 and row.count(0) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == 0

    def get_pw_absent(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of absent pathways (unique or not) : considered unique if
        only this species is not having the pathway among all species

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the absence is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            absent pathways and its set (unique or not)
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        absent_pw_dict = {}
        for sp in species:
            absent_pw = set()
            for pw in self.pathways_list:
                if self.is_absent(sp, pw, unique):
                    absent_pw.add(pw)
            absent_pw_dict[sp] = (len(absent_pw), absent_pw)
        return absent_pw_dict

    # ## Incomplete

    def is_incomplete(self, species: str, pathway: str, unique: bool = False) -> bool:
        """ Indicate if the pathway is incomplete for the species (unique or not) : considered unique if all other
        species have the pathway completed

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool, optional (default=False)
            True if the incomplete pathway is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is incomplete (unique or not),
            False otherwise
        """
        species += self.STR_COMP
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return sum(row) > self.nb_species - 1 and self.data_pathways_float.loc[pathway, species] < 1
        else:
            return 0 < self.data_pathways_float.loc[pathway, species] < 1

    def get_pw_incomplete(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of incomplete pathways (unique or not) : considered unique
        if all other species have the pathway completed

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the incomplete pathway is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            incomplete pathways and its set (unique or not)
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        incomplete_pw_dict = {}
        for sp in species:
            incomplete_pw = set()
            for pw in self.pathways_list:
                if self.is_incomplete(sp, pw, unique):
                    incomplete_pw.add(pw)
            incomplete_pw_dict[sp] = (len(incomplete_pw), incomplete_pw)
        return incomplete_pw_dict

    # Gain

    # ## Present

    def is_present(self, species: str, pathway: str, unique: bool = False) -> bool:
        """ Indicate if the pathway is present for the species (unique or not) : considered unique if only this species
        is having the pathway among all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            pathway to be considered
        unique: bool, optional (default=False)
            True if the presence is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is present (unique or not),
            False otherwise
        """
        species += self.STR_COMP
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return 0 < self.data_pathways_float.loc[pathway, species] == sum(row)
        else:
            return self.data_pathways_float.loc[pathway, species] > 0

    def get_pw_present(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of present pathways (unique or not) : considered unique if
        only this species is having the pathway among all species

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the presence is unique, False otherwise

        Returns
        -------
        present_pw_dict: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            present pathways and its set (unique or not)
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        present_pw_dict = {}
        for sp in species:
            present_pw = set()
            for pw in self.pathways_list:
                if self.is_present(sp, pw, unique):
                    present_pw.add(pw)
            present_pw_dict[sp] = (len(present_pw), present_pw)
        return present_pw_dict

    # ## Max

    def is_max(self, species: str, pathway: str, unique: bool = False) -> bool:
        """ Indicate if the completion of the pathway for the species is maximal (unique or not) between all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool, optional (default=False)
            True if the maximum is unique, False otherwise

        Returns
        -------
        bool
            True if the completion of the pathway for the species is maximal (unique or not) between all species,
            False otherwise
        """
        species += self.STR_COMP
        row = list(self.data_pathways_float.loc[pathway])
        max_comp = max(row)
        if unique:
            return self.data_pathways_float.loc[pathway, species] == max_comp and row.count(max_comp) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == max_comp

    def get_pw_max(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of pathways having the maximal completion (unique or not),

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the maximum is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            pathways and its set having the maximum (unique or not) completion value
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        max_pw_dict = {}
        for sp in species:
            max_pw = set()
            for pw in self.pathways_list:
                if self.is_max(sp, pw, unique):
                    max_pw.add(pw)
            max_pw_dict[sp] = (len(max_pw), max_pw)
        return max_pw_dict

    # ## Complete

    def is_complete(self, species: str, pathway: str, unique: bool = False) -> bool:
        """ Indicate if the pathway is complete for the species (unique or not) : considered unique if all other
        species have not the pathway completed

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool, optional (default=False)
            True if the complete pathway is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is complete (unique or not),
            False otherwise
        """
        species += self.STR_COMP
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return self.data_pathways_float.loc[pathway, species] == 1 and row.count(1) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == 1

    def get_pw_complete(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of complete pathways (unique or not) : considered unique
        if all other species have not the pathway completed

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the complete pathway is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            complete pathways and its set (unique or not)
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        complete_pw_dict = {}
        for sp in species:
            complete_pw = set()
            for pw in self.pathways_list:
                if self.is_complete(sp, pw, unique):
                    complete_pw.add(pw)
            complete_pw_dict[sp] = (len(complete_pw), complete_pw)
        return complete_pw_dict

    def get_pw_over_treshold(self, species: str or List[str], threshold: float, strict: bool = False) \
            -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns the pathways with completion over a given threshold for species

        Parameters
        ----------
        species: str or List[str]
            Species or list of species to be considered
        threshold : float
            Threshold of completion
        strict : bool, optional (default=False)
            Whether the inequality relative to the threshold should be strict or not

        Returns
        -------
        over_t_pw_dict : Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            pathways over the threshold and its set
        """
        if threshold < 0 or threshold > 1:
            raise ValueError("The threshold must be within the range of 0 to 1")
        if type(species) == str:
            species = [species]
        over_t_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            over_t_pw = set()
            for pw in self.pathways_list:
                if strict:
                    if self.data_pathways_float.loc[pw, sp] > threshold:
                        over_t_pw.add(pw)
                else:
                    if self.data_pathways_float.loc[pw, sp] >= threshold:
                        over_t_pw.add(pw)
            over_t_pw_dict[sp[:-len(self.STR_COMP)]] = (len(over_t_pw), over_t_pw)
        return over_t_pw_dict

    def get_pw_under_treshold(self, species: str or List[str], threshold: float, strict: bool = False) \
            -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns the pathways with completion under a given threshold for species

        Parameters
        ----------
        species: str or List[str]
            Species or list of species to be considered
        threshold : float
            Threshold of completion
        strict : bool, optional (default=False)
            Whether the inequality relative to the threshold should be strict or not

        Returns
        -------
        over_t_pw_dict : Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            pathways under the threshold and its set
        """
        if strict:
            over_t_dict = self.get_pw_over_treshold(species, threshold, False)
        else:
            over_t_dict = self.get_pw_over_treshold(species, threshold, True)
        under_t_pw_dict = {}
        for sp in over_t_dict:
            under_t_pw_dict[sp] = (self.nb_pathways - over_t_dict[sp][0],
                                   set(self.pathways_list).difference(over_t_dict[sp][1]))
        return under_t_pw_dict

    # Other

    def get_pw_names(self):
        """ Associate for each pathway in pathways_list, its common name in a dictionary.

        Returns
        -------
        pw_names : Dict[str, str]
            Dictionary Dict[ID of pw, common name of pw] associating common pathways names to their ID.
        """
        pw_set = set(self.pathways_list)
        pw_names = {}
        padmet_file = os.path.join(self.path_runs, self.name, "analysis", "all",
                                   "all_panmetabolism.padmet")
        with open(padmet_file, "r") as f:
            for l in f:
                if "COMMON-NAME" in l and "pathway" in l:
                    l = l.split("\t")
                    if l[1] in pw_set:
                        pw_names[l[1]] = l[3]
        return pw_names

    def convert_df_to_binary(self, threshold: float, strict: bool = False,
                             output_file: bool = False, common_name: bool = False,
                             nb_rxn: bool = False) -> 'pd.DataFrame':
        """ Returns the pathways with completion over a given threshold for species

        Parameters
        ----------
        threshold : float
            Threshold of completion
        strict : bool, optional (default=False)
            Whether the inequality relative to the threshold should be strict or not
        output_file : bool, optional (default=False)
            Whether write the binary data frame in an output file or not
        common_name : bool, optional (default=False)
            Whether write the common name of the pathway in 1st column or not
        nb_rxn : bool, optional (default=False)
            Whether write the total number of reactions of the pathway in 2nd column or not

        Returns
        -------
        df_binary : 'pd.DataFrame'
            Dataframe with binary values according to the threshold indicated :
            value = 1 if value >= (or > if strict) threshold
            value = 0 if value < (or <= if strict) threshold
        """
        if threshold < 0 or threshold > 1:
            raise ValueError("The threshold must be within the range of 0 to 1")
        if common_name:
            common_pw_name = self.get_pw_names()
            if nb_rxn:
                columns = ["common name"] + self.species_list + ["nb rnx pw"]
            else:
                columns = ["common name"] + self.species_list
        else:
            if nb_rxn:
                columns = self.species_list + ["nb rnx pw"]
            else:
                columns = self.species_list
        df_binary = pd.DataFrame(columns=columns, index=self.pathways_list)
        for pw in self.pathways_list:
            if nb_rxn:
                df_binary.loc[pw, "nb rnx pw"] = int(self.nb_rnx_pw[pw])
            if common_name:
                df_binary.loc[pw, "common name"] = common_pw_name[pw]
            for sp in self.species_list:
                sp_c = sp + self.STR_COMP
                if strict:
                    if self.data_pathways_float.loc[pw, sp_c] > threshold:
                        df_binary.loc[pw, sp] = int(1)
                    else:
                        df_binary.loc[pw, sp] = int(0)
                else:
                    if self.data_pathways_float.loc[pw, sp_c] >= threshold:
                        df_binary.loc[pw, sp] = int(1)
                    else:
                        df_binary.loc[pw, sp] = int(0)

        if output_file:
            file_name = f"{self.name}_{threshold}_binary_pw.tsv"
            file_path = os.path.join(self.path_study, "output_data", "pathways_data", "binary_df",
                                     file_name)
            df_binary.to_csv(file_path, sep="\t", index_label="pathway")
        else:
            return df_binary

    def generate_pw_dendrogram(self, a: 'Analysis', df_binary_theshold: float, strict: bool = False,
                               name: str = None, phylo_file: str = None, n_boot: int = 100000):
        if name is None:
            self.nb_dend += 1
            name = f"dendrogram{self.nb_dend}"
        name = f"pw_{int(df_binary_theshold*100)}_{name}"
        df_bin = self.convert_df_to_binary(df_binary_theshold, strict)
        d = analysis_runs.dendrograms.Dendrogram(a, self.path_study, df_bin, self.name, name,
                                                 phylo_file)
        d.get_dendro_pvclust(n_boot)

    def print_completion_pw(self, pathway, species):
        species += self.STR_COMP
        comp_str = self.data_pathways_str.loc[pathway, species]
        comp_float = round(self.data_pathways_float.loc[pathway, species], 3)
        print(f"{comp_str} reactions present in pathway \"{pathway}\" = {comp_float*100}%")

