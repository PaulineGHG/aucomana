"""
Pathways class
"""
from typing import Dict, List, Tuple, Set
import pandas as pd


class Pathways:
    """
    Attributes
    ----------
    name : str
        name of the file
    species_list : List[str]
        List of species studied
    data_pathways_float :
        Dataframe indicating completion rate of each pathway for each species (prop float)
    data_pathways_str :
        Dataframe indicating completion rate of each pathway for each species (prop str)
    data_reacs_assoc :
        Dataframe indicating reactions associated with each pathway for each species
    pathways_list : List[str]
        List of all pathways
    nb_pathways : int
        number of pathways
    nb_species : int
        number of species studied
    """
    STR_COMP = "_completion_rate"
    STR_RNX_ASSOC = "_rxn_assoc (sep=;)"

    def __init__(self, file_pathways_tsv: str, species_list: List[str] = None,
                 out: int = None):
        """ Init the Reactions class

        Parameters
        ----------
        file_pathways_tsv : str
            file pathways.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study.
            If not specified, will contain all the species from the run.
        """
        self.name = file_pathways_tsv.split("/")[-4]
        self.species_list = species_list
        self.data_pathways_str, self.data_reacs_assoc = self.__init_data(file_pathways_tsv)
        self.data_pathways_float = self.data_pathways_str.copy(deep=True)
        self.__convert_data_pathways()
        self.pathways_list = self.__get_filtered_pathways(out)
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
        reac_assoc_list = [x + self.STR_RNX_ASSOC for x in self.species_list]
        data_reacs_assoc = data[reac_assoc_list]
        return data_species_all_pathways, data_reacs_assoc

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
        for sp in self.data_pathways_float.columns:
            for pw in self.data_pathways_float.index:
                val_str = self.data_pathways_float.loc[pw, sp]
                if type(val_str) == str:
                    val_l = val_str.split("/")
                    self.data_pathways_float.loc[pw, sp] = int(val_l[0]) / int(val_l[1])
                else:
                    self.data_pathways_float.loc[pw, sp] = float(0)
                    nb_r = "?"
                    for v in self.data_pathways_str.loc[pw]:
                        if type(v) == str:
                            nb_r = v.split("/")[1]
                            break
                    self.data_pathways_str.loc[pw, sp] = f"0/{nb_r}"

    def __get_filtered_pathways(self, out: int) -> List[str]:
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
            return list(self.data_pathways_float.index)
        filtered_pw = []
        nb_species = len(self.species_list)
        for pw in self.data_pathways_float.index:
            count_pw = 0
            for sp in self.data_pathways_float.columns:
                if self.data_pathways_float.loc[pw, sp] != 0:
                    count_pw += 1
            if count_pw > nb_species - (out + 1):
                filtered_pw.append(pw)
        return filtered_pw

    # Loss

    # ## Min

    def is_min(self, species: str, pathway: str, unique) -> bool:
        """ Indicate if the completion of the pathway for the species is minimal (unique or not) between all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool
            True if the minimum is unique, False otherwise

        Returns
        -------
        bool
            True if the completion of the pathway for the species is minimal (unique or not) between all species,
            False otherwise
        """
        row = list(self.data_pathways_float.loc[pathway])
        min_completion = min(row)
        if unique:
            return self.data_pathways_float.loc[pathway, species] == min_completion and min_completion > 0 and \
                   row.count(min_completion) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == min_completion and min_completion > 0

    def get_pw_min(self, species: str or List[str], unique: bool = True) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of pathways having the minimal completion (unique or not),
        the pathways returned are not absent (at least 1 reaction in the pathway)

        Parameters
        ----------
        species: str or List[str]
            species or list of species to be considered
        unique: bool, optional (default=True)
            True if the minimum is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            pathways and its set having the minimal (unique or not) completion value
        """
        if type(species) == str:
            species = [species]
        min_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            min_pw = set()
            for pw in self.pathways_list:
                if self.is_min(sp, pw, unique):
                    min_pw.add(pw)
            min_pw_dict[sp[:-len(self.STR_COMP)]] = (len(min_pw), min_pw)
        return min_pw_dict

    # ## Absent

    def is_absent(self, species: str, pathway: str, unique: bool) -> bool:
        """ Indicate if the pathway is absent for the species (unique or not) : considered unique if only this species
        is not having the pathway among all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool
            True if the absence is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is absent (unique or not),
            False otherwise
        """
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return self.data_pathways_float.loc[pathway, species] == 0 and row.count(0) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == 0

    def get_pw_absent(self, species: str or List[str], unique: bool = True) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of absent pathways (unique or not) : considered unique if
        only this species is not having the pathway among all species

        Parameters
        ----------
        species: str or List[str]
            species or list of species to be considered
        unique: bool, optional (default=True)
            True if the absence is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            absent pathways and its set (unique or not)
        """
        if type(species) == str:
            species = [species]
        absent_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            absent_pw = set()
            for pw in self.pathways_list:
                if self.is_absent(sp, pw, unique):
                    absent_pw.add(pw)
            absent_pw_dict[sp[:-len(self.STR_COMP)]] = (len(absent_pw), absent_pw)
        return absent_pw_dict

    # ## Incomplete

    def is_incomplete(self, species: str, pathway: str, unique: bool) -> bool:
        """ Indicate if the pathway is incomplete for the species (unique or not) : considered unique if all other
        species have the pathway completed

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool
            True if the incomplete pathway is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is incomplete (unique or not),
            False otherwise
        """
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return sum(row) > self.nb_species - 1 and self.data_pathways_float.loc[pathway, species] < 1
        else:
            return 0 < self.data_pathways_float.loc[pathway, species] < 1

    def get_pw_incomplete(self, species: str or List[str], unique: bool = True) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of incomplete pathways (unique or not) : considered unique
        if all other species have the pathway completed

        Parameters
        ----------
        species: str or List[str]
            species or list of species to be considered
        unique: bool, optional (default=True)
            True if the incomplete pathway is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            incomplete pathways and its set (unique or not)
        """
        if type(species) == str:
            species = [species]
        incomplete_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            incomplete_pw = set()
            for pw in self.pathways_list:
                if self.is_incomplete(sp, pw, unique):
                    incomplete_pw.add(pw)
            incomplete_pw_dict[sp[:-len(self.STR_COMP)]] = (len(incomplete_pw), incomplete_pw)
        return incomplete_pw_dict

    # Gain

    # ## Present

    def is_present(self, species: str, pathway: str, unique: bool) -> bool:
        """ Indicate if the pathway is present for the species (unique or not) : considered unique if only this species
        is having the pathway among all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool
            True if the presence is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is present (unique or not),
            False otherwise
        """
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return 0 < self.data_pathways_float.loc[pathway, species] == sum(row)
        else:
            return self.data_pathways_float.loc[pathway, species] > 0

    def get_pw_present(self, species: str or List[str], unique: bool = True) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of present pathways (unique or not) : considered unique if
        only this species is having the pathway among all species

        Parameters
        ----------
        species: str or List[str]
            species or list of species to be considered
        unique: bool, optional (default=True)
            True if the presence is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            present pathways and its set (unique or not)
        """
        if type(species) == str:
            species = [species]
        present_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            present_pw = set()
            for pw in self.pathways_list:
                if self.is_present(sp, pw, unique):
                    present_pw.add(pw)
            present_pw_dict[sp[:-len(self.STR_COMP)]] = (len(present_pw), present_pw)
        return present_pw_dict

    # ## Max

    def is_max(self, species: str, pathway: str, unique: bool) -> bool:
        """ Indicate if the completion of the pathway for the species is maximal (unique or not) between all species

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool
            True if the maximum is unique, False otherwise

        Returns
        -------
        bool
            True if the completion of the pathway for the species is maximal (unique or not) between all species,
            False otherwise
        """
        row = list(self.data_pathways_float.loc[pathway])
        max_comp = max(row)
        if unique:
            return self.data_pathways_float.loc[pathway, species] == max_comp and row.count(max_comp) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == max_comp

    def get_pw_max(self, species: str or List[str], unique: bool = True) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of pathways having the maximal completion (unique or not),

        Parameters
        ----------
        species: str or List[str]
            species or list of species to be considered
        unique: bool, optional (default=True)
            True if the maximum is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            pathways and its set having the maximum (unique or not) completion value
        """
        if type(species) == str:
            species = [species]
        max_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            max_pw = set()
            for pw in self.pathways_list:
                if self.is_max(sp, pw, unique):
                    max_pw.add(pw)
            max_pw_dict[sp[:-len(self.STR_COMP)]] = (len(max_pw), max_pw)
        return max_pw_dict

    # ## Complete

    def is_complete(self, species: str, pathway: str, unique: bool) -> bool:
        """ Indicate if the pathway is complete for the species (unique or not) : considered unique if all other
        species have not the pathway completed

        Parameters
        ----------
        species: str
            species to be considered
        pathway: str
            species to be considered
        unique: bool
            True if the complete pathway is unique, False otherwise

        Returns
        -------
        bool
            True if the pathway for the species is complete (unique or not),
            False otherwise
        """
        if unique:
            row = list(self.data_pathways_float.loc[pathway])
            return self.data_pathways_float.loc[pathway, species] == 1 and row.count(1) == 1
        else:
            return self.data_pathways_float.loc[pathway, species] == 1

    def get_pw_complete(self, species: str or List[str], unique: bool = True) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of complete pathways (unique or not) : considered unique
        if all other species have not the pathway completed

        Parameters
        ----------
        species: str or List[str]
            species or list of species to be considered
        unique: bool, optional (default=True)
            True if the complete pathway is unique, False otherwise

        Returns
        -------
        min_pw: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_pathways, Set[pathways]]]) dictionary associating for each species the number of
            complete pathways and its set (unique or not)
        """
        if type(species) == str:
            species = [species]
        complete_pw_dict = {}
        for sp in species:
            sp += self.STR_COMP
            complete_pw = set()
            for pw in self.pathways_list:
                if self.is_complete(sp, pw, unique):
                    complete_pw.add(pw)
            complete_pw_dict[sp[:-len(self.STR_COMP)]] = (len(complete_pw), complete_pw)
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

    def print_completion_pw(self, pathway, species):
        species += self.STR_COMP
        comp_str = self.data_pathways_str.loc[pathway, species]
        comp_float = round(self.data_pathways_float.loc[pathway, species], 3)
        print(f"{comp_str} reactions present in pathway \"{pathway}\" = {comp_float*100}%")

