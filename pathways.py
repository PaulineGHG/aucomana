from typing import Dict, List, Tuple, Set
import pandas as pd
from algae_project import get_cat_l
from init_analysis import PATH_STUDY


class Pathways:
    """
    Attributes
    ----------
    name : str
        name of the file
    species_list : List[str]
        List of species studied
    data_pathways :
        Dataframe indicating completion rate of each pathway for each species
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
        self.data_pathways, \
            self.data_reacs_assoc = self.__init_data(file_pathways_tsv)
        self.__convert_data_pathways()
        self.pathways_list = self.__get_filtered_pathways(out)
        self.nb_pathways, self.nb_species = self.data_pathways.shape

    def __init_data(self, file_pathways_tsv: str) \
            -> Tuple['pd.DataFrame', 'pd.DataFrame']:
        """ Generate the data_reactions, data_genes_assoc and reactions_list attributes

        Parameters
        ----------
        file_pathways_tsv : str
            file pathways.tsv

        Returns
        -------
        data_pathways :
            data_pathways attribute
        data_reacs_assoc :
            data_reacs_assoc attribute
        pathways_list : List[str]
            pathways_list attribute
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
        for sp in self.data_pathways.columns:
            for pw in self.data_pathways.index:
                val_str = self.data_pathways.loc[pw, sp]
                if type(val_str) == str:
                    val_l = val_str.split("/")
                    self.data_pathways.loc[pw, sp] = int(val_l[0]) / int(val_l[1])
                else:
                    self.data_pathways.loc[pw, sp] = float(0)

    def __get_filtered_pathways(self, out: int):
        if out is None:
            return list(self.data_pathways.index)
        filtered_pw = []
        nb_species = len(self.species_list)
        for pw in self.data_pathways.index:
            count_pw = 0
            for sp in self.data_pathways.columns:
                if self.data_pathways.loc[pw, sp] != 0:
                    count_pw += 1
            if count_pw > nb_species - (out + 1):
                filtered_pw.append(pw)
        return filtered_pw

    def get_pw_lost_1_species(self, species):
        species += self.STR_COMP
        loss = set()
        for pw in self.pathways_list:
            if self.data_pathways.loc[pw, species] == 0:
                loss.add(pw)
        return len(loss), loss

    def get_pw_min_1_species(self, species):
        species += self.STR_COMP
        loss = set()
        for pw in self.pathways_list:
            min_comp = min(self.data_pathways.loc[pw])
            max_comp = max(self.data_pathways.loc[pw])
            if self.data_pathways.loc[pw, species] == min_comp and min_comp < max_comp:
                loss.add(pw)
        return len(loss), loss

    def get_pw_incomplete_1_species(self, species):
        species += self.STR_COMP
        loss = set()
        for pw in self.pathways_list:
            val = self.data_pathways.loc[pw, species]
            if sum(self.data_pathways.loc[pw]) > self.nb_species - 1 and val != 1:
                loss.add(pw)
        return len(loss), loss


R_FILE_01 = "data/reactions_data/run01b_reactions.tsv"
ORG_TSV = "data/species_group.tsv"
BROWN_01 = get_cat_l(R_FILE_01, ORG_TSV, "brown")

SLAT = "Saccharina_latissima_FEMALE"
PLAC = "Pleurocladia_lacustris"

FILE01 = "data/pathways_data/run01b_pathways.tsv"
P01 = Pathways(FILE01, BROWN_01, out=1)

# print(P01.get_pw_incomplete_1_species(SLAT))
# print(P01.get_pw_min_1_species(SLAT))
# print(P01.get_pw_lost_1_species(SLAT))
# #
# print(P01.get_pw_incomplete_1_species(PLAC))
# print(P01.get_pw_min_1_species(PLAC))
# print(P01.get_pw_lost_1_species(PLAC))

for sp in P01.species_list:
    print(sp, P01.get_pw_lost_1_species(sp))

