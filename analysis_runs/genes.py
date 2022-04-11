"""
Genes class
"""
from typing import Dict, List, Tuple, Set
from analysis_runs.init_analysis import PATH_STUDY, PATH_RUNS
import pandas as pd
import json
import os


class Genes:
    """
    Attributes
    ----------
    name : str
        name of the file
    species_list : List[str]
        List of species studied
    data_reactions :
        Dataframe indicating if filtered reactions are present in each species
    data_genes_assoc :
        Dataframe indicating genes associated with each filtered reaction for each species
    reactions_list : List[str]
        List of all reactions filtered according to the number of species not having the reaction
    nb_reactions : int
        number of reactions (filtered)
    nb_species : int
        number of species studied
    reactions_loss : Dict[str, Tuple[int, Set[str]]]
        Dictionary indicating for each species how many and which reactions are lost
    """
    nb_common_reac = 0
    nb_genes_assoc = 0
    STR_GENE_ASSOC = "_genes_assoc (sep=;)"

    def __init__(self, file_reactions_tsv: str, species_list: List[str] = None, out: int = None):
        """ Init the Reactions class

        Parameters
        ----------
        file_reactions_tsv : str
            file reactions.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study (must correspond to their name in reactions.tsv file).
            If not specified, will contain all the species from reactions.tsv file.
        out : int, optional (default = 1)
            number of species maximum not having the reaction for the reaction to be kept
        """
        self.name = file_reactions_tsv.split("/")[-4]
        self.species_list = species_list
        self.data_reactions, \
            self.data_genes_assoc, \
            self.reactions_list = self.__init_data(file_reactions_tsv, out)
        self.nb_reactions, self.nb_species = self.data_reactions.shape

    def __init_data(self, file_reactions_tsv: str, out: int) \
            -> Tuple['pd.DataFrame', 'pd.DataFrame', List[str]]:
        """ Generate the data_reactions, data_genes_assoc and reactions_list attributes

        Parameters
        ----------
        file_reactions_tsv : str
            file reactions.tsv
        out : int
            number of species maximum not having the reaction for the reaction to be kept

        Returns
        -------
        data_reactions :
            data_reactions attribute
        data_genes_assoc :
            data_genes_assoc attribute
        reactions_list : List[str]
            reactions_list
        """
        data = pd.read_csv(file_reactions_tsv, sep="\t", header=0, index_col='reaction')
        if self.species_list is None:
            self.__generate_species_list(data)
        data_species_all_reactions = data[self.species_list]
        genes_assoc_list = [x + self.STR_GENE_ASSOC for x in self.species_list]
        data_genes_assoc = data[genes_assoc_list]
        filtered_reactions = self.__get_filtered_reactions(data_species_all_reactions, out)
        return data_species_all_reactions.loc[filtered_reactions], \
            data_genes_assoc.loc[filtered_reactions], filtered_reactions

    def __generate_species_list(self, data: 'pd.DataFrame'):
        """ Generate the species_list attribute if is None

        Parameters
        ----------
        data :
            The dataframe created from reactions.tsv file
        """
        self.species_list = []
        for x in data.columns:
            if x[-7:] == '(sep=;)':
                break
            self.species_list.append(x)