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
    data_genes :
        Dataframe indicating if genes are present for each species
    data_rnx_assoc :
        Dataframe indicating reactions associated with each gene for each species
    genes_list : List[str]
        List of all genes
    nb_genes : int
        number of genes
    nb_species : int
        number of species studied
    """
    STR_RNX_ASSOC = "_rxn_assoc (sep=;)"

    def __init__(self, file_genes_tsv: str, species_list: List[str] = None):
        """ Init the Genes class

        Parameters
        ----------
        file_genes_tsv : str
            file genes.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study (must correspond to their name in genes.tsv file).
            If not specified, will contain all the species from genes.tsv file.
        """
        self.name = file_genes_tsv.split("/")[-4]
        self.species_list = species_list
        self.data_genes, \
            self.data_rnx_assoc, \
            self.genes_list = self.__init_data(file_genes_tsv)
        self.nb_genes, self.nb_species = self.data_genes.shape

    def __init_data(self, file_genes_tsv: str) \
            -> Tuple['pd.DataFrame', 'pd.DataFrame', List[str]]:
        """ Generate the data_genes, data_rnx_assoc and genes_list attributes

        Parameters
        ----------
        file_genes_tsv : str
            file genes.tsv

        Returns
        -------
        data_genes :
            data_genes attribute
        data_rnx_assoc :
            data_rnx_assoc attribute
        genes_list : List[str]
            genes_list
        """
        data = pd.read_csv(file_genes_tsv, sep="\t", header=0, index_col='gene')
        if self.species_list is None:
            self.__generate_species_list(data)
        data_species_all_genes = data[self.species_list]
        rnx_assoc_list = [x + self.STR_RNX_ASSOC for x in self.species_list]
        data_rnx_assoc = data[rnx_assoc_list]
        genes_list = list(data_species_all_genes.index)
        return data_species_all_genes.loc[genes_list], \
            data_rnx_assoc.loc[genes_list], genes_list

    def __generate_species_list(self, data: 'pd.DataFrame'):
        """ Generate the species_list attribute if is None

        Parameters
        ----------
        data :
            The dataframe created from genes.tsv file
        """
        self.species_list = []
        for x in data.columns:
            if x[-7:] == '(sep=;)':
                break
            self.species_list.append(x)