"""
Genes class
"""
from typing import Dict, List, Tuple, Set
import pandas as pd


class Genes:
    """
    Attributes
    ----------
    species_list : List[str]
        List of species studied
    data_genes : DataFrame
        The Dataframe indicating if genes are present for each species
    data_rxn_assoc : DataFrame
        The Dataframe indicating reactions associated with each gene for each species
    genes_list : List[str]
        List of all genes
    nb_genes : int
        number of genes
    nb_species : int
        number of species studied
    """
    STR_RXN_ASSOC = "_rxn_assoc (sep=;)"

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
        self.species_list = species_list
        self.data_genes, \
            self.data_rxn_assoc, \
            self.genes_list = self.__init_data(file_genes_tsv)
        self.nb_genes, self.nb_species = self.data_genes.shape
        self.nb_genes_species = self.__generate_nb_genes_sp()

    def __init_data(self, file_genes_tsv: str) \
            -> Tuple['pd.DataFrame', 'pd.DataFrame', List[str]]:
        """ Generate the data_genes, data_rxn_assoc and genes_list attributes

        Parameters
        ----------
        file_genes_tsv : str
            file genes.tsv

        Returns
        -------
        data_genes : pd.DataFrame
            data_genes attribute
        data_rxn_assoc : pd.DataFrame
            data_rxn_assoc attribute
        genes_list : List[str]
            genes_list
        """
        data = pd.read_csv(file_genes_tsv, sep="\t", header=0, index_col='gene')
        if self.species_list is None:
            self.__generate_species_list(data)
        data_species_all_genes = data[self.species_list]
        data_species_all_genes = data_species_all_genes.fillna(int(0))
        rnx_assoc_list = [x + self.STR_RXN_ASSOC for x in self.species_list]
        data_rnx_assoc = data[rnx_assoc_list]
        genes_list = list(data_species_all_genes.index)
        return data_species_all_genes.loc[genes_list], \
            data_rnx_assoc.loc[genes_list], genes_list

    def __generate_species_list(self, data: 'pd.DataFrame'):
        """ Generate the species_list attribute if is None

        Parameters
        ----------
        data : pd.DataFrame
            The dataframe created from genes.tsv file
        """
        self.species_list = []
        for x in data.columns:
            if x[-7:] == '(sep=;)':
                break
            self.species_list.append(x)

    def __generate_nb_genes_sp(self) -> Dict[str, int]:
        """ Generate the nb_genes_sp attribute. Returns a dictionary associating for each species the number of genes
        it has.

        Returns
        -------
        nb_genes_sp_dict : Dict[str, int]
            Dictionary Dict[species name, number of genes] associating for each species the number of genes it has.
        """
        nb_genes_sp_dict = {}
        for sp in self.species_list:
            nb_genes_sp_dict[sp] = int(sum(self.data_genes[sp]))
        return nb_genes_sp_dict

    def __get_genes_1sp(self, species: str) -> Set[str]:
        """ Returns for the species to consider, its set of genes.

        Parameters
        ----------
        species : str
            Name of the species to consider

        Returns
        -------
        genes_set : Set[str]
            Set of genes the species to consider have
        """
        genes_set = set()
        for gene in self.genes_list:
            if self.data_genes.loc[gene, species] == 1:
                genes_set.add(gene)
        return genes_set

    def get_genes_species(self, species_list: str or List[str] = None) -> Dict[str, Set[str]]:
        """ Returns a dictionary associating for each species, its set of genes it has.

        Parameters
        ----------
        species_list : str or List[str], optional (default=None)
            species list to consider, if None will be the species_list attribute (all the species)

        Returns
        -------
        genes_sp_dict : Dict[str, Set[str]]
            Dictionary Dict[species name, Set[genes id]] associating for each species, its set of genes it has.
        """
        if type(species_list) == str:
            species_list = [species_list]
        elif species_list is None:
            species_list = self.species_list
        genes_sp_dict = {}
        for sp in species_list:
            genes_set = self.__get_genes_1sp(sp)
            genes_sp_dict[sp] = genes_set
        return genes_sp_dict
