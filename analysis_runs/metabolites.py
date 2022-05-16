"""
Metabolites class
"""
from typing import Dict, List, Tuple, Set
from analysis_runs.init_analysis import PATH_STUDY, PATH_RUNS
import analysis_runs.dendrograms
import pandas as pd
import json
import os


class Metabolites:
    """
    Attributes
    ----------
    name : str
        name of the run
    species_list : List[str]
        List of species studied
    data_metabolites_consumed :
        Dataframe indicating wich for each species
    data_metabolites_produced :
        Dataframe indicating reactions associated with each gene for each species
    metabolites_list : List[str]
        List of all genes
    nb_metabolites : int
        number of genes
    nb_species : int
        number of species studied
    """
    STR_CONSUME = "_rxn_consume"
    STR_PRODUCE = "_rxn_produce"
    nb_dend = 0

    def __init__(self, file_metabolites_tsv: str, species_list: List[str] = None):
        """ Init the Genes class

        Parameters
        ----------
        file_metabolites_tsv : str
            file metabolites.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study (must correspond to their name in metabolites.tsv file).
            If not specified, will contain all the species from metabolites.tsv file.
        """
        self.name = file_metabolites_tsv.split("/")[-4]
        self.species_list = species_list
        self.data_metabolites_consumed, \
            self.data_metabolites_produced, \
            self.metabolites_list = self.__init_data(file_metabolites_tsv)
        self.data_metabolites = self.__generate_metabolites_df()
        self.nb_metabolites, self.nb_species = self.data_metabolites_consumed.shape

    def __init_data(self, file_metabolites_tsv: str) \
            -> Tuple['pd.DataFrame', 'pd.DataFrame', List[str]]:
        """ Generate the data_metabolites_consumed, data_metabolites_produced and metabolites_list
         attributes

        Parameters
        ----------
        file_metabolites_tsv : str
            file metabolites.tsv

        Returns
        -------
        data_genes :
            data_genes attribute
        data_rnx_assoc :
            data_rnx_assoc attribute
        genes_list : List[str]
            genes_list
        """
        data = pd.read_csv(file_metabolites_tsv, sep="\t", header=0, index_col='metabolite')
        if self.species_list is None:
            self.__generate_species_list(data)
        rnx_consume_list = [x + self.STR_CONSUME for x in self.species_list]
        rnx_produce_list = [x + self.STR_PRODUCE for x in self.species_list]
        data_consume_all_metabolites = data[rnx_consume_list]
        data_produce_all_metabolites = data[rnx_produce_list]
        data_consume_all_metabolites = data_consume_all_metabolites.fillna(int(0))
        data_produce_all_metabolites = data_produce_all_metabolites.fillna(int(0))
        metabolites_list = list(data.index)
        return data_consume_all_metabolites.loc[metabolites_list], \
            data_produce_all_metabolites.loc[metabolites_list], metabolites_list

    def __generate_species_list(self, data: 'pd.DataFrame'):
        """ Generate the species_list attribute if is None

        Parameters
        ----------
        data :
            The dataframe created from genes.tsv file
        """
        self.species_list = []
        for x in data.columns:
            if x[-12:] == self.STR_PRODUCE:
                break
            self.species_list.append(x[:-12])

    def __generate_metabolites_df(self):
        metabolites_df = pd.DataFrame(index=self.metabolites_list, columns=self.species_list)
        metabolites_df.index_label = "metabolite"
        for met in self.metabolites_list:
            for sp in self.species_list:
                sp_prod = sp + self.STR_PRODUCE
                sp_cons = sp + self.STR_CONSUME
                if self.data_metabolites_consumed.loc[met, sp_cons] == 0 and \
                        self.data_metabolites_produced.loc[met, sp_prod] == 0:
                    metabolites_df.loc[met, sp] = 0
                else:
                    metabolites_df.loc[met, sp] = 1
        return metabolites_df

    def generate_met_dendrogram(self, name=None, phylo_file=None, n_boot=100000):
        if name is None:
            self.nb_dend += 1
            name = f"dendrogram{self.nb_dend}"
        name = "met_" + name
        analysis_runs.dendrograms.get_dendro_pvclust(self.data_metabolites, name, self.name,
                                                     phylo_file, n_boot)


