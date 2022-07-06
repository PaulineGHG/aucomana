"""
Metabolites class
"""
from typing import Dict, List, Tuple, Set
import analysis_runs.dendrograms
import pandas as pd
import os


class Metabolites:
    """
    Attributes
    ----------
    path_runs: str
        path of AuCoMe runs results
    path_study: str
        path of outputs_data of the study
    name : str
        name of the run
    species_list : List[str]
        List of species studied
    data_metabolites_consumed : pd.DataFrame
        DataFrame indicating which reaction consumes the metabolite for each species
    data_metabolites_produced : pd.DataFrame
        DataFrame indicating which reaction produces the metabolite for each species
    data_metabolites : pd.DataFrame
        DataFrame indicating which metabolite is produced and/or consumed (1 or 0) by a reaction for each species
    metabolites_list : List[str]
        List of all metabolites
    nb_metabolites : int
        number of metabolites
    nb_species : int
        number of species studied
    """
    STR_CONSUME = "_rxn_consume"
    STR_PRODUCE = "_rxn_produce"
    nb_dend = 0

    def __init__(self, path_runs: str, path_study: str, file_metabolites_tsv: str, species_list: List[str] = None):
        """ Init the Genes class

        Parameters
        ----------
        path_runs: str
            path of AuCoMe runs results
        path_study: str
            path of outputs_data of the study
        file_metabolites_tsv : str
            file metabolites.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study (must correspond to their name in metabolites.tsv file).
            If not specified, will contain all the species from metabolites.tsv file.
        """
        self.path_runs = path_runs
        self.path_study = path_study
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
        data_metabolites_consumed : pd.DataFrame
            data_metabolites_consumed attribute
        data_metabolites_produced : pd.DataFrame
            data_metabolites_produced attribute
        metabolites_list : List[str]
            metabolites_list attribute
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
        data : pd.DataFrame
            The dataframe created from metabolites.tsv file
        """
        self.species_list = []
        for x in data.columns:
            if x[-12:] == self.STR_PRODUCE:
                break
            self.species_list.append(x[:-12])

    def __generate_metabolites_df(self):
        """ Generate the metabolites_data attribute. The dataframe indicate if a metabolite is consumed and/or produced
        by a reaction or not (1 or 0) for each species.

        Returns
        -------
        metabolite_df : pd.DataFrame
            data_metabolites attribute
        """
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
        d = analysis_runs.dendrograms.Dendrogram(self.path_runs, self.path_study, self.data_metabolites, name,
                                                 self.name, phylo_file)
        d.get_dendro_pvclust(n_boot)

    def get_rxn_consuming(self, metabolites_list: str or List[str] = None) -> \
            Dict[str, Dict[str, Dict[str, List[str]]]]:
        """ Returns a dictionary of reactions consuming each metabolite for each species.

        Parameters
        ----------
        metabolites_list : str or List[str], optional (default=None)
            List of metabolites to find reactions consuming them.
            If None, will be metabolites_list attribute.

        Returns
        -------
        rxn_consuming : Dict[str, Dict[str, Dict[str, List[str]]]]
            Dictionary of reactions consuming each metabolite for each species
        """
        rxn_consuming = {}
        if metabolites_list is None:
            metabolites_list = self.metabolites_list
        elif type(metabolites_list) == str:
            metabolites_list = [metabolites_list]
        for mb in metabolites_list:
            rxn_consuming[mb] = {}
        for mb in metabolites_list:
            for species in self.species_list:
                rxn_consuming[mb][species] = str(
                    self.data_metabolites_consumed[species + self.STR_CONSUME][mb]).split(";")
        return rxn_consuming

    def get_rxn_producing(self, metabolites_list: str or List[str] = None) -> \
            Dict[str, Dict[str, Dict[str, List[str]]]]:
        """ Returns a dictionary of reactions producing each metabolite for each species.

        Parameters
        ----------
        metabolites_list : str or List[str], optional (default=None)
            List of metabolites to find reactions producing them.
            If None, will be metabolites_list attribute.

        Returns
        -------
        rxn_producing : Dict[str, Dict[str, Dict[str, List[str]]]]
            Dictionary of reactions producing each metabolite for each species
        """
        rxn_producing = {}
        if metabolites_list is None:
            metabolites_list = self.metabolites_list
        elif type(metabolites_list) == str:
            metabolites_list = [metabolites_list]
        for mb in metabolites_list:
            rxn_producing[mb] = {}
        for mb in metabolites_list:
            for species in self.species_list:
                rxn_producing[mb][species] = str(
                    self.data_metabolites_produced[species + self.STR_PRODUCE][mb]).split(";")
        return rxn_producing

    def get_metabolites_names(self):
        """ Associate for each metabolite in metabolites_list, its common name in a dictionary.

        Returns
        -------
        mb_names : Dict[str, str]
            Dictionary Dict[ID of mb, common name of mb] associating common metabolites names to their ID.
        """
        mb_set = set(self.metabolites_list)
        mb_names = {}
        padmet_file = os.path.join(self.path_runs, self.name, "analysis", "all",
                                   "all_panmetabolism.padmet")
        with open(padmet_file, "r") as f:
            for l in f:
                if "COMMON-NAME" in l and "compound" in l:
                    l = l.split("\t")
                    if l[1] in mb_set:
                        mb_names[l[1]] = l[3]
        return mb_names


