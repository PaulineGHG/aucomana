"""
Metabolites class
"""
from typing import Dict, List, Tuple
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
    reactions_consuming : Dict[str, Dict[str, List[str]]]
        Dictionary of reactions consuming each metabolite for each species
    reactions_producing : Dict[str, Dict[str, List[str]]]
        Dictionary of reactions producing each metabolite for each species
    metabolites_produced : pd.DataFrame.
        DataFrame indicating which metabolite is produced (1 or 0) by a reaction for each species
    metabolites_consumed : pd.DataFrame.
        DataFrame indicating which metabolite is consumed (1 or 0) by a reaction for each species
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
        self.reactions_consuming, \
            self.reactions_producing, \
            self.metabolites_produced, \
            self.metabolites_consumed, \
            self.metabolites_list = self.__init_data(file_metabolites_tsv)
        self.nb_metabolites = len(self.metabolites_list)
        self.nb_species = len(self.species_list)

    def __init_data(self, file_metabolites_tsv: str) \
            -> Tuple[Dict[str, Dict[str, List[str]]], Dict[str, Dict[str, List[str]]],
                     pd.DataFrame, pd.DataFrame,
                     List[str]]:
        """ Generate the data_metabolites_consumed, data_metabolites_produced and metabolites_list
         attributes

        Parameters
        ----------
        file_metabolites_tsv : str
            file metabolites.tsv

        Returns
        -------
        reactions_consuming : Dict[str, Dict[str, List[str]]]
            reactions_consuming attribute
        reactions_producing : Dict[str, Dict[str, List[str]]]
            reactions_producing attribute
        metabolites_produced : pd.DataFrame
            metabolites_produced attribute
        metabolites_consumed : pd.DataFrame
            metabolites_consumed attribute
        metabolites_list : List[str]
            metabolites_list attribute
        """
        data = pd.read_csv(file_metabolites_tsv, sep="\t", header=0, index_col='metabolite')
        if self.species_list is None:
            self.__generate_species_list(data)
        rnx_consume_list = [x + self.STR_CONSUME for x in self.species_list]
        rnx_produce_list = [x + self.STR_PRODUCE for x in self.species_list]
        metabolites_list = list(data.index)
        reactions_consuming = self.__get_rxn(data[rnx_consume_list], metabolites_list, True)
        reactions_producing = self.__get_rxn(data[rnx_produce_list], metabolites_list, False)

        metabolites_produced, metabolites_consumed = self.__generate_metabolites_df(data)

        return reactions_consuming, reactions_producing, metabolites_produced, metabolites_consumed, metabolites_list

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

    def __generate_metabolites_df(self, metabolites_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """ Generate the metabolites_data attribute. The dataframe indicate if a metabolite is consumed and/or produced
        by a reaction or not (1 or 0) for each species.

        Returns
        -------
        metabolites_produced : pd.DataFrame
            metabolites_produced attribute
        metabolites_consumed : pd.DataFrame
            metabolites_consumed attribute
        """
        metabolites_df = metabolites_df.fillna(int(0))
        metabolites_df[metabolites_df != 0] = int(1)
        metabolites_produced = metabolites_df[[sp + self.STR_PRODUCE for sp in self.species_list]]
        metabolites_consumed = metabolites_df[[sp + self.STR_CONSUME for sp in self.species_list]]
        metabolites_produced.columns = metabolites_produced.columns.str.replace(self.STR_PRODUCE, '')
        metabolites_consumed.columns = metabolites_consumed.columns.str.replace(self.STR_CONSUME, '')
        return metabolites_produced, metabolites_consumed

    def __get_rxn(self, df: pd.DataFrame, metabolites_list: List[str], consumed: bool) -> \
            Dict[str, Dict[str, List[str]]]:
        """ Returns a dictionary of reactions consuming or producing each metabolite for each species.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame of data from metabolites.tsv file
        metabolites_list : List[str]
            List of metabolites to find reactions consuming or producing them.
        consumed : bool
            If True will return the reactions consuming the metabolites
            If False will return the reactions producing the metabolites

        Returns
        -------
        rxn_dict : Dict[str, Dict[str, Dict[str, List[str]]]]
            Dictionary of reactions consuming or producing each metabolite for each species
        """
        rxn_dict = {}
        for mb in metabolites_list:
            rxn_dict[mb] = {}
            for species in self.species_list:
                if consumed:
                    species_id = species + self.STR_CONSUME
                else:
                    species_id = species + self.STR_PRODUCE
                if type(df[species_id][mb]) == str:
                    rxn_dict[mb][species] = str(df[species_id][mb]).split(";")
        return rxn_dict

    def generate_met_dendrogram(self, name=None, phylo_file=None, n_boot=10000):
        if name is None:
            self.nb_dend += 1
            name = f"dendrogram{self.nb_dend}"
        name = "met_" + name
        d = analysis_runs.dendrograms.Dendrogram(self.path_runs, self.path_study, self.metabolites_produced, self.name,
                                                 name, phylo_file)
        d.get_dendro_pvclust(n_boot)

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

    def is_produced(self, species: str, metabolite: str, unique=False) -> bool:
        """ Indicate if the metabolite is produced for the species (unique or not) : considered unique if only this
        species is producing the metabolite among all species

        Parameters
        ----------
        species: str
            species to be considered
        metabolite: str
            metabolite to be considered
        unique: bool, optional (default=False)
            True if the production is unique, False otherwise

        Returns
        -------
        bool
            True if the metabolite for the species is produced (unique or not),
            False otherwise
        """
        if metabolite not in self.metabolites_list:
            print(f'Metabolite {metabolite} not in the run results')
            return False
        else:
            if unique:
                row = list(self.metabolites_produced.loc[metabolite])
                return self.metabolites_produced.loc[metabolite, species] == 1 and sum(row) == 1
            else:
                return self.metabolites_produced.loc[metabolite, species] == 1

    def is_consumed(self, species: str, metabolite: str, unique=False) -> bool:
        """ Indicate if the metabolite is consumed for the species (unique or not) : considered unique if only this
        species is consuming the metabolite among all species

        Parameters
        ----------
        species: str
            species to be considered
        metabolite: str
            metabolite to be considered
        unique: bool, optional (default=False)
            True if the consumption is unique, False otherwise

        Returns
        -------
        bool
            True if the metabolite for the species is consumed (unique or not),
            False otherwise
        """
        if metabolite not in self.metabolites_list:
            print(f'Metabolite {metabolite} not in the run results')
            return False
        else:
            if unique:
                row = list(self.metabolites_consumed.loc[metabolite])
                return self.metabolites_consumed.loc[metabolite, species] == 1 and sum(row) == 1
            else:
                return self.metabolites_consumed.loc[metabolite, species] == 1
