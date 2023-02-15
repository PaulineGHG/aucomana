"""
Reactions class
"""
from typing import Dict, List, Tuple, Set
import pandas as pd


class Reactions:
    """
    Attributes
    ----------
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
    """
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
        out : int, optional (default = None)
            number of species maximum not having the reaction for the reaction to be kept
        """
        self.species_list = species_list
        self.data_reactions, \
            self.data_genes_assoc, \
            self.reactions_list = self.__init_data(file_reactions_tsv, out)
        self.nb_reactions, self.nb_species = self.data_reactions.shape
        self.nb_rxn_sp = self.__generate_nb_rxn_sp()

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

    def __get_filtered_reactions(self, data_all_reactions: 'pd.DataFrame', out: int) \
            -> List[str]:
        """ Filter the reactions according to the number of species not having the reaction

        Parameters
        ----------
        data_all_reactions:
            Dataframe with filtered columns and unfiltered reactions (rows)
        out: int
            number of species maximum not having the reaction for the reaction to be kept

        Returns
        -------
        filtered_reactions : List[str]
            List of reactions filtered
        """
        nb_species = len(self.species_list)
        if out is None:
            out = nb_species - 1
        filtered_reactions = []
        for reaction in data_all_reactions.index:
            count = sum(data_all_reactions.loc[reaction])
            if count > nb_species - (out + 1):
                filtered_reactions.append(reaction)
        return filtered_reactions


    def __generate_nb_rxn_sp(self) -> Dict[str, int]:
        """ Generate the nb_rxn_sp attribute. Returns a dictionary associating for each species the number of genes
        it has.

        Returns
        -------
        nb_genes_sp_dict : Dict[str, int]
            Dictionary Dict[species name, number of genes] associating for each species the number of genes it has.
        """
        nb_rxn_sp_dict = {}
        for sp in self.species_list:
            nb_rxn_sp_dict[sp] = int(sum(self.data_reactions[sp]))
        return nb_rxn_sp_dict


    def is_present(self, species: str, reaction: str, unique: bool = False) -> bool:
        """ Indicate if the reaction is present for the species (unique or not) : considered unique if only this species
        is having the reaction among all species

        Parameters
        ----------
        species: str
            species to be considered
        reaction: str
            reaction to be considered
        unique: bool, optional (default=False)
            True if the presence is unique, False otherwise

        Returns
        -------
        bool
            True if the reaction for the species is present (unique or not),
            False otherwise
        """
        if unique:
            row = list(self.data_reactions.loc[reaction])
            return self.data_reactions.loc[reaction, species] == sum(row) == 1
        else:
            return self.data_reactions.loc[reaction, species] == 1

    def get_rxn_present(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of present reactions (unique or not) : considered unique if
        only this species is having the reaction among all species

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None, will be all the species
        unique: bool, optional (default=False)
            True if the presence is unique, False otherwise

        Returns
        -------
        present_rxn_dict: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_reactions, Set[reactions]]]) dictionary associating for each species the number
            of present reactions and its set (unique or not)
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        present_rxn_dict = {}
        for sp in species:
            present_rxn = set()
            for rxn in self.reactions_list:
                if self.is_present(sp, rxn, unique):
                    present_rxn.add(rxn)
            present_rxn_dict[sp] = (len(present_rxn), present_rxn)
        return present_rxn_dict


    def get_rxn_presence(self, reaction: str or List[str] = None, unique: bool = False) -> \
            Dict[str, Tuple[Tuple[int, float], Set[str]]]:
        """ Returns for each reaction the number, the percentage, and the set of species having the reaction present
        (unique or not) : considered unique if only one species is having the reaction among all species

        Parameters
        ----------
        reaction: str or List[str], optional (default=None)
            reaction or list of reactions to be considered
            if None, will be all the reactions
        unique: bool, optional (default=False)
            True if the presence is unique, False otherwise

        Returns
        -------
        rxn_presence_dict: Dict[str, Tuple[Tuple[int, float], Set[str]]]
            (Dict[species, Tuple[Tuple[number_species, percentage_species], Set[reactions]]])
            dictionary associating for each reaction the number and the percentage of species having the reaction and
            its set (unique or not)
        """
        if reaction is None:
            reaction = self.reactions_list
        elif type(reaction) == str:
            reaction = [reaction]
        rxn_presence_dict = {}
        for rxn in reaction:
            rxn_presence = set()
            if rxn not in self.reactions_list:
                rxn_presence_dict[rxn] = (0, 0, rxn_presence)
            else:
                for sp in self.species_list:
                    if self.is_present(sp, rxn, unique):
                        rxn_presence.add(sp)
            rxn_presence_dict[rxn] = ((len(rxn_presence), len(rxn_presence)/self.nb_species), rxn_presence)
        return rxn_presence_dict

    def is_absent(self, species: str, reaction: str, unique: bool = False) -> bool:
        """ Indicate if the reaction is absent for the species (unique or not) : considered unique if only this species
        is not having the reaction among all species

        Parameters
        ----------
        species: str
            species to be considered
        reaction: str
            reaction to be considered
        unique: bool, optional (default=False)
            True if the absence is unique, False otherwise

        Returns
        -------
        bool
            True if the reaction for the species is absent (unique or not),
            False otherwise
        """
        if unique:
            row = list(self.data_reactions.loc[reaction])
            return self.data_reactions.loc[reaction, species] == 0 and sum(row) == self.nb_species - 1
        else:
            return self.data_reactions.loc[reaction, species] == 0

    def get_rxn_absent(self, species: str or List[str] = None, unique: bool = False) -> Dict[str, Tuple[int, Set[str]]]:
        """ Returns for each species the number and the set of absent reactions (unique or not) : considered unique if
        only this species is not having the reaction among all species

        Parameters
        ----------
        species: str or List[str], optional (default=None)
            species or list of species to be considered
            if None will be all the species
        unique: bool, optional (default=False)
            True if the absence is unique, False otherwise

        Returns
        -------
        absent_rxn_dict: Dict[str, Tuple[int, Set[str]]]
            (Dict[species, Tuple[number_reactions, Set[reactions]]]) dictionary associating for each species the number
            of absent reactions and its set (unique or not)
        """
        if species is None:
            species = self.species_list
        elif type(species) == str:
            species = [species]
        absent_rxn_dict = {}
        for sp in species:
            absent_rxn = set()
            for rxn in self.reactions_list:
                if self.is_absent(sp, rxn, unique):
                    absent_rxn.add(rxn)
            absent_rxn_dict[sp] = (len(absent_rxn), absent_rxn)
        return absent_rxn_dict

    def get_genes_assoc(self, reactions_list: str or List[str] = None,) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
        """ Returns a dictionary of genes associated with each reaction for each species. Can write proteins sequences
        associated in fasta files.
        Parameters
        ----------
        reactions_list : str or List[str], optional (default=None)
            List of reactions to find genes associated with.
            If None, will be reactions_list attribute.

        Returns
        -------
        genes_assoc : Dict[str, Dict[str, Dict[str, List[str]]]]
            Dictionary of genes associated with each reaction for each species
        """
        genes_assoc = {}
        if reactions_list is None:
            reactions_list = self.reactions_list
        elif type(reactions_list) == str:
            reactions_list = [reactions_list]

        for reaction in reactions_list:
            genes_assoc[reaction] = {}
        for reaction in reactions_list:
            for species in self.species_list:
                genes_assoc[reaction][species] = \
                    set(str(self.data_genes_assoc[species + self.STR_GENE_ASSOC][reaction]).split(";"))

        return genes_assoc
