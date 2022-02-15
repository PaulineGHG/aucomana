import pandas as pd
from typing import Dict, List, Tuple


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
    reactions_loss : Dict[str, Tuple[int, List[str]]]
        Dictionary indicating for each species how many and which reactions are lost
    """

    def __init__(self, file_reactions_tsv: str, species_list: List[str] = None):
        """ Init the Reactions class

        Parameters
        ----------
        file_reactions_tsv : str
            file reactions.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study (must correspond to their name in reactions.tsv file). If not specified, will
            contain all the species from reactions.tsv file.
        """
        self.species_list = species_list
        self.data_reactions, \
            self.data_genes_assoc, \
            self.reactions_list = self.__init_data(file_reactions_tsv)
        self.nb_reactions, self.nb_species = self.data_reactions.shape
        self.reactions_loss = self.__init_reactions_loss()

    def __init_data(self, file_reactions_tsv: str):
        """ Generate the data_reactions, data_genes_assoc and reactions_list attributes

        Parameters
        ----------
        file_reactions_tsv : str
            file reactions.tsv

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
        genes_assoc_list = [x + "_genes_assoc (sep=;)" for x in self.species_list]
        data_genes_assoc = data[genes_assoc_list]
        filtered_reactions = self.__get_filtered_reactions(data_species_all_reactions)
        return data_species_all_reactions.loc[filtered_reactions], data_genes_assoc.loc[filtered_reactions], \
            filtered_reactions

    def __generate_species_list(self, data):
        """ Generate the species_list attribute if is None

        Parameters
        ----------
        data :
            The dataframe created from reactions.tsv file
        """
        self.species_list = []
        for x in data.columns[1:]:
            if x[-7:] == '(sep=;)':
                break
            self.species_list.append(x)

    def __get_filtered_reactions(self, data_all_reactions, out: int = 1):
        """ Filter the reactions according to the number of species not having the reaction

        Parameters
        ----------
        data_all_reactions:
            Dataframe with filtered columns and unfiltered reactions (rows)
        out: int, optional (default=1)
            number of species maximum not having the reaction for the reaction to be kept

        Returns
        -------
        filtered_reactions : List[str]
            List of reactions filtered
        """
        nb_species = len(self.species_list)
        filtered_reactions = []
        for reaction in data_all_reactions.index:
            count = 0
            for species in data_all_reactions.columns:
                if data_all_reactions[species][reaction] == 1:
                    count += 1
            if count > nb_species - (out + 1):
                filtered_reactions.append(reaction)
        return filtered_reactions

    def __get_reactions_loss_1_species(self, species: str):
        """ Capture the reaction lost for a given species

        Parameters
        ----------
        species : str
            Name of the species

        Returns
        -------
        Tuple[int, List[str]]
            Tuple with the number of reactions lost and the list of these reactions
        """
        loss = []
        for reaction in self.reactions_list:
            if self.data_reactions[species][reaction] == 0:
                loss.append(reaction)
        return len(loss), loss

    def __init_reactions_loss(self):
        """ Init the reactions_loss attribute

        Returns
        -------
        reactions_loss : Dict[str, Tuple[int, List[str]]]
            reactions_loss attribute
        """
        reactions_loss = {}
        for species in self.species_list:
            reactions_loss[species] = self.__get_reactions_loss_1_species(species)
        return reactions_loss

    def get_genes_assoc(self, reactions_list: List[str] = None):
        """
        Parameters
        ----------
        reactions_list : List[str], optional (default=None)
            List of reactions to find genes associated with. If None, will pe reactions_list attribute.

        Returns
        -------
        genes_assoc : Dict[str, Dict[str, List[str]]]
            Dictionary of genes associated with each reaction for each species
        """
        genes_assoc = {}
        for species in self.species_list:
            genes_assoc[species] = {}
        if reactions_list is None:
            reactions_list = self.reactions_list
        for species in self.species_list:
            for reaction in reactions_list:
                genes_assoc[species][reaction] = \
                    str(self.data_genes_assoc[species + "_genes_assoc (sep=;)"][reaction]).split(";")
        return genes_assoc

    def get_common_reactions(self, data_to_compare: 'Reactions', species: str):
        """ Returns the reactions lost in common between 2 Reactions instance for 1 common species

        Parameters
        ----------
        data_to_compare : "Reactions"
            Reactions instance to compare with self instance
        species : str
            Species to compare

        Returns
        -------
        Tuple[int, List[str]]
            Number of reactions in common and their list
        """
        common_reactions = []
        for reaction in self.reactions_loss[species][1]:
            if reaction in data_to_compare.reactions_loss[species][1]:
                common_reactions.append(reaction)
        return len(common_reactions), common_reactions

    @staticmethod
    def print_genes_assoc(dict_genes_assoc: Dict[str, Dict[str, List[str]]]):
        for species, gene_assoc in dict_genes_assoc.items():
            print(f"\n{species} :\n{'=' * 50}")
            for reaction, genes in gene_assoc.items():
                print(f"{reaction} : {genes}")
