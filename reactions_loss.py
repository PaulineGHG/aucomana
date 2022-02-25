"""
Reactions class
"""
from typing import Dict, List, Tuple, Set
import pandas as pd
import json


class Reactions:
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
    nb_analysis = 0

    def __init__(self, file_reactions_tsv: str, species_list: List[str] = None, out: int = 1):
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
        self.name = file_reactions_tsv.split(".")[0].split("/")[-1]
        self.species_list = species_list
        self.data_reactions, \
            self.data_genes_assoc, \
            self.reactions_list = self.__init_data(file_reactions_tsv, out)
        self.nb_reactions, self.nb_species = self.data_reactions.shape
        self.reactions_loss = self.__init_reactions_loss()

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
        genes_assoc_list = [x + "_genes_assoc (sep=;)" for x in self.species_list]
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
        filtered_reactions = []
        for reaction in data_all_reactions.index:
            count = sum(data_all_reactions.loc[reaction])
            if count > nb_species - (out + 1):
                filtered_reactions.append(reaction)
        return filtered_reactions

    def __get_reactions_loss_1_species(self, species: str) -> Tuple[int, Set[str]]:
        """ Capture the reaction lost for a given species

        Parameters
        ----------
        species : str
            Name of the species

        Returns
        -------
        Tuple[int, Set[str]]
            Tuple with the number of reactions lost and the list of these reactions
        """
        loss = set()
        for reaction in self.reactions_list:
            if self.data_reactions[species][reaction] == 0:
                loss.add(reaction)
        return len(loss), loss

    def __init_reactions_loss(self) -> Dict[str, Tuple[int, Set[str]]]:
        """ Init the reactions_loss attribute

        Returns
        -------
        reactions_loss : Dict[str, Tuple[int, Set[str]]]
            reactions_loss attribute
        """
        reactions_loss = {}
        for species in self.species_list:
            reactions_loss[species] = self.__get_reactions_loss_1_species(species)
        return reactions_loss

    def get_genes_assoc(self, species: str,  reactions_list: Set[str]=None) -> Dict[str, Dict[str, List[str]]]:
        """
        Parameters
        ----------
        reactions_list : List[str], optional (default=None)
            List of reactions to find genes associated with.
            If None, will be reactions_list attribute.

        Returns
        -------
        genes_assoc : Dict[str, Dict[str, List[str]]]
            Dictionary of genes associated with each reaction for each species
        """
        genes_assoc = {}
        if reactions_list is None:
            reactions_list = self.reactions_loss[species][1]
        for reaction in reactions_list:
            genes_assoc[reaction] = {}
        for reaction in reactions_list:
            for species in self.species_list:
                genes_assoc[reaction][species] = str(self.data_genes_assoc[species +
                                                                           "_genes_assoc (sep=;)"]
                                                     [reaction]).split(";")
        return genes_assoc

    @classmethod
    def get_common_reactions(cls, datas: List["Reactions"], species: str, output_file=None) \
            -> Tuple[int, Set[str]]:
        """ Returns the reactions lost in common between at least 2 Reactions instance for 1
        common species

        Parameters
        ----------
        datas : List["Reactions"]
            List of Reactions instance to compare
        species : str
            Species of interest compare
        output_file : str, optional (default=None)
            'json' to write output in a .json file
            'txt' to write output in a .txt file
            None to return the output

        Returns
        -------
        Tuple[int, Set[str]]
            Number of reactions in common and their list
        """
        set_reac_list = []
        for data in datas:
            set_reac_list.append(data.reactions_loss[species][1])
        common_reactions = set.intersection(*set_reac_list)
        if output_file is None:
            return len(common_reactions), common_reactions
        elif output_file == "json":
            cls.nb_analysis += 1
            cls.__write_common_reactions_json(datas, common_reactions, species)
        elif output_file == "txt":
            cls.nb_analysis += 1
            cls.__write_common_reactions_txt(datas, common_reactions, species)
        else:
            raise ValueError("output_file value must be 'json' or 'txt'")

    @classmethod
    def __write_common_reactions_txt(cls, datas_list: List["Reactions"],
                                     common_reactions: Set[str], species: str):
        """Write get_common_reactions results in a .txt file

        Parameters
        ----------
        datas_list : List["Reactions"]
            List of Reactions instance compared
        common_reactions : Set[str]
            Set of common reactions between all the datas
        species : str
            Interest species
        """
        outfile_name = f'outputs/analyse_{cls.nb_analysis}.txt'
        with open(outfile_name, 'w') as o:
            o.write("Compared files :\n"
                    "----------------\n")
            for data in datas_list:
                o.write(data.name + "\n")
            o.write(f"\nInterest species :\n"
                    f"------------------\n"
                    f"{species}\n\n")
            o.write(f"Number of common reactions : {len(common_reactions)}\n"
                    f"----------------------------\n\n"
                    f"Reactions :\n"
                    f"-----------\n")
            for reaction in common_reactions:
                o.write(reaction + "\n")

    @classmethod
    def __write_common_reactions_json(cls, datas_list: List["Reactions"],
                                      common_reactions: Set[str], species: str):
        """Write get_common_reactions results in a .json file

       Parameters
       ----------
       datas_list : List["Reactions"]
           List of Reactions instance compared
       common_reactions : Set[str]
           Set of common reactions between all the datas
       species : str
           Interest species
       """
        outfile_name = f'outputs/analyse_{cls.nb_analysis}.json'
        data = {"Compared files": [data.name for data in datas_list],
                "Interest species": species,
                "Number of common reactions": len(common_reactions),
                "Reactions": [reaction for reaction in common_reactions]}
        with open(outfile_name, 'w') as o:
            json.dump(data, o, indent=4)

