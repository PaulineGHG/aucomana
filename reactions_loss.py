"""
Reactions class
"""
from typing import Dict, List, Tuple, Set
import pandas as pd
import json
import os


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
    nb_common_reac = 0
    nb_genes_assoc = 0
    STR_GENE_ASSOC = "_genes_assoc (sep=;)"

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

    def get_genes_assoc(self, interest_species: str,  reactions_set: Set[str] = None,
                        output_file=False) -> Dict[str, Dict[str, List[str]]]:
        """
        Parameters
        ----------
        interest_species : str
            species of interest
        reactions_set : List[str], optional (default=None)
            Set of reactions to find genes associated with.
            If None, will be reactions_list attribute.
        output_file : Bool, optional (default=False)
            If False, return the gene_assoc dict
            If true write the results with prot sequence

        Returns
        -------
        genes_assoc : Dict[str, Dict[str, List[str]]]
            Dictionary of genes associated with each reaction for each species
        """
        genes_assoc = {interest_species: {}}
        if reactions_set is None:
            reactions_set = self.reactions_loss[interest_species][1]
        for reaction in reactions_set:
            genes_assoc[interest_species][reaction] = {}
        for reaction in reactions_set:
            for species in self.species_list:
                genes_assoc[interest_species][reaction][species] = str(
                    self.data_genes_assoc[species + self.STR_GENE_ASSOC][reaction]).split(";")
        if not output_file:
            return genes_assoc
        else:
            self.nb_genes_assoc += 1
            self.write_genes(genes_assoc, interest_species)

    def write_genes(self, genes_assoc, interest_species):
        fa_file_path = f"data/{self.name.split('_')[0]}_studied_organism/"
        os.makedirs(f"outputs/genes_assoc/{self.nb_genes_assoc}_{interest_species}/")
        file_out_json = f"outputs/genes_assoc/{self.nb_genes_assoc}_{interest_species}/" \
                        f"{self.name.split('_')[0]}_genes_assoc.json"
        dir_out_fasta = f"outputs/genes_assoc/{self.nb_genes_assoc}_{interest_species}/"

        with open(file_out_json, 'w') as o:
            json.dump(genes_assoc, o, indent=4)

        genes_dict = genes_assoc[interest_species]
        for reaction, species_dict in genes_dict.items():
            fasta = dir_out_fasta + reaction + ".fa"
            with open(fasta, 'w') as o:
                for species, genes_list in species_dict.items():
                    genes_fa = f"{fa_file_path}{species}/{species}.faa"
                    with open(genes_fa, 'r') as f:
                        write_seq = False
                        for line in f:
                            if line[0] == ">":
                                write_seq = False
                                gene = line[1:][:-1]
                                if gene in genes_list:
                                    o.write(f">{gene}\n")
                                    write_seq = True
                            elif write_seq:
                                o.write(line)

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
            cls.nb_common_reac += 1
            cls.__write_common_reactions_json(datas, common_reactions, species)
        elif output_file == "txt":
            cls.nb_common_reac += 1
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
        outfile_name = f'outputs/common_reac/analyse_{cls.nb_common_reac}.txt'
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
        outfile_name = f'outputs/common_reac/analyse_{cls.nb_common_reac}.json'
        data = {"Compared files": [data.name for data in datas_list],
                "Interest species": species,
                "Number of common reactions": len(common_reactions),
                "Reactions": [reaction for reaction in common_reactions]}
        with open(outfile_name, 'w') as o:
            json.dump(data, o, indent=4)

