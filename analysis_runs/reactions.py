"""
Reactions class
"""
from typing import Dict, List, Tuple, Set
import analysis_runs.dendrograms
import analysis_runs.analysis
import pandas as pd
import json
import os
import datetime


class Reactions:
    """
    Attributes
    ----------
    path_runs: str
        path of AuCoMe runs results
    path_study: str
        path of outputs_data of the study
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
    nb_reactions_sp : Dict[str, int]
        Dictionary indicating for each species how many reactions are found
    """
    nb_common_reac = 0
    nb_genes_assoc = 0
    nb_dend = 0
    STR_GENE_ASSOC = "_genes_assoc (sep=;)"

    def __init__(self, path_runs: str, path_study: str, file_reactions_tsv: str, species_list: List[str] = None,
                 out: int = None):
        """ Init the Reactions class

        Parameters
        ----------
        path_runs: str
            path of AuCoMe runs results
        path_study: str
            path of outputs_data of the study
        file_reactions_tsv : str
            file reactions.tsv output from aucome analysis
        species_list : List[str], optional (default=None)
            List of species to study (must correspond to their name in reactions.tsv file).
            If not specified, will contain all the species from reactions.tsv file.
        out : int, optional (default = None)
            number of species maximum not having the reaction for the reaction to be kept
        """
        self.path_runs = path_runs
        self.path_study = path_study
        self.name = file_reactions_tsv.split("/")[-4]
        self.species_list = species_list
        self.data_reactions, \
            self.data_genes_assoc, \
            self.reactions_list = self.__init_data(file_reactions_tsv, out)
        self.nb_reactions, self.nb_species = self.data_reactions.shape
        self.reactions_loss = self.__init_reactions_loss()
        self.nb_reactions_sp = self.__get_reaction_nb()

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

    def __get_reaction_nb(self):
        reactions_nb_dict = {}
        for species in self.species_list:
            reactions_nb_dict[species] = sum(self.data_reactions[species])
        return reactions_nb_dict

    def get_genes_assoc(self, reactions_set: Set[str] = None,
                        output_file=False) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
        """
        Parameters
        ----------
        reactions_set : List[str], optional (default=None)
            Set of reactions to find genes associated with.
            If None, will be reactions_list attribute.
        output_file : Bool, optional (default=False)
            If False, return the gene_assoc dict
            If true write the results with prot sequence

        Returns
        -------
        genes_assoc : Dict[str, Dict[str, Dict[str, List[str]]]]
            Dictionary of genes associated with each reaction for each species
        """
        genes_assoc = {}
        if reactions_set is None:
            reactions_set = self.reactions_list
        for reaction in reactions_set:
            genes_assoc[reaction] = {}
        for reaction in reactions_set:
            for species in self.species_list:
                genes_assoc[reaction][species] = str(
                    self.data_genes_assoc[species + self.STR_GENE_ASSOC][reaction]).split(";")
        if not output_file:
            return genes_assoc
        else:
            self.nb_genes_assoc += 1
            self.__write_genes(genes_assoc)

    def __write_genes(self, genes_assoc):
        """ Allow writing sequences of genes in fasta out files

        Parameters
        ----------
        genes_assoc : Dict[str, Dict[str, Dict[str, List[str]]]]
            Dictionary of genes associated with each reaction for each species

        Write results in {PATH_STUDY}/output_data/reactions_data/genes_assoc/{now}_{self.nb_genes_assoc}/
        """
        fa_file_path = f"{self.path_runs}/{self.name}/studied_organisms/"
        now = datetime.datetime.now().strftime('%d-%m-%Y_%Hh-%Mmin-%Ss')
        out_file = f"{self.path_study}/output_data/reactions_data/genes_assoc/{now}_{self.nb_genes_assoc}/"
        if not os.path.exists(out_file):
            os.makedirs(out_file)
        file_out_json = f"{out_file}{self.name}_genes_assoc.json"

        with open(file_out_json, 'w') as o:
            json.dump(genes_assoc, o, indent=4)

        for reaction, species_dict in genes_assoc.items():
            fasta = f"{out_file}{reaction}.fa"
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
                                    o.write(f">{species}|{gene}\n")
                                    write_seq = True
                            elif write_seq:
                                o.write(line)

    @classmethod
    def get_common_reactions(cls, datas: List["Reactions"], species: str, output_file=False,
                             union=False) -> Tuple[int, Set[str]]:
        """ Returns the reactions lost in common between at least 2 Reactions instance for 1
        common species

        Parameters
        ----------
        datas : List["Reactions"]
            List of Reactions instance to compare
        species : str
            Species of interest compare
        output_file : bool, optional (default=False)
            True : write json and txt output
            False : returns the output
        union : bool, optional (default=False)
            True : get the union set of reactions between runs
            False : get the intersection set of reactions between runs

        Returns
        -------
        Tuple[int, Set[str]]
            Number of reactions in common and their list
        """
        set_reac_list = []
        for data in datas:
            set_reac_list.append(data.reactions_loss[species][1])
        if not union:
            common_reactions = set.intersection(*set_reac_list)
        else:
            common_reactions = set.union(*set_reac_list)
        if output_file:
            cls.nb_common_reac += 1
            cls.__write_common_reactions_json(datas, common_reactions, species, union)
            cls.__write_common_reactions_txt(datas, common_reactions, species, union)
        else:
            return len(common_reactions), common_reactions

    def __write_common_reactions_txt(self, datas_list: List["Reactions"],
                                     common_reactions: Set[str], species: str, union: bool):
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
        now = datetime.datetime.now().strftime('%d_%m_%Y__%Hh_%Mmin_%Ss')
        if not union:
            outfile_name = f'{self.path_study}/output_data/reactions_data/common_reac/intersection/' \
                           f'{now}_{self.nb_common_reac}.txt'
        else:
            outfile_name = f'{self.path_study}/output_data/reactions_data/common_reac/union/{now}_' \
                           f'{self.nb_common_reac}.txt'
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

    def generate_rnx_dendrogram(self, a, name=None, phylo_file=None, n_boot=10000):
        if name is None:
            self.nb_dend += 1
            name = f"dendrogram{self.nb_dend}"
        name = "rnx_" + name
        d = analysis_runs.dendrograms.Dendrogram(a, self.path_study, self.data_reactions,
                                                 self.name, name, phylo_file)
        d.get_dendro_pvclust(n_boot)

    def __write_common_reactions_json(self, datas_list: List["Reactions"],
                                      common_reactions: Set[str], species: str, union: bool):
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
        now = datetime.datetime.now().strftime('%d_%m_%Y__%Hh_%Mmin_%Ss')
        if not union:
            outfile_name = f'{self.path_study}/output_data/reactions_data/common_reac/intersection/' \
                           f'{now}_{self.nb_common_reac}.json'
        else:
            outfile_name = f'{self.path_study}/output_data/reactions_data/common_reac/union/{now}_' \
                           f'{self.nb_common_reac}.json'
        data = {"Compared files": [data.name for data in datas_list],
                "Interest species": species,
                "Number of common reactions": len(common_reactions),
                "Reactions": [reaction for reaction in common_reactions]}
        with open(outfile_name, 'w') as o:
            json.dump(data, o, indent=4)

