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
from Bio import SeqIO


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

    def get_genes_assoc(self, reactions_list: str or List[str] = None,
                        output_file=False) -> Dict[str, Dict[str, Dict[str, List[str]]]]:
        """ Returns a dictionary of genes associated with each reaction for each species. Can write proteins sequences
        associated in fasta files.

        Parameters
        ----------
        reactions_list : str or List[str], optional (default=None)
            List of reactions to find genes associated with.
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
        if reactions_list is None:
            reactions_list = self.reactions_list
        elif type(reactions_list) == str:
            reactions_list = [reactions_list]
        for reaction in reactions_list:
            genes_assoc[reaction] = {}
        for reaction in reactions_list:
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
        # PATH TO FASTA FILES
        fa_file_path = os.path.join(self.path_runs, self.name, "studied_organisms")

        # GET OUTFILE PATH NOT EXISTING IN DIRECTORY
        out_file = os.path.join(self.path_study, "output_data", "reactions_data", "genes_assoc",
                                str(self.nb_genes_assoc))
        while os.path.exists(out_file):
            self.nb_genes_assoc += 1
            out_file = os.path.join(self.path_study, "output_data", "reactions_data", "genes_assoc",
                                    str(self.nb_genes_assoc))
        os.makedirs(out_file)

        # CREATE JSON FILE
        file_out_json = os.path.join(out_file, f"{self.name}_genes_assoc.json")
        with open(file_out_json, 'w') as o:
            json.dump(genes_assoc, o, indent=4)

        # FILL FASTA FILE
        for reaction, species_dict in genes_assoc.items():
            fasta_file = os.path.join(out_file, f"{reaction}.fa")
            with open(fasta_file, 'w') as o:
                for species, genes_list in species_dict.items():
                    genes_fa_file = os.path.join(fa_file_path, species, f"{species}.faa")
                    seq_dict = SeqIO.to_dict(SeqIO.parse(genes_fa_file, "fasta"))
                    for gene in genes_list:
                        id_gene = f">{species}|{gene}"
                        seq = str(seq_dict[gene].seq)
                        o.write(id_gene + "\n")
                        o.write(seq + "\n")

    def generate_rxn_dendrogram(self, name=None, phylo_file=None, n_boot=10000):
        if name is None:
            self.nb_dend += 1
            name = f"dendrogram{self.nb_dend}"
        name = "Rxn_" + name
        d = analysis_runs.dendrograms.Dendrogram(self.path_runs, self.path_study, self.data_reactions,
                                                 self.name, name, phylo_file)
        d.get_dendro_pvclust(n_boot)

    def get_reactions_names(self):
        """ Associate for each reaction in reactions_list, its common name in a dictionary.

        Returns
        -------
        rxn_names : Dict[str, str]
            Dictionary Dict[ID of rxn, common name of rxn] associating common reactions names to their ID.
        """
        rxn_set = set(self.reactions_list)
        rxn_names = {}
        padmet_file = os.path.join(self.path_runs, self.name, "analysis", "all",
                                   "all_panmetabolism.padmet")
        with open(padmet_file, "r") as f:
            for l in f:
                if "COMMON-NAME" in l and "reaction" in l:
                    l = l.split("\t")
                    if l[1] in rxn_set:
                        name_index = l.index('COMMON-NAME') + 1
                        rxn_names[l[1]] = l[name_index].strip('\n')
        return rxn_names

