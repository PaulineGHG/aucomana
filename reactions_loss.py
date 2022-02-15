import pandas as pd
from typing import Dict, List, Tuple


class Reactions:

    def __init__(self, file_reactions_tsv: str, species_list: List[str] = None):
        self.species_list = species_list
        self.data_reactions, \
            self.data_genes_assoc, \
            self.reactions_list = self.__init_data(file_reactions_tsv)
        self.nb_reactions, self.nb_species = self.data_reactions.shape
        self.reactions_loss = self.__init_reactions_loss()

    def __init_data(self, file_reactions_tsv: str):
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
        self.species_list = []
        for x in data.columns[1:]:
            if x[-7:] == '(sep=;)':
                break
            self.species_list.append(x)

    def __get_filtered_reactions(self, data_all_reactions, out: int = 1):
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
        loss = []
        for reaction in self.reactions_list:
            if self.data_reactions[species][reaction] == 0:
                loss.append(reaction)
        return len(loss), loss

    def __init_reactions_loss(self):
        loss_all = {}
        for species in self.species_list:
            loss_all[species] = self.__get_reactions_loss_1_species(species)
        return loss_all

    def get_genes_assoc(self, reactions_list: List[str] = None):
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

    def get_common_reactions(self, data_to_compare, species: str):
        reaction_list_self = self.reactions_loss[species][1]
        common_reactions = []
        for r in reaction_list_self:
            if r in data_to_compare.reactions_loss[species][1]:
                common_reactions.append(r)
        return len(common_reactions), common_reactions

    @staticmethod
    def print_genes_assoc(dict_genes_assoc: Dict[str, Dict[str, str]]):
        for species, gene_assoc in dict_genes_assoc.items():
            print(f"\n{species} :\n{'=' * 50}")
            for reaction, genes in gene_assoc.items():
                print(f"{reaction} : {genes}")
