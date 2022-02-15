import pandas as pd


class Reactions:

    def __init__(self, file_reactions_tsv, species_list=None):
        self.species_list = species_list
        self.data_reactions, \
        self.data_genes_assoc,\
        self.reactions_list = self.__init_data(file_reactions_tsv)
        self.nb_row, self.nb_col = self.data_reactions.shape
        self.reactions_loss = self.__init_reactions_loss()

    def __init_data(self, file_reactions_tsv):
        data = pd.read_csv(file_reactions_tsv, sep="\t", header=0, index_col='reaction')

        if self.species_list is None:
            self.species_list = []
            for x in data.columns[1:]:
                if x[-7:] == '(sep=;)':
                    break
                self.species_list.append(x)

        data_species_all_reac = data[self.species_list]

        genes_assoc_list = [x + "_genes_assoc (sep=;)" for x in self.species_list]
        data_genes_assoc = data[genes_assoc_list]

        reactions = self.__get_reactions(data_species_all_reac)

        return data_species_all_reac.loc[reactions], data_genes_assoc.loc[reactions], reactions

    @staticmethod
    def __get_reactions(data_all_reac):
        nb_row, nb_col = data_all_reac.shape
        reactions = []
        for r in data_all_reac.index:
            count = 0
            for s in data_all_reac.columns:
                if data_all_reac[s][r] == 1:
                    count += 1
            if count > nb_col - 2:
                reactions.append(r)
        return reactions

    def __get_reactions_loss_1_specie(self, specie):
        loss = []
        for r in self.data_reactions.index:
            if self.data_reactions[specie][r] == 0:
                loss.append(r)
        return len(loss), loss

    def __init_reactions_loss(self):
        loss_all = {}
        for specie in self.species_list:
            loss_all[specie] = self.__get_reactions_loss_1_specie(specie)
        return loss_all

    def get_gene_assos(self, reactions_list=None):
        gene_assoc = {}
        for specie in self.species_list:
            gene_assoc[specie] = {}
        if reactions_list is None:
            reactions_list = self.reactions_list
        for specie in self.species_list:
            for reaction in reactions_list:
                gene_assoc[specie][reaction] = self.data_genes_assoc[specie + "_genes_assoc (sep=;)"][reaction]
        return gene_assoc

    def get_common_reactions(self, data_to_compare, specie):
        reaction_list_self = self.reactions_loss[specie][1]
        common_reac = []
        for r in reaction_list_self:
            if r in data_to_compare.reactions_loss[specie][1]:
                common_reac.append(r)
        return len(common_reac), common_reac

    def print_gene_assoc(self, dict_gene_assoc):
        for specie, gene_assoc in dict_gene_assoc.items():
            print(f"\n{specie} :\n==============================================================")
            for reaction, genes in gene_assoc.items():
                print(f"{reaction} : {genes}")

DATA_FILE_1 = "run1_reactions.tsv"

DATA_FILE_40 = 'run40_reactions.tsv'
BROWN_ALGAE_40 = ['Thalassiosira_pseudonana',
                  'Fragilariopsis_cylindrus',
                  'Phaeodactylum_tricornutum',
                  'Nannochloropsis_gaditana',
                  'Ectocarpus_siliculosus',
                  'Ectocarpus_crouaniorum',
                  'Ectocarpus_subulatus',
                  'Ectocarpus_fasciculatus',
                  'Scytosiphon_lomentaria',
                  'Porterinema_fluviatile',
                  'Nemacystus_decipiens',
                  'Cladosiphon_okamuranus',
                  'Laminarionema_elsbetiae',
                  'Saccharina_japonica',
                  'Undaria_pinnatifida']


# R1 = Reactions(DATA_FILE_1)
R40 = Reactions(DATA_FILE_40, BROWN_ALGAE_40)


R40_loss_lami_e = R40.reactions_loss['Laminarionema_elsbetiae']
# R1_loss_lami_e = R1.reactions_loss['Laminarionema_elsbetiae']
# print(R40_loss_lami_e)
# print(R1_loss_lami_e)
# print(R1.get_common_reactions(R40, 'Laminarionema_elsbetiae'))
# print(R40.reactions_loss)
# print(R40.data_genes_assoc)
R40.print_gene_assoc(R40.get_gene_assos(R40_loss_lami_e[1]))


