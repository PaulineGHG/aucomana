import pandas as pd


class Reactions:

    def __init__(self, file_reactions_tsv, species_list=None):
        self.species_list = species_list
        self.data_reactions, self.reactions_list, self.data_genes_assoc = self.__init_data(file_reactions_tsv)
        self.nb_row, self.nb_col = self.data_reactions.shape
        self.reactions_loss = self.__init_reactions_loss()

    def __init_data(self, file_reactions_tsv):
        data = pd.read_csv(file_reactions_tsv, sep="\t", header=0)

        if self.species_list is None:
            self.species_list = []
            for x in data.columns[1:]:
                if x[-7:] == '(sep=;)':
                    break
                self.species_list.append(x)

        data_species_all_reac = data[self.species_list].to_numpy()

        reactions_list = data['reaction']

        genes_assoc_list = [x + "_genes_assoc (sep=;)" for x in self.species_list]
        data_genes_assoc = data[genes_assoc_list]

        reactions = self.__get_reactions(data_species_all_reac)

        return data_species_all_reac[reactions], list(reactions_list[reactions]), data_genes_assoc

    def mutate_species_list(self):
        new_list = []
        for s in self.species_list:
            print(s)
            keep = input("keep specie ? [y/n]")
            if keep != 'y' and keep != 'n':
                print("Not accepted, must be 'y' or 'n'")
            elif keep == 'y':
                new_list.append(s)
        self.species_list = new_list


    @staticmethod
    def __get_reactions(data_all_reac):
        nb_row, nb_col = data_all_reac.shape
        reactions = []
        for i in range(nb_row):
            count = 0
            for j in range(nb_col):
                if data_all_reac[i, j] == 1:
                    count += 1
            if count > nb_col - 2:
                reactions.append(i)
        return reactions

    def __get_reactions_loss_1_specie(self, specie):
        i_specie = self.species_list.index(specie)
        loss = []
        for i in range(self.nb_row):
            if self.data_reactions[i, i_specie] == 0:
                loss.append(self.reactions_list[i])
        return len(loss), loss

    def __init_reactions_loss(self):
        loss_all = {}
        for specie in self.species_list:
            loss_all[specie] = self.__get_reactions_loss_1_specie(specie)
        return loss_all

    def get_common_reactions(self, data_to_compare, specie):
        reaction_list_self = self.reactions_loss[specie][1]
        common_reac = []
        for r in reaction_list_self:
            if r in data_to_compare.reactions_loss[specie][1]:
                common_reac.append(r)
        return len(common_reac), common_reac


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


R1 = Reactions(DATA_FILE_1)
R40 = Reactions(DATA_FILE_40)


# R40_loss_lami_e = R40.reactions_loss['Laminarionema_elsbetiae']
# R1_loss_lami_e = R1.reactions_loss['Laminarionema_elsbetiae']
# print(R40_loss_lami_e)
# print(R1_loss_lami_e)
# print(R1.get_common_reactions(R40, 'Laminarionema_elsbetiae'))
# print(R1.reactions_loss)
# print(R40.data_genes_assoc)
print(R40.species_list)

