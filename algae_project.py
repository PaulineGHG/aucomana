from reactions_loss import Reactions
from pandas_ods_reader import read_ods


def reactions_from_file(file):
    reactions = []
    d = read_ods(file)
    for i in range(d.shape[0]):
        r = d.loc[i, "Identifiant Metacyc"]
        reactions.append(r)
    return reactions


def get_brown_algae_l(reactions_file, organisms_file):
    with open(organisms_file, "r") as org_f, open(reactions_file, "r") as rea_f:
        species_l = set()
        brown_l = []
        for l in rea_f:
            l = l.split("\t")
            for x in l[1:]:
                if x[-7:] == '(sep=;)':
                    break
                species_l.add(x)
            break
        for l in org_f:
            l = l.split()
            if l[1] == "brown" and l[0] in species_l:
                brown_l.append(l[0])
        return brown_l


def write_cut_reactions_file(original_file, cut_nb):
    name = original_file.split("/")[-1]
    with open(original_file, "r") as f, open(f"outputs/cut_reactions_data/cut{cut_nb}_{name}", "w") as o:
        for line in f:
            l = line.split("\t")
            if l[0] in reac_list or l[0] == "reaction":
                o.write(line)


# ###Description ###

# 40 : run 40+7
# 01 : run Pauline all species
# A0 : run Alexandre
# A1 : run Pauline, Alexandre like exactly
# A2 : run Pauline, Alexandre like improved

# ### FILES ####

DATA_FILE_01 = "data/reactions_data/run01_reactions.tsv"
DATA_FILE_03 = "data/reactions_data/run03_reactions.tsv"
DATA_FILE_40 = 'data/reactions_data/run40_reactions.tsv'
DATA_FILE_A0 = 'data/reactions_data/runA0_reactions.tsv'
DATA_FILE_A1 = 'data/reactions_data/runA1_reactions.tsv'
DATA_FILE_A2 = 'data/reactions_data/runA2_reactions.tsv'
DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"
ORG_TSV = "data/species_group.tsv"

# ### Select species ###

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

LAMINARIONEMA_E = 'Laminarionema_elsbetiae'

# ### Class instances ###

BROWN_ALGAE_01 = get_brown_algae_l(DATA_FILE_01, ORG_TSV)
BROWN_ALGAE_03 = get_brown_algae_l(DATA_FILE_03, ORG_TSV)

# R01 = Reactions(DATA_FILE_01, BROWN_ALGAE_01, out=2, prio=("Dictyota_dichotoma_m", "Desmarestia_herbacea_m"))
# R40 = Reactions(DATA_FILE_40, BROWN_ALGAE_40)
# RA0 = Reactions(DATA_FILE_A0)
# RA1 = Reactions(DATA_FILE_A1)
# RA2 = Reactions(DATA_FILE_A2)
R03 = Reactions(DATA_FILE_03)

reac_lostA = reactions_from_file(DATA_LELSB_LOSSES)


# ### Laminarionema loss ###

# print(R40.reactions_loss[LAMINARIONEMA_E])
# print(R01.reactions_loss[LAMINARIONEMA_E])
# print(RA0.reactions_loss[LAMINARIONEMA_E])
# print(RA1.reactions_loss[LAMINARIONEMA_E])
# print(RA2.reactions_loss[LAMINARIONEMA_E])


# ### Common Reactions ###

# print(Reactions.get_common_reactions([R01, R40, RA0], LAMINARIONEMA_E, output_file=None))
# print(Reactions.get_common_reactions([R01, R40, RA1], LAMINARIONEMA_E, output_file=None))
# print(Reactions.get_common_reactions([R01, R40, RA2], LAMINARIONEMA_E, output_file=None))


# ### reactions lost percentage (based on highly shared reactions)

# for k, v in R01.reactions_loss.items():
#     print(k, " : ", round((v[0]/R01.nb_reactions)*100, 3), "%")

# ### genes assoc ###

# print(R01.get_genes_assoc(LAMINARIONEMA_E,
#                           Reactions.get_common_reactions([R01, R40, RA2], LAMINARIONEMA_E)[1],
#                           output_file=True))

reac_list = R03.reactions_list

# write_cut_reactions_file(DATA_FILE_03, 4)

print(R03.get_reaction_nb(R03.species_list))