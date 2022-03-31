import os
from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.init_analysis import PATH_RUNS


def get_cat_l(reactions_file, organisms_file, cat):
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
            if l[1] == cat and l[0] in species_l:
                brown_l.append(l[0])
        return brown_l


def write_cut_reactions_file(original_file, cut_nb, reac_list):
    name = original_file.split("/")[-1]
    with open(original_file, "r") as f, open(f"outputs/cut_reactions_data/cut{cut_nb}_{name}", "w") as o:
        for line in f:
            l = line.split("\t")
            if l[0] in reac_list or l[0] == "reaction":
                o.write(line)


def get_reactions_inst(path_runs, org_tsv, cat=None, out=None):
    r_dic = {}
    for run in os.listdir(path_runs):
        r_path = os.path.join(path_runs, run, "analysis", "all", "reactions.tsv")
        if os.path.exists(r_path):
            if cat is not None:
                cat = get_cat_l(r_path, org_tsv, cat)
            r_dic[run] = Reactions(r_path, cat, out)
    return r_dic


def get_pathways_inst(path_runs, org_tsv, cat=None, out=None):
    p_dic = {}
    for run in os.listdir(path_runs):
        r_path = os.path.join(path_runs, run, "analysis", "all", "reactions.tsv")
        p_path = os.path.join(path_runs, run, "analysis", "all", "pathways.tsv")
        if os.path.exists(r_path) and os.path.exists(p_path):
            if cat is not None:
                cat = get_cat_l(r_path, org_tsv, cat)
            p_dic[run] = Pathways(p_path, cat, out)
    return p_dic


# create_folders(PATH_STUDY)

# ### FILES #######################################################################################

# DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"
ORG_TSV = "../data/species_group.tsv"

# ### Select species ##############################################################################

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

LELS = 'Laminarionema_elsbetiae'
SLAT = "Saccharina_latissima_FEMALE"
PLAC = "Pleurocladia_lacustris"

# ### Class instances #############################################################################

# REACTIONS = get_reactions_inst(PATH_RUNS, ORG_TSV, "brown", 1)

PATHWAYS = get_pathways_inst(PATH_RUNS, ORG_TSV, "brown")

R01 = "run01"
R02 = "run02"
R03 = "run03"
R04 = "run04"
R40 = "run40"

# ### Common Reactions ############################################################################

# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True)
# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True, union=True)

# ### genes assoc #################################################################################

# REACTIONS[R01].get_genes_assoc(LELS, {"12-OXOPHYTODIENOATE-REDUCTASE-RXN"}, output_file=True)

print(PATHWAYS[R01].get_pw_max([SLAT, PLAC], unique=False))


# with open("pathways_slat_plac.txt", "w") as f:
#     f.write("1) Pathways absents chez Slat et présents chez les autres : ")
#     p_l = PATHWAYS[R01].get_pw_absent(SLAT)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#     f.write("\nPathways absents chez Plac et présents chez les autres : ")
#     p_l = PATHWAYS[R01].get_pw_absent(PLAC)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#
#     f.write("\n\n2) Pathways de completion minimale chez Slat : ")
#     p_l = PATHWAYS[R01].get_pw_min(SLAT)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#     f.write("\nPathways de completion minimale chez Plac : ")
#     p_l = PATHWAYS[R01].get_pw_min(PLAC)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#
#     f.write("\n\n3) Pathways incomplets chez Slat et complets chez autres espèces : ")
#     p_l = PATHWAYS[R01].get_pw_incomplete(SLAT)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
#     f.write("\nPathways incomplets chez Plac et complets chez autres espèces : ")
#     p_l = PATHWAYS[R01].get_pw_incomplete(PLAC)
#     f.write(f"{p_l[0]}\n")
#     for p in p_l[1]:
#         f.write(f"{p}\n")
