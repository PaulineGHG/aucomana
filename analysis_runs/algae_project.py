import os
from analysis_runs.reactions import Reactions
from analysis_runs.pathways import Pathways
from analysis_runs.init_analysis import PATH_RUNS, PATH_STUDY
from pandas_ods_reader import read_ods


def reactions_from_file(file):
    reactions = []
    d = read_ods(file)
    for i in range(d.shape[0]):
        r = d.loc[i, "Identifiant Metacyc"]
        reactions.append(r)
    return reactions


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
        if os.path.exists(r_path):
            if cat is not None:
                cat = get_cat_l(r_path, org_tsv, cat)
            p_dic[run] = Pathways(p_path, cat, out)
    return p_dic


# ### FILES #######################################################################################

# DATA_FILE_01 = "data/reactions_data/run01_reactions.tsv"
# DATA_FILE_03 = "data/reactions_data/run03_reactions.tsv"
# DATA_FILE_40 = 'data/reactions_data/run40_reactions.tsv'
# DATA_FILE_A0 = 'data/reactions_data/runA0_reactions.tsv'
# DATA_FILE_A1 = 'data/reactions_data/runA1_reactions.tsv'
# DATA_FILE_A2 = 'data/reactions_data/runA2_reactions.tsv'
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
# PATHWAYS = get_pathways_inst(PATH_RUNS, ORG_TSV, "brown", 1)

# ### Common Reactions ############################################################################

# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True)
# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True, union=True)

# ### genes assoc #################################################################################

# REACTIONS["run01"].get_genes_assoc(LELS, {"12-OXOPHYTODIENOATE-REDUCTASE-RXN"}, output_file=True)


# for sp in PATHWAYS["run01"].species_list:
#     print(sp, PATHWAYS["run01"].get_pw_lost_1_species(sp))
