from analysis_runs.analysis import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ### FILES #######################################################################################

ORG_FILE = 'data/species_group.tsv'
PATH_STUDY = os.getcwd()
PATH_RUNS = 'data/runs/'
# PATH_RUNS = '/home/phamongi/Documents/Runs/'
# DATA_LELSB_LOSSES = "data/Lelsb_losses.ods"
# ORG_TSV = "data/species_group.tsv"
phylo_f1 = "data/Phaeoexplorer_MLtree_rooted.nex"
phylo_f2 = "data/SpeciesTree_rooted.nex"

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

R01 = "run01"
R02 = "run02"
R03 = "run03"
R04 = "run04"
R40 = "run40"

A = Analysis(PATH_RUNS, PATH_STUDY, ORG_FILE)

# REACTIONS = A.get_reactions_inst(runs=[R04], out=1, group=("brown", 1))
PATHWAYS = A.get_pathways_inst(runs=[R04], group=("brown", 1))
# GENES = A.get_genes_inst(runs=[R04], group=("brown", 1))
# METABOLITES = A.get_metabolites_inst(runs=[R04], group=("brown", 1))

# ### Common Reactions ############################################################################

# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True)
# Reactions.get_common_reactions([R01, R40, RA2, R03], LELS, output_file=True, union=True)

# ### genes assoc #################################################################################

# REACTIONS[R01].get_genes_assoc(LELS, {"12-OXOPHYTODIENOATE-REDUCTASE-RXN"}, output_file=True)
# print(REACTIONS[R04].get_genes_assoc({"12-OXOPHYTODIENOATE-REDUCTASE-RXN",
#                                       "LEUKOTRIENE-C4-SYNTHASE-RXN",
#                                       "PROSTAGLANDIN-E-SYNTHASE-RXN"}, output_file=True))


# ### pathways ##################################################################################

# print(PATHWAYS[R01].get_pw_complete([SLAT, PLAC], unique=True))
# PATHWAYS[R04].convert_df_to_binary(1, output_file=True)
# print(PATHWAYS[R04].get_pw_names())
# print(PATHWAYS)
# PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_allrxn", phylo_f2, 1000)
# print(GENES)
# print(METABOLITES)
# REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out10", phylo_f2, n_boot=1000)
lr = A.get_grp_l(R04, ("LR", 2))
sr = A.get_grp_l(R04, ("SR", 2))

# oxy = ["12-OXOPHYTODIENOATE-REDUCTASE-RXN", "LEUKOTRIENE-C4-SYNTHASE-RXN", "PROSTAGLANDIN-E-SYNTHASE-RXN"]
# loss = REACTIONS[R04].reactions_loss
# loss_l = []
# for sp, val in loss.items():
#     if sp in sr:
#         print(sp, val[0])
#         loss_l.append(val[0])
# print(np.mean(loss_l))
# print(np.sqrt(np.var(loss_l)))
species = PATHWAYS[R04].species_list

abs = PATHWAYS[R04].get_pw_absent(species=species, unique=True)
abs_sr = []
abs_lr = []
abs_all = []
abs_lam = []

for sp, val in abs.items():
    print(sp, val[0])
    abs_all.append(val[0])
    if sp in sr:
        abs_sr.append(val[0])
    if sp in lr:
        abs_lr.append(val[0])
    if sp == LELS:
        abs_lam.append(val[0])

mini = PATHWAYS[R04].get_pw_min(species=species, unique=True)
min_sr = []
min_lr = []
min_all = []
min_lam = []

for sp, val in mini.items():
    print(sp, val[0])
    min_all.append(val[0])
    if sp in sr:
        min_sr.append(val[0])
    if sp in lr:
        min_lr.append(val[0])
    if sp == LELS:
        min_lam.append(val[0])


inc = PATHWAYS[R04].get_pw_incomplete(species=species, unique=True)
inc_sr = []
inc_lr = []
inc_all = []
inc_lam = []

for sp, val in inc.items():
    print(sp, val[0])
    inc_all.append(val[0])
    if sp in sr:
        inc_sr.append(val[0])
    if sp in lr:
        inc_lr.append(val[0])
    if sp == LELS:
        inc_lam.append(val[0])

# print(np.mean(loss_l))
# print(np.sqrt(np.var(loss_l)))

col = ["absent", "minimal", "incomplete"]
ind = ["all", "LR", "SR", "L.elsb"]
df = pd.DataFrame(columns=ind, index=col)
df.loc["absent", "all"] = round(float(np.mean(abs_all)), 2)
df.loc["absent", "LR"] = round(float(np.mean(abs_lr)), 2)
df.loc["absent", "SR"] = round(float(np.mean(abs_sr)), 2)
df.loc["absent", "L.elsb"] = round(float(np.mean(abs_lam)), 2)

df.loc["minimal", "all"] = round(float(np.mean(min_all)), 2)
df.loc["minimal", "LR"] = round(float(np.mean(min_lr)), 2)
df.loc["minimal", "SR"] = round(float(np.mean(min_sr)), 2)
df.loc["minimal", "L.elsb"] = round(float(np.mean(min_lam)), 2)

df.loc["incomplete", "all"] = round(float(np.mean(inc_all)), 2)
df.loc["incomplete", "LR"] = round(float(np.mean(inc_lr)), 2)
df.loc["incomplete", "SR"] = round(float(np.mean(inc_sr)), 2)
df.loc["incomplete", "L.elsb"] = round(float(np.mean(inc_lam)), 2)

print(df)
ax = df.plot.bar()
for container in ax.containers:
    ax.bar_label(container)
ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
ax.set_ylabel("Number of pathways")
ax.set_title("Unique absent/minimal/incomplete pathway for L.eslb compared to all/LR/SR groups mean")
plt.xticks(rotation=0)
plt.tight_layout()
plt.show()
