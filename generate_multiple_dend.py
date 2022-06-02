from analysis_runs.analysis import *


ORG_FILE = 'data/species_group.tsv'
PATH_STUDY = os.getcwd()
PATH_RUNS = 'data/runs/'
phylo_f2 = "data/SpeciesTree_rooted.nex"
R04 = "run04"
A = Analysis(PATH_RUNS, PATH_STUDY, ORG_FILE)


print("REACTIONS : 1")
REACTIONS = A.get_reactions_inst(runs=[R04], out=34)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out05perc", phylo_f2, n_boot=10000)

print("REACTIONS : 2")
REACTIONS = A.get_reactions_inst(runs=[R04], out=32)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out10perc", phylo_f2, n_boot=10000)

print("REACTIONS : 3")
REACTIONS = A.get_reactions_inst(runs=[R04], out=28)
REACTIONS[R04].generate_rnx_dendrogram(A, "all_sp_out20perc", phylo_f2, n_boot=10000)

print("PATHWAYS ALL RXN : 1")
PATHWAYS = A.get_pathways_inst(runs=[R04])
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_allrxn", phylo_f2, 10000)
print("PATHWAYS ALL RXN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_allrxn", phylo_f2, 10000)

PATHWAYS = A.get_pathways_inst(runs=[R04], nb_rnx_px_min=2)
print("PATHWAYS 2 MIN : 1")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_2min", phylo_f2, 10000)
print("PATHWAYS 2 MIN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_2min", phylo_f2, 10000)

PATHWAYS = A.get_pathways_inst(runs=[R04], nb_rnx_px_min=3)
print("PATHWAYS 3 MIN : 1")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_3min", phylo_f2, 10000)
print("PATHWAYS 3 MIN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_3min", phylo_f2, 10000)

PATHWAYS = A.get_pathways_inst(runs=[R04], nb_rnx_px_min=4)
print("PATHWAYS 4 MIN : 1")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.5, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 2")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.6, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 3")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.75, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 4")
PATHWAYS[R04].generate_pw_dendrogram(A, 0.8, False, "all_sp_4min", phylo_f2, 10000)
print("PATHWAYS 4 MIN : 5")
PATHWAYS[R04].generate_pw_dendrogram(A, 1, False, "all_sp_4min", phylo_f2, 10000)