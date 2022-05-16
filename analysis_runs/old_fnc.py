# ## analysis.py #########################################################################################################

from pandas_ods_reader import read_ods

# def reactions_from_file(file):
#     reactions = []
#     d = read_ods(file)
#     for i in range(d.shape[0]):
#         r = d.loc[i, "Identifiant Metacyc"]
#         reactions.append(r)
#     return reactions

# ## init_analysis.py #################################################################################################

# def fill_input_auto(path_study, path_runs, runs=None):
#     for run in os.listdir(path_runs):
#         analysis = os.path.join(path_runs, run, 'analysis')
#         if os.path.exists(analysis):
#             for grp in os.listdir(analysis):
#                 if grp != 'group_template.tsv':
#                     analysis_grp = os.path.join(analysis, grp)
#                     nw_obj_l = ['reactions', 'pathways', 'genes', 'metabolites']
#                     for nw_obj in nw_obj_l:
#                         copy_file_tsv(path_study, run, analysis_grp, grp, nw_obj)
#
#
# def copy_file_tsv(path_study, run, analysis_grp, grp, nw_obj):
#     r_file = os.path.join(analysis_grp, f'{nw_obj}.tsv')
#     r_dest = os.path.join(path_study, 'input_data', f'{nw_obj}_data', f'{run}_{grp}_{nw_obj}.tsv')
#     if os.path.exists(r_file):
#         shutil.copy(r_file, r_dest)
#         print(f'copy {r_file} in {r_dest}')

# ## reactions.py #####################################################################################################

# Reactions methods :

# def __get_filtered_reactions(self, data_all_reactions: 'pd.DataFrame', out: int, prio) \
#         -> Set[str]:
#     """ Filter the reactions according to the number of species not having the reaction
#
#     Parameters
#     ----------
#     data_all_reactions:
#         Dataframe with filtered columns and unfiltered reactions (rows)
#     out: int
#         number of species maximum not having the reaction for the reaction to be kept
#
#     Returns
#     -------
#     filtered_reactions : List[str]
#         List of reactions filtered
#     """
#     nb_species = len(self.species_list)
#     filtered_reactions = set()
#     for reaction in data_all_reactions.index:
#         count = sum(data_all_reactions.loc[reaction])
#         if count > (out + 1):
#             filtered_reactions.add(reaction)
#         else:
#             if prio is not None:
#                 for p_sp in prio:
#                     if data_all_reactions.loc[reaction][p_sp] == 1:
#                         filtered_reactions.add(reaction)
#     return filtered_reactions

# def __get_filtered_reactions(self, data_all_reactions: 'pd.DataFrame', out: int, prio) \
#         -> Set[str]:
#     """ Filter the reactions according to the number of species not having the reaction
#
#     Parameters
#     ----------
#     data_all_reactions:
#         Dataframe with filtered columns and unfiltered reactions (rows)
#     out: int
#         number of species maximum not having the reaction for the reaction to be kept
#
#     Returns
#     -------
#     filtered_reactions : List[str]
#         List of reactions filtered
#     """
#     nb_species = len(self.species_list)
#     filtered_reactions = set()
#     for reaction in data_all_reactions.index:
#         count = sum(data_all_reactions.loc[reaction])
#         if count > out:
#             filtered_reactions.add(reaction)
#         else:
#             if prio is not None:
#                 add = True
#                 for p_sp in prio:
#                     if data_all_reactions.loc[reaction][p_sp] == 1:
#                         add = False
#                 if add:
#                     filtered_reactions.add(reaction)
#     return filtered_reactions
