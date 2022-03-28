import os
import shutil


def create_folders(path_study):
    arbo = [['input_data'], ['output_data'],
            ['output_data', 'dendro_tanglegrams'],
            ['output_data', 'pathways_data'],
            ['output_data', 'reactions_data'],
            ['output_data', 'reactions_data', 'common_reac'],
            ['output_data', 'reactions_data', 'common_reac', 'union'],
            ['output_data', 'reactions_data', 'common_reac', 'intersection'],
            ['output_data', 'reactions_data', 'genes_assoc'],
            ['output_data', 'genes_data'],
            ['output_data', 'metabolites_data'],
            ['output_data', 'pages']]
    for folder in arbo:
        folder_path = os.path.join(path_study, *folder)
        print(f'creating folder {folder_path}')
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)


def fill_input_auto(path_study, path_runs, runs=None):
    for run in os.listdir(path_runs):
        analysis = os.path.join(path_runs, run, 'analysis')
        if os.path.exists(analysis):
            for grp in os.listdir(analysis):
                if grp != 'group_template.tsv':
                    analysis_grp = os.path.join(analysis, grp)
                    nw_obj_l = ['reactions', 'pathways', 'genes', 'metabolites']
                    for nw_obj in nw_obj_l:
                        copy_file_tsv(path_study, run, analysis_grp, grp, nw_obj)


def copy_file_tsv(path_study, run, analysis_grp, grp, nw_obj):
    r_file = os.path.join(analysis_grp, f'{nw_obj}.tsv')
    r_dest = os.path.join(path_study, 'input_data', f'{nw_obj}_data', f'{run}_{grp}_{nw_obj}.tsv')
    if os.path.exists(r_file):
        shutil.copy(r_file, r_dest)
        print(f'copy {r_file} in {r_dest}')


PATH_STUDY = os.getcwd()
PATH_RUNS = '/home/paulinehg/Documents/Aucome_runs'

# create_folders(PATH_STUDY)
# fill_input_auto(PATH_STUDY, PATH_RUNS)


