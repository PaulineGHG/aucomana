import os


def create_folders(path_study):
    arbo = [['output_data'],
            ['output_data', 'compare_groups'],
            ['output_data', 'dendro_tanglegrams'],
            ['output_data', 'pathways_data'],
            ['output_data', 'pathways_data', 'binary_df'],
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


PATH_STUDY = os.getcwd()
PATH_RUNS = 'data/runs/'
# PATH_RUNS = '/home/phamongi/Documents/Runs/'
