from aucomana.aucomana import AuCoMAna

g_file = 'Runs/run62/analysis/group_template.tsv'
comp_dir = 'Runs/run62/analysis/all'
spha = 'Sphacelaria-rigidula_FEMALE'

Auc = AuCoMAna(comp_dir, g_file)
Auc.group_pathway_completion_comparison(('LR', 'SR', spha), 'Phaeophyceae')
