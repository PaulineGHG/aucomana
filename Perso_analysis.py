from aucomana.aucomana import AuCoMAna
from aucomana.utils.reactions import Reactions

g_file = 'Runs/lelsb_run/analysis/group_template.tsv'
comp_dir = 'Runs/run62/analysis/all'
spha = 'Sphacelaria-rigidula_FEMALE'

Auc = AuCoMAna(comp_dir, g_file)
# Auc.group_pathway_completion_comparison(('LR', 'SR', spha), 'Phaeophyceae')
Auc.group_reactions_comparison(('LR', 'SR', spha), 'Phaeophyceae')

# R = Reactions(file_reactions_tsv=comp_dir + '/reactions.tsv')
