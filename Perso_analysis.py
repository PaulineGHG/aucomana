import os.path

from aucomana.aucomana import GroupsAnalysis, SequencesAnalysis
from aucomana.utils.utils import get_grp_set

g_file = 'Runs/run62/analysis/group_template.tsv'
comp_dir = 'Runs/run62/analysis/all'
spha = 'Sphacelaria-rigidula_FEMALE'
orders = ['Fucales', 'Ectocarpales', 'Discosporangiales', 'Desmarestiales', 'Tilopteridales',
          'Ralfsiales', 'Laminariales', 'Sphacelariales', 'Dictyotales', 'Outgroup']

sp_seq_dir = 'Runs/run62/studied_organisms'
brown = get_grp_set(g_file, 'Phaeophyceae')

# GA = GroupsAnalysis(comp_dir, g_file)

# Auc.group_pathway_completion_comparison(('LR', 'SR', spha), 'Phaeophyceae')
# GA.group_reactions_comparison(('LR', 'SR', spha), 'Phaeophyceae')
# GA.group_supervenn_rxn(groups_comp=orders, fig_size=(64, 48), output='supervenn_run62.png')

mannitol_pw_rxn = ['FRUCTOKINASE-RXN', 'MANNITOL-1-PHOSPHATASE-RXN',
                   'MANNITOL-2-DEHYDROGENASE-RXN', 'MANNPDEHYDROG-RXN']
align_output = 'Runs/run62/analysis/Multiple_align/PWY-6531'
SA = SequencesAnalysis(sp_seq_dir, comp_dir)
for rxn in mannitol_pw_rxn:
    out = os.path.join(align_output, rxn)
    os.mkdir(out)
    SA.multiple_alignments(rxn, 'Ectocarpus-subulatus', out, list(brown))

