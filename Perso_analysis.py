import os.path

from aucomana.aucomana import GroupsAnalysis, SequencesAnalysis
from aucomana.utils.reactions import Reactions
from aucomana.utils.pathways import Pathways
from aucomana.utils.utils import get_grp_set

g_file = 'Runs/run62/analysis/group_template.tsv'
comp_dir = 'Runs/run16best/analysis/all'
spha = 'Sphacelaria-rigidula_FEMALE'
orders = ['Fucales', 'Ectocarpales', 'Discosporangiales', 'Desmarestiales', 'Tilopteridales',
          'Ralfsiales', 'Laminariales', 'Sphacelariales', 'Dictyotales', 'Outgroup']

sp_seq_dir = 'Runs/run62/studied_organisms'
brown = get_grp_set(g_file, 'Phaeophyceae')
outgr = get_grp_set(g_file, 'Outgroup')

# GA = GroupsAnalysis(comp_dir, g_file)

#GA.group_pathway_completion_comparison(('LR', 'SR', spha), 'Phaeophyceae')
#GA.group_reactions_comparison(('LR', 'SR', spha), 'Phaeophyceae')
#GA.group_supervenn_rxn(groups_comp=orders, fig_size=(64, 48), output='supervenn_run62_nf.png')
# GA.group_supervenn_rxn(groups_comp=orders, fig_size=(64, 48), output='supervenn_run16best.png')

# mannitol_pw_rxn = ['FRUCTOKINASE-RXN', 'MANNITOL-1-PHOSPHATASE-RXN',
#                    'MANNITOL-2-DEHYDROGENASE-RXN', 'MANNPDEHYDROG-RXN']
# align_output = 'Runs/run62/analysis/Multiple_align/PWY-6531'
# SA = SequencesAnalysis(sp_seq_dir, comp_dir)
# for rxn in mannitol_pw_rxn:
#     out = os.path.join(align_output, rxn)
#     os.mkdir(out)
#     SA.multiple_alignments(rxn, 'Ectocarpus-subulatus', out, list(brown))

# A = GroupsAnalysis(group_template=g_file, compare_path=comp_dir)
# RB = Reactions(A.reactions_tsv, species_list=list(brown))
# RO = Reactions(A.reactions_tsv, species_list=list(outgr))
# RB1 = Reactions(A.reactions_tsv, species_list=list(brown), out=1)
#
# print(RB.get_rxn_absent('Laminarionema-elsbetiae', unique=True))
# print(RB.nb_reactions)
# print(RB1.nb_reactions)

R = Reactions('Runs/gf/5_wf_new_biomass/out/reactions.tsv')
P = Pathways('Runs/gf/5_wf_new_biomass/out/pathways.tsv')
print(P.get_pw_complete())