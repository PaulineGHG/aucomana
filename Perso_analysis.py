from aucomana.aucomana import GroupsAnalysis, SequencesAnalysis

g_file = 'Runs/run62/analysis/group_template.tsv'
comp_dir = 'Runs/run62/analysis/all'
spha = 'Sphacelaria-rigidula_FEMALE'
orders = ['Fucales', 'Ectocarpales', 'Discosporangiales', 'Desmarestiales', 'Tilopteridales',
          'Ralfsiales', 'Laminariales', 'Sphacelariales', 'Dictyotales', 'Outgroup']
sp_seq_dir = 'Runs/run62/studied_organisms'

# GA = GroupsAnalysis(comp_dir, g_file)

# Auc.group_pathway_completion_comparison(('LR', 'SR', spha), 'Phaeophyceae')
# GA.group_reactions_comparison(('LR', 'SR', spha), 'Phaeophyceae')
# GA.group_supervenn_rxn(groups_comp=orders, fig_size=(64, 48), output='supervenn_run62.png')

SA = SequencesAnalysis(sp_seq_dir)

print(SA.sp_genome)
print(SA.sp_proteome)

