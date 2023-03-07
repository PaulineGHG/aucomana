from aucomana.aucomana import AuCoMAna

g_file = 'Runs/run62/analysis/group_template.tsv'
comp_dir = 'Runs/run62/analysis/all'
spha = 'Sphacelaria-rigidula_FEMALE'
orders = ['Fucales', 'Ectocarpales', 'Discosporangiales', 'Desmarestiales', 'Tilopteridales',
          'Ralfsiales', 'Laminariales', 'Sphacelariales', 'Dictyotales', 'Outgroup']

Auc = AuCoMAna(comp_dir, g_file)

# Auc.group_pathway_completion_comparison(('LR', 'SR', spha), 'Phaeophyceae')
# Auc.group_reactions_comparison(('LR', 'SR', spha), 'Phaeophyceae')
Auc.group_supervenn_rxn(orders, fig_size=(64, 48), output='supervenn_run62.png')

