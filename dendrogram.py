from analysis_runs.utils import *



R04 = "run04"
REACTIONS = get_reactions_inst(runs=[R04])
print(REACTIONS[R04].species_list)
