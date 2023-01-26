import os

from aucomana import analysis

pruns = "/home/phamongi/Documents/Runs"
pst = os.getcwd()
orgf = "/home/phamongi/Documents/Perso_analysis_runs/data/species_group.tsv"

A = analysis.Analysis(pruns, pst, orgf)
A.rename_padmet_id("run04")