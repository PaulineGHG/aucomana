# Reactions_loss

## data

Contains data files used as input :

> Lelsb_losses.ods : file of reactions found manualy
> Proj tut 1.ipynb : jupyter file written by L2 students 

### reactions files (AuCoMe output) :

> reactions_phaeo.tsv : original output of run A0 (with "present" and empty fields)
> runA0_reactions.tsv : run Alexandre 2020 (modified with 1 and 0 fields)
> runA1_reactions.tsv : run Pauline like Alexandre's exactly
> runA2_reactions.tsv : run Pauline like Alexandre's improved
> run40_reactions.tsv : run Jeanne 40 + 7
> run01_reactions.tsv : run Pauline best quality data phaeoxplorer

## output

Contains output files from executed commands

## reactions.py

Contains Reactions class

## algae_project.py

Contains command to analyse files :

> Creation of a Reactions object instance : R = Reactions(input_file, [species_list], out(optional))

