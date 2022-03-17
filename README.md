# Reactions_loss

## data

Contains data files used as input :

- Lelsb_losses.ods : file of reactions found manualy
- Proj tut 1.ipynb : jupyter file written by L2 students 
- Ref tree : nexus tree for original phylogeny comparison
- species_group.tsv : file informations of algae (brown, diatoms, etc) and (LR, SR, PUB)

### reactions_data (AuCoMe output) :

- reactions_phaeo.tsv : original output of run A0 (with "present" and empty fields)
- runA0_reactions.tsv : run Alexandre 2020 (modified with 1 and 0 fields)
- runA1_reactions.tsv : run Pauline like Alexandre's exactly
- runA2_reactions.tsv : run Pauline like Alexandre's improved
- run40_reactions.tsv : run Jeanne 40 + 7
- run01_reactions.tsv : run Pauline best quality data phaeoxplorer
- run02_reactions.tsv : run Pauline long read only + outgroup
- run03_reactions.tsv : run Pauline long read + L.Elsbetiae + outgroups
- run01b_reactions.tsv : run Pauline like run01 with corrected ID

## output

### common_reac

Contains output files txt + json indicating common reactions between several runs

### gene_assoc

Contains output files for each reaction in fasta containing each protein sequence for each species

### cut_reactions_data

reactions.tsv files cut according to selected reactions

### dendrogram_comp

Contains directory for runs containing dendrograms (pvclust + dend + dend_comp)

## R_script

Contains R script used to create dendrograms

## reactions.py

Contains Reactions class

## algae_project.py

Contains command to analyse files :

- Creation of a Reactions object instance :

```R = Reactions(input_file, [species_list], [out])```
- Get the reactions lost : for all species : 

```R.reactions_loss / for a specific species : R.reactions_loss[species_name]```
- Get common reactions between different instance : 

```Reactions.get_common_reactions(list_of-instances, species_name, [output_file])```
  - ```output_file=None``` --> returns the result
  - ```output_file=JSON``` --> create .json output file
  - ```output_file=TXT```  --> create .txt  output file

