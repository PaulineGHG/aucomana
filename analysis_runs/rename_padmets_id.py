import os.path
import ahocorasick
from typing import Dict


def make_automaton(genes_to_replace: Dict[str, str]) -> 'ahocorasick.Automaton':
    """ Creates automaton from the dictionary of renamed genes ID and returns it.

    Parameters
    ----------
    genes_to_replace: Dict[str, str]
        Dictionary associated all the original genes ID to their renamed ID for all the species.

    Returns
    -------
    automaton : ahocorasick.Automaton
        automaton created from the dictionary of genes ID replaced.
    """
    automaton = ahocorasick.Automaton()
    for original_id, renamed_id in genes_to_replace.items():
        automaton.add_word(original_id, (original_id, renamed_id))
    automaton.make_automaton()
    return automaton


def apply_automaton(automaton: 'ahocorasick.Automaton', input_padmet_file: str, output_padmet_file: str):
    """ Apply the automaton to the original padmet file and replaces all the renamed IDs by the originals ID. Write the
    result to a new padmet file.

    Parameters
    ----------
    automaton : ahocorasick.Automaton
        automaton created from the dictionary of genes ID replaced.
    input_padmet_file : str
        Path of the original padmet file.
    output_padmet_file: str
        Path for the new padmet file

    """
    with open(input_padmet_file) as infile, open(output_padmet_file, 'w') as outfile:
        for line in infile:
            for e, (renamed_id, original_id) in automaton.iter(line):
                line = line.replace(renamed_id, original_id)
            outfile.write(line)


def get_dict(run: str, species: str, asso_dict: Dict[str, str], path_runs: str) -> Dict[str, str]:
    """ Complete the assodict parameter (Dictionary associated original genes IDs with the renamed IDs) with the IDs
    of genes associated to the species to consider.

    Parameters
    ----------
    run : str
        ID of the run to consider
    species : str
        Name of the species to consider
    asso_dict : Dict[str, str]
        Dictionary associated original genes IDs with the renamed IDs for species having been treated by the function
        previously
    path_runs : str
        Path of AuCoMe runs

    Returns
    -------
    ass_odict : Dict[str, str]
        Dictionary associated original genes IDs with the renamed IDs for species having been treated by the function
        previously and completed by the IDs of genes associated to the species to consider.
    """
    dic = os.path.join(path_runs, run, "studied_organisms", species, f"{species}_dict.csv")
    if os.path.exists(dic):
        with open(dic, 'r') as d:
            for line in d:
                line = line.split()
                asso_dict[line[1]] = line[0]
        return asso_dict
    else:
        return asso_dict




