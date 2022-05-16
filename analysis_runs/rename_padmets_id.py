import os.path
import ahocorasick


def make_automaton(gene_to_replace):
    """Build an Aho-Corasick automaton from a dictionary file and return
    it. The lines in the dictionary file must consist of a key and a
    value separated by a space.
    """
    automaton = ahocorasick.Automaton()
    for key, value in gene_to_replace.items():
        automaton.add_word(key, (key, value))
    automaton.make_automaton()
    return automaton


def apply_automaton(automaton, input_filename, output_filename):
    """Apply an Aho-Corasick automaton to an input file, replacing the
    first occurrence of a key in each line with the corresponding
    value, and writing the result to the output file.

    """
    with open(input_filename) as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            for e, (key, value) in automaton.iter(line):
                line = line.replace(key, value)
            outfile.write(line)


def get_dict(run, sp, assodict, path_runs):
    dic = os.path.join(path_runs, run, "studied_organisms", sp, f"{sp}_dict.csv")
    if os.path.exists(dic):
        with open(dic, 'r') as d:
            for l in d:
                li = l.split()
                assodict[li[1]] = li[0]
        return assodict
    else:
        return assodict




