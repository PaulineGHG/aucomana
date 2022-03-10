import ahocorasick
import time


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


def get_dict(sp):
    dic = f"data/run01_studied_organism/{sp}/{sp}_dict.csv"
    assodict = {}
    with open(dic, 'r') as d:
        for l in d:
            li = l.split()
            assodict[li[1]] = li[0]
    return assodict


spl = ['Ectocarpus_fasciculatus_m',
       'Desmarestia_herbacea_m', 'Ectocarpus_siliculosus_m', 'Chordaria_linearis',
       'Scytosiphon_promiscuus_MALE',
       'Pleurocladia_lacustris',
       'Ectocarpus_crouaniorum_m',
       'Fucus_serratus_MALE', 'Saccharina_latissima_FEMALE',
       'Schizocladia_ischiensis', 'Dictyota_dichotoma_m',
       'Porterinema_fluviatile', 'Laminarionema_elsbetiae']

spl2 = ['Dictyota_dichotoma_m']

start = time.time()

for sp in spl2:
    automaton = make_automaton(get_dict(sp))
    apply_automaton(automaton, f"data/run01_studied_organism/{sp}/{sp}.padmet",
                    f"data/run01_studied_organism/2PADMETs/test_{sp}.padmet")

end = time.time()
print(end - start)

