import pandas as pd

# file = "reactions_phaeo.tsv"
#
# data = pd.read_csv(file, sep="\t", header=0)
# for i in range(16):
#     for j in range(data.shape[0]):
#         if data.loc[j][i+1] == "present":
#             data.loc[j][i + 1] = 1
#         else:
#             data.loc[j][i + 1] = 0
#
# data.to_csv("reactions_phaeo2.tsv", sep="\t"

spl = ['Ectocarpus_fasciculatus_m',
       'Desmarestia_herbacea_m', 'Ectocarpus_siliculosus_m', 'Chordaria_linearis',
       'Scytosiphon_promiscuus_MALE',
       'Pleurocladia_lacustris',
       'Ectocarpus_crouaniorum_m',
       'Fucus_serratus_MALE', 'Saccharina_latissima_FEMALE',
       'Schizocladia_ischiensis', 'Dictyota_dichotoma_m',
       'Porterinema_fluviatile', 'Laminarionema_elsbetiae']

spl2 = ['Chordaria_linearis']


def edit_padmet(spl, path=None):
    for sp in spl:
        dic = f"data/run01_studied_organism/{sp}/{sp}_dict.csv"
        assodict = {}
        with open(dic, 'r') as d:
            for l in d:
                li = l.split()
                assodict[li[1]] = li[0]
        file = f"data/run01_studied_organism/{sp}/{sp}.padmet"
        out = f"data/run01_studied_organism/2PADMETs/2_{sp}.padmet"
        with open(file, 'r') as f:
            with open(out, 'w') as o:
                nb = 0
                for li in f:
                    li = li.split("\n")[0].split("\t")
                    for i in range(len(li)):
                        if li[i] in assodict:
                            li[i] = assodict[li[i]]
                            nb += 1
                        elif li[i].split("_names")[0] in assodict:
                            li[i] = assodict[li[i].split("_names")[0]] + "_names"
                    o.write("\t".join(li) + "\n")
        print(nb)


edit_padmet(spl)
