import pandas as pd
import time

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

spl2 = ['Dictyota_dichotoma_m']


def get_dict(sp):
    dic = f"data/run01_studied_organism/{sp}/{sp}_dict.csv"
    assodict = {}
    with open(dic, 'r') as d:
        for l in d:
            li = l.split()
            assodict[li[1]] = li[0]
    return assodict


def edit_padmet(spl):
    for sp in spl:
        assodict = get_dict(sp)
        file = f"data/run01_studied_organism/{sp}/{sp}.padmet"
        out = f"data/run01_studied_organism/2PADMETs/{sp}.padmet"
        with open(file, 'r') as f:
            with open(out, 'w') as o:
                nb = 0
                for li in f:
                    if "c_organism" in li:
                        for el in assodict:
                            if el in li:
                                li = li.replace(el, assodict[el])
                    # for i in range(len(li)):
                    #     if "c_organism" in li[i]:
                    #         for el in assodict:
                    #             if el in li[i]:
                    #                 li[i] = li[i].replace(el, assodict[el])
                        # if ' or ' in li[i]:
                        #     li2 = li[i].split(' or ')
                        #     for j in range(len(li2)):
                        #         el = li2[j].strip()[1:-1]
                        #         if el in assodict:
                        #             li2[j] = f"({assodict[el]})"
                        #     li[i] = " or ".join(li2)
                        # if li[i] in assodict:
                        #     li[i] = assodict[li[i]]
                        # elif li[i].split("_names")[0] in assodict:
                        #     li[i] = assodict[li[i].split("_names")[0]] + "_names"
                    o.write(li)
        print(nb)


start = time.time()
for sp in spl2:
    dic = get_dict(sp)
    file = f"data/run01_studied_organism/{sp}/{sp}.padmet"
    out = f"data/run01_studied_organism/2PADMETs/{sp}.padmet"
    with open(file, "r") as in_file, open(out, "w") as out_file:
        for l in in_file:
            l2 = l.translate(dic)
            out_file.write(l2)
end = time.time()
print(end - start)

# 1 : 443.7
# 2 : 255.9
# 3 : 235.9
