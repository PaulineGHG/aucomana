from Bio import SeqIO
from collections import Counter
import pandas as pd
import os


gbk_path = ""
file = "Chordaria_linearis.gbk"

file_path = f"{gbk_path}{file}"

df = pd.DataFrame(columns=("id_contig", "#gene", "#CDS", "#b"))
for record in SeqIO.parse(file_path, "genbank"):
    all_features = [feature.type for feature in record.features]
    features_nb = {"gene": 0, "CDS": 0}
    features_count = dict(Counter(all_features))
    features_nb.update(features_count)
    row = pd.DataFrame({"id_contig": record.id, "#gene": features_nb["gene"],
                        "#CDS": features_nb["CDS"], "#b": len(record.seq)}, index=[0])
    df = pd.concat([df, row], ignore_index=True)

df.to_csv(f"{file.split('.')[0]}.tsv", sep="\t")
