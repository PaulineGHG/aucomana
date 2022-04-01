from Bio import SeqIO
from collections import Counter
import pandas as pd
import os


gbk_path = "/media/paulinehg/OS/Users/Pauline/Documents/MyDocs/Gbk/short_read/uncut_genomes/"
out_path = "/media/paulinehg/OS/Users/Pauline/Documents/MyDocs/Gbk/info_gbk/"

for dir in os.listdir(gbk_path):

    print(f"Creating table for {dir}")

    file = f"{dir}/{dir}.gbk"
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

    df.to_csv(f"{out_path}{dir}.tsv", sep="\t")
