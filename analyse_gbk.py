import numpy
from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import os


GBK_PATH = "/media/paulinehg/OS/Users/Pauline/Documents/MyDocs/Gbk/short_read/cut_genomes/"
OUT_PATH = "/media/paulinehg/OS/Users/Pauline/Documents/MyDocs/Gbk/info_gbk/"

# TO DO : adapt create info file for boxplot !


def individual_analysis(gbk_path, out_path):
    for dir in os.listdir(gbk_path):

        print(f"Creating table for {dir}")

        file_path = f"{gbk_path}{dir}/{dir}.gbk"

        df = pd.DataFrame(columns=("id_contig", "#gene", "#CDS", "#b"))
        for record in SeqIO.parse(file_path, "genbank"):
            all_features = [feature.type for feature in record.features]
            features_nb = {"gene": 0, "CDS": 0}
            features_count = dict(Counter(all_features))
            features_nb.update(features_count)
            row = pd.DataFrame({"id_contig": record.id, "#gene": features_nb["gene"],
                                "#CDS": features_nb["CDS"], "#b": len(record.seq)}, index=[0])
            df = pd.concat([df, row], ignore_index=True)

        df.to_csv(f"{out_path}{dir}.tsv", sep="\t", index=False)


def group_analysis(out_path):
    df = pd.DataFrame(columns=("species", "#contigs", "#gene", "#CDS", "#b", "#contigs no gene", "contigs at least 1 gene (%)"))
    for gbk_analysis_file in os.listdir(out_path):
        gbk_analysis_df = pd.read_csv(f"{out_path}{gbk_analysis_file}", sep="\t", header=0)
        species = gbk_analysis_file.split(".")[0]
        sum_contig = gbk_analysis_df.shape[0]
        sum_gene = sum(gbk_analysis_df["#gene"])
        sum_cds = sum(gbk_analysis_df["#CDS"])
        sum_b = sum(gbk_analysis_df["#b"])
        sum_contig_no_gene = list(gbk_analysis_df["#gene"]).count(0)
        cont_at_least_1_gene = round(100 - (sum_contig_no_gene / sum_contig) * 100, 2)
        row = pd.DataFrame({"species": species, "#contigs": sum_contig, "#gene": sum_gene, "#CDS": sum_cds, "#b": sum_b,
                            "#contigs no gene": sum_contig_no_gene,
                            "contigs at least 1 gene (%)": cont_at_least_1_gene}, index=[0])
        df = pd.concat([df, row], ignore_index=True)
    df.to_csv(f"{out_path}{'all_species'}.tsv", sep="\t", index=False)


def boxplot_gbk_groups(group1_file, group2_file, outpath):
    col = ("#contigs", "#genes", "contigs at least 1 gene (%)")
    df_g1 = pd.read_csv(group1_file, sep="\t", index_col=0)
    df_g2 = pd.read_csv(group2_file, sep="\t", index_col=0)
    fig = plt.figure(figsize=[15.8, 7.6])
    i = 1
    titles = ("a. Number of contigs depending\nof short or long read sequencing",
              "b. Number of genes depending\nof short or long read sequencing",
              "c. Percentage of contigs with\nat least 1 genes depending\n"
              "of short or long read sequencing")
    plt.rcParams.update({'font.size': 15})

    for prop in col:
        axs = fig.add_subplot(1, 3, i)
        axs.boxplot([df_g1[prop], df_g2[prop]], labels=("Long Read", "Short Read"))
        axs.set_ylabel(prop)
        axs.set_title(titles[i-1])
        i += 1
        print(prop, numpy.sqrt(numpy.var(df_g1[prop])),
              numpy.sqrt(numpy.var(df_g2[prop])))
        print(prop, numpy.mean(df_g1[prop]),
              numpy.mean(df_g2[prop]))
        print(prop, numpy.median(df_g1[prop]),
              numpy.median(df_g2[prop]))
    fig.tight_layout()
    plt.savefig(outpath, dpi=200)


individual_analysis(GBK_PATH, OUT_PATH)
# group_analysis(OUT_PATH)
# boxplot_gbk_groups("data/info_gbk/long_read.tsv",
#                    "data/info_gbk/short_read.tsv",
#                    "data/info_gbk/boxplot_info_gbk.png")
