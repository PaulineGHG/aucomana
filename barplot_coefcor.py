import pandas as pd
import matplotlib.pyplot as plt


RNX = ["output_data/dendro_tanglegrams/run04/rnx_all_sp/rnx_all_sp_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/rnx_all_sp_out05perc/rnx_all_sp_out05perc_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/rnx_all_sp_out10perc/rnx_all_sp_out10perc_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/rnx_all_sp_out20perc/rnx_all_sp_out20perc_similarity_indicators.tsv"]

col = ["Correlation cophenetic", "Correlation bakers gamma"]
ind = ["C", "5%", "10%", "20%"]
df = pd.DataFrame(columns=col, index=ind)
for i in range(4):
    tb = pd.read_csv(RNX[i], sep="\t")
    df.loc[ind[i], col[0]] = round(tb.loc[0, col[0]], 2)
    df.loc[ind[i], col[1]] = round(tb.loc[0, col[1]], 2)
print(df)
ax = df.plot.bar()

for container in ax.containers:
    ax.bar_label(container)
ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
ax.set_xlabel("Minimal completion percentage of reactions between species")
ax.set_title("a. Correlations depending on minimal completion percentage of reactions between species")
plt.tight_layout()
plt.show()


PWC = ["output_data/dendro_tanglegrams/run04/rnx_all_sp/rnx_all_sp_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_100_all_sp_allrxn/pw_100_all_sp_allrxn_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_80_all_sp_allrxn/pw_80_all_sp_allrxn_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_75_all_sp_allrxn/pw_75_all_sp_allrxn_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_60_all_sp_allrxn/pw_60_all_sp_allrxn_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_50_all_sp_allrxn/pw_50_all_sp_allrxn_similarity_indicators.tsv"]

col = ["Correlation cophenetic", "Correlation bakers gamma"]
ind = ["C", "100%", "80%", "75%", "60%", "50%"]
df = pd.DataFrame(columns=col, index=ind)
for i in range(6):
    tb = pd.read_csv(PWC[i], sep="\t")
    df.loc[ind[i], col[0]] = round(tb.loc[0, col[0]], 2)
    df.loc[ind[i], col[1]] = round(tb.loc[0, col[1]], 2)
print(df)
ax = df.plot.bar()
for container in ax.containers:
    ax.bar_label(container)
ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
ax.set_xlabel("Pathway completion rate")
ax.set_title("b. Correlations depending on pathway completion rate")
plt.tight_layout()
plt.show()

PWM = ["output_data/dendro_tanglegrams/run04/rnx_all_sp/rnx_all_sp_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_100_all_sp_allrxn/pw_100_all_sp_allrxn_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_100_all_sp_2min/pw_100_all_sp_2min_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_100_all_sp_3min/pw_100_all_sp_3min_similarity_indicators.tsv",
       "output_data/dendro_tanglegrams/run04/pw_100_all_sp_4min/pw_100_all_sp_4min_similarity_indicators.tsv"]

col = ["Correlation cophenetic", "Correlation bakers gamma"]
ind = ["C", "2", "3", "4"]
df = pd.DataFrame(columns=col, index=ind)
for i in range(4):
    tb = pd.read_csv(PWM[i], sep="\t")
    df.loc[ind[i], col[0]] = round(tb.loc[0, col[0]], 2)
    df.loc[ind[i], col[1]] = round(tb.loc[0, col[1]], 2)
print(df)
ax = df.plot.bar()
for container in ax.containers:
    ax.bar_label(container)
ax.legend(bbox_to_anchor=(1.0, 1), loc='upper left')
ax.set_xlabel("Minimal total number of reactions in pathway")
ax.set_title("c. Correlations depending on minimal total number of reactions in pathway")
plt.tight_layout()
plt.show()