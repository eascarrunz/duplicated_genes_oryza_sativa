from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("output/seq_divergence.tsv", sep='\t')
tags = pd.read_csv("output/hstags1.tsv", sep='\t')
genes = pd.read_csv("output/master_gene_table.tsv", sep='\t')

Ks = df["dS"]

seq_to_gene = {s: g for s,g in zip(genes["seqid"], genes["geneID"])}

idseq_tags = set()
for x in tags["genes"].values:
    for seq in x.split(','):
        idseq_tags.add(seq)
        
df["intag"] = False
df["intag"] = \
    np.array(list(map(lambda x: seq_to_gene[x] in idseq_tags, df["seqid1"])))
df["intag"] =  df["intag"] & \
    np.array(list(map(lambda x: seq_to_gene[x] in idseq_tags, df["seqid2"])))

fig, axes = plt.subplots(2, 1, sharex=True, layout="constrained")
Ks_in_tag = Ks[df["intag"]]
axes[0].hist(Ks[Ks < 2.0], bins=20, label="All gene pairs")
axes[0].hist(Ks_in_tag[Ks_in_tag < 2.0], bins=20, label="Gene pairs in TAG")
axes[0].legend()
# axes[0].set_xlabel("Ks")
axes[0].set_ylabel("Numbers of pairs of genes")
axes[0].annotate(
        "A",
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', verticalalignment='top', weight="bold")

Ks_notin_tag = Ks[df["intag"] == False]
axes[1].hist(Ks[Ks < 2.0], bins=20, label="All gene pairs")
axes[1].hist(Ks_notin_tag[Ks_notin_tag < 2.0], bins=20, color="cyan", label="Non-TAG gene pairs")
axes[1].legend()
axes[1].set_xlabel("Ks")
axes[1].set_ylabel("Numbers of pairs of genes")
axes[1].annotate(
        "B",
        xy=(0, 1), xycoords='axes fraction',
        xytext=(+0.5, -0.5), textcoords='offset fontsize',
        fontsize='large', verticalalignment='top', weight="bold")

fig.savefig("output/Ks_histograms.png", dpi=200)
