import pandas as pd
from itertools import combinations
import numpy as np

genes = pd.read_csv("output/master_gene_table.tsv", sep='\t')
family_groups = genes.groupby("hsfamily")
# next(family_groups)

with open("output/gene_combinations.tsv", "w") as o:
    o.write(f"family\tseqid1\tsequid2\n")
    
    for f, g in family_groups:
        if f == 0:
            continue
        for c in combinations(g["seqid"], 2):
            o.write(f"{f}\t{c[0]}\t{c[1]}\n")
