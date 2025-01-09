#!/usr/bin/env python3
import pandas as pd

fasta_file_name:str = "data/Oryza_sativa.IRGSP-1.0.pep.all.fa"
df_fasta:pd.DataFrame = pd.DataFrame(columns=["seqid", "assembly", "chromosome", "start", "end", "strand"])

row_list:list[str] = []

with open(fasta_file_name, "r") as file:
    for line in file:
        if line[0] != '>':
            continue
        else:
            fields:list[str] = line.split(' ')
            pos_string = fields[2]
            row_list:list[str] = pos_string.split(':')[1:]
            row_seqid:str = fields[0].lstrip('>')
            row_assembly:str = row_list[0]
            if row_list[1] == "Mt" or row_list[1] == "Pt":
                continue
            row_chromosome:int = int(row_list[1])
            row_start:int = int(row_list[2])
            row_end:int = int(row_list[3])
            row_strand:int = int(row_list[4])
            df_fasta.loc[len(df_fasta.index)] = [row_seqid, row_assembly, row_chromosome, row_start, row_end, row_strand]


df_geneids:pd.DataFrame = pd.read_csv("output/Oryza_sativa_liste", sep='\t')
df_fasta = df_fasta.join(df_geneids.set_index("proteinID"), on="seqid", lsuffix="gene")

df_fasta.to_csv("output/genomic_positions.tsv", sep='\t', index=False)
