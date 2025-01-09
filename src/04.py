import pandas as pd

# Read in protein lengths file
df_lengths = pd.read_csv("output/protein_lengths.txt", header=None, sep=" ")
df_lengths = df_lengths.rename(columns={0:"seqid", 1:"length"})

# Read in genomic accession ids
df_genacc = pd.read_csv("output/genomic_accessions.txt", header=None, sep=" ")
df_genacc = df_genacc.rename(columns={0:"GenomicAccession"})

df_blast = pd.read_csv("data/Oryza_sativa_Blastp_longIsoforme", header=None, sep='\t')
df_blast = df_blast.rename(columns={0:"qseqid", 1:"sseqid", 2:"pident", 3:"length", 4:"mismatch", 5:"gapopen", 6:"qstart", 7:"qend", 8:"sstart", 9:"send", 10:"evalue", 11:"bitscore"})

df_list = pd.read_csv("output/Oryza_sativa_liste", sep='\t')

nrows1 = df_blast.shape[0]

df_blast_filtered = df_blast.loc[df_blast["qseqid"].isin(df_genacc["GenomicAccession"]) & df_blast["sseqid"].isin(df_genacc["GenomicAccession"])]
nrows2 = df_blast_filtered.shape[0]

# Check how many rows were removed
nrows2 - nrows1

df_join = df_blast_filtered.join(df_lengths.set_index("seqid"), on="qseqid", rsuffix="q")
df_join = df_join.join(df_lengths.set_index("seqid"), on="sseqid", rsuffix="s")
df_join = df_join.rename(columns={"lengthq":"gene1len", "lengths":"gene2len"})

df_join = df_join.join(df_list.set_index("proteinID"), on="qseqid", rsuffix="q")
df_join = df_join.join(df_list.set_index("proteinID"), on="sseqid", rsuffix="s")
df_join = df_join.rename(columns={"geneID":"gene1id", "geneIDs":"gene2id"})

df_join.to_csv("output/complete_blast_result.csv", index=False)
