import pandas as pd
import numpy as np

def families_from_family_list_file(df:pd.DataFrame, filename:str, famcol:str) -> pd.DataFrame:
    protein_ids:list[str] = []
    with open(filename) as file:
        for i, row in enumerate(file):
            protein_ids = row.strip().split('\t')
            # Assert that none of the genes in this family had been assigned to another family
            assert (df.loc[protein_ids, famcol] == 0).all()
            
            df.loc[protein_ids, famcol] = i + 1
    
    return df

if __name__ == "__main__":
    df_genes:pd.DataFrame = pd.read_csv("output/genomic_positions.tsv", sep='\t')
    df_genes.loc[:, "hsfamily"] = 0
    df_genes.loc[:, "lsfamily"] = 0
    df_genes.set_index("seqid", inplace=True)

    df_genes = families_from_family_list_file(df_genes, "output/mcl_high_stringency.tabular", "hsfamily")
    df_genes = families_from_family_list_file(df_genes, "output/mcl_low_stringency.tabular", "lsfamily")

    df_genes.sort_values(by=["chromosome", "start"], inplace=True)
    
    # Discard isoforms that (1) do not belong to any family or (2) are not the longest isoform of the cds
    discarded_isoforms = []
    for gene in df_genes.groupby("geneID"):
        isoforms = gene[1]
        if isoforms.shape[0] > 1:
            selected_isoform = isoforms.index[(isoforms["hsfamily"] > 0) | (isoforms["lsfamily"] > 0)]
            if len(selected_isoform) < 1:
                longest_isoform = np.argmax(isoforms["end"] - isoforms["start"])
                selected_isoform = isoforms.index[longest_isoform]
            discarded_isoforms += np.setdiff1d(isoforms.index, selected_isoform).tolist()

    df_genes.drop(discarded_isoforms, inplace=True)

    df_genes.to_csv("output/master_gene_table.tsv", sep='\t')
