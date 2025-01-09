#!/usr/bin/env python3
import pandas as pd
import argparse
import sys

parser:argparse.ArgumentParser = argparse.ArgumentParser(
    prog="tag_finder.py",
    description="Finds TAGs in a genome based on an input data frame"
)
parser.add_argument("filename", help="Genome data frame file. It must include columns called \"chromosome\", \"start\", \"end\", \"geneID\" and another column with the gene family")
parser.add_argument("-f", type=str, help="Name of the family column")
parser.add_argument("-s", type=int, default=1, help="Maximum number of spacers in a TAG")

def create_tag(chromosome, family):
    return {
        "chromosome": chromosome,
        "family": family,
        "start": 0,
        "end": 0,
        "size": 0,
        "genes": []
    }


def add_gene(tag, row):
    tag["end"] = row["end"]
    tag["genes"].append(row["geneID"])
    tag["size"] += 1
    
    return tag


def checkout_tag(df, tag):
    new_row = tag.copy()
    new_row["genes"] = ','.join(tag["genes"])
    if tag["size"] > 1:
        df.loc[len(df.index)] = new_row
        
    return df, create_tag(tag["chromosome"], tag["family"])

def scan_chromosomes_for_tags(df, famcol, spacer_tolerance):
    group_iterator = genes.groupby(["chromosome", famcol])
    tag_df = pd.DataFrame(columns=["chromosome", "family", "start", "end", "size", "genes"])

    for indices, sdf in group_iterator:
        chromosome, family = indices
        if family == 0:
            continue
        current_tag = create_tag(chromosome, family)
        first_row = sdf.iloc[0]
        current_tag["start"] = first_row["start"]
        current_tag = add_gene(current_tag, first_row)
        prev = first_row.name
        for i, row in sdf[1:].iterrows():
            if i - prev <= spacer_tolerance:
                current_tag = add_gene(current_tag, row)
            else:
                tag_df, current_tag = checkout_tag(tag_df, current_tag)
                current_tag["start"] = row["start"]
                current_tag = add_gene(current_tag, row)
            prev = i
        checkout_tag(tag_df, current_tag)

    return tag_df

if __name__ == "__main__":
    argv = parser.parse_args()

    genes = pd.read_csv(argv.filename, sep='\t')
    genes.sort_values(["chromosome", "start"], inplace=True)
    
    df_tags = scan_chromosomes_for_tags(genes, argv.f, argv.s + 1)
    df_tags.sort_values(["chromosome", "start"], inplace=True)
    
    df_tags.to_csv(sys.stdout, sep='\t', index=False)
