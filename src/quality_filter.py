#!/usr/bin/env python3

# query or source (gene) sequence id
# subject or target (reference genome) sequence id
# percentage of identical positions
# alignment length (sequence overlap)
# number of mismatches
# number of gap openings
# start of alignment in query
# end of alignment in query
# start of alignment in subject
# end of alignment in subject
# expect value
# bit score

import argparse

field_definitions:dict[str,int] = {
    "qseqid"  : 0,  
    "sseqid"  : 1,  
    "pident"  : 2,  
    "length"  : 3,  
    "mismatch": 4,  
    "gapopen" : 5,  
    "qstart"  : 6,  
    "qend"    : 7,  
    "sstart"  : 8,  
    "send"    : 9,  
    "evalue"  : 10, 
    "bitscore": 11, 
    "gene1len": 12, 
    "gene2len": 13, 
    "gene1id" : 14, 
    "gene2id" : 15  
}

parser:argparse.ArgumentParser = argparse.ArgumentParser(
    prog="BlastFilter",
)
parser.add_argument("filename")
parser.add_argument("-i", type=int)
parser.add_argument("-c", type=int)

argv = parser.parse_args()

pident_threshold:int = argv.i
coverage_threshold:int = argv.c

blastfile:str = argv.filename

def coverage(alnlen:int, genelen:int) -> float:
    return 100 * alnlen / genelen

if __name__ == "__main__":
    with open(blastfile, "r") as file:
        for line in file:
            print(line.rstrip('\n'))
            break

    with open(blastfile, "r") as file:
        next(file)
        for line in file:
            fields:list[str] = line.split(',')
            pident:float = float(fields[field_definitions["pident"]])
            alnlen:int = int(fields[field_definitions["length"]])
            gene1id:str = fields[field_definitions["gene1id"]]
            gene1len:int = float(field_definitions["gene1id"])
            gene2id:str = fields[field_definitions["gene2id"]]
            gene2len:int = float(field_definitions["gene2id"])

            min_cover:float = min(coverage(alnlen, gene1len), coverage(alnlen, gene2len))

            if pident < pident_threshold or min_cover < coverage_threshold:
                continue

            print(line.rstrip('\n'))
                 

