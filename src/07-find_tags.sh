#!/usr/bin/env
python3 src/tag_finder.py output/master_gene_table.tsv -s 0 -f hsfamily > output/hstags0.tsv
python3 src/tag_finder.py output/master_gene_table.tsv -s 0 -f lsfamily > output/lstags0.tsv
python3 src/tag_finder.py output/master_gene_table.tsv -s 1 -f hsfamily > output/hstags1.tsv
python3 src/tag_finder.py output/master_gene_table.tsv -s 1 -f lsfamily > output/lstags1.tsv
python3 src/tag_finder.py output/master_gene_table.tsv -s 2 -f hsfamily > output/hstags2.tsv
python3 src/tag_finder.py output/master_gene_table.tsv -s 2 -f lsfamily > output/lstags2.tsv
