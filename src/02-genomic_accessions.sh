awk '$1 ~ ">" {split($3, a, ":"); chr=a[3]; if (chr != "Mt" && chr != "Pt") print substr($1, 2)}' data/Oryza_sativa.IRGSP-1.0.pep.all.fa > output/genomic_accessions.txt
