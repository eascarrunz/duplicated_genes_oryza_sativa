# Extract protein IDs with lengths
awk '/^>/{if (seq) print id, length(seq); id=substr($1, 2); seq=""} /^[^>]/ {seq = seq $0} END {if (seq) print id, length(seq)}' data/Oryza_sativa.IRGSP-1.0.pep.all.fa > output/protein_lengths.txt
