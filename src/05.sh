python3 src/quality_filter.py output/complete_blast_result.csv -i 30 -c 30 > output/blast_low_stringency.csv
python3 src/quality_filter.py output/complete_blast_result.csv -i 50 -c 40 > output/blast_high_stringency.csv
