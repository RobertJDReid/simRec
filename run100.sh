# First run — write header
python simRec.py genome_chrII.csv --n-gen 500 --p-rec 1e-7 --events > results.tsv

# Subsequent runs — append without header
for i in $(seq 2 100); do
    python simRec.py genome_chrII.csv --n-gen 500 --p-rec 1e-7 --events-nh >> results.tsv
done