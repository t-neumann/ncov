source activate nextstrain

./scripts/normalize_gisaid_fasta.sh data/gisaid_cov2020_sequences.fasta data/sequences.fasta

snakemake -p

cd results

../scripts/hCov_consensus.py -f masked.fasta -r ../config/reference.fasta
