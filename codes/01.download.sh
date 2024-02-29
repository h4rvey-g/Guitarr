#!/usr/bin/env bash
#!/bin/bash
mamba activate fastq-dl
# download ERR4330671, ERR4330672, ERR4330673 using fastq-dl
SIRV_accession="ERR4330671,ERR4330672,ERR4330673"
# use while read loop to download the fastq files
echo $SIRV_accession | tr ',' '\n' | while read acc; do
    fastq-dl --accession $acc -o ./data/01.Raw_data/
done
# download other 2 data
fastq-dl --accession PRJNA438990 -o ./data/01.Raw_data
awk -F'\t' 'NR>1 {print $1,$2"_"$3}' "data/01.Raw_data/samples.tsv" | while read old new; do
    mv "data/01.Raw_data/${old}_1.fastq.gz" "data/01.Raw_data/${new}_1.fastq.gz"
    mv "data/01.Raw_data/${old}_2.fastq.gz" "data/01.Raw_data/${new}_2.fastq.gz"
done
