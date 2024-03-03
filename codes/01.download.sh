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
seqkit grep -r -p '^chr([0-9][0-9]?|X|Y)$' ./data/reference/GRCh38.fa >./data/reference/GRCh38_chr.fa
seqkit concat ./data/reference/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta_170612a.fasta ./data/reference/GRCh38_chr.fa -o ./data/reference/SIRV_combined.fa -f
samtools faidx ./data/reference/GRCh38_chr.fa
samtools faidx ./data/reference/SIRV_combined.fa
agat_sp_merge_annotations.pl --gff ./data/reference/GRCh38_chr.gtf --gff data/reference/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf -o ./data/reference/GRCh38_SIRV.gtf
