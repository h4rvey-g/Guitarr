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
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz" -O ./data/reference/GRCh37.gtf.gz && gunzip ./data/reference/GRCh37.gtf.gz
wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz" -O ./data/reference/GRCh37.fa.gz && gunzip ./data/reference/GRCh37.fa.gz
mkdir -p ./data/reference/3rd_gen
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192955&format=file&file=GSE192955%5FPC3E%5FGS689%5FHEK293T%5F1D%5FcDNA%5FN2%5FR0%5Fabundance%2Eesp%2Etxt%2Egz" -O ./data/reference/3rd_gen/3rd_gen.txt.gz && gunzip ./data/reference/3rd_gen/3rd_gen.txt.gz
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192955&format=file&file=GSE192955%5FPC3E%5F1D%5FcDNA%5FSIRV%5FN2%5FR0%5Fabundance%2Eesp%2Etxt%2Egz" -O ./data/reference/3rd_gen/SIRV.txt.gz && gunzip ./data/reference/3rd_gen/SIRV.txt.gz
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192955&format=file&file=GSE192955%5FPC3E%5FGS689%5FHEK293T%5F1D%5FcDNA%5FN2%5FR0%5Fupdated%2Egtf%2Egz" -O ./data/reference/3rd_gen/3rd_gen.gtf.gz && gunzip ./data/reference/3rd_gen/3rd_gen.gtf.gz
wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE192955&format=file&file=GSE192955%5FPC3E%5FGS689%5FHEK293T%5F1D%5FcDNA%5FN2%5FR0%5Fupdated%2Egtf%2Egz" -O ./data/reference/3rd_gen/additional.gtf.gz && gunzip ./data/reference/3rd_gen/additional.gtf.gz
seqkit grep -r -p '^chr([0-9][0-9]?|X|Y)$' ./data/reference/GRCh38.fa >./data/reference/GRCh38_chr.fa
seqkit concat ./data/reference/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta_170612a.fasta ./data/reference/GRCh38_chr.fa -o ./data/reference/SIRV_combined.fa -f
samtools faidx ./data/reference/GRCh38_chr.fa
samtools faidx ./data/reference/SIRV_combined.fa
agat_sp_merge_annotations.pl --gff ./data/reference/GRCh38_chr.gtf --gff data/reference/SIRV_Set1_Norm_Sequences_20210507/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf -o ./data/reference/GRCh38_SIRV.gtf
stringtie --merge -G ./data/reference/GRCh37.gtf ./data/reference/3rd_gen/additional.gtf -o ./data/reference/3rd_gen/3rd_gen.gtf -m 0 -F 0 -T 0
