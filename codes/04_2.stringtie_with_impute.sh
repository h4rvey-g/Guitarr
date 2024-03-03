#!/bin/bash
mkdir -p ./data/04.Tx_quantification/02.Stringtie_with_impute
# use parallel to Guitar gtfidx two gtf files
# if gtf.idx files already exist, skip this step
if [ ! -f ./data/reference/GRCh38_chr.gtf.idx ]; then
    find ./data/reference -name "GRCh38*.gtf" | parallel --progress --keep-order --line-buffer Guitar gtfidx {} -c
fi
# define a function to run Guitar to impute bam files in ./data/03.Align
impute() {
    bam=$1
    # if bam contains ""SIRV"" then use SIRV gtf
    if [[ $bam == *"SIRV"* ]]; then
        gtf=data/reference/GRCh38_SIRV.gtf
        genome=data/reference/SIRV_combined.fa
    else
        gtf=./data/reference/GRCh38_chr.gtf
        genome=./data/reference/GRCh38_chr.fa
    fi
    Guitar -i $bam -a $gtf -g $genome -o ./data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute/ -p 5 -c
}
# test Guitar by running it on one bam file
export -f impute
find ./data/03.Align -name "*Aligned.sortedByCoord.out.bam" | parallel --progress --keep-order --line-buffer impute
# run stringtie to get quantification of imputed bam files
get_quantification() {
    bam=$1
    # if bam contains ""SIRV"" then use SIRV gtf
    if [[ $bam == *"SIRV"* ]]; then
        gtf=data/reference/GRCh38_SIRV.gtf
    else
        gtf=./data/reference/GRCh38_chr.gtf
    fi
    stringtie $bam -e -p 5 -G $gtf -o ./data/04.Tx_quantification/02.Stringtie_with_impute/$(basename $bam .bam).gtf -A ./data/04.Tx_quantification/02.Stringtie_with_impute/$(basename $bam .bam).tab -B
}
export -f get_quantification
find ./data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute -name "*.bam" | parallel --progress --keep-order --line-buffer get_quantification
