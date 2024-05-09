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
    Guitar -i $bam -a $gtf -g $genome -o ./data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute/ -p 40 -c -r
}
# test Guitar by running it on one bam file
export -f impute
bams=$(find ./data/03.Align -name "*sortedByCoord.out.bam")
parallel -S xu_cu21,xu_cu23,: --workdir /data0/work/guozhonghao/Guitarr -j1 --progress --line-buffer --verbose --env impute impute ::: $bams
# run stringtie to get quantification of imputed bam files
get_quantification() {
    bam=$1
    # if bam contains ""SIRV"" then use SIRV gtf
    if [[ $bam == *"SIRV"* ]]; then
        gtf=data/reference/GRCh38_SIRV.gtf
    else
        gtf=./data/reference/GRCh38_chr.gtf
    fi
    sample_name=$(basename $bam Aligned.sortedByCoord.out.bam)
    output_dir="./data/04.Tx_quantification/02.Stringtie_with_impute/$sample_name"
    mkdir -p "$output_dir"
    stringtie $bam -e -p 5 -G $gtf -o "$output_dir/$sample_name.gtf" -A "$output_dir/$sample_name.tab" -B
}
export -f get_quantification
find ./data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute -name "*imputed.bam" | parallel --progress --keep-order --line-buffer get_quantification
