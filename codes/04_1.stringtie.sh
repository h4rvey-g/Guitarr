#!/bin/bash
mkdir -p ./data/04.Tx_quantification
mkdir -p ./data/04.Tx_quantification/01.Stringtie
# use stringtie to assemble the transcripts of bam files in ./data/03.Align, use data/reference/GRCh38.gtf

get_quantification() {
    bam=$1
    # if bam contains ""SIRV"" then use SIRV gtf
    if [[ $bam == *"SIRV"* ]]; then
        gtf=data/reference/GRCh38_SIRV.gtf
    else
        gtf=./data/reference/GRCh38_chr.gtf
    fi
    stringtie $bam -e -p 5 -G $gtf -o ./data/04.Tx_quantification/01.Stringtie/$(basename $bam .bam).gtf -A ./data/04.Tx_quantification/01.Stringtie/$(basename $bam .bam).tab -B
}

export -f get_quantification

find ./data/03.Align -name "*Aligned.sortedByCoord.out.bam" | parallel --progress --keep-order --line-buffer get_quantification
