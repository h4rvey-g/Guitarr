#!/bin/bash
index_dir="./data/reference/STARindex"
mkdir -p $index_dir
if [ -z "$(ls -A ${index_dir})" ]; then
    STAR \
        --runMode genomeGenerate \
        --genomeDir ${index_dir} \
        --runThreadN 20 \
        --genomeFastaFiles data/reference/GRCh38.fa \
        --sjdbGTFfile data/reference/GRCh38.gtf \
        --sjdbOverhang 99
fi
index_dir_SIRV="./data/reference/STARindex_SIRV"
mkdir -p $index_dir_SIRV
if [ -z "$(ls -A ${index_dir_SIRV})" ]; then
    STAR \
        --runMode genomeGenerate \
        --genomeDir ${index_dir_SIRV} \
        --runThreadN 20 \
        --genomeFastaFiles data/reference/SIRV_combined.fa \
        --sjdbGTFfile data/reference/GRCh38_SIRV.gtf \
        --sjdbOverhang 97
fi
mkdir -p ./data/03.Align
getAlign() {
    # if sampleID contains SIRV, align to SIRV reference
    if [[ $1 == *"SIRV"* ]]; then
        index_dir="./data/reference/STARindex_SIRV"
    else
        index_dir="./data/reference/STARindex"
    fi
    indir="data/02.Clean_data"
    outdir="data/03.Align"
    sampleID=$(basename $1 _1.fastq.gz)
    # group=${line[2]}
    # echo "Processing: "$sampleID::$group
    STAR \
        --runMode alignReads \
        --twopassMode Basic \
        --genomeDir ${index_dir} \
        --runThreadN 10 \
        --readFilesIn $1 ${1%_1.fastq.gz}_2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix $outdir/$sampleID \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 5 \
        --outFilterType BySJout \
        --quantMode TranscriptomeSAM GeneCounts
}
export -f getAlign
find ./data/02.Clean_data -name "*_1.fastq.gz" | parallel --progress --keep-order --line-buffer getAlign
# find ./data/03.Align -name "*.bam" | parallel --progress --keep-order --line-buffer samtools sort -o {.}.sorted.bam
# date
/data0/apps/anaconda3/bin/multiqc ./data/03.Align -o ./data/03.Align/report
